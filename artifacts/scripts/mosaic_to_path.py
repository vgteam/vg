#!/usr/bin/env python3
"""Expand a mosaic into one path per thread and check it is an exact walk in the graph.

This is the acceptance test the mosaic format exists to pass. Every other property -- segments
chaining, orientation, walk order, patched gaps, nested excursions -- is a way for this to fail, so
checking it directly subsumes checking them individually.

Inputs are what `vg` can already produce, so this needs no new subcommand:
  paths.gaf   vg paths -x GRAPH -A      one line per stored path, node walk in field 6
  graph.gfa   vg view -g GRAPH          L-lines give the edges
  mosaic.tsv  vg call --mosaic-out

Checks, in order of how much they tell you:
  1. every segment names a haplotype the GAF carries, and its start/end nodes lie on that walk
  2. the slice between them is taken in the haplotype's own direction of travel
  3. consecutive segments of a path join -- the next begins where the previous ended
  4. every step of the concatenated path is a real edge in the graph
Exit status is non-zero if any check fails; the report says which and how often.

A path is one (contig, strand, fragment) triple: the format says so, and `fragment` only advances
where a boundary was left open, so with --mosaic-patch-gaps on there is one path per strand.

--gfa is optional. Without it checks 1-3 still run, which is enough to say whether the rows tile;
check 4 needs the graph's edges, and the GFA for a whole chromosome is large.

A version 6 mosaic (`vg call --mosaic-exact`) also has `A` rows. Each gives the allele a strand was
called with at one site, as a node walk between the site's two bound nodes. They are applied in
file order to the paths the `H` rows expand to: every stretch of the strand's paths that runs from
the row's start node to its next visit of the end node becomes the row's walk. A row whose start
node the strand does not visit adds its walk as a path of its own. Check 4 then runs on the paths
with the `A` rows applied, so it also says that every `A` row's walk is a walk in the graph. The
one exception is a row that gives a sequence and whose walk is its two bound nodes alone, as the
row of an allele built from the reads is (`vg call --read-consensus`): the sequence joins the two
nodes, so the step between them need not be an edge, and it is counted apart.

An `A` row's sequence is the strand's bases over the row's walk. With --gfa, which gives the
nodes' bases, the sequences are read into the bases of the paths: `StrandWalk.settle` works out
from all of a strand's rows together which bases each visit of a node has, and which lie between
two visits. The number of sequences read, not read and left out is printed.

--vcf CALLS.vcf checks that the paths take the alleles the VCF calls. It reads each record whose ID
names a snarl by its two bound nodes (`>a>b`, without the `_N` suffix of a block record), whose
INFO/AT gives every allele's walk between those nodes, and whose GT is phased `x|y`. Strand 0 has
allele x and strand 1 allele y, and a `.` is skipped. Between the snarl's two bound nodes, the
strand's paths are compared with the called allele's AT walk, and, where that walk is the two
bound nodes alone, the paths' bases there with the allele's. It counts the alleles checked, and
of those the ones whose walk the paths take (`equal`), the ones they do not (`different`), and the
ones whose snarl the strand's paths do not cross. Two more classes count the alleles whose AT walk
is not the walk the genotyper chose for the strand, which the VCF itself shows; `check_vcf` defines
them. A site written as flat records (`vg call --flatten-records`, INFO/FL) has no AT walks, and
is compared by its bases: `check_site_records`. Any `different` allele or strand makes the exit
status non-zero. Every record, allele and strand that is not compared is counted, by reason: block
records, the records of edit sites inside nodes, records with no phase, and strands with no
allele. The node IDs of AT are compared with those of the mosaic as they are, so the VCF must be
written without a node translation.
"""
import sys, argparse, re
from collections import defaultdict


def read_walks(gaf):
    """path name -> list of (node_id, is_reverse), from the GAF path field."""
    walks = {}
    # Haplotypes share most of their steps, so each distinct step is kept once.
    shared = {}
    for line in open(gaf):
        if line.startswith('@') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 6:
            continue
        steps, i, p = [], 0, f[5]
        while i < len(p):
            if p[i] in '><':
                rev = p[i] == '<'
                j = i + 1
                while j < len(p) and p[j] not in '><':
                    j += 1
                step = (int(p[i + 1:j]), rev)
                steps.append(shared.setdefault(step, step))
                i = j
            else:
                i += 1
        walks[f[0]] = steps
    return walks


def read_edges(gfa, sequences):
    """set of ((from, from_rev), (to, to_rev)), both orientations of each L-line.

    Also fills `sequences` with each S-line's node ID and bases.
    """
    e = set()
    for line in open(gfa):
        if line.startswith('S\t'):
            f = line.rstrip('\n').split('\t')
            if f[1].isdigit():
                sequences[int(f[1])] = f[2]
            continue
        if not line.startswith('L\t'):
            continue
        f = line.rstrip('\n').split('\t')
        a, ao, b, bo = int(f[1]), f[2] == '-', int(f[3]), f[4] == '-'
        e.add(((a, ao), (b, bo)))
        e.add(((b, not bo), (a, not ao)))      # the same edge walked the other way
    return e


def pick_copy(w, starts, ends, srev, erev, resume):
    """Pick which visit to the start and end nodes this segment means.

    The file states the orientation at both ends, so a copy whose orientation disagrees is not the
    one meant. Among those that agree, take the pair nearest the point the walk has reached -- the
    segments of a fragment are given in walk order, so that is the copy the previous segment hands
    over to. Returns (None, None) when no pair qualifies, which is a real failure and is counted.
    """
    best = None
    for i in starts:
        if w[i][1] != srev:
            continue
        for j in ends:
            if w[j][1] != erev:
                continue
            # Distance from where the walk stands, then span, so a segment does not leap across
            # the contig to a distant repeat copy when a near one fits.
            cost = (abs(i - resume), abs(j - i))
            if best is None or cost < best[0]:
                best = (cost, i, j)
    if best is None:
        # Orientation may legitimately be the reverse-complement reading of the whole fragment, in
        # which case both ends flip together. Try that before giving up.
        for i in starts:
            if w[i][1] == srev:
                continue
            for j in ends:
                if w[j][1] == erev:
                    continue
                cost = (abs(i - resume), abs(j - i))
                if best is None or cost < best[0]:
                    best = (cost, i, j)
    return (best[1], best[2]) if best else (None, None)


def haplotype_walk(walks, name, wanted, _cache={}):
    """Fragments of one haplotype, each with a node->offset index of the nodes in `wanted`.

    Indexed once per haplotype and cached: a walk here is two million nodes, and rebuilding the
    index per segment turns a linear check into 16 billion operations. Only the nodes that rows
    start or end on are looked up, so only those are indexed.
    """
    if name in _cache:
        return _cache[name]
    out = [w for n, w in walks.items() if n == name or n.startswith(name + '#')
           or n.split('[')[0] == name]
    out.sort(key=len, reverse=True)
    # EVERY occurrence, not just the last. A haplotype can visit a node more than once -- a
    # tandem repeat, an inversion re-entering -- and `{n: i for ...}` silently keeps the final
    # visit, so a segment anchored on an earlier one was sliced from the wrong copy. That is the
    # checker's own ambiguity, not the file's: the file names the copy, by stating the orientation
    # at both ends and by listing segments in walk order.
    occ = []
    for w in out:
        d = defaultdict(list)
        for i, (n, _) in enumerate(w):
            if n in wanted:
                d[n].append(i)
        occ.append((w, d))
    _cache[name] = occ
    return _cache[name]


STEP = re.compile(r'([<>])(\d+)')
WALK = re.compile(r'(?:[<>]\d+)+$')
SNARL_ID = re.compile(r'([<>]\d+[<>]\d+)(_\d+)?$')


def parse_walk(text):
    """'>1<2' -> [(1, False), (2, True)], or None if the text is not a walk."""
    if not WALK.match(text):
        return None
    return [(int(i), o == '<') for o, i in STEP.findall(text)]


def backwards(walk):
    """The same walk read from its end to its start."""
    return [(n, not r) for n, r in reversed(walk)]


class StrandWalk:
    """One strand's paths, held so that a stretch can be replaced without moving the rest.

    Each visit of a node is a slot that points at the next one. A replaced slot is dead.
    """

    def __init__(self):
        self.node, self.next, self.live = [], [], bytearray()
        self.heads = []                    # (fragment, first slot) of each path
        self.at = defaultdict(list)        # node -> slots, for the nodes asked for
        self.joined = set()                # (node, node) steps that an A row's sequence joins
        # What the A rows' sequences say, which the nodes of a walk do not:
        self.bases = {}                    # slot -> the strand's bases over that visit
        self.after = {}                    # slot -> its bases between that visit and the next
        self.fused = {}                    # slot -> the slot before it that holds its bases
        self.unclear = set()               # slots whose bases two sequences give differently
        self.given = []                    # (slots, sequence) of each A row with a sequence

    def add_path(self, fragment, steps):
        first = len(self.node)
        self.node.extend(steps)
        self.next.extend(range(first + 1, first + len(steps)))
        self.next.append(-1)
        self.live.extend(b'\x01' * len(steps))
        self.heads.append((fragment, first))
        return first

    def index(self, wanted):
        """Index the live slots of the nodes in `wanted`. Slots added later are indexed too."""
        self.at = defaultdict(list)
        for slot, n in enumerate(self.node):
            if self.live[slot] and n in wanted:
                self.at[n].append(slot)

    def crossings(self, start, end):
        """Each stretch from `start` to the next visit of `end`, as (first slot, last slot)."""
        found = []
        for slot in self.at.get(start, ()):
            if not self.live[slot]:
                continue
            at = slot
            while at >= 0 and self.node[at] != end:
                at = self.next[at]
                if at >= 0 and self.node[at][0] == start[0]:
                    at = -1                # back at the start node before reaching the end
            if at >= 0:
                found.append((slot, at))
        return found

    def slots(self, crossing):
        out, at = [crossing[0]], crossing[0]
        while at != crossing[1]:
            at = self.next[at]
            out.append(at)
        return out

    def taken(self, crossing):
        return [self.node[slot] for slot in self.slots(crossing)]

    def give(self, slots, text):
        """Note that the visits in `slots`, the walk of an A row, have the row's sequence `text`
        between them. `settle` then works out which bases each visit has."""
        self.given.append((tuple(slots), text))

    def settle(self, sequences, counts):
        """Work out the bases of each visit from the sequences that `give` noted.

        A sequence is the strand's bases over a whole walk, and does not say which are over
        which node. Rows share nodes, and each has every change the strand makes to its nodes,
        so the rows are settled together: a row is settled once the bases of all its visits but
        one are known, which makes the last one known, and a row of one visit says its bases
        outright. Where no row can be settled that way, the first one left is settled by taking
        each visit not yet known, from the two ends of the walk, to have its node's own bases if
        the sequence has them there.

        Bases left over between two visits that keep theirs lie between the two nodes, as the
        bases of an allele built from the reads lie between its two bound nodes. Visits that
        cannot be told apart share their bases: the first holds them, and the others are fused
        to it. A row that disagrees with what is known makes its visits unclear.
        """
        known = {}
        pending = [row for row in dict.fromkeys(self.given)
                   if all(self.live[slot] for slot in row[0])]
        counts['A row sequences left out, their walk being replaced by a later row'] += \
            len(set(self.given)) - len(pending)
        own = {}
        for slots, text in pending:
            for slot in slots:
                if slot not in own:
                    own[slot] = oriented(self.node[slot], sequences)
        while pending:
            settled = [row for row in pending if self.place(row, known, own, False)]
            if settled:
                pending = [row for row in pending if row not in settled]
            else:
                self.place(pending.pop(0), known, own, True)
        for slot, bases in known.items():
            if bases != own[slot]:
                self.bases[slot] = bases
        self.given = []

    def place(self, row, known, own, guess):
        """Settle one row, if it can be: see `settle`. With `guess`, it always can."""
        slots, text = row
        lo, hi, a, b = 0, len(slots) - 1, 0, len(text)

        def bases(i):
            if slots[i] in known:
                return known[slots[i]]
            return own[slots[i]] if guess and slots[i] not in self.fused else None

        while hi > lo and bases(hi) is not None and text.endswith(bases(hi), a, b):
            b -= len(bases(hi))
            hi -= 1
        while lo <= hi and bases(lo) is not None and text.startswith(bases(lo), a, b):
            a += len(bases(lo))
            lo += 1
        rest = text[a:b]
        if lo > hi:
            # Every visit has bases. What is left lies between two of them, or after the last.
            for slot in slots:
                known.setdefault(slot, own[slot])
            if rest and lo < len(slots):
                if self.after.setdefault(slots[lo - 1], rest) != rest:
                    self.unclear.update(slots)
            elif rest:
                self.unclear.update(slots)
            return True
        if any(slots[i] in known or slots[i] in self.fused for i in range(lo, hi + 1)):
            # A visit between lo and hi has bases that this row does not have there.
            if lo == hi or guess:
                self.unclear.update(slots)
                return True
            return False
        if lo < hi and not guess:
            return False
        for i in list(range(lo)) + list(range(hi + 1, len(slots))):
            known.setdefault(slots[i], own[slots[i]])
        if lo == hi:
            # One visit is left. Bases beside its node's own lie between it and the visit
            # next to it in the walk, if the row has one there.
            mine = own[slots[lo]]
            if rest != mine and rest.endswith(mine) and lo > 0:
                self.after[slots[lo - 1]] = rest[:len(rest) - len(mine)]
                rest = mine
            elif rest != mine and rest.startswith(mine) and lo + 1 < len(slots):
                self.after[slots[lo]] = rest[len(mine):]
                rest = mine
            known[slots[lo]] = rest
            return True
        known[slots[lo]] = rest
        for i in range(lo + 1, hi + 1):
            known[slots[i]] = ''
            self.fused[slots[i]] = slots[lo]
        return True

    def interior(self, crossing, sequences):
        """The strand's bases between the two bound visits of a stretch, with what the A rows'
        sequences say, or None if a node's bases are not known or a sequence does not say which
        of its bases lie between the bounds."""
        slots = self.slots(crossing)
        if len(slots) < 2 or slots[-1] in self.fused or slots[1] in self.fused \
                or not self.unclear.isdisjoint(slots):
            return None
        out = [self.after.get(slots[0], '')]
        try:
            for slot in slots[1:-1]:
                out.append(self.bases[slot] if slot in self.bases
                           else oriented(self.node[slot], sequences))
                out.append(self.after.get(slot, ''))
        except KeyError:
            return None
        return ''.join(out)

    def replace(self, crossing, walk):
        if crossing[0] == crossing[1]:
            return
        at = self.next[crossing[0]]
        while at != crossing[1]:
            self.live[at] = 0
            at = self.next[at]
        before = crossing[0]
        for n in walk[1:-1]:
            slot = len(self.node)
            self.node.append(n)
            self.next.append(-1)
            self.live.append(1)
            self.at[n].append(slot)
            self.next[before] = slot
            before = slot
        self.next[before] = crossing[1]

    def paths(self):
        for fragment, first in self.heads:
            path, at = [], first
            while at >= 0:
                path.append(self.node[at])
                at = self.next[at]
            yield fragment, path


def apply_allele_rows(strand_walk, rows, fail, counts, sequences):
    """Apply one strand's `A` rows, in file order: the walk of each, and its sequence."""
    wanted = set()
    for f in rows:
        wanted.add((int(f[4]) // 2, bool(int(f[4]) & 1)))
    strand_walk.index(wanted)
    for f in rows:
        # A: contig strand fragment start_node end_node walk sequence
        if len(f) != 8:
            fail['A row does not have 8 columns'] += 1
            continue
        walk = parse_walk(f[6])
        start = (int(f[4]) // 2, bool(int(f[4]) & 1))
        end = (int(f[5]) // 2, bool(int(f[5]) & 1))
        if not walk or walk[0] != start or walk[-1] != end:
            fail['A row walk does not run from its start_node to its end_node'] += 1
            continue
        counts['A rows'] += 1
        text = None if f[7] == '.' else f[7].upper()
        counts['A rows with a sequence'] += text is not None
        if text is not None and not re.fullmatch('[ACGTN]*', text):
            fail['A row sequence is not bases'] += 1
            continue
        if text is not None and len(walk) == 2:
            # The sequence gives the bases between the two bound nodes.
            strand_walk.joined.add((walk[0], walk[1]))
        found = strand_walk.crossings(start, end) if len(walk) == 1 or start != end else []
        if not found:
            # The strand has no path here, so the walk is a path of its own.
            first = strand_walk.add_path(f[3], walk)
            for i, n in enumerate(walk):
                strand_walk.at[n].append(first + i)
            counts['A rows that join no H row'] += 1
            found = [(first, first + len(walk) - 1)]
        else:
            for crossing in found:
                strand_walk.replace(crossing, walk)
        if text is None:
            continue
        if not sequences:
            counts['A row sequences not read, for want of --gfa'] += 1
        elif all(n in sequences for n, _ in walk):
            for crossing in found:
                strand_walk.give(strand_walk.slots(crossing), text)
            counts['A row sequences read'] += 1
        else:
            fail['A row sequence is over a node the graph does not have'] += 1
    if sequences:
        strand_walk.settle(sequences, counts)


def outside_nested(walk, nested, own):
    """`walk` with each stretch that crosses a snarl of `nested` cut down to its two bound nodes.

    `nested` maps a start bound to the end bounds of the snarls that start there, and `own` is the
    (start, end) of the snarl `walk` itself crosses, which is left alone. A snarl inside a cut
    stretch goes with it.
    """
    out, i = [], 0
    while i < len(walk):
        out.append(walk[i])
        ends = nested.get(walk[i], ())
        j = next((k for k in range(i + 1, len(walk))
                  if walk[k] in ends and (walk[i], walk[k]) != own), -1) if ends else -1
        if j < 0:
            i += 1
        elif j + 1 < len(walk) and walk[j] in nested:
            i = j                          # the end bound may start the next snarl
        else:
            out.append(walk[j])
            i = j + 1
    return out


def divergence(path, walk, limit=12):
    """Where two walks between the same bound nodes first part and last rejoin, as text."""
    a = next((i for i, (x, y) in enumerate(zip(path, walk)) if x != y), min(len(path), len(walk)))
    b = next((i for i, (x, y) in enumerate(zip(reversed(path), reversed(walk))) if x != y),
             min(len(path), len(walk)))
    text = lambda w: ''.join(('<' if r else '>') + str(n) for n, r in w[:limit]) + \
        ('...' if len(w) > limit else '')
    return (f"path takes {text(path[max(a - 1, 0):len(path) - b + 1])} where the allele has "
            f"{text(walk[max(a - 1, 0):len(walk) - b + 1])}")


COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def oriented(step, sequences):
    """The bases of one visit of a node, in upper case, read the way the visit reads it."""
    node, reverse = step
    bases = sequences[node].upper()
    return bases[::-1].translate(COMPLEMENT) if reverse else bases


def spell(walk, sequences):
    """The bases of a walk, or None if a node's bases are not known."""
    try:
        return ''.join(sequences[n][::-1].translate(COMPLEMENT) if r else sequences[n]
                       for n, r in walk)
    except KeyError:
        return None


def revcomp(bases):
    return bases[::-1].translate(COMPLEMENT)


class Reference:
    """One contig's reference path: where each node it visits once starts on it."""

    def __init__(self, walk, sequences):
        self.walk, self.index, self.offset, self.length = walk, {}, {}, {}
        seen, at = set(), 0
        for i, (node, reverse) in enumerate(walk):
            if node in seen:
                self.index.pop(node, None)
                self.offset.pop(node, None)
            else:
                seen.add(node)
                self.index[node], self.offset[node] = i, at
                self.length[node] = len(sequences.get(node, ''))
            at += len(sequences.get(node, ''))

    def between(self, bounds, sequences):
        """(offset of the first base after the start bound's node, the bases from there to the
        end bound's node), or None unless the path visits both nodes once, forwards, in order."""
        (a, a_reverse), (b, b_reverse) = bounds
        i, j = self.index.get(a), self.index.get(b)
        if i is None or j is None or j <= i or self.walk[i] != (a, a_reverse) \
                or self.walk[j] != (b, b_reverse):
            return None
        bases = spell(self.walk[i + 1:j], sequences)
        return None if bases is None else (self.offset[a] + len(sequences[a]), bases.upper())


def expected_interior(reference_walk, ref, allele, pos, start_offset, sequences):
    """The bases a record's allele has between its site's bound nodes, or None if they cannot
    be told. A list of them, one for each place REF can lie, where the record's POS cannot say.

    A record's alleles are what is left of its walks' bases once the bases all of them start and
    end with are cut off, so the reference walk spells some bases, REF, and some more bases.
    `start_offset` is where the reference walk's first node starts on the reference path, or None.
    """
    reference = spell(reference_walk, sequences)
    if reference is None or not ref or not re.fullmatch('[ACGTNacgtn]+', allele):
        return None
    reference, ref, allele = reference.upper(), ref.upper(), allele.upper()
    first = len(sequences[reference_walk[0][0]])
    last = len(sequences[reference_walk[-1][0]])
    places = []
    if start_offset is not None and reference.startswith(ref, pos - 1 - start_offset) \
            and pos - 1 >= start_offset:
        places.append(pos - 1 - start_offset)
    else:
        at = reference.find(ref)
        while at >= 0:
            places.append(at)
            at = reference.find(ref, at + 1)
    out = []
    for at in places:
        whole = reference[:at] + allele + reference[at + len(ref):]
        # The bound nodes' own bases are not the site's, and stay as they are.
        if len(whole) >= first + last and whole.startswith(reference[:first]) \
                and whole.endswith(reference[len(reference) - last:]):
            out.append(whole[first:len(whole) - last])
    return out or None


def beside_bounds(reference, bounds, inserted):
    """The bases that the records of edit sites put right after the node of a site's start bound,
    and right before the node of its end bound, on one strand. They lie between the two nodes,
    and are no part of the site.

    `inserted` maps a place on the reference path to the bases an edit site's record puts there
    on the strand.
    """
    (a, a_reverse), (b, b_reverse) = bounds
    if not inserted or reference is None or a_reverse or b_reverse \
            or a not in reference.offset or b not in reference.offset:
        return '', ''
    first = reference.offset[a] + reference.length[a]
    before, after = inserted.get(first, ''), inserted.get(reference.offset[b], '')
    return before, ('' if first == reference.offset[b] else after)


def without_beside(interior, beside, counts):
    """`interior`, the bases a strand has between two bound nodes, without the bases `beside`
    that edit sites' records put right beside the two nodes."""
    before, after = beside
    if before and interior.startswith(before):
        interior = interior[len(before):]
        counts["with an edit site's bases beside a bound node"] += 1
    if after and interior.endswith(after):
        interior = interior[:len(interior) - len(after)]
        counts["with an edit site's bases beside a bound node"] += 1
    return interior


def check_site_records(name, records, strand, strand_walk, reference, sequences, counts,
                       inserted):
    """Compare one strand with the flat records one site is written as.

    Each record replaces its REF, at its POS on the reference, with the allele the strand is
    called with. The reference from the first base of the node of the site's start bound to the
    last base of the node of its end bound, with every record's allele put in, is what the
    strand must spell there: a record can lie on the bound nodes, since a gap is written as far
    left as it can go. The strand's bases on the two bound nodes are taken to be the nodes' own.
    Returns a description of a difference, or None.
    """
    bounds = parse_walk(name)
    span = reference.between(bounds, sequences) if reference and sequences else None
    if span is None:
        counts['strands not compared: the reference path does not give the site'] += 1
        return None
    first, last = (oriented(bound, sequences) for bound in bounds)
    start, bases = span[0] - len(first), first + span[1] + last
    changes = []
    for pos, ref, alts, gt in records:
        if '|' not in gt:
            called = gt.split('/')
            if len(set(called)) > 1:
                counts['strands not compared: a record of the site is not phased'] += 1
                return None
            called = called * 2
        else:
            called = gt.split('|')
        g = called[strand] if strand < len(called) else called[0]
        if g == '.':
            counts['strands not compared: a record of the site has no allele for the strand'] += 1
            return None
        if g == '0':
            continue
        allele = alts[int(g) - 1].upper()
        if not re.fullmatch('[ACGTN]+', allele):
            counts['strands not compared: a record of the site has a symbolic allele'] += 1
            return None
        ref = ref.upper()
        k = 0
        while k < len(ref) and k < len(allele) and ref[k] == allele[k]:
            k += 1
        changes.append((pos - 1 + k - start, len(ref) - k, allele[k:]))
    changes.sort()
    out, at = [], 0
    for where, removed, added in changes:
        if where < at or where + removed > len(bases):
            return f"{name}: a record's allele lies outside the site, or on another record's"
        out.append(bases[at:where])
        out.append(added)
        at = where + removed
    out.append(bases[at:])
    expected = ''.join(out)
    beside = beside_bounds(reference, bounds, inserted)
    forwards = [(c, False) for c in strand_walk.crossings(bounds[0], bounds[1])]
    flipped = backwards(bounds)
    crossings = forwards + [(c, True) for c in strand_walk.crossings(flipped[0], flipped[1])]
    if not crossings:
        counts['strand does not cross'] += 1
        return None
    counts['checked'] += 1
    for crossing, turned in crossings:
        spelled = strand_walk.interior(crossing, sequences)
        if spelled is None:
            counts['strands not compared: a sequence does not say which bases lie in the site'] += 1
            counts['checked'] -= 1
            return None
        spelled = without_beside(revcomp(spelled) if turned else spelled, beside, counts)
        if first + spelled + last != expected:
            return (f"{name}: the strand has {len(spelled)} bases between the bounds, and the "
                    f"records give {len(expected) - len(first) - len(last)}")
    counts['equal'] += 1
    return None


def check_vcf(vcf, strands, sequences, references, quiet):
    """Compare the paths with the phased snarl records of a VCF. Returns two sets of counts.

    Each called allele of a record whose ID names a snarl is put in the first of these classes
    that holds at every crossing of its snarl:
      equal                         the path takes the allele's AT walk. Where that walk is the
                                    two bound nodes alone, as the walk of an allele built from
                                    the reads is, the path must also have the allele's bases
                                    between them, when the nodes' bases are known (--gfa).
      equal outside nested records  it does so except inside snarls that other records call an
                                    allele of on this strand; those records are checked
                                    themselves
      same bases                    it takes another walk, and has the allele's bases between
                                    the snarl's bound nodes. This needs the nodes' bases, so
                                    only with --gfa.
      different                     none of these
    AT gives one walk for each allele, and the genotyper's walk for a strand can be another one:
    allele 0 also stands for walks that differ from the reference walk only inside nested sites,
    an allele's bases are spelled with the alleles chosen at the sites nested in it on the first
    strand that carries it, and walks that spell the same bases are one allele. So a path that
    takes a walk of three nodes or more is not held to the record's bases. The bases a path has
    are those of its nodes, with what the A rows' sequences give.

    A site written as flat records (INFO/FL, `_N` after the site's ID) is compared by its bases
    alone (`check_site_records`), since its records give no walk and spell each strand whole.
    The second set of counts is of the strands of those sites. A block record, which has `_N`
    and no FL, gives part of its site's walk, and is counted and not compared.
    """
    records = []
    site_records = defaultdict(list)
    wanted = defaultdict(set)
    counts = defaultdict(int)
    site_counts = defaultdict(int)
    # For each strand, the snarls a record calls an allele of there, read both ways.
    nested = defaultdict(lambda: defaultdict(set))
    # For each strand, the bases the records of edit sites insert, by place on the reference.
    inserted = defaultdict(dict)
    for line in open(vcf):
        if line.startswith('#'):
            continue
        t = line.rstrip('\n').split('\t')
        snarl = SNARL_ID.match(t[2]) if len(t) >= 10 else None
        gt = t[9].split(':')[0] if len(t) >= 10 else ''
        if not snarl:
            # The record of an edit site, whose ID names a node and an offset. Its changes are
            # inside nodes, where no walk tells one allele from another.
            counts['records skipped: the ID names no snarl'] += 1
            alleles = [t[3].upper()] + t[4].upper().split(',')
            for strand, g in enumerate(gt.split('|') if '|' in gt else ()):
                if g not in ('.', '0') and alleles[int(g)].startswith(alleles[0]):
                    inserted[(t[0], str(strand))][int(t[1]) - 1 + len(alleles[0])] = \
                        alleles[int(g)][len(alleles[0]):]
            continue
        at = [x[3:] for x in t[7].split(';') if x.startswith('AT=')]
        bounds = parse_walk(snarl.group(1))
        flipped = backwards(bounds)
        wanted[t[0]].update((bounds[0], flipped[0]))
        if snarl.group(2):
            # One of several records of its site. The site still counts as one that has
            # records, for the sites around it.
            for strand, g in enumerate(gt.replace('/', '|').split('|')):
                if g != '.':
                    nested[(t[0], str(strand))][bounds[0]].add(bounds[1])
                    nested[(t[0], str(strand))][flipped[0]].add(flipped[1])
            if 'FL' in t[7].split(';'):
                site_records[(t[0], snarl.group(1))].append(
                    (int(t[1]), t[3], t[4].split(','), gt))
            else:
                counts['block records not checked'] += 1
            continue
        if '|' not in gt:
            counts['records skipped: the genotype is not phased'] += 1
            continue
        if not at:
            counts['records skipped: no AT'] += 1
            continue
        for strand, g in enumerate(gt.split('|')):
            if g != '.':
                nested[(t[0], str(strand))][bounds[0]].add(bounds[1])
                nested[(t[0], str(strand))][flipped[0]].add(flipped[1])
        alleles = [parse_walk(x) for x in at[0].split(',')]
        if any(w is None or w[0] != bounds[0] or w[-1] != bounds[1] for w in alleles):
            counts['records skipped: an allele has no full walk in AT'] += 1
            continue
        records.append((t[0], int(t[1]), t[2], gt, alleles, [t[3]] + t[4].split(',')))
    for (contig, strand), strand_walk in strands.items():
        strand_walk.index(wanted.get(contig, ()))
    shown = 0
    classes = ['equal', 'equal outside nested records', 'same bases', 'different']
    for contig, pos, name, gt, alleles, strings in records:
        reference = references.get(contig)
        for strand, g in enumerate(gt.split('|')):
            if g == '.':
                counts['alleles skipped: the strand has none'] += 1
                continue
            if (contig, str(strand)) not in strands:
                counts['alleles skipped: the mosaic has no such strand'] += 1
                continue
            strand_walk = strands[(contig, str(strand))]
            called = alleles[int(g)]
            counts['checked'] += 1
            # A strand can cross a snarl from its end to its start, inside an inverted allele of
            # an enclosing snarl, and then it takes the called walk read backwards.
            taken = [(c, called) for c in strand_walk.crossings(called[0], called[-1])]
            flipped = backwards(called)
            taken += [(c, flipped) for c in strand_walk.crossings(flipped[0], flipped[-1])]
            inside = nested[(contig, str(strand))]
            if not taken:
                counts['strand does not cross'] += 1
                continue
            expected = None
            if sequences and int(g) < len(strings):
                expected = expected_interior(
                    alleles[0], strings[0], strings[int(g)], pos,
                    reference.offset.get(alleles[0][0][0]) if reference else None, sequences)
                beside = beside_bounds(reference, (alleles[0][0], alleles[0][-1]),
                                       inserted.get((contig, str(strand))))
            worst, bare = 0, False
            for crossing, walk in taken:
                path = strand_walk.taken(crossing)
                own = (walk[0], walk[-1])
                spelled = strand_walk.interior(crossing, sequences) if expected else None
                if spelled is not None:
                    spelled = without_beside(revcomp(spelled) if walk is not called else spelled,
                                             beside, counts)
                same = None if spelled is None else spelled in expected
                if path == walk and len(walk) == 2:
                    # The walk does not say what lies between the two nodes.
                    found = 3 if same is False else 0
                    bare = bare or same is None
                elif path == walk:
                    found = 0
                elif outside_nested(path, inside, own) == outside_nested(walk, inside, own):
                    found = 1
                elif same:
                    found = 2
                else:
                    found = 3
                if found >= worst:
                    worst = found
                    example = (path, walk, same)
            counts[classes[worst]] += 1
            if len(called) == 2 and worst == 0:
                counts['of the equal, walks of two nodes whose bases were '
                       + ('not compared' if bare else 'compared')] += 1
            if worst == 3 and shown < 10 and not quiet:
                shown += 1
                why = ('takes the walk, with other bases than the allele: '
                       if example[2] is False and example[0] == example[1] else '')
                print(f"  DIFFERENT  {contig}:{pos} {name} GT {gt} strand {strand}: {why}"
                      f"{divergence(example[0], example[1])}", file=sys.stderr)
    for (contig, name), group in site_records.items():
        site_counts['sites'] += 1
        for strand in (0, 1):
            if (contig, str(strand)) not in strands:
                continue
            problem = check_site_records(name, group, strand, strands[(contig, str(strand))],
                                         references.get(contig), sequences, site_counts,
                                         inserted.get((contig, str(strand))))
            if problem:
                site_counts['different'] += 1
                if shown < 10 and not quiet:
                    shown += 1
                    print(f"  DIFFERENT  {contig} strand {strand} {problem}", file=sys.stderr)
    return counts, site_counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mosaic', required=True)
    ap.add_argument('--gaf', required=True)
    ap.add_argument('--gfa', help='optional; without it check 4 is skipped')
    ap.add_argument('--vcf', help='optional; the calls the mosaic was written with, to check '
                                  'that the paths take the called alleles')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args()

    walks = read_walks(a.gaf)
    sequences = {}
    edges = read_edges(a.gfa, sequences) if a.gfa else None
    hapname, threads, reference_names = {}, defaultdict(list), []
    allele_rows = defaultdict(list)
    for line in open(a.mosaic):
        f = line.rstrip('\n').split('\t')
        if f[0] == '#haplotype':
            hapname[f[1]] = f[2]
        elif f[0] == '#reference':
            reference_names += f[1:]
        elif f[0] == '#mosaic-version' and f[1] not in ('5', '6'):
            print(f"  FAIL  mosaic version {f[1]} is not one this script reads (5 or 6)")
            return 1
        elif f[0] == 'H':
            threads[(f[1], f[2], f[3])].append(f)      # contig, strand, fragment
        elif f[0] == 'A':
            allele_rows[(f[1], f[2])].append(f)        # contig, strand, in file order

    fail = defaultdict(int)
    strands = defaultdict(StrandWalk)
    row_ends = {int(f[c]) // 2 for segs in threads.values() for f in segs for c in (6, 7)}
    reversed_slices = [0]
    orient_mismatch = [0]
    total_steps = 0
    for (contig, strand, frag), segs in sorted(threads.items()):
        path = []
        pick_copy.at = 0          # each fragment starts its own walk
        for f in segs:
            # H: contig strand fragment ref_start ref_end start_node end_node
            #     hap_index haplotype sites gbwt_offset
            # start_node and end_node are ORIENTED: id * 2 + is_reverse.
            hap = f[9]
            start, srev = int(f[6]) // 2, bool(int(f[6]) & 1)
            end,   erev = int(f[7]) // 2, bool(int(f[7]) & 1)
            if f[8] == '*':
                fail['segment has no haplotype to walk'] += 1
                continue
            frags = (haplotype_walk(walks, hap, row_ends)
                     or haplotype_walk(walks, hap + '#' + contig, row_ends))
            if not frags:
                fail['haplotype not found in the graph'] += 1
                continue
            piece = None
            # Where the previous segment left off, so the copy chosen here is the one the walk
            # actually reaches rather than whichever happens to be last in the fragment.
            resume = getattr(pick_copy, 'at', 0)
            for w, idx in frags:
                if start in idx and end in idx:
                    i, j = pick_copy(w, idx[start], idx[end], srev, erev, resume)
                    if i is None:
                        continue
                    pick_copy.at = j
                    if i <= j:
                        piece = w[i:j + 1]
                    else:
                        # The segment runs against the haplotype's stored direction. That is not an
                        # error and not rare: on this graph whole haplotypes are stored reversed --
                        # recombination#17 is '<' at every one of its 2.0 M steps -- so a segment
                        # given in reference order is a walk along the reverse complement. Take the
                        # slice and flip both the order and every step's orientation.
                        piece = [(n, not r) for n, r in reversed(w[j:i + 1])]
                        reversed_slices[0] += 1
                    # The file states the orientation at both ends; the haplotype's own walk states
                    # it too. They have to agree, or the row's oriented anchors are not the walk's.
                    if piece and (piece[0][1] != srev or piece[-1][1] != erev):
                        orient_mismatch[0] += 1
                        print(f"  ORIENT  {f}  walk_start={piece[0]} walk_end={piece[-1]} "
                              f"stated=({start},{srev})->({end},{erev})", file=sys.stderr)
                    break
            if piece is None:
                fail['segment endpoints not both on one fragment of its haplotype'] += 1
                continue
            if path and path[-1] != piece[0]:
                if edges is not None and (path[-1], piece[0]) in edges:
                    fail['segments abut by an edge rather than sharing a node'] += 1
                else:
                    fail['segments do not join'] += 1
                    print(f"  NOJOIN  {f}  prev_end={path[-1]} this_start={piece[0]}",
                          file=sys.stderr)
            path.extend(piece if not path else piece[1:])
        if path:
            strands[(contig, strand)].add_path(frag, path)
        if not a.quiet:
            print(f"  {contig} strand {strand} fragment {frag}: "
                  f"{len(segs)} segments -> {len(path)} node steps")

    # The A rows change the paths, so they are applied before the paths are checked against the
    # graph's edges.
    allele_counts = defaultdict(int)
    for key, rows in allele_rows.items():
        apply_allele_rows(strands[key], rows, fail, allele_counts, sequences)
    for (contig, strand), strand_walk in sorted(strands.items()):
        for frag, path in strand_walk.paths():
            for i in range(len(path) - 1):
                total_steps += 1
                if (path[i], path[i + 1]) in strand_walk.joined:
                    allele_counts['steps joined by a sequence'] += 1
                elif edges is not None and (path[i], path[i + 1]) not in edges:
                    fail['step is not an edge in the graph'] += 1
                    print(f"  NOEDGE  {contig} strand {strand} frag {frag} step {i}: "
                          f"{path[i]} -> {path[i+1]}", file=sys.stderr)

    against = (f"against {len(edges)//2} edges" if edges is not None
               else "(no --gfa given: edges unchecked)")
    print(f"\n  {len(threads)} paths, {total_steps} steps {against}")
    print(f"  {reversed_slices[0]} segments walk against their haplotype's stored direction")
    if orient_mismatch[0]:
        print(f"  {orient_mismatch[0]} segments whose stated orientation disagrees with the walk")
    if allele_counts:
        print(f"  {allele_counts['A rows']} A rows applied, of which "
              f"{allele_counts['A rows with a sequence']} give a sequence and "
              f"{allele_counts['A rows that join no H row']} join no H row")
    if allele_counts['steps joined by a sequence']:
        print(f"  {allele_counts['steps joined by a sequence']} steps joined by an A row's "
              f"sequence rather than an edge")
    if allele_counts['A rows with a sequence']:
        left_out = 'A row sequences left out, their walk being replaced by a later row'
        print(f"  {allele_counts['A row sequences read']} A row sequences read into the paths' "
              f"bases, {allele_counts['A row sequences not read, for want of --gfa']} not read "
              f"for want of --gfa, {allele_counts[left_out]} left out, their walk being replaced "
              f"by a later row")
    if a.vcf:
        # The reference path of each contig, as the VCF names the contig.
        references = {}
        for name in reference_names if sequences else ():
            if name in walks:
                references[name.split('#')[-1]] = Reference(walks[name], sequences)
        vcf_counts, site_counts = check_vcf(a.vcf, strands, sequences, references, a.quiet)
        print(f"  called alleles against {a.vcf}: {vcf_counts['checked']} checked, "
              f"{vcf_counts['equal']} equal, "
              f"{vcf_counts['equal outside nested records']} equal outside nested records, "
              f"{vcf_counts['same bases']} same bases by another walk, "
              f"{vcf_counts['different']} different, "
              f"{vcf_counts['strand does not cross']} strand does not cross")
        for k, v in sorted(vcf_counts.items()):
            if k.startswith(('records skipped', 'alleles skipped', 'of the equal',
                             'block records', 'with an edit')):
                print(f"  {v} {k}")
        if vcf_counts['different']:
            fail['called allele is not what the path takes'] = vcf_counts['different']
        if site_counts['sites']:
            print(f"  sites written as flat records: {site_counts['sites']} sites, "
                  f"{site_counts['checked']} strands checked, {site_counts['equal']} spell "
                  f"their records, {site_counts['different']} different, "
                  f"{site_counts['strand does not cross']} strand does not cross")
            for k, v in sorted(site_counts.items()):
                if k.startswith(('strands not compared', 'with an edit')):
                    print(f"  {v} {k}")
            if site_counts['different']:
                fail["strand does not spell its site's records"] = site_counts['different']
    if not fail:
        print("  OK: every path expands to an exact walk in the graph")
        return 0
    for k, v in sorted(fail.items(), key=lambda kv: -kv[1]):
        print(f"  FAIL {v:>7}  {k}")
    return 1


if __name__ == '__main__':
    sys.exit(main())
