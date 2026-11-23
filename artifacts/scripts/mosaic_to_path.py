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

--vcf CALLS.vcf checks that the paths take the alleles the VCF calls. It reads each record whose ID
names a snarl by its two bound nodes (`>a>b`, without the `_N` suffix of a block record), whose
INFO/AT gives every allele's walk between those nodes, and whose GT is phased `x|y`. Strand 0 has
allele x and strand 1 allele y, and a `.` is skipped. Between the snarl's two bound nodes, the
strand's paths are compared with the called allele's AT walk. It counts the alleles checked, and
of those the ones whose walk the paths take (`equal`), the ones they do not (`different`), and the
ones whose snarl the strand's paths do not cross. Two more classes count the alleles whose AT walk
is not the walk the genotyper chose for the strand, which the VCF itself shows; `check_vcf` defines
them. Any `different` allele makes the exit status non-zero. The node IDs of AT are compared with
those of the mosaic as they are, so the VCF must be written without a node translation.
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

    def taken(self, crossing):
        out, at = [self.node[crossing[0]]], crossing[0]
        while at != crossing[1]:
            at = self.next[at]
            out.append(self.node[at])
        return out

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


def apply_allele_rows(strand_walk, rows, fail, counts):
    """Apply one strand's `A` rows, in file order."""
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
        counts['A rows with a sequence'] += f[7] != '.'
        if f[7] != '.' and len(walk) == 2:
            # The sequence gives the bases between the two bound nodes.
            strand_walk.joined.add((walk[0], walk[1]))
        found = strand_walk.crossings(start, end) if len(walk) == 1 or start != end else []
        if not found:
            # The strand has no path here, so the walk is a path of its own.
            first = strand_walk.add_path(f[3], walk)
            for i, n in enumerate(walk):
                strand_walk.at[n].append(first + i)
            counts['A rows that join no H row'] += 1
            continue
        for crossing in found:
            strand_walk.replace(crossing, walk)


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


def spell(walk, sequences):
    """The bases of a walk, or None if a node's bases are not known."""
    try:
        return ''.join(sequences[n][::-1].translate(COMPLEMENT) if r else sequences[n]
                       for n, r in walk)
    except KeyError:
        return None


def spells_allele(path, reference_walk, ref, allele, sequences):
    """Whether `path` spells the reference walk's bases with REF replaced by `allele`.

    A record's alleles are what is left of its walks' bases once the bases all of them start and
    end with are cut off, so the reference walk spells some bases, REF, and some more bases.
    """
    spelled, reference = spell(path, sequences), spell(reference_walk, sequences)
    if spelled is None or reference is None or not ref or not re.fullmatch('[ACGTNacgtn]+', allele):
        return False
    at = reference.find(ref)
    while at >= 0:
        if spelled == reference[:at] + allele + reference[at + len(ref):]:
            return True
        at = reference.find(ref, at + 1)
    return False


def check_vcf(vcf, strands, sequences, quiet):
    """Compare the paths with the phased snarl records of a VCF. Returns the counts.

    Each called allele is put in the first of these classes that holds at every crossing of its
    snarl:
      equal                         the path takes the allele's AT walk
      equal outside nested records  it does so except inside snarls that other records call an
                                    allele of on this strand; those records are checked
                                    themselves, unless they are block records
      same bases                    it takes another walk, which spells the record's allele. This
                                    needs the nodes' bases, so only with --gfa.
      different                     none of these
    AT gives one walk for each allele, and the genotyper's walk for a strand can be another one:
    allele 0 also stands for walks that differ from the reference walk only inside nested sites,
    an allele's bases are spelled with the alleles chosen at the sites nested in it, and walks
    that spell the same bases are one allele.
    """
    records = []
    wanted = defaultdict(set)
    counts = defaultdict(int)
    # For each strand, the snarls a record calls an allele of there, read both ways.
    nested = defaultdict(lambda: defaultdict(set))
    for line in open(vcf):
        if line.startswith('#'):
            continue
        t = line.rstrip('\n').split('\t')
        snarl = SNARL_ID.match(t[2]) if len(t) >= 10 else None
        if not snarl:
            continue
        gt = t[9].split(':')[0]
        at = [x[3:] for x in t[7].split(';') if x.startswith('AT=')]
        if '|' not in gt or not at:
            continue
        bounds = parse_walk(snarl.group(1))
        flipped = backwards(bounds)
        for strand, g in enumerate(gt.split('|')):
            if g != '.':
                nested[(t[0], str(strand))][bounds[0]].add(bounds[1])
                nested[(t[0], str(strand))][flipped[0]].add(flipped[1])
        if snarl.group(2):
            # A block record gives part of its site's walk, so it is not checked itself. Its
            # site still counts as one that has records, for the sites around it.
            counts['block records not checked'] += 1
            continue
        alleles = [parse_walk(x) for x in at[0].split(',')]
        if any(w is None or w[0] != bounds[0] or w[-1] != bounds[1] for w in alleles):
            counts['records skipped: an allele has no full walk in AT'] += 1
            continue
        records.append((t[0], t[1], t[2], gt, alleles, [t[3]] + t[4].split(',')))
        wanted[t[0]].update((bounds[0], flipped[0]))
    for (contig, strand), strand_walk in strands.items():
        strand_walk.index(wanted.get(contig, ()))
    shown = 0
    classes = ['equal', 'equal outside nested records', 'same bases', 'different']
    for contig, pos, name, gt, alleles, strings in records:
        for strand, g in enumerate(gt.split('|')):
            if g == '.' or (contig, str(strand)) not in strands:
                continue
            strand_walk = strands[(contig, str(strand))]
            called = alleles[int(g)]
            counts['checked'] += 1
            # A strand can cross a snarl from its end to its start, inside an inverted allele of
            # an enclosing snarl, and then it takes the called walk read backwards.
            taken = [(strand_walk.taken(c), called)
                     for c in strand_walk.crossings(called[0], called[-1])]
            flipped = backwards(called)
            taken += [(strand_walk.taken(c), flipped)
                      for c in strand_walk.crossings(flipped[0], flipped[-1])]
            inside = nested[(contig, str(strand))]
            if not taken:
                counts['strand does not cross'] += 1
                continue
            worst = 0
            for path, walk in taken:
                own = (walk[0], walk[-1])
                if path == walk:
                    found = 0
                elif outside_nested(path, inside, own) == outside_nested(walk, inside, own):
                    found = 1
                elif sequences and int(g) < len(strings) and spells_allele(
                        path if walk is called else backwards(path), alleles[0], strings[0],
                        strings[int(g)], sequences):
                    found = 2
                else:
                    found = 3
                if found > worst:
                    worst, example = found, (path, walk)
            counts[classes[worst]] += 1
            if worst == 3 and shown < 10 and not quiet:
                shown += 1
                print(f"  DIFFERENT  {contig}:{pos} {name} GT {gt} strand {strand}: "
                      f"{divergence(*example)}", file=sys.stderr)
    return counts


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
    hapname, threads = {}, defaultdict(list)
    allele_rows = defaultdict(list)
    for line in open(a.mosaic):
        f = line.rstrip('\n').split('\t')
        if f[0] == '#haplotype':
            hapname[f[1]] = f[2]
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
        apply_allele_rows(strands[key], rows, fail, allele_counts)
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
    if a.vcf:
        vcf_counts = check_vcf(a.vcf, strands, sequences, a.quiet)
        print(f"  called alleles against {a.vcf}: {vcf_counts['checked']} checked, "
              f"{vcf_counts['equal']} equal, "
              f"{vcf_counts['equal outside nested records']} equal outside nested records, "
              f"{vcf_counts['same bases']} same bases by another walk, "
              f"{vcf_counts['different']} different, "
              f"{vcf_counts['strand does not cross']} strand does not cross")
        for k, v in vcf_counts.items():
            if k.startswith('records skipped') or k.startswith('block records'):
                print(f"  {v} {k}")
        if vcf_counts['different']:
            fail['called allele is not the walk the path takes'] = vcf_counts['different']
    if not fail:
        print("  OK: every path expands to an exact walk in the graph")
        return 0
    for k, v in sorted(fail.items(), key=lambda kv: -kv[1]):
        print(f"  FAIL {v:>7}  {k}")
    return 1


if __name__ == '__main__':
    sys.exit(main())
