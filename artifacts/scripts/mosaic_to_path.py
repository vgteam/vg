#!/usr/bin/env python3
"""Rebuild the strands of a mosaic and check them against the graph and the VCF.

Inputs: the mosaic of `vg call --mosaic-out`, with or without `--mosaic-exact`; the graph's
paths, from `vg paths -x GRAPH -A`; optionally the graph, from `vg view -g GRAPH`, for its edges
and its nodes' bases; and optionally the VCF of the same run, which needs the graph.

A strand is rebuilt as the file says. Each fragment's `H` rows are expanded along their
haplotypes and joined on the nodes they share. Then every `A` row of the fragment is put on that
walk: at each visit of its start node, the visits up to the next visit of its end node become the
row's walk, and each step takes its entry of `bases`. In an entry `*` stands for the bases the
visit already has, and `+` is followed by the bases that lie after the step on no node, which
join it to the next step without an edge.

Checks, any of which makes the exit status non-zero:
  1. every `H` row names a haplotype the GAF carries, with its two end nodes on one of its walks,
     and consecutive rows of a fragment share the node they meet at
  2. every `A` row has a stretch to go on and one entry of bases per step, and no two rows of a
     strand change the same visit
  3. with --gfa, every step of a rebuilt strand is an edge of the graph or is joined by a `+`
  4. with --vcf, every strand has the bases the records give it:
     - a site, whose records' IDs name its two bound nodes (`>a>b`, or `>a>b_N`): the reference
       from bound node to bound node, with the strand's alleles of the site's records put in, is
       what the strand has there. The bound nodes count with their own bases, since a record can
       lie on them and the changes to them are another site's.
     - the nodes that the records of edit sites (`op...`) change: the reference over them, with
       the strand's alleles put in, is what the strand's visits of them have.
     - where an `A` row makes a strand differ from the reference, a record lies there.

What is not compared is counted, by reason. A site with the records of other sites inside it is
one such class: its record spells an allele once, with the alleles one strand has at the sites
inside, so only those are compared. The VCF must be written without a node translation.
"""
import sys, argparse, re, bisect
from collections import defaultdict

STEP = re.compile(r'([<>])(\d+)')
SNARL_ID = re.compile(r'([<>]\d+[<>]\d+)(_\d+)?$')
COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def parse_walk(text):
    """'>1<2' -> [(1, False), (2, True)]."""
    return [(int(i), o == '<') for o, i in STEP.findall(text)]


def read_walks(gaf):
    """path name -> list of (node_id, is_reverse), from the GAF path field."""
    walks, shared = {}, {}
    for line in open(gaf):
        f = line.rstrip('\n').split('\t')
        if line.startswith('@') or len(f) < 6:
            continue
        # Haplotypes share most of their steps, so each distinct step is kept once.
        walks[f[0]] = [shared.setdefault(step, step) for step in parse_walk(f[5])]
    return walks


def read_graph(gfa):
    """(edges, sequences): each L-line in both orientations, and each node's bases."""
    edges, sequences = set(), {}
    for line in open(gfa):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'S' and f[1].isdigit():
            sequences[int(f[1])] = f[2]
        elif f[0] == 'L':
            a, b = (int(f[1]), f[2] == '-'), (int(f[3]), f[4] == '-')
            edges.add((a, b))
            edges.add(((b[0], not b[1]), (a[0], not a[1])))
    return edges, sequences


def haplotype_walks(walks, name, wanted, cache={}):
    """The walks of one haplotype, longest first, each with the places of the nodes in `wanted`."""
    if name not in cache:
        cache[name] = []
        for w in sorted((w for n, w in walks.items() if n == name or n.startswith(name + '#')
                         or n.split('[')[0] == name), key=len, reverse=True):
            places = defaultdict(list)
            for i, (node, _) in enumerate(w):
                if node in wanted:
                    places[node].append(i)
            cache[name].append((w, places))
    return cache[name]


def pick_copy(w, starts, ends, start_reverse, end_reverse, resume):
    """Pick the visits of a row's start and end nodes that the row means: those with the stated
    orientations, or both flipped where the haplotype is stored the other way, nearest to where
    the fragment's walk has reached. Returns (None, None) if there are none."""
    for flip in (False, True):
        pairs = [(abs(i - resume), abs(j - i), i, j)
                 for i in starts if (w[i][1] != start_reverse) == flip
                 for j in ends if (w[j][1] != end_reverse) == flip]
        if pairs:
            return min(pairs)[2:]
    return None, None


def expand_fragment(rows, contig, walks, row_ends, fail, counts):
    """The node walk that the `H` rows of one fragment spell."""
    path, resume = [], 0
    for f in rows:
        # A row gives an oriented node as id * 2 + is_reverse.
        (start, start_reverse), (end, end_reverse) = divmod(int(f[6]), 2), divmod(int(f[7]), 2)
        found = [] if f[8] == '*' else (haplotype_walks(walks, f[9], row_ends)
                                        or haplotype_walks(walks, f[9] + '#' + contig, row_ends))
        if not found:
            fail['row names no haplotype that the graph has'] += 1
            continue
        piece = None
        for w, places in found:
            i, j = pick_copy(w, places.get(start, ()), places.get(end, ()), start_reverse,
                             end_reverse, resume)
            if i is None:
                continue
            resume = j
            if i <= j:
                piece = w[i:j + 1]
            else:
                piece = [(n, not r) for n, r in reversed(w[j:i + 1])]
                counts['H rows that run against the direction their haplotype is stored in'] += 1
            if piece[0][1] != start_reverse or piece[-1][1] != end_reverse:
                fail["row's stated orientation is not the walk's"] += 1
            break
        if piece is None:
            fail['row ends are not both on one walk of its haplotype'] += 1
            continue
        if path and path[-1] != piece[0]:
            fail['rows do not share the node they meet at'] += 1
        path.extend(piece if not path else piece[1:])
    return path


def crossings(places, start, end):
    """Each visit of `start` with the next visit of `end`, the same visit where the two are one
    node, unless the walk comes back to the node of `start` on the way. `places` has the visits
    of each oriented node."""
    again, turned = places.get(end, ()), (start[0], not start[1])
    for i in places.get(start, ()):
        k = bisect.bisect_left(again, i if start == end else i + 1)
        if k < len(again) and not any(i < p < again[k] for step in (start, turned)
                                      for p in places.get(step, ())):
            yield i, again[k]


def oriented(step, sequences):
    """The bases of a node, read the way `step` reads it."""
    return sequences[step[0]][::-1].translate(COMPLEMENT) if step[1] else sequences[step[0]]


class Fragment:
    """One fragment of a strand: the visit of each node, with the bases a row gave it (None for
    the node's own), the bases a row put after it on no node (None unless a row joins it to the
    next visit), and whether a row made it."""

    def __init__(self, path, rows, fail, counts):
        places = defaultdict(list)
        for i, step in enumerate(path):
            places[step].append(i)
        # Each row's stretches of the walk, as (first visit, last visit, walk, entries).
        stretches = []
        for f in rows:
            walk = parse_walk(f[6])
            entries = ['*'] * len(walk) if f[7] == '.' else f[7].split(',')
            counts['A rows of one node' if len(walk) == 1 else 'A rows of sites'] += 1
            if not walk or [walk[0], walk[-1]] != [divmod(int(f[4]), 2), divmod(int(f[5]), 2)]:
                fail['A row whose walk does not run from its start node to its end node'] += 1
                continue
            if len(entries) != len(walk) or not all(
                    re.fullmatch(r'(\*|[ACGTNacgtn]*)(\+[ACGTNacgtn]*)?', e) for e in entries):
                fail['A row whose bases are not one entry for each step of its walk'] += 1
                continue
            found = [(i, j, walk, entries) for i, j in crossings(places, walk[0], walk[-1])]
            if not found:
                fail['A row with no stretch of its fragment to go on'] += 1
            stretches += found
        self.nodes, self.bases, self.after, self.made = [], [], [], []
        at, site_end, node_row = 0, -1, -1
        for first, last, walk, entries in sorted(stretches, key=lambda s: s[:2]):
            if first < site_end or (first == last and first == node_row):
                fail['A rows that change the same visit'] += 1
                continue
            self.extend(path[at:first + 1])
            at = max(at, first + 1)
            # The visits the row's entries are for.
            steps = [len(self.nodes) - 1]
            if first == last:
                node_row = first
            else:
                steps += range(len(self.nodes), len(self.nodes) + len(walk) - 2)
                self.extend(walk[1:-1], True)
                self.extend(path[last:last + 1])
                steps.append(len(self.nodes) - 1)
                at, site_end = last + 1, last
            for k, entry in zip(steps, entries):
                own, joined, after = entry.partition('+')
                if (own != '*' and self.bases[k] is not None) or (joined
                                                                  and self.after[k] is not None):
                    fail['A rows that change the same visit'] += 1
                if own != '*':
                    self.bases[k] = own
                if joined:
                    self.after[k] = after
                self.made[k] = self.made[k] or entry != '*'
        self.extend(path[at:])
        self.places = defaultdict(list)
        for i, step in enumerate(self.nodes):
            self.places[step].append(i)

    def extend(self, steps, made=False):
        self.nodes += steps
        self.bases += [None] * len(steps)
        self.after += [None] * len(steps)
        self.made += [made] * len(steps)

    def spell(self, i, sequences):
        """The bases of visit `i`."""
        return self.bases[i] if self.bases[i] is not None else oriented(self.nodes[i], sequences)

    def between(self, i, j, sequences):
        """The bases after visit `i` and before visit `j`."""
        return (self.after[i] or '') + ''.join(self.spell(k, sequences) + (self.after[k] or '')
                                               for k in range(i + 1, j))


class Reference:
    """One contig's reference path: where each of its steps starts, and the step of each node
    it visits once."""

    def __init__(self, walk, sequences):
        self.walk, self.sequences, self.starts, self.step = walk, sequences, [], {}
        at = 0
        for i, (node, _) in enumerate(walk):
            self.starts.append(at)
            at += len(sequences[node])
            self.step[node] = None if node in self.step else i

    def on(self, step):
        """The index of the path's step at `step`, if it visits that node once, the same way."""
        i = self.step.get(step[0])
        return i if i is not None and self.walk[i] == step else None

    def spell(self, i, j):
        """The bases of the path's steps `i` up to `j`."""
        return ''.join(oriented(s, self.sequences) for s in self.walk[i:j]).upper()


def read_vcf(vcf, counts):
    """contig -> (the records of each snarl by its bound nodes, the records of edit sites, the
    place, REF length and site of each record). A record is (place of first base, REF, ALTs, GT)."""
    out = defaultdict(lambda: (defaultdict(list), [], []))
    for line in open(vcf):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        record = (int(f[1]) - 1, f[3].upper(), f[4].upper().split(','), f[9].split(':')[0])
        sites, edits, footprints = out[f[0]]
        named = SNARL_ID.match(f[2])
        if named:
            sites[named.group(1)].append(record)
        elif f[2].startswith('op'):
            edits.append(record)
        else:
            counts['records not compared: the ID names neither a snarl nor an edit site'] += 1
        footprints.append((record[0], len(record[1]), f[2].split('_')[0]))
    return out


def change_of(record, strand):
    """The change a record makes on one strand, as (place, bases removed, bases put in, None),
    which removes and puts in nothing for the reference allele. Where the record does not say
    what the strand has, the last is the reason, and the change is the record's whole REF."""
    place, ref, alts, gt = record
    called = gt.split('|') if '|' in gt else gt.split('/')
    allele = called[strand] if strand < len(called) else called[0]
    if '|' not in gt and len(set(called)) > 1:
        return place, len(ref), '', 'a record there is not phased'
    if allele == '.':
        return place, len(ref), '', 'a record there has no allele for the strand'
    alt, k = (alts[int(allele) - 1] if allele != '0' else ref), 0
    if not re.fullmatch('[ACGTN]+', alt):
        return place, len(ref), '', 'a record there has a symbolic allele'
    while k < len(ref) and k < len(alt) and ref[k] == alt[k]:
        k += 1
    return place + k, len(ref) - k, alt[k:], None


def with_changes(bases, start, made):
    """`bases`, which start at place `start` of the reference, with the changes made, or None
    if two of them take the same base or one lies outside the bases."""
    out, at = [], start
    for place, removed, added, _ in sorted(change for change in made if change[1] or change[2]):
        if place < at or place + removed > start + len(bases):
            return None
        out += [bases[at - start:place - start], added]
        at = place + removed
    return ''.join(out) + bases[at - start:]


def check_strand(fragments, strand, vcf, reference, sequences, counts):
    """Compare one strand with the records of its contig. Returns the differences."""
    sites, edits, footprints = vcf
    different = []
    differs = different.append

    # Sites, from bound node to bound node.
    places = sorted((place, name) for place, _, name in footprints)
    for name, records in sites.items():
        first, last = bounds = parse_walk(name)
        i, j = reference.on(first), reference.on(last)
        if i is None or j is None or j <= i:
            counts['site strands not compared: the reference path does not give the site'] += 1
            continue
        inside = places[bisect.bisect_left(places, (reference.starts[i + 1], '')):
                        bisect.bisect_left(places, (reference.starts[j], ''))]
        if any(other != name for _, other in inside):
            counts['site strands not compared: records of other sites lie inside the site'] += 1
            continue
        made = [change_of(record, strand) for record in records]
        reasons = [change[3] for change in made if change[3]]
        if reasons:
            counts[f'site strands not compared: {reasons[0]}'] += 1
            continue
        expected = with_changes(reference.spell(i, j + 1), reference.starts[i], made)
        found = [(f, c) for f in fragments for c in crossings(f.places, first, last)]
        if not found:
            counts['site strands not compared: the strand does not cross the site from start to '
                   'end'] += 1
            continue
        counts['site strands checked'] += 1
        ends = [oriented(bound, sequences).upper() for bound in bounds]
        # Records that take the same base, or lie outside the site, give no bases to have.
        if expected is not None and all(
                ends[0] + f.between(a, b, sequences).upper() + ends[1] == expected
                for f, (a, b) in found):
            counts['site strands equal'] += 1
        else:
            differs(f"{name}: the strand does not have the bases its records give")

    # The nodes that the records of edit sites change, taken together where they share a node.
    groups = []
    # A change that takes and puts in nothing is the reference allele.
    for change in sorted(filter(lambda c: c[1] or c[2], (change_of(r, strand) for r in edits))):
        place, removed = change[:2]
        # The steps of the reference with the first and the last base the change takes. Inserted
        # bases are with the node of the base before them.
        steps = [bisect.bisect_right(reference.starts, at) - 1
                 for at in ((place, place + removed - 1) if removed else (place - 1, place - 1))]
        if groups and steps[0] <= groups[-1][1]:
            groups[-1][1] = max(groups[-1][1], steps[1])
            groups[-1][2].append(change)
        else:
            groups.append([steps[0], steps[1], [change]])
    for i, j, group in groups:
        visits = [[(f, v) for f in fragments for v in f.places.get(step, ())]
                  for step in reference.walk[i:j + 1]]
        expected = with_changes(reference.spell(i, j + 1), reference.starts[i], group)
        reasons = [change[3] for change in group if change[3]]
        if reasons:
            counts[f'edit site nodes not compared: {reasons[0]}'] += 1
        elif any(len(visit) != 1 for visit in visits):
            counts['edit site nodes not compared: the strand does not visit each node once'] += 1
        elif expected is None:
            counts['edit site nodes not compared: two records there take the same base'] += 1
        else:
            counts['edit site nodes checked'] += 1
            # Bases put after the last node can be the edit site's, written as lying after it.
            (f, v), = visits[-1]
            spelled = ''.join(f.spell(v, sequences) for (f, v), in visits).upper()
            if expected in (spelled, spelled + (f.after[v] or '').upper()):
                counts['edit site nodes equal'] += 1
            else:
                differs(f"node {reference.walk[i][0]}: the strand does not have the bases the "
                        f"records of edit sites give")

    # Where an A row makes the strand differ from the reference, a record lies: its REF is
    # there, past its first base, which an insertion keeps, or its site is, between the bounds.
    spans = []
    for place, length, name in footprints:
        ends = [reference.on(bound) for bound in parse_walk(name)]
        if len(ends) == 2 and None not in ends and ends[0] < ends[1]:
            spans.append((reference.starts[ends[0] + 1], reference.starts[ends[1]]))
        spans.append((place + 1, place + length))
    covered = []
    for lo, hi in sorted(spans):
        if covered and lo <= covered[-1][1]:
            covered[-1][1] = max(covered[-1][1], hi)
        else:
            covered.append([lo, hi])

    def recorded(lo, hi):
        k = bisect.bisect_right(covered, [hi, float('inf')]) - 1
        return k >= 0 and covered[k][1] >= lo

    for f in fragments:
        before = None
        for v, step in enumerate(f.nodes):
            r = reference.on(step)
            if r is None:
                continue
            lo, hi = reference.starts[r], reference.starts[r] + len(sequences[step[0]])
            if f.bases[v] is not None and f.bases[v].upper() != reference.spell(r, r + 1) \
                    and not recorded(lo, hi):
                differs(f"node {step[0]}: an A row changes it where no record lies")
            if before is not None and (v > before[0] + 1 or r != before[1] + 1
                                       or f.after[before[0]] is not None):
                if r <= before[1]:
                    counts['stretches not compared: the strand does not follow the reference '
                           'path there'] += 1
                elif f.between(before[0], v, sequences).upper() \
                        != reference.spell(before[1] + 1, r) and not recorded(before[2], lo):
                    if any(f.made[before[0]:v + 1]):
                        differs(f"between nodes {f.nodes[before[0]][0]} and {step[0]}: an A row "
                                f"makes the strand differ from the reference where no record lies")
                    else:
                        counts['stretches not compared: the H rows leave the reference where no '
                               'record lies, and no A row changes them'] += 1
            before = (v, r, hi)
    return different


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mosaic', required=True)
    ap.add_argument('--gaf', required=True)
    ap.add_argument('--gfa', help='optional; without it the steps are not checked to be edges')
    ap.add_argument('--vcf', help='optional; the calls the mosaic was written with')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args()

    walks = read_walks(a.gaf)
    edges, sequences = read_graph(a.gfa) if a.gfa else (None, None)
    h_rows, a_rows, reference_names = defaultdict(list), defaultdict(list), []
    for line in open(a.mosaic):
        f = line.rstrip('\n').split('\t')
        if f[0] == '#reference':
            reference_names += f[1:]
        elif f[0] == '#mosaic-version' and f[1] not in ('5', '6'):
            print(f"  FAIL  mosaic version {f[1]} is not one this script reads (5 or 6)")
            return 1
        elif f[0] in ('H', 'A'):
            (h_rows if f[0] == 'H' else a_rows)[(f[1], f[2], f[3])].append(f)

    fail, counts, strands, steps = defaultdict(int), defaultdict(int), defaultdict(list), 0
    row_ends = {int(f[c]) // 2 for rows in h_rows.values() for f in rows for c in (6, 7)}
    for key in sorted(set(h_rows) | set(a_rows)):
        contig, strand, fragment = key
        path = expand_fragment(h_rows[key], contig, walks, row_ends, fail, counts)
        built = Fragment(path, a_rows[key], fail, counts)
        strands[(contig, int(strand))].append(built)
        if not a.quiet:
            print(f"  {contig} strand {strand} fragment {fragment}: {len(h_rows[key])} H rows, "
                  f"{len(a_rows[key])} A rows, {len(built.nodes)} node steps")
        for i in range(len(built.nodes) - 1):
            steps += 1
            if built.after[i] is not None:
                counts['steps joined by an A row rather than an edge'] += 1
            elif edges is not None and (built.nodes[i], built.nodes[i + 1]) not in edges:
                fail['step is not an edge in the graph'] += 1

    print(f"\n  {len(h_rows)} fragments, {steps} steps "
          + (f"against {len(edges) // 2} edges" if edges is not None else "(no --gfa: unchecked)"))
    if a.vcf and sequences is None:
        ap.error("--vcf needs --gfa, for the nodes' bases")
    if a.vcf:
        records = read_vcf(a.vcf, counts)
        different = []
        for (contig, strand), fragments in sorted(strands.items()):
            name = next((n for n in reference_names if n.split('#')[-1] == contig and n in walks),
                        None)
            if name is None or contig not in records:
                counts['strands not compared: no reference path, or no records'] += 1
                continue
            different += [f"{contig} strand {strand} {what}" for what in check_strand(
                fragments, strand, records[contig], Reference(walks[name], sequences), sequences,
                counts)]
        for what in different if not a.quiet else ():
            print(f"  DIFFERENT  {what}", file=sys.stderr)
        print(f"  against {a.vcf}: {counts['site strands checked']} site strands checked, "
              f"{counts['site strands equal']} equal; {counts['edit site nodes checked']} "
              f"stretches of nodes of edit sites checked, {counts['edit site nodes equal']} "
              f"equal; {len(different)} different")
        if different:
            fail['strand does not have the bases the records give'] = len(different)
    for k, v in sorted(counts.items()):
        if not k.endswith((' checked', ' equal')):
            print(f"  {v} {k}")
    if not fail:
        print("  OK: every strand is a walk in the graph"
              + (", with the bases the records give" if a.vcf else ""))
        return 0
    for k, v in sorted(fail.items(), key=lambda kv: -kv[1]):
        print(f"  FAIL {v:>7}  {k}")
    return 1


if __name__ == '__main__':
    sys.exit(main())
