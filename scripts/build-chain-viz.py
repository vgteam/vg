#!/usr/bin/env python3
"""
Build an interactive SVG visualization of a chaining problem using D3.js.

Usage: build-chain-viz.py <data_directory> <chain_number> <output.svg>
"""

import sys
import os
import glob
import json
import gzip
import base64
from collections import defaultdict, Counter


def parse_seeds_file(filepath):
    """
    Parse a chain seeds file.
    Returns list of dicts with keys: read_pos, ref_name, ref_pos, strand, seed_num, seed_id
    """
    seeds = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            seeds.append({
                'read_pos': int(parts[0]),
                'ref_name': parts[1],
                'ref_pos': int(parts[2]),
                'strand': parts[3],
                'seed_num': int(parts[4]),
                'seed_id': parts[5]
            })
    return seeds


def parse_chaindump_file(filepath):
    """
    Parse a chaindump file.
    Returns list of dicts with keys: source_id, dest_id, score (signed)
    """
    transitions = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            transitions.append({
                'source_id': parts[0],
                'dest_id': parts[1],
                'score': int(parts[2])
            })
    return transitions


def find_seeds_file(data_dir, chain_num):
    """Find the seeds file for the given chain number."""
    pattern = os.path.join(data_dir, f'chain{chain_num}-seeds*.tsv')
    matches = glob.glob(pattern)
    if not matches:
        return None
    return matches[0]


def find_all_chaindump_files(data_dir):
    """Find all chaindump files in the directory."""
    pattern = os.path.join(data_dir, 'chaindump*.tsv')
    return glob.glob(pattern)


def generate_svg(seeds, transitions, output_path):
    """Generate an SVG file with embedded D3.js visualization."""

    # Filter seeds to the most common reference path
    ref_name_counts = Counter(s['ref_name'] for s in seeds)
    if ref_name_counts:
        best_ref_name = ref_name_counts.most_common(1)[0][0]
        seeds = [s for s in seeds if s['ref_name'] == best_ref_name]

    # Build seed lookup by ID
    seed_by_id = {s['seed_id']: s for s in seeds}
    seed_ids = set(seed_by_id.keys())

    # Index seeds by ID
    seed_id_to_index = {s['seed_id']: i for i, s in enumerate(seeds)}

    # Filter transitions to those where both endpoints are in our seeds
    relevant_transitions = [t for t in transitions
                            if t['source_id'] in seed_ids and t['dest_id'] in seed_ids]

    # Compute max positive score per destination
    max_score_by_dest = defaultdict(lambda: float('-inf'))
    for t in relevant_transitions:
        if t['score'] > 0 and t['score'] > max_score_by_dest[t['dest_id']]:
            max_score_by_dest[t['dest_id']] = t['score']

    # Build seeds data with score info
    seeds_data = []
    for i, s in enumerate(seeds):
        max_score = max_score_by_dest.get(s['seed_id'], 0)
        seeds_data.append({
            'index': i,
            'seed_num': s['seed_num'],
            'read_pos': s['read_pos'],
            'ref_pos': s['ref_pos'],
            'strand': s['strand'],
            'seed_id': s['seed_id'],
            'max_score': max_score if max_score > float('-inf') else 0
        })

    # Build transitions data with fraction of max
    transitions_data = []
    for t in relevant_transitions:
        max_score = max_score_by_dest[t['dest_id']]
        fraction = max(0, t['score'] / max_score) if max_score > 0 else 0
        transitions_data.append({
            'source_id': t['source_id'],
            'source_index': seed_id_to_index.get(t['source_id'], -1),
            'dest_id': t['dest_id'],
            'dest_index': seed_id_to_index.get(t['dest_id'], -1),
            'score': t['score'],
            'fraction': fraction,
            'is_max': (t['score'] == max_score)
        })

    # Compress JSON data
    seeds_compressed = base64.b64encode(gzip.compress(
        json.dumps(seeds_data).encode('utf-8'))).decode('ascii')
    transitions_compressed = base64.b64encode(gzip.compress(
        json.dumps(transitions_data).encode('utf-8'))).decode('ascii')

    # CSS and JS are plain strings (no f-string brace escaping needed).
    # Only the data script uses f-string interpolation.
    css = '''  <style>
    .seed {
      fill: steelblue;
      stroke: #333;
      stroke-width: 1px;
      cursor: pointer;
    }
    .seed.on-best-chain {
      fill: #90EE90;
      stroke: darkgreen;
      stroke-width: 2px;
    }
    .seed:hover, .seed.hovered, .seed.on-best-chain:hover, .seed.on-best-chain.hovered {
      r: 8;
    }
    .seed.selected {
      fill: orange;
      r: 8;
    }
    .seed.on-traceback {
      fill: #ff6666;
      stroke: darkred;
      stroke-width: 2px;
    }
    .transition {
      stroke-width: 1.5px;
      fill: none;
      opacity: 0.6;
    }
    .transition.highlighted {
      stroke-width: 3px;
      opacity: 1;
    }
    .transition.to-selected {
      stroke: black !important;
      stroke-width: 4px !important;
      opacity: 1 !important;
    }
    .transition.secondary {
      stroke-width: 1px;
      opacity: 0.4;
    }
    .transition.traceback {
      stroke: red;
      stroke-width: 3px;
      opacity: 1;
    }
    .axis text {
      font-family: sans-serif;
      font-size: 12px;
    }
    .axis path, .axis line {
      stroke: #333;
    }
    .axis-label {
      font-family: sans-serif;
      font-size: 14px;
      font-weight: bold;
    }
  </style>'''

    data_script = (
        '  <script type="text/javascript">\n'
        f'    var SEEDS_COMPRESSED = "{seeds_compressed}";\n'
        f'    var TRANSITIONS_COMPRESSED = "{transitions_compressed}";\n'
        '  </script>'
    )

    js = '''  <script type="text/javascript">
    <![CDATA[

    /** Decompress a gzip+base64-encoded JSON blob using the Compression Streams API. */
    async function decompressData(base64Data) {
      const binaryString = atob(base64Data);
      const bytes = new Uint8Array(binaryString.length);
      for (let i = 0; i < binaryString.length; i++) {
        bytes[i] = binaryString.charCodeAt(i);
      }
      const stream = new Blob([bytes]).stream();
      const decompressedStream = stream.pipeThrough(new DecompressionStream('gzip'));
      const decompressedBlob = await new Response(decompressedStream).blob();
      const text = await decompressedBlob.text();
      return JSON.parse(text);
    }

    document.addEventListener('DOMContentLoaded', async function() {

      // ---------------------------------------------------------------
      // Data Loading and Indexing
      //
      // Decompress embedded data and build maps for fast lookup of seeds
      // by ID, transitions by endpoint index, and the best-scoring
      // transition arriving at each seed.
      // ---------------------------------------------------------------

      const seeds = await decompressData(SEEDS_COMPRESSED);
      const transitions = await decompressData(TRANSITIONS_COMPRESSED);

      const seedById = new Map();
      seeds.forEach(s => seedById.set(s.seed_id, s));

      const transitionsByDestIndex = new Map();
      const transitionsBySourceIndex = new Map();
      transitions.forEach(t => {
        if (!transitionsByDestIndex.has(t.dest_index))
          transitionsByDestIndex.set(t.dest_index, []);
        transitionsByDestIndex.get(t.dest_index).push(t);
        if (!transitionsBySourceIndex.has(t.source_index))
          transitionsBySourceIndex.set(t.source_index, []);
        transitionsBySourceIndex.get(t.source_index).push(t);
      });

      const bestTransitionToSeed = new Map();
      transitions.forEach(t => {
        if (t.is_max) bestTransitionToSeed.set(t.dest_index, t);
      });

      /**
       * Look up a specific transition between two seeds. Used to find
       * connections from a hovered seed to the selected seed, including
       * negative-score transitions that aren't shown in secondary display.
       */
      function getTransitionFromTo(sourceIndex, destIndex) {
        const fromSource = transitionsBySourceIndex.get(sourceIndex) || [];
        return fromSource.find(t => t.dest_index === destIndex);
      }

      // ---------------------------------------------------------------
      // Scales and SVG Layout
      //
      // Set up 1:1 aspect ratio coordinate scales, axes, labels, clip
      // region, and the layered SVG structure. Layers are ordered
      // bottom-to-top: traceback, secondary transitions, hovered-to-
      // selected highlight, seeds.
      // ---------------------------------------------------------------

      const margin = {top: 40, right: 40, bottom: 60, left: 80};
      const svgWidth = 1200;
      const svgHeight = 800;
      const width = svgWidth - margin.left - margin.right;
      const height = svgHeight - margin.top - margin.bottom;
      const svg = d3.select('svg');

      let xExtent = d3.extent(seeds, d => d.ref_pos);
      let yExtent = d3.extent(seeds, d => d.read_pos);
      if (xExtent[0] === undefined) xExtent = [0, 1000];
      if (yExtent[0] === undefined) yExtent = [0, 1000];

      // Use the same data-per-pixel rate on both axes for 1:1 aspect ratio,
      // chosen so the data fits in whichever dimension is tighter.
      const xDataRange = xExtent[1] - xExtent[0] || 1000;
      const yDataRange = yExtent[1] - yExtent[0] || 1000;
      const dataPerPixel = Math.max(xDataRange / width, yDataRange / height);
      const xMid = (xExtent[0] + xExtent[1]) / 2;
      const yMid = (yExtent[0] + yExtent[1]) / 2;
      const initialXDomain = [xMid - width * dataPerPixel / 2, xMid + width * dataPerPixel / 2];
      const initialYDomain = [yMid - height * dataPerPixel / 2, yMid + height * dataPerPixel / 2];

      const xScale = d3.scaleLinear().domain(initialXDomain).range([0, width]);
      const yScale = d3.scaleLinear().domain(initialYDomain).range([height, 0]);
      const colorScale = d3.scaleSequential(d3.interpolateViridis).domain([0, 1]);

      svg.append('defs').append('clipPath')
        .attr('id', 'plot-clip')
        .append('rect')
        .attr('width', width)
        .attr('height', height);

      const g = svg.append('g')
        .attr('transform', `translate(${margin.left},${margin.top})`);

      const xAxis = d3.axisBottom(xScale).ticks(10);
      const yAxis = d3.axisLeft(yScale).ticks(10);

      const xAxisG = g.append('g')
        .attr('class', 'axis x-axis')
        .attr('transform', `translate(0,${height})`)
        .call(xAxis);

      const yAxisG = g.append('g')
        .attr('class', 'axis y-axis')
        .call(yAxis);

      g.append('text')
        .attr('class', 'axis-label')
        .attr('x', width / 2)
        .attr('y', height + 45)
        .attr('text-anchor', 'middle')
        .text('Reference Position');

      g.append('text')
        .attr('class', 'axis-label')
        .attr('transform', 'rotate(-90)')
        .attr('x', -height / 2)
        .attr('y', -60)
        .attr('text-anchor', 'middle')
        .text('Read Position');

      const plotArea = g.append('g')
        .attr('clip-path', 'url(#plot-clip)');

      plotArea.append('rect')
        .attr('width', width)
        .attr('height', height)
        .attr('fill', 'white');

      const zoomG = plotArea.append('g');
      const tracebackLayer = zoomG.append('g');
      const secondaryTransitionLayer = zoomG.append('g');
      const hoveredToSelectedLayer = zoomG.append('g');
      const seedLayer = zoomG.append('g');

      // Current zoom-adjusted scales; start as the base scales and get
      // replaced on each zoom event.
      let currentXScale = xScale;
      let currentYScale = yScale;

      // ---------------------------------------------------------------
      // Helpers
      //
      // Shared utilities for positioning transition lines, formatting
      // tooltips, and computing zoom transforms.
      // ---------------------------------------------------------------

      /**
       * Set x1/y1/x2/y2 on a D3 selection of <line> elements whose data
       * objects have source_id and dest_id fields. Uses the current zoom-
       * adjusted scales, so call this both when creating lines and when
       * the zoom changes. Use as: selection.call(positionLines).
       */
      function positionLines(selection) {
        selection.each(function(d) {
          const src = seedById.get(d.source_id);
          const dst = seedById.get(d.dest_id);
          if (src && dst) {
            d3.select(this)
              .attr('x1', currentXScale(src.ref_pos))
              .attr('y1', currentYScale(src.read_pos))
              .attr('x2', currentXScale(dst.ref_pos))
              .attr('y2', currentYScale(dst.read_pos));
          }
        });
      }

      /**
       * Base tooltip text for a seed point, showing its number, strand,
       * ID, positions, and best achievable score. Extended by
       * updateSeedTooltip() when there is an active selection.
       */
      function seedTooltipText(d) {
        return `Seed #${d.seed_num} (${d.strand})\\nSeed ID: ${d.seed_id}\\nRead: ${d.read_pos}\\nRef: ${d.ref_pos}\\nMax score: ${d.max_score}`;
      }

      /**
       * Tooltip text for a transition line, showing its DP score and what
       * fraction of the best score to this destination it represents.
       */
      function transitionTooltipText(d) {
        return `Score: ${d.score}\\nFraction of max: ${(d.fraction * 100).toFixed(1)}%`;
      }

      /**
       * Apply a D3 zoom transform that centers and scales the viewport
       * to show seedsToFit, with padding as a fraction of the data range
       * added on each side (e.g. 0.5 adds 50% of the range as margin).
       * Picks whichever axis is tighter to determine the zoom level, so
       * all seeds are visible regardless of aspect ratio.
       */
      function zoomToFit(seedsToFit, padding) {
        if (seedsToFit.length === 0) return;
        const [xMin, xMax] = d3.extent(seedsToFit, s => s.ref_pos);
        const [yMin, yMax] = d3.extent(seedsToFit, s => s.read_pos);
        const xPad = (xMax - xMin) * padding || 100;
        const yPad = (yMax - yMin) * padding || 100;
        const targetXRange = (xMax - xMin) + 2 * xPad;
        const targetYRange = (yMax - yMin) + 2 * yPad;
        const kx = (initialXDomain[1] - initialXDomain[0]) / targetXRange;
        const ky = (initialYDomain[1] - initialYDomain[0]) / targetYRange;
        const k = Math.min(kx, ky);
        const cx = (xMin + xMax) / 2;
        const cy = (yMin + yMax) / 2;
        const tx = width / 2 - k * xScale(cx);
        const ty = height / 2 - k * yScale(cy);
        svg.call(zoom.transform, d3.zoomIdentity.translate(tx, ty).scale(k));
      }

      // ---------------------------------------------------------------
      // Interaction
      //
      // State tracking and functions for the DP traceback path, secondary
      // transition display on hover, hovered-to-selected highlighting,
      // and dynamic seed tooltips.
      // ---------------------------------------------------------------

      let selectedSeed = null;
      let hoveredSeed = null;
      let tracebackSeedIndices = new Set();

      /**
       * Follow the chain of best-scoring transitions backward from
       * seedIndex. Returns the DP traceback path as an array of
       * transition objects from the destination back to the chain origin.
       */
      function computeTraceback(seedIndex) {
        const path = [];
        const visited = new Set();
        let currentIndex = seedIndex;
        while (currentIndex !== undefined && !visited.has(currentIndex)) {
          visited.add(currentIndex);
          const bestTrans = bestTransitionToSeed.get(currentIndex);
          if (bestTrans) {
            path.push(bestTrans);
            currentIndex = bestTrans.source_index;
          } else {
            break;
          }
        }
        return path;
      }

      /**
       * Display the optimal traceback path from seedIndex as red lines,
       * and mark the intermediate seeds with on-traceback styling. Also
       * populates tracebackSeedIndices for use by other code (e.g. to
       * identify the best chain for persistent highlighting).
       */
      function showTraceback(seedIndex) {
        const path = computeTraceback(seedIndex);
        tracebackSeedIndices.clear();
        path.forEach(t => {
          tracebackSeedIndices.add(t.source_index);
          tracebackSeedIndices.add(t.dest_index);
        });

        seedCircles.classed('on-traceback', d => tracebackSeedIndices.has(d.index) && d.index !== seedIndex);

        const lines = tracebackLayer.selectAll('.traceback')
          .data(path, d => d.source_id + '->' + d.dest_id);

        lines.enter()
          .append('line')
          .attr('class', 'transition traceback')
          .call(positionLines)
          .append('title')
          .text(transitionTooltipText);

        lines.exit().remove();
      }

      function hideTraceback() {
        tracebackLayer.selectAll('.traceback').remove();
        tracebackSeedIndices.clear();
        seedCircles.classed('on-traceback', false);
      }

      /**
       * Show all positive-score transitions arriving at seedIndex as thin
       * lines colored by their fraction of the best score (viridis scale,
       * red for the best). Negative-score transitions are excluded here
       * but can still appear via highlightTransitionToSelected.
       */
      function showSecondaryTransitions(seedIndex) {
        const allTrans = transitionsByDestIndex.get(seedIndex) || [];
        const positiveTrans = allTrans.filter(t => t.score > 0);
        const className = 'secondary-dest-' + seedIndex;

        secondaryTransitionLayer.selectAll('.' + className)
          .data(positiveTrans, d => d.source_id + '->' + d.dest_id)
          .enter()
          .append('line')
          .attr('class', 'transition secondary ' + className)
          .call(positionLines)
          .attr('stroke', d => d.is_max ? 'red' : colorScale(d.fraction))
          .append('title')
          .text(transitionTooltipText);
      }

      function hideSecondaryTransitions(seedIndex) {
        secondaryTransitionLayer.selectAll('.secondary-dest-' + seedIndex).remove();
      }

      /**
       * If a transition exists from hoveredIndex to the currently selected
       * seed, draw it as a thick black line. This includes negative-score
       * transitions that don't appear in the secondary display, so the
       * user can discover connections the DP considered unfavorable.
       */
      function highlightTransitionToSelected(hoveredIndex) {
        if (!selectedSeed || hoveredIndex === selectedSeed.index) return;
        const trans = getTransitionFromTo(hoveredIndex, selectedSeed.index);
        if (!trans) return;

        hoveredToSelectedLayer.append('line')
          .datum(trans)
          .attr('class', 'transition to-selected')
          .call(positionLines)
          .append('title')
          .text(transitionTooltipText);
      }

      function unhighlightTransitionToSelected() {
        hoveredToSelectedLayer.selectAll('*').remove();
      }

      /**
       * Rebuild a seed's SVG <title> tooltip. When another seed is
       * selected, appends the transition score from this seed to the
       * selected seed (if one exists) and the insertion/deletion offset
       * between them, so the user can evaluate potential connections.
       */
      function updateSeedTooltip(circle, d) {
        let text = seedTooltipText(d);
        if (selectedSeed && selectedSeed !== d) {
          const trans = getTransitionFromTo(d.index, selectedSeed.index);
          if (trans) {
            text += `\\nScore to selected: ${trans.score}`;
          }
          const refDiff = Math.abs(d.ref_pos - selectedSeed.ref_pos);
          const readDiff = Math.abs(d.read_pos - selectedSeed.read_pos);
          const offset = Math.abs(refDiff - readDiff);
          if (readDiff > refDiff) {
            text += `\\nOffset: ${offset}bp INS`;
          } else if (refDiff > readDiff) {
            text += `\\nOffset: ${offset}bp DEL`;
          } else {
            text += '\\nOffset: 0bp';
          }
        }
        circle.select('title').text(text);
      }

      // ---------------------------------------------------------------
      // Seeds
      //
      // Create seed circles and attach mouse/click event handlers.
      // Hover shows secondary transitions and highlights connection to
      // selected seed. Click selects/deselects and shows traceback.
      // ---------------------------------------------------------------

      const seedCircles = seedLayer.selectAll('.seed')
        .data(seeds)
        .enter()
        .append('circle')
        .attr('class', 'seed')
        .attr('cx', d => xScale(d.ref_pos))
        .attr('cy', d => yScale(d.read_pos))
        .attr('r', 5)
        .on('mouseover', function(event, d) {
          hoveredSeed = d;
          if (selectedSeed !== d) {
            d3.select(this).classed('hovered', true);
          }
          updateSeedTooltip(d3.select(this), d);
          showSecondaryTransitions(d.index);
          highlightTransitionToSelected(d.index);
        })
        .on('mouseout', function(event, d) {
          hoveredSeed = null;
          if (selectedSeed !== d) {
            d3.select(this).classed('hovered', false);
            hideSecondaryTransitions(d.index);
          }
          unhighlightTransitionToSelected();
        })
        .on('click', function(event, d) {
          const prevSelected = selectedSeed;
          if (selectedSeed === d) {
            selectedSeed = null;
            d3.select(this).classed('selected', false);
            hideTraceback();
            if (hoveredSeed !== d) {
              hideSecondaryTransitions(d.index);
            }
          } else {
            if (prevSelected) {
              seedCircles.filter(s => s.index === prevSelected.index)
                .classed('selected', false)
                .classed('hovered', false);
              if (hoveredSeed === null || hoveredSeed.index !== prevSelected.index) {
                hideSecondaryTransitions(prevSelected.index);
              }
            }
            selectedSeed = d;
            d3.select(this).classed('selected', true);
            showTraceback(d.index);
            showSecondaryTransitions(d.index);
          }
        });

      seedCircles.append('title').text(seedTooltipText);

      // ---------------------------------------------------------------
      // Zoom
      //
      // Pan/zoom with a single scale factor for both axes, maintaining
      // 1:1 aspect ratio. On each zoom event, update the scales, reposition
      // axes, seeds, and all transition lines.
      // ---------------------------------------------------------------

      const zoom = d3.zoom()
        .scaleExtent([0.1, 10000])
        .on('zoom', function(event) {
          currentXScale = event.transform.rescaleX(xScale);
          currentYScale = event.transform.rescaleY(yScale);
          xAxisG.call(xAxis.scale(currentXScale));
          yAxisG.call(yAxis.scale(currentYScale));
          seedCircles
            .attr('cx', d => currentXScale(d.ref_pos))
            .attr('cy', d => currentYScale(d.read_pos));
          zoomG.selectAll('.transition').call(positionLines);
        });

      svg.call(zoom);

      // ---------------------------------------------------------------
      // Initialization
      //
      // Auto-select the highest-scoring seed, show its traceback, mark
      // the optimal chain with persistent green styling, and zoom to fit.
      // If no seed has a positive score, just zoom to show all seeds.
      // ---------------------------------------------------------------

      const bestSeed = seeds.reduce((best, s) => (s.max_score > best.max_score ? s : best), seeds[0]);
      if (bestSeed && bestSeed.max_score > 0) {
        selectedSeed = bestSeed;
        seedCircles.filter(s => s.index === bestSeed.index).classed('selected', true);
        showTraceback(bestSeed.index);
        const bestChainSeedIndices = new Set(tracebackSeedIndices);
        seedCircles.classed('on-best-chain', d => bestChainSeedIndices.has(d.index));
        zoomToFit(seeds.filter(s => tracebackSeedIndices.has(s.index)), 0.05);
      } else {
        zoomToFit(seeds, 0.05);
      }
    });
    ]]>
  </script>'''

    svg_content = '\n'.join([
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<svg xmlns="http://www.w3.org/2000/svg"'
        ' xmlns:xlink="http://www.w3.org/1999/xlink"'
        ' width="1200" height="800">',
        css,
        '  <script xlink:href="https://d3js.org/d3.v7.js"></script>',
        data_script,
        js,
        '</svg>',
        ''
    ])

    with open(output_path, 'w') as f:
        f.write(svg_content)


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <data_directory> <chain_number> <output.svg>", file=sys.stderr)
        sys.exit(1)

    data_dir = sys.argv[1]
    chain_num = int(sys.argv[2])
    output_path = sys.argv[3]

    if not os.path.isdir(data_dir):
        print(f"Error: {data_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    # Find and parse seeds file
    seeds_file = find_seeds_file(data_dir, chain_num)
    if not seeds_file:
        print(f"Error: No seeds file found for chain {chain_num}", file=sys.stderr)
        sys.exit(1)

    seeds = parse_seeds_file(seeds_file)
    if not seeds:
        print(f"Warning: No seeds found in {seeds_file}", file=sys.stderr)

    print(f"Loaded {len(seeds)} seeds from {seeds_file}")

    # Find and parse all chaindump files
    chaindump_files = find_all_chaindump_files(data_dir)
    all_transitions = []
    for cf in chaindump_files:
        transitions = parse_chaindump_file(cf)
        all_transitions.extend(transitions)
        print(f"Loaded {len(transitions)} transitions from {cf}")

    print(f"Total transitions: {len(all_transitions)}")

    # Generate SVG
    generate_svg(seeds, all_transitions, output_path)
    print(f"Generated {output_path}")


if __name__ == '__main__':
    main()
