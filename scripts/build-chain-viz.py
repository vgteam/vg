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

    /**
     * Data model and query interface for the chaining problem.
     *
     * Holds the seed and transition arrays plus precomputed indexes for
     * fast lookup. Seeds are points in (ref_pos, read_pos) space;
     * transitions represent DP edges between them with cumulative scores.
     * The `fraction` and `is_max` fields on transitions are relative to
     * the best positive score arriving at each destination.
     */
    class ChainData {
      /**
       * Decode and gunzip a base64 blob via Compression Streams API.
       */
      static async decompress(base64Data) {
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

      /**
       * Async factory: decompresses both blobs, returns a new ChainData.
       */
      static async fromCompressed(seedsB64, transitionsB64) {
        const seeds = await ChainData.decompress(seedsB64);
        const transitions = await ChainData.decompress(transitionsB64);
        return new ChainData(seeds, transitions);
      }

      /**
       * Stores arrays and builds four index Maps.
       */
      constructor(seeds, transitions) {
        /// The raw seed array, for iteration and extent computation.
        this.seeds = seeds;
        /// The raw transition array, for iteration and extent computation.
        this.transitions = transitions;

        /// Map from seed_id to seed object.
        this._seedById = new Map();
        seeds.forEach(s => this._seedById.set(s.seed_id, s));

        /// Map from dest_index to array of transitions arriving there.
        this._transitionsByDestIndex = new Map();
        /// Map from source_index to array of transitions departing there.
        this._transitionsBySourceIndex = new Map();
        transitions.forEach(t => {
          if (!this._transitionsByDestIndex.has(t.dest_index))
            this._transitionsByDestIndex.set(t.dest_index, []);
          this._transitionsByDestIndex.get(t.dest_index).push(t);
          if (!this._transitionsBySourceIndex.has(t.source_index))
            this._transitionsBySourceIndex.set(t.source_index, []);
          this._transitionsBySourceIndex.get(t.source_index).push(t);
        });

        /// Map from dest_index to the is_max transition arriving there.
        this._bestTransitionToSeed = new Map();
        transitions.forEach(t => {
          if (t.is_max) this._bestTransitionToSeed.set(t.dest_index, t);
        });
      }

      /**
       * Look up a seed by its seed_id.
       */
      seedById(id) {
        return this._seedById.get(id);
      }

      /**
       * The is_max transition arriving at seedIndex, or undefined.
       */
      bestTransitionTo(seedIndex) {
        return this._bestTransitionToSeed.get(seedIndex);
      }

      /**
       * A specific transition from sourceIndex to destIndex, or undefined.
       */
      transitionFromTo(sourceIndex, destIndex) {
        const fromSource = this._transitionsBySourceIndex.get(sourceIndex) || [];
        return fromSource.find(t => t.dest_index === destIndex);
      }

      /**
       * All transitions arriving at seedIndex.
       */
      transitionsTo(seedIndex) {
        return this._transitionsByDestIndex.get(seedIndex) || [];
      }

      /**
       * All transitions departing sourceIndex.
       */
      transitionsFrom(sourceIndex) {
        return this._transitionsBySourceIndex.get(sourceIndex) || [];
      }
    }

    /**
     * Renders the chaining problem as an interactive SVG scatterplot.
     *
     * Seeds are circles; transitions are lines drawn in four z-ordered
     * layers (traceback, secondary, hovered-to-selected, seeds).
     *
     * Interaction model: hover a seed to see its incoming positive-score
     * transitions as thin colored lines; click to select and show the DP
     * traceback in red; hover another seed while one is selected to see
     * the direct connection (including negative scores) as a thick black
     * line. Secondary transitions are given a CSS class
     * `secondary-dest-{index}` so they can be efficiently removed when
     * the hover ends. Zoom maintains 1:1 aspect ratio via a single D3
     * zoom transform.
     */
    class ChainViz {
      constructor(data) {
        /// The ChainData instance providing seed/transition queries.
        this.data = data;
        /// The currently click-selected seed object, or null. Selecting a
        /// seed shows its traceback and enables hovered-to-selected highlighting.
        this.selectedSeed = null;
        /// The seed currently under the mouse, or null. Used to manage
        /// secondary transition display and avoid redundant hide/show.
        this.hoveredSeed = null;
        /// Set of seed indices along the current traceback path. Used both
        /// for styling and as the seed set for #zoomToFit on startup.
        this.tracebackSeedIndices = new Set();

        this.#setupScales();
        this.#setupSVG();
        this.#setupSeeds();
        this.#setupZoom();
        this.#initialize();
      }

      // -----------------------------------------------------------
      // Private setup methods
      // -----------------------------------------------------------

      /**
       * Compute 1:1 aspect ratio domains, create xScale/yScale/colorScale.
       */
      #setupScales() {
        const seeds = this.data.seeds;
        const margin = {top: 40, right: 40, bottom: 60, left: 80};
        const svgWidth = 1200;
        const svgHeight = 800;
        this._margin = margin;
        this._width = svgWidth - margin.left - margin.right;
        this._height = svgHeight - margin.top - margin.bottom;
        const width = this._width;
        const height = this._height;

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

        /// The domains of the base scales; needed by #zoomToFit to compute
        /// the zoom transform's scale factor from the ratio of domain widths.
        this.initialXDomain = [xMid - width * dataPerPixel / 2, xMid + width * dataPerPixel / 2];
        this.initialYDomain = [yMid - height * dataPerPixel / 2, yMid + height * dataPerPixel / 2];

        /// The unzoomed 1:1 aspect ratio scales; the D3 zoom transform is
        /// computed relative to these, so they must not be mutated.
        this.baseXScale = d3.scaleLinear().domain(this.initialXDomain).range([0, width]);
        this.baseYScale = d3.scaleLinear().domain(this.initialYDomain).range([height, 0]);

        /// Zoom-adjusted scales; updated on every zoom event and read by
        /// #positionLines and seed repositioning.
        this.currentXScale = this.baseXScale;
        this.currentYScale = this.baseYScale;

        /// Viridis color scale mapping transition fraction [0,1] to color.
        this.colorScale = d3.scaleSequential(d3.interpolateViridis).domain([0, 1]);
      }

      /**
       * Create clip path, axes, labels, layers.
       */
      #setupSVG() {
        const margin = this._margin;
        const width = this._width;
        const height = this._height;

        /// D3 selection of the root <svg> element; needed to attach and
        /// programmatically trigger the zoom behavior.
        this.svg = d3.select('svg');

        this.svg.append('defs').append('clipPath')
          .attr('id', 'plot-clip')
          .append('rect')
          .attr('width', width)
          .attr('height', height);

        const g = this.svg.append('g')
          .attr('transform', `translate(${margin.left},${margin.top})`);

        /// D3 axis generators; passed to xAxisG/yAxisG.call() with the
        /// current scale on zoom.
        this.xAxis = d3.axisBottom(this.baseXScale).ticks(10);
        this.yAxis = d3.axisLeft(this.baseYScale).ticks(10);

        /// D3 selections of the axis <g> elements; updated on zoom to
        /// redraw tick marks.
        this.xAxisG = g.append('g')
          .attr('class', 'axis x-axis')
          .attr('transform', `translate(0,${height})`)
          .call(this.xAxis);

        this.yAxisG = g.append('g')
          .attr('class', 'axis y-axis')
          .call(this.yAxis);

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

        /// The <g> inside the clip region that contains all layers; used
        /// in the zoom handler to reposition all transition lines at once.
        this.zoomG = plotArea.append('g');

        /// D3 selections of the <g> layers for each transition type.
        this.tracebackLayer = this.zoomG.append('g');
        this.secondaryTransitionLayer = this.zoomG.append('g');
        this.hoveredToSelectedLayer = this.zoomG.append('g');
        this._seedLayer = this.zoomG.append('g');
      }

      /**
       * Create seed circles, attach event handlers, set initial tooltips.
       */
      #setupSeeds() {
        const viz = this;
        const seeds = this.data.seeds;

        /// D3 selection of all seed <circle> elements; stored so interaction
        /// methods can restyle seeds (e.g. classed 'on-traceback') without
        /// re-selecting.
        this.seedCircles = this._seedLayer.selectAll('.seed')
          .data(seeds)
          .enter()
          .append('circle')
          .attr('class', 'seed')
          .attr('cx', d => viz.baseXScale(d.ref_pos))
          .attr('cy', d => viz.baseYScale(d.read_pos))
          .attr('r', 5)
          .on('mouseover', function(event, d) {
            viz.hoveredSeed = d;
            if (viz.selectedSeed !== d) {
              d3.select(this).classed('hovered', true);
            }
            viz.#updateSeedTooltip(d3.select(this), d);
            viz.#showSecondaryTransitions(d.index);
            viz.#highlightTransitionToSelected(d.index);
          })
          .on('mouseout', function(event, d) {
            viz.hoveredSeed = null;
            if (viz.selectedSeed !== d) {
              d3.select(this).classed('hovered', false);
              viz.#hideSecondaryTransitions(d.index);
            }
            viz.#unhighlightTransitionToSelected();
          })
          .on('click', function(event, d) {
            const prevSelected = viz.selectedSeed;
            if (viz.selectedSeed === d) {
              viz.selectedSeed = null;
              d3.select(this).classed('selected', false);
              viz.#hideTraceback();
              if (viz.hoveredSeed !== d) {
                viz.#hideSecondaryTransitions(d.index);
              }
            } else {
              if (prevSelected) {
                viz.seedCircles.filter(s => s.index === prevSelected.index)
                  .classed('selected', false)
                  .classed('hovered', false);
                if (viz.hoveredSeed === null || viz.hoveredSeed.index !== prevSelected.index) {
                  viz.#hideSecondaryTransitions(prevSelected.index);
                }
              }
              viz.selectedSeed = d;
              d3.select(this).classed('selected', true);
              viz.#showTraceback(d.index);
              viz.#showSecondaryTransitions(d.index);
            }
          });

        this.seedCircles.append('title').text(d => viz.#seedTooltipText(d));
      }

      /**
       * Create zoom behavior, attach to svg.
       */
      #setupZoom() {
        const viz = this;

        /// The D3 zoom behavior object; stored so #zoomToFit can call
        /// svg.call(zoom.transform, ...).
        this.zoom = d3.zoom()
          .scaleExtent([0.1, 10000])
          .on('zoom', function(event) {
            viz.currentXScale = event.transform.rescaleX(viz.baseXScale);
            viz.currentYScale = event.transform.rescaleY(viz.baseYScale);
            viz.xAxisG.call(viz.xAxis.scale(viz.currentXScale));
            viz.yAxisG.call(viz.yAxis.scale(viz.currentYScale));
            viz.seedCircles
              .attr('cx', d => viz.currentXScale(d.ref_pos))
              .attr('cy', d => viz.currentYScale(d.read_pos));
            viz.zoomG.selectAll('.transition').call(viz.#positionLinesFn());
          });

        this.svg.call(this.zoom);
      }

      /**
       * Find best seed, select it, show traceback, mark best chain, zoom to fit.
       */
      #initialize() {
        const seeds = this.data.seeds;
        const bestSeed = seeds.reduce((best, s) => (s.max_score > best.max_score ? s : best), seeds[0]);
        if (bestSeed && bestSeed.max_score > 0) {
          this.selectedSeed = bestSeed;
          this.seedCircles.filter(s => s.index === bestSeed.index).classed('selected', true);
          this.#showTraceback(bestSeed.index);
          const bestChainSeedIndices = new Set(this.tracebackSeedIndices);
          this.seedCircles.classed('on-best-chain', d => bestChainSeedIndices.has(d.index));
          this.#zoomToFit(seeds.filter(s => this.tracebackSeedIndices.has(s.index)), 0.05);
        } else {
          this.#zoomToFit(seeds, 0.05);
        }
      }

      // -----------------------------------------------------------
      // Formatting methods
      // -----------------------------------------------------------

      /**
       * Base tooltip text for a seed.
       */
      #seedTooltipText(d) {
        return `Seed #${d.seed_num} (${d.strand})\\nSeed ID: ${d.seed_id}\\nRead: ${d.read_pos}\\nRef: ${d.ref_pos}\\nMax score: ${d.max_score}`;
      }

      /**
       * Tooltip text for a transition line.
       */
      #transitionTooltipText(d) {
        return `Score: ${d.score}\\nFraction of max: ${(d.fraction * 100).toFixed(1)}%`;
      }

      /**
       * Set x1/y1/x2/y2 on a D3 selection of <line> elements using
       * current scales.
       */
      #positionLines(selection) {
        const viz = this;
        selection.each(function(d) {
          const src = viz.data.seedById(d.source_id);
          const dst = viz.data.seedById(d.dest_id);
          if (src && dst) {
            d3.select(this)
              .attr('x1', viz.currentXScale(src.ref_pos))
              .attr('y1', viz.currentYScale(src.read_pos))
              .attr('x2', viz.currentXScale(dst.ref_pos))
              .attr('y2', viz.currentYScale(dst.read_pos));
          }
        });
      }

      /**
       * Return a closure suitable for D3 selection.call() that positions
       * <line> elements using the current zoom-adjusted scales.
       */
      #positionLinesFn() {
        const viz = this;
        return function(selection) { viz.#positionLines(selection); };
      }

      // -----------------------------------------------------------
      // Interaction methods
      // -----------------------------------------------------------

      /**
       * Follow best-scoring transitions backward from seedIndex.
       */
      #computeTraceback(seedIndex) {
        const path = [];
        const visited = new Set();
        let currentIndex = seedIndex;
        while (currentIndex !== undefined && !visited.has(currentIndex)) {
          visited.add(currentIndex);
          const bestTrans = this.data.bestTransitionTo(currentIndex);
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
       * Display the traceback path from seedIndex as red lines.
       */
      #showTraceback(seedIndex) {
        const path = this.#computeTraceback(seedIndex);
        this.tracebackSeedIndices.clear();
        path.forEach(t => {
          this.tracebackSeedIndices.add(t.source_index);
          this.tracebackSeedIndices.add(t.dest_index);
        });

        const tracebackSet = this.tracebackSeedIndices;
        this.seedCircles.classed('on-traceback', d => tracebackSet.has(d.index) && d.index !== seedIndex);

        const lines = this.tracebackLayer.selectAll('.traceback')
          .data(path, d => d.source_id + '->' + d.dest_id);

        lines.enter()
          .append('line')
          .attr('class', 'transition traceback')
          .call(this.#positionLinesFn())
          .append('title')
          .text(d => this.#transitionTooltipText(d));

        lines.exit().remove();
      }

      #hideTraceback() {
        this.tracebackLayer.selectAll('.traceback').remove();
        this.tracebackSeedIndices.clear();
        this.seedCircles.classed('on-traceback', false);
      }

      /**
       * Show positive-score transitions arriving at seedIndex as colored lines.
       */
      #showSecondaryTransitions(seedIndex) {
        const allTrans = this.data.transitionsTo(seedIndex);
        const positiveTrans = allTrans.filter(t => t.score > 0);
        const className = 'secondary-dest-' + seedIndex;

        this.secondaryTransitionLayer.selectAll('.' + className)
          .data(positiveTrans, d => d.source_id + '->' + d.dest_id)
          .enter()
          .append('line')
          .attr('class', 'transition secondary ' + className)
          .call(this.#positionLinesFn())
          .attr('stroke', d => d.is_max ? 'red' : this.colorScale(d.fraction))
          .append('title')
          .text(d => this.#transitionTooltipText(d));
      }

      #hideSecondaryTransitions(seedIndex) {
        this.secondaryTransitionLayer.selectAll('.secondary-dest-' + seedIndex).remove();
      }

      /**
       * Draw a thick black line from hoveredIndex to the selected seed.
       */
      #highlightTransitionToSelected(hoveredIndex) {
        if (!this.selectedSeed || hoveredIndex === this.selectedSeed.index) return;
        const trans = this.data.transitionFromTo(hoveredIndex, this.selectedSeed.index);
        if (!trans) return;

        this.hoveredToSelectedLayer.append('line')
          .datum(trans)
          .attr('class', 'transition to-selected')
          .call(this.#positionLinesFn())
          .append('title')
          .text(d => this.#transitionTooltipText(d));
      }

      #unhighlightTransitionToSelected() {
        this.hoveredToSelectedLayer.selectAll('*').remove();
      }

      /**
       * Rebuild a seed's tooltip, adding transition info when a seed is selected.
       */
      #updateSeedTooltip(circle, d) {
        let text = this.#seedTooltipText(d);
        if (this.selectedSeed && this.selectedSeed !== d) {
          const trans = this.data.transitionFromTo(d.index, this.selectedSeed.index);
          if (trans) {
            text += `\\nScore to selected: ${trans.score}`;
          }
          const refDiff = Math.abs(d.ref_pos - this.selectedSeed.ref_pos);
          const readDiff = Math.abs(d.read_pos - this.selectedSeed.read_pos);
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

      /**
       * Zoom to fit the given seeds with padding as a fraction of data range.
       */
      #zoomToFit(seedsToFit, padding) {
        if (seedsToFit.length === 0) return;
        const [xMin, xMax] = d3.extent(seedsToFit, s => s.ref_pos);
        const [yMin, yMax] = d3.extent(seedsToFit, s => s.read_pos);
        const xPad = (xMax - xMin) * padding || 100;
        const yPad = (yMax - yMin) * padding || 100;
        const targetXRange = (xMax - xMin) + 2 * xPad;
        const targetYRange = (yMax - yMin) + 2 * yPad;
        const kx = (this.initialXDomain[1] - this.initialXDomain[0]) / targetXRange;
        const ky = (this.initialYDomain[1] - this.initialYDomain[0]) / targetYRange;
        const k = Math.min(kx, ky);
        const cx = (xMin + xMax) / 2;
        const cy = (yMin + yMax) / 2;
        const tx = this._width / 2 - k * this.baseXScale(cx);
        const ty = this._height / 2 - k * this.baseYScale(cy);
        this.svg.call(this.zoom.transform, d3.zoomIdentity.translate(tx, ty).scale(k));
      }
    }

    document.addEventListener('DOMContentLoaded', async function() {
      const data = await ChainData.fromCompressed(SEEDS_COMPRESSED, TRANSITIONS_COMPRESSED);
      new ChainViz(data);
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
