#!/usr/bin/env python3
"""
Build an interactive SVG visualization of a chaining problem using D3.js.

Usage: build-chain-viz.py <data_directory> <chain_number> <output.svg>
"""

import sys
import os
import glob
import json
import struct
from collections import defaultdict


def uint64_to_int64(val):
    """Reinterpret a 64-bit unsigned integer as a signed integer."""
    return struct.unpack('q', struct.pack('Q', val))[0]


def parse_seeds_file(filepath):
    """
    Parse a chain seeds file.
    Returns list of dicts with keys: read_pos, ref_name, ref_pos, seed_num, seed_id
    """
    seeds = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 5:
                continue
            seeds.append({
                'read_pos': int(parts[0]),
                'ref_name': parts[1],
                'ref_pos': int(parts[2]),
                'seed_num': int(parts[3]),
                'seed_id': parts[4]
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
            score_uint64 = int(parts[2])
            score_int64 = uint64_to_int64(score_uint64)
            transitions.append({
                'source_id': parts[0],
                'dest_id': parts[1],
                'score': score_int64
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

    # Build seed lookup by ID
    seed_by_id = {s['seed_id']: s for s in seeds}
    seed_ids = set(seed_by_id.keys())

    # Assign numeric indices to seeds for use in CSS class names
    seed_id_to_index = {s['seed_id']: i for i, s in enumerate(seeds)}

    # Filter transitions to only include those where both endpoints are in our seeds
    # and score is positive
    relevant_transitions = []
    for t in transitions:
        if t['source_id'] in seed_ids and t['dest_id'] in seed_ids and t['score'] > 0:
            relevant_transitions.append(t)

    # Compute max score per destination
    max_score_by_dest = defaultdict(lambda: float('-inf'))
    for t in relevant_transitions:
        if t['score'] > max_score_by_dest[t['dest_id']]:
            max_score_by_dest[t['dest_id']] = t['score']

    # Add score info to seeds
    seeds_data = []
    for i, s in enumerate(seeds):
        max_score = max_score_by_dest.get(s['seed_id'], 0)
        seeds_data.append({
            'index': i,
            'seed_num': s['seed_num'],
            'read_pos': s['read_pos'],
            'ref_pos': s['ref_pos'],
            'seed_id': s['seed_id'],
            'max_score': max_score if max_score > float('-inf') else 0
        })

    # Prepare transitions data with fraction of max
    transitions_data = []
    for t in relevant_transitions:
        max_score = max_score_by_dest[t['dest_id']]
        fraction = t['score'] / max_score if max_score > 0 else 0
        is_max = (t['score'] == max_score)
        dest_index = seed_id_to_index.get(t['dest_id'], -1)
        source_index = seed_id_to_index.get(t['source_id'], -1)
        transitions_data.append({
            'source_id': t['source_id'],
            'source_index': source_index,
            'dest_id': t['dest_id'],
            'dest_index': dest_index,
            'score': t['score'],
            'fraction': fraction,
            'is_max': is_max
        })

    seeds_json = json.dumps(seeds_data)
    transitions_json = json.dumps(transitions_data)

    svg_content = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="1200" height="800">
  <style>
    .seed {{
      fill: steelblue;
      stroke: #333;
      stroke-width: 1px;
      cursor: pointer;
    }}
    .seed.on-best-chain {{
      fill: #90EE90;
      stroke: darkgreen;
      stroke-width: 2px;
    }}
    .seed:hover, .seed.hovered, .seed.on-best-chain:hover, .seed.on-best-chain.hovered {{
      r: 8;
    }}
    .seed.selected {{
      fill: orange;
      r: 8;
    }}
    .seed.on-traceback {{
      fill: #ff6666;
      stroke: darkred;
      stroke-width: 2px;
    }}
    .transition {{
      stroke-width: 1.5px;
      fill: none;
      opacity: 0.6;
    }}
    .transition.highlighted {{
      stroke-width: 3px;
      opacity: 1;
    }}
    .transition.to-selected {{
      stroke: black !important;
      stroke-width: 4px !important;
      opacity: 1 !important;
    }}
    .transition.secondary {{
      stroke-width: 1px;
      opacity: 0.4;
    }}
    .transition.traceback {{
      stroke: red;
      stroke-width: 3px;
      opacity: 1;
    }}
    .axis text {{
      font-family: sans-serif;
      font-size: 12px;
    }}
    .axis path, .axis line {{
      stroke: #333;
    }}
    .axis-label {{
      font-family: sans-serif;
      font-size: 14px;
      font-weight: bold;
    }}
  </style>
  <script xlink:href="https://d3js.org/d3.v7.js"></script>
  <script type="text/javascript">
    <![CDATA[
    document.addEventListener('DOMContentLoaded', function() {{
      const seeds = {seeds_json};
      const transitions = {transitions_json};

      const margin = {{top: 40, right: 40, bottom: 60, left: 80}};
      const svgWidth = 1200;
      const svgHeight = 800;
      const width = svgWidth - margin.left - margin.right;
      const height = svgHeight - margin.top - margin.bottom;

      const svg = d3.select('svg');

      // Build lookup
      const seedById = new Map();
      seeds.forEach(s => seedById.set(s.seed_id, s));

      // Group transitions by destination index for efficient lookup
      const transitionsByDestIndex = new Map();
      transitions.forEach(t => {{
        if (!transitionsByDestIndex.has(t.dest_index)) {{
          transitionsByDestIndex.set(t.dest_index, []);
        }}
        transitionsByDestIndex.get(t.dest_index).push(t);
      }});

      // Group transitions by source index for efficient lookup
      const transitionsBySourceIndex = new Map();
      transitions.forEach(t => {{
        if (!transitionsBySourceIndex.has(t.source_index)) {{
          transitionsBySourceIndex.set(t.source_index, []);
        }}
        transitionsBySourceIndex.get(t.source_index).push(t);
      }});

      // Find best (max score) transition to each seed
      const bestTransitionToSeed = new Map();
      transitions.forEach(t => {{
        if (t.is_max) {{
          bestTransitionToSeed.set(t.dest_index, t);
        }}
      }});

      // Get transition from source to dest (if exists)
      function getTransitionFromTo(sourceIndex, destIndex) {{
        const fromSource = transitionsBySourceIndex.get(sourceIndex) || [];
        return fromSource.find(t => t.dest_index === destIndex);
      }}

      // Compute data extents (handle empty data)
      let xExtent = d3.extent(seeds, d => d.ref_pos);
      let yExtent = d3.extent(seeds, d => d.read_pos);

      // Handle empty data case
      if (xExtent[0] === undefined) xExtent = [0, 1000];
      if (yExtent[0] === undefined) yExtent = [0, 1000];

      // Compute ranges
      const xRange = xExtent[1] - xExtent[0] || 1000;
      const yRange = yExtent[1] - yExtent[0] || 1000;

      // For 1:1 aspect ratio, we need pixels per data unit to be equal for both axes
      // Compute the scale factor that fits both dimensions
      const xPixelsPerUnit = width / xRange;
      const yPixelsPerUnit = height / yRange;
      const pixelsPerUnit = Math.min(xPixelsPerUnit, yPixelsPerUnit);

      // Compute centered domains with 1:1 aspect ratio
      const xCenter = (xExtent[0] + xExtent[1]) / 2;
      const yCenter = (yExtent[0] + yExtent[1]) / 2;
      const xHalfRange = (width / pixelsPerUnit) / 2;
      const yHalfRange = (height / pixelsPerUnit) / 2;

      // Add 5% padding
      const xPadding = xHalfRange * 0.05;
      const yPadding = yHalfRange * 0.05;

      // Initial scale (store for zoom reset)
      const initialXDomain = [xCenter - xHalfRange - xPadding, xCenter + xHalfRange + xPadding];
      const initialYDomain = [yCenter - yHalfRange - yPadding, yCenter + yHalfRange + yPadding];

      // Scales - note: both use the same pixels-per-unit ratio
      const xScale = d3.scaleLinear()
        .domain(initialXDomain)
        .range([0, width]);

      const yScale = d3.scaleLinear()
        .domain(initialYDomain)
        .range([height, 0]);  // Inverted so higher read pos is at top

      // Color scale for transitions (fraction of max score)
      const colorScale = d3.scaleSequential(d3.interpolateViridis)
        .domain([0, 1]);

      // Create main group with clipping
      const clipId = 'clip-' + Math.random().toString(36).substr(2, 9);
      svg.append('defs').append('clipPath')
        .attr('id', clipId)
        .append('rect')
        .attr('width', width)
        .attr('height', height);

      const g = svg.append('g')
        .attr('transform', `translate(${{margin.left}},${{margin.top}})`);

      // Axes
      const xAxis = d3.axisBottom(xScale).ticks(10);
      const yAxis = d3.axisLeft(yScale).ticks(10);

      const xAxisG = g.append('g')
        .attr('class', 'axis x-axis')
        .attr('transform', `translate(0,${{height}})`)
        .call(xAxis);

      const yAxisG = g.append('g')
        .attr('class', 'axis y-axis')
        .call(yAxis);

      // Axis labels
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

      // Create container for zoomable content
      const plotArea = g.append('g')
        .attr('clip-path', `url(#${{clipId}})`);

      // Add a background rect for capturing zoom events
      plotArea.append('rect')
        .attr('width', width)
        .attr('height', height)
        .attr('fill', 'white');

      const zoomG = plotArea.append('g');

      // Layer for traceback transitions (red path)
      const tracebackLayer = zoomG.append('g').attr('class', 'traceback-layer');
      // Layer for dynamically added secondary transitions (on hover)
      const secondaryTransitionLayer = zoomG.append('g').attr('class', 'secondary-transition-layer');
      // Layer for seeds (on top)
      const seedLayer = zoomG.append('g').attr('class', 'seed-layer');

      // Current scale references (updated by zoom)
      let currentXScale = xScale;
      let currentYScale = yScale;

      // Helper to compute line coordinates
      function getLineCoords(d, xS, yS) {{
        const src = seedById.get(d.source_id);
        const dest = seedById.get(d.dest_id);
        return {{
          x1: src ? xS(src.ref_pos) : 0,
          y1: src ? yS(src.read_pos) : 0,
          x2: dest ? xS(dest.ref_pos) : 0,
          y2: dest ? yS(dest.read_pos) : 0
        }};
      }}

      // Track state
      let selectedSeed = null;
      let hoveredSeed = null;
      let tracebackSeedIndices = new Set();

      // Compute traceback path from a seed
      function computeTraceback(seedIndex) {{
        const path = [];
        const visited = new Set();
        let currentIndex = seedIndex;

        while (currentIndex !== undefined && !visited.has(currentIndex)) {{
          visited.add(currentIndex);
          const bestTrans = bestTransitionToSeed.get(currentIndex);
          if (bestTrans) {{
            path.push(bestTrans);
            currentIndex = bestTrans.source_index;
          }} else {{
            break;
          }}
        }}
        return path;
      }}

      // Show traceback path for selected seed
      function showTraceback(seedIndex) {{
        const path = computeTraceback(seedIndex);
        tracebackSeedIndices.clear();

        // Collect all seed indices on the traceback
        if (path.length > 0) {{
          path.forEach(t => {{
            tracebackSeedIndices.add(t.source_index);
            tracebackSeedIndices.add(t.dest_index);
          }});
        }}

        // Update seed styling for traceback
        seedCircles.classed('on-traceback', d => tracebackSeedIndices.has(d.index) && d.index !== seedIndex);

        // Draw traceback transitions
        const lines = tracebackLayer.selectAll('.traceback')
          .data(path, d => d.source_id + '->' + d.dest_id);

        lines.enter()
          .append('line')
          .attr('class', 'transition traceback')
          .each(function(d) {{
            const coords = getLineCoords(d, currentXScale, currentYScale);
            d3.select(this)
              .attr('x1', coords.x1)
              .attr('y1', coords.y1)
              .attr('x2', coords.x2)
              .attr('y2', coords.y2);
          }})
          .append('title')
          .text(d => `Score: ${{d.score}}\\nFraction of max: ${{(d.fraction * 100).toFixed(1)}}%`);

        lines.exit().remove();
      }}

      // Hide traceback path
      function hideTraceback() {{
        tracebackLayer.selectAll('.traceback').remove();
        tracebackSeedIndices.clear();
        seedCircles.classed('on-traceback', false);
      }}

      // Function to show all transitions to a seed (by index) on hover
      function showSecondaryTransitions(seedIndex) {{
        const allTrans = transitionsByDestIndex.get(seedIndex) || [];

        // Use numeric index for reliable class selection
        const className = 'secondary-dest-' + seedIndex;

        const lines = secondaryTransitionLayer.selectAll('.' + className)
          .data(allTrans, d => d.source_id + '->' + d.dest_id);

        lines.enter()
          .append('line')
          .attr('class', 'transition secondary ' + className)
          .each(function(d) {{
            const coords = getLineCoords(d, currentXScale, currentYScale);
            d3.select(this)
              .attr('x1', coords.x1)
              .attr('y1', coords.y1)
              .attr('x2', coords.x2)
              .attr('y2', coords.y2);
          }})
          .attr('stroke', d => d.is_max ? 'red' : colorScale(d.fraction))
          .append('title')
          .text(d => `Score: ${{d.score}}\\nFraction of max: ${{(d.fraction * 100).toFixed(1)}}%`);
      }}

      // Function to hide secondary transitions to a seed (by index)
      function hideSecondaryTransitions(seedIndex) {{
        const className = 'secondary-dest-' + seedIndex;
        secondaryTransitionLayer.selectAll('.' + className).remove();
      }}

      // Function to update seed tooltip based on current selection
      function updateSeedTooltip(circle, d) {{
        let tooltipText = `Seed #${{d.seed_num}}\\nSeed ID: ${{d.seed_id}}\\nRead: ${{d.read_pos}}\\nRef: ${{d.ref_pos}}\\nMax score: ${{d.max_score}}`;
        if (selectedSeed && selectedSeed !== d) {{
          const trans = getTransitionFromTo(d.index, selectedSeed.index);
          if (trans) {{
            tooltipText += `\\nScore to selected: ${{trans.score}}`;
          }}
          // Compute insertion/deletion offset
          const refDiff = Math.abs(d.ref_pos - selectedSeed.ref_pos);
          const readDiff = Math.abs(d.read_pos - selectedSeed.read_pos);
          const offset = Math.abs(refDiff - readDiff);
          let offsetStr;
          if (readDiff > refDiff) {{
            offsetStr = `${{offset}}bp INS`;
          }} else if (refDiff > readDiff) {{
            offsetStr = `${{offset}}bp DEL`;
          }} else {{
            offsetStr = `0bp`;
          }}
          tooltipText += `\\nOffset: ${{offsetStr}}`;
        }}
        circle.select('title').text(tooltipText);
      }}

      // Function to highlight transition from hovered to selected
      function highlightTransitionToSelected(hoveredIndex) {{
        if (!selectedSeed || hoveredIndex === selectedSeed.index) return;
        // Find and bold the transition line from hovered to selected in secondary layer
        secondaryTransitionLayer.selectAll('.secondary-dest-' + selectedSeed.index)
          .filter(t => t.source_index === hoveredIndex)
          .classed('to-selected', true);
      }}

      // Function to unhighlight transition from hovered to selected
      function unhighlightTransitionToSelected(hoveredIndex) {{
        if (!selectedSeed) return;
        secondaryTransitionLayer.selectAll('.secondary-dest-' + selectedSeed.index)
          .filter(t => t.source_index === hoveredIndex)
          .classed('to-selected', false);
      }}

      // Draw seeds with SVG title tooltips
      const seedCircles = seedLayer.selectAll('.seed')
        .data(seeds)
        .enter()
        .append('circle')
        .attr('class', 'seed')
        .attr('cx', d => xScale(d.ref_pos))
        .attr('cy', d => yScale(d.read_pos))
        .attr('r', 5)
        .on('mouseover', function(event, d) {{
          hoveredSeed = d;
          if (selectedSeed !== d) {{
            d3.select(this).classed('hovered', true);
          }}
          // Update tooltip to show transition score to selected if applicable
          updateSeedTooltip(d3.select(this), d);
          // Show all transitions to this seed
          showSecondaryTransitions(d.index);
          // Highlight transition from this seed to selected seed
          highlightTransitionToSelected(d.index);
        }})
        .on('mouseout', function(event, d) {{
          hoveredSeed = null;
          if (selectedSeed !== d) {{
            d3.select(this).classed('hovered', false);
            // Hide secondary transitions
            hideSecondaryTransitions(d.index);
          }}
          // Unhighlight transition to selected
          unhighlightTransitionToSelected(d.index);
        }})
        .on('click', function(event, d) {{
          const prevSelected = selectedSeed;

          // If clicking same seed, deselect it
          if (selectedSeed === d) {{
            selectedSeed = null;
            d3.select(this).classed('selected', false);
            hideTraceback();
            // If not still hovering, hide transitions
            if (hoveredSeed !== d) {{
              hideSecondaryTransitions(d.index);
            }}
          }} else {{
            // Deselect previous if any
            if (prevSelected) {{
              seedCircles.filter(s => s.index === prevSelected.index)
                .classed('selected', false)
                .classed('hovered', false);
              // Hide previous seed's transitions if not hovering it
              if (hoveredSeed === null || hoveredSeed.index !== prevSelected.index) {{
                hideSecondaryTransitions(prevSelected.index);
              }}
            }}

            // Select new seed
            selectedSeed = d;
            d3.select(this).classed('selected', true);
            showTraceback(d.index);
            showSecondaryTransitions(d.index);
          }}
        }});

      // Add SVG title tooltips to seeds (initial)
      seedCircles.append('title')
        .text(d => `Seed #${{d.seed_num}}\\nSeed ID: ${{d.seed_id}}\\nRead: ${{d.read_pos}}\\nRef: ${{d.ref_pos}}\\nMax score: ${{d.max_score}}`);

      // Zoom behavior with 1:1 aspect ratio constraint
      const zoom = d3.zoom()
        .scaleExtent([0.1, 10000])
        .on('zoom', function(event) {{
          const transform = event.transform;

          // Update scales (both use same transform to maintain 1:1 ratio)
          currentXScale = transform.rescaleX(xScale);
          currentYScale = transform.rescaleY(yScale);

          // Update axes
          xAxisG.call(xAxis.scale(currentXScale));
          yAxisG.call(yAxis.scale(currentYScale));

          // Update seed positions
          seedCircles
            .attr('cx', d => currentXScale(d.ref_pos))
            .attr('cy', d => currentYScale(d.read_pos));

          // Update traceback transition lines
          tracebackLayer.selectAll('.transition').each(function(d) {{
            const coords = getLineCoords(d, currentXScale, currentYScale);
            d3.select(this)
              .attr('x1', coords.x1)
              .attr('y1', coords.y1)
              .attr('x2', coords.x2)
              .attr('y2', coords.y2);
          }});

          // Update secondary transition lines
          secondaryTransitionLayer.selectAll('.transition').each(function(d) {{
            const coords = getLineCoords(d, currentXScale, currentYScale);
            d3.select(this)
              .attr('x1', coords.x1)
              .attr('y1', coords.y1)
              .attr('x2', coords.x2)
              .attr('y2', coords.y2);
          }});
        }});

      svg.call(zoom);

      // Function to zoom to fit a set of seed indices
      function zoomToFitSeeds(seedIndices) {{
        if (seedIndices.size === 0) return;

        const seedsToFit = seeds.filter(s => seedIndices.has(s.index));
        if (seedsToFit.length === 0) return;

        const xMin = d3.min(seedsToFit, s => s.ref_pos);
        const xMax = d3.max(seedsToFit, s => s.ref_pos);
        const yMin = d3.min(seedsToFit, s => s.read_pos);
        const yMax = d3.max(seedsToFit, s => s.read_pos);

        // Add padding (50% on each side for a more zoomed-out view)
        const xPad = (xMax - xMin) * 0.5 || 100;
        const yPad = (yMax - yMin) * 0.5 || 100;

        const targetXMin = xMin - xPad;
        const targetXMax = xMax + xPad;
        const targetYMin = yMin - yPad;
        const targetYMax = yMax + yPad;

        const targetXRange = targetXMax - targetXMin;
        const targetYRange = targetYMax - targetYMin;

        // Compute scale to fit while maintaining 1:1 aspect ratio
        const scaleX = width / targetXRange;
        const scaleY = height / targetYRange;
        const scale = Math.min(scaleX, scaleY);

        // Compute the center of the target region
        const targetXCenter = (targetXMin + targetXMax) / 2;
        const targetYCenter = (targetYMin + targetYMax) / 2;

        // Compute translation to center the target region
        // The base scales map initialXDomain/initialYDomain to [0,width]/[height,0]
        // We need to find the transform that maps our target to fill the view

        // In the base scale, what pixel coords would the target center be at?
        const baseCenterX = xScale(targetXCenter);
        const baseCenterY = yScale(targetYCenter);

        // We want the target center to be at the view center (width/2, height/2)
        // After applying transform: newX = transform.k * baseX + transform.x
        // We want: transform.k * baseCenterX + transform.x = width/2
        // And: transform.k * baseCenterY + transform.y = height/2

        // First compute the scale factor relative to base scale
        const baseScaleXRange = initialXDomain[1] - initialXDomain[0];
        const k = baseScaleXRange / targetXRange * (width / width);  // simplified: scale relative to base

        // Actually, let's compute it more directly
        // The base xScale maps initialXDomain to [0, width]
        // We want to map targetX range to [0, width] with scale k
        // k = (initialXDomain range) / (targetX range) when both map to same pixel range
        const k2 = (initialXDomain[1] - initialXDomain[0]) / targetXRange;

        // Translation: we want targetXCenter to map to width/2
        // xScale(targetXCenter) gives us the pixel position in base coords
        // After transform: k * xScale(targetXCenter) + tx = width/2
        const tx = width / 2 - k2 * xScale(targetXCenter);
        const ty = height / 2 - k2 * yScale(targetYCenter);

        const transform = d3.zoomIdentity.translate(tx, ty).scale(k2);
        svg.call(zoom.transform, transform);
      }}

      // Auto-select the seed with highest max_score and zoom to traceback
      const bestSeed = seeds.reduce((best, s) => (s.max_score > best.max_score ? s : best), seeds[0]);
      const bestChainSeedIndices = new Set();
      if (bestSeed && bestSeed.max_score > 0) {{
        selectedSeed = bestSeed;
        seedCircles.filter(s => s.index === bestSeed.index).classed('selected', true);
        showTraceback(bestSeed.index);

        // Store the best chain seed indices permanently
        tracebackSeedIndices.forEach(idx => bestChainSeedIndices.add(idx));

        // Mark seeds on the best chain with persistent styling
        seedCircles.classed('on-best-chain', d => bestChainSeedIndices.has(d.index));

        // Zoom to fit the traceback
        if (tracebackSeedIndices.size > 0) {{
          zoomToFitSeeds(tracebackSeedIndices);
        }}
      }}
    }});
    ]]>
  </script>
</svg>
'''

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
