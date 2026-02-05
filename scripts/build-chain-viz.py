#!/usr/bin/env python3
"""
Build an interactive SVG visualization of a chaining problem using D3.js.

Usage: build-chain-viz.py <data_directory> <chain_number> <output.svg>
"""

import sys
import os
import glob
import re
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
    for s in seeds:
        max_score = max_score_by_dest.get(s['seed_id'], 0)
        seeds_data.append({
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
        transitions_data.append({
            'source_id': t['source_id'],
            'dest_id': t['dest_id'],
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
    .seed:hover, .seed.hovered {{
      r: 8;
    }}
    .seed.selected {{
      fill: orange;
      r: 8;
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
    .axis text {{
      font-family: sans-serif;
      font-size: 12px;
    }}
    .axis path, .axis line {{
      stroke: #333;
    }}
    .tooltip {{
      position: absolute;
      background: rgba(0, 0, 0, 0.8);
      color: white;
      padding: 8px;
      border-radius: 4px;
      font-family: monospace;
      font-size: 11px;
      pointer-events: none;
      max-width: 400px;
      word-wrap: break-word;
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
      const width = 1200 - margin.left - margin.right;
      const height = 800 - margin.top - margin.bottom;

      const svg = d3.select('svg');

      // Create tooltip div
      const tooltip = d3.select('body').append('div')
        .attr('class', 'tooltip')
        .style('opacity', 0)
        .style('position', 'absolute');

      // Build lookup
      const seedById = new Map();
      seeds.forEach(s => seedById.set(s.seed_id, s));

      // Compute data extents (handle empty data)
      let xExtent = d3.extent(seeds, d => d.ref_pos);
      let yExtent = d3.extent(seeds, d => d.read_pos);

      // Handle empty data case
      if (xExtent[0] === undefined) xExtent = [0, 1000];
      if (yExtent[0] === undefined) yExtent = [0, 1000];

      // Add some padding
      const xPadding = (xExtent[1] - xExtent[0]) * 0.05 || 100;
      const yPadding = (yExtent[1] - yExtent[0]) * 0.05 || 100;

      // Scales
      const xScale = d3.scaleLinear()
        .domain([xExtent[0] - xPadding, xExtent[1] + xPadding])
        .range([0, width]);

      const yScale = d3.scaleLinear()
        .domain([yExtent[0] - yPadding, yExtent[1] + yPadding])
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

      const zoomG = plotArea.append('g');

      // Draw transitions first (so seeds are on top)
      const transitionLines = zoomG.selectAll('.transition')
        .data(transitions)
        .enter()
        .append('line')
        .attr('class', 'transition')
        .attr('x1', d => {{
          const src = seedById.get(d.source_id);
          return src ? xScale(src.ref_pos) : 0;
        }})
        .attr('y1', d => {{
          const src = seedById.get(d.source_id);
          return src ? yScale(src.read_pos) : 0;
        }})
        .attr('x2', d => {{
          const dest = seedById.get(d.dest_id);
          return dest ? xScale(dest.ref_pos) : 0;
        }})
        .attr('y2', d => {{
          const dest = seedById.get(d.dest_id);
          return dest ? yScale(dest.read_pos) : 0;
        }})
        .attr('stroke', d => d.is_max ? 'red' : colorScale(d.fraction))
        .on('mouseover', function(event, d) {{
          d3.select(this).classed('highlighted', true);
          tooltip.transition().duration(200).style('opacity', 0.9);
          tooltip.html(`Score: ${{d.score}}<br>Fraction of max: ${{(d.fraction * 100).toFixed(1)}}%`)
            .style('left', (event.pageX + 10) + 'px')
            .style('top', (event.pageY - 28) + 'px');
        }})
        .on('mouseout', function() {{
          d3.select(this).classed('highlighted', false);
          tooltip.transition().duration(500).style('opacity', 0);
        }});

      // Draw seeds
      let selectedSeed = null;

      const seedCircles = zoomG.selectAll('.seed')
        .data(seeds)
        .enter()
        .append('circle')
        .attr('class', 'seed')
        .attr('cx', d => xScale(d.ref_pos))
        .attr('cy', d => yScale(d.read_pos))
        .attr('r', 5)
        .on('mouseover', function(event, d) {{
          if (selectedSeed !== d) {{
            d3.select(this).classed('hovered', true);
          }}
          // Highlight incoming transitions
          transitionLines.filter(t => t.dest_id === d.seed_id)
            .classed('highlighted', true);
          tooltip.transition().duration(200).style('opacity', 0.9);
          tooltip.html(`Seed ID: ${{d.seed_id}}<br>Read: ${{d.read_pos}}<br>Ref: ${{d.ref_pos}}<br>Max score: ${{d.max_score}}`)
            .style('left', (event.pageX + 10) + 'px')
            .style('top', (event.pageY - 28) + 'px');
        }})
        .on('mouseout', function(event, d) {{
          if (selectedSeed !== d) {{
            d3.select(this).classed('hovered', false);
          }}
          if (selectedSeed === null || selectedSeed.seed_id !== d.seed_id) {{
            transitionLines.filter(t => t.dest_id === d.seed_id)
              .classed('highlighted', false);
          }}
          tooltip.transition().duration(500).style('opacity', 0);
        }})
        .on('click', function(event, d) {{
          // Deselect previous
          if (selectedSeed) {{
            seedCircles.filter(s => s.seed_id === selectedSeed.seed_id)
              .classed('selected', false)
              .classed('hovered', false);
            transitionLines.filter(t => t.dest_id === selectedSeed.seed_id)
              .classed('highlighted', false);
          }}

          if (selectedSeed === d) {{
            // Clicking same seed deselects
            selectedSeed = null;
          }} else {{
            // Select new seed
            selectedSeed = d;
            d3.select(this).classed('selected', true);
            transitionLines.filter(t => t.dest_id === d.seed_id)
              .classed('highlighted', true);
          }}
        }});

      // Zoom behavior
      const zoom = d3.zoom()
        .scaleExtent([0.1, 50])
        .on('zoom', function(event) {{
          const transform = event.transform;

          // Update scales
          const newXScale = transform.rescaleX(xScale);
          const newYScale = transform.rescaleY(yScale);

          // Update axes
          xAxisG.call(xAxis.scale(newXScale));
          yAxisG.call(yAxis.scale(newYScale));

          // Update seed positions
          seedCircles
            .attr('cx', d => newXScale(d.ref_pos))
            .attr('cy', d => newYScale(d.read_pos));

          // Update transition lines
          transitionLines
            .attr('x1', d => {{
              const src = seedById.get(d.source_id);
              return src ? newXScale(src.ref_pos) : 0;
            }})
            .attr('y1', d => {{
              const src = seedById.get(d.source_id);
              return src ? newYScale(src.read_pos) : 0;
            }})
            .attr('x2', d => {{
              const dest = seedById.get(d.dest_id);
              return dest ? newXScale(dest.ref_pos) : 0;
            }})
            .attr('y2', d => {{
              const dest = seedById.get(d.dest_id);
              return dest ? newYScale(dest.read_pos) : 0;
            }});
        }});

      svg.call(zoom);

      // Add a background rect for capturing zoom events
      plotArea.insert('rect', ':first-child')
        .attr('width', width)
        .attr('height', height)
        .attr('fill', 'white');
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
