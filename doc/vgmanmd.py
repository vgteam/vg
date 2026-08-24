#!/usr/bin/env python3
import subprocess


# commands to include
######### If you add to this, also add it to the intro section of vgmandmd.desc.md
cmds = ['align', 'annotate', 'augment', 'autoindex', 'bench-dist-query',
        'benchmark', 'call', 'chains', 'chunk', 'circularize', 'clip',
        'cluster', 'combine', 'construct', 'convert', 'deconstruct', 'depth',
        'describe', 'dotplot', 'filter', 'find', 'gamcompare', 'gampcompare',
        'gamsort', 'gbwt', 'genotype', 'giraffe', 'haplotypes', 'help', 'ids',
        'index', 'inject', 'map', 'mask', 'minimizer', 'mod', 'mpmap',
        'pack', 'paths', 'primers', 'prune', 'rna', 'sim', 'simplify', 'snarls',
        'stats', 'surject', 'test', 'trace', 'validate', 'vectorize', 'version',
        'view', 'viz', 'zipcode'
]
cmds.sort()


# parse short descriptions
try:
    desc_inf = open('./doc/vgmanmd.desc.md', 'rt')
except FileNotFoundError:
    desc_inf = open('vgmanmd.desc.md', 'rt')

desc = {}
cur_desc = ''
cur_header = ''
for line in desc_inf:
    if line[0] == '#':
        # new header
        if cur_header != '':
            desc[cur_header] = cur_desc
        cur_header = line.rstrip().replace('# ', '')
        cur_desc = ''
    else:
        cur_desc += line
desc[cur_header] = cur_desc
desc_inf.close()

# start page
#print('# vg manpage')

# get vg version
ret = subprocess.run(['vg', 'version'], capture_output=True)
vg_v = ret.stdout.decode().split('\n')[0]

#Metadata
print("% vg(1)  | Variation Graph Toolkit\n\n")

#Start with the name
print("NAME")
print("====")
print('vg - variation graph tool, ' + vg_v + '.\n\n') 

if 'description' in desc:
    print("DESCRIPTION")
    print("====")
    print(desc['description'])
    print('\n\n')

if 'synopsis' in desc:
    print("SYNOPSIS")
    print("====")
    print(desc['synopsis'])
    print('\n\n')

print("COMMANDS")
print("====")

# help for each cmd
vg_help = subprocess.run(['vg', 'help'], capture_output=True)
cmd_desc = dict()
for line in vg_help.stderr.decode().split('\n'):
    if '-' in line:
        parts = line.split()
        cmd_desc[parts[1]] = ' '.join(parts[2:])
        cmd = parts[1]

for cmd in cmds:
    print('<a id="{cmd}"/>\n\n## {cmd}: {blurb}\n\n'.format(cmd=cmd, blurb=cmd_desc[cmd]))
    ret = subprocess.run(['vg', cmd, '--help'], capture_output=True)
    print('```')
    print(ret.stderr.decode())
    print('```\n\n')

if 'bugs' in desc:
    print("BUGS")
    print("====")
    print(desc['bugs'])
    print('\n\n')
