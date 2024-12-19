#!/usr/bin/env python3
import subprocess


# commands to include
######### If you add to this, also add it to the intro section of vgmandmd.desc.md
cmds = ['index', 'view', 'autoindex', 'pack', 'giraffe', 'map', 'call',
        'mpmap', 'rna', 'chunk', 'stats', 'gbwt', 'paths', 'find', 'filter',
        'construct', 'minimizer', 'haplotypes', 'deconstruct', 'convert',
        'gamsort', 'inject', 'surject', 'mod', 'prune', 'ids', 'sim', 'annotate']
cmds.sort()


# parse short descriptions
desc_inf = open('./doc/vgmanmd.desc.md', 'rt')
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

# table of contents
#for cmd in cmds:
#    print(' - [vg {cmd}](#{cmd})'.format(cmd=cmd))
#
#print('\n\n')

print("COMMANDS")
print("====")

# help for each cmd
for cmd in cmds:
    print('## {cmd}\n\n'.format(cmd=cmd))
    # run subcommand without -h because not everything has -h
    ret = subprocess.run(['vg', cmd], capture_output=True)
    print('```')
    if cmd in desc:
        print(desc[cmd])
        print('\n\n')
    print(ret.stderr.decode())
    print('```\n\n')

if 'bugs' in desc:
    print("BUGS")
    print("====")
    print(desc['bugs'])
    print('\n\n')
