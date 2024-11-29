import subprocess


# commands to include
######### If you add to this, also add it to the intro section of vgmandmd.desc.md
cmds = ['index', 'view', 'autoindex', 'pack', 'giraffe', 'map', 'call',
        'mpmap', 'rna', 'chunk', 'stats', 'gbwt', 'paths', 'find', 'filter',
        'construct', 'minimizer', 'haplotypes', 'deconstruct', 'convert',
        'gamsort', 'inject', 'surject', 'mod', 'prune', 'ids', 'sim', 'annotate']
cmds.sort()

# parse short descriptions
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

print('\n*Automatically made for ' + vg_v + '.*\n\n')

# add intro text
if 'intro' in desc:
    print(desc['intro'])
    print('\n\n')

# table of contents
#for cmd in cmds:
#    print(' - [vg {cmd}](#{cmd})'.format(cmd=cmd))
#
#print('\n\n')

# help for each cmd
for cmd in cmds:
    print('## {cmd}\n\n'.format(cmd=cmd))
    if cmd in desc:
        print(desc[cmd])
        print('\n\n')
    # run subcommand with -h
    ret = subprocess.run(['vg', cmd, '-h'], capture_output=True)
    print('```')
    print(ret.stderr.decode())
    print('```\n\n')
