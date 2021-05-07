"""
Note that this only works because I don't have a working max_handle_size involved. Because
if I did, I wouldn't be able to assume that snarls in before_norm always share source/sink
with some snarl(s) in after_norm 

TODO: I could avoid double counting if I always skip the source + sink sequence when counting sequence in snarls (in 0_analyze_snarls.py). That would at least keep size counting to a standard.
"""

#%%
import pandas as pd
import numpy as np
import collections as col
import matplotlib.pyplot as plt

# load sizes
names = ["source", "sink", "size"]
## the following snarl size calculations include the sources and sinks of each snarl. This 
## allows significant double-counting of sequence in handles that are shared in the
## source/sink of two adjacent snarls.
before_norm = pd.read_csv("~/paten_lab/vg/robin-graphs/yeast_subset/yeast_subset.snarl_new_sizes.txt", sep="\t", names=names) #using updated snarl calculation, since vg snarls output has changed in past year.
after_norm = pd.read_csv("~/paten_lab/vg/robin-graphs/yeast_subset/yeast_subset.normalized.snarl_sizes.txt", sep="\t", names=names)
## the following snarl size calculations ignore sequence in source and sink:
# before_norm = pd.read_csv("~/paten_lab/vg/robin-graphs/yeast_subset/yeast_subset.snarls_new.skip_source_sink_seq.snarl_sizes.txt", sep="\t", names=names) #using updated snarl calculation, since vg snarls output has changed in past year.
# after_norm = pd.read_csv("~/paten_lab/vg/robin-graphs/yeast_subset/yeast_subset.normalized.skip_source_sink_seq.snarl_sizes.txt", sep="\t", names=names)


#%%
"""
Simple stats: how did the mean and mode snarl sizes change?
As a reality check, do I get the same overall change in sequence size for the snarls?

Interesting problem: before_norm/after_norm total snarl size calculated here doesn't match
the snarl sizes according to normalize_snarls:

according to normalize_snarls:
amount of sequence in normalized snarls before normalization: 880131
amount of sequence in normalized snarls after normalization: 879739
(difference: -392 bases)

Simple Stats output:
before_norm_num_snarls 54822
after_norm_num_snarls 58350
before_norm_mean 278.5043230819744
after_norm_mean 261.632853470437
before_norm_mode 9
after_norm_mode 9
before_norm_range (3, 14380261)
after_norm_range (3, 14380261)
before_norm_total_snarl_size 15268164
after_norm_total_snarl_size 15266277
(difference: -1887 bases)

## NOTE: {the following turned out to not be true:} In both cases, we're double counting a significant quantity of sequence. Any two 
##     snarls sharing source/sink will double count the sequence in that shared handle.

Total sequence in before_norm graph:
14926050 
Total sequence in after_norm graph:
14925658
(difference: -392 bases)

NOTE: Huh. previous note is wrong. No double counting inside normalize_snarls. 

NOTE: maybe the difference is because of snarls that were never normalized, + double 
    counting somehow. There should be more double counting in norm graph, since we get 
    multiple snarl "pinches". 
    There were also two cyclic snarls skipped over in the normalization, and
    so aren't included on the count. They're really big.

    Some problems remain, however:
    It is weird how there's no double counting in the normalize_snarls proper. How?

    And there's slightly more sequence in before_norm_total_snarl_size than actually 
    exists in the total graph. which suggests that counts are substantially off for in-snarl
    counting.

    For now, I think I should proceed with my data analysis. I'm short on time, and the
    few tests I checked in on seemed accurate. Analysis should still give interesting data
    on what snarls are interesting.
    
     Maybe once I find specifically interesting snarls, I'll have better ammo for debugging
    the counts. 
"""
before_norm_num_snarls = len(before_norm["size"])
print("before_norm_num_snarls", before_norm_num_snarls)
after_norm_num_snarls = len(after_norm["size"])
print("after_norm_num_snarls", after_norm_num_snarls)

before_norm_mean = before_norm["size"][1:].mean()
print("before_norm_mean", before_norm_mean)
after_norm_mean = after_norm["size"][1:].mean()
print("after_norm_mean", after_norm_mean)

before_norm_mode = int(before_norm["size"].mode())
print("before_norm_mode", before_norm_mode)
after_norm_mode = int(after_norm["size"].mode())
print("after_norm_mode", after_norm_mode)

before_norm_range = (before_norm["size"].min(), int(before_norm["size"].max()))
print("before_norm_range", before_norm_range)
after_norm_range = (after_norm["size"].min(), int(after_norm["size"].max()))
print("after_norm_range", after_norm_range)

before_norm_total_snarl_size = before_norm["size"].sum()
print("before_norm_total_snarl_size", before_norm_total_snarl_size)
after_norm_total_snarl_size = after_norm["size"].sum()
print("after_norm_total_snarl_size", after_norm_total_snarl_size)
print("(difference:", after_norm_total_snarl_size - before_norm_total_snarl_size, "bases)")
#%%
print("hi")
#%%
"""
Fairly advanced vis: how did specific snarls change? This is for plots for identifying
interesting snarls for case studies, and for visualizing which general regions of graphs
had significant changes.

Two main plots:
* show the size change of each snarl with position along x and size change along y.
* show the snarls that are split into multiple pieces, with position along x and size change along y.

"""

snarl_size_before_after = col.OrderedDict() #key: source of original snarl, value: (before, after) in size (negative if reduced in size)
snarls_that_split = col.defaultdict(int) #key: source of original snarl, value: number of snarls it is split into.

after_i = 0
for before_i in range(len(before_norm["source"])):
    before_source = before_norm["source"][before_i]
    before_sink = before_norm["sink"][before_i]
    snarl_size_before_after[before_source] = [before_norm["size"][before_i], 0]
    while before_sink != after_norm["sink"][after_i]:
        #we haven't reached a new snarl in before_norm, because the snarl has been split
        # into mulitple parts in after_norm.
        snarl_size_before_after[before_source][1] += after_norm["size"][after_i]
        snarls_that_split[before_source] += 1
        after_i += 1
    
    #we've reached the last snarl in after_norm that corresponds to the the current snarl in before_norm.
    snarl_size_before_after[before_source][1] += after_norm["size"][after_i]
    if before_source in snarls_that_split:
        snarls_that_split[before_source] += 1
    after_i += 1


# %%
# reality check: do we see the number of snarls that we'd expect? A: we do!
print(len(snarl_size_before_after))
additional_splits = int()
for splits in snarls_that_split.values():
    additional_splits+=splits-1
print(additional_splits)
print("total number of snarls in after_norm (reality check):", len(snarl_size_before_after) + additional_splits)

#%%
#* show the size change of each snarl with position along x and size change along y.
snarl_size_change = col.OrderedDict()
shrinking_snarl_num = int()
shrinking_snarl_before_norm_size = list()
shrinking_snarl_after_norm_size = list()
growing_snarl_num = int()
growing_snarl_before_norm_size = list()
growing_snarl_after_norm_size = list()
unchanging_snarl_num = int()
for source, before_after in snarl_size_before_after.items():
    snarl_size_change[source] = before_after[1] - before_after[0]
    if (before_after[1] - before_after[0]) < 0:
        shrinking_snarl_num += 1
        shrinking_snarl_before_norm_size.append(before_after[0])
        shrinking_snarl_after_norm_size.append(before_after[1])
    if (before_after[1] - before_after[0]) == 0:
        unchanging_snarl_num += 1
    if (before_after[1] - before_after[0]) > 0:
        growing_snarl_num += 1
        growing_snarl_before_norm_size.append(before_after[0])
        growing_snarl_after_norm_size.append(before_after[1])

print("shrinking_snarl_before_norm_size mean", np.mean(shrinking_snarl_before_norm_size))
print("shrinking_snarl_after_norm_size mean", np.mean(shrinking_snarl_after_norm_size))
print("growing_snarl_before_norm_size mean", np.mean(growing_snarl_before_norm_size))
print("growing_snarl_after_norm_size mean", np.mean(growing_snarl_after_norm_size))

print("shrinking_snarl_before_norm_size range", np.min(shrinking_snarl_before_norm_size), np.max(shrinking_snarl_before_norm_size))
print("shrinking_snarl_after_norm_size range", np.min(shrinking_snarl_after_norm_size), np.max(shrinking_snarl_after_norm_size))
print("growing_snarl_before_norm_size range", np.min(growing_snarl_before_norm_size), np.max(growing_snarl_before_norm_size))
print("growing_snarl_after_norm_size range", np.min(growing_snarl_after_norm_size), np.max(growing_snarl_after_norm_size))




print("shrinking_snarl_num", shrinking_snarl_num)
print("growing_snarl_num", growing_snarl_num)
print("unchanging_snarl_num", unchanging_snarl_num)

snarl_size_change_without_splits = col.OrderedDict()
snarl_size_change_only_splits = col.OrderedDict()
for x, y in snarl_size_change.items():
    if x not in snarls_that_split:
        snarl_size_change_without_splits[x] = y
    else:
        snarl_size_change_only_splits[x] = y
print("all snarls:")
plt.scatter(x=snarl_size_change.keys(), y=snarl_size_change.values(), s=8)
plt.show()
print("all snarls without splits:")
plt.scatter(x=snarl_size_change_without_splits.keys(), y=snarl_size_change_without_splits.values(), s=8)
plt.show()

print("only snarls with splits:")
plt.scatter(x=snarl_size_change_only_splits.keys(), y=snarl_size_change_only_splits.values(), s=8)
plt.show()

print("all snarls, dif colors")
plt.scatter(x=snarl_size_change_without_splits.keys(), y=snarl_size_change_without_splits.values(), s=8)
plt.scatter(x=snarl_size_change_only_splits.keys(), y=snarl_size_change_only_splits.values(), s=8)
plt.show()
#%%
shrink = int()
same = int()
grow = int()
for x in snarl_size_change_without_splits.values():
    if x < 0:
        shrink += 1
    if x == 0:
        same += 1
    if x > 0:
        grow += 1
print(shrink, same, grow)
# print(len([x if x < 0 for x in snarl_size_change_without_splits.values()]))
#%%
shrink_splits = list(snarl_size_change_only_splits.values())
most_shrink_split = min(shrink_splits)
most_shrink_split_i = shrink_splits.index(most_shrink_split)
most_shrink_snarl = list(snarl_size_change_only_splits.keys())[most_shrink_split_i]
print(most_shrink_snarl)
print(snarl_size_change_only_splits[most_shrink_snarl])
print(list(snarl_size_change.keys()).index(most_shrink_snarl))
#%%
"""
This is for visualizing a good example of a large snarl in yeast_subset that split 
into five subsnarls. Original snarl coordinates: 3881503-3881729
"""

print(before_norm.loc[list(snarl_size_change.keys()).index(most_shrink_snarl)])
# df.loc[df['favorite_color'] == 'yellow']
print(after_norm.loc[after_norm["source"] == most_shrink_snarl])
print(after_norm.loc[58342])
print(after_norm.loc[58343])
print(after_norm.loc[58344])
print(after_norm.loc[58345])
print(after_norm.loc[58346])
#%%
print(len(snarl_size_change_only_splits))

#%%
"""
plot arranged with x axis is size of snarl, y axis is change of snarl size by percentage.
Idea is to find relatively small snarls with large change of snarl size, either favorably or unfavorably.
"""
after_i = int()
snarl_num = list()
size = list()
size_change = list()
debug_stop = int()
for before_i in range(len(before_norm)):
    ## skip the one massive snarl:
    # if before_norm["size"][before_i] == 14380261:
    #     after_i += 1
    #     continue
    # if debug_stop == 10:
    #     break
    # debug_stop += 1
    # print('before_norm["source"][before_i]', before_norm["source"][before_i])
    # print('after_norm["source"][after_i]', after_norm["source"][after_i])
    if before_norm["source"][before_i] == after_norm["source"][after_i]:
        if before_norm["sink"][before_i] == after_norm["sink"][after_i]:
            # we have a snarl that is represented fully in both graphs.
            size.append(before_norm["size"][before_i])
            snarl_num.append((before_i, after_i))
            size_change.append(after_norm["size"][after_i] - before_norm["size"][before_i])
            after_i += 1
        else:
            # print("in else statement, while bool: ", (before_norm["source"][before_i + 1] != after_norm["source"][after_i]))
            # we have a snarl that is split during normalization. Scan for next before_norm source value in after_norm.
            while before_norm["source"][before_i + 1] != after_norm["source"][after_i]:
                # print("next source in before_norm", before_norm["source"][before_i + 1], "current source in after_norm:", after_norm["source"][after_i])
                after_i += 1

size_change_percent = [(change/size)*100 for (change, size) in zip(size_change, size)]
print(size_change_percent[:5])


# size_change_sorted = [x for _, x in sorted(zip(size, size_change))]
# size_sorted = sorted(size)
# sort_indices = np.argsort(size)
plt.xscale("log")
plt.scatter(x=size, y= size_change_percent, s=8)
#%%

plt.xscale("log")
plt.scatter(x=size, y= size_change_percent, s=8)
# plt.xscale("linear")
#%%
max_change = max(size_change_percent)
max_change_i = size_change_percent.index(max_change)
print(max_change_i)
snarl_num[max_change_i]
print(before_norm.loc[snarl_num[max_change_i][0]])
print(after_norm.loc[snarl_num[max_change_i][1]])
#%%
"""
Yeast_subset normalization info:
normalized 54820 snarls, skipped 2 snarls because. . .
they exceeded the size limit (0 snarls),
had haplotypes starting/ending in the middle of the snarl (0),
there were handles not connected by the gbwt info (0 snarls),
the snarl was cyclic (2 snarls),
or the snarl was trivial - composed of only one or two nodes (0 snarls).
amount of sequence in normalized snarls before normalization: 880131
amount of sequence in normalized snarls after normalization: 879739
Elapsed time: 29.4443 s
"""

#%%

# snarl_size_change = dict() #key: source of original snarl, value: change in size (negative if reduced in size)
# snarls_that_split = dict() #key: source of original snarl, value: number of snarls it is split into.

# prev_source = int()
# after_i = 0
# prev_i = 0
# while (prev_i < prev_norm["source"].size() and after_i < after_norm["source"].size())
#     if source == after_norm[after_i]:
#         #we've reached a new snarl in before_norm, and it's still in after_norm. 
#         snarl_size_change[source] = after_norm[after_i]
#         prev_source = source
#         prev_i
#         after_i += 1
#     else:
#         #we haven't reached a new snarl in before_norm, because we're still in a snarl in after_norm.
#         snarl_size_change[prev_source]
    