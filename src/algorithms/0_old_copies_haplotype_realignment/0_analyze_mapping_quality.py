#%%

import json
import numpy as np
import matplotlib.pyplot as plt
import collections as col


#%%

def measure_map_qualities(json_file):
    mapping_qualities = list()
    mapping_qualities_not_60 = list()
    no_quality = int()
    more_than_1_qual = int()
    with open(mapped_reads_file) as file:
        line_ct = 0
        for line in file:
            line_ct +=1
            # print(line_ct)
            line_split = (line.split('"mapping_quality":'))
            if len(line_split) == 2:
                # print(line_split[1])
                # print(([line_split[1].split(',"')]))
                # print(([line_split[1].split(',"')])[0])
                # print()
                map_qual = int(line_split[1].split(',"')[0])
                if map_qual != 60:
                    mapping_qualities_not_60.append(map_qual)
                mapping_qualities.append(int(line_split[1].split(',"')[0]))
            elif len(line_split) == 1:
                no_quality += 1 
            else:
                more_than_1_qual += 1

    print("mean of mapping_qualities", np.mean(mapping_qualities))
    print("mean of mapping_qualities that aren't perfect (!= 60)", np.mean(mapping_qualities_not_60))
    print("no_quality", no_quality)
    print("more_than_1_qual", more_than_1_qual)

    qual_count = col.Counter(mapping_qualities)
    print(qual_count)
    # plt.plot(qual_count)
    plt.plot(mapping_qualities_not_60)
    # plt.bar(mapping_qualities_not_60, 60)

#%%


# mapped_reads_file = "~/test/robin_tests/full_chr10/test_map/hgsvc_chr10_construct_test_map.json"
mapped_reads_file = "hgsvc_chr10_construct_test_map.json"
# mapped_reads_file = "0_test.txt"
measure_map_qualities(mapped_reads_file)


#%%
"""
Current output with no normalization:
map_qual avg: 59.059475806451616
22.68
no_quality 8
more_than_1_qual 0
"""














"""
import json
import numpy as np
mapped_reads_file = "../../../test/robin_tests/full_chr10/test_map/hgsvc_chr10_construct_test_map.gam"

mapping_qualities = list()
mapping_qualities_not_60 = list()
no_quality = int()
more_than_1_qual = int()
with open(mapped_reads_file) as file:
    line_ct = 0
    for line in file:
        line_ct +=1
        # print(line_ct)
        line_split = (line.split('"mapping_quality":'))
        if len(line_split) == 2:
            # print(line_split[1])
            # print(([line_split[1].split(',"')]))
            # print(([line_split[1].split(',"')])[0])
            # print()
            map_qual = int(line_split[1].split(',"')[0])
            if map_qual != 60:
                mapping_qualities_not_60.append(map_qual)
            mapping_qualities.append(int(line_split[1].split(',"')[0]))
        elif len(line_split) == 1:
            no_quality += 1 
        else:
            more_than_1_qual += 1
# print(mapping_qualities)
print(np.mean(mapping_qualities))
print(np.mean(mapping_qualities_not_60))
print("no_quality", no_quality)
print("more_than_1_qual", more_than_1_qual)


    # mapped_reads = json.load(file)
    # print(mapped_reads)



"""