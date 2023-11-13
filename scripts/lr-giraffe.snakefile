REFERENCES=["chm13"]
INDEX_PARAM_SETS=["k31.w50.W"]
SAMPLES=["HG002"]
REALNESSES=["real", "sim"]
TECHS=["r9", "r10", "hifi"]

GRAPHS_DIR="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs"
READS_DIR="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/reads"
WORK_DIR="trash/exp"

wildcard_constraints:
    trimmedness="\\.trimmed|",
    sample=".+(?<!\\.trimmed)"

def graph_base(wildcards):
    """
    Find the base name for a collection fo graph files from reference.
    """
    return os.path.join(GRAPHS_DIR, "hprc-v1.1-mc-" + wildcards["reference"] + ".d9")

def gbz(wildcards):
    """
    Find a graph GBZ file from reference.
    """
    return graph_base(wildcards) + ".gbz"

def dist_indexed_graph(wildcards):
    """
    Find a GBZ and its dist index from reference.
    """
    base = graph_base(wildcards)
    return {
        "gbz": gbz(wildcards),
        "dist": base + ".dist"
    }

def indexed_graph(wildcards):
    """
    Find an indexed graph and all its indexes from reference and minparams.
    """
    base = graph_base(wildcards)
    indexes = dist_indexed_graph(wildcards)
    new_indexes = {
        "minfile": base + "." + wildcards["minparams"] + ".withzip.min",
        "zipfile": base + "." + wildcards["minparams"] + ".zipcodes"
    }
    new_indexes.update(indexes)
    return new_indexes

def fastq(wildcards):
    """
    Find a FASTQ from realness, tech, sample, trimmedness, and subset, even if there is extra stuff in the name besides sample.
    """
    import glob
    pattern = os.path.join(READS_DIR, "{realness}/{tech}/*{sample}*{trimmedness}[._]{subset}.f*q".format(**wildcards))
    results = glob.glob(pattern)
    if len(results) == 0:
        raise FileNotFoundError("Nothing matched " + pattern)
    if len(results) > 1:
        raise AmbiguousRuleException("Multiple files matched " + pattern)
    return results[0]

def all_experiment_conditions(expname):
    """
    Yield dictionaries of all conditions for the given experiment.
    
    The config file should have a dict in "experiments", of which the given
    expname should be a key. THe value is the experiment dict.

    The experiment dict should have a "control" dict, listing names and values
    of variables to keep constant.

    The experiment dict should have a "vary" dict, listing names and values
    lists of variables to vary. All combinations will be generated.

    The experiment dict should have a "constrain" list. Each item is a dict of
    variable names and values. A condition must match *at least* one of these
    dicts on *all* values in the dict in order to pass.

    Yields variable name to value dicts for all passing conditions for the
    given experiment.
    """

    exp_dict = config.get("experiments", {}).get(expname, {})

    # Make a base dict of all controlled variables.
    base_condition = exp_dict.get("control", {})

    to_vary = exp_dict.get("vary", {})

    to_constrain = exp_dict.get("constrain", [])

    for condition in augmented_with_all(base_condition, to_vary):
        # For each combination of independent variables on top of the base condition

        # We need to see if this is a combination we want to do
        
        if len(to_constrain) == 0 or matches_any_constraint(condition, to_constrain):
            yield condition
    

def augmented_with_each(base_dict, new_key, possible_values):
    """
    Yield copies of base_dict with each value from possible_values under new_key.
    """

    for value in possible_values:
        clone = dict(base_dict)
        clone[new_key] = value
        yield clone

def augmented_with_all(base_dict, keys_and_values):
    """
    Yield copies of base_dict augmented with all combinations of values from
    keys_and_values, under the corresponding keys.
    """

    if len(keys_and_values) == 0:
        # Base case: nothing to add
        yield base_dict
    else:
        # Break off one facet
        first_key = next(iter(keys_and_values.keys()))
        first_values = keys_and_values[first_key]
        rest = dict(keys_and_values)
        del rest[first_key]
        for with_rest in augmented_with_all(base_dict, rest):
            # Augment with the rest
            for with_first in augmented_with_each(with_rest, first_key, first_values):
                # And augment with this key
                yield with_first


def matches_constraint(condition, constraint):
    """
    Returns True if all keys in constraint are in condition with the same
    values.
    """
    for k, v in constraint.items():
        if k not in condition or condition[k] != v:
            return False
    return True

def matches_any_constraint(condition, constraints):
    """
    Return True if, for some constraint dict, the condition dict matches all
    values in the constraint dict.
    """

    for constraint in constraints:
        if matches_constraint(condition, constraint):
            return True
    return False

def wildcards_to_condition(all_wildcards):
    """
    Filter dowen wildcards to just the condition parameters for the experiment in expname.
    
    Raises an error if any variable in the experiment cannot be determined.
    """

    exp_dict = config.get("experiments", {}).get(all_wildcards["expname"], {})
    base_condition = exp_dict.get("control", {})
    to_vary = exp_dict.get("vary", {})
    all_vars = list(base_condition.keys()) + list(to_vary.keys())

    condition = {}

    for var in all_vars:
        condition[var] = all_wildcards[var]

    return condition

def condition_name(wildcards):
    """
    Determine a human-readable condition name from expname, reference, minparams, realness, tech, sample, trimmedness, and subset.
    """
    
    # Get what changes in the experiment
    exp_dict = config.get("experiments", {}).get(wildcards["expname"], {})
    to_vary = exp_dict.get("vary", {})

    # Get the condition dict in use here
    condition = wildcards_to_condition(wildcards)
    
    # Paste together all the varied variable values from the condition.
    varied = list(to_vary.keys())
    varied_values = [condition[v] for v in varied]
    return ",".join(varied_values)

def all_experiment_mapping_rate_stats(wildcards):
    """
    Produce the names of all mapping rate stats files for the current experiment, form expname and root.
    """
    
    for condition in all_experiment_conditions(wildcards["expname"]):
        filename = wildcards["root"] + "/experiments/" + wildcards["expname"] + "/{reference}/{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_rate.tsv".format(**condition)
        yield filename


rule minimizer_index_graph:
    input:
        unpack(dist_indexed_graph)
    output:
        minfile="{graphs_dir}/hprc-v1.1-mc-{reference}.d9.k{k}.w{w}{weightedness}.withzip.min",
        zipfile="{graphs_dir}/hprc-v1.1-mc-{reference}.d9.k{k}.w{w}{weightedness}.zipcodes"
    wildcard_constraints:
        weightedness="\\.W|",
        k="[0-9]+",
        w="[0-9]+"
    threads: 16
    resources:
        mem_mb=80000,
        runtime=240
    shell:
        "vg minimizer --progress -k {wildcards.k} -w {wildcards.w} -t {threads} -p -d {input.dist} -z {output.zipfile} -o {output.minfile} {input.gbz}"


rule align_real_reads:
    input:
        unpack(indexed_graph),
        fastq=fastq,
    output:
        gam="{root}/aligned/{reference}/{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        realness="real"
    threads: 16
    resources:
        mem_mb=300000,
        runtime=240
    shell:
        "vg giraffe -t{threads} --parameter-preset lr --progress --track-provenance -Z {input.gbz} -d {input.dist} -m {input.minfile} -z {input.zipfile} -f {input.fastq} >{output.gam}"

rule align_sim_reads:
    input:
        unpack(indexed_graph),
        gam=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"),
    output:
        gam="{root}/aligned/{reference}/{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        realness="sim"
    threads: 16
    resources:
        mem_mb=300000,
        runtime=60
    shell:
        "vg giraffe -t{threads} --parameter-preset lr --progress --track-provenance --track-correctness -Z {input.gbz} -d {input.dist} -m {input.minfile} -z {input.zipfile} -G {input.gam} >{output.gam}"

rule annotate_and_compare_alignments:
    input:
        gbz,
        gam="{root}/aligned/{reference}/{minparams}/sim/{tech}/{sample}{trimmedness}.{subset}.gam",
        truth_gam="{READS_DIR}/sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam",
    output:
        gam="{root}/annotated/{reference}/{minparams}/sim/{tech}/{sample}{trimmedness}.{subset}.gam",
        tsv="{root}/compared/{reference}/{minparams}/sim/{tech}/{sample}{trimmedness}.{subset}.compared.tsv",
        report="{root}/compared/{reference}/{minparams}/sim/{tech}/{sample}{trimmedness}.{subset}.compare.txt"
    threads: 8
    resources:
        mem_mb=25000,
        runtime=60
    shell:
        "vg annotate -t{threads - 1} -a {input.gam} -x {input.gbz} -m | tee >{output.gam} | vg gamcompare --range 200 - {input.truth_gam} -T > {output.tsv} 2>{output.report}"

rule stats_alignments:
    input:
        gam="{root}/aligned/{reference}/{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        stats="{root}/stats/{reference}/{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gamstats.txt"
    threads: 16
    resources:
        mem_mb=10000,
        runtime=30
    shell:
        "vg stats -p {threads} -a {input.gam} >{output.stats}"

rule mapping_rate_stats:
    input:
        stats="{root}/stats/{reference}/{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gamstats.txt"
    params:
        condition_name=condition_name
    output:
        rate="{root}/experiments/{expname}/{reference}/{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_rate.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5
    shell:
        "printf '{params.condition_name}\\t' >{output.rate} && cat {input.stats} | grep 'Total aligned:' | cut -f2 -d':' | tr -d ' ' >>{output.rate}"

rule experiment_mapping_rate_table:
    input:
        all_experiment_mapping_rate_stats
    output:
        table="{root}/experiments/{expname}/results/mapping_rate.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5
    shell:
        "cat {input} >{output.table}"

rule experiment_mapping_rate_plot:
    input:
        tsv="{root}/experiments/{expname}/results/mapping_rate.tsv"
    output:
        "{root}/experiments/{expname}/plots/mapping_rate.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5
    shell:
        "barchart.py {input.tsv} --title '{wildcards.expname} Mapping Rate' --y_label 'Mapped Reads' --x_label 'Condition' --no_n --save {output}"

rule chain_coverage_alignments:
    input:
        gam="{root}/aligned/{reference}/{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_coverage.tsv"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=120
    shell:
        "vg view -aj {input.gam} | jq -r '.annotation.best_chain_coverage' >{output}"

rule chain_coverage_histogram:
    input:
        tsv="{root}/stats/{reference}/{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_coverage.tsv"
    output:
        "{root}/plots/{reference}/{minparams}/best_chain_coverage-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=10
    shell:
        "histogram.py {input.tsv} --bins 100 --title '{wildcards.tech} {wildcards.realness} Fraction Covered' --y_label 'Items' --x_label 'Coverage' --no_n --save {output}"

rule read_length_alignments:
    input:
        gam="{root}/aligned/{reference}/{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_mapping.tsv"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=120
    shell:
        "vg view -aj {input.gam} | jq -r '[if (.path.mapping // []) == [] then \"unmapped\" else \"mapped\" end, (.sequence | length)] | @tsv' >{output}"

rule read_length_histogram:
    input:
        tsv="{root}/stats/{reference}/{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_mapping.tsv"
    output:
        "{root}/plots/{reference}/{minparams}/length_by_mapping-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=10
    shell:
        "histogram.py {input.tsv} --bins 100 --title '{wildcards.tech} {wildcards.realness} Read Length' --y_label 'Items' --x_label 'Length (bp)' --no_n --legend_overlay best --save {output}"





