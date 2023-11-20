GRAPHS_DIR="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/graphs"
READS_DIR="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/reads"
REFS_DIR="/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/references"
WORK_DIR="trash/exp"

wildcard_constraints:
    trimmedness="\\.trimmed|",
    sample=".+(?<!\\.trimmed)"

def repetitive_kmers(wildcards):
    """
    Find the Winnowmap repetitive kmers file from a reference.
    """
    return os.path.join(REFS_DIR, wildcards["reference"] + "-pansn.repetitive_k15.txt")

def reference_fasta(wildcards):
    """
    Find the linear reference FASTA from a reference.
    """
    return os.path.join(REFS_DIR, wildcards["reference"] + "-pansn.fa")

def graph_base(wildcards):
    """
    Find the base name for a collection of graph files from reference.
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
    Find a FASTQ from realness, tech, sample, trimmedness, and subset.

    Works even if there is extra stuff in the name besides sample. Accounts for
    being able to make a FASTQ from a GAM.
    """
    import glob
    fastq_pattern = os.path.join(READS_DIR, "{realness}/{tech}/*{sample}*{trimmedness}[._-]{subset}.f*q".format(**wildcards))
    fastq_by_sample_pattern = os.path.join(READS_DIR, "{realness}/{tech}/{sample}/*{sample}*{trimmedness}[._-]{subset}.f*q".format(**wildcards))
    results = glob.glob(fastq_pattern) + glob.glob(fastq_by_sample_pattern)
    if len(results) == 0:
        # Maybe there's a GAM to extract from? GAMs are always under per-sample directories.
        gam_pattern = os.path.join(READS_DIR, "{realness}/{tech}/{sample}/*{sample}*{trimmedness}[._-]{subset}.gam".format(**wildcards))
        results = glob.glob(gam_pattern)
        if len(results) == 0 and wildcards["realness"] == "sim":
            # TODO: We give up and assume we can make this subset.
            results = [os.path.join(READS_DIR, "{realness}/{tech}/{sample}/{sample}-{realness}-{tech}{trimmedness}-{subset}.gam".format(**wildcards))]
        if len(results) > 1:
            raise AmbiguousRuleException("Multiple files matched " + gam_pattern)
        # Replace the extension
        return results[0][:-3] + "fq"
    if len(results) > 1:
        raise AmbiguousRuleException("Multiple files matched " + fastq_pattern + " and " + fastq_by_sample_pattern)
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

    if "experiments" not in config:
        raise RuntimeError(f"No experiments section in configuration; cannot run experiment {expname}")
    all_experiments = config["experiments"]
    
    if expname not in all_experiments:
        raise RuntimeError(f"Experiment {expname} not in configuration")
    exp_dict = all_experiments[expname]

    # Make a base dict of all controlled variables.
    base_condition = exp_dict.get("control", {})

    to_vary = exp_dict.get("vary", {})

    to_constrain = exp_dict.get("constrain", [])

    total_conditions = 0
    for condition in augmented_with_all(base_condition, to_vary):
        # For each combination of independent variables on top of the base condition

        # We need to see if this is a combination we want to do
        
        if len(to_constrain) == 0 or matches_any_constraint(condition, to_constrain):
            total_conditions += 1
            yield condition
        else:
            print(f"Condition {condition} does not match a constraint")
    print(f"Experiment {expname} has {total_conditions} conditions")
    

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


def matches_constraint(condition, constraint, debug=False):
    """
    Returns True if all keys in constraint are in condition with the same
    values.
    """
    for k, v in constraint.items():
        if k not in condition or condition[k] != v:
            if debug:
                print(f"Condition {condition} mismatched constraint {constraint} on {k}")
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
    Determine a human-readable condition name from expname and the experiment's variable values.
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

def all_experiment(wildcard_values, pattern, debug=False):
    """
    Produce all values of pattern substituted with the wildcards and the experiment conditions' values, from expname.
    
    Needs to be used like:
        lambda w: all_experiment(w, "your pattern")
    """

    for condition in all_experiment_conditions(wildcard_values["expname"]):
        merged = dict(wildcard_values)
        merged.update(condition)
        if debug:
            print(f"Evaluate {pattern} in {merged} from {wildcard_values} and {condition}")
        filename = pattern.format(**merged)
        yield filename

def winnowmap_mode(wildcards):
    """
    Determine the right Winnowmap preset (map-pb, etc.) from tech.
    """

    return {
        "r9": "map-ont",
        "r10": "map-ont",
        "hifi": "map-pb"
    }[wildcards["tech"]]

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

rule alias_gam_k:
    input:
        gam="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{part_subset}000.gam"
    output:
        gam="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{part_subset}k.gam"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5
    shell:
        "ln {input.gam} {output.gam}"

rule alias_gam_m:
    input:
        gam="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{part_subset}000000.gam"
    output:
        gam="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{part_subset}m.gam" 
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5
    shell:
        "ln {input.gam} {output.gam}"

rule extract_fastq:
    input:
        gam="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"
    output:
        fastq="{reads_dir}/sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.fq"
    threads: 16
    resources:
        mem_mb=10000,
        runtime=60
    shell:
        "vg view --threads {threads} {input.gam} >{output.fastq}"

rule giraffe_real_reads:
    input:
        unpack(indexed_graph),
        fastq=fastq,
    output:
        gam="{root}/aligned/{reference}/giraffe-{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        realness="real"
    threads: 16
    resources:
        mem_mb=300000,
        runtime=240
    shell:
        "vg giraffe -t{threads} --parameter-preset lr --progress --track-provenance -Z {input.gbz} -d {input.dist} -m {input.minfile} -z {input.zipfile} -f {input.fastq} >{output.gam}"

rule giraffe_sim_reads:
    input:
        unpack(indexed_graph),
        gam=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"),
    output:
        gam="{root}/aligned/{reference}/giraffe-{minparams}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    wildcard_constraints:
        realness="sim"
    threads: 16
    resources:
        mem_mb=300000,
        runtime=60
    shell:
        "vg giraffe -t{threads} --parameter-preset lr --progress --track-provenance --track-correctness -Z {input.gbz} -d {input.dist} -m {input.minfile} -z {input.zipfile} -G {input.gam} >{output.gam}"

rule winnowmap_reads:
    input:
        reference_fasta=reference_fasta,
        repetitive_kmers=repetitive_kmers,
        fastq=fastq
    params:
        winnowmap_mode=winnowmap_mode
    output:
        bam="{root}/aligned/{reference}/winnowmap/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    threads: 16
    resources:
        mem_mb=300000,
        runtime=120
    shell:
        "winnowmap -t 15 -W {input.repetitive_kmers} -ax {params.winnowmap_mode} {input.reference_fasta} {input.fastq} | samtools view -h -F 2048 -F 256 --bam - >{output.bam}"

rule inject_bam:
    input:
        gbz=gbz,
        bam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.bam"
    output:
        gam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam"
    threads: 16
    resources:
        mem_mb=300000,
        runtime=120
    shell:
        "vg inject --threads {threads} -x {input.gbz} {input.bam} >{output.gam}"

rule annotate_and_compare_alignments:
    input:
        gbz=gbz,
        gam="{root}/aligned/{reference}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.gam",
        truth_gam=os.path.join(READS_DIR, "sim/{tech}/{sample}/{sample}-sim-{tech}-{subset}.gam"),
    output:
        gam="{root}/annotated/{reference}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.gam",
        tsv="{root}/compared/{reference}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.compared.tsv",
        report="{root}/compared/{reference}/{mapper}/sim/{tech}/{sample}{trimmedness}.{subset}.compare.txt"
    threads: 8
    resources:
        mem_mb=25000,
        runtime=60
    shell:
        "vg annotate -t7 -a {input.gam} -x {input.gbz} -m | tee >{output.gam} | vg gamcompare --range 200 - {input.truth_gam} -T > {output.tsv} 2>{output.report}"

rule correctness_from_comparison:
    input:
        report="{root}/compared/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.compare.txt"
    params:
        condition_name=condition_name
    output:
        correct="{root}/experiments/{expname}/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.correct.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
    shell:
        "printf '{params.condition_name}\\t' >{output.correct} && cat {input.report} | grep ' reads correct$' | cut -f1 -d' ' >>{output.correct}"

rule experiment_correctness_table:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.correct.tsv")
    output:
        table="{root}/experiments/{expname}/results/correct.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10
    shell:
        "cat {input} >{output.table}"

rule experiment_correctness_plot:
    input:
        tsv="{root}/experiments/{expname}/results/correct.tsv"
    output:
        "{root}/experiments/{expname}/plots/correct.{ext}"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5
    shell:
        "barchart.py {input.tsv} --title '{wildcards.expname} Correctness' --y_label 'Correct Reads' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule stats_from_alignments:
    input:
        gam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        stats="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gamstats.txt"
    threads: 16
    resources:
        mem_mb=10000,
        runtime=30
    shell:
        "vg stats -p {threads} -a {input.gam} >{output.stats}"

rule mapping_rate_from_stats:
    input:
        stats="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gamstats.txt"
    params:
        condition_name=condition_name
    output:
        rate="{root}/experiments/{expname}/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_rate.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5
    shell:
        "printf '{params.condition_name}\\t' >{output.rate} && cat {input.stats} | grep 'Total aligned:' | cut -f2 -d':' | tr -d ' ' >>{output.rate}"

rule experiment_mapping_rate_table:
    input:
        lambda w: all_experiment(w, "{root}/experiments/{expname}/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.mapping_rate.tsv")
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
        "barchart.py {input.tsv} --title '{wildcards.expname} Mapping Rate' --y_label 'Mapped Reads' --x_label 'Condition' --x_sideways --no_n --save {output}"

rule chain_coverage_alignments:
    input:
        gam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_coverage.tsv"
    threads: 2
    wildcard_constraints:
        mapper="giraffe"
    resources:
        mem_mb=2000,
        runtime=120
    shell:
        "vg view -aj {input.gam} | jq -r '.annotation.best_chain_coverage' >{output}"

rule chain_coverage_histogram:
    input:
        tsv="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.best_chain_coverage.tsv"
    output:
        "{root}/plots/{reference}/{mapper}/best_chain_coverage-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    wildcard_constraints:
        mapper="giraffe"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=10
    shell:
        "histogram.py {input.tsv} --bins 100 --title '{wildcards.tech} {wildcards.realness} Fraction Covered' --y_label 'Items' --x_label 'Coverage' --no_n --save {output}"

rule read_length_alignments:
    input:
        gam="{root}/aligned/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.gam",
    output:
        "{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_mapping.tsv"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=120
    shell:
        "vg view -aj {input.gam} | jq -r '[if (.path.mapping // []) == [] then \"unmapped\" else \"mapped\" end, (.sequence | length)] | @tsv' >{output}"

rule read_length_histogram:
    input:
        tsv="{root}/stats/{reference}/{mapper}/{realness}/{tech}/{sample}{trimmedness}.{subset}.length_by_mapping.tsv"
    output:
        "{root}/plots/{reference}/{mapper}/length_by_mapping-{realness}-{tech}-{sample}{trimmedness}.{subset}.{ext}"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=10
    shell:
        "histogram.py {input.tsv} --bins 100 --title '{wildcards.tech} {wildcards.realness} Read Length' --y_label 'Items' --x_label 'Length (bp)' --no_n --legend_overlay best --save {output}"





