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

rule minimizer_index_graph:
    input:
        unpack(dist_indexed_graph)
    output:
        minfile="{graphs_dir}/hprc-v1.1-mc-{reference}.d9.k{k}.w{w}{weightedness}.withzip.min",
        zipfile="{graphs_dir}/hprc-v1.1-mc-{reference}.d9.k{k}.w{w}{weightedness}.zipcodes.min"
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





