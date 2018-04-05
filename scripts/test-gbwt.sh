#!/usr/bin/env bash
# test-gbwt.sh: Plot the effect of haplotype information on mapping performance

# We are going to compare 5 mapping regimes:
# 1. The snp1kg graph without the variants in the sample (negative control)
# 2. The full snp1kg graph (Heng Effect positive control)
# 3. The full snp1kg graph with GBWT haplotype information (under test)
# 4. The frequency-filtered minaf snp1kg graph (current best)
# 5. The primary path graph (Heng Effect negative control)

# We want to know if the GBWT reduces the Heng Effect (poor mapping due to
# spurious alignments to rare variants)

# OK so on to the actual code.
set -ex

# Define constants

# Define a region name to process. This sets the name that the graphs and
# indexes will be saved/looked for under.
REGION_NAME="CHR21"
# Define the VCF and FASTA basenames. We assume the VCF has a TBI.
VCF_BASENAME="1kg_hg19-CHR21.vcf.gz"
FASTA_BASENAME="CHR21.fa"
# Define where to get them
SOURCE_BASE_URL="s3://cgl-pipeline-inputs/vg_cgl/bakeoff"
# Define the contig we are using
GRAPH_CONTIG="21"
# Define the region to build the graph on, as contig[:start-end]
GRAPH_REGION="${GRAPH_CONTIG}"


# Set a FASTQ to model reads after
TRAINING_FASTQ="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_005_BH814YADXX/Project_RM8398/Sample_U5a/U5a_AGTCAA_L002_R1_007.fastq.gz"
# And a read simulation seed
READ_SEED="90"
# And a read count
READ_COUNT="10000000"
# Chunks to simulate in (which affects results)
READ_CHUNKS="32"

MODE="21"
if [[ "${MODE}" == "mhc" ]]; then
    # Actually do a smaller test
    READ_COUNT="100000"
    REGION_NAME="MHC"
    GRAPH_CONTIG="6"
    GRAPH_REGION="${GRAPH_CONTIG}:28510119-33480577"
    FASTA_BASENAME="chr6.fa.gz"
    VCF_BASENAME="1kg_hg38-MHC.vcf.gz"
elif [[ "${MODE}" == "tiny" ]]; then
    # Do just 20 kb of MHC and a very few reads
    READ_COUNT="1000"
    REGION_NAME="MHC"
    GRAPH_CONTIG="6"
    GRAPH_REGION="${GRAPH_CONTIG}:28510119-28520119"
    FASTA_BASENAME="chr6.fa.gz"
    VCF_BASENAME="1kg_hg38-MHC.vcf.gz"
fi

# Define the sample to use for synthesizing reads
SAMPLE_NAME="HG00096"

# What min allele frequency limit do we use?
MIN_AF="0.0335570469"

# Do we actually want to run the mapeval jobs? Or just do plotting?
RUN_JOBS="1"


# Now we need to parse our arguments
usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [Options] TREE_PATH GRAPHS_PATH OUTPUT_PATH\n"
    printf "\tTREE_PATH\ta Toil job tree location\n"
    printf "\tGRAPHS_PATH\ta directory of graphs (which, if extant, allows graph construction to be skipped)\n"
    printf "\tOUTPUT_PATH\ta results directory\n"
    exit 1
}

while getopts "" o; do
    case "${o}" in
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

if [[ "$#" -lt "3" ]]; then
    # Too few arguments
    usage
fi

TREE_PATH="${1}"
shift
GRAPHS_PATH="${1}"
shift
OUTPUT_PATH="${1}"
shift

if [[ -e "${TREE_PATH}" ]]; then
    # Make sure we don't clobber the first arg on accident
    echo "ERROR: Tree path ${TREE_PATH} already exists and needs to be removed" 1>&2
    exit 1
fi

mkdir -p "${TREE_PATH}"

# Generate a toil-vg config
toil-vg generate-config > "${TREE_PATH}/toil-vg.conf"
sed -i "s/alignment-cores:.*/alignment-cores: 64/" "${TREE_PATH}/toil-vg.conf"
sed -i "s/alignment-disk:.*/alignment-disk: '20G'/" "${TREE_PATH}/toil-vg.conf"

# Now we need to make sure our graphs exist and are downloaded

if [[ ! -d "${GRAPHS_PATH}" ]]; then
    # Graphs need to be gotten
    
    if [[ -e "${GRAPHS_PATH}" ]]; then
        # It needs to not exist at all
        echo "ERROR: Graph path ${GRAPHS_PATH} is not a directory" 1>&2
        exit 1
    fi
    
    # Make the directory
    mkdir -p "${GRAPHS_PATH}"
    
    # Construct the graphs
    # Hardcoded constants here go with the snp1kg URLs above.
    toil-vg construct "${TREE_PATH}/construct" "${GRAPHS_PATH}" \
        --config "${TREE_PATH}/toil-vg.conf" \
        --vcf "${SOURCE_BASE_URL}/${VCF_BASENAME}" \
        --fasta "${SOURCE_BASE_URL}/${FASTA_BASENAME}" \
        --out_name "snp1kg-${REGION_NAME}" \
        --alt_paths \
        --realTimeLogging \
        --control_sample "${SAMPLE_NAME}" \
        --haplo_sample "${SAMPLE_NAME}" \
        --filter_samples "${SAMPLE_NAME}" \
        --regions "${GRAPH_REGION}" \
        --min_af "${MIN_AF}" \
        --primary \
        --gcsa_index \
        --xg_index \
        --gbwt_index \
        --snarls_index
fi

READS_DIR="${GRAPHS_PATH}/sim-${READ_SEED}-${READ_COUNT}-${READ_CHUNKS}"

if [[ ! -e "${READS_DIR}" ]]; then 
    # Now we need to simulate reads from the two haplotypes
    # This will make a "sim.gam"
    toil-vg sim "${TREE_PATH}/sim" \
        "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo_thread_0.xg" \
        "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo_thread_1.xg" \
        "${READ_COUNT}" \
        "${READS_DIR}" \
        --config "${TREE_PATH}/toil-vg.conf" \
        --annotate_xg "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo.xg" \
        --gam \
        --fastq_out \
        --seed "${READ_SEED}" \
        --sim_chunks "${READ_CHUNKS}" \
        --fastq "${TRAINING_FASTQ}"
fi

# Make sure we have the SLLS linear index file
# We will work on the filtered VCF
FILTERED_VCF_BASENAME="${VCF_BASENAME%.vcf.gz}_filter.vcf.gz"
SLLS_INDEX="${GRAPHS_PATH}/slls/${FILTERED_VCF_BASENAME}.slls"
if [[ ! -e "${SLLS_INDEX}" ]]; then
    # We need to make the SLLS index
    mkdir "${GRAPHS_PATH}/slls"
    cp "${GRAPHS_PATH}/${FILTERED_VCF_BASENAME}" "${GRAPHS_PATH}/slls/${FILTERED_VCF_BASENAME}"
    cd deps/sublinear-Li-Stephens && make && cd ../..
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`/deps/sublinear-Li-Stephens/deps/htslib/ ./deps/sublinear-Li-Stephens/bin/serializer "${GRAPHS_PATH}/slls/${FILTERED_VCF_BASENAME}"
fi

 # Now we do a bunch of stuff in parallel
JOB_ARRAY=()

# What do they return?
JOB_RETURNS=()

# What condition names have we run
CONDITIONS=()

# How many jobs should we let run at once
MAX_JOBS=2

# We have a function to wait for the parallel jobs to finish, if MAX_JOBS or
# more are running
function wait_on_jobs() {
    # How many jobs are running?
    CURRENT_JOBS="${#JOB_ARRAY[@]}"

    if [[ "${CURRENT_JOBS}" -lt "${MAX_JOBS}" ]]; then
        # If we haven't hit the cap, don't do anything.
        return
    fi

    # Otherwise we have to collect some jobs
    COLLECTED_JOBS=0

    # Now wait for all the jobs and fail if any failed
    for JOB in "${JOB_ARRAY[@]}"; do
        if [[ -z "${JOB}" ]]; then
            # Drop empty strings that get in here. This happens if we forget
            # the trailing & on something intended to be in parallel.
            continue
        fi
        wait "${JOB}"
        RETURN_CODE="$?"
        
        if [[ "${RETURN_CODE}" != "0" ]]; then
            # A job has failed.
            # Collect up all the jobs now, actually.
            echo "Job PID ${JOB} failed with return code ${RETURN_CODE}; flushing queue" 1>&2
            MAX_JOBS=0
        fi
        
        JOB_RETURNS+=("${RETURN_CODE}")
        ((COLLECTED_JOBS+=1))
        
        if [[ "$((CURRENT_JOBS-COLLECTED_JOBS))" -lt "${MAX_JOBS}" ]]; then
            # No need to clean up any more jobs
            break
        fi
    done

    JOB_NUMBER=1
    for JOB_RETURN in "${JOB_RETURNS[@]}"; do
        echo "Job ${JOB_NUMBER} exit status: ${JOB_RETURN}"
        ((JOB_NUMBER+=1))
        if [[ "${JOB_RETURN}" != "0" ]]; then
            echo "Job failed!" 1>&2
            exit 1
        fi
    done
    
    # Pop off the finished jobs
    JOB_ARRAY=("${JOB_ARRAY[@]:${COLLECTED_JOBS}}")
    # Delete all the return codes
    JOB_RETURNS=()
}

if [[ "${RUN_JOBS}" == "1" ]]; then
    # We actually want to run the toil-vg jobs

    #CONDITIONS+=("snp1kg-mp")
    #if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp" ]]; then
    #    # Do the full snp1kg graph multipath
    #    toil-vg mapeval "${TREE_PATH}/snp1kg-mp" "${OUTPUT_PATH}/snp1kg-mp" \
    #        --single_reads_chunk \
    #        --config "${TREE_PATH}/toil-vg.conf" \
    #        --maxDisk 100G \
    #        --multipath-only \
    #        --fastq "${READS_DIR}/sim.fq.gz" \
    #        --truth "${READS_DIR}/true.pos" \
    #        --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter" \
    #        --gam-names snp1kg 2>&1 & 
    #    JOB_ARRAY+=("$!")
    #fi
    #wait_on_jobs
    
    CONDITIONS+=("snp1kg-mp-snarlcut")
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-snarlcut" ]]; then
        # Do the full snp1kg graph multipath
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp-snarlcut" "${OUTPUT_PATH}/snp1kg-mp-snarlcut" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --use-snarls \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter" \
            --gam-names snp1kg-snarlcut 2>&1 & 
        JOB_ARRAY+=("$!")
    fi
    wait_on_jobs
    
    #CONDITIONS+=("snp1kg-mp-gbwt")
    #if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-gbwt" ]]; then
    #    # Do the full snp1kg graph multipath with gbwt
    #    toil-vg mapeval "${TREE_PATH}/snp1kg-mp-gbwt" "${OUTPUT_PATH}/snp1kg-mp-gbwt" \
    #        --single_reads_chunk \
    #        --config "${TREE_PATH}/toil-vg.conf" \
    #        --maxDisk 100G \
    #        --multipath-only \
    #        --use-gbwt \
    #        --fastq "${READS_DIR}/sim.fq.gz" \
    #        --truth "${READS_DIR}/true.pos" \
    #        --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter" \
    #        --gam-names snp1kg-gbwt 2>&1 &
    #    JOB_ARRAY+=("$!")
    #fi
    #wait_on_jobs
    
    #CONDITIONS+=("snp1kg-mp-gbwt-traceback")
    #if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-gbwt-traceback" ]]; then
    #    # Do the full snp1kg graph multipath with gbwt
    #    toil-vg mapeval "${TREE_PATH}/snp1kg-mp-gbwt-traceback" "${OUTPUT_PATH}/snp1kg-mp-gbwt-traceback" \
    #        --single_reads_chunk \
    #        --config "${TREE_PATH}/toil-vg.conf" \
    #        --maxDisk 100G \
    #        --multipath-only \
    #        --use-gbwt \
    #        --mpmap_opts "--max-paths 10" \
    #        --fastq "${READS_DIR}/sim.fq.gz" \
    #        --truth "${READS_DIR}/true.pos" \
    #        --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter" \
    #        --gam-names snp1kg-gbwt-traceback 2>&1 &
    #    JOB_ARRAY+=("$!")
    #fi
    #wait_on_jobs
    
    CONDITIONS+=("snp1kg-mp-gbwt-snarlcut")
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-gbwt-snarlcut" ]]; then
        # Do the full snp1kg graph multipath with snarl cutting and gbwt but single path traceback
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp-gbwt-snarlcut" "${OUTPUT_PATH}/snp1kg-mp-gbwt-snarlcut" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --use-gbwt \
            --use-snarls \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter" \
            --gam-names snp1kg-gbwt-snarlcut 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    wait_on_jobs
    
    CONDITIONS+=("snp1kg-mp-gbwt-traceback-snarlcut")
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-gbwt-traceback-snarlcut" ]]; then
        # Do the full snp1kg graph multipath with gbwt
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp-gbwt-traceback-snarlcut" "${OUTPUT_PATH}/snp1kg-mp-gbwt-traceback-snarlcut" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --use-gbwt \
            --use-snarls \
            --mpmap_opts "--max-paths 10" \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter" \
            --gam-names snp1kg-gbwt-traceback-snarlcut 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    wait_on_jobs
    
    # This should replicate our best known performance...
    CONDITIONS+=("snp1kg-mp-gbwt-traceback-snarlcut-unfiltered")
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-gbwt-traceback-snarlcut-unfiltered" ]]; then
        # Do the full snp1kg graph multipath with gbwt
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp-gbwt-traceback-snarlcut-unfiltered" "${OUTPUT_PATH}/snp1kg-mp-gbwt-traceback-snarlcut-unfiltered" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --use-gbwt \
            --use-snarls \
            --mpmap_opts "--max-paths 10" \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}" \
            --gam-names snp1kg-gbwt-traceback-snarlcut-unfiltered 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    wait_on_jobs
    
    CONDITIONS+=("snp1kg-mp-gbwt-traceback")
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-gbwt-traceback" ]]; then
        # Do the full snp1kg graph multipath with gbwt
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp-gbwt-traceback" "${OUTPUT_PATH}/snp1kg-mp-gbwt-traceback" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --use-gbwt \
            --use-snarls \
            --mpmap_opts "--max-paths 10" \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter" \
            --gam-names snp1kg-gbwt-traceback 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    wait_on_jobs
    
    CONDITIONS+=("snp1kg-mp-minaf-snarlcut")
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-minaf-snarlcut" ]]; then
        # And with the min allele frequency
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp-minaf-snarlcut" "${OUTPUT_PATH}/snp1kg-mp-minaf-snarlcut" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --use-snarls \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_minaf_${MIN_AF}" \
            --gam-names snp1kg-minaf-snarlcut 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    wait_on_jobs
    
    CONDITIONS+=("snp1kg-mp-minaf")
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-minaf" ]]; then
        # And with the min allele frequency
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp-minaf" "${OUTPUT_PATH}/snp1kg-mp-minaf" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --use-snarls \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_minaf_${MIN_AF}" \
            --gam-names snp1kg-minaf 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    wait_on_jobs
    
    CONDITIONS+=("snp1kg-mp-positive-snarlcut")
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-positive-snarlcut" ]]; then    
        # And the positive control with only real variants
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp-positive-snarlcut" "${OUTPUT_PATH}/snp1kg-mp-positive-snarlcut" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --use-snarls \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}" \
            --gam-names snp1kg-positive-snarlcut 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    wait_on_jobs
    
    CONDITIONS+=("snp1kg-mp-positive")
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-positive" ]]; then    
        # And the positive control with only real variants
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp-positive" "${OUTPUT_PATH}/snp1kg-mp-positive" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --use-snarls \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}" \
            --gam-names snp1kg-positive 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    wait_on_jobs
    
    CONDITIONS+=("primary-mp-snarlcut")
    if [[ ! -e "${OUTPUT_PATH}/primary-mp-snarlcut" ]]; then
        # And the primary path only
        toil-vg mapeval "${TREE_PATH}/primary-mp-snarlcut" "${OUTPUT_PATH}/primary-mp-snarlcut" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --use-snarls \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/primary" \
            --gam-names primary-snarlcut 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    wait_on_jobs
    
    CONDITIONS+=("primary-mp")
    if [[ ! -e "${OUTPUT_PATH}/primary-mp" ]]; then
        # And the primary path only
        toil-vg mapeval "${TREE_PATH}/primary-mp" "${OUTPUT_PATH}/primary-mp" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --use-snarls \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/primary" \
            --gam-names primary 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    wait_on_jobs
    
    #CONDITIONS+=("snp1kg-mp-slls")
    #if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-slls/position.results.tsv" ]]; then
    #    # This one's a bit different since we need to manually do the mapping
    #    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-slls/aligned-snp1kg-slls-pe_default.gam" ]]; then
    #        mkdir -p "${OUTPUT_PATH}/snp1kg-mp-slls"
    #        vg mpmap --linear-index "${SLLS_INDEX}" \
    #            --linear-path "${GRAPH_CONTIG}" \
    #            -x "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter.xg" \
    #            -g "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter.gcsa" \
    #            --fastq "${READS_DIR}/sim.fq.gz" \
    #            -i \
    #            -S \
    #            -t 32 \
    #            >"${OUTPUT_PATH}/snp1kg-mp-slls/aligned-snp1kg-slls-pe_default.gam"
    #    fi
    #    # Then do the mapeval
    #    toil-vg mapeval "${TREE_PATH}/snp1kg-mp-slls" "${OUTPUT_PATH}/snp1kg-mp-slls" \
    #        --gams "${OUTPUT_PATH}/snp1kg-mp-slls/aligned-snp1kg-slls-pe_default.gam" \
    #        --config "${TREE_PATH}/toil-vg.conf" \
    #        --maxDisk 100G \
    #        --truth "${READS_DIR}/true.pos" \
    #        --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter" \
    #        --gam-names snp1kg-slls-mp-pe 2>&1 &
    #    JOB_ARRAY+=("$!")
    #fi
    #wait_on_jobs
    
    #CONDITIONS+=("snp1kg-mp-slls-snarlcut")
    #if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-slls-snarlcut/position.results.tsv" ]]; then
    #    # This one's a bit different since we need to manually do the mapping
    #    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-slls-snarlcut/aligned-snp1kg-slls-snarlcut-pe_default.gam" ]]; then
    #        mkdir -p "${OUTPUT_PATH}/snp1kg-mp-slls-snarlcut"
    #        vg mpmap --linear-index "${SLLS_INDEX}" \
    #            --linear-path "${GRAPH_CONTIG}" \
    #            -x "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter.xg" \
    #            -g "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter.gcsa" \
    #            --snarls "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter.snarls" \
    #            --fastq "${READS_DIR}/sim.fq.gz" \
    #            -i \
    #            -S \
    #            -t 32 \
    #            >"${OUTPUT_PATH}/snp1kg-mp-slls-snarlcut/aligned-snp1kg-slls-snarlcut-pe_default.gam"
    #    fi
    #    # Then do the mapeval
    #    toil-vg mapeval "${TREE_PATH}/snp1kg-mp-slls-snarlcut" "${OUTPUT_PATH}/snp1kg-mp-slls-snarlcut" \
    #        --gams "${OUTPUT_PATH}/snp1kg-mp-slls-snarlcut/aligned-snp1kg-slls-snarlcut-pe_default.gam" \
    #        --config "${TREE_PATH}/toil-vg.conf" \
    #        --maxDisk 100G \
    #        --truth "${READS_DIR}/true.pos" \
    #        --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_filter" \
    #        --gam-names snp1kg-slls-snarlcut-mp-pe 2>&1 &
    #    JOB_ARRAY+=("$!")
    #fi
    #wait_on_jobs

fi

# Make all the jobs finish
MAX_JOBS=0
wait_on_jobs

if [[ ! -e "${OUTPUT_PATH}/position.results.tsv" ]]; then
    # Concatenate all the conditions' position results files
    
    FIRST_CONDITION=1
    for CONDITION in "${CONDITIONS[@]}"; do
        if [[ "${FIRST_CONDITION}" == "1" ]]; then
            # Keep the header on the first condition we process
            cat "${OUTPUT_PATH}/${CONDITION}/position.results.tsv" > "${OUTPUT_PATH}/position.results.tsv"
        else
            # Drop the header
            cat "${OUTPUT_PATH}/${CONDITION}/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
        fi
        FIRST_CONDITION=0
    done
fi

if [[ ! -e "${OUTPUT_PATH}/qq" ]]; then
    # Collect all the qq plots to one directory to page through easily
    mkdir "${OUTPUT_PATH}/qq"
    
    CONDITION_NUMBER=0
    for CONDITION in "${CONDITIONS[@]}"; do
        cp "${OUTPUT_PATH}/${CONDITION}/plot-qq.svg" "${OUTPUT_PATH}/qq/qq-${CONDITION_NUMBER}-${CONDITION}.svg"
        ((CONDITION_NUMBER+=1))
    done

fi

if [[ ! -e "${OUTPUT_PATH}/table.tsv" ]]; then
    # Make a table of wrong reads

    # First we need a baseline for comparing against
    BASELINE_CONDITION="snp1kg-mp-snarlcut"

    # Make header for table of wrong read counts
    printf "Condition\tWrong reads total\tAt MAPQ 60\tAt MAPQ 0\tAt MAPQ >0\tNew vs. ${BASELINE_CONDITION}\tFixed vs. ${BASELINE_CONDITION}\tAvg. Correct MAPQ\tCorrect MAPQ 0\n" > "${OUTPUT_PATH}/table.tsv"

    # Pull out the baseline wrong reads
    cat "${OUTPUT_PATH}/${BASELINE_CONDITION}/position.results.tsv" | sed 1d | grep -- "-pe" | grep -v "^1" | cut -f4 | sort > "${OUTPUT_PATH}/baseline-wrong-names.tsv"

    for CONDITION in "${CONDITIONS[@]}"; do
        
        # We want a table like
        # Condition 	Wrong reads total 	At MAPQ 60 	At MAPQ 0 	At MAPQ >0 	New vs. baseline 	Fixed vs. baseline
        printf "${CONDITION}\t" >> "${OUTPUT_PATH}/table.tsv"
        
        
        # Get the wrong reads for this condition
        cat "${OUTPUT_PATH}/${CONDITION}/position.results.tsv" | sed 1d | grep -- "-pe" | grep -v "^1" > "${OUTPUT_PATH}/${CONDITION}/wrong.tsv"
        
        # Wrong total
        cat "${OUTPUT_PATH}/${CONDITION}/wrong.tsv" | wc -l | tr -d '\n' >> "${OUTPUT_PATH}/table.tsv"
        printf "\t" >> "${OUTPUT_PATH}/table.tsv"
        
        # Wrong at MAPQ 60
        cat "${OUTPUT_PATH}/${CONDITION}/wrong.tsv" | grep -P "\t60\t" | wc -l | tr -d '\n' >> "${OUTPUT_PATH}/table.tsv"
        printf "\t" >> "${OUTPUT_PATH}/table.tsv"
        
        # Wrong at MAPQ 0
        cat "${OUTPUT_PATH}/${CONDITION}/wrong.tsv" | grep -P "\t0\t" | wc -l | tr -d '\n' >> "${OUTPUT_PATH}/table.tsv"
        printf "\t" >> "${OUTPUT_PATH}/table.tsv"
        
        # Wrong at MAPQ >0
        cat "${OUTPUT_PATH}/${CONDITION}/wrong.tsv" | grep -v -P "\t0\t" | wc -l | tr -d '\n' >> "${OUTPUT_PATH}/table.tsv"
        printf "\t" >> "${OUTPUT_PATH}/table.tsv"
        
        # Get the wrong read names for this condition
        cat "${OUTPUT_PATH}/${CONDITION}/wrong.tsv" | cut -f4 | sort > "${OUTPUT_PATH}/${CONDITION}/wrong-names.tsv"
       
        # Count newly wrong names (not in file 1 or both)
        comm -1 -3 "${OUTPUT_PATH}/baseline-wrong-names.tsv" "${OUTPUT_PATH}/${CONDITION}/wrong-names.tsv" | wc -l | tr -d '\n' >> "${OUTPUT_PATH}/table.tsv"
        printf "\t" >> "${OUTPUT_PATH}/table.tsv"
        
        # Count newly right names (not in file 2 or both)
        comm -2 -3 "${OUTPUT_PATH}/baseline-wrong-names.tsv" "${OUTPUT_PATH}/${CONDITION}/wrong-names.tsv" | wc -l | tr -d '\n' >> "${OUTPUT_PATH}/table.tsv"
        printf "\t" >> "${OUTPUT_PATH}/table.tsv"
        
        # Get the right reads
        cat "${OUTPUT_PATH}/${CONDITION}/position.results.tsv" | sed 1d | grep -- "-pe" | grep "^1" > "${OUTPUT_PATH}/${CONDITION}/right.tsv"
        
        # Compute average MAPQ for correct reads
        # See <https://stackoverflow.com/a/19149931>
        cat "${OUTPUT_PATH}/${CONDITION}/right.tsv" | cut -f2 | awk '{total += $1} END {if (NR > 0) {print total / NR} else {print "N/A"}}' | tr -d '\n' >> "${OUTPUT_PATH}/table.tsv"
        printf "\t" >> "${OUTPUT_PATH}/table.tsv"
        
        # Count correct reads at MAPQ 0
        cat "${OUTPUT_PATH}/${CONDITION}/right.tsv" | grep -P "\t0\t" | wc -l | tr -d '\n' >> "${OUTPUT_PATH}/table.tsv"
        printf "\n" >> "${OUTPUT_PATH}/table.tsv"
    done
    
    rm "${OUTPUT_PATH}/baseline-wrong-names.tsv"

fi

# Determine our source directory, where the ROC plotting script also lives
# See <https://stackoverflow.com/a/246128>
SCRIPT_DIRECTORY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Do the R plots

if [[ ! -e "${OUTPUT_PATH}/roc.svg" ]]; then
    Rscript "${SCRIPT_DIRECTORY}/plot-roc.R" "${OUTPUT_PATH}/position.results.tsv" "${OUTPUT_PATH}/roc.svg"
fi
if [[ ! -e "${OUTPUT_PATH}/pr.svg" ]]; then
    Rscript "${SCRIPT_DIRECTORY}/plot-pr.R" "${OUTPUT_PATH}/position.results.tsv" "${OUTPUT_PATH}/pr.svg"
fi

rm "${TREE_PATH}/toil-vg.conf"
rmdir "${TREE_PATH}"
echo "Results available as plots ${OUTPUT_PATH}/roc.svg and ${OUTPUT_PATH}/pr.svg and table ${OUTPUT_PATH}/table.tsv"






























