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
# Define the region to build the graph on, as contig[:start-end]
GRAPH_REGION="21"

# Set a FASTQ to model reads after
TRAINING_FASTQ="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_005_BH814YADXX/Project_RM8398/Sample_U5a/U5a_AGTCAA_L002_R1_007.fastq.gz"
# And a read simulation seed
READ_SEED="75"
# And a read count
READ_COUNT="3000000"
    
# Actually do a smaller test
#REGION_NAME="MHC"
#GRAPH_REGION="6:28510119-33480577"
#FASTA_BASENAME="chr6.fa.gz"
#VCF_BASENAME="1kg_hg38-MHC.vcf.gz"

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
        --vcf "${SOURCE_BASE_URL}/${VCF_BASENAME}" \
        --fasta "${SOURCE_BASE_URL}/${FASTA_BASENAME}" \
        --out_name "snp1kg-${REGION_NAME}" \
        --alt_paths \
        --realTimeLogging \
        --control_sample "${SAMPLE_NAME}" \
        --haplo_sample "${SAMPLE_NAME}" \
        --regions "${GRAPH_REGION}" \
        --min_af "${MIN_AF}" \
        --primary \
        --gcsa_index \
        --xg_index \
        --gbwt_index
fi

READS_DIR="${GRAPHS_PATH}/sim-${READ_SEED}-${READ_COUNT}"

if [[ ! -e "${READS_DIR}" ]]; then 
    # Now we need to simulate reads from the two haplotypes
    # This will make a "sim.gam"
    toil-vg sim "${TREE_PATH}/sim" \
        "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo_thread_0.xg" \
        "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo_thread_1.xg" \
        "${READ_COUNT}" \
        "${READS_DIR}" \
        --annotate_xg "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo.xg" \
        --gam \
        --seed "${READ_SEED}" \
        --fastq "${TRAINING_FASTQ}"
fi

if [[ "${RUN_JOBS}" == "1" ]]; then
    # We actually want to run the toil-vg jobs

    # Now we do a bunch of stuff in parallel
    JOB_ARRAY=()
    
    # What do they return?
    JOB_RETURNS=()

    if [[ ! -e "${OUTPUT_PATH}/snp1kg" ]]; then
        # Do the full snp1kg graph without the GBWT
        toil-vg mapeval "${TREE_PATH}/snp1kg" "${OUTPUT_PATH}/snp1kg" \
            --gam_input_reads "${READS_DIR}/sim.gam" \
            --gam-input-xg "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo.xg" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}" \
            --gam-names snp1kg 2>&1 &
        JOB_ARRAY+=("$!")
    fi

    if [[ ! -e "${OUTPUT_PATH}/snp1kg-gbwt" ]]; then
        # Do the full snp1kg graph with GBWT
        toil-vg mapeval "${TREE_PATH}/snp1kg-gbwt" "${OUTPUT_PATH}/snp1kg-gbwt" \
            --use-gbwt \
            --map_opts "--hap-exp 1" \
            --gam_input_reads "${READS_DIR}/sim.gam" \
            --gam-input-xg "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo.xg" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}" \
            --gam-names snp1kg-gbwt 2>&1 &
        JOB_ARRAY+=("$!")
    fi

     if [[ ! -e "${OUTPUT_PATH}/snp1kg-fullgbwt" ]]; then
        # Do the full snp1kg graph with GBWT
        toil-vg mapeval "${TREE_PATH}/snp1kg-fullgbwt" "${OUTPUT_PATH}/snp1kg-fullgbwt" \
            --use-gbwt \
            --map_opts "--hap-exp 1" \
            --gam_input_reads "${READS_DIR}/sim.gam" \
            --gam-input-xg "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo.xg" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_fullgbwt" \
            --gam-names snp1kg-fullgbwt 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    
    #if [[ ! -e "${OUTPUT_PATH}/snp1kg-gbwt-slight" ]]; then
    #    # Do the full snp1kg graph with GBWT
    #    toil-vg mapeval "${TREE_PATH}/snp1kg-gbwt" "${OUTPUT_PATH}/snp1kg-gbwt-slight" \
    #        --use-gbwt \
    #        --map_opts "--hap-exp 0.5" \
    #        --gam_input_reads "${READS_DIR}/sim.gam" \
    #        --gam-input-xg "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo.xg" \
    #        --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}" \
    #        --gam-names snp1kg-gbwt-slight 2>&1 &
    #    JOB_ARRAY+=("$!")
    #fi
        
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-minaf" ]]; then
        # And with the min allele frequency
        toil-vg mapeval "${TREE_PATH}/snp1kg-minaf" "${OUTPUT_PATH}/snp1kg-minaf" \
            --gam_input_reads "${READS_DIR}/sim.gam" \
            --gam-input-xg "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo.xg" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_minaf_${MIN_AF}" \
            --gam-names snp1kg-minaf 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-negative" ]]; then    
        # And the negative control with correct variants removed
        toil-vg mapeval "${TREE_PATH}/snp1kg-negative" "${OUTPUT_PATH}/snp1kg-negative" \
            --gam_input_reads "${READS_DIR}/sim.gam" \
            --gam-input-xg "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo.xg" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_minus_${SAMPLE_NAME}" \
            --gam-names snp1kg-negative 2>&1 &
        JOB_ARRAY+=("$!")
    fi
        
    if [[ ! -e "${OUTPUT_PATH}/primary" ]]; then
        # And the primary path only
        toil-vg mapeval "${TREE_PATH}/primary" "${OUTPUT_PATH}/primary" \
            --gam_input_reads "${READS_DIR}/sim.gam" \
            --gam-input-xg "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}_haplo.xg" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_primary" \
            --gam-names primary 2>&1 &
        JOB_ARRAY+=("$!")
    fi

    # Now wait for all the jobs and fail if any failed
    for JOB in "${JOB_ARRAY[@]}"; do
        wait "${JOB}"
        JOB_RETURNS+=("$?")
    done

    JOB_NUMBER=1
    for JOB_RETURN in "${JOB_RETURNS[@]}"; do
        echo "Job ${JOB_NUMBER} exit status: ${JOB_RETURN}"
        ((JOB_NUMBER=JOB_NUMBER+1))
        if [[ "${JOB_RETURN}" != "0" ]]; then
            echo "Job failed!" 1>&2
            exit 1
        fi
    done

fi
    
# Combine all the position.results.tsv files into one
cat "${OUTPUT_PATH}/snp1kg/position.results.tsv" > "${OUTPUT_PATH}/position.results.tsv"
cat "${OUTPUT_PATH}/snp1kg-gbwt/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
cat "${OUTPUT_PATH}/snp1kg-fullgbwt/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
#cat "${OUTPUT_PATH}/snp1kg-gbwt-slight/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
cat "${OUTPUT_PATH}/snp1kg-minaf/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
cat "${OUTPUT_PATH}/snp1kg-negative/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
cat "${OUTPUT_PATH}/primary/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"

# Determine our source directory, where the ROC plotting script also lives
# See <https://stackoverflow.com/a/246128>
SCRIPT_DIRECTORY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Do the R plot
Rscript "${SCRIPT_DIRECTORY}/plot-roc.R" "${OUTPUT_PATH}/position.results.tsv" "${OUTPUT_PATH}/roc.svg"

rmdir "${TREE_PATH}"
echo "Successfully produced ROC plot ${OUTPUT_PATH}/roc.svg"




























