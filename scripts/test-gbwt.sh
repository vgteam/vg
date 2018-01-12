#!/usr/bin/env bash
# test-gbwt.sh: Plot the effect of haplotype information on mapping performance

# We are going to compare 5 mapping regimes:
# 1. The snp1kg graph without the variants in the sample (negative control)
# 2. The full snp1kg graph (Heng Effect positive control)
# 3. The full snp1kg graph with GBWT haplotype information (under test)
# 4. The frequency-filtered snp1kg graph (current best)
# 5. The primary path graph (Heng Effect negative control)

# We want to know if the GBWT reduces the Heng Effect (poor mapping due to
# spurious alignments to rare variants)

# OK so on to the actual code.
set -e

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
    
# Actually do a smaller test
REGION_NAME="BRCA1"
GRAPH_REGION="17:43044293-43125482"
FASTA_BASENAME="chr17.fa.gz"
VCF_BASENAME="1kg_hg38-BRCA1.vcf.gz"

# Define the sample to use for synthesizing reads
SAMPLE_NAME="HG00096"


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
    toil-vg construct "${TREE_PATH}" "${GRAPHS_PATH}" \
        --vcf "${SOURCE_BASE_URL}/${VCF_BASENAME}" \
        --fasta "${SOURCE_BASE_URL}/${FASTA_BASENAME}" \
        --out_name "snp1kg-${REGION_NAME}" \
        --alt_paths \
        --realTimeLogging \
        --control_sample "${SAMPLE_NAME}" \
        --haplo_sample "${SAMPLE_NAME}" \
        --regions "${GRAPH_REGION}" \
        --min_af 0.0335570469 \
        --primary \
        --gcsa_index \
        --xg_index \
        --gbwt_index
        
fi

exit

# Now we run the actual mappings for the experiment

# Where do we get our simulated reads for our sample? Last time I got them by running the Jenkins tests.

# Do the full snp1kg graph with GBWT
toil-vg mapeval "${TREE_PATH}" "${OUTPUT_PATH}/snp1kg+gbwt" \
    --use-gbwt \
    --map_opts "--hap-exp 1" \
    --gam_input_reads sim-mhc.gam \
    --gam-input-xg index-gpbwt-mhc.xg \
    --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}" \
    --gam-names snp1kg



























