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
READ_SEED="90"
# And a read count
READ_COUNT="200"
# Chunks to simulate in (which affects results)
READ_CHUNKS="32"

MODE="21"
if [[ "${MODE}" == "mhc" ]]; then
    # Actually do a smaller test
    READ_COUNT="100000"
    REGION_NAME="MHC"
    GRAPH_REGION="6:28510119-33480577"
    FASTA_BASENAME="chr6.fa.gz"
    VCF_BASENAME="1kg_hg38-MHC.vcf.gz"
elif [[ "${MODE}" == "tiny" ]]; then
    # Do just 20 kb of MHC and a very few reads
    READ_COUNT="1000"
    REGION_NAME="MHC"
    GRAPH_REGION="6:28510119-28520119"
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

 # Now we do a bunch of stuff in parallel
JOB_ARRAY=()

# What do they return?
JOB_RETURNS=()

# We have a function to wait for the parallel jobs to finish
function wait_on_jobs() {
    # Now wait for all the jobs and fail if any failed
    for JOB in "${JOB_ARRAY[@]}"; do
        if [[ -z "${JOB}" ]]; then
            # Drop empty strings that get in here. This happens if we forget
            # the trailing & on something intended to be in parallel.
            continue
        fi
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
    
    JOB_ARRAY=()
    JOB_RETURNS=()
}

if [[ "${RUN_JOBS}" == "1" ]]; then
    # We actually want to run the toil-vg jobs

    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp" ]]; then
        # Do the full snp1kg graph multipath
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp" "${OUTPUT_PATH}/snp1kg-mp" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}" \
            --gam-names snp1kg 2>&1 & 
        JOB_ARRAY+=("$!")
    fi
    
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-gbwt" ]]; then
        # Do the full snp1kg graph multipath with gbwt
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp-gbwt" "${OUTPUT_PATH}/snp1kg-mp-gbwt" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --use-gbwt \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}" \
            --gam-names snp1kg-gbwt 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    
    wait_on_jobs
    
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-gbwt-traceback" ]]; then
        # Do the full snp1kg graph multipath with gbwt
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp-gbwt-traceback" "${OUTPUT_PATH}/snp1kg-mp-gbwt-traceback" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --use-gbwt \
            --mpmap_opts "--max-paths 10" \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}" \
            --gam-names snp1kg-gbwt-traceback 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    
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
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}" \
            --gam-names snp1kg-gbwt-traceback 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    
    wait_on_jobs
    
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-minaf" ]]; then
        # And with the min allele frequency
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp-minaf" "${OUTPUT_PATH}/snp1kg-mp-minaf" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_minaf_${MIN_AF}" \
            --gam-names snp1kg-minaf 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    
    if [[ ! -e "${OUTPUT_PATH}/snp1kg-mp-positive" ]]; then    
        # And the positive control with only real variants
        toil-vg mapeval "${TREE_PATH}/snp1kg-mp-positive" "${OUTPUT_PATH}/snp1kg-mp-positive" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/snp1kg-${REGION_NAME}_${SAMPLE_NAME}" \
            --gam-names snp1kg-positive 2>&1 &
        JOB_ARRAY+=("$!")
    fi
    
    wait_on_jobs
        
    if [[ ! -e "${OUTPUT_PATH}/primary-mp" ]]; then
        # And the primary path only
        toil-vg mapeval "${TREE_PATH}/primary-mp" "${OUTPUT_PATH}/primary-mp" \
            --single_reads_chunk \
            --config "${TREE_PATH}/toil-vg.conf" \
            --maxDisk 100G \
            --multipath-only \
            --fastq "${READS_DIR}/sim.fq.gz" \
            --truth "${READS_DIR}/true.pos" \
            --index-bases "${GRAPHS_PATH}/primary" \
            --gam-names primary 2>&1 &
        JOB_ARRAY+=("$!")
    fi

    wait_on_jobs

fi

if [[ ! -e "${OUTPUT_PATH}/position.results.tsv" ]]; then

    # Combine all the position.results.tsv files into one
    #cat "${OUTPUT_PATH}/snp1kg/position.results.tsv" > "${OUTPUT_PATH}/position.results.tsv"
    #cat "${OUTPUT_PATH}/snp1kg-gbwt/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
    cat "${OUTPUT_PATH}/snp1kg-mp/position.results.tsv" > "${OUTPUT_PATH}/position.results.tsv"
    cat "${OUTPUT_PATH}/snp1kg-mp-gbwt/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
    cat "${OUTPUT_PATH}/snp1kg-mp-gbwt-traceback/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
    cat "${OUTPUT_PATH}/snp1kg-mp-gbwt-traceback-snarlcut/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
    #cat "${OUTPUT_PATH}/snp1kg-fullgbwt/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
    #cat "${OUTPUT_PATH}/snp1kg-gbwt-slight/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
    cat "${OUTPUT_PATH}/snp1kg-mp-minaf/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
    #cat "${OUTPUT_PATH}/snp1kg-negative/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
    cat "${OUTPUT_PATH}/snp1kg-mp-positive/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
    cat "${OUTPUT_PATH}/primary-mp/position.results.tsv" | sed 1d >> "${OUTPUT_PATH}/position.results.tsv"
    
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

if [[ ! -e "${OUTPUT_PATH}/table.tsv" ]]; then

    # Generate a table of wrong read counts
    printf "Condition\tWrong reads total\tAt MAPQ 60\tAt MAPQ 0\tAt MAPQ >0\tNew vs. mpmap\tFixed vs. mpmap\n" > "${OUTPUT_PATH}/table.tsv"

    # First we need a baseline of snp1kg-mp for comparing against

    cat "${OUTPUT_PATH}/snp1kg-mp/position.results.tsv" | sed 1d | grep -- "-pe" | grep -v "^1" | cut -f4 | sort > "${OUTPUT_PATH}/baseline-wrong-names.tsv"

    for CONDITION in snp1kg-mp snp1kg-mp-gbwt snp1kg-mp-gbwt-traceback snp1kg-mp-gbwt-traceback-snarlcut snp1kg-mp-minaf snp1kg-mp-positive primary-mp; do
        
        # We want a table like
        # Condition 	Wrong reads total 	At MAPQ 60 	At MAPQ 0 	At MAPQ >0 	New vs. mpmap 	Fixed vs. mpmap
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
        printf "\n" >> "${OUTPUT_PATH}/table.tsv"
    done
    
    rm "${OUTPUT_PATH}/baseline-wrong-names.tsv"

fi

rm "${TREE_PATH}/toil-vg.conf"
rmdir "${TREE_PATH}"
echo "Results available as plots ${OUTPUT_PATH}/roc.svg and ${OUTPUT_PATH}/pr.svg and table ${OUTPUT_PATH}/table.tsv"






























