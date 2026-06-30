#!/usr/bin/env bash

# Usage: ./script.sh <fasta_path> <model_name> <model_link> <base_output_dir> <repr_layers>

# (0) Argument handling
FASTA=$1
MODEL_NAME=$2
MODEL_LINK=$3
BASE_OUT=$4
REPR_LAYERS=$5

if [[ -z "$FASTA" || -z "$MODEL_NAME" || -z "$MODEL_LINK" || -z "$BASE_OUT" ]]; then
    echo "Usage: $0 <fasta_path> <model_name> <model_link> <output_dir> [repr_layers]" >&2
    exit 1
fi

# Path logic
FILENAME="${FASTA##*/}"
PROTEIN_NAME="${FILENAME%.*}"
PROT_OUT_DIR="${BASE_OUT}/${PROTEIN_NAME}"
EMB_DIR="${PROT_OUT_DIR}/${MODEL_NAME}_embeddings"

echo "--------------------------------------------------"
echo "Processing Protein: ${PROTEIN_NAME}"
echo "Model: ${MODEL_NAME}"
echo "--------------------------------------------------"

mkdir -p "${PROT_OUT_DIR}"

### (1) Saturation Mutagenesis
echo "Step 1: Running Saturation Mutagenesis..."
MUT_DIR="${PROT_OUT_DIR}/mutagenesis"
mkdir -p "${MUT_DIR}"

python scripts/saturation_mutagenesis.py \
    --fasta "${FASTA}" \
    --output "${MUT_DIR}/${PROTEIN_NAME}" 2>&2

if [ $? -ne 0 ]; then echo "Error in Mutagenesis" >&2; exit 1; fi

### (2) Generate Embeddings
echo "Step 2: Generating Embeddings..."
mkdir -p "${EMB_DIR}"

tasks=(
    "WT|${FASTA}|wt"
    "SM|${MUT_DIR}/${PROTEIN_NAME}_SM.fa|sm"
    "Trunc|${MUT_DIR}/${PROTEIN_NAME}_trunc.fa|trunc"
)

for task in "${tasks[@]}"; do
    IFS="|" read -r label infile suffix <<< "$task"
    echo "  -> Embedding ${label}..."
    
    # 1. ESMC (EvolutionaryScale) Logic
    if [[ "$MODEL_NAME" == *"ESMC"* ]]; then
        TASK_OUT="${EMB_DIR}/${suffix}"
        mkdir -p "${TASK_OUT}"
        
        # Using the new compatibility script we adapted
        python /home/laura.cano/projects/plms/esm/esm-scripts/extract_esmc.py \
            "${MODEL_LINK}" \
            "${infile}" \
            "${TASK_OUT}" \
            --repr_layers "${REPR_LAYERS}" \
            --include mean \
            --toks_per_batch 40000 2>&1

    # 2. VESM Logic
    elif [[ "$MODEL_NAME" == *"VESM"* ]]; then
        TASK_OUT="${EMB_DIR}/${suffix}"
        mkdir -p "${TASK_OUT}"
        
        python /home/laura.cano/projects/plms/esm/esm-scripts/extract_dist.py \
            "${MODEL_LINK}" \
            "${infile}" \
            "${TASK_OUT}" \
            --repr_layers "${REPR_LAYERS}" \
            --include mean \
            --toks_per_batch 16400 2>&1

    # 3. Legacy ESM (ESM-1b, ESM-2) Logic
    elif [[ "$MODEL_NAME" == *"esm"* || "$MODEL_NAME" == *"ESM"* ]]; then
        TASK_OUT="${EMB_DIR}/${suffix}"
        mkdir -p "${TASK_OUT}"
        
        python /home/laura.cano/projects/plms/esm/esm-scripts/extract.py \
            "${MODEL_LINK}" \
            "${infile}" \
            "${TASK_OUT}" \
            --repr_layers "${REPR_LAYERS}" \
            --include mean \
            --toks_per_batch 16400 2>&1

    # 4. ProtT5 Logic
    else
        # ProtT5 Logic: Single .npz file
        python /home/laura.cano/projects/plms/prostt5/scripts/prostt5_embedder_np.py \
            "${infile}" \
            --output "${EMB_DIR}/${PROTEIN_NAME}_${suffix}_${MODEL_NAME}.npz" \
            --model "${MODEL_LINK}" \
            --max-residues 32000 \
            --max-batch 36000 2>&1
    fi
done
### (3) Embeddings' Analysis
echo "Step 3: Analyzing Embedding Deviation..."
DEV_DIR="${PROT_OUT_DIR}/${MODEL_NAME}_deviation_analysis"
mkdir -p "${DEV_DIR}"

# Define input paths based on model type
if [[ "$MODEL_NAME" == *"esm"* || "$MODEL_NAME" == *"ESM"* ]]; then
    WT_IN="${EMB_DIR}/wt"
    SM_IN="${EMB_DIR}/sm"
    TR_IN="${EMB_DIR}/trunc"
else
    WT_IN="${EMB_DIR}/${PROTEIN_NAME}_wt_${MODEL_NAME}.npz"
    SM_IN="${EMB_DIR}/${PROTEIN_NAME}_sm_${MODEL_NAME}.npz"
    TR_IN="${EMB_DIR}/${PROTEIN_NAME}_trunc_${MODEL_NAME}.npz"
fi

python scripts/analyse_embedding_deviation.py \
    --wt_fasta "${FASTA}" \
    --wt_embedding "${WT_IN}" \
    --SM_embeddings "${SM_IN}" \
    --trunc_embeddings "${TR_IN}" \
    --output "${DEV_DIR}/${MODEL_NAME}" \
    --model "${MODEL_NAME}" \
    --sigma 1 2>&2

echo "Done. Results available in: ${PROT_OUT_DIR}"
