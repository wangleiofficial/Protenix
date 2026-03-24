# Copyright 2024 ByteDance and/or its affiliates.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# RNA inference modes:
# - basic: RNA sequence only, with optional RNA MSA/template toggles
# - ss: same as basic, but also enables the contact constraint branch for
#   same-chain RNA secondary-structure constraints already embedded in the JSON
RNA_INFER_MODE="${RNA_INFER_MODE:-basic}"

# Required environment:
# export PROTENIX_ROOT_DIR="/modify/to/your/data_root_dir"

MODEL_NAME="${MODEL_NAME:-protenix_base_default_v1.0.0}"
INPUT_JSON="${INPUT_JSON:-examples/examples_with_rna_msa/example_9gmw_2.json}"
OUT_DIR="${OUT_DIR:-./test_outputs/rna_inference}"
SEEDS="${SEEDS:-101}"
N_CYCLE="${N_CYCLE:-10}"
N_STEP="${N_STEP:-200}"
N_SAMPLE="${N_SAMPLE:-5}"
DTYPE="${DTYPE:-bf16}"
USE_TEMPLATE="${USE_TEMPLATE:-false}"
USE_RNA_MSA="${USE_RNA_MSA:-true}"
TRIANGLE_ATTENTION="${TRIANGLE_ATTENTION:-cuequivariance}"
TRIANGLE_MULTIPLICATIVE="${TRIANGLE_MULTIPLICATIVE:-cuequivariance}"

case "${RNA_INFER_MODE}" in
  basic)
    exec protenix pred \
      -i "${INPUT_JSON}" \
      -o "${OUT_DIR}" \
      -s "${SEEDS}" \
      -c "${N_CYCLE}" \
      -p "${N_STEP}" \
      -e "${N_SAMPLE}" \
      -d "${DTYPE}" \
      -n "${MODEL_NAME}" \
      --use_default_params false \
      --use_template "${USE_TEMPLATE}" \
      --use_rna_msa "${USE_RNA_MSA}" \
      --triatt_kernel "${TRIANGLE_ATTENTION}" \
      --trimul_kernel "${TRIANGLE_MULTIPLICATIVE}"
    ;;
  ss)
    exec python3 runner/inference.py \
      --model_name "${MODEL_NAME}" \
      --seeds "${SEEDS}" \
      --dump_dir "${OUT_DIR}" \
      --input_json_path "${INPUT_JSON}" \
      --dtype "${DTYPE}" \
      --model.N_cycle "${N_CYCLE}" \
      --sample_diffusion.N_sample "${N_SAMPLE}" \
      --sample_diffusion.N_step "${N_STEP}" \
      --triangle_attention "${TRIANGLE_ATTENTION}" \
      --triangle_multiplicative "${TRIANGLE_MULTIPLICATIVE}" \
      --use_template "${USE_TEMPLATE}" \
      --use_rna_msa "${USE_RNA_MSA}" \
      --load_strict false \
      --model.constraint_embedder.contact_embedder.enable true
    ;;
  *)
    echo "Unknown RNA_INFER_MODE: ${RNA_INFER_MODE}" >&2
    echo "Expected one of: basic, ss" >&2
    exit 1
    ;;
esac
