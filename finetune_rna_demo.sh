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
# fast_layernorm is used by default, no need to set explicitly. Set LAYERNORM_TYPE=torch to disable.
# export LAYERNORM_TYPE=fast_layernorm
# Kernel options:
# - triangle_attention: supports 'triattention', 'cuequivariance', 'deepspeed', 'torch'
# - triangle_multiplicative: supports 'cuequivariance', 'torch'

# Specify your data root directory by uncommenting the following line.
# export PROTENIX_ROOT_DIR="/modify/to/your/data_root_dir"
# wget -P $PROTENIX_ROOT_DIR/checkpoint/ https://protenix.tos-cn-beijing.volces.com/checkpoint/protenix_base_default_v1.0.0.pt
checkpoint_path="${PROTENIX_ROOT_DIR}/checkpoint/protenix_base_default_v1.0.0.pt"

# RNA finetuning modes:
# - full: conservative full-model finetuning
# - ss_head: prioritize the newly added RNA secondary-structure contact branch
RNA_FINETUNE_MODE="${RNA_FINETUNE_MODE:-full}"

common_args=(
  python3 ./runner/train.py
  --model_name protenix_base_default_v1.0.0
  --seed 42
  --base_dir ./output
  --dtype bf16
  --project protenix
  --use_wandb false
  --eval_interval 400
  --log_interval 50
  --checkpoint_interval 400
  --ema_decay 0.999
  --train_crop_size 384
  --diffusion_batch_size 48
  --max_steps 20000
  --warmup_steps 1000
  --lr_scheduler cosine_annealing
  --model.N_cycle 4
  --sample_diffusion.N_step 20
  --triangle_attention cuequivariance
  --triangle_multiplicative cuequivariance
  --load_checkpoint_path "${checkpoint_path}"
  --load_ema_checkpoint_path "${checkpoint_path}"
  --data.train_sets rna_monomer_train
  --data.test_sets rna_monomer_val
  --data.msa.enable_rna_msa true
  --data.template.enable_rna_template true
  --data.rna_monomer_train.constraint.enable true
  --data.rna_monomer_train.constraint.rna_ss.enable true
  --data.rna_monomer_train.constraint.contact.prob 0.0
  --model.constraint_embedder.contact_embedder.enable true
  --load_strict false
)

case "${RNA_FINETUNE_MODE}" in
  full)
    exec "${common_args[@]}" \
      --run_name rna_monomer_finetune_full \
      --lr 3e-4
    ;;
  ss_head)
    exec "${common_args[@]}" \
      --run_name rna_monomer_finetune_ss_head \
      --lr 1e-4 \
      --finetune.lr 1e-3 \
      --finetune.lr_scheduler cosine_annealing \
      --finetune.warmup_steps 1000 \
      --finetune.max_steps 20000 \
      --finetune_params_with_substring constraint_embedder.contact_z_embedder
    ;;
  *)
    echo "Unknown RNA_FINETUNE_MODE: ${RNA_FINETUNE_MODE}" >&2
    echo "Expected one of: full, ss_head" >&2
    exit 1
    ;;
esac
