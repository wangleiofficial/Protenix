# RNA Training README

This document describes the data layout and configuration required to train **monomer RNA structure** models in the current Protenix codebase.

It is written for the workflow that is already implemented in this repository:

- single-chain RNA `Bioassembly Dict`
- optional RNA MSA
- optional RNA template hits (`.csv` / `.m8`) + template mmCIF
- optional predicted RNA secondary structure converted to `constraint.contact`

## Scope

This guide assumes:

- you are training **single-chain RNA** structures
- your training samples are stored as Protenix `Bioassembly Dict` `.pkl.gz` files
- you use `prepare_training_data.py --rna-mode rna_monomer` or an equivalent pipeline
  to build **single-chain RNA samples**

Important notes:

- If you use `sequence_uid`-based mapping, the current convention is:
  - `sequence_uid = <pdb_id>_<asym_id_int>`
- If you split chains out of complexes, keep the intended `asym_id_int` stable when
  building mapping files. The current `rna_monomer` mode already does this by naming
  each saved sample as `<pdb_id>_<asym_id_int>.pkl.gz`.
- External predicted RNA secondary-structure constraints are currently interpreted in the **token order of the training sample**. In practice, this is easiest when a training sample contains the **full RNA chain** and fits inside the crop.

## Required Data

The minimum required data for RNA training is:

1. `Bioassembly Dict` directory
2. `indices` CSV

Optional but supported:

1. RNA MSA
2. RNA templates
3. predicted RNA secondary-structure constraints

## Optional: Build a Raw RNA PDB Corpus

If you do **not** already have a local RNA structure collection, you can first build one
from the PDB directly.

The helper script below will:

1. search RCSB for entries with RNA polymers
2. download the corresponding `mmCIF` files into one directory

Example:

```bash
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
export PROTENIX_ROOT_DIR=/path/to/your/data_root

python scripts/database/download_rna_pdb_dataset.py \
  --mmcif-dir $PROTENIX_ROOT_DIR/mmcif \
  --entry-id-list-out $PROTENIX_ROOT_DIR/indices/rna_pdb_ids.txt
```

Then preprocess the downloaded `mmCIF` files with an explicit RNA mode:

```bash
python scripts/prepare_training_data.py \
  -i $PROTENIX_ROOT_DIR/mmcif \
  -o $PROTENIX_ROOT_DIR/indices/rna_monomer_train.csv \
  -b $PROTENIX_ROOT_DIR/rna_monomer_bioassembly \
  -c /path/to/clusters-by-entity-40.txt \
  -n 8 \
  --rna-mode rna_monomer
```

Notes:

- By default the search uses `rcsb_entry_info.polymer_entity_count_RNA > 0`
  with `resolution <= 4.5`, which is closer to the current Protenix training filters.
- This gives you a **raw RNA-related corpus**.
- RNA-specific preprocessing is now **separated** from download. Use
  `prepare_training_data.py --rna-mode ...` to choose the subset you want.
- Supported `--rna-mode` values are:
  - `rna_monomer`: extract every RNA chain and save it as its own single-chain
    sample, including RNA chains originating from RNA complexes, RNA-protein
    complexes, or DNA-RNA complexes
  - `rna_only`: RNA-only assemblies, including RNA monomers and RNA-RNA complexes
  - `rna_all`: RNA monomers, RNA-RNA complexes, and RNA-protein complexes
- DNA and DNA-RNA hybrid chains are **not** emitted as `rna_monomer` samples.
- `rna_only` and `rna_all` still exclude assemblies containing DNA / DNA-RNA hybrid
  polymer chains.
- Useful helper:
  [download_rna_pdb_dataset.py](../scripts/database/download_rna_pdb_dataset.py)
  and [prepare_training_data.py](../scripts/prepare_training_data.py)

## Recommended Folder Layout

```text
$PROTENIX_ROOT_DIR/
  mmcif/
    1abc.cif
    2xyz.cif
    4v4q.cif
    9gmw.cif

  rna_monomer_bioassembly/
    1abc_2.pkl.gz
    1abc_5.pkl.gz
    2xyz_0.pkl.gz

  indices/
    rna_monomer_train.csv
    rna_monomer_val.csv

  rna_msa/
    mapping_sequence.json
    msas/
      rna_000001/
        rna_000001_all.a3m
      rna_000002/
        rna_000002_all.a3m

  rna_template_hits/
    mapping_sequence_uid.json
    1abc_2.csv
    2xyz_0.m8

  rna_ss_constraints/
    mapping_sequence_uid.json
    1abc_2.json
    2xyz_0.json
```

Notes:

- One shared `mmcif/` directory is usually enough. It can serve both:
  - `data.<dataset>.base_info.mmcif_dir`
  - `data.template.rna_template_mmcif_dir`
- The only requirement is that this directory must contain every `<pdb_id>.cif` needed by your RNA template hit tables, not just the training samples themselves.
- If you prefer, you can still keep a separate template-specific mmCIF directory, but it is optional.
- `rna_template_hits/` and `rna_ss_constraints/` are both keyed by `sequence_uid`.
- `rna_msa/` is keyed by **sequence**, not `sequence_uid`.

## 1. Bioassembly Dict

The `.pkl.gz` schema does **not** need to be changed for RNA template or RNA secondary-structure support. The format is the same as described in [prepare_training_data.md](./prepare_training_data.md).

For monomer RNA training, the practical recommendation is:

- one sample contains one RNA chain
- `atom_array` only contains that RNA chain
- `token_array` is regenerated from that RNA chain
- `entity_poly_type` contains only the RNA entity
- `num_prot_chains = 0`
- `num_assembly_polymer_chains = 1`

When using `prepare_training_data.py --rna-mode rna_monomer`, the saved sample name is
`<pdb_id>_<asym_id_int>.pkl.gz`, and the stored `atom_array` keeps that original
`asym_id_int`. This is the recommended basis for `sequence_uid`-keyed RNA template and
predicted secondary-structure mappings.

## 2. Indices CSV

The indices CSV is still required by the dataset loader. The schema follows [prepare_training_data.md](./prepare_training_data.md).

For monomer RNA chain rows, the important fields are:

- `type = chain`
- `pdb_id`
- `bioassembly_name`
- `entity_1_id`
- `chain_1_id`
- `mol_1_type = nuc`
- `cluster_1_id`

For `type = chain` rows, the chain-2-related columns can remain empty.

When `rna_monomer` is produced automatically from mixed assemblies:

- `pdb_id` remains the original PDB code
- `bioassembly_name = <pdb_id>_<asym_id_int>` is used to locate the saved `.pkl.gz`
- `sequence_uid` for templates and RNA secondary structure still follows
  `<pdb_id>_<asym_id_int>`

## 2.1 Time-Based Test Split

For monomer RNA evaluation, a common setup is:

- use structures released **before** a cutoff date for training
- use structures released **on or after** that cutoff date for testing
- remove test samples whose canonical RNA sequence is **100% identical** to any
  training sample

The helper script below does exactly this for RNA chain rows:

```bash
python scripts/build_rna_time_split.py \
  --indices $PROTENIX_ROOT_DIR/indices/rna_monomer_all.csv \
  --bioassembly-dir $PROTENIX_ROOT_DIR/rna_monomer_bioassembly \
  --test-release-date-from 2021-09-30 \
  --train-out $PROTENIX_ROOT_DIR/indices/rna_monomer_train.csv \
  --test-out $PROTENIX_ROOT_DIR/indices/rna_monomer_val.csv \
  --dropped-test-out $PROTENIX_ROOT_DIR/indices/rna_monomer_test_dropped_exact_train_match.csv \
  --test-pdb-list-out $PROTENIX_ROOT_DIR/indices/rna_monomer_test_pdb_id.txt
```

Notes:

- The script only considers `type = chain` rows with `sub_mol_1_type` in
  `{rna, modified_rna}`.
- Exact-match filtering is done on the canonical single-chain sequence loaded from the
  saved `Bioassembly Dict` files.
- The output test set is therefore a **time-split RNA monomer set with 100% sequence
  matches to train removed**.

## 3. RNA MSA Format

RNA MSA support in training is sequence-based.

Current expectations:

- mapping key: **RNA sequence**
- mapping value: a list containing one directory ID
- final file path:
  - `<rna_msadir_raw_paths>/<eid>/<eid>_all.a3m`

Example:

```json
{
  "ACCAGGAUGGCCGAGUGGUUAAGGCGUUGGACUUAAGAUCCAAUGGACAUAUGUCCGCGUGGGUUCGAACCCCACUCCUGGUACCA": ["rna_000001"]
}
```

Corresponding file:

```text
$PROTENIX_ROOT_DIR/rna_msa/msas/rna_000001/rna_000001_all.a3m
```

Notes:

- RNA MSA currently uses `sequence`, not `sequence_uid`
- format must be `.a3m`

## 4. RNA Template Format

RNA template support in training is `sequence_uid`-based.

Current expectations:

- mapping key: `sequence_uid`
- mapping value:
  - either a direct `.csv` / `.m8` file
  - or a directory containing `result.csv` or `result.m8`
- template structure source:
  - local mmCIF file named `<pdb_id>.cif`

Recommended mapping example:

```json
{
  "1abc_2": "/abs/path/to/rna_template_hits/1abc_2.csv",
  "2xyz_0": "/abs/path/to/rna_template_hits/2xyz_0.m8"
}
```

Supported hit table columns:

```text
query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
```

Recommended target naming:

```text
<pdb_id>_<assembly_id>_<entity_id>_<label_asym_id>
```

Example:

```text
9gmw_1_2_C
```

Useful helper:

- [split_rna_template_hits.py](../scripts/split_rna_template_hits.py)

Example:

```bash
python scripts/split_rna_template_hits.py \
  --input /path/to/result.csv \
  --output-dir $PROTENIX_ROOT_DIR/rna_template_hits \
  --mapping-out $PROTENIX_ROOT_DIR/rna_template_hits/mapping_sequence_uid.json \
  --query-map-json /path/to/query_to_sequence_uid.json \
  --drop-self-hit \
  --drop-same-pdb
```

The `--drop-self-hit` and `--drop-same-pdb` options are recommended for training to avoid trivial leakage from the query structure into the template branch.

## 5. Predicted RNA Secondary Structure Format

Predicted RNA secondary structure is supported as external **same-chain RNA contact constraints**.

Current expectations:

- mapping key: `sequence_uid`
- mapping value:
  - direct JSON file
  - or directory containing `constraint.json` or `contact.json`

Supported JSON payloads:

```json
{
  "contact": [
    {
      "entity1": 1,
      "copy1": 1,
      "position1": 5,
      "entity2": 1,
      "copy2": 1,
      "position2": 22,
      "max_distance": 10.0,
      "min_distance": 0.0
    }
  ]
}
```

or:

```json
{
  "constraint": {
    "contact": [
      {
        "entity1": 1,
        "copy1": 1,
        "position1": 5,
        "entity2": 1,
        "copy2": 1,
        "position2": 22,
        "max_distance": 10.0,
        "min_distance": 0.0
      }
    ]
  }
}
```

Useful helper:

- [rna_ss_to_constraint_json.py](../scripts/rna_ss_to_constraint_json.py)

Examples:

```bash
python scripts/rna_ss_to_constraint_json.py \
  --input /path/to/ss.dbn \
  --wrap constraint \
  --output $PROTENIX_ROOT_DIR/rna_ss_constraints/1abc_2.json
```

```bash
python scripts/rna_ss_to_constraint_json.py \
  --input /path/to/pairs.csv \
  --input-format pairs \
  --pair-index-base 0 \
  --wrap constraint \
  --output $PROTENIX_ROOT_DIR/rna_ss_constraints/1abc_2.json
```

Practical recommendations:

- use **predicted** secondary structure for training if you plan to use predicted secondary structure at inference
- prefer high-confidence pairs only
- start with token-level contact and keep:
  - `max_distance = 8` to `10`
  - `min_distance = 0`

## 6. sequence_uid Convention

The current RNA template and external RNA secondary-structure mappings both support:

- `sequence_uid = <pdb_id>_<asym_id_int>`

Example:

```text
1abc_2
```

If you enable `reassign_continuous_chain_ids`, rebuild your `sequence_uid` mappings after that choice is finalized. Otherwise the lookup key may not match the final `asym_id_int` used during training.

## 7. Config Example

The repository now ships with built-in dataset entries:

- `rna_monomer_train`
- `rna_monomer_val`

They are defined in [configs_data.py](../configs/configs_data.py) and expect the following default layout under `PROTENIX_ROOT_DIR`:

- `mmcif/`
- `rna_monomer_bioassembly/`
- `indices/rna_monomer_train.csv`
- `indices/rna_monomer_val.csv`

If your files follow that layout, you can use these dataset names directly without editing `configs/configs_data.py`.

For RNA feature loading, the built-in config now also auto-points to the standard paths:

- RNA MSA mapping: `rna_msa/mapping_sequence.json`
  fallback: `rna_msa/rna_sequence_to_pdb_chains.json`
- RNA MSA directory: `rna_msa/msas/`
- RNA template mapping: `rna_template_hits/mapping_sequence_uid.json`
- RNA template root: `rna_template_hits/`
- predicted RNA SS mapping: `rna_ss_constraints/mapping_sequence_uid.json`
- predicted RNA SS root: `rna_ss_constraints/`

So in the standard layout, you usually only need to:

1. set `PROTENIX_ROOT_DIR`
2. place files in the default directories
3. enable the corresponding feature flags

### Dataset entry

```python
data_configs["rna_monomer_train"] = {
    "base_info": {
        "mmcif_dir": os.path.join(PROTENIX_ROOT_DIR, "mmcif"),
        "bioassembly_dict_dir": os.path.join(PROTENIX_ROOT_DIR, "rna_monomer_bioassembly"),
        "indices_fpath": os.path.join(PROTENIX_ROOT_DIR, "indices/rna_monomer_train.csv"),
        "pdb_list": "",
        "random_sample_if_failed": True,
        "max_n_token": -1,
        "use_reference_chains_only": False,
        "exclusion": {},
    },
    **deepcopy(default_weighted_pdb_configs),
}

data_configs["rna_monomer_val"] = {
    "base_info": {
        "mmcif_dir": os.path.join(PROTENIX_ROOT_DIR, "mmcif"),
        "bioassembly_dict_dir": os.path.join(PROTENIX_ROOT_DIR, "rna_monomer_bioassembly"),
        "indices_fpath": os.path.join(PROTENIX_ROOT_DIR, "indices/rna_monomer_val.csv"),
        "pdb_list": "",
        "max_n_token": GlobalConfigValue("test_max_n_token"),
        "sort_by_n_token": False,
        "group_by_pdb_id": True,
        "find_eval_chain_interface": False,
    },
    **deepcopy(default_test_configs),
}
```

### RNA MSA

```python
data_configs["msa"]["enable_rna_msa"] = True
data_configs["msa"]["rna_seq_or_filename_to_msadir_jsons"] = [
    os.path.join(PROTENIX_ROOT_DIR, "rna_msa/mapping_sequence.json")
]
data_configs["msa"]["rna_msadir_raw_paths"] = [
    os.path.join(PROTENIX_ROOT_DIR, "rna_msa/msas")
]
data_configs["msa"]["rna_indexing_methods"] = ["sequence"]
```

### RNA template

```python
data_configs["template"]["enable_rna_template"] = True
data_configs["template"]["rna_template_mmcif_dir"] = os.path.join(
    PROTENIX_ROOT_DIR, "mmcif"
)
data_configs["template"]["rna_template_raw_paths"] = [
    os.path.join(PROTENIX_ROOT_DIR, "rna_template_hits")
]
data_configs["template"]["rna_seq_or_filename_to_templatedir_jsons"] = [
    os.path.join(PROTENIX_ROOT_DIR, "rna_template_hits/mapping_sequence_uid.json")
]
data_configs["template"]["rna_indexing_methods"] = ["sequence_uid"]
```

### Predicted RNA secondary structure

```python
data_configs["rna_monomer_train"]["constraint"]["enable"] = True
data_configs["rna_monomer_train"]["constraint"]["contact"]["prob"] = 0.0
data_configs["rna_monomer_train"]["constraint"]["rna_ss"]["enable"] = True
data_configs["rna_monomer_train"]["constraint"]["rna_ss"]["raw_paths"] = [
    os.path.join(PROTENIX_ROOT_DIR, "rna_ss_constraints")
]
data_configs["rna_monomer_train"]["constraint"]["rna_ss"][
    "seq_or_filename_to_ss_jsons"
] = [
    os.path.join(PROTENIX_ROOT_DIR, "rna_ss_constraints/mapping_sequence_uid.json")
]
data_configs["rna_monomer_train"]["constraint"]["rna_ss"]["indexing_methods"] = [
    "sequence_uid"
]
```

### Model-side contact embedder

Enable the contact constraint branch in the model config:

```python
configs.model.constraint_embedder.contact_embedder.enable = True
```

## 8. Fine-tuning Command

In practice, the recommended RNA training workflow is to **fine-tune from a released
Protenix checkpoint**, not to start from random initialization.

If you use:

- RNA MSA only
- RNA template only

then the released `protenix_base_default_v1.0.0` checkpoint already matches the base
architecture.

If you additionally enable predicted RNA secondary structure through
`model.constraint_embedder.contact_embedder.enable = True`, you are adding a branch
that is **disabled in the released base checkpoint config**. In that case, load with
`load_strict = false` so the new constraint-embedder parameters can be initialized
instead of causing a strict load failure.

Also note:

- if EMA is enabled (`ema_decay > 0`), set both `load_checkpoint_path` and
  `load_ema_checkpoint_path`
- using only `load_checkpoint_path` is not enough for a clean fine-tuning start when EMA
  is active
- the generic `finetune_demo.sh` values (`100000` steps, `2000` warmup, `0.001` lr)
  are too aggressive for most RNA fine-tuning runs

Recommended starting points:

- full-model RNA fine-tuning:
  - `lr = 3e-4`
  - `warmup_steps = 1000`
  - `max_steps = 20000`
  - `lr_scheduler = cosine_annealing`
- RNA secondary-structure + template prioritized fine-tuning:
  - backbone `lr = 1e-4`
  - `constraint_embedder.contact_z_embedder` `finetune.lr = 1e-3`
  - `template_embedder` `finetune.lr = 1e-3`
  - `warmup_steps = 1000`
  - `max_steps = 20000`
  - `lr_scheduler = cosine_annealing`

If you want a conservative full-model RNA fine-tuning command, use:

```bash
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
export PROTENIX_ROOT_DIR=/path/to/your/data_root
checkpoint_path="${PROTENIX_ROOT_DIR}/checkpoint/protenix_base_default_v1.0.0.pt"

python3 runner/train.py \
  --run_name rna_monomer_train \
  --model_name protenix_base_default_v1.0.0 \
  --base_dir ./output \
  --dtype bf16 \
  --use_wandb false \
  --train_crop_size 384 \
  --diffusion_batch_size 48 \
  --max_steps 20000 \
  --warmup_steps 1000 \
  --lr_scheduler cosine_annealing \
  --lr 3e-4 \
  --model.N_cycle 4 \
  --sample_diffusion.N_step 20 \
  --triangle_attention cuequivariance \
  --triangle_multiplicative cuequivariance \
  --load_checkpoint_path ${checkpoint_path} \
  --load_ema_checkpoint_path ${checkpoint_path} \
  --load_strict false \
  --data.train_sets rna_monomer_train \
  --data.test_sets rna_monomer_val \
  --data.msa.enable_rna_msa true \
  --data.template.enable_rna_template true \
  --data.rna_monomer_train.constraint.enable true \
  --data.rna_monomer_train.constraint.rna_ss.enable true \
  --data.rna_monomer_train.constraint.contact.prob 0.0 \
  --model.constraint_embedder.contact_embedder.enable true
```

If you want to prioritize the RNA secondary-structure contact branch and RNA template
adaptation while still updating the backbone more conservatively, use:

```bash
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
export PROTENIX_ROOT_DIR=/path/to/your/data_root
checkpoint_path="${PROTENIX_ROOT_DIR}/checkpoint/protenix_base_default_v1.0.0.pt"

python3 runner/train.py \
  --run_name rna_monomer_train_ss_template_head \
  --model_name protenix_base_default_v1.0.0 \
  --base_dir ./output \
  --dtype bf16 \
  --use_wandb false \
  --train_crop_size 384 \
  --diffusion_batch_size 48 \
  --max_steps 20000 \
  --warmup_steps 1000 \
  --lr_scheduler cosine_annealing \
  --lr 1e-4 \
  --finetune.lr 1e-3 \
  --finetune.lr_scheduler cosine_annealing \
  --finetune.warmup_steps 1000 \
  --finetune.max_steps 20000 \
  --finetune_params_with_substring constraint_embedder.contact_z_embedder,template_embedder \
  --model.N_cycle 4 \
  --sample_diffusion.N_step 20 \
  --triangle_attention cuequivariance \
  --triangle_multiplicative cuequivariance \
  --load_checkpoint_path ${checkpoint_path} \
  --load_ema_checkpoint_path ${checkpoint_path} \
  --load_strict false \
  --data.train_sets rna_monomer_train \
  --data.test_sets rna_monomer_val \
  --data.msa.enable_rna_msa true \
  --data.template.enable_rna_template true \
  --data.rna_monomer_train.constraint.enable true \
  --data.rna_monomer_train.constraint.rna_ss.enable true \
  --data.rna_monomer_train.constraint.contact.prob 0.0 \
  --model.constraint_embedder.contact_embedder.enable true
```

If you are **not** using predicted RNA secondary structure and therefore keep
`model.constraint_embedder.contact_embedder.enable = false`, you can usually keep
`load_strict = true`.

If you prefer a runnable demo script, use [finetune_rna_demo.sh](../finetune_rna_demo.sh)
for RNA fine-tuning and [finetune_demo.sh](../finetune_demo.sh) for the general
non-RNA example. The RNA demo supports:

- `RNA_FINETUNE_MODE=full`
- `RNA_FINETUNE_MODE=ss_template_head`
- `RNA_FINETUNE_MODE=ss_head` (backward-compatible alias)

Practical note:

- `contact_z_embedder` is a newly enabled branch, so giving it a higher learning rate is
  important when you use predicted RNA secondary structure
- `template_embedder` already exists in the released v1.0.0 checkpoint, but if your goal
  is to improve **RNA template utilization**, it is usually better to include it in
  `finetune_params_with_substring` as well, rather than leaving it only on the smaller
  backbone learning rate

## 9. Checklist

Before starting training, verify the following:

- every RNA training sample is a single-chain `Bioassembly Dict`
- your `indices` CSV points to the correct chain rows
- RNA MSA mapping is keyed by **sequence**
- RNA template mapping is keyed by **sequence_uid**
- predicted RNA secondary-structure mapping is keyed by **sequence_uid**
- template self-hits were removed for training
- `constraint.contact.prob = 0.0` if you only want external predicted RNA SS, not random GT contact
- `model.constraint_embedder.contact_embedder.enable = True`
- your monomer RNA chains fit inside the training crop, or you otherwise control cropping so token positions remain aligned with the external RNA SS contacts

## 10. Related Files

- [prepare_training_data.md](./prepare_training_data.md)
- [training_inference_instructions.md](./training_inference_instructions.md)
- [infer_json_format.md](./infer_json_format.md)
- [download_rna_pdb_dataset.py](../scripts/database/download_rna_pdb_dataset.py)
- [build_rna_time_split.py](../scripts/build_rna_time_split.py)
- [split_rna_template_hits.py](../scripts/split_rna_template_hits.py)
- [rna_ss_to_constraint_json.py](../scripts/rna_ss_to_constraint_json.py)
