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
- if your RNA chains were extracted from complexes, you already split them into **single-chain samples**

Important notes:

- If you use `sequence_uid`-based mapping, the current convention is:
  - `sequence_uid = <pdb_id>_<asym_id_int>`
- If you split chains out of complexes, keep the intended `asym_id_int` stable when building mapping files.
- External predicted RNA secondary-structure constraints are currently interpreted in the **token order of the training sample**. In practice, this is easiest when a training sample contains the **full RNA chain** and fits inside the crop.

## Required Data

The minimum required data for RNA training is:

1. `Bioassembly Dict` directory
2. `indices` CSV

Optional but supported:

1. RNA MSA
2. RNA templates
3. predicted RNA secondary-structure constraints

## Recommended Folder Layout

```text
$PROTENIX_ROOT_DIR/
  mmcif/
    1abc.cif
    2xyz.cif

  rna_monomer_bioassembly/
    1abc.pkl.gz
    2xyz.pkl.gz

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

  rna_template_mmcif/
    1abc.cif
    4v4q.cif
    9gmw.cif

  rna_ss_constraints/
    mapping_sequence_uid.json
    1abc_2.json
    2xyz_0.json
```

Notes:

- `rna_template_mmcif/` can be the same directory as your main mmCIF library if it already contains all template structures.
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

If the chain came from a complex, keep the intended `asym_id_int` consistent with the `sequence_uid` mapping you plan to use for templates and predicted secondary structure.

## 2. Indices CSV

The indices CSV is still required by the dataset loader. The schema follows [prepare_training_data.md](./prepare_training_data.md).

For monomer RNA chain rows, the important fields are:

- `type = chain`
- `pdb_id`
- `entity_1_id`
- `chain_1_id`
- `mol_1_type = nuc`
- `cluster_1_id`

For `type = chain` rows, the chain-2-related columns can remain empty.

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

Below is a minimal example of how to wire a custom monomer RNA dataset into `configs/configs_data.py`.

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
    PROTENIX_ROOT_DIR, "rna_template_mmcif"
)
data_configs["template"]["rna_template_raw_paths"] = [""]
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
data_configs["rna_monomer_train"]["constraint"]["rna_ss"]["raw_paths"] = [""]
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

## 8. Training Command

After editing your dataset/config entries, a typical launch command looks like:

```bash
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
export PROTENIX_ROOT_DIR=/path/to/your/data_root

python3 runner/train.py \
  --run_name rna_monomer_train \
  --model_name protenix_base_default_v1.0.0 \
  --base_dir ./output \
  --dtype bf16 \
  --use_wandb false \
  --train_crop_size 384 \
  --diffusion_batch_size 48 \
  --max_steps 100000 \
  --warmup_steps 2000 \
  --lr 0.001 \
  --model.N_cycle 4 \
  --sample_diffusion.N_step 20 \
  --triangle_attention cuequivariance \
  --triangle_multiplicative cuequivariance \
  --data.train_sets rna_monomer_train \
  --data.test_sets rna_monomer_val \
  --data.msa.enable_rna_msa true \
  --data.template.enable_rna_template true \
  --data.rna_monomer_train.constraint.enable true \
  --data.rna_monomer_train.constraint.rna_ss.enable true \
  --data.rna_monomer_train.constraint.contact.prob 0.0 \
  --model.constraint_embedder.contact_embedder.enable true
```

If you prefer the demo script, you can also start from [train_demo.sh](../train_demo.sh) and replace the dataset/config arguments there.

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
- [split_rna_template_hits.py](../scripts/split_rna_template_hits.py)
- [rna_ss_to_constraint_json.py](../scripts/rna_ss_to_constraint_json.py)
