# RNA Inference README

This document describes the current **RNA inference** workflow implemented in this
repository.

It covers:

- single-chain RNA inference
- optional RNA MSA
- optional RNA template hits (`.csv` / `.m8`) + template mmCIF
- optional predicted RNA secondary-structure constraints through `constraint.contact`

## Scope

This guide assumes:

- you are running inference with `protenix_base_default_v1.0.0` or
  `protenix_base_20250630_v1.0.0`
- your RNA input is represented as one or more `rnaSequence` entries in the input JSON
- optional RNA MSA, RNA template, and predicted RNA secondary structure are prepared
  ahead of time

Notes:

- RNA MSA can be generated automatically by `protenix prep` if `--use_rna_msa true`.
- RNA templates are currently **manual** at inference time. You must provide a
  `.csv` or `.m8` hit table through `rnaSequence.templatesPath`.
- Predicted RNA secondary structure is provided as same-chain `constraint.contact`
  entries in the input JSON.

## Recommended Model

For RNA inference, use:

- `protenix_base_default_v1.0.0`

Use `protenix_base_20250630_v1.0.0` if you want the more recent practical checkpoint.

If you want fair public benchmarking against models trained on the AlphaFold3-aligned
2021-09-30 cutoff, prefer `protenix_base_default_v1.0.0`.

## Folder Layout

Recommended layout:

```text
$PROTENIX_ROOT_DIR/
  checkpoint/
    protenix_base_default_v1.0.0.pt

  mmcif/
    4v4q.cif
    7lvk.cif
    9gmw.cif

  rna_inference/
    example_9gmw.json
    9gmw_2_all.a3m
    9gmw_2_templates.csv
```

Notes:

- `checkpoint/` is used automatically by `runner/inference.py` and `protenix pred`.
- `mmcif/` is needed if you use RNA templates and want local template structure files.
- The JSON file can use absolute paths, or paths relative to the current working
  directory.

## Input JSON

The top-level format is still a list:

```json
[
  {
    "name": "rna_9gmw",
    "modelSeeds": [101],
    "sequences": [
      {
        "rnaSequence": {
          "sequence": "GGCGCGUN",
          "count": 1,
          "unpairedMsaPath": "/abs/path/to/9gmw_2_all.a3m",
          "templatesPath": "/abs/path/to/9gmw_2_templates.csv",
          "templateQueryId": "9gmw_1_2_C"
        }
      }
    ],
    "constraint": {
      "contact": [
        {
          "entity1": 1,
          "copy1": 1,
          "position1": 5,
          "entity2": 1,
          "copy2": 1,
          "position2": 22,
          "min_distance": 0.0,
          "max_distance": 10.0
        }
      ]
    }
  }
]
```

Key RNA-specific fields:

- `rnaSequence.sequence`: RNA sequence
- `rnaSequence.unpairedMsaPath`: RNA MSA `.a3m`
- `rnaSequence.templatesPath`: RNA template hit table `.csv` or `.m8`
- `rnaSequence.templateQueryId`: optional query selector for multi-query template tables
- `constraint.contact`: optional predicted RNA secondary-structure contacts

For the full JSON schema, see [infer_json_format.md](./infer_json_format.md).

## RNA MSA

You have two options:

1. provide a precomputed RNA MSA file through `unpairedMsaPath`
2. let `protenix prep` generate RNA MSA automatically

Example auto-preprocessing:

```bash
protenix prep \
  --input /path/to/rna_input.json \
  --out_dir ./output \
  --nhmmer_binary_path /usr/bin/nhmmer \
  --hmmalign_binary_path /usr/bin/hmmalign \
  --hmmbuild_rna_binary_path /usr/bin/hmmbuild
```

This updates the JSON with `unpairedMsaPath`.

## RNA Template

RNA templates are inference-only and manual.

Expected input:

- `rnaSequence.templatesPath = /abs/path/to/result.csv`
  or
- `rnaSequence.templatesPath = /abs/path/to/result.m8`

The table format is:

```text
query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
```

The template structure itself is loaded from local or remote mmCIF using the `target`
PDB ID and chain information.

Notes:

- if your template table has exactly one query, `templateQueryId` can be omitted
- if your template table contains multiple queries, set `templateQueryId`
- for fair evaluation, remove self-hits and same-PDB hits before inference if needed

## Predicted RNA Secondary Structure

Predicted RNA secondary structure should be converted into `constraint.contact`.

Use:

```bash
python scripts/rna_ss_to_constraint_json.py \
  --input /path/to/ss.dbn \
  --wrap constraint \
  --output /path/to/rna_constraint.json
```

Then copy the resulting `constraint.contact` block into your inference JSON.

Important:

- same-chain RNA contacts are supported in inference
- same-chain protein contacts are still not the intended path here
- for RNA secondary-structure constraints, you should run `runner/inference.py` with:
  - `--load_strict false`
  - `--model.constraint_embedder.contact_embedder.enable true`

## RNA Inference Script

Use [inference_rna_demo.sh](../inference_rna_demo.sh).

It supports two modes:

- `RNA_INFER_MODE=basic`
  - uses `protenix pred`
  - suitable for RNA sequence only, RNA MSA, or RNA template inference
- `RNA_INFER_MODE=ss`
  - uses `runner/inference.py`
  - enables the contact constraint embedder for predicted RNA secondary structure

Examples:

```bash
export PROTENIX_ROOT_DIR=/path/to/your/data_root
INPUT_JSON=examples/examples_with_rna_msa/example_9gmw_2.json \
USE_RNA_MSA=true \
USE_TEMPLATE=false \
bash inference_rna_demo.sh
```

```bash
export PROTENIX_ROOT_DIR=/path/to/your/data_root
INPUT_JSON=/path/to/rna_with_template.json \
USE_RNA_MSA=true \
USE_TEMPLATE=true \
bash inference_rna_demo.sh
```

```bash
export PROTENIX_ROOT_DIR=/path/to/your/data_root
RNA_INFER_MODE=ss \
INPUT_JSON=/path/to/rna_with_template_and_ss.json \
USE_RNA_MSA=true \
USE_TEMPLATE=true \
bash inference_rna_demo.sh
```

## Direct Commands

### RNA MSA or RNA Template Only

```bash
protenix pred \
  -i /path/to/rna_input.json \
  -o ./output \
  -n protenix_base_default_v1.0.0 \
  --use_default_params true \
  --use_rna_msa true \
  --use_template true
```

### RNA Secondary Structure Constraint

```bash
python runner/inference.py \
  --model_name protenix_base_default_v1.0.0 \
  --input_json_path /path/to/rna_with_ss.json \
  --dump_dir ./output \
  --seeds 101 \
  --model.N_cycle 10 \
  --sample_diffusion.N_sample 5 \
  --sample_diffusion.N_step 200 \
  --triangle_attention cuequivariance \
  --triangle_multiplicative cuequivariance \
  --use_rna_msa true \
  --use_template true \
  --load_strict false \
  --model.constraint_embedder.contact_embedder.enable true
```

## Practical Notes

- `protenix pred` is the simplest route for RNA sequence, RNA MSA, and RNA template
  inference.
- If you want predicted RNA secondary-structure constraints, use
  `runner/inference.py` or `inference_rna_demo.sh` with `RNA_INFER_MODE=ss`.
- RNA template inference with local mmCIF is easiest if
  `data.template.prot_template_mmcif_dir` / `data.template.rna_template_mmcif_dir`
  both point to `$PROTENIX_ROOT_DIR/mmcif`.
- For RNA template realignment, `kalign` is required unless query and template
  sequences are identical.

## Related Files

- [infer_json_format.md](./infer_json_format.md)
- [training_inference_instructions.md](./training_inference_instructions.md)
- [rna_training_readme.md](./rna_training_readme.md)
- [inference_demo.sh](../inference_demo.sh)
- [inference_rna_demo.sh](../inference_rna_demo.sh)
- [rna_ss_to_constraint_json.py](../scripts/rna_ss_to_constraint_json.py)
