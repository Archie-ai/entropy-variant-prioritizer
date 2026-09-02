# Entropy-Based Cross-Species Variant Prioritization

## Project continuity record

This file is the version-controlled source of truth for the research project. Update it after each meaningful analysis milestone so the project does not depend on a chat transcript.

Last updated: 2026-09-01  
Repository: `https://github.com/Archie-ai/entropy-variant-prioritizer.git`  
Current branch: `main`  
Latest verified milestone commit: `fb96cd6` (`Add strict ClinVar benchmark class counts`)

## 1. Project identity

Working title: **Entropy-Based Cross-Species Variant Prioritization Integrating Functional Prediction and RNA-Seq Validation**

Project type: Individual, disease-agnostic research publication and open-source bioinformatics tool. This is separate from the group capstone project, *Microbiome Abundance Consensus Machine*.

Primary goal: Determine whether Shannon entropy calculated from cross-species protein multiple-sequence alignments contributes useful evolutionary-constraint information beyond established variant-effect predictors, and whether integrating entropy, predictor scores, and RNA-seq evidence improves pathogenic-versus-benign variant prioritization.

## 2. Central hypothesis

Cross-species protein-MSA-derived Shannon entropy captures evolutionary constraint that is not fully represented by conventional functional predictors. A combined model using entropy, established predictor scores, and appropriately designed RNA-seq evidence will prioritize clinically classified missense variants more effectively than conventional predictors alone across diverse disease contexts.

Low entropy represents a highly conserved alignment column. High entropy represents a more variable column. Entropy is evidence of evolutionary constraint; it is not itself a clinical pathogenicity label.

## 3. Planned model comparison

- **M1 - Predictors:** established functional-prediction scores without the new entropy feature.
- **M2 - Entropy:** cross-species Shannon entropy and alignment-quality features.
- **M3 - Predictors + entropy:** tests whether entropy adds information beyond conventional predictors.
- **M4 - Predictors + entropy + RNA evidence:** evaluates the additional contribution of transcriptomic evidence.

Planned primary evaluation measures:

- ROC-AUC
- PR-AUC
- Matthews correlation coefficient (MCC)
- Calibration
- Sensitivity and specificity at prespecified or transparently selected thresholds

Evaluation must use gene-held-out splitting so variants from the same gene do not appear in both training and test sets. This reduces gene-level information leakage and better tests generalization to unseen genes.

## 4. Benchmark 1A definition

Benchmark 1A is the strict primary human clinical-label benchmark.

- Reference assembly: GRCh38
- Species: Human
- Variant context: Germline ClinVar records
- Initial molecular class: Missense single-nucleotide variants
- Positive label (`1`): Pathogenic or Likely pathogenic
- Negative label (`0`): Benign or Likely benign
- Excluded labels: uncertain significance, conflicting classifications, association, risk factor, drug response, and other non-binary interpretations
- Strict evidence target: ClinVar 2- to 4-star review categories
- ClinVar source file: `data/raw/clinvar_20260822.vcf.gz`

Accepted clinical-significance strings in the current parser:

```text
Pathogenic
Likely_pathogenic
Pathogenic/Likely_pathogenic
Benign
Likely_benign
Benign/Likely_benign
```

Accepted strict review-status strings:

```text
criteria_provided,_multiple_submitters,_no_conflicts
reviewed_by_expert_panel
practice_guideline
```

These correspond to the intended 2-, 3-, and 4-star evidence categories. Exact values in the frozen VCF must still be audited before the final dataset is declared complete.

## 5. Data-source roles

- **ClinVar:** supplies the human clinical labels for Benchmark 1A.
- **Orthologous protein sequences:** supply the cross-species information used to construct MSAs and calculate Shannon entropy.
- **Established predictors:** planned features include SIFT, PolyPhen-2, CADD, and REVEL, subject to licensing, accessibility, versioning, coverage, and missing-data rules.
- **RNA-seq:** planned as a later evidence layer for M4. GEUVADIS has been identified as a candidate resource, but tissue relevance, expression, read depth, allele-specific analysis, and independence from model development must be handled explicitly.
- **EVA:** considered for later multi-species/input expansion; it is not the clinical label source for Benchmark 1A.

## 6. Repository and environment

Observed local repository path in Git Bash:

```text
/d/2. Research Work - Entropy Based App Tool/entropy-variant-prioritizer
```

Observed environment:

- Windows 10
- Git Bash (`MINGW64`)
- Python 3 virtual environment: `.venv`
- Activation command: `source .venv/Scripts/activate`

An earlier failure occurred because `python` resolved to Python 2.7.11 from MGLTools. The project was corrected to use Python 3.13.7 in a project-specific virtual environment.

Important scripts:

- `scripts/inspect_clinvar.py`: exploratory inspection of ClinVar VCF records and INFO fields.
- `scripts/build_clinvar_benchmark.py`: dataset-wide filtering counters for Benchmark 1A.

Important documentation:

- `README.md`
- `data/raw/README.md`
- `PROJECT_CONTEXT.md` (this file)
- `Entropy_Variant_Prioritization_Study_Notes_Through_ClinVar_Filtering.docx`

Raw ClinVar data must remain outside Git history. Its frozen filename, source, date, checksum when available, and retrieval instructions should be documented.

## 7. Completed technical workflow

### 7.1 ClinVar download verification

The first attempted file, `clinvar_20260816.vcf.gz`, was only 990 bytes and contained XHTML rather than a valid ClinVar VCF. It was rejected.

The corrected frozen file, `clinvar_20260822.vcf.gz`, was verified as a compressed VCF. The header included the expected columns:

```text
#CHROM POS ID REF ALT QUAL FILTER INFO
```

### 7.2 Exploratory VCF inspection

The compressed VCF was read without decompressing it permanently. The INFO column was parsed into key-value fields. An early inspected example was at chromosome 1, position 66926, ClinVar variation ID 3385321, with an uncertain-significance label and intron consequence. This confirmed the parser worked, but that example was not eligible for the missense binary benchmark.

### 7.3 Dataset-wide staged filtering

`scripts/build_clinvar_benchmark.py` processes the frozen ClinVar VCF sequentially and applies these filters:

1. Skip VCF metadata/header lines beginning with `#`.
2. Retain records whose `MC` annotation contains `missense_variant`.
3. Require `CLNVC=single_nucleotide_variant`.
4. Require both `REF` and `ALT` to contain exactly one nucleotide.
5. Retain only accepted pathogenic/likely-pathogenic or benign/likely-benign clinical-significance strings.
6. Retain only accepted strict 2- to 4-star review-status strings.
7. Count the binary classes after the strict review filter and verify that their sum equals the strict total.

The dual SNV check (`CLNVC` plus one-base `REF`/`ALT`) is intentional defensive validation.

## 8. Verified ClinVar filtering results

The following results were produced by:

```bash
python scripts/build_clinvar_benchmark.py
```

### 8.1 Attrition through the filters

| Filtering stage | Records |
|---|---:|
| Total ClinVar variants | 4,467,990 |
| Missense variants | 2,536,550 |
| Missense SNVs | 2,527,633 |
| Valid one-base REF/ALT SNVs | 2,527,633 |
| Accepted binary clinical significance | 223,911 |
| Strict 2- to 4-star candidates | 65,225 |

### 8.2 Class counts before strict review filtering

| Class | Records | Percentage |
|---|---:|---:|
| Pathogenic/Likely pathogenic | 71,362 | 31.87% |
| Benign/Likely benign | 152,549 | 68.13% |
| Total | 223,911 | 100.00% |

### 8.3 Class counts after strict review filtering

| Class | Records | Percentage |
|---|---:|---:|
| Pathogenic/Likely pathogenic | 24,904 | 38.18% |
| Benign/Likely benign | 40,321 | 61.82% |
| Total | 65,225 | 100.00% |

Verification:

```text
24,904 + 40,321 = 65,225
```

The benign-to-pathogenic ratio is approximately 1.62:1. The strict review filter retained approximately 34.90% of the pre-review pathogenic group, 26.43% of the pre-review benign group, and 29.13% of all binary-labeled records. The pathogenic proportion consequently increased from 31.87% to 38.18%.

This is a moderate class imbalance. The benchmark should retain the available records rather than randomly undersampling the benign group. PR-AUC, MCC, calibration, sensitivity, and specificity will complement ROC-AUC during evaluation.

## 9. Implementation completed on 2026-09-01

Three changes were made to `scripts/build_clinvar_benchmark.py`:

1. Added `high_review_pathogenic_variants` and `high_review_benign_variants` counters.
2. Incremented the appropriate counter only after a record passed the strict review-status filter, using the existing binary `label`.
3. Printed the two post-review class counts and a verification sum.

The verified script was committed and pushed:

```bash
git add scripts/build_clinvar_benchmark.py
git commit -m "Add strict ClinVar benchmark class counts"
git push origin main
```

Commit:

```text
fb96cd6 Add strict ClinVar benchmark class counts
```

Only `scripts/build_clinvar_benchmark.py` was included in that commit. At the last observed `git status`, these separate local changes remained intentionally uncommitted and must be reviewed independently:

```text
modified: README.md
modified: data/raw/README.md
modified: scripts/inspect_clinvar.py
untracked: .gitignore
```

Do not discard or overwrite these files without reviewing their contents.

## 10. Current scientific status

The 65,225 records are **strict candidates**, not yet the final publication-ready Benchmark 1A dataset. The class-count checkpoint is complete. Dataset-integrity QC must be completed before exporting the analytical table.

## 11. Exact next steps

### Immediate next step: integrity audit

Extend the benchmark workflow or add a dedicated QC script to:

1. Count unique variants using `(CHROM, POS, REF, ALT)`.
2. Detect duplicate genomic representations.
3. Test whether an identical genomic variant has both pathogenic and benign labels.
4. Count candidates separately by exact 2-star, 3-star, and 4-star review status.
5. Inspect exact `CLNSIG` and `CLNREVSTAT` values in the frozen release for edge cases.

### Following QC stages

1. Decide and document normalization rules for genomic representation.
2. Extract and normalize gene identifiers.
3. Select transcript policy and map genomic variants to transcripts/proteins.
4. Establish protein position and reference/alternate amino-acid consistency.
5. Define ortholog inclusion, taxonomic breadth, sequence-identity thresholds, and minimum usable sequence count.
6. Construct protein MSAs and calculate entropy plus alignment-quality features.
7. Add conventional functional-predictor scores with explicit version and missingness policies.
8. Create development and publication datasets without contaminating held-out evaluation.
9. Apply gene-held-out validation for M1-M4.
10. Add RNA-seq evidence only after its role and independence are clearly defined.

## 12. Methodological guardrails

- Do not call the 65,225 candidates the final dataset until integrity, mapping, and representation QC are complete.
- Do not use arbitrary hand-selected weights for the integrated score; derive and validate model parameters empirically.
- Do not randomly split variants across training and test sets when they share genes.
- Do not treat lack of an RNA effect in an irrelevant tissue as evidence of benignity.
- Do not use the same RNA-derived evidence for model construction and then describe it as independent validation.
- Do not silently update ClinVar or other databases; freeze and record data versions.
- Record stage counts and exclusion reasons so the manuscript can report dataset attrition transparently.
- Retain raw public data outside Git while tracking scripts, provenance documentation, configuration, and small derived summaries.

## 13. Update protocol for future sessions

After each completed milestone:

1. Run the relevant script and save the complete terminal output.
2. Verify numerical invariants, such as class totals equaling the filtered total.
3. Record the result and interpretation in this file.
4. Update the detailed Word study notes with commands, code logic, results, scientific purpose, and the new stopping point.
5. Review `git diff` and `git status` before staging.
6. Commit only logically related files with a descriptive message.
7. Push the verified checkpoint.
8. Record the commit hash and next action here.

## Benchmark 1A: ClinVar High-Confidence Missense SNVs

### Status

Benchmark construction and record-level dataset-integrity auditing are complete.

Final benchmark:

- 65,222 unique canonical missense SNVs
- 24,904 pathogenic/likely pathogenic variants
- 40,318 benign/likely benign variants
- Reference assembly: GRCh38
- ClinVar input release: 2026-08-22

### Frozen Input

```text
data/raw/clinvar_20260822.vcf.gz



