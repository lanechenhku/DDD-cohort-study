# HK Spine Cohort: lumbar disc degeneration (LDD) analyses

Analysis code for two manuscripts from the Hong Kong Spine Cohort (2026 re-analysis):

- **Paper 1** – *Onset, progression and correlates of lumbar disc degeneration* (longitudinal MRI cohort, three lumbar scans; plus the whole-spine MRI section, Figures 6–8 and Table 4).
- **Paper 2** – *GWAS of LDD phenotypes and the ITGA9 locus (rs1549596) for Modic change.*

The repository holds the scripts exactly as they were run to produce the numbers, tables and figures in the
2026-11-07 drafts, together with **aggregate** result tables and the rendered figures.

> **No participant-level data are in this repository.** MRI gradings, questionnaire data, genotypes, imputed
> onset ages and any file with one row per participant (or per participant × scan) stay on the study's
> secure storage. `.gitignore` blocks the known file names and folders; see *Data* below.

## Layout

The folder layout mirrors the project working folders the scripts were written against, so relative imports
(`style.py`, `common_v2.py`, `common.py`) keep working.

```
analyses_v2/            Papers 1 and 2: re-analyses (v2 2026-10-23, v3 2026-11-06/07)
  code/                 p1_* = Paper 1, p2_* = Paper 2, shared helpers, build/ (docx/pdf build)
  results/              aggregate result CSV/JSON (model coefficients, curves, summaries)
  results/tables_md/    supplementary/main tables rendered as Markdown
analyses_ws/            Paper 1 whole-spine section (2026-11-06)
  code/
  results/              aggregate results; _tv/ = sensitivity run excluding transitional vertebrae
  results/tables_md/
figures_v2/             Redrawn figures (2026-10-20, updated 2026-11-06)
  code/                 figure scripts
  data/                 digitised / derived aggregate curve tables
  Paper1/  Paper2/      PNG + PDF figures (TIFF masters kept locally only)
```

## Running the code

Scripts were run in a sandbox where the inputs were staged at fixed paths. To re-run, point these at your
local copies:

| Path in scripts | What it is | Local equivalent |
|---|---|---|
| `/mnt/user-data/uploads/DDD_Paper` | Lane's `DDD_Paper` folder (raw cohort files) | `~/Desktop/DDD_Paper` |
| `/mnt/project-files/DDD_papers/analyses_v2/data_clean/DDD202008` | cleaned grading files (69 undated T2 records removed) | `Analyses_v2_2026-10-23/data_clean/DDD202008` |
| `/mnt/project-files/DDD_papers/analyses_v2/results` | results folder (also holds participant-level intermediates) | `Analyses_v2_2026-10-23/results` |
| `/mnt/project-files/DDD_papers/...` (other) | this repository's folders | repo root |

`figures_v2/code/common.py` resolves paths relative to the repo and honours `DDD_ROOT` (the `DDD_Paper` folder)
and `DDD_CLEAN` environment variables.

Software: Python 3.13 (pandas, numpy, scipy, statsmodels, scikit-learn, matplotlib, Pillow, pgmpy, pysam,
scanpy + harmonypy for the public scRNA-seq); R 4.x (geepack, SISE from github.com/tubbsjd/SISE, coloc, lavaan, lme4).

### Order of execution (Paper 1)

1. `analyses_v2/code/p1_exclude_undated_t2.py` → cleaned grading files (`data_clean/`, participant-level, not in repo)
2. `p1_build_long.py` → `results/p1_long.csv`, `p1_person.csv` (participant-level, not in repo)
3. `p1_descriptive.py`, `p1_turnbull.py`, `p1_sise_intervals.py` → `p1_sise.R`, `p1_interrater.py`, `p1_levelwise_transitions.py`, `p1_dbn.py`
4. `p1_wgee.R`, `p1_ipw_mri.py` → `p1_wgee_ipwmri.R`, `p1_wgee_m23_ageheight.R`, `p1_pain.py`, `p1_bmi_height_model_check.R`
5. Whole spine: `analyses_ws/code/ws_clean.py` → `ws_link_prev.py` → `ws_corr_decay.py`, `ws_factor.R`, `ws_models.py`, `ws_reader_agreement.py`, `ws_sensitivity_tv.py` → `ws_tables.py`, `ws_figs.py`
6. `make_tables.py` (Markdown tables), `figs_v2.py` and `figures_v2/code/*` (figures)

## Script reference

### analyses_v2/code – shared

| Script | Purpose | Main outputs |
|---|---|---|
| `common_v2.py` | Input/output paths; prefers cleaned grading files | – |
| `style.py` | Shared figure style and feature palette | – |
| `figs_v2.py` | v2 figures from re-analysis results (Paper 1 Fig 4, 5, S5, S6; Paper 2 Fig S3) | `figures_v2/Paper1`, `Paper2` |
| `make_tables.py` | Renders result CSVs as Markdown tables | `results/tables_md/*.md` |
| `refcheck.py` | Checks numbered references against PubMed (E-utilities) | console / `refcheck/` |
| `renumber_refs.py` | Renumbers Vancouver citations by first appearance | edited .md |
| `build/build_doc.py` | Builds docx + pdf from the living .md drafts, TODOs highlighted | `docx/`, `pdf/` |

### analyses_v2/code – Paper 1

| Script | Purpose | Main outputs (aggregate unless noted) |
|---|---|---|
| `p1_exclude_undated_t2.py` | Drops T2 gradings of participants with no T2 scan date (v3 decision) | `data_clean/` *(participant-level)* |
| `p1_build_long.py` | Builds person × scan file; SF-36, LBP; attendance (IPW) weights | `p1_long.csv`, `p1_person.csv` *(participant-level)* |
| `p1_descriptive.py` | Attrition, regression among affected, scanner effects, level-wise summary, Turnbull onset | `p1_attrition_*.csv`, `p1_transitions.csv`, `p1_scanner_effects.csv`, `p1_tableS2_levelwise.csv`, `p1_misc.txt` |
| `p1_turnbull.py` | Turnbull NPMLE of age at onset per feature (check on SISE) | `p1_turnbull_onset.csv`, `p1_turnbull_curves.csv` |
| `p1_sise_intervals.py` | Onset intervals per participant for SISE | `p1_sise_intervals.csv` *(participant-level)* |
| `p1_sise.R` | SISE onset curves (logNe main; logNm, logN robustness) | `p1_sise_curves.csv`, `p1_sise_summary.csv`; `p1_sise_*.rds`, `p1_sise_imputed_*.csv` *(participant-level)* |
| `p1_figS2_onset_hist.py` | Fig S2: distribution of imputed onset ages | `FigS2_imputed_onset_ages` |
| `p1_interrater.py` | Inter-rater agreement (Table S1), raters anonymised | `p1_interrater_*.csv` |
| `p1_levelwise_transitions.py` | Disc-level transitions between scans (Table S7) | `p1_levelwise_transitions.csv` |
| `p1_dbn.py` | Dynamic Bayesian networks (feature and disc level), bootstrap arc strengths | `p1_dbn_feature.csv`, `p1_dbn_level.csv` |
| `p1_wgee.R` | WGEE risk-factor Models 1–3 + sensitivity (geepack) | `p1_wgee_model1/2/3.csv` |
| `p1_ipw_mri.py` | Alternative attendance weights using sex and prior MRI | weights *(participant-level)* |
| `p1_wgee_ipwmri.R` | Model 1 with MRI-informed weights (appended to model1) | `p1_wgee_model1.csv` |
| `p1_wgee_m23_ageheight.R` | Models 2–3 keeping sport and age × height | `p1_wgee_m23_keep_sport_ageheight.csv` |
| `p1_pain.py` | MRI vs pain, disability, mental health (cross-sectional and prospective) | `p1_pain_*.csv` |
| `p1_bmi_height_model_check.R` | Re-check of the 2020 BMI/height-at-40 predictions | `p1_bmi_height_model_check.csv` |

### analyses_v2/code – Paper 2

| Script | Purpose | Main outputs |
|---|---|---|
| `p2_itga9_sensitivity.py` | rs1549596 sensitivity analyses on hard-called genotypes; phenotype descriptives | `p2_rs1549596_sensitivity.csv`, `_by_level.csv`, `_genotype_means.csv`, `p2_phenotype_*.csv`, `p2_misc.txt` |
| `p2_itga9_perm.py` | Freedman–Lane permutation p (1e6) and robust count models | `p2_misc.txt` |
| `p2_figS_gwas_pca.py` | Fig S1b: genotype PCA of 2,121 QC-passed participants | `FigS1b_genotype_PCA`, `p2_gwas_pca_summary.txt` |
| `p2_fig2_data.py` | MC prevalence by rs1549596 carrier status (Fig 2) | `p2_fig2_mc_prevalence_by_carrier.csv` |
| `p2_figs_v3_redesign.py` | Redesigned Figures 2–3 and S7 | `Fig2_rs1549596_MC_in_cohort`, `Fig3_ITGA9_context_mechanism`, `FigS7_ITGA9_eQTL_coloc` |
| `p2_itga9_gtex.py` | ITGA9 median expression across GTEx v8 tissues (Fig 2A) | `p2_ITGA9_GTEx_v8_median.json`, `Fig2A_ITGA9_GTEx` |
| `p2_itga9_coloc.R` | coloc of MC GWAS signal with ITGA9 eQTLs (eQTL Catalogue) | `p2_itga9_coloc.csv` |
| `p2_itga9_expr.py` | ITGA9 single-cell expression (HPA v24) | `p2_ITGA9_HPA_singlecell.tsv` |
| `p2_np_scrna_load.py`, `p2_np_scrna_process.py` | Public NP scRNA-seq (Tu et al. 2022, GSE165722): load, QC, Harmony, UMAP | cell table (not in repo, public source) |
| `p2_figS4_S6_public.py` | Figs S4 and S6 from public data only (expression, chromatin) | `FigS4_ITGA9_expression_public`, `FigS6_rs1549596_chromatin_public`, `p2_s6_ucsc/` |
| `p2_phewas.py` | rs1549596 PheWAS in FinnGen R12 and BioBank Japan | `p2_rs1549596_phewas.csv`, `fg_phewas.json`, `bbj_phewas.json`, `FigS5` |
| `p2_go2_lookup.py` | Lookup in GO2 (Genetics of Osteoarthritis) summary stats | `p2_go2_rs1549596.csv`, `p2_go2_itga9_region.tsv` |
| `p2_spine_lookup_finngen.py`, `p2_spine_lookup_gwascatalog.py` | Spine-trait lookups for rs1549596 | `p2_rs1549596_spine_lookups.csv` |
| `p2_nearest_genes.py` | Nearest genes for suggestive loci | `p2_suggestive_nearest_genes.json` |

### analyses_ws/code – whole spine (Paper 1, Figures 6–8, Table 4, Tables S12–S17)

| Script | Purpose | Main outputs |
|---|---|---|
| `ws_clean.py` | Cleans and reshapes the whole-spine gradings | `ws_long.csv`, `ws_wide.csv` *(participant-level)*, `ws_clean_log.csv` |
| `ws_link_prev.py` | Links to lumbar cohort; lumbar grading agreement; prevalence atlas | `ws_linkage_log.csv`, `ws_lumbar_agreement.csv`, `ws_prevalence_by_*.csv` |
| `ws_corr_decay.py` | Between-level partial correlations by distance and region | `ws_corr_pairs.csv`, `ws_corr_by_lag.csv`, `ws_corr_by_block.csv` |
| `ws_factor.R` | Ordinal CFA (one factor, three region factors, bifactor; MIMIC) | `ws_factor_*.csv`, `ws_factor.log` |
| `ws_models.py` | Risk factors by region (GEE); lumbar course vs C/T degeneration | `ws_riskfactors_*.csv`, `ws_lumbar_course_vs_regions.csv`, `ws_symptoms_by_region.csv`; `ws_regional_scores.csv` *(participant-level)* |
| `ws_reader_agreement.py` | Whole-spine reader vs own earlier lumbar grading | `ws_reader_agreement.csv` |
| `ws_sensitivity_tv.py` | Re-runs excluding transitional lumbosacral vertebrae | `results/_tv/`, `ws_sens_tv_*.csv` |
| `ws_tables.py` | Markdown tables (Table 4, S12–S17) | `results/tables_md/` |
| `ws_figs.py` | Figures 6–8 | `Fig6_whole_spine_atlas`, `Fig7_spatial_structure`, `Fig8_regional_risk_factors` |

### figures_v2/code

| Script | Figure |
|---|---|
| `common.py`, `style.py`, `derive.py` | paths, style, derived composite-state tables |
| `p1_fig1_design.py` | Paper 1 Fig 1 cohort design |
| `p1_fig2_features.py` | Paper 1 Fig 2 LDD features (A–D) (`_origDBN_2026-10-20` = earlier version) |
| `p1_fig3_states.py` | Paper 1 Fig 3 states and transitions (+ `Fig3A/B_*.csv`) |
| `p1_fig4_wgee.py` | Paper 1 Fig 4 WGEE (first version from digitised coefficients; superseded by `figs_v2.py`) |
| `p1_supp.py` | Paper 1 Figs S1, S3, S4 |
| `p1_figS7_mri_examples.py` (`p1_figS9_…` = same, earlier numbering) | Paper 1 Fig S7 MRI examples (figure itself not in repo) |
| `p2_gwas.py`, `prep_sumstats.py` | Paper 2 Fig 1 Manhattan/locus, Figs S1 QQ, S2 Manhattan (thinned sumstats) |

`figures_v2/README_2026-10-20.md` documents the original redraw and the sources of digitised values.

## Data

Not in this repository (participant-level; held on study storage):
`data_clean/`, `analyses_v2/data/` (regional GWAS summary statistics and locus genotypes), `p1_long.csv`,
`p1_person.csv`, `p1_imputed_onset_ages.csv`, `p1_sise_imputed_*.csv`, `p1_sise_intervals.csv`, `p1_sise_*.rds`,
`ws_long.csv`, `ws_wide.csv`, `ws_regional_scores.csv`, genotype files (`.bed/.bim/.fam`, `.bgen`, `.vcf`),
raw grading spreadsheets and questionnaires, MRI images.

Also left out: TIFF masters, the MRI-example figure (Fig S7), superseded figures, and the cell-level table
from the public Tu et al. 2022 scRNA-seq data (re-creatable with `p2_np_scrna_*.py`).

Access to HK Spine Cohort data is by application to the study team.
