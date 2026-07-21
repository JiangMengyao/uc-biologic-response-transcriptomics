# Analysis execution contract — cross-mechanism prognostic state

This contract supersedes the anti-TNF-specific branch. It was locked before
running longitudinal or single-cell analyses. The primary claim is deliberately
prognostic: a high baseline mucosal inflammatory state is associated with lower
response across biologic mechanisms. A drug-specific predictive claim requires
a treatment-by-score interaction and is never inferred from a treated arm alone.

## Locked primary score and endpoints

- Primary score: `HALLMARK_INFLAMMATORY_RESPONSE`, scored by GSVA and
  standardized to each accession's baseline UC distribution.
- Primary binary endpoints: each study's prespecified clinical response or
  remission endpoint. For GSE206285 the primary endpoint is Week 8 clinical
  remission; mucosal healing is secondary.
- Primary baseline estimand: response odds ratio per 1-SD higher score within
  each treated cohort, using Firth logistic regression.
- Primary longitudinal estimand: follow-up-minus-baseline score difference in
  responders versus nonresponders among paired treated participants. Negative
  values indicate greater molecular resolution in responders.
- TNF–NFκB and IL6–JAK–STAT3 are supportive positive-control modules;
  oxidative phosphorylation and the locked random set are specificity controls.

## Locked evidence gates

1. **Cross-mechanism baseline gate.** The primary score must have an odds ratio
   below 1 in the anti-TNF synthesis, the UST active arm, and the VDZ synthesis.
   Statistical significance is not required in every small cohort. The overall
   REML/Hartung–Knapp synthesis is supportive because endpoints differ.
2. **Prognostic-versus-predictive gate.** RCTs are shown as arm-specific
   prognostic effects plus treatment-by-score interactions. No interaction is
   required for the common-prognostic-state hypothesis. A predictive or
   drug-specific claim is allowed only if the relevant interaction supports it.
3. **Longitudinal direction gate.** The responder-minus-nonresponder change must
   be below 0 in at least two independent treated cohorts and may not show a
   coherent opposite direction in a mechanism class. Failure stops the workflow
   before single-cell localization and final figure production.
4. **High-state transition analysis.** “High” is locked as the top baseline
   tertile within each cohort; “molecular reversal” means the follow-up score is
   below that cohort's baseline median. This is descriptive and is not used to
   optimize a clinical cut-point.
5. **Single-cell localization.** GSE282122 colonic epithelial, myeloid and
   fibroblast/pericyte h5ad files are analysed together. This is orthogonal
   localization, not causal proof. UC donors are the inferential unit;
   cell-level p-values are prohibited. Cell-state localization is reported as
   observed without changing the frozen bulk signature.

The single-cell implementation is further locked by
`docs/analysis_amendment_2026-07-20_singlecell_sample_first.md`. Samples are the
first aggregation unit and donors are the inferential unit. Primary localization
and the baseline clinical audit use baseline inflamed biopsies. Longitudinal
audits require publisher `Match=Yes` and exact donor-site Pre/Post matching.
The pre-amendment state/cell-weighted epithelial result is retained as a labelled
sensitivity analysis; it is not used to reselect genes or states.

## Multiplicity and robustness

- The inflammatory-response score and clinical endpoint are primary.
- Secondary modules, nonlinear splines, mucosal healing, late time points and
  transition rates are explicitly supportive/exploratory.
- Profile-likelihood Firth intervals are retained for sparse binary outcomes.
- Meta-analysis reports cohort effects, REML/Hartung–Knapp intervals,
  heterogeneity, and leave-one-out sensitivity; mechanism subgroup summaries
  are descriptive when fewer than three independent datasets are available.
- Samples are collapsed to one subject-time record before paired analysis.

## Figure contract

- Backend: R only for plotting, rendering, export and visual QA.
- Archetypes: schematic evidence map; forest/ridgeline/nonlinear quantitative
  grid; prognostic-versus-predictive RCT grid; longitudinal triptych with raw
  pairs, effect forest and high-state transitions; donor-aware single-cell
  cross-compartment localization composite. Figure 5 is localization-only;
  single-cell clinical outcome and longitudinal audits are supplementary.
- Final size: 183 mm wide; SVG/PDF with editable text; TIFF at 600 dpi; PNG
  previews; one source-data CSV per main figure.
- Visual stance: dense but legible layering, mechanism-specific gradients,
  ribbons and direct labels. Raw observations, cohort n, effect sizes, 95%
  intervals and provenance remain visible.
