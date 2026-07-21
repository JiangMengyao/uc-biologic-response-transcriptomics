# Supplementary figure contracts

These contracts define the scientific job, evidence hierarchy and review risks
before plotting. All panels use the locked inflammatory-response signature and
the existing analysis outputs; no gene, cohort, endpoint or threshold is
reselected for display.

## Figure S1 | Cohort inclusion and evidence flow

- **Core conclusion:** Five public accessions yield seven active-treatment
  baseline cohorts, three randomized interaction comparisons and six active
  longitudinal cohorts after transparent sample and participant filtering.
- **Archetype:** schematic-led quantitative composite.
- **Final size:** 183 mm × 154 mm; editable SVG and PDF, 600-dpi TIFF and PNG
  preview.
- **Panel a:** accession-to-analysis flow, with cohort, mechanism and final
  participant counts.
- **Panel b:** mutually exclusive sample-disposition counts by cohort.
- **Panel c:** evidence-branch summary for baseline prognosis, randomized
  effect modification and active longitudinal change.
- **Statistics:** descriptive counts only; no P value or FDR.
- **Reviewer risk:** accession-level expression samples must not be confused
  with independent participants; repeated baseline biopsies are collapsed.

## Figure S2 | Baseline robustness and non-linearity audit

- **Core conclusion:** The adverse baseline association is directionally
  stable across cohorts and leave-one-out analyses, while cohort-specific
  curvature is heterogeneous and does not define a shared clinical threshold.
- **Archetype:** asymmetric quantitative grid with the cohort forest as hero.
- **Final size:** 183 mm × 240 mm.
- **Panel a:** complete seven-cohort Firth forest plot.
- **Panel b:** cohort-specific spline curves across the observed score range.
- **Panel c:** likelihood-ratio non-linearity P values and spline-minus-linear
  AIC differences.
- **Panel d:** leave-one-cohort-out pooled estimates.
- **Panel e:** complete prespecified module meta-analysis, including negative
  controls.
- **Statistics:** profile-likelihood 95% CIs for Firth estimates;
  REML/Hartung–Knapp summaries; likelihood-ratio tests for non-linearity.
- **Reviewer risk:** statistically detectable curvature in small cohorts must
  not be converted into a common threshold; non-linearity tests are shown in
  full and not FDR-selected.

## Figure S3 | Longitudinal sensitivity and module audit

- **Core conclusion:** Deeper reduction among responders remains
  directionally stable under bootstrap and mixed-model formulations and is
  concentrated in inflammatory rather than locked negative-control modules.
- **Archetype:** quantitative grid.
- **Final size:** 183 mm × 184 mm.
- **Panel a:** cohort-level responder-minus-non-responder change with analytic
  and 4,000-resample bootstrap intervals.
- **Panel b:** participant mixed-model time-by-response estimates.
- **Panel c:** cohort-by-module change-effect heatmap with all estimates.
- **Panel d:** module-level distribution and sign summary across cohorts.
- **Statistics:** analytic and bootstrap 95% intervals; Wald intervals for
  mixed models; all module estimates shown without discovery selection.
- **Reviewer risk:** longitudinal associations are observational and do not
  establish mediation or causal molecular normalization.

## Figure S4 | High-state transition and randomized paired-change audit

- **Core conclusion:** Threshold crossing from the high baseline tertile is
  common in both outcome groups and therefore not response-specific, while
  paired active–placebo contrasts remain imprecise.
- **Archetype:** schematic-led quantitative composite.
- **Final size:** 183 mm × 154 mm.
- **Panel a:** locked high-state and molecular-reversal definitions with total
  participant flow.
- **Panel b:** cohort- and response-specific reversal counts and rates.
- **Panel c:** exploratory reversal-response ORs, including non-estimable
  cohorts.
- **Panel d:** active-minus-placebo paired-change estimates and arm-specific
  mean changes in the two eligible RCT subsets.
- **Statistics:** Firth ORs where estimable and conventional 95% CIs for paired
  active–placebo change contrasts.
- **Reviewer risk:** the upper-tertile/median thresholds are descriptive and
  must not be presented as a validated classifier; null randomized contrasts
  are underpowered and do not establish treatment equivalence.

## Shared visual and export contract

- R-only rendering with `ggplot2`, `patchwork`, `svglite` and `ragg`.
- White background, Helvetica/Arial-compatible sans-serif text, lowercase bold
  panel labels and a restrained mechanism/outcome palette reused from Figures
  1–5.
- Main text remains readable at final double-column size; no oversized overall
  title is embedded in the figure.
- Every quantitative panel is written to a panel-tagged `FigureS*_source_data.csv`.
- Figure legends report analysis unit, cohort/participant counts, interval
  definition, model and interpretation boundary.
