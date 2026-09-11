# Stage 24 computational acceleration — 2026-09-10

This is a computational amendment, not a new statistical method. It preserves
seed 42, the CPM filter, gene order, all samples, design matrices, priors, both
voom passes and the five DE methods. The official limma 3.68.0 implementation
remains the reference. Other limma versions or unsupported inputs fall back to
the official implementation. No package namespace or installed library is edited.

## Bottleneck and equivalent algebra

The default gene-by-gene sample-weight update uses a sum-contrast matrix Z2 with
n-1 columns and repeatedly constructs and solves dense systems. For the default
sample-specific variance model, let h be one minus the leverage, b = h[-n]-h[n],
and H = sum(h). Its information increment is exactly
diag(h[-n]) + h[n] * 11' - bb'/H. Its score is z[-n]-z[n], and weights are
exp(c(-gamma, sum(gamma))). These identities avoid multiplication by dense
contrast matrices. The system is solved with preconditioned conjugate gradients,
using the accumulated diagonal-plus-rank-one term and the Woodbury identity.
The true residual is checked against 5e-12 * max(1, ||score||2); failure falls
back to the original dense solve. Statistical convergence criteria are unchanged.

Reference: [limma manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf),
arrayWeights and voomWithQualityWeights; installed limma 3.68.0 is authoritative
for this implementation comparison.

## Acceptance gates fixed before production use

- Compare against official gene-by-gene weights on deterministic synthetic
  datasets with 118, 564 and 842 samples, and paired/composition designs.
- Compare the full completed TCGA cache (22,787 genes, 564 samples): maximum
  absolute log weight ratio < 1e-7; absolute logFC and t differences < 1e-6;
  absolute p/FDR differences < 1e-7; zero changed DEG classifications and identical
  gene ranking. A failed gate prohibits enabling the accelerated production run.
- Validate input sample/gene identifiers and record hashes, timing and environment.
- Preserve the completed reference cache. Any migration must match its original
  complete input/design/metadata/package/seed/threshold key, never just its name.
- Save subsequent methods atomically so interruption does not discard completed
  methods. Record progress every 1,000 genes and per-method elapsed time.

Synthetic pilot (60 genes, seed 42; current laptop while the original run remains
active): 564 samples 34.58 s versus 0.69 s (~50x); 842 samples 90.19 s versus
1.51 s (~60x). These are pilot timings, not promises of whole-pipeline runtime.
The complete real-data equivalence check is still pending at this amendment.
