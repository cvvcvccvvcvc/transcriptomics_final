## Main submission files

- Concise report: `report/assignment_report.md`
- Full reproducible R entry point: `scripts/run_all.R`

## Analysis scripts

- Raw import, normalization, PCA/QC, gene-level matrix construction: `scripts/run_analysis.R`
- Final batch-adjusted DEG model, enrichment, pathway figures: `scripts/finalize_model_and_enrichment.R`

## Key result files

- Final DEG table: `results/deg_all.tsv`
- Upregulated DEGs: `results/deg_up.tsv`
- Downregulated DEGs: `results/deg_down.tsv`
- Model sensitivity summary: `results/model_sensitivity.tsv`
- Hallmark enrichment: `results/enrichment_hallmark.tsv`
- Reactome enrichment: `results/enrichment_reactome.tsv`
- KEGG enrichment: `results/enrichment_kegg_legacy.tsv`
- GO BP enrichment: `results/enrichment_go_bp.tsv`

## Key figures

- PCA before normalization: `figures/pca_raw_gene_level.pdf`
- PCA after normalization: `figures/pca_normalized_gene_level.pdf`
- Volcano plot: `figures/volcano_deg.pdf`
- Raw-to-processed correlation heatmap: `figures/raw_to_geo_processed_correlation_heatmap.pdf`
- Model robustness heatmaps: `figures/model_logFC_correlation_heatmap.pdf`, `figures/model_t_correlation_heatmap.pdf`, `figures/model_deg_jaccard_heatmap.pdf`
- Pathway stability across models: `figures/pathway_stability_across_models.pdf`
- Pathway barcode plot: `figures/pathway_barcode_translation_initiation.pdf`
- Pathway heatmap: `figures/pathway_heatmap_translation_initiation.pdf`

## Re-run

From the project root:

```bash
Rscript scripts/run_all.R
```
