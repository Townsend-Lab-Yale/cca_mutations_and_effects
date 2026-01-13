# Run from project directory
mkdir -p final_figures

cp data_prep/beautified_data_source_summary.pdf final_figures/Table1.pdf

cp figures/prepped_figures/signature_attributions.png final_figures/Figure1.png
cp figures/prepped_figures/variant_effects.png final_figures/Figure2.png
cp figures/prepped_figures/signature_effects.png final_figures/Figure3.png
cp figures/prepped_figures/gene_epistasis_plots.png final_figures/Figure4.png
cp figures/prepped_figures/pathway_effects.png final_figures/Figure5.png
cp figures/prepped_figures/MSK_landscape.png final_figures/Figure6.png

cp figures/prepped_figures/combined_APOBEC_SBS12_signatures.png final_figures/SuppFigureS1.png
cp output/survival/subtype_survival/cancer_gene_pw_cox_ihc.pdf final_figures/SuppTableS1.pdf
cp output/survival/subtype_survival/wxs_pw_lasso_cox_table_ihc.pdf final_figures/SuppTableS2.pdf
cp output/survival/cancer_gene_pw_cox.pdf final_figures/SuppTableS3.pdf
cp output/survival/wxs_pw_lasso_cox_table.pdf final_figures/SuppTableS4.pdf

