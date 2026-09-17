# Run from project directory
mkdir -p final_figures

cp figures/prepped_figures/signature_attributions.png final_figures/Figure1.png
cp figures/prepped_figures/variant_effects.png final_figures/Figure2.png
cp figures/prepped_figures/signature_effects.png final_figures/Figure3.png
cp figures/prepped_figures/gene_epistasis_plots.png final_figures/Figure4.png
cp figures/prepped_figures/pathway_effects.png final_figures/Figure5.png
cp figures/prepped_figures/MSK_landscape.png final_figures/Figure6.png

# Convert files as requested by journal
cd final_figures/
for file in $(ls Figure*png); do converted=${file/%png/tif}; converted=${converted/Figure/Fig}; tiffutil -lzw $file -o $converted; done
cd -


mkdir -p supplementary_materials
cp data_source_summary.txt supplementary_materials/S1_file.txt
cp figures/prepped_figures/combined_APOBEC_SBS12_signatures.png supplementary_materials/S2_fig.png
cp output/pathway/pathway_gene_membership_for_S3.pdf supplementary_materials/S3_text.pdf
cp output/survival/supp_survival_tables.pdf supplementary_materials/S4_tables.pdf

