library(data.table)
library(ggplot2)
library(glmnet)
library(survival)
library(gt)
library(gtsummary)
library(survminer)

# Need recent features
stopifnot(packageVersion('gtsummary') >= as.package_version('2.2'))

# Load pathway info and beautify for display
pw_info = fread('output/pathway/redefined_pathway_info.txt')
path_to_gene = fread('output/pathway/redefined_path_to_gene.txt')
needs_ellipses = path_to_gene[, .N, by = 'path_id'][N > 3, path_id]
pw_info[, display_dots := fcase(path_id %in% needs_ellipses, '…', default = '')]
pw_info[, genes := paste0(top_genes, display_dots)]
pw_info[, table_display_name := paste0(stringr::str_pad(paste0('P', pw_rank),
                                                        pad = ' ', width = 8, side = 'right'),
                                       '(', genes, ')')]
pw_info[is_cancer_gene_path == TRUE, table_display_name := sub('(^P\\d+)', '\\1\u2605', table_display_name)]

sample_key = fread('combined_sample_key.txt')
setnames(sample_key, 'Unique_Patient_Identifier', 'patient_id')

# Fit Cox models on given features/pathways and return formatted table. Adds pM (from sample key) to the feature table.
run_cox = function(features, path_ids, table_name = '', 
                   caption = 'Cox model of patient survival', signif_threshold = .1) {
  stopifnot(is.data.table(features),
            all(path_ids %in% names(features)))
  all_samples = features$patient_id
  path_ids = unique(path_ids)
  melted = melt(features[, .SD, .SDcols = c('patient_id', path_ids)], 
                id.vars = 'patient_id', variable.name = 'path_id')
  
  melted = melted[value != ''][, .SD[1], by = c('patient_id', 'path_id')] # not double-counting indel/SBS in same path
  melted[, value := 1]
  melted$patient_id = factor(melted$patient_id, levels = all_samples)
  fm = dcast(melted, patient_id ~ path_id, value.var = 'value', fill = 0)
  missing_samples = data.table(patient_id = setdiff(all_samples, fm$patient_id))
  missing_samples[, (path_ids) := 0]
  fm = rbind(fm, missing_samples)
  fm[, names(.SD) := lapply(.SD, as.numeric), .SDcols = setdiff(names(fm), 'patient_id')]
  fm[sample_key, let(surv_month = surv_month, surv_status = surv_status, pM = pM), on = 'patient_id']
  
  fm = fm[complete.cases(fm)]
  affected_counts = as.list(fm[, colSums(.SD), .SDcols = c(path_ids, 'pM')])
  cox1 = coxph(Surv(surv_month, surv_status) ~ ., data = fm[, -'patient_id'])
  
  stopifnot(cox1$n == fm[, .N])
  
  cox1_coef = as.data.table(summary(cox1)$coefficients, keep.rownames = 'feature')[order(`Pr(>|z|)`)]
  
  not_signif = cox1_coef[`Pr(>|z|)` > signif_threshold, feature]
  cox1_coef = cox1_coef[`Pr(>|z|)` <= signif_threshold]
  
  ordered_path_ids = setdiff(cox1_coef[order(-coef), feature], 'pM') # Order by desc HR
  
  table_labels = c(setNames(pw_info$table_display_name, pw_info$path_id),
                   list(pM = 'Presence of\nmetastatic disease'))
  
  # gene_labels = c(setNames(pw_info$genes, pw_info$path_id),
  #                 list(pM = '—'))
  
  curr_table_labels = table_labels[names(table_labels) %in% cox1_coef$feature]
  
  footer = paste0('<i>n</i> = ', cox1$n, ' patients (', cox1$nevent, ' events)')
  
  # Credit to https://stackoverflow.com/questions/65665465/grouping-rows-in-gtsummary
  # for how to get "Pathway mutation" as a group label
  
  table_name = paste0('<b>', table_name, '</b>')
  caption = paste0(caption, '<br><br>')
  cox_table = tbl_regression(cox1, label = curr_table_labels, exponentiate = TRUE, 
                             include = cox1_coef$feature, 
                             pvalue_fun = function(x) {
                               # Custom P value formatting function
                               formatted <- style_pvalue(x, digits = 2)
                               formatted <- gsub("<", " < ", formatted, fixed = TRUE)
                               # Add significance stars
                               stars <- fcase(
                                 x < 0.001, "***",
                                 x < 0.01, "**", 
                                 x < 0.05, "*",
                                 x >= 0.05, ""
                               )
                               return(paste0(formatted, stars))
                             }) |>
    modify_table_body(
      ~ .x |> 
        mutate(#genes = gene_labels[variable],
          num_affected = affected_counts[variable]) |>
        dplyr::relocate(c(num_affected), .before = estimate) |>
        dplyr::bind_rows(
          tibble::tibble(variable="mut_pw", var_label = "mut_pw", row_type="label",
                         label="Subpath mutation in...")) |> 
        dplyr::arrange(factor(variable, levels=c('pM', 'mut_pw', ordered_path_ids)))
    ) |> 
    # add_significance_stars(pattern = '{p.value}{stars}', hide_se = T, hide_ci = F, hide_p = F) |>
    modify_column_indent(columns=label, rows = variable %in% path_ids) |> 
    
    modify_header(label = '**Predictor**',
                  #genes = "**Frequently altered genes**",
                  estimate = '**Hazard ratio**',
                  num_affected = "**Number affected**",
                  p.value = "***P* value**"
    ) |>
    #modify_column_alignment(columns = 'genes', align = 'left') |> 
    remove_abbreviation() |> 
    as_gt() |>
    gt::tab_options(table.font.names = "Times New Roman",
                    footnotes.marks = c("a", "b", "c")) |>
    # Add sample size footnote
    gt::tab_footnote(
      footnote = gt::html(footer),
      locations = gt::cells_column_labels(columns = num_affected)
    ) |>
    # Add custom P value footnote
    gt::tab_footnote(
      footnote = gt::html("*<i>P</i> < 0.05; **<i>P</i> < 0.01; ***<i>P</i> < 0.001"),
      locations = gt::cells_column_labels(columns = p.value)
    ) |>
    gt::tab_header(title = gt::html(table_name),
                   subtitle = gt::html(caption)) |>
    gt::opt_align_table_header(align = 'left') |>
    gt::tab_options(heading.border.bottom.style = 'none',
                    heading.title.font.size = '120%',
                    heading.subtitle.font.size = '120%',
                    heading.border.lr.style = 'none',
                    table.border.top.style = 'none')
  return(list(gt = cox_table, fit = cox1))
}


# Prep feature table for LASSO 
prep_for_lasso = function(features, path_ids) {
  all_samples = features$patient_id
  melted = melt(features[, .SD, .SDcols = c('patient_id', path_ids)], 
                id.vars = 'patient_id', variable.name = 'path_id')
  melted = melted[value != ''][, .SD[1], by = c('patient_id', 'path_id')] # not double-counting indel/SBS in same path
  melted[, value := 1]
  melted$patient_id = factor(melted$patient_id, levels = all_samples)
  fm = dcast(melted, patient_id ~ path_id, value.var = 'value', fill = 0)
  missing_samples = data.table(patient_id = setdiff(all_samples, fm$patient_id))
  missing_samples[, (path_ids) := 0]
  fm = rbind(fm, missing_samples)
  fm[, names(.SD) := lapply(.SD, as.numeric), .SDcols = setdiff(names(fm), 'patient_id')]
  fm[sample_key, let(surv_month = surv_month, surv_status = surv_status, pM = pM), on = 'patient_id']
  fm = fm[complete.cases(fm)][surv_month > 0][]
  return(fm)
}

# Runs LASSO Cox
do_lasso = function(fm, seed = 12345, ...) {
  x = as.matrix(fm[, .SD, .SDcols = setdiff(names(fm), c('surv_month', 'surv_status', 'patient_id'))])
  y = as.matrix(fm[, .(time = surv_month, status = surv_status)])
  set.seed(seed)
  cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "C", cox.ties = 'breslow', ...)
  return(cvfit)
}

# Function to create a table for LASSO Cox regression results
plot_lasso_table = function(lasso_fit, features_data, pathway_info, sample_key, 
                            table_name = 'Table A', caption = '', lambda_choice = "lambda.1se") {
  
  table_name = paste0('<b>', table_name, '</b>')
  caption = paste0(caption, '<br><br>')
  
  # Extract coefficients at chosen lambda
  lasso_coef = as.data.table(
    as.matrix(coef(lasso_fit$glmnet.fit, s = lasso_fit[[lambda_choice]])),
    keep.rownames = TRUE
  )[`1` != 0]  # Only non-zero coefficients
  
  setnames(lasso_coef, c("feature", "coefficient"))
  lasso_coef[, hr := exp(coefficient)]
  
  # Get sample counts for each feature
  # Recreate the processed data to count affected samples
  all_path_ids = pathway_info$path_id
  processed_data = prep_for_lasso(features_data, path_ids = all_path_ids)
  
  # Calculate affected counts for each selected feature
  affected_counts_list = list()
  for(feat in lasso_coef$feature) {
    if(feat %in% names(processed_data)) {
      affected_counts_list[[feat]] = sum(processed_data[[feat]], na.rm = TRUE)
    } else {
      affected_counts_list[[feat]] = 0
    }
  }
  
  # Create table labels
  pathway_labels = setNames(pathway_info$table_display_name, pathway_info$path_id)
  clinical_labels = list(pM = 'Presence of\nmetastatic disease')
  all_labels = c(pathway_labels, clinical_labels)
  
  # Prepare data for table
  table_data = copy(lasso_coef)
  table_data[, `:=`(
    label = ifelse(feature %in% names(all_labels), all_labels[feature], feature),
    num_affected = as.numeric(affected_counts_list[feature]),
    hr_rounded = round(hr, 2) # 3 significant figures
  )]
  
  # Separate pathway mutations (order by HR) from clinical variables
  pathway_features = table_data[feature %in% pathway_info$path_id][order(-hr)]
  clinical_features = table_data[!feature %in% pathway_info$path_id]
  
  # Create final table structure
  # First clinical features, then pathway header, then pathway features
  final_data = list()
  
  # Add clinical features first
  if(nrow(clinical_features) > 0) {
    final_data = c(final_data, list(clinical_features))
  }
  
  # Add pathway header if there are pathway features
  if(nrow(pathway_features) > 0) {
    pathway_header = data.table(
      feature = "mut_pw",
      coefficient = NA_real_,
      hr = NA_real_,
      label = "Supbath mutation in...",
      num_affected = "",
      hr_rounded = ""
    )
    final_data = c(final_data, list(pathway_header), list(pathway_features))
  }
  
  # Combine all data
  final_table_data = rbindlist(final_data, fill = TRUE)
  
  # Convert to data frame for gt
  gt_input = final_table_data[, .(
    label = label,
    num_affected = num_affected,
    hr = hr_rounded
  )]
  
  # Replace NAs with empty strings for display
  gt_input[is.na(num_affected), num_affected := ""]
  gt_input[is.na(hr), hr := ""]
  # gt_input[label == "Subpath mutation in..." & hr_ci == "", hr_ci := ""]
  
  # Convert to data frame
  gt_df = as.data.frame(gt_input)
  
  # Create footer
  processed_data = prep_for_lasso(features_data, path_ids = pathway_info$path_id)
  n_patients = nrow(processed_data)
  n_events = sum(processed_data$surv_status, na.rm = TRUE)
  footer = paste0('<i>n</i> = ', n_patients, ' patients (', n_events, ' events)')
  
  # Create gt table
  gt_table = gt_df |>
    gt::gt() |>
    gt::cols_label(
      label = gt::html("<b>Predictor</b>"),
      num_affected = gt::html("<b>Number affected</b>"),
      hr = gt::html("<b>Hazard ratio</b>")
    ) |>
    # Left-align the first column (label)
    gt::cols_align(
      align = "left",
      columns = label
    ) |>
    # Style the pathway header row
    gt::tab_style(
      style = list(gt::cell_text(weight = "bold")),
      locations = gt::cells_body(
        columns = label,
        rows = label == "Subpath mutation in..."
      )
    ) |>
    # Indent pathway features
    gt::tab_style(
      style = list(gt::cell_text(indent = gt::px(20))),
      locations = gt::cells_body(
        columns = label,
        rows = !gt_df$label %in% c("Subpath mutation in...", "Presence of\nmetastatic disease") & gt_df$hr != ""
      )
    ) |>
    gt::tab_options(
      table.font.names = "Times New Roman",
      footnotes.marks = c("a", "b", "c")
    ) |>
    gt::tab_footnote(
      footnote = gt::html(footer),
      locations = gt::cells_column_labels(columns = num_affected)
    ) |>
    gt::tab_header(title = gt::html(table_name),
                   subtitle = gt::html(caption)) |>
    gt::opt_align_table_header(align = 'left') |>
    gt::tab_options(heading.border.bottom.style = 'none',
                    heading.title.font.size = '120%',
                    heading.subtitle.font.size = '120%',
                    heading.border.lr.style = 'none',
                    table.border.top.style = 'none')
  
  return(list(gt = gt_table, coefficients = lasso_coef))
}

# Run Cox model on full data set.
# Features are indel/SBS status (binarized to mutated/not mutated) in pathways, and pM status.
features = fread('output/landscape/pan_landscape_pw_plotted_features.txt')
cancer_pw = pw_info[is_cancer_gene_path == T, path_id]
cancer_pw_cox = run_cox(features, cancer_pw, table_name = 'Table C', 
                        caption = 'Cox model of pan-CCA patient survival using cancer-gene subpaths')

gtsave(cancer_pw_cox$gt, 'output/survival/cancer_gene_pw_cox_(supp3C).pdf')

summary(cancer_pw_cox$fit)$concordance
# C     se(C) 
# 0.6395650 0.0133303 


# Same model, using just the MSK samples.
# Unlike all other feature tables, non-neutral CNA status is included and counts as mutated.
## Probably a supplementary table (if used).
# msk_features = fread('output/landscape/msk_landscape_plotted_features.txt')
# msk_cox = run_cox(msk_features, cancer_pw)
# gtsave(msk_cox$gt, 'output/survival/ihc_msk_cox.pdf')
# summary(msk_cox$fit)$concordance
# C      se(C) 
# 0.72057295 0.01555276 


# Do LASSO Cox using exome pathways (which means, can only use WXS/WGS samples)
all_path_ids = pw_info$path_id

# Note that the feature table gives indel, SBS, or SBS+indel for each pathway (no CNA info).
# The lasso prep function binarizes the table to mutated/not mutated.
wxm = fread('output/landscape/nontarget_pan_landscape_plotted_features.txt') 

wxs_model = prep_for_lasso(wxm, path_ids = all_path_ids)
wxs_fit = do_lasso(wxs_model)

# View nonzero predictors
as.data.table(as.matrix(coef(wxs_fit$glmnet.fit, s= wxs_fit$lambda.1se)),
              keep.rownames = T)[`1` != 0][, .(predictor = rn, hr = exp(`1`))][pw_info, 
                                                                               pw_name := path_display_name,
                                                                               on = c(predictor = 'path_id')][]
# predictor       hr                pw_name
# <char>    <num>                 <char>
# 1: pw.2698.tgs 1.040716    KRAS|NRAS|RAC1 (P1)
# 2: pw.2227.tgs 1.144544 TP53|CTNNB1|EGFR (P14)
# 3:      pw.444 1.169877   FBN1|FN1|ITGB8 (P15)
# 4:          pM 1.375146                   <NA>

wxs_fit
# Call:  cv.glmnet(x = x, y = y, type.measure = "C", family = "cox") 
# 
# Measure: C-index 
# 
# Lambda Index Measure      SE Nonzero
# min 0.06947    10  0.5992 0.02458      11
# 1se 0.12139     4  0.5789 0.02155       4

pdf('output/survival/wxs_pw_lasso_path.pdf')
plot(wxs_fit)
dev.off()

wxs_table_result = plot_lasso_table(
  lasso_fit = wxs_fit,
  features_data = wxm,
  pathway_info = pw_info,
  sample_key = sample_key,
  table_name = 'Table D',
  caption = 'LASSO-regularized Cox model of pan-CCA patient survival using whole-exome and cancer-gene subpaths',
  lambda_choice = "lambda.1se"  # or "lambda.min"
)
gt::gtsave(wxs_table_result$gt, "./output/survival/wxs_pw_lasso_cox_table_(supp3D).pdf")

# LASSO on cancer gene pathways (not currently in use)
# for_cancer_pw_model = fread('output/landscape/pan_landscape_pw_plotted_features.txt')
# cancer_pw = pw_info[is_cancer_gene_path == T, path_id]
# cancer_pw_model = prep_for_lasso(for_cancer_pw_model, cancer_pw)
# cancer_fit = do_lasso(cancer_pw_model)
# 
# as.data.table(as.matrix(coef(cancer_fit$glmnet.fit, 
#                              s = cancer_fit$lambda.1se)),
#               keep.rownames = T)[`1` != 0][, .(predictor = rn, hr = exp(`1`))][pw_info, 
#                                                                                pw_name := path_display_name,
#                                                                                on = c(predictor = 'path_id')][]
# pdf('output/survival/cancer_pw_lasso_path.pdf')
# plot(cancer_fit)
# dev.off()

# predictor       hr                pw_name
# <char>    <num>                 <char>
# 1: pw.2227.tgs 1.180028 TP53|CTNNB1|EGFR (P14)
# 2:          pM 1.574628                   <NA>

### Subtype-specific survival analysis ###
ihc_pids = sample_key[cca_type == "IHC", patient_id]
ehc_pids = sample_key[cca_type %in% c("DCC", "PHC", "EHC"), patient_id]

# For combined eCCA (dCCA + pCCA + unspecified) only 19 samples with survival data, so no subtype modeling.
wxm_ehc = wxm[patient_id %in% ehc_pids]
wxs_model_ehc = prep_for_lasso(wxm_ehc, path_ids = all_path_ids)
stopifnot(wxs_model_ehc[, .N] == 19)

# IHC whole-exome model
wxm_ihc = wxm[patient_id %in% ihc_pids]
wxs_model_ihc = prep_for_lasso(wxm_ihc, path_ids = all_path_ids)
wxs_fit_ihc = do_lasso(wxs_model_ihc)
# View nonzero predictors
as.data.table(as.matrix(coef(wxs_fit_ihc$glmnet.fit, s= wxs_fit_ihc$lambda.1se)),
              keep.rownames = T)[`1` != 0][, .(predictor = rn, hr = exp(`1`))][pw_info, 
                                                                               pw_name := path_display_name,
                                                                               on = c(predictor = 'path_id')][]

# 1se concordance ("Measure") stated in manuscript
wxs_fit_ihc
# Call:  cv.glmnet(x = x, y = y, type.measure = "C", family = "cox") 
# 
# Measure: C-index 
# 
# Lambda Index Measure      SE Nonzero
# min 0.04967    14  0.5981 0.01453      19
# 1se 0.07908     9  0.5895 0.01693      11

pdf('output/survival/wxs_pw_lasso_path_ihc.pdf')
plot(wxs_fit_ihc)
dev.off()

wxs_table_result_ihc = plot_lasso_table(
  lasso_fit = wxs_fit_ihc,
  features_data = wxm_ihc,
  pathway_info = pw_info,
  sample_key = sample_key,
  table_name = 'Table B',
  caption = 'LASSO-regularized Cox model of iCCA patient survival using whole-exome and cancer-gene subpaths',
  lambda_choice = "lambda.1se"
)
gt::gtsave(wxs_table_result_ihc$gt, "./output/survival/wxs_pw_lasso_cox_table_ihc_(supp3B).pdf")


cancer_pw_cox_ihc = run_cox(features[cca_subtype == 'Intrahepatic'], path_ids = cancer_pw,
                            table_name = 'Table A', caption = 'Cox model of iCCA patient survival using cancer-gene subpaths')
gtsave(cancer_pw_cox_ihc$gt, 'output/survival/cancer_gene_pw_cox_ihc_(supp3A).pdf')

# Concordance stated in manuscript
summary(cancer_pw_cox_ihc$fit)$concordance
# C     se(C) 
# 0.6487932 0.0141930 


model_tables = c('output/survival/cancer_gene_pw_cox_ihc_(supp3A).pdf',
                 'output/survival/wxs_pw_lasso_cox_table_ihc_(supp3B).pdf',
                 'output/survival/cancer_gene_pw_cox_(supp3C).pdf',
                 'output/survival/wxs_pw_lasso_cox_table_(supp3D).pdf')

# Make one PDF
qpdf::pdf_combine(input = model_tables, output = 'output/survival/supp_survival_tables.pdf')

