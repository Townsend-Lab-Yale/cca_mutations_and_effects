# Use transcript definitions (i.e., ces.refset.hg38::transcripts) to
# come up with some information for the given noncoding variants.
annotate_noncoding = function(cesa, variant_ids, transcripts) {
  transcripts = copy(transcripts)
  variant_ids = unique(variant_ids)
  coding_ends = transcripts[type == 'CDS', .(coding_end = max(start, end),
                                             coding_start = min(start, end)), by = 'transcript_id']
  transcripts[coding_ends, let(coding_start = coding_start, coding_end = coding_end), on = 'transcript_id']
  transcripts[type == 'UTR' & strand == '+', utr_type := fcase(end < coding_start, "5'UTR", default = "3'UTR")]
  transcripts[type == 'UTR' & strand == '-', utr_type := fcase(start > coding_end, "5'UTR", default = "3'UTR")]
  transcripts[type == 'UTR', type := utr_type]
  
  variants = cesa$variants[variant_ids, on = 'variant_id', nomatch = NULL][, .(variant_id, chr, start, end)]
  
  if(variants[, .N] != length(variant_ids)) {
    stop('Not all variants found.')
  }
  if(! all(variants$variant_type == 'snv')) {
    stop('Variant IDs should be noncoding')
  }
  setkey(transcripts, chr, start, end)
  setkey(variants, chr, start, end)
  ol = foverlaps(variants, transcripts)
  no_match = ol[is.na(gene_name), .(variant_id = variant_id)]
  no_match[, anno := 'intergenic']
  ol = ol[! is.na(gene_name)]
  if(any(ol$type == 'CDS')) {
    stop('Found coding annotation.')
  }
  top_anno = ol[order(! type %in% c("5'UTR", "3'UTR"), -is_mane), .SD[1], by = 'variant_id']
  top_anno[, anno := fcase(gene_type != 'protein_coding', paste0(gene_name, ' (', gene_type, ')'),
                           type == "5'UTR", paste0(gene_name, ":5'UTR"),
                           type == "3'UTR", paste0(gene_name, ":3'UTR"),
                           type == 'transcript', paste0(gene_name, ':intron'))]
  
  
  top_anno[, bare_anno := fcase(gene_type != 'protein_coding', gene_type,
                                type == "5'UTR", "5'UTR",
                                type == "3'UTR", "3'UTR",
                                type == 'transcript', 'intronic')]
  output = rbind(top_anno[, .(variant_id, anno, gene = gene_name, bare_anno)],
                 no_match[, .(variant_id, anno, gene = 'intergenic', bare_anno = '')])
  output = output[variants, .(variant_id, variant_name = sub('_', ' ', variant_id), 
                              gene, anno, bare_anno), on = 'variant_id']
  return(output)
}