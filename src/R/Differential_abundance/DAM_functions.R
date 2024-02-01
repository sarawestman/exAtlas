# Load functions 
diff_ttest <- function(met, group_factors, test_gr) {
  new_meta <- met@colData  %>% as.data.frame() %>% filter(get(group_factors) %in% test_gr)
  new_meta[[group_factors]] <- factor(new_meta[[group_factors]])
  new_assay <- assays(met)[["norm_imputed"]] %>% as.data.frame() %>% select(all_of(rownames(new_meta))) %>% as.matrix()
  metadata(met) = vector("list", 0)
  for (i in 1:length(group_factors)) {
    if (length(levels(as.factor(new_meta[[group_factors[i]]]))) > 2) {
      coeff = sapply(1:nrow(new_assay), function(x) aov(new_assay[x, ] ~ as.factor(new_meta[[group_factors[i]]]))$coefficient)
      xlevels = unlist(aov(new_assay[1, ] ~ as.factor(new_meta[[group_factors[i]]]))$xlevels)
      res = sapply(1:nrow(new_assay), function(x) summary(aov(new_assay[x, ] ~ as.factor(new_meta[[group_factors[i]]]))))
      res_df = data.frame(pval = as.vector(sapply(sapply(res, "[", i = 5), "[", i = 1)), 
                          adj_pval = p.adjust(as.vector(sapply(sapply(res, "[", i = 5), "[", i = 1)), method = "fdr"), 
                          dm = coeff[2, ])
      metadata(met)[[paste0("anova_", group_factors[i], "_", paste(xlevels, collapse = "_vs_"))]] = res_df
    }
    else {
      xlev = levels(as.factor(new_meta[[group_factors[i]]]))
      df = genefilter::rowttests(new_assay, 
                                 fac = as.factor(new_meta[[group_factors[i]]]))
      res_df = data.frame(metabolite = rownames(df), pval = df$p.value, 
                          adj_pval = p.adjust(df$p.value, method = "fdr"), 
                          dm = df$dm)
      metadata(met)[[paste0("ttest_", group_factors[i], "_", xlev[2], "_vs_", xlev[1])]] = res_df
    }
  }
  metadata(met)
}