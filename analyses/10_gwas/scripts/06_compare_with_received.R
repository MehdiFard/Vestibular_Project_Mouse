#!/usr/bin/env Rscript

### 06 Check against the GWAS received from UCLA (input/gwas_ucla_results, copied from ../../data/gwas_results; used by 02/04): SNPs matched by chr:bp (mm10)
# Out: results/compare_with_received.tsv (per trait and model: SNPs, correlation of -log10 P, inflation, top hits),
#      figures/compare_with_received.jpeg (-log10 P, ours vs received; both models)
# Usage (from analyses/10_gwas): Rscript scripts/06_compare_with_received.R <number of PCs in withPCs>
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

received <- c(threshold = "Young_VsEP_threshold.gwas", latency = "Young_VsEP_P2_latency.gwas",
              amplitude = "Young_VsEP_P1_Amplitude.rsID.gwas", beamtest = "Young_beamtest.gwas")
gif <- function(p) median(qchisq(p, 1, lower.tail = F)) / qchisq(0.5, 1)
top_n <- 100
n_pcs <- commandArgs(trailingOnly = TRUE)[1]

rows <- list(); pts <- list()
for (t in names(received)) {
  rec <- read.table(file.path("input/gwas_ucla_results", received[[t]]), header = T, sep = "\t") %>%
    transmute(key = paste(sub("^X$", "20", CHR), BP, sep = ":"), P_rec = P)   # chr X is "X" or 20 in the received files
  for (model in c("noPCs", "withPCs")) {
    ours <- read.table(paste0("results/", model, "/", t, ".gwas.loco.tsv"), header = T, sep = "\t") %>%
      transmute(key = paste(Chr, ChrPos, sep = ":"), P_ours = PValue)
    m <- inner_join(ours, rec, by = "key")
    top_ours <- m$key[order(m$P_ours)][1:top_n]; top_rec <- m$key[order(m$P_rec)][1:top_n]
    rows[[length(rows) + 1]] <- tibble(
      trait = t, model = ifelse(model == "noPCs", model, paste0("withPCs_PC1-", n_pcs)), received = received[[t]],
      snps_ours = nrow(ours), snps_received = nrow(rec), snps_shared = nrow(m),
      r_log10p = cor(-log10(m$P_ours), -log10(m$P_rec)),
      spearman_p = cor(m$P_ours, m$P_rec, method = "spearman"),
      gif_ours = gif(m$P_ours), gif_received = gif(m$P_rec),
      min_p_ours = min(m$P_ours), min_p_received = min(m$P_rec),
      top_snp_ours = top_ours[1], top_snp_received = top_rec[1],
      top100_overlap = length(intersect(top_ours, top_rec)))
    pts[[paste(t, model)]] <- m %>% mutate(trait = t, model = model)
  }
}
res <- bind_rows(rows)
write.table(res, "results/compare_with_received.tsv", sep = "\t", quote = F, row.names = F)
print(as.data.frame(res %>% select(trait, model, snps_shared, r_log10p, spearman_p, gif_ours, gif_received, top100_overlap)))

pl <- ggplot(bind_rows(pts), aes(-log10(P_rec), -log10(P_ours))) +
  geom_point(size = 0.4, alpha = 0.3) + geom_abline(linetype = "dotted") +
  facet_grid(trait ~ model, scales = "free") +
  labs(x = "-log10 P, received GWAS", y = "-log10 P, FaST-LMM LOCO (ours)") + theme_minimal()
ggsave("figures/compare_with_received.jpeg", pl, width = 7, height = 10, bg = "white")

### END
