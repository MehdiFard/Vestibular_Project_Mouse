# Preparing edit sites

library(data.table)
library(dplyr)

edsites <- fread("Pantry_RNA_editing/phenotyping/input/ref/TABLE1_mm10_v3.txt")

edsites_cln <- edsites %>%
            select(Region,Position,Ref,Ed,Strand) %>%
            filter(Region != "chrX_GL456233_random") %>%
            mutate(end = Position + 1,
            chr = gsub("chr", "", Region)) %>%
            select(chr,Position,end,Ref,Ed,Strand)

write.table(edsites_cln, "Pantry_RNA_editing/phenotyping/input/ref/edit_sites_clean_mm10_v3.bed",
            col.names = F, row.names = F, quote = F, sep = "\t")

### END