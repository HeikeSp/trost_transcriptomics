
merge_mapman_pgsc <- read.table("data/merge_mapman_pgsc.txt", header=T, sep = "\t")
pgsc <- read.table("data/PGSC_DM_v3.4_g2t2c2p2func_edit.txt", header=F, sep = "\t")

colnames(pgsc) <- c("PGSC_DMG", "PGSC_DMT", "PGSC_DMC", "PGSC_DMP", "function")
head(pgsc)


##############################

func_supp_table_overlap <- function(supp_table, external_up, external_down, col_supp = "pgsc_dmg", col_external = "gene.ID"){
  
  common_up <- which(supp_table[, col_supp] %in% intersect(external_up[, col_external], supp_table[, col_supp]))
  common_down <- which(supp_table[, col_supp] %in% intersect(external_down[, col_external], supp_table[, col_supp]))
  
  common <- rep(NA, dim(supp_table)[1])
  common[common_up] <- "up"
  common[common_down] <- "down"
  
  return(common)
}

###############################

supp_table_s12 <- read.table("data/Supp_Table_S12.txt", header=T, sep = "\t")
dim(supp_table_s12)
# 248 13
length(unique(supp_table_s12$pgsc_dmg)) # 248

# merge s12 with mapman annotation
supp_table_s12_merged <- merge(supp_table_s12, merge_mapman_pgsc, by.x = "pgsc_dmg", by.y = "pgsc_dmg")
dim(supp_table_s12_merged)
# 371 21

write.table(supp_table_s12_merged, "data/Supp_Table_S12_merged.txt", row.names = F, sep = "\t")

###############################

supp_table_s12_new <- read.table("data/Supp_Table_S12_new.txt", header=T, sep = "\t")

# Bengtsson 2014 with Supp Table S12

bengtsson2014a_common_S12 <- func_supp_table_overlap(supp_table = supp_table_s12_new, 
                                                     external_up = bengtsson2014a_S3_up_dmg,
                                                     external_down = bengtsson2014a_S3_down_dmg,
                                                     col_external = "PGSC_DMG")

############ 
# create final table with additional information

supp_table_s12_final <- data.frame(
  pgsc_dmg = supp_table_s12_new$pgsc_dmg, 
  bengtsson2014 = bengtsson2014a_common_S12)

write.table(supp_table_s12_final, "output/Supp_Table_S12_final.txt", row.names = F, sep = "\t")


##############################
# Load Supp Table 10

supp_table_s10 <- read.table("data/Supp_Table_S10_extended.txt", header=T, sep = "\t")
dim(supp_table_s10)
# 1908 13
dim(supp_table_s10)[1]


############
# Zhang 2014

zhang2014_common <- func_supp_table_overlap(supp_table = supp_table_s10, 
                                            external_up = zhang2014_up,
                                            external_down = zhang2014_down,
                                            col_external = "Gene_id")

############
# Bengtsson 2014a

bengtsson2014a_common <- func_supp_table_overlap(supp_table = supp_table_s10, 
                                                 external_up = bengtsson2014a_S3_up_dmg,
                                                 external_down = bengtsson2014a_S3_down_dmg,
                                                 col_external = "PGSC_DMG")

############
# Wiesel 2015

wiesel2015_S1_ABA_common <- func_supp_table_overlap(supp_table = supp_table_s10, 
                                                    external_up = wiesel2015_S1_ABA_up,
                                                    external_down = wiesel2015_S1_ABA_down)

wiesel2015_S2_ABA_common <- func_supp_table_overlap(supp_table = supp_table_s10, 
                                                    external_up = wiesel2015_S2_ABA_up,
                                                    external_down = wiesel2015_S2_ABA_down)

wiesel2015_S1_MeJA_common <- func_supp_table_overlap(supp_table = supp_table_s10, 
                                                     external_up = wiesel2015_S1_MeJA_up,
                                                     external_down = wiesel2015_S1_MeJA_down)

wiesel2015_S2_MeJA_common <- func_supp_table_overlap(supp_table = supp_table_s10, 
                                                     external_up = wiesel2015_S2_MeJA_up,
                                                     external_down = wiesel2015_S2_MeJA_down)

wiesel2015_S1_SA_common <- func_supp_table_overlap(supp_table = supp_table_s10, 
                                                     external_up = wiesel2015_S1_SA_up,
                                                     external_down = wiesel2015_S1_SA_down)

wiesel2015_S2_SA_common <- func_supp_table_overlap(supp_table = supp_table_s10, 
                                                   external_up = wiesel2015_S2_SA_up,
                                                   external_down = wiesel2015_S2_SA_down)

############ 
# create final table with additional information

supp_table_s10_final <- data.frame(
  pgsc_dmg = supp_table_s10$pgsc_dmg, 
  zhang2014 = zhang2014_common,
  bengtsson2014 = bengtsson2014a_common,
  wiesel2015_ABA_1h = wiesel2015_S1_ABA_common,
  wiesel2015_ABA_6h = wiesel2015_S2_ABA_common,
  wiesel2015_MeJA_1h = wiesel2015_S1_MeJA_common,
  wiesel2015_MeJA_6h = wiesel2015_S2_MeJA_common,
  wiesel2015_SA_1h = wiesel2015_S1_SA_common,
  wiesel2015_SA_6h = wiesel2015_S2_SA_common)

write.table(supp_table_s10_final, "output/Supp_Table_S10_final.txt", row.names = F, sep = "\t")
# write.table(supp_table_s10_merged, "output/Supp_Table_S10_merged.txt", row.names = F, sep = "\t")



