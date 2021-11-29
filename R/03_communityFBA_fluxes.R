library(ggplot2)
library(ggtext)

# Prepare foodgrp data
dt_dietGRP <- copy(DietaryInfo_mainGrp)
setkey(dt_dietGRP, "sample")

# prepare diet component data
dt_dietComp <- copy(DietaryInfo_components)
totalenergy <- dt_dietComp[item_code == "GJ",.(sample, kJ = item_value)]
dt_dietComp <- merge(dt_dietComp, totalenergy)
dt_dietComp <- dt_dietComp[item_code != "GJ"]
setkey(dt_dietComp, "sample")

dt_dietGRP <- merge(dt_dietGRP, totalenergy)

# coupling parameters
parc <- 40
paru <- 0.01

relspls <- unique(relCountsDT$sample)
cFBA_out <- list()
k <- 1
for(ispl in relspls) {
  cat("\n",k, "/", length(relspls),"\n")
  dt_tmp <- relCountsDT[sample == ispl & relCount.norm > 0]
  
  commFBA <- communityFBA_FB(models = mods[dt_tmp$repr_MAG],
                             model.prop = dt_tmp$relCount.norm,
                             cpl_c = parc, cpl_u = paru)
  commFBA$met.interchange[, gr := commFBA$community.growth]
  cFBA_out[[ispl]] <- commFBA
  k <- k + 1
}

names(cFBA_out) <- relspls

outflow <- lapply(cFBA_out, function(x) x$met.interchange)
outflow <- rbindlist(outflow, idcol = "sample")
outflow[, o.flux.n := o.flux / gr]
#outflow[, o.flux.n := o.flux]
outflow <- dcast(outflow, sample ~ rxn, value.var = "o.flux.n", fill = 0)
outflow <- melt(outflow, id.vars = "sample", variable.name = "compound", value.name = "flux")
outflow <- outflow[compound != "EX_cpd00001_e0"]
outflow[compound == "EX_cpd00221_e0", compound := "EX_cpd00159_e0"] # merge both lactate forms
outflow <- outflow[, .(flux = sum(flux)), by = c("sample","compound")]
outflow[abs(flux)< 1e-6 & flux != 0, flux := 0]
# Filter compounds:
# * remove compound were SD of flux is in the first quartile (0-25%)
# * whose SD is lower than 1e-5
# * contains a zero flux in more than 50% of the samples -> rm sample
tmpsd <- outflow[, .(sd = sd(flux), n0 = sum(abs(flux)<1e-5)), by = compound]
tmpsd <- as.character(tmpsd[sd > quantile(sd, prob = 0.25) & sd > 1e-5 & n0 <= .N/2, compound])
outflow <- outflow[compound %in% tmpsd]
setkey(outflow, "sample")
relcpds <- unique(as.character(outflow$compound))

#-----------------------------------------------------------#
# spearman correlation of outflow fluxes with food groups   #
#-----------------------------------------------------------#
outflow_dietGrp_stat <- list()
relspls <- intersect(as.character(outflow$sample),
                     unique(dt_dietGRP$sample))
k <- 1
for(idiet in unique(dt_dietGRP$Group)) {
  x <- dt_dietGRP[relspls][Group == idiet, gramm_mainGrp]
  for(ipheno in relcpds) {
    y <- outflow[compound == ipheno][relspls][,flux]
    
    corStat <- cor.test(x,y, method = "spearman")
    
    Q2 <- quantile(y, prob = 0.5)
    if(Q2 == 0)
      Q2 <- mean(y)
    
    outflow_dietGrp_stat[[k]] <- data.table(item    = idiet,
                                            pheno   = ipheno,
                                            pval    = corStat$p.value,
                                            est     = corStat$estimate,
                                            Q2      = Q2)
    
    k <- k + 1
  }
}
outflow_dietGrp_stat <- rbindlist(outflow_dietGrp_stat)
outflow_dietGrp_stat[, padj := p.adjust(pval, method = "fdr"), by = "item"]
#outflow_dietGrp_stat[padj < 0.05]
outflow_dietGrp_stat[, category := "<b>Food groups</b>"]


#-----------------------------------------------------------#
# spearman correlation of outflow fluxes with food compon.  #
#-----------------------------------------------------------#
outflow_dietComp_stat <- list()
relspls <- intersect(as.character(outflow$sample),
                     unique(dt_dietComp$sample))
k <- 1
for(idiet in unique(dt_dietComp$item_code)) {
  x <- dt_dietComp[relspls][item_code == idiet, item_value / kJ]
  for(ipheno in relcpds) {
    y <- outflow[compound == ipheno][relspls][,flux]
    
    corStat <- cor.test(x,y, method = "spearman")
    
    Q2 <- quantile(y, prob = 0.5)
    if(Q2 == 0)
      Q2 <- mean(y)
    
    outflow_dietComp_stat[[k]] <- data.table(item    = idiet,
                                            pheno   = ipheno,
                                            pval    = corStat$p.value,
                                            est     = corStat$estimate,
                                            Q2      = Q2)
    
    k <- k + 1
  }
}
outflow_dietComp_stat <- rbindlist(outflow_dietComp_stat)
outflow_dietComp_stat[, padj := p.adjust(pval, method = "fdr"), by = "item"]
#outflow_dietComp_stat[padj < 0.05]
outflow_dietComp_stat[, category := "<b>Nutrients</b> <i>(%E)</i>"]



outflow_statsum <- rbind(outflow_dietGrp_stat,outflow_dietComp_stat)
rel_items <- outflow_statsum[padj < 0.05, unique(item)]
rel_cpds  <- outflow_statsum[padj < 0.05, unique(pheno)]
outflow_statsum[padj < 0.05, p.plot := "P[\"adj\"] < 0.05"]
# tmp <- copy(keyec[,.(pheno = keyRea, keyRea.name, pwy = name)])
# tmp[, keyRea.name := keyRea.name[1], by = pheno]
# tmp <- tmp[, .(pwy = paste(pwy, collapse = ";")), by = .(pheno, keyRea.name)]
# #tmp <- tmp[!duplicated(pheno)]
# tmp[, pheno.plot := paste0("[<b>",pwy,"</b>] ",keyRea.name)]
# 
# outflow_statsum <- merge(outflow_statsum, tmp,
#                        by = "pheno")
# outflow_statsum[, pheno.plot := gsub("\\(S\\)-","",pheno.plot)]
# outflow_statsum[, pheno.plot := gsub("\\(R\\)-","",pheno.plot)]

outflow_statsum_plot <- outflow_statsum[item %in% rel_items & pheno %in% rel_cpds]
outflow_statsum_plot[, pheno := gsub("^EX_|_e0$","",pheno)]
outflow_statsum_plot <- merge(outflow_statsum_plot, gs_cpds, by.x = "pheno",
                              by.y = "compound")

outflow_statsum_plot[, direction := "Consumption"]
outflow_statsum_plot[Q2 > 0, direction := "Production"]
outflow_statsum_plot[direction == "Consumption", est := -est]

outflow_statsum_plot <- merge(outflow_statsum_plot, nutr_names, by.x = "item", by.y = "item_code", all.x = T)
outflow_statsum_plot[is.na(item_name), item_name := item]

cluster_factors <- function(dt, by.x, by.y, value) {
  # dt <- keyec_statsum_plot[category == "<b>Food groups</b>"]
  # by.x <- "item"
  # by.y <- "pheno.plot"
  # value <- "est"
  f <- as.formula(paste(by.x, "~ ",by.y))
  tmp_wide <- dcast(dt, formula = f, value.var = value)
  
  rn <- tmp_wide[[by.x]]
  tmp_wide <- as.matrix(tmp_wide[,-1])
  
  rclust <- hclust(dist(tmp_wide))
  return(rn[rclust$order])
}

fgrp_order <- cluster_factors(outflow_statsum_plot[category == "<b>Food groups</b>"],
                              by.x = "item_name", by.y = "name",
                              value = "est")
fcomp_order <- cluster_factors(outflow_statsum_plot[category == "<b>Nutrients</b> <i>(%E)</i>"],
                               by.x = "item_name", by.y = "name",
                               value = "est")

outflow_statsum_plot$item_name <- factor(outflow_statsum_plot$item_name,
                                  levels = c(fgrp_order, fcomp_order))

cons_pwy_order <- cluster_factors(outflow_statsum_plot[direction == "Consumption"], 
                             by.x = "name",
                             by.y = "item_name", value = "est")
prod_pwy_order <- cluster_factors(outflow_statsum_plot[direction == "Production"], 
                                  by.x = "name",
                                  by.y = "item_name", value = "est")

outflow_statsum_plot$name <- factor(outflow_statsum_plot$name,
                                        levels = c(cons_pwy_order, prod_pwy_order))

p_commFBA_scor <- ggplot(outflow_statsum_plot,
                         aes(item_name, name, fill = est)) +
  geom_tile() +
  geom_point(aes(shape = p.plot)) +
  scale_shape_manual(values = c(19), na.translate = F, labels = expression(P['FDR'] < 0.05)) +
  scale_fill_gradient2(high = "#ca0020", mid = "#f7f7f7", low = "#0571b0") +
  scale_y_discrete(position = "right", expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  labs(x = "Dietary intake", y = "Compound", shape = "", fill = "Spearman's rho") +
  facet_grid(direction~category, space = "free", scales = "free", switch = "y") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_markdown(),
        strip.text.y = element_text(face = "bold", color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y.right = element_markdown(color = "black"),
        legend.position = "right") +
  guides(fill = guide_colourbar(title.position = "top"))

ggsave("output/plots/commFBA_scor.pdf", plot = p_commFBA_scor,
       width = 6.55, height = 6)
fwrite(outflow_statsum_plot, "output/tables/commFBA_spearman.csv")

#--------------------------#
# Calculate change in diet #
# food groups              #
#--------------------------#
dt_dietGRP[TimePoint == "F1", tmp := gramm_mainGrp]
dt_dietGRP[, tmp := min(tmp, na.rm = T), by = .(subject, Group)]
dt_dietGRP[, gramm_mainGrp_chg := gramm_mainGrp - tmp]
#dt_dietGRP[, tmp := NULL]

dt_dgrps <- merge(dt_dietGRP, dt_meta, by = c("subject","TimePoint"))


foodGrpChange_stat <- list()
k <- 1
for(idiet in unique(dt_dietGRP$Group)) {
  zut <- copy(dt_dgrps[Group == idiet])
  setkey(zut, "sample.x")
  relspls <- intersect(zut$sample.x,
                       as.character(outflow$sample))
  zut <- zut[relspls]
  for(ipheno in relcpds) {
    zut$rxnfreq <- outflow[compound == ipheno][zut$sample.x][,flux]
    
    if(sd(zut$rxnfreq)> 0) {
      lme0 <- lmer(rxnfreq ~ TimePoint + tmp + gramm_mainGrp + BMI.scaled + BD.scaled + sex + Depth.log10 + (1|subject),
                   data = zut, REML = F)
      lme1 <- lmer(rxnfreq ~ TimePoint + tmp +                 BMI.scaled + BD.scaled + sex + Depth.log10 + (1|subject),
                   data = zut, REML = F)
      lmeaov <- anova(lme0, lme1)
      
      foodGrpChange_stat[[k]] <- data.table(item      = idiet,
                                            pheno     = ipheno,
                                            pval      = lmeaov$`Pr(>Chisq)`[2],
                                            tval_item = summary(lme0)$coefficients["gramm_mainGrp","t value"],
                                            Q1 = quantile(zut$rxnfreq, prob = 0.25),
                                            Q2 = quantile(zut$rxnfreq, prob = 0.5),
                                            Q3 = quantile(zut$rxnfreq, prob = 0.75))
      
      k <- k + 1
    }
  }
}
foodGrpChange_stat <- rbindlist(foodGrpChange_stat)
foodGrpChange_stat[, padj := p.adjust(pval, method = "fdr"), by = "item"]
foodGrpChange_stat[padj < 0.05]
foodGrpChange_stat[, category := "<b>Food groups</b>"]

#--------------------------#
# Calculate change in diet #
# food components          #
#--------------------------#
dt_dietComp[TimePoint == "F1", tmp := item_value / kJ]
dt_dietComp[, tmp := min(tmp, na.rm = T), by = .(subject, item_code)]
dt_dietComp[, gramm_itemE := item_value/kJ]
#dt_dietGRP[, tmp := NULL]

dt_dcomp <- merge(dt_dietComp, dt_meta, by = c("subject","TimePoint"))


foodCompChange_stat <- list()
k <- 1
jn <- length(unique(dt_dcomp$item_code))
j <- 1
for(idiet in unique(dt_dcomp$item_code)) {
  cat("\r",j,"/",jn)
  zut <- copy(dt_dcomp[item_code == idiet])
  setkey(zut, "sample.x")
  relspls <- intersect(zut$sample.x,
                       as.character(outflow$sample))
  zut <- zut[relspls]
  for(ipheno in relcpds) {
    zut$rxnfreq <- outflow[compound == ipheno][zut$sample.x][, flux]
    
    if(sd(zut$rxnfreq)> 0) {
      lme0 <- lmer(rxnfreq ~ TimePoint + tmp + gramm_itemE + BMI.scaled + BD.scaled + sex + Depth.log10 + (1|subject),
                   data = zut, REML = F)
      lme1 <- lmer(rxnfreq ~ TimePoint + tmp +               BMI.scaled + BD.scaled + sex + Depth.log10 + (1|subject),
                   data = zut, REML = F)
      lmeaov <- anova(lme0, lme1)
      
      foodCompChange_stat[[k]] <- data.table(item      = idiet,
                                             pheno     = ipheno,
                                             pval      = lmeaov$`Pr(>Chisq)`[2],
                                             tval_item = summary(lme0)$coefficients["gramm_itemE","t value"],
                                             Q1 = quantile(zut$rxnfreq, prob = 0.25),
                                             Q2 = quantile(zut$rxnfreq, prob = 0.5),
                                             Q3 = quantile(zut$rxnfreq, prob = 0.75))
      
      k <- k + 1
    }
  }
  j <- j + 1
}
foodCompChange_stat <- rbindlist(foodCompChange_stat)
foodCompChange_stat[, padj := p.adjust(pval, method = "fdr"), by = "item"]
foodCompChange_stat[padj < 0.05]
foodCompChange_stat[, category := "<b>Nutrients</b> <i>(%E)</i>"]



outflowchg_statsum <- rbind(foodGrpChange_stat,foodCompChange_stat)
rel_items <- outflowchg_statsum[padj < 0.05, unique(item)]
rel_cpds  <- outflowchg_statsum[padj < 0.05, unique(pheno)]
outflowchg_statsum[padj < 0.05, p.plot := "P[\"adj\"] < 0.05"]

outflowchg_statsum_plot <- outflowchg_statsum[item %in% rel_items & pheno %in% rel_cpds]
outflowchg_statsum_plot[, pheno := gsub("^EX_|_e0$","",pheno)]
outflowchg_statsum_plot <- merge(outflowchg_statsum_plot, gs_cpds, by.x = "pheno",
                              by.y = "compound")

outflowchg_statsum_plot[, direction := "Consumption"]
outflowchg_statsum_plot[Q2 > 0, direction := "Production"]
outflowchg_statsum_plot[direction == "Consumption", tval_item := -tval_item]

outflowchg_statsum_plot <- merge(outflowchg_statsum_plot, nutr_names, by.x = "item", by.y = "item_code", all.x = T)
outflowchg_statsum_plot[is.na(item_name), item_name := item]

fgrp_order <- cluster_factors(outflowchg_statsum_plot[category == "<b>Food groups</b>"],
                              by.x = "item_name", by.y = "name",
                              value = "tval_item")
fcomp_order <- cluster_factors(outflowchg_statsum_plot[category == "<b>Nutrients</b> <i>(%E)</i>"],
                               by.x = "item_name", by.y = "name",
                               value = "tval_item")

outflowchg_statsum_plot$item_name <- factor(outflowchg_statsum_plot$item_name,
                                    levels = c(fgrp_order, fcomp_order))

cons_pwy_order <- cluster_factors(outflowchg_statsum_plot[direction == "Consumption"], 
                                  by.x = "name",
                                  by.y = "item_name", value = "tval_item")
prod_pwy_order <- cluster_factors(outflowchg_statsum_plot[direction == "Production"], 
                                  by.x = "name",
                                  by.y = "item_name", value = "tval_item")

outflowchg_statsum_plot$name <- factor(outflowchg_statsum_plot$name,
                                    levels = c(cons_pwy_order, prod_pwy_order))

p_commfba_chg <- ggplot(outflowchg_statsum_plot,
                        aes(item_name, name, fill = tval_item)) +
  geom_tile() +
  geom_point(aes(shape = p.plot)) +
  scale_shape_manual(values = c(19), na.translate = F, labels = expression(P['FDR'] < 0.05)) +
  scale_fill_gradient2(high = "#ca0020", mid = "#f7f7f7", low = "#0571b0") +
  scale_y_discrete(position = "right", expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  labs(x = "Dietary intake", y = "Compound", shape = "", fill = "t value") +
  facet_grid(direction~category, space = "free", scales = "free", switch = "y") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_markdown(),
        strip.text.y = element_text(face = "bold", color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y.right = element_markdown(color = "black"),
        legend.position = "right") +
  guides(fill = guide_colourbar(title.position = "top"))

ggsave("output/plots/commFBA_change.pdf", plot = p_commfba_chg,
       width = 8.25, height = 7.7)
fwrite(outflowchg_statsum_plot, "output/tables/commFBA_LME.csv")

