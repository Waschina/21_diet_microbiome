library(data.table)
library(sybil)
library(MicrobiomeGS2) # https://www.nutrinf.uni-kiel.de/MicrobiomeGS2_0.1.0.tar.gz
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

# get compound names
gs_cpds <- fread("data/mgx/gapseq_files/seed_metabolites_edited.tsv")
gs_cpds <- gs_cpds[,.(id, name)]
setnames(gs_cpds, "id","compound")
gs_cpds[compound == "cpd90003", name := "Starch"]
gs_cpds[compound == "cpd00159", name := "Lactate"]

#-------------------------------#
# Preprocess MGX data from Alba #
#-------------------------------#

MxS <- fread("data/mgx/mag_descriptions_update/MAGs.in.all.samples.prevalence.f1.f22.tsv")
setnames(MxS, c("SampleName", "Sample"), c("subject","sample"))
MxS[, sample := substr(sample, 1, 11)]
MxG <- fread("data/mgx/mag_descriptions_update/MAGs.unique.prevalence.f1.f22.tsv")
setnames(MxG, c("SampleName", "Sample"), c("subject","sample"))
MxG[, sample := substr(sample, 1, 11)]
MxS <- merge(MxS, MxG[,.(Genome, repr_MAG = Bins)], by = "Genome")

bifimods <- MxG[grepl("Bifidobacterium", Genus), Bins]

relCountTab <- dcast(MxS, repr_MAG ~ sample, value.var = "RelativeCounts",
                     fill = 0)
relCountsDT <- melt(relCountTab, id.vars = "repr_MAG", variable.name = "sample",
                    value.name = "relCount")
relCountsDT[, sumRelCount := sum(relCount), by = "sample"]
relCountsDT[, relCount.norm := relCount / sumRelCount]

#------------------------------#
# Load models and analyse them #
#------------------------------#

mods <- dir("/mnt/nuuk/2021/Diet_microbiome/models/", recursive = T)
mods <- mods[grepl("\\.RDS$", mods) & !grepl("draft\\.RDS$|rxnXgenes\\.RDS|rxnWeights\\.RDS",mods)]

mods <- lapply(mods, FUN = function(x) {
  mod <- readRDS(paste0("/mnt/nuuk/2021/Diet_microbiome/models/",x))
  
  if(mod@mod_id %in% bifimods) {
    #print("BIFI!")
    ec2_7_1_11 <- mod@react_id[which(mod@react_attr$ec == "2.7.1.11")]
    mod <- rmReact(mod, react = ec2_7_1_11)
    
    mod <- changeBounds(mod, react = "EX_cpd00071_e0", ub = 0) # no acetaldehyde prod
    mod <- changeBounds(mod, react = "EX_cpd00363_e0", ub = 0) # no ethanol prod
    mod <- changeBounds(mod, react = "EX_cpd00159_e0", lb = 0) # no external lactate supply
    mod <- changeBounds(mod, react = "EX_cpd00221_e0", lb = 0) # no external lactate supply
    mod <- changeBounds(mod, react = "EX_cpd00029_e0", lb = 0) # no external acetate supply
    #print(optimizeProb(mod)@lp_obj)
    
  }
  return(mod)
})
mods.ids <- unlist(lapply(mods, function(x) x@mod_id))
names(mods) <- mods.ids

# growth <- lapply(mods, function(x) get_growth(x))
# growth <- data.table(model = names(growth),
#                      gr = unlist(growth))
# 
# # PRODUCED METABOLITES
# prod <- lapply(mods, function(x) get_produced_metabolites(x))
# prod <- rbindlist(prod, idcol = "model")
# prod <- merge(prod, growth, by = "model")
# prod[, cflux := mtf.flux / gr]
# prod <- prod[!grepl("cpd11416", rxn.name)]
# prod[cflux<0, cflux := 0]
# 
# # UTILIZED METABOLITES
# util <- lapply(mods, function(x) get_utilized_metabolites(x))
# util <- rbindlist(util, idcol = "model")
# util <- merge(util, growth, by = "model")
# util[, cflux := mtf.flux / gr]
# util[cflux>0, cflux := 0]
# 
# # FLUXES
# flux <- lapply(mods, function(x) get_flux_distribution(x, exclude.unused = F))
# flux <- rbindlist(flux, idcol = "model")
# flux <- merge(flux, growth, by = "model")
# flux[, cflux := flux / gr]
# 
# tmpwys <- unique(flux$pathway)
# tmpwys <- gsub("\\|","",tmpwys)
# keyec <- keyec[id %in% tmpwys]# filters out non-bacterial pathways
# 
# # Concat Predicted-Phenotype data
# MAG_phenoDT <- copy(MxG[,.(Genome, repr_MAG = Bins)])
# # prod: but
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      prod[ex == "EX_cpd00211_e0", .(repr_MAG = model, prod_but = cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(prod_but), prod_but := 0]
# # prod: ac
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      prod[ex == "EX_cpd00029_e0", .(repr_MAG = model, prod_ac = cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(prod_ac), prod_ac := 0]
# # prod: lac
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      prod[ex %in% c("EX_cpd00221_e0","EX_cpd00159_e0"), .(prod_lac = sum(cflux)), by = "model"][,.(repr_MAG = model, prod_lac)],
#                      all.x = TRUE, by = "repr_MAG")
# MAG_phenoDT[is.na(prod_lac), prod_lac := 0]
# # prod: succ
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      prod[ex == "EX_cpd00036_e0", .(repr_MAG = model, prod_succ = cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(prod_succ), prod_succ := 0]
# # prod: prop
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      prod[ex == "EX_cpd00141_e0", .(repr_MAG = model, prod_ppa = cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(prod_ppa), prod_ppa := 0]
# # prod: for
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      prod[ex == "EX_cpd00047_e0", .(repr_MAG = model, prod_for = cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(prod_for), prod_for := 0]
# # prod: h2
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      prod[ex == "EX_cpd11640_e0", .(repr_MAG = model, prod_h2 = cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(prod_h2), prod_h2 := 0]
# # prod: h2s
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      prod[ex == "EX_cpd00239_e0", .(repr_MAG = model, prod_h2s = cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(prod_h2s), prod_h2s := 0]
# # prod: eth
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      prod[ex == "EX_cpd00363_e0", .(repr_MAG = model, prod_eth = cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(prod_eth), prod_eth := 0]
# # prod: DL-ala
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      prod[ex %in% c("EX_cpd00035_e0","EX_cpd00117_e0"), .(prod_ala = sum(cflux)), by = "model"][,.(repr_MAG = model, prod_ala)],
#                      all.x = TRUE, by = "repr_MAG")
# MAG_phenoDT[is.na(prod_ala), prod_ala := 0]
# # prod: indole
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      prod[ex == "EX_cpd00359_e0", .(repr_MAG = model, prod_indole = cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(prod_indole), prod_indole := 0]
# 
# # util: ac
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex == "EX_cpd00029_e0", .(repr_MAG = model, util_ac = -cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(util_ac), util_ac := 0]
# # util: lac
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex %in% c("EX_cpd00221_e0","EX_cpd00159_e0"), .(util_lac = -sum(cflux)), by = "model"][,.(repr_MAG = model, util_lac)],
#                      all.x = TRUE, by = "repr_MAG")
# MAG_phenoDT[is.na(util_lac), util_lac := 0]
# # util: succ
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex == "EX_cpd00036_e0", .(repr_MAG = model, util_succ = -cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(util_succ), util_succ := 0]
# # util: DL-ala
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex %in% c("EX_cpd00035_e0","EX_cpd00117_e0"), .(util_ala = -sum(cflux)), by = "model"][,.(repr_MAG = model, util_ala)],
#                      all.x = TRUE, by = "repr_MAG")
# MAG_phenoDT[is.na(util_ala), util_ala := 0]
# 
# # Carbs
# # util: glc
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex == "EX_cpd00027_e0", .(repr_MAG = model, util_glc = -cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(util_glc), util_glc := 0]
# # util: frc
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex == "EX_cpd00082_e0", .(repr_MAG = model, util_fru = -cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(util_fru), util_fru := 0]
# # util: lcts
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex == "EX_cpd00208_e0", .(repr_MAG = model, util_lcts = -cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(util_lcts), util_lcts := 0]
# # util: starch
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex == "EX_cpd90003_e0", .(repr_MAG = model, util_starch = -cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(util_starch), util_starch := 0]
# # util: malt
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex == "EX_cpd00179_e0", .(repr_MAG = model, util_malt = -cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(util_malt), util_malt := 0]
# # util: NeuNAC
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex == "EX_cpd00232_e0", .(repr_MAG = model, util_neunac = -cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(util_neunac), util_neunac := 0]
# # util: glcNAC
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex == "EX_cpd00122_e0", .(repr_MAG = model, util_glcnac = -cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(util_glcnac), util_glcnac := 0]
# # util: inulin
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex == "EX_cpd28763_e0", .(repr_MAG = model, util_inulin = -cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(util_inulin), util_inulin := 0]
# # util: sucr
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex == "EX_cpd00076_e0", .(repr_MAG = model, util_sucr = -cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(util_sucr), util_sucr := 0]
# # util: trp
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex == "EX_cpd00065_e0", .(repr_MAG = model, util_trp = -cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(util_trp), util_trp := 0]
# # util: his
# MAG_phenoDT <- merge(MAG_phenoDT,
#                      util[ex == "EX_cpd00119_e0", .(repr_MAG = model, util_his = -cflux)],all.x = TRUE)
# MAG_phenoDT[is.na(util_his), util_his := 0]
# 
# # fluxes through key reactions
# MAG_keyECFlux <- copy(MxG[,.(Genome, repr_MAG = Bins)])
# 
# tmp <- merge(flux, keyec, by.x = "rxn_2",by.y = "keyRea")
# MAG_keyECFlux <- merge(MAG_keyECFlux,
#                        tmp[!duplicated(paste(rxn_2, model)),
#                            .(keyRea = rxn_2, repr_MAG = model, flux = abs(cflux))])
# MAG_keyECFlux <- dcast(MAG_keyECFlux, formula = repr_MAG + Genome ~ keyRea,
#                        value.var = "flux",fill = -1)
# keyecflux_sd <- apply(MAG_keyECFlux[,-(1:2)],2,sd)
# keyec_rm <- names(which(apply(MAG_keyECFlux[,-(1:2)],2,sd) < quantile(keyecflux_sd, prob = 0.25)))
# MAG_keyECFlux <- MAG_keyECFlux[,-keyec_rm, with = F]
