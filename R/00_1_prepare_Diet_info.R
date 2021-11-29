library(data.table)
library(stringr)

DietaryInfo <- fread("data/mgx/Dietary_data_food_groups_nutrients/11.food.groups.nutrients.f1.f2.tsv")
DietaryInfo[, sample := paste0(Samples,"_",TimePoint)]
DietaryInfo[, subject := paste0(Samples)]
DietaryInfo[, `:=`(Samples = NULL, new_id = NULL)]

# Total diet gramm per day
DietaryInfo_DietMass <- DietaryInfo[, .(daily_gramm = sum(gramm)), by = sample]

# Individual dietary components
comp_idxrange <- which(colnames(DietaryInfo) == "EALA") :
                   which(colnames(DietaryInfo) == "ZW")
DietaryInfo_components <- melt(DietaryInfo, id.vars = c("sample",
                                                        "subject","TimePoint"),
                               measure.vars = colnames(DietaryInfo)[comp_idxrange],
                               variable.name = "item_code", value.name = "item_value")
DietaryInfo_components <- DietaryInfo_components[, .(item_value = sum(item_value)), by = c("sample","subject","TimePoint",
                                                 "item_code")]
DietaryInfo_components <- merge(DietaryInfo_components, DietaryInfo_DietMass,
                                by = "sample")

# summed main groups
DietaryInfo_mainGrp <- DietaryInfo[, .(gramm_mainGrp = sum(gramm)),
                                   by = c("sample","subject","TimePoint","Group")]
DietaryInfo_mainGrp <- merge(DietaryInfo_mainGrp, DietaryInfo_DietMass, 
                             by = "sample")

# EPIX nutrient code translation
nutr_names <- fread("data/EPIC_nutrients_codes.tsv")
nutr_names[, unit := str_match(item_name, "\\[.*\\]")]
nutr_names[, item_name := gsub(" \\[.*\\]","", item_name)]
