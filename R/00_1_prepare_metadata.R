library(zoo)
library(lme4)

dtmeta1 <- fread("data/mgx/meta/1.metadata1.tsv")
dtmeta1[, TimePoint := "F1"]
dtmeta2 <- fread("data/mgx/meta/1.metadata2.tsv")
dtmeta2[, TimePoint := "F2"]
colnames(dtmeta2) <- gsub("\\.x$|_F2$|_f2$","",colnames(dtmeta2))
colnames(dtmeta2) <- gsub("climbstairs","stairclimb",colnames(dtmeta2))

dt_meta <- rbind(dtmeta1, dtmeta2, fill = T)
setnames(dt_meta, "SampleID","sample")
dt_meta[, subject := substr(new_id, 10, 17)]
dt_meta[, new_id := NULL]

dt_meta[, birthday := as.Date(format(as.Date(gebdat, format = "%d/%m/%y"), "19%y-%m-%d"))]
dt_meta[, birthday := min(birthday, na.rm = T), by = "subject"]
dt_meta[, BD.scaled := scale(birthday)]
dt_meta[, sex := min(sex, na.rm = T), by = "subject"]
dt_meta[, Depth.scaled := scale(Depth)]
dt_meta[, Depth.log10 := log10(Depth)]
dt_meta[, BMI.scaled := scale(BMI)]

rm(dtmeta1)
rm(dtmeta2)