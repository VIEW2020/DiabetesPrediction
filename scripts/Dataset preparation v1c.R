# Diabetes Prediction for Primary Prevention v3.0
# Create a diabetes-free population
# 
# Last updated April 2020 - created by Billy Wu
# 
# Exclude IF
# •	Prior admission for diabetes; OR
# •	Treated with antidiabetic drugs in last 6 months; OR
# •	Noted as diabetic in PREDICT; OR
# • Non-existing HbA1c test in prior 2 years; OR
# •	Any HbA1c in prior 2 years >= 50mmol/mol

# Study Index Date = first available HbA1c test date in TestSafe*

# Further refinement:
# There is a chance that individual is not in TestSafe catchment at index time point.
# 

library(fst)
library(data.table)

# ---- A. Find first Predict in which inclusion is met ----

# All DM admissions
files.path  <- "V:/source_data/R/PREDICT/2018/National Collection/"
files.list <- list.files(files.path, pattern = "PUBLIC")

VIEW_DM_ICD <- readxl::read_xlsx("V:/common_lookups/Definitions/VIEW_CVD_ICD10_II_28JUN18_Clean042.xlsx", sheet = 2)
hx.diabetes <- VIEW_DM_ICD$CLINICALCODE[which(VIEW_DM_ICD$hx_diabetes=="Y")]

for(file in files.list){
 
 DAT <- read.fst(paste0(files.path, file), as.data.table = T)
 DAT <- DAT[CLIN_CD_10 %in% hx.diabetes]
 
 if(file == files.list[1]){
  ADM_DM <- DAT
 } else {
  ADM_DM <- rbind(ADM_DM, DAT)
 }
 
 print(paste0(file, " completed"))
 
}

# All anti-DM Dispensing
files.list <- list.files(files.path, pattern = "PHH")

DIM_FORM_PACK <- read.fst("V:/common_lookups/CURRENT_PHARMS_LOOKUP_VIEW.fst", 
                          as.data.table = T)

all.diabetes <- DIM_FORM_PACK$DIM_FORM_PACK_SUBSIDY_KEY[which(DIM_FORM_PACK$antidiabetes == 1)]

for(file in files.list){
 
 DAT <- read.fst(paste0(files.path, file), as.data.table = T)
 DAT <- DAT[DIM_FORM_PACK_SUBSIDY_KEY %in% all.diabetes]
 
 if(file == files.list[1]){
  MX_DM <- DAT
 } else {
  MX_DM <- rbind(MX_DM, DAT)
 }
 
 print(paste0(file, " completed"))
 
}

# All HbA1c 
files.path  <- "V:/source_data/R/PREDICT/2018/TestSafe/"
files.list <- list.files(files.path, pattern = "HBA1C")

for(file in files.list){
 
 DAT <- read.fst(paste0(files.path, file), as.data.table = T)
 
 if(file == files.list[1]){
  TS_DM <- DAT
 } else {
  TS_DM <- rbind(TS_DM, DAT)
 }
 
 print(paste0(file, " completed"))
 
}

# Save
write.fst(ADM_DM, "V:/rpyl001/DM Risk Prediction/Working/ADM_DM.fst", 75)
write.fst(MX_DM, "V:/rpyl001/DM Risk Prediction/Working/MX_DM.fst", 75)
write.fst(TS_DM, "V:/rpyl001/DM Risk Prediction/Working/TS_DM.fst", 75)

#--- Import Data --
PREDICT <- read.fst("V:/source_data/R/PREDICT/2018/Cleaned_PREDICT_2018_All_Records_v1.fst",
                    as.data.table = T)

ADM_DM <- read.fst("V:/rpyl001/DM Risk Prediction/Working/ADM_DM.fst", 
                   as.data.table = T)

MX_DM <- read.fst("V:/rpyl001/DM Risk Prediction/Working/MX_DM.fst", 
                  as.data.table = T)

TS_DM <- read.fst("V:/rpyl001/DM Risk Prediction/Working/TS_DM.fst", 
                  as.data.table = T)

#--- FUN: FindExclusion ---
FindExclusion <- function(x){
 
 HX_DM  <- copy(ADM_DM)[VSIMPLE_INDEX_MASTER %in% x$VSIMPLE_INDEX_MASTER
 ][, visit_date := x$view_visit_date[match(VSIMPLE_INDEX_MASTER, x$VSIMPLE_INDEX_MASTER)]]
 
 MX6_DM <- copy(MX_DM)[VSIMPLE_INDEX_MASTER %in% x$VSIMPLE_INDEX_MASTER
 ][, visit_date := x$view_visit_date[match(VSIMPLE_INDEX_MASTER, x$VSIMPLE_INDEX_MASTER)]]
 
 TS_DM <- copy(TS_DM)[VSIMPLE_INDEX_MASTER %in% x$VSIMPLE_INDEX_MASTER
 ][, visit_date := x$view_visit_date[match(VSIMPLE_INDEX_MASTER, x$VSIMPLE_INDEX_MASTER)]]
 
 PT_DM <- copy(PREDICT)[VSIMPLE_INDEX_MASTER %in% x$VSIMPLE_INDEX_MASTER
 ][, .(VSIMPLE_INDEX_MASTER, view_visit_date, pt_en_hba1c_mm, pt_en_hba1c_mmdate)]
 
 HX_DM  <- HX_DM[EVSTDATE < visit_date]
 MX6_DM <- MX6_DM[(DATE_DISPENSED - visit_date >= -182) & (DATE_DISPENSED <= visit_date)]
 TS_DM <- TS_DM[RESULT_DATE - visit_date <= 30]
 PT_DM <- PT_DM[pt_en_hba1c_mmdate - view_visit_date <= 30]
 
 setnames(TS_DM, "EN_OBSR_RESULT_NUM", "RESULT_MM")
 setnames(PT_DM, "pt_en_hba1c_mm", "RESULT_MM")
 
 LAB_DM <- rbind(TS_DM[,.(VSIMPLE_INDEX_MASTER, RESULT_MM)], 
                 PT_DM[,.(VSIMPLE_INDEX_MASTER, RESULT_MM)])
 
 A1C_50 <- LAB_DM[LAB_DM[, .I[any(RESULT_MM >= 50)], by = VSIMPLE_INDEX_MASTER]$V1] # Filter by group with one or more tests >= 50mmol/mol
 
 x$exclude <- +(x$VSIMPLE_INDEX_MASTER %in% HX_DM$VSIMPLE_INDEX_MASTER |    # DM Admissions
                 x$VSIMPLE_INDEX_MASTER %in% MX6_DM$VSIMPLE_INDEX_MASTER |  # DM Treatment
                 !x$VSIMPLE_INDEX_MASTER %in% LAB_DM$VSIMPLE_INDEX_MASTER | # Absent HbA1c value  
                 x$VSIMPLE_INDEX_MASTER %in% A1C_50$VSIMPLE_INDEX_MASTER |  # HbA1c value >= 50mmol/mol
                 x$pt_diabetes %in% 1:3)                                    # DM in Predict
 
 return(x$exclude)
 
}

# --- Apply FUN as List Operation ---
# Examine each visit and apply exclusion criteria. Continue with each nth visit to find those who could meet the inclusion criteria in subsequent visits.
# Nb: Starting with 564180 participants, 394265 met criteria after first visit; In the remaining 169915, 9386 meet criteria in subsequent visits.
#     Total 403651 met the criteria

# Limit PREDICT to 2006 - 2018
PREDICT <- PREDICT[year(PREDICT$view_visit_date) >= 2006]

PREDICT_LS <- lapply(split(PREDICT, by = "view_visit_seq"), 
                     function(x){
                      
                      vis.index <- unique(x$view_visit_seq)
                      
                      if(vis.index > 1){
                       x <- x[!VSIMPLE_INDEX_MASTER %in% x_last$VSIMPLE_INDEX_MASTER]
                      }
                      
                      x <- x[, exclude := FindExclusion(x)][exclude == 0]
                      
                      assign("x_last", x, envir = parent.frame(2L))
                      
                      print(paste("Visit", vis.index, "completed."))
                      
                      return(x)
                      
                     })

PREDICT_BL <- do.call("rbind", PREDICT_LS)
PREDICT_BL <- PREDICT_BL[PREDICT_BL[, .I[which.min(view_visit_seq)], 
                                    by = VSIMPLE_INDEX_MASTER]$V1] # 405576

# Apply Age Limits: 25-74
PREDICT_BL <- PREDICT_BL[view_ag_age %in% 25:74] # -16186

# 389390 met the criteria
# write.fst(PREDICT_BL, "V:/rpyl001/DM Risk Prediction/Working/PREDICT_BL_v1c.fst", 75)

# --- Index HbA1c ---

PREDICT_BL <- read.fst("V:/rpyl001/DM Risk Prediction/Working/PREDICT_BL_v1c.fst", as.data.table = T)

# Earliest available HbA1c date: Draws from both PREDICT and TestSafe
TS2_DM <- copy(TS_DM)[VSIMPLE_INDEX_MASTER %in% PREDICT_BL$VSIMPLE_INDEX_MASTER
][, visit_date := PREDICT_BL$view_visit_date[match(VSIMPLE_INDEX_MASTER, PREDICT_BL$VSIMPLE_INDEX_MASTER)]]

PT2_DM <- copy(PREDICT)[VSIMPLE_INDEX_MASTER %in% PREDICT_BL$VSIMPLE_INDEX_MASTER
][, visit_date := PREDICT_BL$view_visit_date[match(VSIMPLE_INDEX_MASTER, PREDICT_BL$VSIMPLE_INDEX_MASTER)]]

setnames(TS2_DM, "EN_OBSR_RESULT_NUM", "RESULT_MM")
setnames(PT2_DM, "pt_en_hba1c_mm", "RESULT_MM")
setnames(PT2_DM, "pt_en_hba1c_mmdate", "RESULT_DATE")

LAB_DM <- rbind(TS2_DM[!is.na(RESULT_MM),.(VSIMPLE_INDEX_MASTER, visit_date, RESULT_MM, RESULT_DATE)], 
                PT2_DM[!is.na(RESULT_MM),.(VSIMPLE_INDEX_MASTER, visit_date, RESULT_MM, RESULT_DATE)])

LAB_DM <- LAB_DM[RESULT_DATE - visit_date <= 30] # -79915 removed due to no result prior to predict visit (nb: they may have entered cohort due to "no evidence of meeting exclusion criteria")

summary(LAB_DM$RESULT_MM) # QC - ensure all results (of cohort members) are <50mmol/mol

# Capture nearest result to baseline
# Evalutes all available HbA1c tests within 30 days of Predict (no limit on earliest record)
library(Allie)
library(dplyr)

LAB_DM[, by = VSIMPLE_INDEX_MASTER
       , flag := FindNearestDate(index = visit_date, 
                                 comparison = RESULT_DATE, 
                                 mode = "n",
                                 from = -6000,
                                 to = 30)]
DM_INDEX <- LAB_DM[flag == 1, 
                   unique(.SD), 
                   .SDcols = c("VSIMPLE_INDEX_MASTER", "visit_date", "RESULT_DATE")]

DM_INDEX$days_diff <- as.numeric(DM_INDEX$RESULT_DATE - DM_INDEX$visit_date)


DM_INDEX$hba1c_lag <- 
 factor(dplyr::case_when(
  DM_INDEX$days_diff %in% 30:0 ~ "+30-0 days",
  DM_INDEX$days_diff %in% -1:-365 ~ "1-12 mths",
  DM_INDEX$days_diff %in% -366:-730 ~ "13-24 mths",
  DM_INDEX$days_diff %in% -731:-913 ~ "25-30 mths",
  DM_INDEX$days_diff %in% -914:-1096 ~ "31-36 mths",
  DM_INDEX$days_diff %in% -1097:-1279 ~ "37-42 mths",
  DM_INDEX$days_diff %in% -1280:-1462 ~ "43-48 mths",
  TRUE ~ ">48 mths"),
  levels = c("+30-0 days", "1-12 mths", "13-24 mths", "25-30 mths", "31-36 mths", "37-42 mths", "43-48 mths", ">48 mths"))

table(DM_INDEX$hba1c_lag)

LAB_DM2$days_diff_label <- LAB_DM2



DM_INDEX  <- copy(LAB_DM)[LAB_DM[, .I[which.min(RESULT_DATE)], 
                                 by = VSIMPLE_INDEX_MASTER]$V1]

setnames(DM_INDEX, "RESULT_DATE", "study_index_date")


# Criteria: must have a result, ALL results ANYTIME prior to Predict visit must be <50mmol 
# Merge to PREDICT Baseline Data
PREDICT_BL <- merge(DM_INDEX[,.(VSIMPLE_INDEX_MASTER, study_index_date, study_index_age)],
                    PREDICT_BL,
                    by = "VSIMPLE_INDEX_MASTER",
                    all.x = F)

# Remove anyone whose index year (earliest HbA1c) was prior to 2006
PREDICT_BL2 <- PREDICT_BL[year(study_index_date) >= 2006] #-10914

write.fst(PREDICT_BL, "V:/rpyl001/DM Risk Prediction/Working/PREDICT_BL2_v1c.fst", 75) #298561 remaining

294014 - 298561
# ---- B.   History / Outcomes ----

PREDICT <- read.fst("V:/rpyl001/DM Risk Prediction/Working/PartA.fst",
                    as.data.table = T)

# Loaded pre-prepared national collection 
files.path <- "V:/source_data/R/PREDICT/2018/National Collection/"
files.list <- list.files(files.path, pattern="PUBLIC")

ALL_PUBLIC <- rbindlist(lapply(files.list, 
                               function(x)
                                read.fst(paste0(files.path, x), 
                                         as.data.table = T)))

ALL_DEATHS <- read.fst(paste0(files.path, "VSIMPLE_DEATHS_v1.fst"), as.data.table = T)


# -- 1.   Definitions --

library(readxl)

# Index Date
ALL_PUBLIC <- ALL_PUBLIC[VSIMPLE_INDEX_MASTER %in% PREDICT$VSIMPLE_INDEX_MASTER]
ALL_PUBLIC <- ALL_PUBLIC[, study_index_date := PREDICT$study_index_date[match(VSIMPLE_INDEX_MASTER, PREDICT$VSIMPLE_INDEX_MASTER)]] 

# - CVD - 
VIEW_CVD_ICD  <- read_xlsx("V:/common_lookups/Definitions/VIEW_CVD_ICD10_II_28JUN18_Clean042.xlsx")

# Retrieve ICD10 Codes for Events
# Nb: PREDICT datasets must always have broad CVD, diabetes, AF, and heart failure
cvd.var.names <- c("hx_broad_cvd", "hx_heart_failure", "mortality_broad_cvd_with_other", "out_broad_cvd")

for(var in cvd.var.names){
 
 var.codes <- VIEW_CVD_ICD$CLINICALCODE[which(eval(VIEW_CVD_ICD[,var])=="Y")]
 
 assign(gsub("_", ".", var), var.codes)
 
}

# - DM - 
VIEW_DM_ICD   <- read_xlsx("V:/common_lookups/Definitions/VIEW_CVD_ICD10_II_28JUN18_Clean042.xlsx", sheet = 2)

out.any.diabetes  <- VIEW_DM_ICD$CLINICALCODE[which(VIEW_DM_ICD$out_diabetes=="Y")]

# Type 2 
VIEW_DM_ICD9     <- read_xlsx("V:/bwu009/Definitions/2017_Diabetes_ICD_Typ1_2_U.xlsx", sheet = "ICD_9_Diab")
VIEW_DM_ICD10v1  <- read_xlsx("V:/bwu009/Definitions/2017_Diabetes_ICD_Typ1_2_U.xlsx", sheet = "ICD_10v1_Diab")
VIEW_DM_ICD10v2  <- read_xlsx("V:/bwu009/Definitions/2017_Diabetes_ICD_Typ1_2_U.xlsx", sheet = "ICD_10v2_Diab")
VIEW_DM_ICD10v3  <- read_xlsx("V:/bwu009/Definitions/2017_Diabetes_ICD_Typ1_2_U.xlsx", sheet = "ICD_10v3_Diab")

# Convert ICD9 to ICD10
ICD9_10 <- read_xlsx("V:/common_lookups/ICD 10 Mapping/ICD9 - 10/ICD9_to_10_forward_key.xlsx")

icd9.t2dm <- VIEW_DM_ICD9$clinical_code[which(VIEW_DM_ICD9$diabetes_type == 2)]

icd100.t2dm <- ICD9_10$ICD10[which(ICD9_10$ICD9 %in% icd9.t2dm)]
icd101.t2dm <- VIEW_DM_ICD10v1$clinical_code[which(VIEW_DM_ICD10v1$diabetes_type == 2)]
icd102.t2dm <- VIEW_DM_ICD10v2$clinical_code[which(VIEW_DM_ICD10v2$diabetes_type == 2)]
icd103.t2dm <- VIEW_DM_ICD10v3$clinical_code[which(VIEW_DM_ICD10v3$diabetes_type == 2)]

out.t2.diabetes   <- unique(c(icd100.t2dm, icd101.t2dm, icd102.t2dm, icd103.t2dm))

# - Other -
# Incorporates all codes within subcategory
hx.af            <- "^I48"
hx.gest.diabetes <- "^O2441"
hx.polycys.ovary <- "^E282"
hx.rheuma.arthri <- "^M06"
hx.asthma        <- "^J45"
hx.hyperthyroid  <- "^E059"
hx.hep.c         <- "^B182" 

subcat.codes <- c("hx.af", "hx.gest.diabetes", "hx.polycys.ovary", "hx.rheuma.arthri", "hx.asthma", "hx.hyperthyroid", "hx.hep.c")


# -- 2.  All-time History --
HX_EVENTS <- copy(ALL_PUBLIC)[EVSTDATE <= study_index_date]

history.vars <- ls(pattern = "^hx.")

source("V:/bwu009/Predict/Linkage/Any History v2.R")

# -- 3.  Outcomes --
OUT_EVENTS <- copy(ALL_PUBLIC)[EVSTDATE > study_index_date]

outcome.vars <- c("out.broad.cvd", "out.any.diabetes", "out.t2.diabetes")

source("V:/bwu009/Predict/Linkage/Any Outcomes v1.R")

# -- 4.  Admissions 28 days prior to death --
# Captures CVD hospitalisation within 28 days of death occurring
# - required for improved fatal CVD / improved fatal DM

# Add procedures to Broad CVD Definition
out.pvd.procs <- VIEW_CVD_ICD$CLINICALCODE[which(VIEW_CVD_ICD$out_pvd_procs == "Y")]
out.pci.cabg  <- VIEW_CVD_ICD$CLINICALCODE[which(VIEW_CVD_ICD$out_pci_cabg == "Y")]

fatal.broad.cvd    <- c(out.broad.cvd, out.pvd.procs, out.pci.cabg)
fatal.any.diabetes <- out.any.diabetes

out.28d.vars <- c("fatal.broad.cvd", "fatal.any.diabetes")

# Add info - date of death
ALL_PUBLIC[, "en_dod" := PREDICT$view_ag_dod[match(VSIMPLE_INDEX_MASTER, PREDICT$VSIMPLE_INDEX_MASTER)]]

# Find events
# The following conditions are extremely important to include - ALL must be included!
#   - DOD less than "2016-12-31 (ensure all deaths occur within study period and not after)
#   - Event admission greater than index date
#   - Event admission is within 28 days of death
EVENTS_28D_PRIOR <- copy(ALL_PUBLIC)[en_dod <= "2018-12-31" & EVSTDATE > study_index_date
][(EVSTDATE-en_dod >= -28) & (EVSTDATE <= en_dod)]

source("V:/bwu009/Predict/Linkage/28 Day Rule v1.R")

# -- 5.  Mortality --

# CVD/non-CVD & DM/non-DM causes of death
# Define death codes and categories: Capture all ICD codes for CVD and DM Deaths
mortality.any.diabetes  <- out.any.diabetes

death.vars  <- c("mortality.broad.cvd.with.other", "mortality.any.diabetes")

source("V:/bwu009/Predict/Linkage/Cause Specific Death v1.R")

# -- 6.  Merge all to Baseline --
Baseline <- Reduce(function(...)
 merge(..., 
       by="VSIMPLE_INDEX_MASTER",
       all.x = T),
 list(PREDICT, ALL_HX_VARS, ALL_OUT_VARS, ALL_28D_OUT_VARS[,.(VSIMPLE_INDEX_MASTER, fatal_broad_cvd, fatal_any_diabetes)], ALL_FATAL_VARS))


write.fst(Baseline, "V:/rpyl001/DM Risk Prediction/Working/PartB.fst", 75)


# ---- C. History of Mental Illness ----

# 1. Capture Activites

ALL_PRIMHD <- read.fst("V:/source_data/R/PREDICT/2018/Working/PP_Link_ALL_PRIMHD_v1.fst", 
                       as.data.table = T)

ALL_PRIMHD[, study_index_date := Baseline$study_index_date[match(VSIMPLE_INDEX_MASTER, Baseline$VSIMPLE_INDEX_MASTER)]]

# Definitions
inpatient   <- c("T02","T03","T04","T11","T12","T13","T14","T21")
resident    <- c("T05", "T16", "T20", "T27", "T28", "T29", "T30", "T48")
treatment   <- c("T01", "T09", "T22", "T36", "T38", "T39", "T40", "T41", "T42", "T17", "T18", "T19")
support     <- c("T07", "T22", "T23", "T24", "T43", "T44", "T45")

Activity.types <- c(inpatient, resident, treatment, support)

HX_ACTIVITY <- copy(ALL_PRIMHD)[ACTIVITY_TYPE_CODE %in% Activity.types & ACTIVITY_START_DATE <= study_index_date]

Baseline[, hx_primhd_activity := +(VSIMPLE_INDEX_MASTER %in% HX_ACTIVITY$VSIMPLE_INDEX_MASTER)]


# 2.  Mental disorders
PRIMHD_CLASS <- read.fst("source_data/R/PRIMHD/Classification Lookup.fst",
                         as.data.table = T)

PRIMHD_CLASS <- PRIMHD_CLASS[VSIMPLE_INDEX_MASTER %in% Baseline$VSIMPLE_INDEX_MASTER
][, study_index_date := Baseline$study_index_date[match(VSIMPLE_INDEX_MASTER, Baseline$VSIMPLE_INDEX_MASTER)]
][CLASSIFICATION_START_DATE <= study_index_date]

# Definitions
schiz_non_effect_dsm  <- paste0("^", c(2950:2959, 2970:2979, 2980:2989))
bipolar_effective_dsm <- paste0("^", c(2960:2961, 2964:2969))
dementia_org_dis_dsm  <- paste0("^", c(2900:2909, 2930:2949, 3100:3102, 3310:3333))
substance_use_dis_dsm <- paste0("^", c(2910:2929, 3030:3059))
depress_anxiety_dsm   <- paste0("^", c(2962:2963, 3000:3009, 3080:3099, "3110"))
dev_child_menret_dsm  <- paste0("^", c(2990:2998, 3120:3190))
other_mental_dis_dsm  <- paste0("^", c(3010:3019, 3020:3029, 3060:3078))

# ICD Codes:
schiz_non_effect_icd  <- paste0("^", paste0("F",20:29))
bipolar_effective_icd <- paste0("^", paste0("F",30:31))
dementia_org_dis_icd  <- paste0("^", c(paste0("F0",0:9), "F59", "G30"))
substance_use_dis_icd <- paste0("^", c(paste0("F",10:19), "F55"))
depress_anxiety_icd   <- paste0("^", paste0("F",32:49))
dev_child_menret_icd  <- paste0("^", paste0("F",70:98))
other_mental_dis_icd  <- paste0("^", c(paste0("F",50:54), paste0("F",60:69), "F99"))

Disorder.classes <- c("schiz_non_effect", "bipolar_effective", "dementia_org_dis", "substance_use_dis", 
                      "depress_anxiety", "dev_child_menret", "other_mental_dis")

for(i in Disorder.classes){
 
 dsm.codes <- get(paste0(i, "_dsm"))
 icd.codes <- get(paste0(i, "_icd"))
 
 CLASS_DATA_DSM <- PRIMHD_CLASS[grepl(paste0(dsm.codes, collapse = "|"), CLINICAL_CODE) & CLINICAL_CODING_SYSTEM_ID==7]
 CLASS_DATA_ICD <- PRIMHD_CLASS[grepl(paste0(icd.codes, collapse = "|"), CLINICAL_CODE) & CLINICAL_CODING_SYSTEM_ID!=7]
 
 CLASS_DATA <- as.data.table(rbind(CLASS_DATA_DSM, CLASS_DATA_ICD))
 
 CLASS_DATA <- CLASS_DATA[order(CLASSIFICATION_START_DATE),
                          index := seq_len(.N), 
                          by=VSIMPLE_INDEX_MASTER][index == 1, .(VSIMPLE_INDEX_MASTER, index)]
 
 setnames(CLASS_DATA, "index", i)
 
 ALL_CLASS <- if(i == Disorder.classes[1]){
  merge(Baseline[,.(VSIMPLE_INDEX_MASTER)], CLASS_DATA,
        by="VSIMPLE_INDEX_MASTER",
        all.x=T)
 } else {
  merge(ALL_CLASS, CLASS_DATA,
        by="VSIMPLE_INDEX_MASTER",
        all.x=T)
 }
 
 print(paste0(i, " completed")); rm(dsm.codes, icd.codes, CLASS_DATA_DSM, CLASS_DATA_ICD, CLASS_DATA, i)
 
}

ALL_CLASS[, hx_primhd_mental_disorder := apply(.SD, 1, 
                                               function(x)
                                                +(any(x == 1, na.rm = T))),
          .SDcols = Disorder.classes]

Baseline <- merge(Baseline, ALL_CLASS[,.(VSIMPLE_INDEX_MASTER, hx_primhd_mental_disorder)],
                  by = "VSIMPLE_INDEX_MASTER",
                  all.x = T)

write.fst(Baseline, "V:/rpyl001/DM Risk Prediction/Working/PartC.fst", 75)


# ---- D.   Pharms ----

PREDICT <- read.fst("V:/rpyl001/DM Risk Prediction/Working/PartA.fst",
                    as.data.table = T)

# files.path <- "V:/source_data/R/PREDICT/2018/National Collection/"
# files.list <- list.files(files.path, pattern = "PHH")
# 
# ALL_PHARMS <- rbindlist(lapply(files.list,
#                                function(x){
#                                  DATA <- read.fst(paste0(files.path, x), as.data.table = TRUE)
#                                  DATA[VSIMPLE_INDEX_MASTER %in% PREDICT$VSIMPLE_INDEX_MASTER]
#                                }))
# 
# ALL_PHARMS <- ALL_PHARMS[year(DATE_DISPENSED) %in% 2005:2018]
# 
# write.fst(ALL_PHARMS, "V:/rpyl001/DM Risk Prediction/Working/ALL_PHARMS_v1.fst", 75)

ALL_PHARMS <- read.fst("V:/rpyl001/DM Risk Prediction/Working/ALL_PHARMS_v1.fst", 
                       as.data.table = T)

# 1.  Definition

# VIEW standard
DIM_FORM_PACK <- read.fst("V:/common_lookups/CURRENT_PHARMS_LOOKUP_VIEW.fst", 
                          as.data.table = T)

vars.needed <- c("antiplatelets", "anticoagulants", "lipid_lowering", "bp_lowering", 
                 "antidiabetes", "insulin", "metformin", "other_oralhypos",
                 "metolazone",  "loopdiuretics", "antianginals", "corticosteroid")

view.varname <- c("all.antipla", "all.anticoag", "all.llds", "all.bplds", 
                  "all.diabetes", "insulin", "metformin", "other_oralhypos", 
                  "metolazone", "loop.diuretics", "antianginals", "corticosteroid")

setnames(DIM_FORM_PACK, vars.needed, view.varname)

# Capture Form Pack IDs (Lookup method)
for(class in view.varname){
 
 class.codes <- eval(substitute(
  DIM_FORM_PACK$DIM_FORM_PACK_SUBSIDY_KEY[which(DIM_FORM_PACK[, class, with = F]==1)]
 ))
 
 assign(class, class.codes)  
 
}

# Antipsychotics
ANTI_PSYCHOTICS <- readxl::read_xlsx("V:/common_lookups/Definitions/antipsychotics_ valproate_clozapine_olanzapine.xlsx", col_names = F)

var.names <- c("DIM_FORM_PACK_SUBSIDY_KEY", "CHEMICAL_ID", "CHEMICAL_NAME", "TG_NAME1", "TG_NAME2","TG_NAME3", 
               "FORMULATION_ID", "FORMULATION_NAME", "NA1", "NA2", "type")

setnames(ANTI_PSYCHOTICS, names(ANTI_PSYCHOTICS), var.names)

atypical_oral   <- unique(ANTI_PSYCHOTICS$DIM_FORM_PACK_SUBSIDY_KEY[which(ANTI_PSYCHOTICS$type==1)])
typical_oral    <- unique(ANTI_PSYCHOTICS$DIM_FORM_PACK_SUBSIDY_KEY[which(ANTI_PSYCHOTICS$type==2)])
atypical_depot  <- unique(ANTI_PSYCHOTICS$DIM_FORM_PACK_SUBSIDY_KEY[which(ANTI_PSYCHOTICS$type==3)])
typical_depot   <- unique(ANTI_PSYCHOTICS$DIM_FORM_PACK_SUBSIDY_KEY[which(ANTI_PSYCHOTICS$type==4)])
lithium         <- unique(ANTI_PSYCHOTICS$DIM_FORM_PACK_SUBSIDY_KEY[which(ANTI_PSYCHOTICS$type==5)])
navalproate     <- unique(ANTI_PSYCHOTICS$DIM_FORM_PACK_SUBSIDY_KEY[which(ANTI_PSYCHOTICS$type==6)]) # Antiepilepsy Drugs
clozapine       <- unique(ANTI_PSYCHOTICS$DIM_FORM_PACK_SUBSIDY_KEY[which(ANTI_PSYCHOTICS$type==11)])
olanzapine      <- unique(ANTI_PSYCHOTICS$DIM_FORM_PACK_SUBSIDY_KEY[which(ANTI_PSYCHOTICS$type==12)])
olanzapine_dept <- unique(ANTI_PSYCHOTICS$DIM_FORM_PACK_SUBSIDY_KEY[which(ANTI_PSYCHOTICS$type==32)])
ruth_exclusion  <- unique(ANTI_PSYCHOTICS$DIM_FORM_PACK_SUBSIDY_KEY[which(ANTI_PSYCHOTICS$type==99)])

atypical_any    <- c(atypical_oral, atypical_depot)
typical_any     <- c(typical_oral, typical_depot) 
olanzapine_any  <- c(olanzapine, olanzapine_dept) 

antipsy.mx <- c("atypical_any", "typical_any", "lithium", "navalproate", "clozapine", "olanzapine_any", "ruth_exclusion")

# Finalise / tidy
all.mx.groups <- c(view.varname, antipsy.mx)

rm(list=setdiff(ls(), c("Baseline","ALL_PHARMS", "all.mx.groups", all.mx.groups)))


# 2. Drug Capture 
ALL_PHARMS[, study_index_date := Baseline$study_index_date[match(VSIMPLE_INDEX_MASTER, Baseline$VSIMPLE_INDEX_MASTER)]]

# Prior dispensing
PRIOR_DISPENSING <- copy(ALL_PHARMS)[DATE_DISPENSED <= study_index_date]

prior.mx.groups <- c("all.antipla", "all.anticoag", "all.llds", "all.bplds",
                     "metolazone", "loop.diuretics", "antianginals", "corticosteroid",
                     "atypical_any", "typical_any", "lithium", "navalproate", "clozapine", 
                     "olanzapine_any", "ruth_exclusion")

source("V:/bwu009/Predict/Linkage/Prior Dispensing v1.R")

# Post dispensing
# First dispensing of diabetes medication
POST_DISPENSING <- copy(ALL_PHARMS)[DATE_DISPENSED > study_index_date]

post.mx.groups <- c("all.diabetes", "insulin", "metformin", "other_oralhypos")

source("V:/bwu009/Predict/Linkage/Post Dispensing v1.R")


# 3. Merge
Baseline <- Reduce(function(...)
 merge(..., 
       by="VSIMPLE_INDEX_MASTER",
       all.x = T),
 list(Baseline, ALL_PRIOR_DRUGS, ALL_POST_DRUGS))

# Quality Check
# Ensure no record of dispensing of any drug > 28 days after date of death 
ALL_PHARMS[, dod := Baseline$view_ag_dod[match(VSIMPLE_INDEX_MASTER, Baseline$VSIMPLE_INDEX_MASTER)]]

DOD_PHARMS <- ALL_PHARMS[!is.na(dod)]

# Capture last dispensing
DOD_PHARMS <- DOD_PHARMS[order(DATE_DISPENSED, decreasing = T), 
                         index := seq_len(.N), 
                         by="VSIMPLE_INDEX_MASTER"
][index == 1]

DOD_PHARMS <- DOD_PHARMS[, exclude := +(DATE_DISPENSED > dod+28)][exclude == 1]

# Merge flag to Baseline & exclude
# Removes 110 people
Baseline[, phh_post_dod := DOD_PHARMS$exclude[match(VSIMPLE_INDEX_MASTER, DOD_PHARMS$VSIMPLE_INDEX_MASTER)]]

Baseline <- Baseline[phh_post_dod==0 | is.na(phh_post_dod)][, -"phh_post_dod", with=F] # -110 


write.fst(Baseline, "V:/rpyl001/DM Risk Prediction/Working/PartD.fst", 75)


# ---- E.  Lab Values ----








