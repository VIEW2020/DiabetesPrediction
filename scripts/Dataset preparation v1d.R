# Diabetes Prediction for Primary Prevention v1d
# Create a diabetes-free population
# 
# Last updated April 2020 - created by Billy Wu
# 
# Exclude IF
# Prior admission for diabetes; OR
# Treated with antidiabetic drugs in last 6 months; OR
# Noted as diabetic in PREDICT; OR
# Non-existing HbA1c test in prior 2 years; OR
#	Any HbA1c in prior 2 years / post 30 days >= 50mmol/mol

# Further refinement:
# There is a chance that individual is not in TestSafe catchment at index time point.
# 

library(fst)
library(data.table)

# ---- A. Cohort Development ----
# Find first Predict in which inclusion is met 

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
write.fst(ADM_DM, "V:/rpyl001/DiabetesPrediction/temp/ADM_DM.fst", 75)
write.fst(MX_DM, "V:/rpyl001/DiabetesPrediction/temp/MX_DM.fst", 75)
write.fst(TS_DM, "V:/rpyl001/DiabetesPrediction/temp/TS_DM.fst", 75)

#--- Import Data --
ALL_PREDICT <- read.fst("V:/source_data/R/PREDICT/2018/Cleaned_PREDICT_2018_All_Records_v1.fst",
                        as.data.table = T)

ADM_DM <- read.fst("V:/rpyl001/DiabetesPrediction/temp/ADM_DM.fst", 
                   as.data.table = T)

MX_DM <- read.fst("V:/rpyl001/DiabetesPrediction/temp/MX_DM.fst", 
                  as.data.table = T)

TS_DM <- read.fst("V:/rpyl001/DiabetesPrediction/temp/TS_DM.fst", 
                  as.data.table = T)

#--- FUN: FindExclusion ---
FindExclusion <- function(x){
 
 HX  <- copy(ADM_DM)[VSIMPLE_INDEX_MASTER %in% x$VSIMPLE_INDEX_MASTER
 ][, visit_date := x$view_visit_date[match(VSIMPLE_INDEX_MASTER, x$VSIMPLE_INDEX_MASTER)]]
 
 MX <- copy(MX_DM)[VSIMPLE_INDEX_MASTER %in% x$VSIMPLE_INDEX_MASTER
 ][, visit_date := x$view_visit_date[match(VSIMPLE_INDEX_MASTER, x$VSIMPLE_INDEX_MASTER)]]
 
 TS <- copy(TS_DM)[VSIMPLE_INDEX_MASTER %in% x$VSIMPLE_INDEX_MASTER
 ][, visit_date := x$view_visit_date[match(VSIMPLE_INDEX_MASTER, x$VSIMPLE_INDEX_MASTER)]]
 
 PT <- copy(ALL_PREDICT)[VSIMPLE_INDEX_MASTER %in% x$VSIMPLE_INDEX_MASTER
 ][, visit_date := view_visit_date]

 HX <- HX[EVSTDATE < visit_date]
 MX <- MX[(DATE_DISPENSED - visit_date >= -182) & (DATE_DISPENSED <= visit_date)]
 
 # Lab values - 1) Ensure that all prior values are normal; 2) ensure that there is a value within the last 2 years (used as index)
 TS <- TS[RESULT_DATE - visit_date <= 30]
 PT <- PT[pt_en_hba1c_mmdate - visit_date <= 30]
 
 setnames(TS, "EN_OBSR_RESULT_NUM", "RESULT_MM")
 setnames(PT, "pt_en_hba1c_mm", "RESULT_MM")
 setnames(PT, "pt_en_hba1c_mmdate", "RESULT_DATE")
 
 LAB <- rbind(TS[,.(VSIMPLE_INDEX_MASTER, visit_date, RESULT_MM, RESULT_DATE)], 
              PT[,.(VSIMPLE_INDEX_MASTER, visit_date, RESULT_MM, RESULT_DATE)])
 
 LAB <- LAB[LAB[, .I[all(RESULT_MM < 50)], by = VSIMPLE_INDEX_MASTER]$V1] # removes anyone with 1 or more high result
 LAB <- LAB[RESULT_DATE - visit_date >= -730]
 
 # In short: creates an index of people who qualify for inclusion
 summary(LAB$RESULT_MM) # all normal
 summary(as.numeric(LAB$RESULT_DATE - LAB$visit_date)) # all within 2 years / 30 days
 
 # Remove anyone how had DM in admission, treatment, or Predict
 LAB <- LAB[!VSIMPLE_INDEX_MASTER %in% 
             c(HX$VSIMPLE_INDEX_MASTER, MX$VSIMPLE_INDEX_MASTER, x$VSIMPLE_INDEX_MASTER[which(x$pt_diabetes %in% 1:3)])]
 
 x$exclude <- +(!x$VSIMPLE_INDEX_MASTER %in% LAB$VSIMPLE_INDEX_MASTER)
 
 return(x$exclude)
 
}

# --- Apply FUN as List Operation ---
# Examine each visit and apply exclusion criteria. Continue with each nth visit to find those who could meet the inclusion criteria in subsequent visits.
# Nb: Starting with 564180 participants, 371362 met criteria after first visit; 12191 meet criteria in subsequent visits.
#     Total 383553 met the criteria

# Limit PREDICT to 2006 - 2018
ALL_PREDICT <- ALL_PREDICT[year(ALL_PREDICT$view_visit_date) >= 2006]

PREDICT_LS <- lapply(split(ALL_PREDICT, by = "view_visit_seq"), 
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
                                    by = VSIMPLE_INDEX_MASTER]$V1] 

# Apply Age Limits: 25-74
PREDICT_BL <- PREDICT_BL[view_ag_age %in% 25:74] # -14806

# 368747 met the criteria
write.fst(PREDICT_BL, "V:/rpyl001/DiabetesPrediction/temp/PREDICT_BL_v1d.fst", 75)


# --- Index HbA1c ---

# PREDICT_BL <- read.fst("V:/rpyl001/DiabetesPrediction/temp/PREDICT_BL_v1d.fst", as.data.table = T)

# Draws from both PREDICT and TestSafe
# Ensure that there is a value within the last 2 years / 30 days
# Capture the nearest HbA1c result as index

TS <- copy(TS_DM)[VSIMPLE_INDEX_MASTER %in% PREDICT_BL$VSIMPLE_INDEX_MASTER
][, visit_date := PREDICT_BL$view_visit_date[match(VSIMPLE_INDEX_MASTER, PREDICT_BL$VSIMPLE_INDEX_MASTER)]]

PT <- copy(PREDICT)[VSIMPLE_INDEX_MASTER %in% PREDICT_BL$VSIMPLE_INDEX_MASTER
][, visit_date := PREDICT_BL$view_visit_date[match(VSIMPLE_INDEX_MASTER, PREDICT_BL$VSIMPLE_INDEX_MASTER)]]

setnames(TS, "EN_OBSR_RESULT_NUM", "RESULT_MM")
setnames(PT, "pt_en_hba1c_mm", "RESULT_MM")
setnames(PT, "pt_en_hba1c_mmdate", "RESULT_DATE")

TS <- TS[RESULT_DATE - visit_date >= -730 & RESULT_DATE - visit_date <= 30]
PT <- PT[RESULT_DATE - visit_date >= -730 & RESULT_DATE - visit_date <= 30]

INDEX <- rbind(TS[,.(VSIMPLE_INDEX_MASTER, visit_date, RESULT_MM, RESULT_DATE)], 
               PT[,.(VSIMPLE_INDEX_MASTER, visit_date, RESULT_MM, RESULT_DATE)])

# In short: creates an index of people who qualify for inclusion
summary(INDEX$RESULT_MM) # all normal
summary(as.numeric(INDEX$RESULT_DATE - INDEX$visit_date)) # all within 2 years / 30 days

# Find Nearest HbA1c value as index
library(Allie)

INDEX[, by = VSIMPLE_INDEX_MASTER
      , flag := FindNearestDate(index = visit_date, 
                                comparison = RESULT_DATE, 
                                mode = "n",
                                from = -730,
                                to = 30)]

INDEX <- INDEX[flag == 1]
INDEX <- INDEX[INDEX[, .I[which.min(RESULT_MM)], by = VSIMPLE_INDEX_MASTER]$V1]

setnames(INDEX, "RESULT_MM", "imp_hba1c_2y30d_result")
setnames(INDEX, "RESULT_DATE", "imp_hba1c_2y30d_resultdate")

# Merge to PREDICT Baseline Data
PREDICT_BL <- merge(INDEX[,.(VSIMPLE_INDEX_MASTER, index_hba1c_mm, index_hba1c_mmdate)],
                    PREDICT_BL,
                    by = "VSIMPLE_INDEX_MASTER",
                    all.x = F)

summary(PREDICT_BL$index_hba1c_mm) # all normal
summary(as.numeric(PREDICT_BL$index_hba1c_mmdate - PREDICT_BL$view_visit_date)) # all within 2 years / 30 days

write.fst(PREDICT_BL, "V:/rpyl001/DiabetesPrediction/temp/PREDICT_BL_v1d.fst", 75) #278177 remaining

rm(list = ls())

# ---- B.   History / Outcomes ----

PREDICT <- read.fst("V:/rpyl001/DiabetesPrediction/temp/PREDICT_BL_v1d.fst",
                    as.data.table = T)

# Loaded pre-prepared national collection 
files.path <- "V:/source_data/R/PREDICT/2018/National Collection/"
files.list <- list.files(files.path, pattern = "PUBLIC")

ALL_PUBLIC <- rbindlist(lapply(files.list, 
                               function(x)
                                read.fst(paste0(files.path, x), 
                                         as.data.table = T)))

ALL_DEATHS <- read.fst(paste0(files.path, "VSIMPLE_DEATHS_v1.fst"), as.data.table = T)


# -- 1.   Definitions --

library(readxl)

# Index Date
ALL_PUBLIC <- ALL_PUBLIC[VSIMPLE_INDEX_MASTER %in% PREDICT$VSIMPLE_INDEX_MASTER]
ALL_PUBLIC <- ALL_PUBLIC[, study_index_date := PREDICT$view_visit_date[match(VSIMPLE_INDEX_MASTER, PREDICT$VSIMPLE_INDEX_MASTER)]] 

# - CVD - 
VIEW_CVD_ICD  <- read_xlsx("V:/common_lookups/Definitions/VIEW_CVD_ICD10_II_28JUN18_Clean042.xlsx")

# Retrieve ICD10 Codes for Events
# Nb: PREDICT datasets must always have broad CVD, diabetes, AF, and heart failure
cvd.var.names <- c("hx_broad_cvd", "hx_heart_failure", "hx_unst_angina", "mortality_broad_cvd_with_other", "out_broad_cvd")

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

# - Cancer - 
hx.cancer <- readxl::read_excel("V:/common_lookups/Definitions/cancer_code_table_8th_STATUS_A_ONLY.xlsx")$CLINICAL_CODE

# - Transplatation / Dialysis- 
VIEW_RENAL    <- readxl::read_xlsx("V:/common_lookups/Definitions/VIEW_CVD_ICD10_II_28JUN18_Clean042.xlsx", sheet = "Renal")

hx.renal.dialysis.transp <- VIEW_RENAL$CLINICALCODE[which(VIEW_RENAL$hx_renal_dialysis_transp == "Y")]

# - Other -
# Incorporates all codes within subcategory
hx.af            <- "^I48"
hx.gest.diabetes <- "^O2441"
hx.polycys.ovary <- "^E282"
hx.rheuma.arthri <- "^M06"
hx.gout          <- "^M10"
hx.copd          <- "^J44" 
hx.asthma        <- "^J45"
hx.hyperthyroid  <- "^E05"
hx.hep.c         <- "^B182"

subcat.codes <- c("hx.af", "hx.gest.diabetes", "hx.polycys.ovary", "hx.rheuma.arthri", "hx.gout", "hx.copd", "hx.asthma", "hx.hyperthyroid", "hx.hep.c")

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

# -- 6. History of Mental Illness --

# 1. Capture Activites
ALL_PRIMHD <- read.fst("V:/source_data/R/PREDICT/2018/Working/PP_Link_ALL_PRIMHD_v1.fst", 
                       as.data.table = T)

ALL_PRIMHD[, study_index_date := PREDICT$view_visit_date[match(VSIMPLE_INDEX_MASTER, PREDICT$VSIMPLE_INDEX_MASTER)]]

# Definitions
inpatient   <- c("T02","T03","T04","T11","T12","T13","T14","T21")
resident    <- c("T05", "T16", "T20", "T27", "T28", "T29", "T30", "T48")
treatment   <- c("T01", "T09", "T22", "T36", "T38", "T39", "T40", "T41", "T42", "T17", "T18", "T19")
support     <- c("T07", "T22", "T23", "T24", "T43", "T44", "T45")

Activity.types <- c(inpatient, resident, treatment, support)

HX_ACTIVITY <- copy(ALL_PRIMHD)[ACTIVITY_TYPE_CODE %in% Activity.types & ACTIVITY_START_DATE <= study_index_date]

PREDICT[, hx_primhd_activity := +(VSIMPLE_INDEX_MASTER %in% HX_ACTIVITY$VSIMPLE_INDEX_MASTER)]


# 2.  Mental disorders
PRIMHD_CLASS <- read.fst("source_data/R/PRIMHD/Classification Lookup.fst",
                         as.data.table = T)

PRIMHD_CLASS <- PRIMHD_CLASS[VSIMPLE_INDEX_MASTER %in% PREDICT$VSIMPLE_INDEX_MASTER
][, study_index_date := PREDICT$view_visit_date[match(VSIMPLE_INDEX_MASTER, PREDICT$VSIMPLE_INDEX_MASTER)]
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

source("V:/bwu009/Predict/Linkage/Mental Disorder v1.R")

# -- 7. Cancer Register --
# NB: Improvement required with NMDS
CANCER <- read.fst("V:/source_data/R/CANCER/CANCER_REGISTRY_1995_2014.fst",
                   as.data.table = T)

# Registry
HX_CANCER_REG <- CANCER[, study_index_date := PREDICT$view_visit_date[match(VSIMPLE_INDEX_MASTER, PREDICT$VSIMPLE_INDEX_MASTER)]
                        ][diag_date <= study_index_date]

HX_CANCER_REG <- HX_CANCER_REG[, hx_cancer_reg := 1][, unique(.SD), .SDcols = c("VSIMPLE_INDEX_MASTER", "hx_cancer_reg")]

ALL_CANCER_REG <- merge(PREDICT[,.(VSIMPLE_INDEX_MASTER)], HX_CANCER_REG,
                        by = "VSIMPLE_INDEX_MASTER",
                        all.x = T)

ALL_CANCER_REG[, hx_cancer_reg := +(!is.na(hx_cancer_reg))]


# ---- C.   Pharms ----

# PREDICT <- read.fst("V:/rpyl001/DiabetesPrediction/temp/PREDICT_BL_v1d.fst",
#                     as.data.table = T)
# 
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
# write.fst(ALL_PHARMS, "V:/rpyl001/DiabetesPrediction/temp/ALL_PHARMS_v1d.fst", 75)

ALL_PHARMS <- read.fst("V:/rpyl001/DiabetesPrediction/temp/ALL_PHARMS_v1.fst", 
                       as.data.table = T)

# 1.  Definition

# VIEW standard
DIM_FORM_PACK <- read.fst("V:/common_lookups/CURRENT_PHARMS_LOOKUP_VIEW.fst", 
                          as.data.table = T)

vars.needed <- c("antiplatelets", "anticoagulants", "lipid_lowering", "bp_lowering", 
                 "antidiabetes", "insulin", "metformin", "other_oralhypos",
                 "metolazone",  "loopdiuretics", "antianginals", "corticosteroid",
                 "nonasp_nsaids")

view.varname <- c("all.antipla", "all.anticoag", "all.llds", "all.bplds", 
                  "all.diabetes", "insulin", "metformin", "other_oralhypos", 
                  "metolazone", "loop.diuretics", "antianginals", "corticosteroid",
                  "nonasp.nsaids")

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

antipsy.mx <- c("atypical_any", "typical_any", "lithium", "navalproate", "clozapine", "olanzapine_any")

# Rheumatoid arthritis
rheuma.chems <- c("Auranofin", "Hydroxychloroquine", "Leflunomide", "Penicillamine")

# Gout Specific (Chemical Name Method)
gout.chems <- c("allopurinol", "febuxostat", "benzbromarone", "colchicine")

for(chem in c(gout.chems, rheuma.chems)){
 
 form.packs <- DIM_FORM_PACK$DIM_FORM_PACK_SUBSIDY_KEY[which(tolower(DIM_FORM_PACK$CHEMICAL_NAME) %in% tolower(chem))]
 
 assign(tolower(chem), form.packs)
 
}

gout.specific  <- c(allopurinol, febuxostat, benzbromarone, colchicine)
antirheumatoid <- c(auranofin, hydroxychloroquine, leflunomide, penicillamine)

# COPD/Asthma Specific
# LAMA / LABA Definitions
Tiotropium <- 3805
Glycopyrronium <- 4043
Umeclidinium <- 4057

Salmeterol <- 1066
Eformoterol <- 1083
Indacaterol <- 4042

Salmeterol.fluticasone <- 3858
Eformoterol.budesonide <- 3758
Vilanterol.fluticasone <- 4056

Tiotropium.olodaterol <- 4059
Glycopyrronium.indacaterol <- 4058
Umeclidinium.vilanterol <- 4060

lama <- c(Tiotropium, Glycopyrronium, Umeclidinium)
laba <- c(Salmeterol, Eformoterol, Indacaterol)
labaics   <- c(Salmeterol.fluticasone, Eformoterol.budesonide, Vilanterol.fluticasone)
lamalaba  <- c(Tiotropium.olodaterol, Glycopyrronium.indacaterol, Umeclidinium.vilanterol)

copd.chems <- c("lama", "laba", "labaics", "lamalaba")

for(chem in copd.chems){
 
 form.packs <- DIM_FORM_PACK$DIM_FORM_PACK_SUBSIDY_KEY[which(DIM_FORM_PACK$CHEMICAL_ID %in% get(chem))]
 
 assign(tolower(chem), form.packs)
 
}

# Finalise / tidy
all.mx.groups <- c(view.varname, antipsy.mx, "gout.specific", "antirheumatoid", copd.chems)

# 2. Drug Capture 
ALL_PHARMS[, study_index_date := PREDICT$view_visit_date[match(VSIMPLE_INDEX_MASTER, PREDICT$VSIMPLE_INDEX_MASTER)]]

# Prior dispensing
PRIOR_DISPENSING <- copy(ALL_PHARMS)[DATE_DISPENSED <= study_index_date]

prior.mx.groups <- all.mx.groups

source("V:/bwu009/Predict/Linkage/Prior Dispensing v1.R")

# Post dispensing
# First dispensing of diabetes medication
POST_DISPENSING <- copy(ALL_PHARMS)[DATE_DISPENSED > study_index_date]

post.mx.groups <- c("all.diabetes", "insulin", "metformin", "other_oralhypos")

source("V:/bwu009/Predict/Linkage/Post Dispensing v1.R")


# -- 6.  Merge all to Baseline --
Baseline <- Reduce(function(...)
 merge(..., 
       by="VSIMPLE_INDEX_MASTER",
       all.x = T),
 list(PREDICT, ALL_HX_VARS, ALL_CANCER_REG, ALL_CLASS, ALL_OUT_VARS, ALL_28D_OUT_VARS[,.(VSIMPLE_INDEX_MASTER, fatal_broad_cvd, fatal_any_diabetes)], 
      ALL_FATAL_VARS, ALL_PRIOR_DRUGS, ALL_POST_DRUGS))


# Quality Check
# 1. Ensure no record of dispensing of any drug > 28 days after date of death

DOD_PHARMS <- ALL_PHARMS[, dod := Baseline$view_ag_dod[match(VSIMPLE_INDEX_MASTER, Baseline$VSIMPLE_INDEX_MASTER)]
                         ][!is.na(dod)]

DOD_PHARMS <- DOD_PHARMS[DOD_PHARMS[, .I[which.max(DATE_DISPENSED)], by = "VSIMPLE_INDEX_MASTER"]$V1]
DOD_PHARMS <- DOD_PHARMS[, exclude := +(DATE_DISPENSED > dod + 28)
                         ][exclude == 1] #97

Baseline <- Baseline[!VSIMPLE_INDEX_MASTER %in% DOD_PHARMS$VSIMPLE_INDEX_MASTER] #Removes 107 people

# 2. Make sure DODs do not exceed start or end of study!
sum(Baseline$view_ag_dod > "2018-12-31", na.rm = T)

# 3. Remove people with renal dialysis & transplantation: -826
Baseline <- Baseline[hx_renal_dialysis_transp == 0] 


write.fst(Baseline, "V:/rpyl001/DiabetesPrediction/temp/PartC_v1d.fst", 75)

rm(list = setdiff(ls(), "Baseline"))

# ---- D.  Lab Values ----

Baseline <- read.fst("V:/rpyl001/DiabetesPrediction/temp/PartC_v1d.fst",
                     as.data.table = T)

# Baseline eGFR, ACR, TC/HDL values
# Evaluate all values within 2yr + 30 day (as per index HbA1c)

ALL_PREDICT <- read.fst("V:/source_data/R/PREDICT/2018/Cleaned_PREDICT_2018_All_Records_v1.fst",
                        as.data.table = T)

# ALL_PREDICT <- ALL_PREDICT[VSIMPLE_INDEX_MASTER %in% Baseline$VSIMPLE_INDEX_MASTER]
# 
# # - Prepare Data -
# files.path     <- "V:/source_data/R/PREDICT/2018/TestSafe/"
# standard.names <- c("RESULT", "RESULT_DATE")
# 
# # i) SCR/eGFR
# TS_SCR <- rbindlist(lapply(list.files(files.path, pattern = "SCR"),
#                            function(x){
#                             DATA <- read.fst(paste0(files.path, x), as.data.table = TRUE)
#                             DATA[VSIMPLE_INDEX_MASTER %in% Baseline$VSIMPLE_INDEX_MASTER]
#                            }))
# 
# PT_SCR <- copy(ALL_PREDICT)[!is.na(pt_en_serum_creatinine),
#                             .(VSIMPLE_INDEX_MASTER, view_ag_sex, view_ag_dob, pt_en_serum_creatinine, pt_en_serum_creatininedate)]
# 
# # Capture eGFR-ckdepi
# egfr.demo <- c("SEX", "TEST_AGE", "NULL_ETH")
# 
# PT_SCR[, (egfr.demo) := list(
#  +(view_ag_sex == "M"),
#  as.numeric(floor((pt_en_serum_creatininedate - view_ag_dob) / 365.25)),
#  0
# )]
# 
# PT_SCR <- PT_SCR[, EN_OBSR_RESULT_EGFR := round(
#  nephro::CKDEpi.creat(pt_en_serum_creatinine * 0.0113, SEX, TEST_AGE, NULL_ETH),
#  2)]
# 
# setnames(PT_SCR, c("EN_OBSR_RESULT_EGFR", "pt_en_serum_creatininedate"), standard.names)
# setnames(TS_SCR, "EN_OBSR_RESULT_EGFR", "RESULT")
# 
# ALL_EGFR <- unique(rbind(PT_SCR[, c("VSIMPLE_INDEX_MASTER", standard.names), with = F],
#                          TS_SCR[, c("VSIMPLE_INDEX_MASTER", standard.names), with = F]))
# 
# # ii) TCHDL
# TS_TCHDL <- rbindlist(lapply(list.files(files.path, pattern = "TCHDL"),
#                              function(x){
#                               DATA <- read.fst(paste0(files.path, x), as.data.table = TRUE)
#                               DATA[VSIMPLE_INDEX_MASTER %in% Baseline$VSIMPLE_INDEX_MASTER]
#                              }))
# 
# PT_TCHDL <- copy(ALL_PREDICT)[!is.na(pt_en_tchdl_ratio),
#                               .(VSIMPLE_INDEX_MASTER, pt_en_tchdl_ratio, pt_en_tchdl_ratiodate)]
# 
# setnames(PT_TCHDL, c("pt_en_tchdl_ratio", "pt_en_tchdl_ratiodate"), standard.names)
# setnames(TS_TCHDL, "EN_OBSR_RESULT_NUM", "RESULT")
# 
# ALL_TCHDL <- unique(rbind(PT_TCHDL[, c("VSIMPLE_INDEX_MASTER", standard.names), with = F],
#                           TS_TCHDL[, c("VSIMPLE_INDEX_MASTER", standard.names), with = F]))
# 
# # iii) ACR
# TS_ACR <- read.fst("V:/source_data/R/TESTSAFE/ACR/VSIMP_ACR_2004_2018.fst",
#                    as.data.table = T)[VSIMPLE_INDEX_MASTER %in% Baseline$VSIMPLE_INDEX_MASTER]
# 
# PT_ACR <- copy(ALL_PREDICT)[!is.na(pt_en_diab_acr),
#                             .(VSIMPLE_INDEX_MASTER, pt_en_diab_acr, pt_en_diab_acrdate)]
# 
# setnames(PT_ACR, c("pt_en_diab_acr", "pt_en_diab_acrdate"), standard.names)
# setnames(TS_ACR, "EN_OBSR_RESULT_NUM", "RESULT")
# 
# ALL_ACR <- unique(rbind(PT_ACR[, c("VSIMPLE_INDEX_MASTER", standard.names), with = F],
#                         TS_ACR[, c("VSIMPLE_INDEX_MASTER", standard.names), with = F]))
# 
# # -- Save
# write.fst(ALL_TCHDL, "V:/rpyl001/DiabetesPrediction/temp/ALL_TCHDL_v1d.fst", 75)
# write.fst(ALL_ACR, "V:/rpyl001/DiabetesPrediction/temp/ALL_ACR_v1d.fst", 75)
# write.fst(ALL_EGFR, "V:/rpyl001/DiabetesPrediction/temp/ALL_EGFR_v1d.fst", 75)

#--  Fetch Data

ALL_TCHDL <- read.fst("V:/rpyl001/DiabetesPrediction/temp/ALL_TCHDL_v1d.fst", as.data.table = T)
ALL_ACR <- read.fst("V:/rpyl001/DiabetesPrediction/temp/ALL_ACR_v1d.fst", as.data.table = T)
ALL_EGFR <- read.fst("V:/rpyl001/DiabetesPrediction/temp/ALL_EGFR_v1d.fst", as.data.table = T)


# -- Find nearest index value -- 
# Uses -2yr / +30d as per index HbA1c

standard.names <- c("RESULT", "RESULT_DATE")
test.types <- c("TCHDL", "EGFR", "ACR")
 
for(test in test.types){
 
 DATA <- get(paste0("ALL_", test))
 DATA <- DATA[, visit_date := Baseline$view_visit_date[match(VSIMPLE_INDEX_MASTER, Baseline$VSIMPLE_INDEX_MASTER)]
              ][RESULT_DATE - visit_date >= -730 & RESULT_DATE - visit_date <= 30]
 
 DATA[, by = VSIMPLE_INDEX_MASTER
      , flag := FindNearestDate(index = visit_date, 
                                comparison = RESULT_DATE, 
                                mode = "n",
                                from = -730,
                                to = 30)]
 
 DATA <- DATA[flag == 1
              ][,.(VSIMPLE_INDEX_MASTER, RESULT, RESULT_DATE)]
 
 DATA <- DATA[DATA[, .I[which.min(RESULT)], by = VSIMPLE_INDEX_MASTER]$V1]
 
 test.names <- c(paste("imp", tolower(test), "2yr30d_result", sep = "_"),
                 paste("imp", tolower(test), "2yr30d_resultdate", sep = "_"))
 
 setnames(DATA, standard.names, test.names)
 
 assign(test, DATA)
 
 print(paste(test, "completed"))
 
} 
 
# - ACR special 18 month look forward - 
ACR18 <- ALL_ACR[, visit_date := Baseline$view_visit_date[match(VSIMPLE_INDEX_MASTER, Baseline$VSIMPLE_INDEX_MASTER)]
                ][RESULT_DATE - visit_date >= -730 & RESULT_DATE - visit_date <= 547]

ACR18[, by = VSIMPLE_INDEX_MASTER
     , flag := FindNearestDate(index = visit_date, 
                               comparison = RESULT_DATE, 
                               mode = "d",
                               from = -730,
                               to = 547)]

ACR18 <- ACR18[flag == 1
             ][,.(VSIMPLE_INDEX_MASTER, RESULT, RESULT_DATE)]

ACR18 <- ACR18[ACR18[, .I[which.min(RESULT)], by = VSIMPLE_INDEX_MASTER]$V1]

test.names <- c("imp_acr_2yr18m_result", "imp_acr_2yr18m_resultdate")

setnames(ACR18, standard.names, test.names)


# -- HbA1c Outcome --

# # Data Prep
# TS_HBA1C <- read.fst("V:/rpyl001/DiabetesPrediction/temp/TS_DM.fst", 
#                      as.data.table = T)
# 
# PT_HBA1C <- copy(ALL_PREDICT)[!is.na(pt_en_hba1c_mm),
#                               .(VSIMPLE_INDEX_MASTER, pt_en_hba1c_mm, pt_en_hba1c_mmdate)]
# 
# setnames(PT_HBA1C, c("pt_en_hba1c_mm", "pt_en_hba1c_mmdate"), standard.names)
# setnames(TS_HBA1C, "EN_OBSR_RESULT_NUM", "RESULT")
# 
# ALL_HBA1C <- unique(rbind(PT_HBA1C[, c("VSIMPLE_INDEX_MASTER", standard.names), with = F],
#                           TS_HBA1C[, c("VSIMPLE_INDEX_MASTER", standard.names), with = F]))
# 
# write.fst(ALL_HBA1C, "V:/rpyl001/DiabetesPrediction/temp/ALL_HBA1C_v1d.fst", 75)

ALL_HBA1C <- read.fst("V:/rpyl001/DiabetesPrediction/temp/ALL_HBA1C_v1d.fst", 
                      as.data.table = T)

ALL_HBA1C <- ALL_HBA1C[VSIMPLE_INDEX_MASTER %in% Baseline$VSIMPLE_INDEX_MASTER]

POST_HBA1C <- ALL_HBA1C[, visit_date := Baseline$view_visit_date[match(VSIMPLE_INDEX_MASTER, Baseline$VSIMPLE_INDEX_MASTER)]
                        ][RESULT >= 50 & RESULT_DATE - visit_date > 30]

POST_HBA1C <- POST_HBA1C[POST_HBA1C[, .I[which.min(RESULT_DATE)], by = VSIMPLE_INDEX_MASTER]$V1]
POST_HBA1C <- POST_HBA1C[, -c("visit_date")]

post.names <- c("imp_hba1c_post_50mm", "imp_hba1c_post_50mmdate")

setnames(POST_HBA1C, standard.names, post.names)

# Merge all
Baseline <- Reduce(function(...)
 merge(..., 
       by="VSIMPLE_INDEX_MASTER",
       all.x = T),
 list(Baseline, TCHDL, EGFR, ACR, ACR18, POST_HBA1C))

# Binary
summary(Baseline$imp_hba1c_post_50mm)
summary(as.numeric(Baseline$imp_hba1c_post_50mmdate - Baseline$view_visit_date))

Baseline[, imp_hba1c_post_50mm := +(!is.na(imp_hba1c_post_50mm))]


# save
write.fst(Baseline, "V:/rpyl001/DiabetesPrediction/temp/PartD_v1d.fst", 75)

rm(list = setdiff(ls(), "Baseline"))

# --- E. Improvements ----
Baseline <- read.fst("V:/rpyl001/DiabetesPrediction/temp/PartD_v1d.fst",
                     as.data.table = T)

ALL_PREDICT <- read.fst("V:/source_data/R/PREDICT/2018/Cleaned_PREDICT_2018_All_Records_v1.fst",
                        as.data.table = T)

# Predict DM outcome
ALL_PREDICT[, study_visit_date := Baseline$view_visit_date[match(VSIMPLE_INDEX_MASTER, Baseline$VSIMPLE_INDEX_MASTER)]]

PT_DIABETES <- copy(ALL_PREDICT)[!is.na(pt_diabetes) & view_visit_date > study_visit_date & (pt_diabetes == 2 | pt_diabetes == 3)
                                 ][, .(VSIMPLE_INDEX_MASTER, study_visit_date, view_visit_date, pt_diabetes)]

PT_DIABETES <- PT_DIABETES[PT_DIABETES[, .I[which.min(view_visit_date)], by = VSIMPLE_INDEX_MASTER] $V1]
PT_DIABETES[, c("pt_diabetes_post",
                "pt_diabetes_post_date") := .(1, view_visit_date)]

Baseline <- merge(Baseline, copy(PT_DIABETES[,.(VSIMPLE_INDEX_MASTER, pt_diabetes_post, pt_diabetes_post_date)]),
                  by = "VSIMPLE_INDEX_MASTER",
                  all.x =T)

Baseline[is.na(pt_diabetes_post), 
         pt_diabetes_post := 0]

# # QC: 23 people have a Predict diabetes status within 30 days of study index
# table(as.numeric(Baseline$pt_diabetes_post_date - Baseline$view_visit_date) <= 30)
# 
# # 149 people have a diabetes drug dispensing within 30 days of study index
# table(as.numeric(Baseline$ph_all_diabetes_post_date - Baseline$view_visit_date) <= 30)
# 
# table(as.numeric(Baseline$ph_insulin_post_date - Baseline$view_visit_date) <= 30) # insulin = 1
# table(as.numeric(Baseline$ph_metformin_post_date - Baseline$view_visit_date) <= 30) # metformin = 144
# table(as.numeric(Baseline$ph_other_oralhypos_post_date - Baseline$view_visit_date) <= 30) # other orals = 4
# 
# # 14 people have a DM admissions within 30 days of study index
# table(as.numeric(Baseline$out_any_diabetes_adm_date - Baseline$view_visit_date) <= 30)
# table(as.numeric(Baseline$out_t2_diabetes_adm_date - Baseline$view_visit_date) <= 30)
# 
# # 179 people have any of the above
# table(as.numeric(Baseline$out_any_diabetes_adm_date - Baseline$view_visit_date <30 |
#        Baseline$ph_all_diabetes_post_date - Baseline$view_visit_date <30 |
#        Baseline$pt_diabetes_post_date - Baseline$view_visit_date <30))

# - VDR Outcomes -
ALL_VDR <- read.fst("source_data/R/PREDICT/2018/Working/PP_Link_ALL_VDR_v1.fst",
                    as.data.table = T)

ALL_VDR[, study_index_date := Baseline$view_visit_date[match(VSIMPLE_INDEX_MASTER, Baseline$VSIMPLE_INDEX_MASTER)]]

POST_VDR <- copy(ALL_VDR)[first_identification_date - study_index_date >= 30]
POST_VDR <- POST_VDR[POST_VDR[, .I[which.min(first_identification_date)], by = VSIMPLE_INDEX_MASTER]$V1]
POST_VDR <- POST_VDR[, vdr_out_diabetes := 1
                     ][, .(VSIMPLE_INDEX_MASTER, vdr_out_diabetes, first_identification_date)]

setnames(POST_VDR, "first_identification_date", "vdr_out_diabetes_date")

Baseline <- merge(Baseline, POST_VDR,
              by = "VSIMPLE_INDEX_MASTER",
              all.x = T)

Baseline[is.na(vdr_out_diabetes),
         vdr_out_diabetes := 0]

# -- QC Exclusions-- 
Baseline[, exclude_lt30d := +(out_any_diabetes_adm_date - view_visit_date <30 | 
                               ph_all_diabetes_post_date - view_visit_date <30 |
                               pt_diabetes_post_date - view_visit_date <30 |
                               (imp_hba1c_post_50mmdate - view_visit_date <30 & imp_hba1c_post_50mm == 1))]

Baseline <- Baseline[exclude_lt30d == 0 | is.na(exclude_lt30d)] #-179

# Improved outcome of diabetes / time to diabetes event

t2dm.vars <- c("out_t2_diabetes", "ph_metformin_post", "ph_other_oralhypos_post", "pt_diabetes_post", "imp_hba1c_post_50mm")
anymx.vars <- c("out_t2_diabetes", "ph_all_diabetes_post", "pt_diabetes_post", "imp_hba1c_post_50mm")
anydm.vars <- c("out_any_diabetes", "ph_all_diabetes_post", "pt_diabetes_post", "imp_hba1c_post_50mm")
anyvdr.vars <- c("out_any_diabetes", "ph_all_diabetes_post", "pt_diabetes_post", "imp_hba1c_post_50mm", "vdr_out_diabetes")

t2dm.datevars <- c("out_t2_diabetes_adm_date", "ph_metformin_post_date", "ph_other_oralhypos_post_date", "pt_diabetes_post_date", "imp_hba1c_post_50mmdate")
anymx.datevars <- c("out_t2_diabetes_adm_date", "ph_all_diabetes_post_date", "pt_diabetes_post_date", "imp_hba1c_post_50mmdate")
anydm.datevars <- c("out_any_diabetes_adm_date", "ph_all_diabetes_post_date", "pt_diabetes_post_date", "imp_hba1c_post_50mmdate")
anyvdr.datevars <- c("out_any_diabetes_adm_date", "ph_all_diabetes_post_date", "pt_diabetes_post_date", "imp_hba1c_post_50mmdate", "vdr_out_diabetes_date")

imp.dm.flag <- c("imp_out_t2_diabetes", "imp_out_t2_diabetes_anymx", "imp_out_any_diabetes", "imp_out_any_diabetes_vdr")

Baseline[, (imp.dm.flag) := list(
 apply(.SD[,t2dm.vars, with = F], 1, function(x) +(any(x == 1, na.rm = T))),
 apply(.SD[,anymx.vars, with = F], 1, function(x) +(any(x == 1, na.rm = T))),
 apply(.SD[,anydm.vars, with = F], 1, function(x) +(any(x == 1, na.rm = T))),
 apply(.SD[,anyvdr.vars, with = F], 1, function(x) +(any(x == 1, na.rm = T)))
 )]

Baseline[imp_out_t2_diabetes == 1, 
     imp_out_t2_diabetes_date := as.Date(apply(.SD, 1, function(x) min(x, na.rm = T))),
     .SDcols = t2dm.datevars]

Baseline[imp_out_t2_diabetes_anymx == 1, 
     imp_out_t2_diabetes_anymx_date := as.Date(apply(.SD, 1, function(x) min(x, na.rm = T))),
     .SDcols = anymx.datevars]

Baseline[imp_out_any_diabetes == 1, 
     imp_out_any_diabetes_date := as.Date(apply(.SD, 1, function(x) min(x, na.rm = T))),
     .SDcols = anydm.datevars]

Baseline[imp_out_any_diabetes_vdr == 1, 
     imp_out_any_diabetes_vdr_date := as.Date(apply(.SD, 1, function(x) min(x, na.rm = T))),
     .SDcols = anyvdr.datevars]


# Improved Fatal CVD/DM
Baseline[, imp_fatal_cvd := +(fatal_broad_cvd == 1 | mortality_broad_cvd_with_other == 1)]
Baseline[, imp_fatal_diabetes := +(fatal_any_diabetes == 1 | mortality_any_diabetes == 1)]

# Improved hx CVD
imp.hx.vars <- c("imp_hx_af", "imp_hx_cvd_with_hf")

Baseline[, (imp.hx.vars) := list(+(pt_atrial_fibrillation==1 | hx_af==1),
                                 
                                 +(hx_heart_failure==1 | ph_loop_diuretics_prior_5yrs_3evts==1 | ph_metolazone_prior_6mths==1 | 
                                    pt_angina==1 | hx_unst_angina==1 | ph_antianginals_prior_5yrs_3evts==1 |
                                    hx_broad_cvd==1 | pt_mi==1 | pt_ihd==1 | pt_ptca_cabg==1 | pt_stroke==1 | pt_tia==1 | pt_pvd==1))]

# Improved hx CVD treatment 
imp.mx.vars <- c("imp_hx_antihypertensives", "imp_hx_antithrombotics", "imp_hx_lipidlowering")

Baseline[, (imp.mx.vars) := list(+(pt_en_ace_inhibitor==1 | pt_en_at2==1 | pt_en_beta_blocker==1 | pt_en_calcium_antagonist==1 | 
                                pt_en_other_hyp_drugs==1 | pt_en_thiazide==1 | ph_all_bplds_prior_6mths==1),
                             
                             +(ph_all_antipla_prior_6mths==1 | pt_en_aspirin==1 | pt_en_clopidogrel==1 |
                                ph_all_anticoag_prior_6mths==1 | pt_en_warfarin==1),
                             
                             +(ph_all_llds_prior_6mths==1 | pt_en_statin==1 | pt_en_fibrate==1 | pt_en_other_lipid_drugs==1))]

# Other improved vars
# LABA and LABA/ICS products are also used by people with asthma (from Otago HRC protocol)

imp.other.vars <- c("imp_hx_copd", "imp_hx_asthma", "imp_hx_gout", "hx_severe_mental_illness", "imp_hx_cancer", "imp_hx_rheumatoid")

Baseline[, (imp.other.vars) := list(+(hx_copd==1 | ph_lama_prior_6mths==1 |  ph_lamalaba_prior_6mths==1),
                                    +(hx_asthma==1 | ph_laba_prior_6mths==1 | ph_labaics_prior_6mths==1),
                                    +(hx_gout==1 | ph_gout_specific_prior_6mths==1),
                                    +(hx_primhd_activity==1 | hx_primhd_mental_disorder==1),
                                    +(hx_cancer_reg==1 | hx_cancer==1),
                                    +(hx_rheuma_arthri==1 | ph_antirheumatoid_prior_6mths==1))] 

# Binary
all.imps <- c(imp.hx.vars, imp.mx.vars, imp.other.vars)

Baseline[, (all.imps) := lapply(.SD, function(x)
 ifelse(is.na(x), 0, x)),
 .SDcols = all.imps]

# save
write.fst(Baseline, "V:/rpyl001/DiabetesPrediction/temp/PartE_v1d.fst", 75) #277075

rm(list = setdiff(ls(), "Baseline"))

# ---- F. Final Tidy up ----

library(fst)
library(dplyr)
library(data.table)

Baseline <- read.fst("V:/rpyl001/DiabetesPrediction/temp/PartE_v1d.fst",
                     as.data.table = T)

# Rename
setnames(Baseline, c("index_hba1c_mm", "index_hba1c_mmdate"), 
         c("imp_hba1c_2y30d_result", "imp_hba1c_2y30d_resultdate"))
setnames(Baseline, c("pt_diabetes_post", "pt_diabetes_post_date"),
         c("sub_pt_diabetes23_post", "sub_pt_diabetes23_post_date"))

# Remove redudant vars
ptrun.vars  <- names(Baseline)[startsWith(names(Baseline), "pt_run")]
ptdiab.vars  <- names(Baseline)[startsWith(names(Baseline), "pt_diab_")]
admin.vars  <- c("pt_current_hba1cdate", "pt_submissions_timestamp")
mort.vars   <- c("mortality_broad_cvd_with_other", "fatal_broad_cvd", "fatal_any_diabetes", "mortality_any_diabetes")

Baseline <- Baseline %>% 
 select(
  setdiff(
   names(Baseline), 
   c(mort.vars, ptrun.vars, ptdiab.vars, admin.vars)))

Baseline <- Baseline %>%
 select(VSIMPLE_INDEX_MASTER, view_visit_date,
        starts_with("view"),
        starts_with("pt"),
        hx_gest_diabetes, hx_hep_c, hx_hyperthyroid, hx_polycys_ovary, hx_severe_mental_illness,
        starts_with("out"),
        "sub_pt_diabetes23_post", "sub_pt_diabetes23_post_date",
        starts_with("vdr"),
        starts_with("ph"),
        starts_with("imp"))

# Remove all components used for improvement (i.e. keep only improved variable)
imp.components <- c("pt_atrial_fibrillation" , "ph_loop_diuretics_prior_5yrs_3evts" , "ph_metolazone_prior_6mths",
                    "pt_angina" , "ph_antianginals_prior_5yrs_3evts", "pt_mi" , "pt_ihd" , "pt_ptca_cabg" , "pt_stroke" , "pt_tia", "pt_pvd" , "pt_stroke_tia",
                    "pt_en_ace_inhibitor" , "pt_en_at2" , "pt_en_beta_blocker" , "pt_en_calcium_antagonist" , "pt_en_other_hyp_drugs" , 
                    "pt_en_thiazide" , "ph_all_bplds_prior_6mths", "ph_all_antipla_prior_6mths" , "pt_en_aspirin" , "pt_en_clopidogrel" ,
                    "ph_all_anticoag_prior_6mths" , "pt_en_warfarin", "ph_all_llds_prior_6mths" , "pt_en_statin" , "pt_en_fibrate" , "pt_en_other_lipid_drugs", 
                    "pt_diabetes", "pt_diabetes_yr", "pt_en_hba1c_mmdate", "pt_en_hba1cdate","pt_en_hba1c_mm", "pt_en_diab_acr","pt_en_diab_acrdate",
                    "pt_en_diab_diet_referral", "pt_en_diab_edu_referral", "pt_en_diab_eye_lastret", "pt_en_diab_feet_date_last_check",
                    "pt_en_serum_creatinine", "pt_en_serum_creatininedate",  "ph_all_diabetes_prior_6mths", "ph_insulin_prior_6mths",
                    "ph_metformin_prior_6mths", "ph_other_oralhypos_prior_6mths",
                    "ph_gout_specific_prior_6mths", "ph_antirheumatoid_prior_6mths",
                    "ph_lama_prior_6mths", "ph_laba_prior_6mths", "ph_labaics_prior_6mths", "ph_lamalaba_prior_6mths")

Baseline <- Baseline[, -imp.components, with=F]

names(Baseline)

# Save
write.fst(Baseline, "V:/rpyl001/DiabetesPrediction/data//PREDICT2018_DM_PP_v3.fst", 75)

# Save to Stata
library(haven)

write_dta(setDT(Baseline), "V:/rpyl001/DiabetesPrediction/data/PREDICT2018_DM_PP_v3.dta", version = 14)

