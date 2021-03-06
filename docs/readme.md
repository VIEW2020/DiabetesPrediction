### Approved DAP
DAP diabetes risk prediction.docx

### Definitions
- 2017_Diabetes_ICD_Typ1_2_U.xlsx (reviewed by Sue Weels 20 July 2017)
- VIEW_CVD_ICD10_II_28JUN18_Clean042.xlsx (latest version at 15 April 2020)
- cancer_code_table_8th_STATUS_A_ONLY.xlsx (reviewed by Vanessa Selak *date required*)
- ICD9_to_10_forward_key.xlsx (current MoH ICD9 to ICD10.v1 mapping key)

### Non-standard Classes
```r
#Incorporates all codes within subcategory
hx.af            <- "^I48"
hx.gest.diabetes <- "^O2441"
hx.polycys.ovary <- "^E282"
hx.rheuma.arthri <- "^M06"
hx.gout          <- "^M10"
hx.copd          <- "^J44" 
hx.asthma        <- "^J45"
hx.hyperthyroid  <- "^E05"
hx.hep.c         <- "^B182"
```

### Medication
```r
#COPD
lama <- c(Tiotropium, Glycopyrronium, Umeclidinium)
laba <- c(Salmeterol, Eformoterol, Indacaterol)
labaics   <- c(Salmeterol.fluticasone, Eformoterol.budesonide, Vilanterol.fluticasone)
lamalaba  <- c(Tiotropium.olodaterol, Glycopyrronium.indacaterol, Umeclidinium.vilanterol)
```

