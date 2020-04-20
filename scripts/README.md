### Exclude IF
- Prior admission for diabetes; OR
-	Treated with antidiabetic drugs in last 6 months; OR
-	Noted as diabetic in PREDICT; OR
- Non-existing HbA1c test in prior 2 years; OR
-	Any HbA1c in prior 2 years >= 50mmol/mol
- Exclusion critiera met within 30 days of study index date

### HbA1c
HbA1c data pulled from all repeated PREDICT assessments & all available TestSafe records. Must ensure filtering by group with one or more tests >= 50mmol/mol
```r
LAB_DM[LAB_DM[, .I[any(RESULT_MM >= 50)], by = VSIMPLE_INDEX_MASTER]$V1]
```

### Roll to next
Examine each visit and apply exclusion criteria. Continue with each nth visit to find those who could meet the inclusion criteria in subsequent visits.
- Starting with 564180 participants
- 371362 met criteria after first visit
- 12191 meet criteria in subsequent visit

### Further Removals
- Remove out of age Limits 25-74: -14806
- Remove non-avaliable HbA1c within +30 days of Predict: -75764
- Remove dispensing records 1 month beyond DOD :-107
- Remove exclusions detected within 30 days of study index: 190
- Total Remaining = 277880

### Preparation Scripts
Linkage preparation scripts <a href="https://github.com/VIEW2020/Predict/tree/master/Linkage" target="_blank">
located in Predict repository
</a>

