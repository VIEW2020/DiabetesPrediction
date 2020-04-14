### Exclude IF
- Prior admission for diabetes; OR
-	Treated with antidiabetic drugs in last 6 months; OR
-	Noted as diabetic in PREDICT; OR
- Non-existing HbA1c test in prior 2 years; OR
-	Any HbA1c in prior 2 years / post 30 days >= 50mmol/mol


### HbA1c
HbA1c data pulled from all repeated PREDICT assessments & all available TestSafe records. Must ensure filtering by group with one or more tests >= 50mmol/mol
```r
LAB_DM[LAB_DM[, .I[any(RESULT_MM >= 50)], by = VSIMPLE_INDEX_MASTER]$V1]
```

### Rollaround
Examine each visit and apply exclusion criteria. Continue with each nth visit to find those who could meet the inclusion criteria in subsequent visits.
- Starting with 564180 participants
- 394265 met criteria after first visit
- in remaining 169915, 9386 meet criteria in subsequent visit
- Initial criteria met by 403651

### Further Removals
- Age Limits 25-74: -16186
- Remove non-avaliable HbA1c within +30 days of Predict: -79915
- Remaining = 309475

## To Do
see issue #1
