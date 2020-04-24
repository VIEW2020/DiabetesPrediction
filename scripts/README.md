### HbA1c
HbA1c data pulled from all repeated PREDICT assessments & all available TestSafe records. Must ensure filtering by group with one or more tests >= 50mmol/mol
```r
LAB_DM[LAB_DM[, .I[any(RESULT_MM >= 50)], by = VSIMPLE_INDEX_MASTER]$V1]
```

### Preparation Scripts
Linkage preparation scripts <a href="https://github.com/VIEW2020/Predict/tree/master/Linkage" target="_blank">
located in Predict repository
</a>

