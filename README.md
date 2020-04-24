# Type-II Diabetes Prediction Model

This study will provide validated T2D prediction models, relevant in contemporary New Zealand’s primary care. 

## Core Data
<a href="https://github.com/VIEW2020/DiabetesPrediction/wiki" target="_blank">See Wiki</a> for details regarding the core dataset including variable descriptions.<br></br>
The core dataset is a PREDICT population from 2006 - 2018. Individuals are eligeable if they are diabetes-free at time of PREDICT risk assessment and where an index HbA1c value is avaiable.
Inclusion does not limit the number of prior PREDICT assessments or prior CVD. Where eligeability is met at multiple PREDICT assessments, the earliest PREDICT record is used as baseline.

Information from the National Health Collection are linked to provide demographic, hospitalised history, hospitalised outcomes, death-specific outcomes, and baseline treatment. 
To ensure consistency, the exclusion criteria (see Wiki) have been applied in data management. 

### Exclude IF
- Prior admission for diabetes; OR
-	Treated with antidiabetic drugs in last 6 months; OR
-	Noted as diabetic in PREDICT; OR
- Non-existing HbA1c test in prior 2 years; OR
-	Any HbA1c in prior 2 years >= 50mmol/mol
- Exclusion critiera met within 30 days of study index date
- QC Conflicts

### Roll to next
Examine each visit and apply exclusion criteria. Continue with each nth visit to find those who could meet the inclusion criteria in subsequent visits.
- Starting with 564180 participants
- 371362 met criteria after first visit
- 12191 meet criteria in subsequent visit

### Further Removals
- Remove out of age Limits 25-74: -14806
- Remove non-avaliable HbA1c within +30 days of Predict: -75764
- Remove dispensing records 1 month beyond DOD :-107
- Remove people with renal dialysis & transplantation: -826
- Remove exclusions detected within 30 days of study index: -179
- Total Remaining = 277075
