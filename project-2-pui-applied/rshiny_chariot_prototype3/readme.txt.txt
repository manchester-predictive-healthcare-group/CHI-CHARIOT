The application currently consider pre-defined interventions statins, antihypertensives and smoking cessation. 

Other pre-defined interventions can be added, if their effect on SBP, Non-HDL cholesterol, BMI and smoking status is known or can be estimated.

It also allows undefined interventions, however we assume that the expected effect on each of SBP, Non-HDL cholesterol, BMI and smoking status is known. 

For this reason, we do not currently apply the knock-on' effect of changes in the modifiable risk factors. For example, based on the DAG, we expect a change in BMI
to also have a result on Non-HDL cholesterol and blood pressure, however we do apply these changes if undefined intervention is selected and only BMI is changed.
Future iterations of the front end may 'suggest' these as expected values, but then allow them to be further changed by the user, however this will require interactive input/output
which this rshiny does not have.