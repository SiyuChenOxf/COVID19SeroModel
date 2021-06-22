# COVID19SeroModel
We have developed a method that combines data on daily death, seroprevalence and case positivity rates in England with a mechanistic mathematical model to infer the temporal trends of exposure and seroprevalence during the COVID-19 epidemic. We fit the mathematical model jointly to serological survey data available in seven regions in England namely, London, North West, North East, South East, South West, Midlands and East of England using a statistical observation model.  

## Constant infection fatality ratio (IFR)

### Input:
* regional death data
* regional seroprevalence data adjusted by sensitivity and specificity of the antibody test

### Output:
* Model predicted time-varying seroprevalence by regions 
* Model predicted time-varying exposure by regions 
* Model predicted constant infection fatality ratio by regions

## Time-varying infection fatality ratio (IFR)
### Input:
* regional death data
* regional virus positivity data
* regional seroprevalence data adjusted by sensitivity and specificity of the antibody test

### Output:
* Model predicted time-varying seroprevalence by regions 
* Model predicted time-varying exposure by regions 
* Model predicted time-varying infection fatality rate by regions

