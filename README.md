# COVID19SeroModel
We have developed a method that combines data on daily death, seroprevalence and virus test positivity rates in England with a mechanistic mathematical model to infer the temporal trends of exposure and seroprevalence during the COVID-19 epidemic. We fit the mathematical model jointly to serological survey data available in seven regions in England namely, London, North West, North East, South East, South West, Midlands and East of England using a Bayesian statistical observation model.  

## Constant infection fatality ratio (IFR) Model

### Input:
* regional death data
* regional seroprevalence data adjusted by sensitivity and specificity of the antibody test

### Output:
* Model predicted time-varying seroprevalence by regions 
* Model predicted time-varying exposure by regions 
* Model predicted constant infection fatality ratio by regions

## Time-varying infection fatality ratio (IFR) Model

### Input:
* regional death data
* regional virus test positivity data
* regional seroprevalence data adjusted by sensitivity and specificity of the antibody test

### Output:
* Model predicted time-varying seroprevalence by regions 
* Model predicted time-varying exposure by regions 
* Model predicted time-varying infection fatality rate by regions

## Sensitivity analysis

To ascertain the robustness of our main results, we explore how our estimates change as we:
* use different values for the delay between testing PCR positive and death , and for the delay between infection and death
* use a different prior distribution for seroreversion rate
* use a different set of mortality data 