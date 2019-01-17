
## Estimate separate Mixed Membership models


We make use of the <tt>mixedMem</tt> library to estimate the models, and on the
<tt>boot</tt> one for bootstrap confidence interval.


```r
library(mixedMem)
library(boot)
```
The actual code to initialize object need for estimation can be found in the <tt>MM_macro.R</tt>, which is called here via <tt>source()</tt> function


## Scenario 1


```r
## load data
load('simulated_data/Sim1_data_and_true_parameters.Rdata')
## call script that estimate the model
source('MM_macro.R')

## save results
Scenario1 = list(mixedMem1     = out1,
                 lambda.point1 = lambda.point1,
                 mixedMem2     = out2,
                 lambda.point2 = lambda.point2,
                 boot_res1     = boot_res1,
                 boot_res2     = boot_res2)


save(Scenario1, file = 'results/Sim1_MM_results.Rdata')

cat('SCENARIO 1 --- DONE \n\n')
rm(list = ls());gc()
```

## Scenario 2


```r
## load data
load('simulated_data/Sim2_data_and_true_parameters.Rdata')
## call script that estimate the model
source('MM_macro.R')

## save results
Scenario2 = list(mixedMem1     = out1,
                 lambda.point1 = lambda.point1,
                 mixedMem2     = out2,
                 lambda.point2 = lambda.point2,
                 boot_res1     = boot_res1,
                 boot_res2     = boot_res2)
                 
                 


save(Scenario2, file = 'results/Sim2_MM_results.Rdata')

cat('SCENARIO 2 --- DONE \n\n')
rm(list = ls());gc()
```


## Scenario 3


```r
## load data
load('simulated_data/Sim3_data_and_true_parameters.Rdata')
## call script that estimate the model
source('MM_macro.R')

## save results
Scenario3 = list(mixedMem1     = out1,
                 lambda.point1 = lambda.point1,
                 mixedMem2     = out2,
                 lambda.point2 = lambda.point2,
                 boot_res1     = boot_res1,
                 boot_res2     = boot_res2)
                 


save(Scenario3, file = 'results/Sim3_MM_results.Rdata')

cat('SCENARIO 3 --- DONE \n\n')
rm(list = ls());gc()
```


## Scenario 4


```r
## load data
load('simulated_data/Sim4_data_and_true_parameters.Rdata')
## call script that estimate the model
source('MM_macro.R')

## save results
Scenario4 = list(mixedMem1     = out1,
                 lambda.point1 = lambda.point1,
                 mixedMem2     = out2,
                 lambda.point2 = lambda.point2,
                 boot_res1     = boot_res1,
                 boot_res2     = boot_res2)


save(Scenario4, file = 'results/Sim4_MM_results.Rdata')

cat('SCENARIO 4 --- DONE \n\n')
rm(list = ls());gc()
```

## Misspecification: more than two pure types

The code macro <tt>MM_macro_miss.R</tt>, acts as <tt>MM_macro4.R</tt>, but it also returns the correct specified MM model with H=4.


```r
load('simulated_data/data_and_true_parameters_misspecified.Rdata')


source('MM_macro_miss.R')

## save results
Scenario_miss = list(mixedMem1     = out1,
                     lambda.point1 = lambda.point1,
                     mixedMem2     = out2,
                     lambda.point2 = lambda.point2,
                     mixedMem4     = out4,
                     lambda.point4 = lambda.point4)

save(Scenario_miss, file = 'results/Sim_miss_MM_results.Rdata')
```


