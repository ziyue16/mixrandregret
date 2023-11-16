## *```mixrandregret```*: A command for fitting mixed random regret minimization models using Stata 

```
*! randregret 1.1.0 18Mar2022
*! [aut & dev] Ziyue Zhu  &  Álvaro A. Gutiérrez-Vargas

/*********************************************************************************
  ____  __  _   _    ___   ____          __    ___   ___  __   ___   ___  ____  
 / / /  /    \ /    /__/  ____/  /\  /  /  \  /__/  /__  / _  /__/  /__    /  
/ / / _/_   _/ \_  /  \  /___/  /  \/  /___/ /  \  /__  /__/ /  \  /__    /   
 
 version 1.1.0:  gf0 ml evaluator that run mixed RRM model
		
	-> mixed Random Regret Model (Hensher et al., 2016)
	
*********************************************************************************/
```

```mixrandregret``` command utilizes the mixed random regret minimization model described in Hensher et al. (2016), which is a mixed version of the classic random regret minimization model introduced in Chorus. C. (2010). ```mixrandregret``` extends the ```randregret```  (Gutiérrez-Vargas et al, 2021) and allows to specify normal and log-normally distributed taste parameters inside the regret function. The command uses maximum simulated likelihood for estimation (Train. K., 2003). Users can obtain the predicted probabilities from the model using the ```mixrpred``` command. Individual-level parameters can also be obtained using the ```mixrbeta``` command.


### UPDATE: ```mixrandregret``` is uploaded to the SSC Archive

``` 
ssc install mixrandregret
```


### Install ```mixrandregret``` 

``` 
*Describe mixrandregret
net describe mixrandregret, from("https://raw.githubusercontent.com/ziyue16/mixrandregret/master/src/")


*Install mixrandregret
cap ado uninstall mixrandregret
net install mixrandregret, from("https://raw.githubusercontent.com/ziyue16/mixrandregret/master/src/")
```


### Examples 

Consider the following toy example that contains the data for two individuals who makes two choices. On each choice occasion, the individual has three alternatives. choice is the dependent variable, and x_rnd and x_fix are the independent variables or alternative attributes:

```
        id_ind id_cs choice_set alt  x_rnd   x_fix  choice
          1      1       1       1   -4.17   -0.83    0
          1      1       1       2    1.39   -0.87    0
          1      1       1       3    2.63   -0.37    1
          1      2       2       1   -0.45    1.12    1
          1      2       2       2    2.98    1.07    0
          1      2       2       3    1.37    0.75    0
          2      3       1       1    2.88   -0.25    0
          2      3       1       2   -3.12   -0.52    1
          2      3       1       3   -1.44    0.39    0
          2      4       2       1   -0.56    0.51    0
          2      4       2       2   -1.03    0.02    0
          2      4       2       3    0.99   -2.22    1
 ```     

A mixed random regret model where x_rnd has a normally distributed coefficient and x_fix has a fixed coefficient, and ASCs are suppressed can be specified as follows:
```
  . mixrandregret choice x_fix, group(id_cs) id(id_ind) rand(x_rnd) alternatives(alt) nocons
```
A model where x_rnd has a lognormally distributed coefficient can be specified as follows:
```
  . mixrandregret choice x_fix, group(id_cs) id(id_ind) rand(x_rnd) alternatives(alt) ln(1) nocons
```
Use alternative 1 as base alternative to calculate ASC can be specified as follows:
```
  . mixrandregret choice x_fix, group(id_cs) id(id_ind) rand(x_rnd) alternatives(alt) ln(1) basealternative(1)
```

```mixrandregret``` does't estimate random regret model to get inital values in mixed version. Users are advised to estimate on their own to make comparsion.


### References 

> * Chorus. C. 2010.  A New Model of Random Regret Minimization.  European Journal of Transport and Infrastructure Research 10: pp. 181-196.

> * Hensher, David A and Greene, William H and Ho, Chinh Q. 2016.  Random regret minimization and random utility maximization in the presence of preference heterogeneity: an empirical contrast.  Journal of Transportation Engineering Volume 142 Number 4: pp. 04016009.

> * Chorus. C. 2010.  A New Model of Random Regret Minimization.  European Journal of Transport and Infrastructure Research 10: pp. 181-196.

> * Gutiérrez-Vargas. Á., Meulders. M., and Vandebroek. M. 2021.  randregret: A command for fitting random regret minimization models using Stata.
        The Stata Journal Volume 21 Number 3: pp. 626-658.

> * Train. K. 2003.  Discrete Choice Methods with Simulation.  Cambridge University Press.

