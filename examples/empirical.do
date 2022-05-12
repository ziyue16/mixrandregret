*! 1.1.0 10May 2022
* Author: Ziyue ZHU (ziyue.zhu@student.kuleuven.be) 

clear all
global route = "D:\mixrandregret\src"
cd "$route"
findfile mixrandregret.ado

* use the data file from: Oluoch-Aridi J, Adam MB, Wafula F, et al. Understanding what women  want: eliciting preference for delivery health facility in a rural subcounty in Kenya, a discrete choice experiment. BMJ Open 2020;10:e038865. doi:10.1136/bmjopen-2020-038865

import excel "D:\mixrandregret\certification\empirical test\data\BMJOPEN_Datafor_DRYAD.xlsx", sheet("Naivasha database") firstrow clear

* clean the dataset
drop AM-DN

* fixed variable: cost
* random variable: The five variables that described the attributes of place of delivery
** QualityClinicalcare 
** attitudeofhealthworkers
** Medicalequipmentandsupplies
** distance
** referral services

timer clear 
* mixed random regret *
timer on 1
mixrandregret choice cost, nrep(2000) group(ncs) id(respondent_id) rand(badqualityclinicalservices distancetofacilitylong attitudeunkindnotsupportive medicalequipmentdrugsnotavail referralservicesnotavailable) alt(alt) base(3) tech(bfgs)

timer off 1

mixrpred pregret, xb
mixrpred p_r, proba

mixrbeta badqualityclinicalservices distancetofacilitylong attitudeunkindnotsupportive medicalequipmentdrugsnotavail referralservicesnotavailable, plot saving(betas_r) replace

* mixed logit *
timer on 2

mixlogit choice cost, nrep(2000) group(ncs) id(respondent_id) rand( badqualityclinicalservices distancetofacilitylong attitudeunkindnotsupportive medicalequipmentdrugsnotavail referralservicesnotavailable)

timer off 2
timer list

mixlpred p_l
mixlbeta badqualityclinicalservices distancetofacilitylong attitudeunkindnotsupportive medicalequipmentdrugsnotavail referralservicesnotavailable, saving(betas_l) replace

save predicted_data.dta, replace
