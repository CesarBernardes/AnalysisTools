The macro in https://github.com/CesarBernardes/AnalysisTools/blob/main/HBT-ATLAS-Method/singleRatio_MC_ATLASmethod.C
will create the correcpondence between Non-femtoscopic parameters from opposite-charge pairs to same-charge pairs.

The values we extracted are:
1) For HIJING
log(B_A)_same-sign = (0.996 +- 0.017)*log(B_A)_opp-sign + (-0.247 +- 0.034)
sigma_B,A_same-sign = (1.211 +- 0.219)*sigma_B,A_opp-sign + (-0.062 +- 0.062)

2) For Pythia
log(B_A)_same-sign = (1.209 +- 0.041)*log(B_A)_opp-sign + (0.190 +- 0.098)
sigma_B,A_same-sign = (0.690 +- 0.163)*sigma_B,A_opp-sign + (0.079 +- 0.043)

Please, follow the sequence in data:

1) build opposite-sign single ratio and fit with the function here:
https://github.com/CesarBernardes/AnalysisTools/blob/main/HBT-ATLAS-Method/singleRatio_MC_ATLASmethod.C#L5-L23
* Parameter [1] is B_A (note that is not log and the relation above is log, so you need to convert it).
* Parameter [2] is sigma_B,A
Convert these parameters to same sign pairs by using the functions above. You can have results for Pythia and Hijing.

2) build same-sign single ratio and fit with a function of the form f = func_BE * func_NonFemto:
* func_BE is your Bose-Einstein function used in other methods
* func_NonFemto is same function as in https://github.com/CesarBernardes/AnalysisTools/blob/main/HBT-ATLAS-Method/singleRatio_MC_ATLASmethod.C#L5-L23 . The parameters of this function should be fixed by the values extracted from the "Step 1". 



