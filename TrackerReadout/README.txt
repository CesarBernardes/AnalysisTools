###For 2023 data

ssh -XY caber@lxplus8.cern.ch

export SCRAM_ARCH=el8_amd64_gcc11

voms-proxy-init -voms cms

###only once

cmsrel CMSSW_13_2_4

cd CMSSW_13_2_4/src

cmsenv

git cms-addpkg RecoLocalTracker/SiStripZeroSuppression

git cms-addpkg RecoLocalTracker/SiStripClusterizer

###Get the files "SiStripCMNanalyzer.cc", "master_RECO.py" and "plot_CMN_dependence_on_LumiSec.C"

###"master_RECO.py" is used to produce the Tree (please, use proper GT, etc..)

###"plot_CMN_dependence_on_LumiSec.C" is used to plot CMN vs LumiSec for a given FED (the one blocking the run)  

###copy "SiStripCMNanalyzer.cc" in "RecoLocalTracker/SiStripZeroSuppression/test/"


### To calculate APV CM shifts - 2024 data

ssh -XY caber@lxplus9.cern.ch
export SCRAM_ARCH=el9_amd64_gcc12
cmsrel CMSSW_14_1_4
cd CMSSW_14_1_4/src
cmsenv

### The code below calculates the CMN shifts per APV and also performs hybrid ZS emulation on VR data to calculate number of Bad APVs, for example

git cms-addpkg RecoLocalTracker/SiStripZeroSuppression

## Copy the module of this repo "SiStripHybridFormatAnalyzer.cc" in "RecoLocalTracker/SiStripZeroSuppression/test/"

scram b -j8

##first produce tree with APV CM for each APV. Use the config of this repo "config_for_apv_cmn_shifts_FEDId_FEDCh_APVId_AverageAPVCMN.py". The input file is in tracker VR mode.

cmsRun config_for_apv_cmn_shifts_FEDId_FEDCh_APVId_AverageAPVCMN.py >& OutPut.txt &

##From the output of the config above, print the per-APV CM shifts to be applied (shifts by 2-sigmas). Use the macro from this repo "print_list_average_APV_CMN.C". With this macro, you can check the value of the sigma of the CM of the distributions.

root -l -b -q print_list_average_APV_CMN.C











