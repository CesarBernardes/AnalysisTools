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
