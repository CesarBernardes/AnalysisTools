ssh -XY caber@lxplus.cern.ch
cmssw-el5
export SCRAM_ARCH=slc5_amd64_gcc434
cmsrel CMSSW_4_4_7
cd CMSSW_4_4_7/src
cmsenv

###Copy the whole folder "HiForest" in the "src" folder and compile
cp -rf HiForest .
scram b -j4

###to run ntuple producer do:
cmsRun hiforestanalyzer_UCC_cfg.py >& OutPut.txt &

###Notes: 1) event selection done in config level, but trigger selection was done in the c++ level. Somehow the filters at config level were not working. PU filter as defined in 2011 is also applied in the c++ level

PDs used are:
/HIMinBiasUPC/HIRun2011-22Feb2012-v2/RECO (available in disk at T2_SPRACE)
/HIMinBiasUPC/HIRun2011-22Feb2012-v1/RECO (on tape only)

MC simulation (available in CMS Open Data):
https://opendata.cern.ch/record/25001

