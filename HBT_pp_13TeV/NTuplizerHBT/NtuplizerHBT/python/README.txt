===>>>RUNNING ON RECO
==============
***Important files:

-->Configuration files to set the analysis (in NTuplizerHBT/NtuplizerHBT/python):

ConfFile_cfg.py, ConfFile_forCrab_MC_cfg.py, ConfFile_forCrab_Data_cfg.py, ConfFile_forCrab_DataHM_ntrackoffline*_cfg.py


You just need to change the "first section" (sure, when running directly in cmssw you need to choose the files):
  --you can choose if run on data, which selection, which mixing and if save tree.

For crab:
 crabConfig_HBTanalysis13TeV_MC.py, crabConfig_HBTanalysis13TeV_Data.py
(these two files will use ConfFile_forCrab_MC_cfg.py, ConfFile_forCrab_Data_cfg.py, ConfFile_forCrab_DataHM_ntrackoffline*_cfg.py)
There we have the samples we are currently using for the analysis.
In general you need to change the parameters: 
config.General.requestName
config.Data.inputDataset
config.Data.unitsPerJob      
config.Data.totalUnits       
config.Data.outLFNDirBase 
 

--Analysis module (in NTuplizerHBT/NtuplizerHBT/plugins):

NtuplizerHBT.cc

Here are defined the selections, signal and reference samples, and the histograms 
(please check histograms names in function "NtuplizerHBT::initHistos(const edm::Service<TFileService> & fs)")


==============
***How to run:

Three Forms (in  NTuplizerHBT/NtuplizerHBT/python):

1. stand-alone (No crab)
cmsRun ConfFile_cfg.py


2. crab (for data)
crab submit -c crabConfig_HBTanalysis13TeV_Data.py


3. crab (for mc)
crab submit -c crabConfig_HBTanalysis13TeV_MC.py


==============
***Output:

--Three directories with histograms:
1. /demo/TrackHistograms
 --variables for control plots

2. /demo/MixedHistograms
 --1D/2D/3D histos for q_inv,nCh,kT for mixing reference samples 

3. /demo/qinvHistograms
 --1D/2D/3D histos for q_inv,nCh,kT for signal, inverted, rotated, opposite sign reference samples


--If you chose to save tree in the config. file, in addition to the 
histograms you will have a tree






===>>>RUNNING ON GEN-ONLY
==============
***Important files:





