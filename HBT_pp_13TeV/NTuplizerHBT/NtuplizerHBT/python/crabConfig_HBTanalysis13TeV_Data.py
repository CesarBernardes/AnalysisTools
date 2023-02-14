###For a description of the crabConfig.py parameters. See:
###https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile

from WMCore.Configuration import Configuration
config = Configuration() ###create a Configuration object

config.section_('General')###add a new section of type "General"
###General: In this section, the user specifies generic parameters about the request (e.g. request name).
config.General.workArea     = 'HBTanalysispp13TeV_3D_Data_Projects' ###fixed name for projects dir in my area

#config.General.requestName  = 'hbtAnalysispp13TeV_3D_data_L1MinimumBiasHF1_SPRACEcuts_EtaMix_V1_RecoCorr' #sub dir with prefix "crab_". Change it for each task
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_data_L1MinimumBiasHF2_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_data_L1MinimumBiasHF3_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_data_L1MinimumBiasHF4_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_data_L1MinimumBiasHF5_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_data_L1MinimumBiasHF6_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_data_L1MinimumBiasHF7_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_data_L1MinimumBiasHF8_SPRACEcuts_EtaMix_V1_RecoCorr'
config.General.requestName  = 'hbtAnalysispp13TeV_3D_data_HighMultiplicity_ntrkoffline80to104_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_data_HighMultiplicity_ntrkoff105to129_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_data_HighMultiplicity_ntrkoff130_SPRACEcuts_EtaMix_V1_RecoCorr'


config.General.transferLogs = True 
config.General.transferOutputs = True

################################

config.section_('JobType')###add a new section of type "JobType"
###JobType: This section aims to contain all the parameters of the user job type and 
###related configurables (e.g. CMSSW parameter-set configuration file, additional input files, etc.).
config.JobType.pluginName     = 'Analysis'
#config.JobType.psetName       = 'ConfFile_forCrab_DataMinBias_cfg.py'
config.JobType.psetName       = 'ConfFile_forCrab_DataHM_ntrackoffline80to104_cfg.py'
#config.JobType.psetName       = 'ConfFile_forCrab_DataHM_ntrackoffline105to129_cfg.py'
#config.JobType.psetName       = 'ConfFile_forCrab_DataHM_ntrackoffline130_cfg.py'

#config.JobType.allowNonProductionCMSSW = True
config.JobType.maxMemoryMB    = 2500

config.JobType.inputFiles     = ['/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_ntrkOffline_Pythia8/results/trkCorr_forHBTAnalysis_tot_Pythia8.root','/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_ntrkOffline_EPOS/results/trkCorr_forHBTAnalysis_tot_EPOS.root','/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_spraceMinusChi2NdofCuts_Pythia8_V2/results/trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_V2.root','/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_spraceMinusChi2NdofCuts_EPOS/results/trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_EPOS.root','/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_6_4/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_spraceMinusChi2NdofCuts_Pythia6/results/trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia6.root']

#config.JobType.inputFiles     = ['/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_spraceMinusChi2NdofCuts_Pythia8_Tight_v2/results/trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Tight.root','/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_spraceMinusChi2NdofCuts_Pythia8_Loose/results/trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Loose.root','/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_spraceMinusChi2NdofCuts_Pythia8_Medium/results/trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Medium.root','/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_spraceMinusChi2NdofCuts_Pythia8_Zvtx3cm/results/trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Zvtx3cm.root','/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_spraceMinusChi2NdofCuts_Pythia8_3Zvtx15cm_v2/results/trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_3Zvtx15cm.root']

#config.JobType.inputFiles     = ['/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_spraceMinusChi2NdofCuts_Pythia8_Loose2/results/trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Loose2.root','/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_spraceMinusChi2NdofCuts_Pythia8_Medium2/results/trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Medium2.root','/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_spraceMinusChi2NdofCuts_Pythia8_Tight2/results/trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Tight2.root']



#config.JobType.inputFiles     = ['']##put new selections here

################################

config.section_('Data')###add a new section of type "Data"
###Data: This section contains all the parameters related to the data to be analyzed, 
###including the splitting parameters.
#config.Data.inputDataset      = '/L1MinimumBiasHF1/davidlw-RecoSkim2015_2015DLowPU_ReTracking_v4-6ca56d9a50b153298a3c6a7ed0fc5558/USER'
#config.Data.inputDataset      = '/L1MinimumBiasHF2/davidlw-RecoSkim2015_2015DLowPU_ReTracking_v4-6ca56d9a50b153298a3c6a7ed0fc5558/USER'
#config.Data.inputDataset      = '/L1MinimumBiasHF3/davidlw-RecoSkim2015_2015DLowPU_ReTracking_v4-6ca56d9a50b153298a3c6a7ed0fc5558/USER'
#config.Data.inputDataset      = '/L1MinimumBiasHF4/davidlw-RecoSkim2015_2015DLowPU_ReTracking_v4-6ca56d9a50b153298a3c6a7ed0fc5558/USER'
#config.Data.inputDataset      = '/L1MinimumBiasHF5/davidlw-RecoSkim2015_2015DLowPU_ReTracking_v4-6ca56d9a50b153298a3c6a7ed0fc5558/USER'
#config.Data.inputDataset      = '/L1MinimumBiasHF6/davidlw-RecoSkim2015_2015DLowPU_ReTracking_v4-6ca56d9a50b153298a3c6a7ed0fc5558/USER'
#config.Data.inputDataset      = '/L1MinimumBiasHF7/davidlw-RecoSkim2015_2015DLowPU_ReTracking_v4-6ca56d9a50b153298a3c6a7ed0fc5558/USER'
#config.Data.inputDataset      = '/L1MinimumBiasHF8/davidlw-RecoSkim2015_2015DLowPU_ReTracking_v4-6ca56d9a50b153298a3c6a7ed0fc5558/USER'
config.Data.inputDataset      = '/HighMultiplicity/davidlw-RecoSkim2015_2015DLowPU_ReTracking_v4-266f47bcc90a343001055b934437977e/USER'
config.Data.splitting         = 'FileBased'
#config.Data.unitsPerJob       = X ###files per job (but not impose)
#config.Data.totalUnits        = Y ###how many files to analyze
#config.Data.unitsPerJob       = 60  ##L1MinimumBiasHF1
#config.Data.totalUnits        = 1331
#config.Data.unitsPerJob       = 60 ##L1MinimumBiasHF2
#config.Data.totalUnits        = 1339
#config.Data.unitsPerJob       = 30 ##L1MinimumBiasHF3
#config.Data.totalUnits        = 670
#config.Data.unitsPerJob       = 30 ##L1MinimumBiasHF4
#config.Data.totalUnits        = 670
#config.Data.unitsPerJob       = 30 ##L1MinimumBiasHF5
#config.Data.totalUnits        = 670
#config.Data.unitsPerJob       = 30 ##L1MinimumBiasHF6
#config.Data.totalUnits        = 670
#config.Data.unitsPerJob       = 30 ##L1MinimumBiasHF7
#config.Data.totalUnits        = 662
#config.Data.unitsPerJob       = 30 ##L1MinimumBiasHF8
#config.Data.totalUnits        = 665
config.Data.unitsPerJob       = 5 ##High Multiplicity
config.Data.totalUnits        = 1330
config.Data.inputDBS          = 'phys03'
#config.Data.inputDBS          = 'global'
#config.Data.publishDBS = 'phys03'
#config.Data.publishDataName = 'TestingJobs_RecoWithNewPixelTracks'
#config.Data.outputDatasetTag = 'HBTAnalysispp13TeV_3D_MC_Test01' ###change for each sample(task) --only if publish

#config.Data.outLFN            = '/store/user/caber/TESTJOBS_RecoWithNewPixelTracks'

#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_Data_L1MinimumBiasHF1_SPRACEcuts_EtaMix_V1_RecoCorr' ###change for each sample(task)
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_Data_L1MinimumBiasHF2_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_Data_L1MinimumBiasHF3_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_Data_L1MinimumBiasHF4_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_Data_L1MinimumBiasHF5_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_Data_L1MinimumBiasHF6_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_Data_L1MinimumBiasHF7_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_Data_L1MinimumBiasHF8_SPRACEcuts_EtaMix_V1_RecoCorr'
config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_Data_HighMultiplicity_ntrkoff80to104_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_Data_HighMultiplicity_ntrkoff105to129_SPRACEcuts_EtaMix_V1_RecoCorr'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_Data_HighMultiplicity_ntrkoffline130_SPRACEcuts_EtaMix_V1_RecoCorr'



################################

config.section_('Site')###add a new section of type "Site"
###Site: Grid site parameters are defined in this section, including the stage out information 
###(e.g. stage out destination site, white/black lists, etc.).
config.Site.storageSite       = 'T2_BR_SPRACE'
#config.Site.whitelist         = ['T2_US_MIT']
