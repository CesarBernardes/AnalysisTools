###For a description of the crabConfig.py parameters. See:
###https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile

from WMCore.Configuration import Configuration
config = Configuration() ###create a Configuration object

config.section_('General')###add a new section of type "General"
###General: In this section, the user specifies generic parameters about the request (e.g. request name).
config.General.workArea     = 'HBTanalysispp13TeV_MC_Projects' ###fixed name for projects dir in my area

config.General.requestName  = 'hbtAnalysispp13TeV_mc_EPOS_MB_SPRACEcuts_EtaMix_NoTree_V10_RecoCorr' #sub dir with prefix "crab_". Change it for each task
#config.General.requestName  = 'hbtAnalysispp13TeV_mc_pythia8_MB_SPRACEcuts_EtaMix_NoTree_V10_RecoCorr'



config.General.transferLogs = True 
config.General.transferOutputs = True

################################

config.section_('JobType')###add a new section of type "JobType"
###JobType: This section aims to contain all the parameters of the user job type and 
###related configurables (e.g. CMSSW parameter-set configuration file, additional input files, etc.).
config.JobType.pluginName     = 'Analysis'
config.JobType.psetName       = 'ConfFile_forCrab_MC_cfg.py'

#config.JobType.allowNonProductionCMSSW = True
config.JobType.maxMemoryMB    = 2500

config.JobType.inputFiles     = ['/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_ntrkOffline_Pythia8/results/trkCorr_forHBTAnalysis_tot_Pythia8.root','/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_ntrkOffline_EPOS/results/trkCorr_forHBTAnalysis_tot_EPOS.root','/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_spraceMinusChi2NdofCuts_Pythia8_V2/results/trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_V2.root','/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_spraceMinusChi2NdofCuts_EPOS/results/trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_EPOS.root']

################################

config.section_('Data')###add a new section of type "Data"
###Data: This section contains all the parameters related to the data to be analyzed, 
###including the splitting parameters.
config.Data.inputDataset      = '/ReggeGribovPartonMC_13TeV-EPOS/RunIISpring15DR74-NoPURealisticRecodebug_castor_741_p1_mcRun2_Realistic_50ns_v0-v1/GEN-SIM-RECODEBUG'
#config.Data.inputDataset      = '/MinBias_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-NoPURealisticRecodebug_741_p1_mcRun2_Realistic_50ns_v0-v1/GEN-SIM-RECODEBUG'
config.Data.splitting         = 'FileBased'
#config.Data.unitsPerJob       = X ###files per job (but not impose)
#config.Data.totalUnits        = Y ###how many files to analyze
config.Data.unitsPerJob       = 6  ##EPOS
config.Data.totalUnits        = 674
#config.Data.unitsPerJob       = 6 ##Pythia8
#config.Data.totalUnits        = 1456
#config.Data.inputDBS          = 'phys03'
config.Data.inputDBS          = 'global'
#config.Data.publishDBS = 'phys03'
#config.Data.publishDataName = 'TestingJobs_RecoWithNewPixelTracks'
#config.Data.outputDatasetTag = 'HBTAnalysispp13TeV_MC_Test01' ###change for each sample(task) --only if publish

#config.Data.outLFN            = '/store/user/caber/TESTJOBS_RecoWithNewPixelTracks'
config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_MC_EPOS_MB_SPRACEcuts_EtaMix_NoTree_V10_RecoCorr' ###change for each sample(task)
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_MC_pythia8_MB_SPRACEcuts_EtaMix_NoTree_V10_RecoCorr'

###IMPORTANT:only for pythia8 that suddenly becomes INVALID in DAS.. ????
config.Data.allowNonValidInputDataset = True

################################

config.section_('Site')###add a new section of type "Site"
###Site: Grid site parameters are defined in this section, including the stage out information 
###(e.g. stage out destination site, white/black lists, etc.).
config.Site.storageSite       = 'T2_BR_SPRACE'
#config.Site.whitelist         = ['T2_US_MIT']
