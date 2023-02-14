###For a description of the crabConfig.py parameters. See:
###https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile

from CRABAPI.RawCommand import crabCommand

from dbs.apis.dbsClient import DbsApi
dbs = DbsApi('https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader')

dataset = '/HighMultiplicity_13TeV_pythia6_TuneZ2star/caber-HM95_pythia6_gen_batch2-bdb556d1414f7c9e53e113504f19a65b/USER'
fileDictList=dbs.listFiles(dataset=dataset)

print ("dataset %s has %d files" % (dataset, len(fileDictList)))

lfnList = [ dic['logical_file_name'] for dic in fileDictList ]


from WMCore.Configuration import Configuration
config = Configuration() ###create a Configuration object

config.section_('General')###add a new section of type "General"
###General: In this section, the user specifies generic parameters about the request (e.g. request name).
config.General.workArea     = 'HBTanalysispp13TeV_3D_MC_Projects' ###fixed name for projects dir in my area

#config.General.requestName  = 'hbtAnalysispp13TeV_3D_mc_EPOS_MB_SPRACEcuts_EtaMix_NoTree_V1_Gen_V2' #sub dir with prefix "crab_". Change it for each task
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_mc_pythia8_MB_SPRACEcuts_EtaMix_NoTree_V1_Gen'
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_mc_pythia6_HM95_SPRACEcuts_RandomMix_NoTree_V1_Gen_batch1'
config.General.requestName  = 'hbtAnalysispp13TeV_3D_mc_pythia6_HM95_SPRACEcuts_RandomMix_NoTree_V1_Gen_batch2'

#config.General.requestName  = 'hbtAnalysispp13TeV_3D_mc_pythia8Tune4C_HM95_SPRACEcuts_EtaMix_V1_Gen_batch1'
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_mc_pythia8TuneHMDPS_SPRACEcuts_EtaMix_V1_Gen_batch1'

#config.General.requestName  = 'hbtAnalysispp13TeV_3D_mc_pythia8Tune4C_HM95_SPRACEcuts_RandomMix_V1_Gen_batch1_V2'
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_mc_pythia6_HM95_SPRACEcuts_RandomMix_V1_Gen_batch1_V2'
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_mc_pythia6_HM95_SPRACEcuts_RandomMix_V1_Gen_batch1_V2'
#config.General.requestName  = 'hbtAnalysispp13TeV_3D_mc_pythia8Tune4C_HM95_SPRACEcuts_RandomMix_V1_Gen_batch1_V2'

#config.General.requestName  = 'hbtAnalysispp13TeV_3D_mc_EPOS_MB_SPRACEcuts_RandomMix100_NoTree_V1_Gen'


config.General.transferLogs = False
config.General.transferOutputs = True

################################

config.section_('JobType')###add a new section of type "JobType"
###JobType: This section aims to contain all the parameters of the user job type and 
###related configurables (e.g. CMSSW parameter-set configuration file, additional input files, etc.).
config.JobType.pluginName     = 'Analysis'
config.JobType.psetName       = 'ConfFile_cfg.py'

#config.JobType.allowNonProductionCMSSW = True
config.JobType.maxMemoryMB    = 2500


config.JobType.inputFiles     = ['/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/NTuplizerHBT_GenOnly/NTuplizerHBT_GenOnly/python/trkCorr_forHBTAnalysis_tot_Pythia8.root']

################################

config.section_('Data')###add a new section of type "Data"
###Data: This section contains all the parameters related to the data to be analyzed, 
###including the splitting parameters.
#config.Data.inputDataset      = '/ReggeGribovPartonMC_13TeV-EPOS/RunIISpring15DR74-NoPURealisticRecodebug_castor_741_p1_mcRun2_Realistic_50ns_v0-v1/GEN-SIM-RECODEBUG'
#config.Data.inputDataset      = '/MinBias_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-NoPURealisticRecodebug_741_p1_mcRun2_Realistic_50ns_v0-v1/GEN-SIM-RECODEBUG'
#config.Data.inputDataset      = '/MinBias_TuneCUETP8M1_13TeV_pythia8/davidlw-HM95_gen_batch1-1712b3406596258c915a8d1c953f3e82/USER'
#config.Data.inputDataset      = '/MinBias_TuneCUETP8M1_13TeV_pythia8/davidlw-HM95_gen_batch2-1712b3406596258c915a8d1c953f3e82/USER'
#config.Data.inputDataset      = '/HighMultiplicity_13TeV_pythia6_TuneZ2star/caber-HM95_pythia6_gen_batch1-bdb556d1414f7c9e53e113504f19a65b/USER'
#config.Data.inputDataset      = '/HighMultiplicity_13TeV_pythia6_TuneZ2star/caber-HM95_pythia6_gen_batch2-bdb556d1414f7c9e53e113504f19a65b/USER'
#config.Data.inputDataset      = '/HighMultiplicity_13TeV_pythia8_Tune4C/caber-HM95_pythia8_Tune4C_gen_batch1-9fcce8d38a50a9ab5d0b13cb5029d21d/USER' ###Monash
#config.Data.inputDataset      = '/HighMultiplicity_13TeV_pythia8_Tune4C_NEW/caber-HM95_pythia8_Tune4C_NEW_gen_batch1-36c36a4568dd8afd1abbbd3bea5174f8/USER'###Really 4C 
#config.Data.inputDataset      = '/HighMultiplicity_13TeV_pythia8_Tune4C_NEW_HMtune/caber-HM95_pythia8_Tune4C_NEW_gen_HMtune_batch1-375bdab058651cc9a1fc5b14783f1015/USER' ###DPS tune

config.Data.userInputFiles = lfnList

config.Data.splitting         = 'FileBased'
#config.Data.unitsPerJob       = X ###files per job (but not impose)
#config.Data.totalUnits        = Y ###how many files to analyze
#config.Data.unitsPerJob       = 6  ##EPOS
#config.Data.totalUnits        = 674
#config.Data.unitsPerJob       = 6 ##Pythia8
#config.Data.totalUnits        = 1456
#config.Data.unitsPerJob       = 50  ##HM batch1
#config.Data.totalUnits        = 4872
config.Data.unitsPerJob       = 50  ##HM batch2
##config.Data.totalUnits        = 4763
#config.Data.unitsPerJob       = 50  ##HM batch1 - pythia6
#config.Data.totalUnits        = 4850
#config.Data.unitsPerJob       = 50  ##HM batch2 - pythia6
#config.Data.totalUnits        = 4662
#config.Data.unitsPerJob       = 50  ##HM batch1 - pythia8Tune4C -> Monash
#config.Data.totalUnits        = 4993
#config.Data.unitsPerJob       = 50  ##HM batch1 - pythia8Tune4C - Really 4C tune
#config.Data.totalUnits        = 4609
#config.Data.unitsPerJob       = 50  ##HM batch1 -  HM DPS tune
#config.Data.totalUnits        = 4860
##config.Data.inputDBS          = 'phys03'##for HM
#config.Data.inputDBS          = 'global' ##for MB
#config.Data.publishDBS = 'phys03'
#config.Data.publishDataName = 'TestingJobs_RecoWithNewPixelTracks'
#config.Data.outputDatasetTag = 'HBTAnalysispp13TeV_3D_MC_Test01' ###change for each sample(task) --only if publish

#config.Data.outLFN            = '/store/user/caber/TESTJOBS_RecoWithNewPixelTracks'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_MC_EPOS_MB_SPRACEcuts_EtaMix_NoTree_V1_Gen_V2' ###change for each sample(task)
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_MC_pythia8_MB_SPRACEcuts_EtaMix_NoTree_V1_Gen'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_MC_pythia6_HM95_SPRACEcuts_RandomMix_NoTree_V1_Gen_batch1'
config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_MC_pythia6_HM95_SPRACEcuts_RandomMix_NoTree_V1_Gen_batch2'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_MC_pythia8Tune4C_HM95_SPRACEcuts_EtaMix_V1_Gen_batch1'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_MC_pythia8TuneHMDPS_SPRACEcuts_EtaMix_V1_Gen_batch1'

#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_MC_pythia8Tune4C_HM95_SPRACEcuts_RandomMix_V1_Gen_batch1_V2'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_MC_pythia6_HM95_SPRACEcuts_RandomMix_V1_Gen_batch1_V2'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_MC_pythia6_HM95_SPRACEcuts_RandomMix_V1_Gen_batch1_V2'
#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_MC_pythia8Tune4C_HM95_SPRACEcuts_RandomMix_V1_Gen_batch1_V2'

#config.Data.outLFNDirBase            = '/store/user/caber/HBTAnalysispp13TeV_3D_MC_EPOS_MB_SPRACEcuts_RandomMix100_NoTree_V1_Gen'

###IMPORTANT:only for pythia8 that suddenly becomes INVALID in DAS.. ????
config.Data.allowNonValidInputDataset = True

config.Data.publication = False

###IMPORTANT: comment this after
#config.Data.ignoreLocality = True


################################

config.section_('Site')###add a new section of type "Site"
###Site: Grid site parameters are defined in this section, including the stage out information 
###(e.g. stage out destination site, white/black lists, etc.).
config.Site.storageSite       = 'T2_BR_SPRACE'
config.Site.whitelist         = ['T2_BR_SPRACE']

###IMPORTANT: comment this after
#config.Site.blacklist = ['T2_BR_SPRACE']



result = crabCommand('submit', config = config)
print (result)



