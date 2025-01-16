#Loading necessary libraries
import FWCore.ParameterSet.Config as cms
###from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
import PhysicsTools.PythonAnalysis.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
process = cms.Process('HiForest')
process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

#Number of events: put '-1' unless testing
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#HiForest script init
process.load("HiForest_cff")
process.HiForest.inputLines = cms.vstring("HiForest V3",)
version = 'no git info'
process.HiForest.HiForestVersion = cms.string(version)

goodJSON = 'Cert_181530-183126_HI7TeV_PromptReco_Collisions11_JSON.txt'
myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')
import FWCore.Utilities.FileUtils as FileUtils
files2011data = FileUtils.loadListFromFile ('CMS_HIRun2011_HIMinBiasUPC_SPRACE_part2.txt')
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*files2011data)
###    fileNames = cms.untracked.vstring("file:file2011PbPb.root")
#######    fileNames = cms.untracked.vstring(
#######"root://cmsxrootd.fnal.gov///store/data/HIRun2011/HIMinBiasUPC/RECO/22Feb2012-v2/10000/A0A5FD9A-937D-E111-AF9B-842B2B688E9B.root",
#######"root://cmsxrootd.fnal.gov///store/data/HIRun2011/HIMinBiasUPC/RECO/22Feb2012-v1/0013/BC7CB2E8-495E-E111-A087-003048FFCBA8.root",
#######"root://cmsxrootd.fnal.gov///store/data/HIRun2011/HIMinBiasUPC/RECO/22Feb2012-v1/0000/06A74613-445E-E111-97DA-003048678FE6.root"
#######    )
)
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)

#Global Tag: change the name according to the instructions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
##process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/GR_R_44_V15.db')
process.GlobalTag.globaltag = 'GR_R_44_V15::All'
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
###process.load("Configuration.StandardSequences.MagneticField_cff")
process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag

# Common stuff if you haven't already loaded it
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

# Common offline event selection
process.load("HeavyIonsAnalysis.Configuration.collisionEventSelection_cff")
###process.primaryVertexFilter.src = "offlinePrimaryVertices"

#Define the output root file (change each run not to overwrite previous output)
process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string("HiForestAOD_DATA2011_UCC_part2.root"))

#Init Trigger Analyzer
process.hltanalysis = cms.EDAnalyzer('TriggerInfoAnalyzer',
                              processName = cms.string("HLT"),
                              triggerName = cms.string("@"),         
                              datasetName = cms.string("HIMinBiasUPC"),  #'HICorePhysics' to look at Core Physics only
                              ###datasetName = cms.string("MinimumBias"),       
                              triggerResults = cms.InputTag("TriggerResults","","HLT"),
                              triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")                             
                              )

##Select triggered events here
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltUCC = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
###process.hltUCC.HLTPaths = ["HLT_HIUCC015_v*","HLT_HIUCC010_v*","HLT_HIMinBiasHfOrBSC_v*"] ###UCC or MB
process.hltUCC.HLTPaths = ["HLT_HIMinBiasHfOrBSC_v*"] ###MB
process.hltUCC.andOr = cms.bool(True)  # True = OR, False = AND between the HLT paths
process.hltUCC.throw = cms.bool(False) # throw exception on unknown path names



#Collect event data
process.demo = cms.EDAnalyzer('AnalyzerHBT',
                              processName = cms.string("HLT"),
                              triggerName = cms.string("@"),
                              datasetName = cms.string("HIMinBiasUPC"),  #'HICorePhysics' to look at Core Physics only
                              ###datasetName = cms.string("MinimumBias"),
                              triggerResults = cms.InputTag("TriggerResults","","HLT"),
                              triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")
                             ) 
process.dump=cms.EDAnalyzer('EventContentAnalyzer') #easy check of Event structure and names without using the TBrowser

###process.ana_step = cms.Path(process.hltUCC*process.collisionEventSelection*process.demo)
process.ana_step = cms.Path(process.collisionEventSelection*process.demo)

###process.ana_step = cms.Path(process.hltanalysis*
###                            #process.hltUCC*
###		  	    #process.dump+  #uncomment if necessary to check the name. Do not forget to change the number of events to '1'
###			    process.demo*
###                            process.HiForest 
###)

