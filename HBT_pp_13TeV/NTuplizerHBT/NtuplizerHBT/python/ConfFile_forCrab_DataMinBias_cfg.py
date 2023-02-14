import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("Demo")

process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')

#####################################################################################################################
##########  Run over data or mc? Which Selection? 
#####################################################################################################################

run_over_data = True     
#run_over_data = False ###Run on MC samples

#isHighMultiplicitySample = True ###Only use "True" for DATA HM skim 
isHighMultiplicitySample = False
multiplicity = 0 ###MinBias 
#multiplicity = 1 ###use only: 80=< ntrkoffline < 105
#multiplicity = 2 ###use only: 105=< ntrkoffline < 130
#multiplicity = 3 ###use only: 130=< ntrkoffline
if (isHighMultiplicitySample == False) :
   multiplicity = 0 ###MinBias Samples

##which mixing procedure???
#mix_procedure = 1 ###random
mix_procedure = 2 ###eta-mix

#tk_selection = "Baseline" 
#tk_selection = "HBT_Padova"     
#tk_selection = "Ridge"     
tk_selection = "HBT_Sprace"

fillTree = False
#fillTree = True ###WARNING: if running in all data, we expect trees of order of TB for all types of selection!!!

#####################################################################################################################
##########  Services, conditions, ...
#####################################################################################################################
process.load('Configuration/StandardSequences/Services_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

out_file_name = 'mcNewSkim_TEST.root'
if run_over_data :
   out_file_name = 'dataNewSkim_TEST.root'
process.TFileService = cms.Service("TFileService",fileName = cms.string(out_file_name))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
if run_over_data :
   #process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
   process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v4', '')

#####################################################################################################################


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )

#mylist = FileUtils.loadListFromFile('list_of_files/TuneCUETP8M1_new_bs_1.list')
##mylist = FileUtils.loadListFromFile('list_of_files/ReggeGribovPartonMC_13TeV_EPOS.list')
#if run_over_data :
#   mylist = FileUtils.loadListFromFile('list_of_files/test_MinBias1_TOTEM_ToCompareWithPadovaCode.list')
#   #mylist = FileUtils.loadListFromFile('list_of_files/dataHMfiles.list')
#readFiles = cms.untracked.vstring(*mylist)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:October_TOTEM_PU0p1_pPb_HM_10.root'
        ######
        #'/store/user/davidlw/L1MinimumBiasHF1/RecoSkim2015_2015DLowPU_ReTracking_v4/151109_223122/0000/pPb_HM_1.root',
        #'/store/user/davidlw/L1MinimumBiasHF1/RecoSkim2015_2015DLowPU_ReTracking_v4/151109_223122/0000/pPb_HM_2.root'   
        ######
        #'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM-RECODEBUG/NoPURealisticRecodebug_741_p1_mcRun2_Realistic_50ns_v0-v1/00000/0266184C-3053-E511-9582-0025905B85E8.root'
        ######
        #'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/ReggeGribovPartonMC_13TeV-EPOS/GEN-SIM-RECODEBUG/NoPURealisticRecodebug_castor_741_p1_mcRun2_Realistic_50ns_v0-v1/00000/00AA1827-FF52-E511-80E5-1CC1DE1D0AD4.root'
        ######  
        ###'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/MinBias_TuneCUETHS1_13TeV-herwigpp/GEN-SIM-RECODEBUG/NoPURealisticRecodebug_castor_741_p1_mcRun2_Realistic_50ns_v0-v1/00000/0020F55E-D452-E511-B305-000F5327349C.root'
        ###'root://xrootd.unl.edu//store/mc/RunIISpring15DR74/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM-DIGI-RECO/NoPURecodebug_TrackingParticle_MCRUN2_74_V8B-v2/30000/04BB7B1D-1B2D-E511-91C5-001E67A4069F.root'
        ###'file:/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CreateMCSkimsForHighMultiplicity/CMSSW_7_4_15_patch1/src/MCSkimsForHighMultiplicity/MCSkimForHM/python/pp13TeV_HBTanalysis_HM.root'
    ),
    #fileNames = cms.untracked.vstring(*mylist),
    secondaryFileNames=cms.untracked.vstring( 
    #    '/store/data/Run2015D/L1MinimumBiasHF1/RECO/PromptReco-v4/000/259/152/00000/08A9BC20-9275-E511-8EAF-02163E01436D.root',
    #    '/store/data/Run2015D/L1MinimumBiasHF1/RECO/PromptReco-v4/000/259/152/00000/E43801EE-9175-E511-9FB2-02163E014606.root',
    )
)

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 15 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True),   # otherwise it won't filter the events
)

process.PAprimaryVertexFilterOnSkimData = process.PAprimaryVertexFilter.clone(cut = cms.string("abs(z) <= 15"))

#Reject beam scraping events standard pp configuration
process.NoScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

process.PAcollisionEventSelection = cms.Sequence(process.hfCoincFilter *
                                                 process.PAprimaryVertexFilter *
                                                 process.NoScraping
)
###no need to apply all the cuts in data since this is already skim
if run_over_data : 
   process.PAcollisionEventSelection = cms.Sequence(process.PAprimaryVertexFilterOnSkimData)   


if isHighMultiplicitySample :
   import HLTrigger.HLTfilters.hltHighLevel_cfi
   process.hltHighMultiplicity = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
   process.hltHighMultiplicity.HLTPaths = ["HLT_L1MinimumBiasHF1OR_part*_v*"]
   if multiplicity==1 :
      process.hltHighMultiplicity.HLTPaths = ["HLT_PixelTracks_Multiplicity60_v*"]
   if multiplicity==2 :
      process.hltHighMultiplicity.HLTPaths = ["HLT_PixelTracks_Multiplicity85_v*","HLT_PixelTracks_Multiplicity60_v*"]
   if multiplicity==3 :
      process.hltHighMultiplicity.HLTPaths = ["HLT_PixelTracks_Multiplicity110_v*","HLT_PixelTracks_Multiplicity85_v*","HLT_PixelTracks_Multiplicity60_v*"] 
   process.hltHighMultiplicity.andOr = cms.bool(True)  # True = OR, False = AND between the HLT paths
   process.hltHighMultiplicity.throw = cms.bool(False) # throw exception on unknown path names


process.demo = cms.EDAnalyzer('NtuplizerHBT',
               generalTrackTag   = cms.InputTag("generalTracks"),
               offlinePrimaryVerticesTag = cms.InputTag("offlinePrimaryVertices"), 
               tpSrc = cms.InputTag('mix','MergedTrackTruth'),
               RealData          = cms.bool(run_over_data),
               Selection         = cms.string(tk_selection),
               FillTree          = cms.bool(fillTree),
               Multiplicity      = cms.uint32(multiplicity),
               Mix_procedure     = cms.uint32(mix_procedure)
)


#process.p = cms.Path(process.demo)
process.p = cms.Path(process.PAcollisionEventSelection*process.demo)
if isHighMultiplicitySample :
   process.p = cms.Path(process.hltHighMultiplicity*process.PAcollisionEventSelection*process.demo)   
# Schedule definition
process.schedule = cms.Schedule(process.p)
