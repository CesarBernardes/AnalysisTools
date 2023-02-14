import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")


###choose if will run on MB or HM samples
#multiplicity = 0 ###MinBias --> ntrkoff in [0,79]
multiplicity = 1 ###High-multiplicity --> ntrkoff in [80,250]

##which mixing procedure???
mix_procedure = 1 ###random
#mix_procedure = 2 ###eta-mix


process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration/StandardSequences/Services_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.TFileService = cms.Service("TFileService",fileName = cms.string('mcGenOnly.root'))


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/MinBias_13TeV_pythia8_TuneCUETP8M1_cfi_py_GEN.root'
        #####'file:/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/MinBias_13TeV_pythia6_TuneZ2star_cfi_py_GEN.root' 
        ####'/store/mc/RunIIFall15DR76/ChMulti85_ReggeGribov_PartonMC_13TeV-EPOS/GEN-SIM-RECODEBUG/PU25nsData2015v1FSQ_76X_mcRun2_asymptotic_v12-v1/00000/00F85EFD-C3DC-E511-B7E6-A0369F7FC070.root'
        #'file:/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/MinBias_13TeV_pythia8_Tune4C_cfi_py_GEN.root'
        #'file:/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/MinBias_13TeV_pythia8_pdfCTEQ6L1_cfi_py_GEN.root'  
        #'/store/mc/RunIISpring15DR74/MinBias_TuneCUETP8M1_13TeV-pythia8/GEN-SIM-RECODEBUG/NoPURealisticRecodebug_741_p1_mcRun2_Realistic_50ns_v0-v1/00000/0266184C-3053-E511-9582-0025905B85E8.root'
      #'/store/user/davidlw/MinBias_TuneCUETP8M1_13TeV_pythia8/HM95_gen_batch1/160717_143730/0000/MinBias_13TeV_pythia8_TuneCUETP8M1_cfi_py_GEN_1.root',
      #'/store/user/davidlw/MinBias_TuneCUETP8M1_13TeV_pythia8/HM95_gen_batch1/160717_143730/0000/MinBias_13TeV_pythia8_TuneCUETP8M1_cfi_py_GEN_10.root',
      #'/store/user/davidlw/MinBias_TuneCUETP8M1_13TeV_pythia8/HM95_gen_batch1/160717_143730/0000/MinBias_13TeV_pythia8_TuneCUETP8M1_cfi_py_GEN_100.root',
      #'/store/user/davidlw/MinBias_TuneCUETP8M1_13TeV_pythia8/HM95_gen_batch1/160717_143730/0000/MinBias_13TeV_pythia8_TuneCUETP8M1_cfi_py_GEN_101.root',
      #'/store/user/davidlw/MinBias_TuneCUETP8M1_13TeV_pythia8/HM95_gen_batch1/160717_143730/0000/MinBias_13TeV_pythia8_TuneCUETP8M1_cfi_py_GEN_102.root',
      #'/store/user/davidlw/MinBias_TuneCUETP8M1_13TeV_pythia8/HM95_gen_batch1/160717_143730/0000/MinBias_13TeV_pythia8_TuneCUETP8M1_cfi_py_GEN_103.root',
      #'/store/user/davidlw/MinBias_TuneCUETP8M1_13TeV_pythia8/HM95_gen_batch1/160717_143730/0000/MinBias_13TeV_pythia8_TuneCUETP8M1_cfi_py_GEN_104.root',
      #'/store/user/davidlw/MinBias_TuneCUETP8M1_13TeV_pythia8/HM95_gen_batch1/160717_143730/0000/MinBias_13TeV_pythia8_TuneCUETP8M1_cfi_py_GEN_105.root',
      #'/store/user/davidlw/MinBias_TuneCUETP8M1_13TeV_pythia8/HM95_gen_batch1/160717_143730/0000/MinBias_13TeV_pythia8_TuneCUETP8M1_cfi_py_GEN_106.root',
    )
)

process.demo = cms.EDAnalyzer('NTuplizerHBT_GenOnly',
   genPTag           = cms.InputTag("genParticles"),
   Multiplicity      = cms.uint32(multiplicity),
   Mix_procedure     = cms.uint32(mix_procedure)
)


process.p = cms.Path(process.demo)
