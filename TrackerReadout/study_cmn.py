# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --eventcontent RAW --datatier RAW --conditions auto:run2_hlt_hi --step RAW2DIGI,DIGI2RAW --scenario HeavyIons --data --era Run2_HI --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('DIGI2RAW',eras.Run2_2018_pp_on_AA)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

##process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1))

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
###      'file:/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/StreamerForests/HIPhysicsMinimumBias0/AOD/v0/000/326/381/step3_3.root',    
###      'file:/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/StreamerForests/HIPhysicsMinimumBias0/AOD/v0/000/326/381/step3_13.root',
###      'file:/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/StreamerForests/HIPhysicsMinimumBias0/AOD/v0/000/326/381/step3_16.root',
###      'file:/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/StreamerForests/HIPhysicsMinimumBias0/AOD/v0/000/326/381/step3_21.root',
###      'file:/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/StreamerForests/HIPhysicsMinimumBias0/AOD/v0/000/326/381/step3_4.root',
###      'file:/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/StreamerForests/HIPhysicsMinimumBias0/AOD/v0/000/326/381/step3_5.root'
##       'file:/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/StreamerForests/HIPhysicsMinimumBias0/AOD/000/326/389/step3_20.root',
##       'file:/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/StreamerForests/HIPhysicsMinimumBias0/AOD/000/326/389/step3_27.root',
##       'file:/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/StreamerForests/HIPhysicsMinimumBias0/AOD/000/326/389/step3_4.root'
#      'file:/eos/cms/store/group/phys_heavyions/bdiab/PbPb2018/Forests/StreamerForests/HIPhysicsMinimumBias0/AOD/000/326/392/step3_29.root',
###      'file:out_326483_600events_NoAPVshotCleaning.root'
      'file:out_326483_600events_Default.root'
    ),
    ##eventsToSkip = cms.untracked.VEventRange('323416:111297492','323416:111120231','323416:110857784')
    ##skipEvents=cms.untracked.uint32(0)
) 

# Output definition

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')



process.TFileService = cms.Service("TFileService",
        fileName=cms.string("fed_cmn.root")
)

process.siStripCMNanalyzer = cms.EDAnalyzer("SiStripCMNanalyzer",
        srcAPVCM =  cms.InputTag("siStripDigis","CommonMode")
)

process.clusterAnaZS = cms.EDAnalyzer("SiStripClusterAnalyzer",
        srcClusters =  cms.InputTag('siStripClusters','')
)

process.trackingAnaZS = cms.EDAnalyzer("TrackingAnalyzer",
        srcTracks =  cms.InputTag('generalTracks')
)


# Path and EndPath definitions
###process.raw2digi_step = cms.Path(process.siStripCMNanalyzer*process.clusterAnaZS)
###process.raw2digi_step = cms.Path(process.clusterAnaZS*process.trackingAnaZS)
process.raw2digi_step = cms.Path(process.siStripCMNanalyzer)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step)
