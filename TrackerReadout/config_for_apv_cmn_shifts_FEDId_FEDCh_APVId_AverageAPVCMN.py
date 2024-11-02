# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --eventcontent RAW --datatier RAW --conditions auto:run2_hlt_hi --step RAW2DIGI,DIGI2RAW --scenario HeavyIons --data --era Run2_HI --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

from Configuration.Eras.Era_Run3_pp_on_PbPb_approxSiStripClusters_2024_cff import Run3_pp_on_PbPb_approxSiStripClusters_2024
process = cms.Process('DIGI2RAW',Run3_pp_on_PbPb_approxSiStripClusters_2024)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_DataMapper_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
###    input = cms.untracked.int32(5000)
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
###process.source = cms.Source("NewEventStreamFileReader",
    fileNames = cms.untracked.vstring(
            '/store/hidata/HIRun2015/HITrackerVirginRaw/RAW/v1/000/263/400/00000/40322926-4AA3-E511-95F7-02163E0146A8.root'
################'root://osg-se.sprace.org.br:1094//store/hidata/HIRun2023A/HITrackerNZS/RAW/v1/000/374/354/00000/6b6f4b36-6036-44bf-b662-bf745a743dcf.root'
    ),
###    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition


# Additional output definition
## Offline Silicon Tracker Zero Suppression
process.load('RecoLocalTracker/SiStripZeroSuppression/SiStripZeroSuppression_cfi')
process.siStripZeroSuppression.produceRawDigis = False
process.siStripZeroSuppression.produceHybridFormat = True
#####process.siStripZeroSuppression.produceHybridFormat = False
process.siStripZeroSuppression.Algorithms.APVInspectMode = 'HybridEmulation'
######process.siStripZeroSuppression.Algorithms.APVInspectMode = 'BaselineFollower'
process.siStripZeroSuppression.Algorithms.APVRestoreMode = ''
process.siStripZeroSuppression.Algorithms.CommonModeNoiseSubtractionMode = 'Median'
process.siStripZeroSuppression.Algorithms.MeanCM = 0
process.siStripZeroSuppression.Algorithms.DeltaCMThreshold = 20
process.siStripZeroSuppression.Algorithms.Use10bitsTruncation = True
process.siStripZeroSuppression.RawDigiProducersList = cms.VInputTag(cms.InputTag("siStripDigis","VirginRaw"))

process.siStripRepackZS = cms.EDProducer("SiStripDigiToRawModule",
    InputDigis       = cms.InputTag("siStripZeroSuppression", "VirginRaw"),
    FedReadoutMode = cms.string('ZERO_SUPPRESSED'),
    PacketCode = cms.string("ZERO_SUPPRESSED10"),
    UseFedKey = cms.bool(False),
    UseWrongDigiType = cms.bool(False),
    CopyBufferHeader = cms.bool(True),
    RawDataTag = cms.InputTag('rawDataCollector')
)

process.hybridRawDataRepacker = cms.EDProducer("RawDataCollectorByLabel",
    verbose = cms.untracked.int32(0),     # 0 = quiet, 1 = collection list, 2 = FED list
    RawCollectionList = cms.VInputTag( 
                                       cms.InputTag('siStripRepackZS'),
                                       cms.InputTag('rawDataCollector')
    ),
)


#this is just to check if is using the Hybrid Format properly
process.siStripUnpackRepackedZS = process.siStripDigis.clone(
        ProductLabel     = cms.InputTag('siStripRepackZS'),
)

process.TFileService = cms.Service("TFileService",
        fileName=cms.string("output_CalculateCMNonTheFly.root")
)

###for this module below should pass only the type : edm::DetSetVector<SiStripDigi>
process.hybridAna = cms.EDAnalyzer("SiStripHybridFormatAnalyzer",
    ###srcDigis =  cms.InputTag('siStripZeroSuppression','VirginRaw'), 
    srcDigis =  cms.InputTag("siStripUnpackRepackedZS","ZeroSuppressed"), 
    srcAPVCM =  cms.InputTag('siStripZeroSuppression','APVCMVirginRaw'),
    nModuletoDisplay = cms.uint32(10000),
    plotAPVCM   = cms.bool(True), ###use "False" when printing CMN calculated from pedestals on-the-fly
)


# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_dataRun3_Express_v3', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.siStripDigis)
process.DigiToRaw2 = cms.Sequence(process.siStripZeroSuppression+process.siStripRepackZS+process.hybridRawDataRepacker+process.siStripUnpackRepackedZS+process.hybridAna)
process.digi2raw_step = cms.Path(process.DigiToRaw2)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.digi2raw_step,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
