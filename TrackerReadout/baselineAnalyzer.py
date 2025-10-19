import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras


from Configuration.Eras.Era_Run3_pp_on_PbPb_2023_cff import Run3_pp_on_PbPb_2023
process = cms.Process('DIGI2RAW',Run3_pp_on_PbPb_2023)

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
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
            'file:input.root'
    ),
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
process.TFileService = cms.Service("TFileService",
        fileName=cms.string("baselines.root")
)

###for this module below should pass only the type : edm::DetSetVector<SiStripDigi>
process.hybridAna = cms.EDAnalyzer("SiStripHybridFormatAnalyzer",
    ###srcDigis =  cms.InputTag('siStripZeroSuppression','VirginRaw'), 
    srcDigis =  cms.InputTag("siStripUnpackRepackedZS","ZeroSuppressed"),
    srcAPVCM =  cms.InputTag('siStripZeroSuppression','APVCMVirginRaw'),
    nModuletoDisplay = cms.uint32(10000),
    plotAPVCM   = cms.bool(True)
)


# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Express_v4', '')


# Path and EndPath definitions
process.DigiToRaw2 = cms.Sequence(process.hybridAna)
process.digi2raw_step = cms.Path(process.DigiToRaw2)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.digi2raw_step,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
