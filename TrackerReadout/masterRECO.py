# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: myMasterRecoFile --conditions 103X_dataRun2_Express_v1 -s RAW2DIGI,L1Reco,RECO --data --era Run2_2018_pp_on_AA --eventcontent AOD --runUnscheduled --datatier RECO --repacked --no_exec --nThreads=8
import FWCore.ParameterSet.Config as cms

###from Configuration.StandardSequences.Eras import eras

###process = cms.Process('RECO',eras.Run2_2018_pp_on_AA)

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
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

####process.source = cms.Source("PoolSource",
####    fileNames = cms.untracked.vstring(
process.source = cms.Source("NewEventStreamFileReader",
    fileNames = cms.untracked.vstring(
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0051_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0052_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0053_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0054_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0056_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0055_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0057_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0058_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0059_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0060_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0061_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0062_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0063_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0064_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0065_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0066_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0067_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0068_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0069_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0070_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0071_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0122_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0123_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0124_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0125_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0126_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0127_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0128_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0129_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0130_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0131_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0132_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0133_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0134_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0135_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0136_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0137_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0138_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0139_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0140_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0142_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/288/run374288_ls0141_streamPhysicsHITrackerNZS_StorageManager.dat',
####
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0070_streamPhysicsHITrackerNZS_StorageManager.dat', 
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0071_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0083_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0084_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0096_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0097_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0110_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0111_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0129_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0128_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0140_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0141_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0157_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0158_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0175_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0176_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0197_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0196_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0216_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0217_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0233_streamPhysicsHITrackerNZS_StorageManager.dat', 
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0235_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0255_streamPhysicsHITrackerNZS_StorageManager.dat',
######'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/289/run374289_ls0257_streamPhysicsHITrackerNZS_StorageManager.dat',
######
'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0052_streamPhysicsHITrackerNZS_StorageManager.dat', 
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0050_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0051_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0053_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0054_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0055_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0057_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0056_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0058_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0059_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0060_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0062_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0063_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0061_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0064_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0065_streamPhysicsHITrackerNZS_StorageManager.dat',
'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0066_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0067_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0068_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0069_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0070_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0071_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0072_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0073_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0076_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0079_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0080_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0082_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0083_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0084_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0085_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0086_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0087_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0088_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0089_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0090_streamPhysicsHITrackerNZS_StorageManager.dat',
'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0091_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0092_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0093_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0094_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0095_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0097_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0096_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0100_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0101_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0078_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0081_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0103_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0075_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0074_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0106_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0077_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0107_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0108_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0102_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0109_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0111_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0110_streamPhysicsHITrackerNZS_StorageManager.dat',
'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0112_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0113_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0114_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0115_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0116_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0117_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0118_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0120_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0121_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0119_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0098_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0099_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0104_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0105_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0123_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0122_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0124_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0126_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0125_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0127_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0128_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0129_streamPhysicsHITrackerNZS_StorageManager.dat',
'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0130_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0131_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0132_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0133_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0134_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0135_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0136_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0137_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0138_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0139_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0141_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0140_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0142_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0143_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0144_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0145_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0146_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0147_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0148_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0149_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0150_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0151_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0153_streamPhysicsHITrackerNZS_StorageManager.dat',
'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0152_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0154_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0156_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0155_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0158_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0157_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0159_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0160_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0161_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0162_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0163_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0164_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0165_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0166_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0167_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0168_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0169_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0170_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0171_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0173_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0172_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0174_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0175_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0176_streamPhysicsHITrackerNZS_StorageManager.dat',
'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0177_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0178_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0180_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0179_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0182_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0181_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0183_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0185_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0184_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0186_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0187_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0188_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0190_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0189_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0191_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0192_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0193_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0194_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0196_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0195_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0197_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0198_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0199_streamPhysicsHITrackerNZS_StorageManager.dat',
'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0201_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0200_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0202_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0203_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0204_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0205_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0206_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0208_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0207_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0209_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0210_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0211_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0212_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0213_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0215_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0214_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0216_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0217_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0219_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0218_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0220_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0221_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0222_streamPhysicsHITrackerNZS_StorageManager.dat',
'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0223_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0224_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0226_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0225_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0227_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0228_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0229_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0230_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0231_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0232_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0233_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0234_streamPhysicsHITrackerNZS_StorageManager.dat',
####'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0236_streamPhysicsHITrackerNZS_StorageManager.dat',
'file:/eos/cms/store/t0streamer/Data/PhysicsHITrackerNZS/000/374/307/run374307_ls0235_streamPhysicsHITrackerNZS_StorageManager.dat',
),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('myMasterRecoFile nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RECO'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string('out.root'),
    outputCommands = process.AODEventContent.outputCommands
)

process.siStripDigis.UnpackCommonModeValues = cms.bool(True)
process.AODEventContent.outputCommands.extend(['keep *_siStripClusters_*_*','keep *_siStripDigis_CommonMode_*'])

# Additional output definition
process.TFileService = cms.Service("TFileService",
      fileName = cms.string("file_cmn_.root"),
      closeFileFast = cms.untracked.bool(True)
)

process.siStripCMNanalyzer = cms.EDAnalyzer("SiStripCMNanalyzer",
        srcAPVCM =  cms.InputTag("siStripDigis","CommonMode")
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
###############process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_Express_v2', '')
##process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
##process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_hlt_hi', '')
##process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_hlt_GRun', '')
###process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_Prompt_v3', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Express_v4', '')

###process.siStripClusters.Clusterizer.RemoveApvShots = cms.bool(False)

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi*process.siStripCMNanalyzer)
###process.raw2digi_step = cms.Path(process.siStripDigis)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODoutput_step = cms.EndPath(process.AODoutput)


# Schedule definition
###process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.AODoutput_step)
process.schedule = cms.Schedule(process.raw2digi_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(1)
process.options.numberOfStreams=cms.untracked.uint32(1)

##############from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
##############MassReplaceInputTag(process, new="rawDataMapperByLabel", old="rawDataCollector")
##############MassReplaceInputTag(process, new="rawDataRepacker", old="rawDataCollector")
##############MassReplaceInputTag(process, new="hltFEDSelector", old="rawDataCollector")

#do not add changes to your config after this point (unless you know what you are doing)
#from FWCore.ParameterSet.Utilities import convertToUnscheduled
#process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
