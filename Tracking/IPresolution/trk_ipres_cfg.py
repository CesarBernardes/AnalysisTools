import FWCore.ParameterSet.Config as cms

process = cms.Process('TRACKKINEFF')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('analysis')

options.register ('n',
                  10, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, bool or float
                  "n")
options.parseArguments()


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.n)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)

### Input source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    ###fileNames =  cms.untracked.vstring("/store/relval/CMSSW_14_0_0_pre1/RelValMinBias_14TeV/MINIAODSIM/133X_mcRun3_2024_realistic_v5_STD_2024_noPU-v1/2590000/97e613f3-7482-4052-8066-4e1d608fc9c4.root"),
    fileNames =  cms.untracked.vstring(
##"/store/relval/CMSSW_14_0_0_pre1/RelValMinBias_14TeV/GEN-SIM-RECO/133X_mcRun3_2024_realistic_v5_STD_2024_noPU-v1/2590000/ed31c075-65f8-4797-a2a0-2e3929cdde7d.root",
##"/store/relval/CMSSW_14_0_0_pre1/RelValMinBias_14TeV/GEN-SIM-RECO/133X_mcRun3_2024_realistic_v5_STD_2024_noPU-v1/2590000/f9879a50-8f78-43d9-96a1-7ea92b010a1a.root",
##"/store/relval/CMSSW_14_0_0_pre1/RelValMinBias_14TeV/GEN-SIM-RECO/133X_mcRun3_2024_realistic_v5_STD_2024_noPU-v1/2590000/06425b6f-acfc-4dbf-99ba-93ad90aac397.root",
##"/store/relval/CMSSW_14_0_0_pre1/RelValMinBias_14TeV/GEN-SIM-RECO/133X_mcRun3_2024_realistic_v5_STD_2024_noPU-v1/2590000/249b0881-09dd-4e09-a2ca-8578059e2e99.root",
##"/store/relval/CMSSW_14_0_0_pre1/RelValMinBias_14TeV/GEN-SIM-RECO/133X_mcRun3_2024_realistic_v5_STD_2024_noPU-v1/2590000/477029d1-c01c-46db-b079-c873fd2503de.root",
##"/store/relval/CMSSW_14_0_0_pre1/RelValMinBias_14TeV/GEN-SIM-RECO/133X_mcRun3_2024_realistic_v5_STD_2024_noPU-v1/2590000/4a791de1-8f0d-4c8a-8524-fa2f91034e74.root",
##"/store/relval/CMSSW_14_0_0_pre1/RelValMinBias_14TeV/GEN-SIM-RECO/133X_mcRun3_2024_realistic_v5_STD_2024_noPU-v1/2590000/6e252477-a9c2-48c4-a879-c99e1cd77042.root",
##" /store/relval/CMSSW_14_0_0_pre1/RelValMinBias_14TeV/GEN-SIM-RECO/133X_mcRun3_2024_realistic_v5_STD_2024_noPU-v1/2590000/ab6b290c-04d5-45cc-91bd-a9e3c4c27fa8.root",
##"/store/relval/CMSSW_14_0_0_pre1/RelValMinBias_14TeV/GEN-SIM-RECO/133X_mcRun3_2024_realistic_v5_STD_2024_noPU-v1/2590000/adc47eb4-778f-4329-8c24-0f6aa134cdc9.root",
##"/store/relval/CMSSW_14_0_0_pre1/RelValMinBias_14TeV/GEN-SIM-RECO/133X_mcRun3_2024_realistic_v5_STD_2024_noPU-v1/2590000/af40b2c4-c855-439b-97d2-eb252a3b1551.root"
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/03e07c88-c621-4f7b-a2e9-04aa7fac50e2.root", 
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/04066584-0fb7-4f91-b269-4f02056ca30d.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/0581cd5a-3592-4528-8cd8-abb862be8ae0.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/0cae5cd5-1cea-433a-9071-58b812631a15.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/0d216794-8620-4252-8aad-a4cc4dde613e.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/1a669cb4-06fc-4fa4-8522-6c14a1c02978.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/1bafbcde-3b37-4d3c-af4b-09f25e28ff1a.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/1d53c50c-5dc4-4b31-bdb3-6856c490b818.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/1ffd5798-3f2e-4abc-bc21-2f8c64cb1234.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/20df3250-05ec-4a39-8331-58bb96a78c06.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/21777830-a28a-45b9-89fb-6eaaf927902c.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/258193c5-3102-4c60-b782-e86bbd8f7447.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/2a87cb01-146b-49a6-8598-fd902bc0e183.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/2bf7d541-6f8e-410d-be9d-9856565ae51f.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/2cf783c3-d347-4a5d-a75c-89a9fe8669d8.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/2d3bd611-bf03-43c9-bc0e-85b185c34b6b.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/2f84a548-60d3-408a-a045-369d2431087e.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/317592f2-22a0-4c60-a4e9-81e63a796795.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/3251e4c9-91bd-48f6-8fe6-58193a8104cf.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/34d81639-be16-4512-805a-3471da1ce0f9.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/35059ab1-e0a8-4fd3-901e-84aa727dd691.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/36448494-29f9-472d-92b7-36febe631ea3.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/3805cfbc-33fd-4872-9e28-1af19c6d0ba3.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/3e0112c2-7bb6-47f7-8bcf-be9533bbec4f.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/3e2180de-770b-49e3-becd-3b5a6485dc11.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/3f78213d-fe9e-49ea-a318-597d672dced2.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/3fb95f3e-9008-4be5-97db-a803fe28f8ee.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/4565bdfa-c752-47d2-8f71-b38ee672a98f.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/463553fa-3061-497f-80c1-e3769fea23c9.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/47a4cbd1-ea02-4e51-9418-0df61d0dead5.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/4a22ec5b-0a5f-415a-86a3-24bd1eb1ad02.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/4a42c172-0c22-44a4-b9b9-ff37977c87b5.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/5009af34-cfab-4e49-83fc-c801a3b61a24.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/503de3b9-8af6-46c0-b83a-e1a1b90eb17b.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/52a555bd-8ac1-40dd-912c-f8c9522b0d82.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/597613c2-7e42-4057-ba2b-948c979165e9.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/5aa92fa7-7891-42bf-bf6f-ba314a0a4cf2.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/5d84286a-574d-4991-90ca-3bd42a08b148.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/5e1c926d-dfdb-4751-a917-c9f06b429a98.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/6215843d-5edf-4396-9c46-355468d221a5.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/66818a77-ead7-40aa-b1a5-e66e2da7f410.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/6827a810-f181-4f23-be9c-8825119b1010.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/6a54e1e6-6e08-4895-861e-5e923117e71c.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/6b3c023a-69bd-4732-97b0-44af9fe9f1e5.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/6b75d4a7-5058-45f3-bcd2-f82eac50acca.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/6bbc5d37-4e53-4f7e-a513-7daf7b3722b0.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/6e994ade-1c13-4785-a669-3a845b324907.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/71fbd0a9-ad4a-416a-8d39-77ac1bd8910f.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/75425330-7919-4a45-8273-4ea07b7dd600.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/81e31d3d-9a87-4e9b-a93e-4e9ccb687b94.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/8ce1d573-3f2d-4f0f-a62f-964b1f245f9a.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/8f6bb72f-db13-4483-a0c4-3d42b23d3837.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/8fbd72b9-0336-43f3-a70a-2517e39a0163.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/90c22c83-4cec-4f32-a1c2-d3ccd3ce9247.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/9227a0d2-0c12-49ce-9677-673080d97ee2.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/928ed43c-bab6-4062-9946-da8fbd8e5ca7.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/94460266-9bcc-49e6-b749-15a25eb946cd.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/94611f6d-e7ff-437b-b5c2-7dfc0de40cde.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/9b1e89f8-a3ba-4549-bca6-0b6976dee758.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/9b977459-ef48-4b65-b1ea-6ddefa0d6d06.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/9f637b4b-4d5f-42d9-94fe-f984093f5eaf.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/9fd11550-1cf0-4e82-b352-2de25c03a0df.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/a1f9b755-a688-471b-bc33-4472a568a6a3.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/a3146d3d-6ead-4b3d-8bbb-0e5dd7fd8172.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/a40fe1b4-e5b0-4fbe-99db-49af4a944da0.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/adec9457-70b3-4270-bf8f-aa4dbe5730c8.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/af0d84af-c355-4cc1-ba76-6ea1b0508979.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/b142ae55-aed7-4b39-a11f-ca541a4e0c67.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/b1e1c997-bbad-4cea-b68b-ecb3336941ab.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/b25dab97-99bf-495a-8e29-29d2127048f2.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/b34fce17-f30b-4216-a27f-91b124e33f13.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/b44f7e1f-5375-4597-9bdd-78f8795af9a7.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/b8845663-552b-4ea9-96f5-aa3c7441e8f3.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/bd032ebb-d6ce-4565-b340-10fbdca5b93f.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/bfc787d3-d6bc-43d6-8192-32eb17685baf.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/c68aacbc-bd7c-4159-905c-1d2ad2173018.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/cbcb29fd-5bf1-42cc-9800-246252555389.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/cdcf5c0d-f950-463f-a853-97ee6c856680.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/d67ec9be-391b-41c2-99cb-5bd0881d7bc5.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/d8e583cf-d383-47fd-a6e3-b4716f64e8d8.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/d9943273-0225-40ee-ae1c-5c47963faf06.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/daa665d1-85b8-4175-bb14-3e11aabac8e7.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/daccb06d-5ef7-4069-8308-0467c1b697d7.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/db1e3d8b-4efd-4ffd-aa86-526b41006b6f.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/e74327b4-9b6d-4563-aa15-0c95cb1e7dc1.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/e946ee97-ccc0-4c47-bc3f-c375535ab362.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/ed3695bd-e3c5-4538-89aa-d9e289c2745a.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/f08c3478-7035-4558-8dc6-c256d15789b2.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/f0fca09e-b670-4cfd-93df-16438c43c3fb.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/f2704a58-b274-4d3b-a8fe-383511358a8b.root",
"file:/eos/cms/store/relval/CMSSW_14_0_0/RelValMinBias_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v1_STD_2026D98_noPU-v1/2580000/f297f61d-d455-4098-82ff-76dd891eb9a6.root"
       ),
    skipEvents = cms.untracked.uint32(0),
    secondaryFileNames = cms.untracked.vstring()
)

### output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('trk_ipres.root')
)

###Global tag
from Configuration.AlCa.GlobalTag import GlobalTag
####process.GlobalTag = GlobalTag(process.GlobalTag, '133X_mcRun3_2024_realistic_v5', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_mcRun4_realistic_v1', '')

###EDAnalyzer for ip resolution calculation
process.trk_ipres = cms.EDAnalyzer("TrackingIPres",
        vertices  = cms.InputTag("offlinePrimaryVertices"),
        tracks   = cms.InputTag("generalTracks"),
        genparticles = cms.InputTag("genParticles")
        #vertices  = cms.InputTag("offlineSlimmedPrimaryVertices"),
        #tracks   = cms.InputTag("packedPFCandidates"),
        #genparticles = cms.InputTag("packedGenParticles")
)

process.p = cms.Path(process.trk_ipres)
