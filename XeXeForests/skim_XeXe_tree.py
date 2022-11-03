###Run in lxplus8 and using python3
### <<ssh -XY user@lxplus8.cern.ch>>
### <<voms-proxy-init -voms cms>>
### <<python3 skim_XeXe_tree.py>>

import ROOT

names = ROOT.std.vector('string')()
listOfFiles = open('list_of_files.txt', 'r')
Lines = listOfFiles.readlines()
for line in Lines:
    names.push_back(line)


#fileName = "root://osg-se.sprace.org.br//store/user/caber/HIMinimumBias1/HeavyIon_Forest_XeXe_5p44TeV_MB_PD_1_out/221019_004207/0000/HiForestAOD_1.root"
treeName1 = "hiEvtAnalyzer/HiTree"
treeName2 = "skimanalysis/HltTree"
treeName3 = "ppTrack/trackTree"


outFileName = "tree_XeXe_MB_data_skim_HBT.root"

getVector_code ='''

std::vector<float> getVectorFloat_Aux(int Ntrks, ROOT::VecOps::RVec<Float_t> vec_trkDxy1, ROOT::VecOps::RVec<Float_t> vec_trkDxyError1, ROOT::VecOps::RVec<Float_t> vec_trkDz1, ROOT::VecOps::RVec<Float_t> vec_trkDzError1, ROOT::VecOps::RVec<Float_t> vec_trkPt, ROOT::VecOps::RVec<Float_t> vec_trkPtError, ROOT::VecOps::RVec<Float_t> vec_ToReturn){

   std::vector<float> v;
   for(int i=0; i<Ntrks; i++){
      if(std::abs(vec_trkDxy1[i]/vec_trkDxyError1[i])<5 && std::abs(vec_trkDz1[i]/vec_trkDzError1[i])<5 && vec_trkPtError[i]/vec_trkPt[i]<0.1)
         v.push_back(vec_ToReturn[i]);
   }
   return v;
}

std::vector<int> getVectorInt_Aux(int Ntrks, ROOT::VecOps::RVec<Float_t> vec_trkDxy1, ROOT::VecOps::RVec<Float_t> vec_trkDxyError1, ROOT::VecOps::RVec<Float_t> vec_trkDz1, ROOT::VecOps::RVec<Float_t> vec_trkDzError1, ROOT::VecOps::RVec<Float_t> vec_trkPt, ROOT::VecOps::RVec<Float_t> vec_trkPtError, ROOT::VecOps::RVec<Int_t> vec_ToReturn){

   std::vector<int> v;
   for(int i=0; i<Ntrks; i++){
      if(std::abs(vec_trkDxy1[i]/vec_trkDxyError1[i])<5 && std::abs(vec_trkDz1[i]/vec_trkDzError1[i])<5 && vec_trkPtError[i]/vec_trkPt[i]<0.1)
         v.push_back(vec_ToReturn[i]);
   }
   return v;
}


'''

ROOT.gInterpreter.Declare(getVector_code)

#d1_ = ROOT.RDataFrame(treeName1, fileName)
#d1_ = ROOT.RDataFrame(treeName1, names)
#d1 = d1_.Range(10) ##only to test
d1 = ROOT.RDataFrame(treeName1, names)
#d2_ = ROOT.RDataFrame(treeName2, fileName)
#d2_ = ROOT.RDataFrame(treeName2, names)
#d2 = d2_.Range(10) ##only to test
d2 = ROOT.RDataFrame(treeName2, names)
#d3_ = ROOT.RDataFrame(treeName3, fileName)
d3_ = ROOT.RDataFrame(treeName3, names)
d3 = d3_.Redefine("trkPt","getVectorFloat_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkPt.trkPt)") \
        .Redefine("trkEta","getVectorFloat_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkEta.trkEta)") \
        .Redefine("trkPhi","getVectorFloat_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkPhi.trkPhi)") \
        .Redefine("trkNHit","getVectorInt_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkNHit.trkNHit)") \
        .Redefine("trkNPixelHit","getVectorInt_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkNPixelHit.trkNPixelHit)") \
        .Redefine("trkNlayer","getVectorInt_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkNlayer.trkNlayer)") \
        .Redefine("highPurity","getVectorInt_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,highPurity.highPurity)") \
        .Redefine("trkChi2","getVectorFloat_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkChi2.trkChi2)") \
        .Redefine("trkNdof","getVectorInt_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkNdof.trkNdof)") \
        .Redefine("pfEcal","getVectorFloat_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,pfEcal.pfEcal)") \
        .Redefine("pfHcal","getVectorFloat_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,pfHcal.pfHcal)") \
        .Redefine("trkDxy1","getVectorFloat_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkDxy1.trkDxy1)") \
        .Redefine("trkDxyError1","getVectorFloat_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkDxyError1.trkDxyError1)") \
        .Redefine("trkDz1","getVectorFloat_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkDz1.trkDz1)") \
        .Redefine("trkDzError1","getVectorFloat_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkDzError1.trkDzError1)") \
        .Redefine("trkPtError","getVectorFloat_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkPtError.trkPtError)")
#        .Redefine("trkPtError","getVectorFloat_Aux(nTrk.nTrk,trkDxy1.trkDxy1,trkDxyError1.trkDxyError1,trkDz1.trkDz1,trkDzError1.trkDzError1,trkPt.trkPt,trkPtError.trkPtError,trkPtError.trkPtError)") \
#        .Range(10) ##only to test

opts = ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"

d1.Snapshot(treeName1,outFileName,{"run","evt","lumi","vx","vy","vz","hiBin"})
d2.Snapshot(treeName2,outFileName,{"HBHENoiseFilterResultRun2Loose","HBHENoiseFilterResultRun2Tight","pPAprimaryVertexFilter","pBeamScrapingFilter","phfCoincFilter1","phfCoincFilter2","phfCoincFilter3","phfCoincFilter4","phfCoincFilter5"},opts)
d3.Snapshot(treeName3,outFileName,{"nTrk","trkPt","trkEta","trkPhi","trkPtError","trkNHit","trkNPixelHit","trkNlayer","highPurity","trkChi2","trkNdof","trkDxy1","trkDxyError1","trkDz1","trkDzError1","pfEcal","pfHcal"},opts)
