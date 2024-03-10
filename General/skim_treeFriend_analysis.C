///To run it, please, do: <<<root -l -b -q skim_treeFriend_analysis.C>>>
{

ROOT::EnableImplicitMT();

auto fileName = "HiForestAOD_DATA2011UCC_test.root";
auto outFileName = "skim_treeUCCAnalysis.root";
auto treeName = "demo/HBT";
auto treeFriendName = "hltanalysis/HltTree";

auto f = TFile::Open(fileName);
auto t = f->Get<TTree>(treeName);
auto ft = f->Get<TTree>(treeFriendName);

t->AddFriend(ft, "hlt");
ROOT::RDataFrame df (*t);

auto entries = df.Filter("hlt.HLT_HIMinBiasHfOrBSC_v1==1").Count();
std::cout<<" entries : "<< *entries <<std::endl;

auto df_cut  = df.Filter("hlt.HLT_HIMinBiasHfOrBSC_v1==1"); //select UCC events

// We are not forced to write the full set of column names. We can also
// specify a regular expression for that. In case nothing is specified, all
// columns are persistified.
df_cut.Snapshot(treeName, outFileName, {"evRunNumber","evEventNumber","Npv","pvNDOF","pvZ","pvRho","HFsumETPlus","HFsumETMinus","HFsumET","zdcSum","Ntrk","trkPt","trkEta","trkPhi","trkPtRes","trkDzSig","trkDxySig","trkNpixLayers","NSSpair","NOSpair","qinvSigSS","coulombWSS","qinvSigOS","coulombWOS"});

}
