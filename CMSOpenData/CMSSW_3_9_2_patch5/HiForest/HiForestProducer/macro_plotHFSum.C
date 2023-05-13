
void macro_plotHFSum(){

auto fileName = "HiForestAOD_DATAtest.root";
auto outFileName = "HFSumEt.root";
TFile *f = TFile::Open(fileName);
TFile *output = new TFile(outFileName,"recreate");


auto treeName_1 = "demo/Muons";
auto treeName_2 = "hltanalysis/HltTree";

auto t_1 = f->Get<TTree>(treeName_1);
auto t_2 = f->Get<TTree>(treeName_2);

t_1->AddFriend(t_2, "trg_sel");

ROOT::RDataFrame d(*t_1);

std::string trg_selection = "HLT_HIMinBiasHfOrBSC_Core==1";
auto d_select = d.Filter(trg_selection,"Trigger Selection")
	         .Define("HFsumET_norm","HFsumET/1000.0");


auto h_hfsum = d_select.Histo1D(TH1D("hfSum","CMS: PbPb at #sqrt{s_{NN}} = 2.76 TeV",100,0,5), "HFsumET_norm");
h_hfsum->GetXaxis()->SetTitle("#Sum{E_{T}} [TeV]");
h_hfsum->GetYaxis()->SetTitle("Number of events");
h_hfsum->GetXaxis()->SetTitleOffset(1.15);
h_hfsum->GetYaxis()->SetTitleOffset(1.15);
auto c = new TCanvas("c", "", 500, 500);
h_hfsum->DrawClone("pError"); 
gPad->SetLogy(1);
h_hfsum->Write();
c->Update(); 



}
