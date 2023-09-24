

void plot_CMN_dependence_on_LumiSec(){


////TFile *file = TFile::Open("file_cmn_Run326854_tot.root");
////TFile *file = TFile::Open("file_cmn_Run326855_tot.root");
////TFile *file = TFile::Open("file_cmn_Run326941_tot.root");
TFile *file = TFile::Open("file_cmn_Run326942_tot.root");


TTree *tree = (TTree*)file->Get("siStripCMNanalyzer/Miscellanea/track_tree");

tree->SetBranchStatus("*",0);  // disable all branches
tree->SetBranchStatus("vec_APVCMN",1); //activate this
tree->SetBranchStatus("vec_fedID",1); //activate this
tree->SetBranchStatus("vec_lumiS",1); //activate this


//TCanvas *c1 = new TCanvas();
//c1->cd();
//tree->Draw("vec_APVCMN:vec_lumiS","vec_APVCMN!=-1024 && vec_fedID==135"); ///IMPORTANT

TCanvas *c2 = new TCanvas();
c2->cd();
TH2F *hist_APVCMN_vs_lumiS = new TH2F("hist_APVCMN_vs_lumiS","hist_APVCMN_vs_lumiS",501,-0.5,500.5,100,-200,200);
//TH2F *hist_APVCMN_vs_lumiS = new TH2F("hist_APVCMN_vs_lumiS","hist_APVCMN_vs_lumiS",1,0,1,100,-200,200);
//hist_APVCMN_vs_lumiS->SetCanExtend(TH2::kXaxis);
////tree->Draw("vec_APVCMN:vec_lumiS>>hist_APVCMN_vs_lumiS","vec_APVCMN!=-1024 && vec_fedID==135","goff"); ///IMPORTANT
////tree->Draw("vec_APVCMN:vec_lumiS>>hist_APVCMN_vs_lumiS","vec_APVCMN!=-1024 && vec_fedID==137","goff"); ///IMPORTANT
tree->Draw("vec_APVCMN:vec_lumiS>>hist_APVCMN_vs_lumiS","vec_APVCMN!=-1024 && vec_fedID==80","goff"); ///IMPORTANT
TProfile * prof = (TProfile *) hist_APVCMN_vs_lumiS->ProfileX("profX",1,-1,"s");
prof->Draw("pError");
prof->SetMarkerStyle(21);
prof->SetMarkerColor(1);
prof->SetLineColor(1);
prof->SetLineWidth(2);
////prof->SetTitle("Run 326942 - FED135"); ///IMPORTANT 
////prof->SetTitle("Run 326855 - FED137"); ///IMPORTANT
prof->SetTitle("Run 326942 - FED80"); ///IMPORTANT
prof->GetXaxis()->SetTitle("LumiSection");
prof->GetYaxis()->SetTitle("Average APV CMN [adc]");

//TCanvas *c3 = new TCanvas();
//c3->cd();
//tree->Draw("vec_APVCMN","vec_APVCMN!=-1024 && vec_fedID==135 && vec_lumiS==236"); ///IMPORTANT -- ADD OTHER LUMISECTIONS??



}
