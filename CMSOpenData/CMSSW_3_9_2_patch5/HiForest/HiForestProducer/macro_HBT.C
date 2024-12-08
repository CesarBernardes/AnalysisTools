

//to reject a range in the fit
Double_t reject_range_min = 0.0;
Double_t reject_range_max = 0.00001;


Double_t func1_exp(Double_t* x, Double_t* par){
Double_t v = 0;
if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
else{v= par[0]*(1 + par[1]*exp(-par[2]*x[0]/0.1973))*(1+par[3]*x[0]);}
return v;
}

Double_t func2_gauss(Double_t* x, Double_t* par){
Double_t v = 0;
if(reject_range_min<x[0] && x[0]<reject_range_max){TF1::RejectPoint();}
else{v= par[0]*(1 + par[1]*exp(-pow(par[2]*x[0]/0.1973,2.)))*(1+par[3]*x[0]);}
return v;
}

void macro_HBT(){

auto fileName = "HiForestAOD_DATAtest_500k.root";
auto outFileName = "HBT_histos.root";
TFile *f = TFile::Open(fileName);
TFile *output = new TFile(outFileName,"recreate");

auto treeName_1 = "demo/HBT";
auto treeName_2 = "hltanalysis/HltTree";

auto t_1 = f->Get<TTree>(treeName_1);
auto t_2 = f->Get<TTree>(treeName_2);

t_1->AddFriend(t_2, "trg_sel");

TH1D * h_qinv_sig_SS = new TH1D("hqinvSigSS","q_{inv} SameEvt SS tracks",50,0,1);
TH1D * h_qinv_sig_OS = new TH1D("hqinvSigOS","q_{inv} SameEvt OS tracks",50,0,1);
TH1D * h_qinv_sig_SS_Ccorr = new TH1D("hqinvSigSS_Ccorr","q_{inv} SameEvt SS tracks - Coulomb Corr",50,0,1);
TH1D * h_qinv_sig_OS_Ccorr = new TH1D("hqinvSigOS_Ccorr","q_{inv} SameEvt OS tracks - Coulomb Corr",50,0,1);

auto c = new TCanvas("c", "", 500, 500);
c->cd();
///See corresponding centrality "HFsumET-->centrality" from Fig.1 of https://arxiv.org/pdf/1107.4800.pdf
//t_1->Draw("qinvSigSS>>hqinvSigSS","HLT_HIMinBiasHfOrBSC_Core==1","goff");
///t_1->Draw("qinvSigSS>>hqinvSigSS","HFsumET>2300","goff"); //grupo4 --> centrality: centrality: ~0-15%
////t_1->Draw("qinvSigSS>>hqinvSigSS","HFsumET>1000 && HFsumET<1850","goff"); //grupo2 --> centrality: ~20-35%
t_1->Draw("qinvSigSS>>hqinvSigSS","HFsumET<800","goff"); //grupo3 --> centrality: ~40-100%
////t_1->Draw("qinvSigSS>>hqinvSigSS","HFsumET>3400","goff"); //grupo1 --> centrality: ~0-5%
h_qinv_sig_SS->Write();
h_qinv_sig_SS->GetXaxis()->SetTitle("q_{inv} [GeV]");
h_qinv_sig_SS->GetYaxis()->SetTitle("Number of pairs / bin");
h_qinv_sig_SS->GetXaxis()->SetTitleOffset(1.15);
h_qinv_sig_SS->GetYaxis()->SetTitleOffset(1.15);
gPad->SetLogy(1);
c->Update();

//t_1->Draw("qinvSigOS>>hqinvSigOS","HLT_HIMinBiasHfOrBSC_Core==1","goff");
///t_1->Draw("qinvSigOS>>hqinvSigOS","HFsumET>2300","goff"); //grupo4
////t_1->Draw("qinvSigOS>>hqinvSigOS","HFsumET>1000 && HFsumET<1850","goff"); //grupo2
t_1->Draw("qinvSigOS>>hqinvSigOS","HFsumET<800","goff");	  //grupo3
////t_1->Draw("qinvSigOS>>hqinvSigOS","HFsumET>3400","goff"); //grupo1
h_qinv_sig_OS->Write();
c->Update();

//t_1->Draw("qinvSigSS>>hqinvSigSS_Ccorr","coulombWSS*(HLT_HIMinBiasHfOrBSC_Core==1)","goff");
///t_1->Draw("qinvSigSS>>hqinvSigSS_Ccorr","coulombWSS*(HFsumET>2300)","goff"); //grupo4
////t_1->Draw("qinvSigSS>>hqinvSigSS_Ccorr","coulombWSS*(HFsumET>1000 && HFsumET<1850)","goff"); //grupo2
t_1->Draw("qinvSigSS>>hqinvSigSS_Ccorr","coulombWSS*(HFsumET<800)","goff"); //grupo3
////t_1->Draw("qinvSigSS>>hqinvSigSS_Ccorr","coulombWSS*(HFsumET>3400)","goff"); //grupo1
h_qinv_sig_SS_Ccorr->Write();
c->Update();

//t_1->Draw("qinvSigOS>>hqinvSigOS_Ccorr","coulombWOS*(HLT_HIMinBiasHfOrBSC_Core==1)","goff");
///t_1->Draw("qinvSigOS>>hqinvSigOS_Ccorr","coulombWOS*(HFsumET>2300)","goff"); //grupo4
////t_1->Draw("qinvSigOS>>hqinvSigOS_Ccorr","coulombWOS*(HFsumET>1000 && HFsumET<1850)","goff"); //grupo2
t_1->Draw("qinvSigOS>>hqinvSigOS_Ccorr","coulombWOS*(HFsumET<800)","goff"); //grupo3
////t_1->Draw("qinvSigOS>>hqinvSigOS_Ccorr","coulombWOS*(HFsumET>3400)","goff"); //grupo1
h_qinv_sig_OS_Ccorr->Write();
c->Update();

TH1D * h_qinv_sig_SS_Ccorr_clone = (TH1D *)h_qinv_sig_SS_Ccorr->Clone("h_qinv_sig_SS_Ccorr_clone");
TH1D * h_qinv_sig_OS_Ccorr_clone = (TH1D *)h_qinv_sig_OS_Ccorr->Clone("h_qinv_sig_OS_Ccorr_clone");
Int_t bin_for_normInt_min = h_qinv_sig_SS_Ccorr_clone->GetXaxis()->FindBin(0.6001);
Int_t bin_for_normInt_max = h_qinv_sig_SS_Ccorr_clone->GetXaxis()->FindBin(0.9999);
Double_t int_num_controlRegion = h_qinv_sig_SS_Ccorr_clone->Integral(bin_for_normInt_min,bin_for_normInt_max);
Double_t int_den_controlRegion = h_qinv_sig_OS_Ccorr_clone->Integral(bin_for_normInt_min,bin_for_normInt_max);
h_qinv_sig_OS_Ccorr_clone->Scale(int_num_controlRegion/int_den_controlRegion);
h_qinv_sig_SS_Ccorr_clone->Divide(h_qinv_sig_OS_Ccorr_clone);
h_qinv_sig_SS_Ccorr_clone->SetLineColor(1);
h_qinv_sig_SS_Ccorr_clone->SetMarkerColor(1);

//Setting function for double ratio
TF1 *f_exp = new TF1("f_exp",func1_exp,0.02, 0.5, 4); //try 0.3 after
f_exp->SetParameters(1.0,1.0,4.0,0.0);
f_exp->SetParName(0,"Const");
f_exp->SetParLimits(0,0.0,2.0);//The specified limits will be used in a fit operation when the option "B" is specified (Bounds). To fix a parameter, use TF1::FixParameter.
f_exp->SetParName(1,"#lambda");
f_exp->SetParLimits(1,0.0,2.0);
f_exp->SetParName(2,"R (fm)");
f_exp->SetParName(3,"#delta");
f_exp->SetLineColor(kRed); 
f_exp->SetLineWidth(2);  
h_qinv_sig_SS_Ccorr_clone->Draw("pError");
h_qinv_sig_SS_Ccorr_clone->GetYaxis()->SetTitle("Single Ratio (SR)");
h_qinv_sig_SS_Ccorr_clone->GetXaxis()->SetTitle("q_{inv} [GeV]");
TFitResultPtr res_exp;
ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
res_exp = h_qinv_sig_SS_Ccorr_clone->Fit(f_exp, "S R"); 
Double_t chi2_BC_exp =  2.* res_exp->MinFcnValue(); 
cout<<"chi2_BC Exp  :  "<<chi2_BC_exp<<endl;
f_exp->Draw("SameL");

TF1 *f_gauss = new TF1("f_gauss",func2_gauss,0.02, 0.5, 4); //try 0.3 after
f_gauss->SetParameters(1.0,1.0,4.0,0.0);
f_gauss->SetParName(0,"Const");
f_gauss->SetParLimits(0,0.0,2.0);
f_gauss->SetParName(1,"#lambda");
f_gauss->SetParLimits(1,0.0,2.0);
f_gauss->SetParName(2,"R (fm)");
f_gauss->SetParName(3,"#delta");
f_gauss->SetLineColor(kBlue); 
f_gauss->SetLineWidth(2);
TFitResultPtr res_gauss;
res_gauss = h_qinv_sig_SS_Ccorr_clone->Fit(f_gauss, "S R");
Double_t chi2_BC_gauss =  2.* res_gauss->MinFcnValue();
cout<<"chi2_BC Gauss  :  "<<chi2_BC_gauss<<endl;
f_gauss->Draw("SameL");
gPad->SetLogy(0);
//h_qinv_sig_SS_Ccorr_clone->Write();
//f_exp->Write();
//f_gauss->Write();
c->Update();

}	
