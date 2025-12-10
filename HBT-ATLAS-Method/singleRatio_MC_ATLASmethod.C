using namespace std;

#include "CMS_lumi.C"

Double_t reject_range_min1 = 0.55; //old (pp) was 0.40
Double_t reject_range_max1 = 0.90; //old (pp) was the same
Double_t reject_range_min2 = 0.20; //old (pp) was the same
Double_t reject_range_max2 = 0.30; //old (pp) was the same
Double_t reject_range_min3 = 0.35; //old (pp) was 0.95
Double_t reject_range_max3 = 0.50; //old (pp) was 1.20

Double_t ref_levy(Double_t* x, Double_t* par){
Double_t v = 0;
if( (reject_range_min1<x[0] && x[0]<reject_range_max1) || (reject_range_min2<x[0] && x[0]<reject_range_max2) || (reject_range_min3<x[0] && x[0]<reject_range_max3) )
{
TF1::RejectPoint();
v= par[0]*(1 + par[1]*exp(-pow(fabs(par[2]*x[0])/0.1973,2)))*(1+par[3]*x[0]);
}
else{
v= par[0]*(1 + par[1]*exp(-pow(fabs(par[2]*x[0])/0.1973,2)))*(1+par[3]*x[0]);
}
return v;
}


void singleRatio_MC_ATLASmethod(){

  gStyle->SetOptStat("n");
  gStyle->SetOptTitle(0);
  
  ///TFile *file = TFile::Open("../ROOTFilesForChecks/pO_Data_5440GeV_norefsample_20250701_HiForest_dpmjetlll_Op.root");
  TFile *file = TFile::Open("../ROOTFilesForChecks/merged_pO_MC_Hijing.root"); int MC = 1;
//  TFile *file = TFile::Open("../ROOTFilesForChecks/merged_MC_Pythia_v9.root"); int MC = 2;

  //Chose variable for numerator and denominator of SR
  const unsigned int nvar=2;
  string vec_var_num_name[nvar] = {"sch","dch"};
  string vec_var_num[nvar] = {"hist_qschncVsKtVsNch_ntrkoff","hist_qdchncVsKtVsNch_ntrkoff"};
  string vec_var_den_name[nvar] = {"schmix","dchmix"};
  string vec_var_den[nvar] = {"hs_qinv_emixschVsKtVsntrkoff","hs_qinv_emixdchVsKtVsntrkoff"};

  //kT range
  //const unsigned int nktbins=2;
  //Double_t kT_min[nktbins] = {0.00,0.35};
  //Double_t kT_max[nktbins] = {0.35,1.00};
  //string kTrange[nktbins] = {", 0.00<k_{T}<0.35", ", 0.35<k_{T}<1.00"};

  const unsigned int nktbins=4;
  Double_t kT_min[nktbins] = {0.20,0.30,0.40,0.60};
  Double_t kT_max[nktbins] = {0.30,0.40,0.60,1.00};
  string kTrange[nktbins] = {", 0.20<k_{T}<0.30", ", 0.30<k_{T}<0.40", ", 0.40<k_{T}<0.60", ", 0.60<k_{T}<1.00"};

  //multiplcity range
  const unsigned int n_NtrkoffBins = 3;
  Double_t Ntrkoff_min[n_NtrkoffBins] = {0,30,60};
  Double_t Ntrkoff_max[n_NtrkoffBins] = {30,60,250};
  string Ntrkoffrange[n_NtrkoffBins] = {"0to30","30to60","60to250"};

  TH1D *matrix_sch_histos[n_NtrkoffBins][nktbins];
  TH1D *matrix_dch_histos[n_NtrkoffBins][nktbins];

  THnSparseD *aux_hsparse_sch_num; THnSparseD *aux_hsparse_sch_den;
  THnSparseD *aux_hsparse_dch_num; THnSparseD *aux_hsparse_dch_den;
  TH1D *aux_histo_sch_num; TH1D *aux_histo_sch_den;
  TH1D *aux_histo_dch_num; TH1D *aux_histo_dch_den;

  for(unsigned int i_ntrkoff=0; i_ntrkoff<n_NtrkoffBins; i_ntrkoff++){

     for(unsigned int i_kt=0; i_kt<nktbins; i_kt++){

	aux_hsparse_sch_num = dynamic_cast<THnSparseD*>(file->Get(("SameEventHistograms/"+vec_var_num[0]).c_str()));
        aux_hsparse_sch_den = dynamic_cast<THnSparseD*>(file->Get(("MixedEventHistograms/"+vec_var_den[0]).c_str()));
        aux_hsparse_dch_num = dynamic_cast<THnSparseD*>(file->Get(("SameEventHistograms/"+vec_var_num[1]).c_str()));
        aux_hsparse_dch_den = dynamic_cast<THnSparseD*>(file->Get(("MixedEventHistograms/"+vec_var_den[1]).c_str()));  	

        TAxis *aux_axisNtrk_sch_num = aux_hsparse_sch_num->GetAxis(2);
	TAxis *aux_axisNtrk_sch_den = aux_hsparse_sch_den->GetAxis(2);
        int binNtrkmin_sch = aux_axisNtrk_sch_num->FindBin(Ntrkoff_min[i_ntrkoff]+0.001);
        int binNtrkmax_sch = aux_axisNtrk_sch_num->FindBin(Ntrkoff_max[i_ntrkoff]-0.001);
        aux_axisNtrk_sch_num->SetRange(binNtrkmin_sch, binNtrkmax_sch);
	aux_axisNtrk_sch_den->SetRange(binNtrkmin_sch, binNtrkmax_sch);

        TAxis *aux_axisNtrk_dch_num = aux_hsparse_dch_num->GetAxis(2);
	TAxis *aux_axisNtrk_dch_den = aux_hsparse_dch_den->GetAxis(2);
        int binNtrkmin_dch = aux_axisNtrk_dch_num->FindBin(Ntrkoff_min[i_ntrkoff]+0.001);
        int binNtrkmax_dch = aux_axisNtrk_dch_num->FindBin(Ntrkoff_max[i_ntrkoff]-0.001);
        aux_axisNtrk_dch_num->SetRange(binNtrkmin_dch, binNtrkmax_dch);
	aux_axisNtrk_dch_den->SetRange(binNtrkmin_dch, binNtrkmax_dch);

        TAxis *aux_axiskT_sch_num = aux_hsparse_sch_num->GetAxis(1);
        TAxis *aux_axiskT_sch_den = aux_hsparse_sch_den->GetAxis(1);
        int binkTmin_sch = aux_axiskT_sch_num->FindBin(kT_min[i_kt]+0.001);
        int binkTmax_sch = aux_axiskT_sch_num->FindBin(kT_max[i_kt]-0.001);
        aux_axiskT_sch_num->SetRange(binkTmin_sch, binkTmax_sch);
        aux_axiskT_sch_den->SetRange(binkTmin_sch, binkTmax_sch);

        TAxis *aux_axiskT_dch_num = aux_hsparse_dch_num->GetAxis(1);
        TAxis *aux_axiskT_dch_den = aux_hsparse_dch_den->GetAxis(1);
        int binkTmin_dch = aux_axiskT_dch_num->FindBin(kT_min[i_kt]+0.001);
        int binkTmax_dch = aux_axiskT_dch_num->FindBin(kT_max[i_kt]-0.001);
        aux_axiskT_dch_num->SetRange(binkTmin_dch, binkTmax_dch);
        aux_axiskT_dch_den->SetRange(binkTmin_dch, binkTmax_dch);


        aux_histo_sch_num = aux_hsparse_sch_num->Projection(0); 
	aux_histo_sch_num->SetName(Form("SR_sch_num_NtrkoffRange_%dto%d_ktRange_%fto%f",(Int_t)Ntrkoff_min[i_ntrkoff],(Int_t)Ntrkoff_max[i_ntrkoff],kT_min[i_kt],kT_max[i_kt]));
	aux_histo_sch_num->SetTitle(Form("SR_sch_num_NtrkoffRange_%dto%d_ktRange_%fto%f",(Int_t)Ntrkoff_min[i_ntrkoff],(Int_t)Ntrkoff_max[i_ntrkoff],kT_min[i_kt],kT_max[i_kt]));

	aux_histo_sch_den = aux_hsparse_sch_den->Projection(0); 
        aux_histo_sch_den->SetName(Form("SR_sch_den_NtrkoffRange_%dto%d_ktRange_%fto%f",(Int_t)Ntrkoff_min[i_ntrkoff],(Int_t)Ntrkoff_max[i_ntrkoff],kT_min[i_kt],kT_max[i_kt]));
	aux_histo_sch_den->SetTitle(Form("SR_sch_den_NtrkoffRange_%dto%d_ktRange_%fto%f",(Int_t)Ntrkoff_min[i_ntrkoff],(Int_t)Ntrkoff_max[i_ntrkoff],kT_min[i_kt],kT_max[i_kt]));

	aux_histo_dch_num = aux_hsparse_dch_num->Projection(0); 
	aux_histo_dch_num->SetName(Form("SR_dch_num_NtrkoffRange_%dto%d_ktRange_%fto%f",(Int_t)Ntrkoff_min[i_ntrkoff],(Int_t)Ntrkoff_max[i_ntrkoff],kT_min[i_kt],kT_max[i_kt]));
	aux_histo_dch_num->SetTitle(Form("SR_dch_num_NtrkoffRange_%dto%d_ktRange_%fto%f",(Int_t)Ntrkoff_min[i_ntrkoff],(Int_t)Ntrkoff_max[i_ntrkoff],kT_min[i_kt],kT_max[i_kt]));

	aux_histo_dch_den = aux_hsparse_dch_den->Projection(0); 
        aux_histo_dch_den->SetName(Form("SR_dch_den_NtrkoffRange_%dto%d_ktRange_%fto%f",(Int_t)Ntrkoff_min[i_ntrkoff],(Int_t)Ntrkoff_max[i_ntrkoff],kT_min[i_kt],kT_max[i_kt]));
	aux_histo_dch_den->SetTitle(Form("SR_dch_den_NtrkoffRange_%dto%d_ktRange_%fto%f",(Int_t)Ntrkoff_min[i_ntrkoff],(Int_t)Ntrkoff_max[i_ntrkoff],kT_min[i_kt],kT_max[i_kt]));


        Int_t bin_for_normInt_min = aux_histo_sch_num->GetXaxis()->FindBin(1.2001);
        Int_t bin_for_normInt_max = aux_histo_sch_num->GetXaxis()->FindBin(1.9999);
        Double_t int_sch_num_controlRegion = aux_histo_sch_num->Integral(bin_for_normInt_min,bin_for_normInt_max);
        Double_t int_sch_den_controlRegion = aux_histo_sch_den->Integral(bin_for_normInt_min,bin_for_normInt_max);
        aux_histo_sch_den->Scale(int_sch_num_controlRegion/int_sch_den_controlRegion);
        Double_t int_dch_num_controlRegion = aux_histo_dch_num->Integral(bin_for_normInt_min,bin_for_normInt_max);
        Double_t int_dch_den_controlRegion = aux_histo_dch_den->Integral(bin_for_normInt_min,bin_for_normInt_max);
        aux_histo_dch_den->Scale(int_dch_num_controlRegion/int_dch_den_controlRegion);
	  
        aux_histo_sch_num->Divide(aux_histo_sch_den);
        aux_histo_dch_num->Divide(aux_histo_dch_den);

	matrix_sch_histos[i_ntrkoff][i_kt] = aux_histo_sch_num;
	matrix_dch_histos[i_ntrkoff][i_kt] = aux_histo_dch_num;
  
     }//kT bins loop

  }//Ntrkoff bins loop

TCanvas* tcanvas[n_NtrkoffBins][nktbins];
TF1* tf1_sch[n_NtrkoffBins][nktbins];
TFitResultPtr res_sch[n_NtrkoffBins][nktbins];
TF1* tf1_dch[n_NtrkoffBins][nktbins];
TFitResultPtr res_dch[n_NtrkoffBins][nktbins];
auto tg_logBA = new TGraphErrors(n_NtrkoffBins*nktbins);
auto tg_sigmaBA = new TGraphErrors(n_NtrkoffBins*nktbins);
int ii=0;
for(unsigned int j_ntrkoff=0; j_ntrkoff<n_NtrkoffBins; j_ntrkoff++){
   for(unsigned int j_kt=0; j_kt<nktbins; j_kt++){
      tcanvas[j_ntrkoff][j_kt] = new TCanvas();
      tcanvas[j_ntrkoff][j_kt]->cd();
      tf1_sch[j_ntrkoff][j_kt] = new TF1(Form("FitFunc_SR_sch_NtrkoffRange%d_ktRange_%d",(Int_t)j_ntrkoff,(Int_t)j_kt),ref_levy,0.1,2.0,4);
      tf1_sch[j_ntrkoff][j_kt]->SetParameters(1.0,0.2,.3,0.0);
      tf1_sch[j_ntrkoff][j_kt]->SetParLimits(0,0.5,1.5);
      tf1_sch[j_ntrkoff][j_kt]->SetParLimits(1,0.0,1.0);
      tf1_sch[j_ntrkoff][j_kt]->SetParLimits(2,0.0,1.0);
      tf1_sch[j_ntrkoff][j_kt]->SetLineColor(kBlack); 
      tf1_sch[j_ntrkoff][j_kt]->SetLineWidth(3);      
      tf1_dch[j_ntrkoff][j_kt] = new TF1(Form("FitFunc_SR_dch_NtrkoffRange%d_ktRange_%d",(Int_t)j_ntrkoff,(Int_t)j_kt),ref_levy,0.1,2.0,4);
      tf1_dch[j_ntrkoff][j_kt]->SetParameters(1.0,0.2,.3,0.0);
      tf1_dch[j_ntrkoff][j_kt]->SetParLimits(0,0.5,1.5);
      tf1_dch[j_ntrkoff][j_kt]->SetParLimits(1,0.0,1.0);
      tf1_dch[j_ntrkoff][j_kt]->SetParLimits(2,0.0,1.0);
      tf1_dch[j_ntrkoff][j_kt]->SetLineColor(kRed);
      tf1_dch[j_ntrkoff][j_kt]->SetLineWidth(3);
      matrix_dch_histos[j_ntrkoff][j_kt]->Draw("pE");
      matrix_dch_histos[j_ntrkoff][j_kt]->SetTitle("");
      matrix_dch_histos[j_ntrkoff][j_kt]->GetXaxis()->SetTitle("q_{inv} [GeV]");
      matrix_dch_histos[j_ntrkoff][j_kt]->GetYaxis()->SetTitle("Single Ratio");
      matrix_dch_histos[j_ntrkoff][j_kt]->SetStats(0);
      matrix_dch_histos[j_ntrkoff][j_kt]->SetLineColor(2);
      matrix_dch_histos[j_ntrkoff][j_kt]->SetMarkerColor(2);
      matrix_dch_histos[j_ntrkoff][j_kt]->GetXaxis()->SetRangeUser(0.0,2.0);
      matrix_dch_histos[j_ntrkoff][j_kt]->GetYaxis()->SetRangeUser(0.96,1.50);
      matrix_sch_histos[j_ntrkoff][j_kt]->Draw("pEsame");
      matrix_sch_histos[j_ntrkoff][j_kt]->SetLineColor(1);
      matrix_sch_histos[j_ntrkoff][j_kt]->SetMarkerColor(1);
      ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
      res_sch[j_ntrkoff][j_kt] = matrix_sch_histos[j_ntrkoff][j_kt]->Fit(tf1_sch[j_ntrkoff][j_kt], "S R");
      tf1_sch[j_ntrkoff][j_kt]->Draw("Lsame");
      res_dch[j_ntrkoff][j_kt] = matrix_dch_histos[j_ntrkoff][j_kt]->Fit(tf1_dch[j_ntrkoff][j_kt], "S R");
      tf1_dch[j_ntrkoff][j_kt]->Draw("Lsame");
      auto leg = new TLegend(0.6,0.75,0.85,0.85);
      leg->SetTextSize(0.04);
      leg->SetFillColor(0);
      leg->SetBorderSize(0);
      leg->AddEntry(matrix_dch_histos[j_ntrkoff][j_kt],"OS","PL");
      leg->AddEntry(matrix_sch_histos[j_ntrkoff][j_kt],"SS","PL");
      leg->Draw("same");
      auto latex = new TLatex();
      latex->SetNDC();
      latex->SetTextSize(0.04);
      latex->SetTextFont(42);
      latex->DrawLatex(0.15,0.83,Form("%.0f < N_{trk}^{off} < %.0f, %.2f < k_{T} < %.2f",Ntrkoff_min[j_ntrkoff],Ntrkoff_max[j_ntrkoff],kT_min[j_kt],kT_max[j_kt]));

      double aux_sch_BA = tf1_sch[j_ntrkoff][j_kt]->GetParameter(1);
      double aux_sch_BA_err = tf1_sch[j_ntrkoff][j_kt]->GetParError(1);
      double aux_sch_BA_log = TMath::Log(aux_sch_BA);
      double aux_sch_BA_log_err = aux_sch_BA_err/aux_sch_BA;
      double aux_dch_BA = tf1_dch[j_ntrkoff][j_kt]->GetParameter(1);
      double aux_dch_BA_err = tf1_dch[j_ntrkoff][j_kt]->GetParError(1);
      double aux_dch_BA_log = TMath::Log(aux_dch_BA);
      double aux_dch_BA_log_err = aux_dch_BA_err/aux_dch_BA;
      tg_logBA->SetPoint(ii,aux_dch_BA_log,aux_sch_BA_log);
      tg_logBA->SetPointError(ii,aux_dch_BA_log_err,aux_sch_BA_log_err);

      double aux_sch_sigmaBA = tf1_sch[j_ntrkoff][j_kt]->GetParameter(2);      
      double aux_sch_sigmaBA_err = tf1_sch[j_ntrkoff][j_kt]->GetParError(2);
      double aux_dch_sigmaBA = tf1_dch[j_ntrkoff][j_kt]->GetParameter(2);
      double aux_dch_sigmaBA_err = tf1_dch[j_ntrkoff][j_kt]->GetParError(2);
      tg_sigmaBA->SetPoint(ii,aux_dch_sigmaBA,aux_sch_sigmaBA);
      tg_sigmaBA->SetPointError(ii,aux_dch_sigmaBA_err,aux_sch_sigmaBA_err);     
 
      ii++;
   }
}   

auto canvas_BA = new TCanvas();
canvas_BA->cd();
tg_logBA->Draw("Ap");
tg_logBA->SetMarkerStyle(20);
tg_logBA->GetXaxis()->SetTitle("log(B_{A})^{+-}");
tg_logBA->GetYaxis()->SetTitle("log(B_{A})^{#pm#pm}");
TF1* fitFunc_logBA = new TF1("fitFunc_logBA", "[0] + [1]*x", -3.5, -1.0);
fitFunc_logBA->SetParNames("Intercept_logBA", "Slope_logBA");
fitFunc_logBA->SetLineColor(kRed);
fitFunc_logBA->SetLineWidth(2);
TFitResultPtr fitResult_logBA = tg_logBA->Fit(fitFunc_logBA, "RS");  // R = use range, S = save result
double intercept_logBA = fitFunc_logBA->GetParameter(0);
double slope_logBA = fitFunc_logBA->GetParameter(1);
double intercept_err_logBA = fitFunc_logBA->GetParError(0);
double slope_err_logBA = fitFunc_logBA->GetParError(1);
double chi2_logBA = fitFunc_logBA->GetChisquare();
int ndf_logBA = fitFunc_logBA->GetNDF();
TLegend* leg_logBA = new TLegend(0.15, 0.65, 0.60, 0.85);
if(MC==1) leg_logBA->AddEntry(tg_logBA, "MC HIJING", "pe");
if(MC==2) leg_logBA->AddEntry(tg_logBA, "MC PYTHIA", "pe");
leg_logBA->AddEntry(fitFunc_logBA, Form("Fit: y = (%.3f #pm %.3f) + (%.3f #pm %.3f)x", intercept_logBA, intercept_err_logBA, slope_logBA, slope_err_logBA), "l");
leg_logBA->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f/%d", chi2_logBA, ndf_logBA), "");
leg_logBA->Draw();



auto canvas_sigBA = new TCanvas();
canvas_sigBA->cd();
tg_sigmaBA->Draw("Ap");
tg_sigmaBA->SetMarkerStyle(20);
tg_sigmaBA->GetXaxis()->SetTitle("#sigma_{B,A}^{+-}");
tg_sigmaBA->GetYaxis()->SetTitle("#sigma_{B,A}^{#pm#pm}");
TF1* fitFunc_sigmaBA = new TF1("fitFunc_sigmaBA", "[0] + [1]*x", 0.23, 0.32);
fitFunc_sigmaBA->SetParNames("Intercept_sigmaBA", "Slope_sigmaBA");
fitFunc_sigmaBA->SetLineColor(kRed);
fitFunc_sigmaBA->SetLineWidth(2);
TFitResultPtr fitResult_sigmaBA = tg_sigmaBA->Fit(fitFunc_sigmaBA, "RS");  // R = use range, S = save result
double intercept_sigmaBA = fitFunc_sigmaBA->GetParameter(0);
double slope_sigmaBA = fitFunc_sigmaBA->GetParameter(1);
double intercept_err_sigmaBA = fitFunc_sigmaBA->GetParError(0);
double slope_err_sigmaBA = fitFunc_sigmaBA->GetParError(1);
double chi2_sigmaBA = fitFunc_sigmaBA->GetChisquare();
int ndf_sigmaBA = fitFunc_sigmaBA->GetNDF();
TLegend* leg_sigmaBA = new TLegend(0.15, 0.65, 0.60, 0.85);
if(MC==1) leg_sigmaBA->AddEntry(tg_sigmaBA, "MC HIJING", "pe");
if(MC==2) leg_sigmaBA->AddEntry(tg_sigmaBA, "MC PYTHIA", "pe");
leg_sigmaBA->AddEntry(fitFunc_sigmaBA, Form("Fit: y = (%.3f #pm %.3f) + (%.3f #pm %.3f)x", intercept_sigmaBA, intercept_err_sigmaBA, slope_sigmaBA, slope_err_sigmaBA), "l");
leg_sigmaBA->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f/%d", chi2_sigmaBA, ndf_sigmaBA), "");
leg_sigmaBA->Draw();


}
