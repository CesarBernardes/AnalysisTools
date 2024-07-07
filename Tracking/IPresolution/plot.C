{

auto frun3 = TFile::Open("trk_ipres_Run3_eta3_noSel.root");
auto frun4 = TFile::Open("trk_ipres_Run4_eta1_and3_noSel.root");

auto h2_trkIPresXY_PtReco_barrel_run3 = (TH2F*) frun3->Get("trk_ipres/Histograms/trkIPresXY_PtReco_barrel");
auto h2_trkIPresZ_PtReco_barrel_run3 = (TH2F*) frun3->Get("trk_ipres/Histograms/trkIPresZ_PtReco_barrel");
auto h2_trkIPresXY_PtReco_forward_run3 = (TH2F*) frun3->Get("trk_ipres/Histograms/trkIPresXY_PtReco_forward");
auto h2_trkIPresZ_PtReco_forward_run3 = (TH2F*) frun3->Get("trk_ipres/Histograms/trkIPresZ_PtReco_forward");
auto h1_pt_run3 = (TH1F*) frun3->Get("trk_ipres/Histograms/TrkPt");

auto h2_trkIPresXY_PtReco_barrel_run4 = (TH2F*) frun4->Get("trk_ipres/Histograms/trkIPresXY_PtReco_barrel");
auto h2_trkIPresZ_PtReco_barrel_run4 = (TH2F*) frun4->Get("trk_ipres/Histograms/trkIPresZ_PtReco_barrel");
auto h2_trkIPresXY_PtReco_forward_run4 = (TH2F*) frun4->Get("trk_ipres/Histograms/trkIPresXY_PtReco_forward");
auto h2_trkIPresZ_PtReco_forward_run4 = (TH2F*) frun4->Get("trk_ipres/Histograms/trkIPresZ_PtReco_forward");
auto h1_pt_run4 = (TH1F*) frun4->Get("trk_ipres/Histograms/TrkPt");

int n=8;
float lowerbin[8]={0.20009,0.40009,0.60009,0.80009,1.0009,1.5009,2.0009,4.0009};
float upperbin[8]={0.39999,0.59999,0.79999,0.99999,1.4999,1.9999,3.9999,19.999};
auto tg_trkIPresXY_PtReco_barrel_run3 = new TGraphErrors(8);
auto tg_trkIPresXY_PtReco_barrel_run4 = new TGraphErrors(8);
auto tg_trkIPresZ_PtReco_barrel_run3 = new TGraphErrors(8);
auto tg_trkIPresZ_PtReco_barrel_run4 = new TGraphErrors(8);
auto tg_trkIPresXY_PtReco_forward_run3 = new TGraphErrors(8);
auto tg_trkIPresXY_PtReco_forward_run4 = new TGraphErrors(8);
auto tg_trkIPresZ_PtReco_forward_run3 = new TGraphErrors(8);
auto tg_trkIPresZ_PtReco_forward_run4 = new TGraphErrors(8);


///
auto c = new TCanvas("c","multipads",900,700);
gStyle->SetOptStat(0);
c->Divide(2,4,0,0);

auto c2 = new TCanvas("c2","multipads2",900,700);
c2->Divide(2,4,0,0);

auto c3 = new TCanvas("c3","multipads3",900,700);
c3->Divide(2,4,0,0);

auto c4 = new TCanvas("c4","multipads4",900,700);
c4->Divide(2,4,0,0);

//auto c5 = new TCanvas("c5","multipads5",900,700);
//c5->Divide(2,4,0,0);

auto c6 = new TCanvas("c6","multipads6",900,700);
c6->Divide(2,4,0,0);

//auto c7 = new TCanvas("c7","multipads7",900,700);
//c7->Divide(2,4,0,0);

auto c8 = new TCanvas("c8","multipads8",900,700);
c8->Divide(2,4,0,0);

for(int i=0; i<n; i++){	
   c->cd(i+1);
   int lbin = h1_pt_run3->FindBin(lowerbin[i]);
   int ubin = h1_pt_run3->FindBin(upperbin[i]);
   //std::cout<<" lbin : "<<lbin<<" ubin : "<<ubin<<std::endl;
   /*auto proj_trkIPresXY_PtReco_barrel_run3 = (TH1F*)h2_trkIPresXY_PtReco_barrel_run3->ProjectionY(TString::Format("trkIPresXY_PtReco_barrel_run3_%d",i),lbin,ubin);
   proj_trkIPresXY_PtReco_barrel_run3->Draw();
   proj_trkIPresXY_PtReco_barrel_run3->Fit("gaus","","",-0.04,0.04);
   TF1 *fr = (TF1*)proj_trkIPresXY_PtReco_barrel_run3->GetListOfFunctions()->FindObject("gaus");
   float sigma = fr->GetParameter(2);
   float sigmaErr = fr->GetParError(2);
   h1_pt_run3->GetXaxis()->SetRange(lbin,ubin);
   float pt = h1_pt_run3->GetMean();
   std::cout<<"ptRun3 barrel : "<<pt<<" sigma : "<<sigma<<" sigmaErr : "<<sigmaErr<<std::endl;
   tg_trkIPresXY_PtReco_barrel_run3->SetPoint(i,pt,sigma);
   tg_trkIPresXY_PtReco_barrel_run3->SetPointError(i,0.,sigmaErr);
   */
   c2->cd(i+1);
   auto proj_trkIPresXY_PtReco_barrel_run4 = (TH1F*)h2_trkIPresXY_PtReco_barrel_run4->ProjectionY(TString::Format("trkIPresXY_PtReco_barrel_run4_%d",i),lbin,ubin);
   proj_trkIPresXY_PtReco_barrel_run4->Draw();
   proj_trkIPresXY_PtReco_barrel_run4->Fit("gaus","","",-0.03,0.03);
   TF1 *fr2 = (TF1*)proj_trkIPresXY_PtReco_barrel_run4->GetListOfFunctions()->FindObject("gaus");
   float sigma2 = fr2->GetParameter(2);
   float sigmaErr2 = fr2->GetParError(2);
   h1_pt_run4->GetXaxis()->SetRange(lbin,ubin);
   float pt2 = h1_pt_run4->GetMean();
   std::cout<<"ptRun4 barrel : "<<pt2<<" sigma : "<<sigma2<<" sigmaErr : "<<sigmaErr2<<std::endl;
   tg_trkIPresXY_PtReco_barrel_run4->SetPoint(i,pt2,sigma2);
   tg_trkIPresXY_PtReco_barrel_run4->SetPointError(i,0.,sigmaErr2);

   /*c3->cd(i+1);
   auto proj_trkIPresZ_PtReco_barrel_run3 = (TH1F*)h2_trkIPresZ_PtReco_barrel_run3->ProjectionY(TString::Format("trkIPresZ_PtReco_barrel_run3_%d",i),lbin,ubin);
   proj_trkIPresZ_PtReco_barrel_run3->Draw();
   proj_trkIPresZ_PtReco_barrel_run3->Fit("gaus","","",-0.04,0.04);
   TF1 *fr3 = (TF1*)proj_trkIPresZ_PtReco_barrel_run3->GetListOfFunctions()->FindObject("gaus");
   float sigma3 = fr3->GetParameter(2);
   float sigmaErr3 = fr3->GetParError(2);
   std::cout<<"ptRun3 barrel : "<<pt<<" sigma : "<<sigma3<<" sigmaErr : "<<sigmaErr3<<std::endl;
   tg_trkIPresZ_PtReco_barrel_run3->SetPoint(i,pt,sigma3);
   tg_trkIPresZ_PtReco_barrel_run3->SetPointError(i,0.,sigmaErr3);
   */
   c4->cd(i+1);
   auto proj_trkIPresZ_PtReco_barrel_run4 = (TH1F*)h2_trkIPresZ_PtReco_barrel_run4->ProjectionY(TString::Format("trkIPresZ_PtReco_barrel_run4_%d",i),lbin,ubin);
   proj_trkIPresZ_PtReco_barrel_run4->Draw();
   proj_trkIPresZ_PtReco_barrel_run4->Fit("gaus","","",-0.03,0.03);
   TF1 *fr4 = (TF1*)proj_trkIPresZ_PtReco_barrel_run4->GetListOfFunctions()->FindObject("gaus");
   float sigma4 = fr4->GetParameter(2);
   float sigmaErr4 = fr4->GetParError(2);
   std::cout<<"ptRun4 barrel : "<<pt2<<" sigma : "<<sigma4<<" sigmaErr : "<<sigmaErr4<<std::endl;
   tg_trkIPresZ_PtReco_barrel_run4->SetPoint(i,pt2,sigma4);
   tg_trkIPresZ_PtReco_barrel_run4->SetPointError(i,0.,sigmaErr4);

  /* c5->cd(i+1);
   auto proj_trkIPresXY_PtReco_forward_run3 = (TH1F*)h2_trkIPresXY_PtReco_forward_run3->ProjectionY(TString::Format("trkIPresXY_PtReco_forward_run3_%d",i),lbin,ubin);
   proj_trkIPresXY_PtReco_forward_run3->Draw();
   proj_trkIPresXY_PtReco_forward_run3->Fit("gaus","","",-0.1,0.1);
   TF1 *fr5 = (TF1*)proj_trkIPresXY_PtReco_forward_run3->GetListOfFunctions()->FindObject("gaus");
   float sigma5 = fr5->GetParameter(2);
   float sigmaErr5 = fr5->GetParError(2);
   h1_pt_run3->GetXaxis()->SetRange(lbin,ubin);
   std::cout<<"ptRun3 forward : "<<pt<<" sigma : "<<sigma5<<" sigmaErr : "<<sigmaErr5<<std::endl;
   tg_trkIPresXY_PtReco_forward_run3->SetPoint(i,pt,sigma5);
   tg_trkIPresXY_PtReco_forward_run3->SetPointError(i,0.,sigmaErr5);
*/
   c6->cd(i+1);
   auto proj_trkIPresXY_PtReco_forward_run4 = (TH1F*)h2_trkIPresXY_PtReco_forward_run4->ProjectionY(TString::Format("trkIPresXY_PtReco_forward_run4_%d",i),lbin,ubin);
   proj_trkIPresXY_PtReco_forward_run4->Draw();
   proj_trkIPresXY_PtReco_forward_run4->Fit("gaus","","",-0.1,0.1);
   TF1 *fr6 = (TF1*)proj_trkIPresXY_PtReco_forward_run4->GetListOfFunctions()->FindObject("gaus");
   float sigma6 = fr6->GetParameter(2);
   float sigmaErr6 = fr6->GetParError(2);
   std::cout<<"ptRun4 forward : "<<pt2<<" sigma : "<<sigma6<<" sigmaErr : "<<sigmaErr6<<std::endl;
   tg_trkIPresXY_PtReco_forward_run4->SetPoint(i,pt2,sigma6);
   tg_trkIPresXY_PtReco_forward_run4->SetPointError(i,0.,sigmaErr6);
/*
   c7->cd(i+1);
   auto proj_trkIPresZ_PtReco_forward_run3 = (TH1F*)h2_trkIPresZ_PtReco_forward_run3->ProjectionY(TString::Format("trkIPresZ_PtReco_forward_run3_%d",i),lbin,ubin);
   proj_trkIPresZ_PtReco_forward_run3->Draw();
   proj_trkIPresZ_PtReco_forward_run3->Fit("gaus","","",-0.1,0.1);
   TF1 *fr7 = (TF1*)proj_trkIPresZ_PtReco_forward_run3->GetListOfFunctions()->FindObject("gaus");
   float sigma7 = fr7->GetParameter(2);
   float sigmaErr7 = fr7->GetParError(2);
   std::cout<<"ptRun3 forward : "<<pt<<" sigma : "<<sigma7<<" sigmaErr : "<<sigmaErr7<<std::endl;
   tg_trkIPresZ_PtReco_forward_run3->SetPoint(i,pt,sigma7);
   tg_trkIPresZ_PtReco_forward_run3->SetPointError(i,0.,sigmaErr7);
*/
   c8->cd(i+1);
   auto proj_trkIPresZ_PtReco_forward_run4 = (TH1F*)h2_trkIPresZ_PtReco_forward_run4->ProjectionY(TString::Format("trkIPresZ_PtReco_forward_run4_%d",i),lbin,ubin);
   proj_trkIPresZ_PtReco_forward_run4->Draw();
   proj_trkIPresZ_PtReco_forward_run4->Fit("gaus","","",-0.04,0.04);
   TF1 *fr8 = (TF1*)proj_trkIPresZ_PtReco_forward_run4->GetListOfFunctions()->FindObject("gaus");
   float sigma8 = fr8->GetParameter(2);
   float sigmaErr8 = fr8->GetParError(2);
   std::cout<<"ptRun4 forward : "<<pt2<<" sigma : "<<sigma8<<" sigmaErr : "<<sigmaErr8<<std::endl;
   tg_trkIPresZ_PtReco_forward_run4->SetPoint(i,pt2,sigma8);
   tg_trkIPresZ_PtReco_forward_run4->SetPointError(i,0.,sigmaErr8);

}

/*
auto c_tg_barrel = new TCanvas("c_tg_barrel","c_tg_barrel",750,500);
c_tg_barrel->cd();
c_tg_barrel->SetGrid();
tg_trkIPresXY_PtReco_barrel_run3->Draw("AplE");
tg_trkIPresXY_PtReco_barrel_run3->SetLineColor(4);
tg_trkIPresXY_PtReco_barrel_run3->SetMarkerColor(4);
tg_trkIPresXY_PtReco_barrel_run3->SetMarkerStyle(20);
tg_trkIPresXY_PtReco_barrel_run3->GetXaxis()->SetTitle("p_{T} (GeV)");
tg_trkIPresXY_PtReco_barrel_run3->GetYaxis()->SetTitle("Pointing Resolution (cm)");
tg_trkIPresXY_PtReco_barrel_run3->SetTitle();
tg_trkIPresXY_PtReco_barrel_run3->GetHistogram()->SetMaximum(0.025);   // along
tg_trkIPresXY_PtReco_barrel_run3->GetHistogram()->SetMinimum(0.0);  //   Y
tg_trkIPresXY_PtReco_barrel_run4->Draw("plEsame");
tg_trkIPresXY_PtReco_barrel_run4->SetLineColor(2);
tg_trkIPresXY_PtReco_barrel_run4->SetMarkerColor(2);
tg_trkIPresXY_PtReco_barrel_run4->SetMarkerStyle(21);
tg_trkIPresZ_PtReco_barrel_run3->Draw("plEsame");
tg_trkIPresZ_PtReco_barrel_run3->SetLineColor(4);
tg_trkIPresZ_PtReco_barrel_run3->SetMarkerColor(4);
tg_trkIPresZ_PtReco_barrel_run3->SetMarkerStyle(24);
tg_trkIPresZ_PtReco_barrel_run3->SetLineStyle(2);
tg_trkIPresZ_PtReco_barrel_run4->Draw("plEsame");
tg_trkIPresZ_PtReco_barrel_run4->SetLineColor(2);
tg_trkIPresZ_PtReco_barrel_run4->SetMarkerColor(2);
tg_trkIPresZ_PtReco_barrel_run4->SetMarkerStyle(25);
tg_trkIPresZ_PtReco_barrel_run4->SetLineStyle(2);
*/

auto c_tg_barrel = new TCanvas("c_tg_barrel","c_tg_barrel",750,500);
c_tg_barrel->cd();
c_tg_barrel->SetGrid();
tg_trkIPresXY_PtReco_barrel_run4->Draw("AplE");
tg_trkIPresXY_PtReco_barrel_run4->SetLineColor(4);
tg_trkIPresXY_PtReco_barrel_run4->SetMarkerColor(4);
tg_trkIPresXY_PtReco_barrel_run4->SetMarkerStyle(20);
tg_trkIPresXY_PtReco_barrel_run4->GetXaxis()->SetTitle("p_{T} (GeV)");
tg_trkIPresXY_PtReco_barrel_run4->GetYaxis()->SetTitle("Pointing Resolution (cm)");
tg_trkIPresXY_PtReco_barrel_run4->SetTitle();
tg_trkIPresXY_PtReco_barrel_run4->GetHistogram()->SetMaximum(0.025);   // along
tg_trkIPresXY_PtReco_barrel_run4->GetHistogram()->SetMinimum(0.0);  //   Y
tg_trkIPresZ_PtReco_barrel_run4->Draw("plEsame");
tg_trkIPresZ_PtReco_barrel_run4->SetLineColor(2);
tg_trkIPresZ_PtReco_barrel_run4->SetMarkerColor(2);
tg_trkIPresZ_PtReco_barrel_run4->SetMarkerStyle(21);
tg_trkIPresXY_PtReco_forward_run4->Draw("plEsame");
tg_trkIPresXY_PtReco_forward_run4->SetLineColor(4);
tg_trkIPresXY_PtReco_forward_run4->SetMarkerColor(4);
tg_trkIPresXY_PtReco_forward_run4->SetMarkerStyle(24);
tg_trkIPresXY_PtReco_forward_run4->SetLineStyle(2);
tg_trkIPresZ_PtReco_forward_run4->Draw("plEsame");
tg_trkIPresZ_PtReco_forward_run4->SetLineColor(2);
tg_trkIPresZ_PtReco_forward_run4->SetMarkerColor(2);
tg_trkIPresZ_PtReco_forward_run4->SetMarkerStyle(25);
tg_trkIPresZ_PtReco_forward_run4->SetLineStyle(2);


/*auto c_tg_forward = new TCanvas("c_tg_forward","c_tg_forward",750,500);
c_tg_forward->cd();
c_tg_forward->SetGrid();
tg_trkIPresXY_PtReco_forward_run3->Draw("AplE");
tg_trkIPresXY_PtReco_forward_run3->SetLineColor(4);
tg_trkIPresXY_PtReco_forward_run3->SetMarkerColor(4);
tg_trkIPresXY_PtReco_forward_run3->SetMarkerStyle(20);
tg_trkIPresXY_PtReco_forward_run3->GetXaxis()->SetTitle("p_{T} (GeV)");
tg_trkIPresXY_PtReco_forward_run3->GetYaxis()->SetTitle("Pointing Resolution (cm)");
tg_trkIPresXY_PtReco_forward_run3->SetTitle();
tg_trkIPresXY_PtReco_forward_run3->GetHistogram()->SetMaximum(0.025);   // along
tg_trkIPresXY_PtReco_forward_run3->GetHistogram()->SetMinimum(0.0);  //   Y
tg_trkIPresXY_PtReco_forward_run4->Draw("plEsame");
tg_trkIPresXY_PtReco_forward_run4->SetLineColor(2);
tg_trkIPresXY_PtReco_forward_run4->SetMarkerColor(2);
tg_trkIPresXY_PtReco_forward_run4->SetMarkerStyle(21);
tg_trkIPresZ_PtReco_forward_run3->Draw("plEsame");
tg_trkIPresZ_PtReco_forward_run3->SetLineColor(4);
tg_trkIPresZ_PtReco_forward_run3->SetMarkerColor(4);
tg_trkIPresZ_PtReco_forward_run3->SetMarkerStyle(24);
tg_trkIPresZ_PtReco_forward_run3->SetLineStyle(2);
tg_trkIPresZ_PtReco_forward_run4->Draw("plEsame");
tg_trkIPresZ_PtReco_forward_run4->SetLineColor(2);
tg_trkIPresZ_PtReco_forward_run4->SetMarkerColor(2);
tg_trkIPresZ_PtReco_forward_run4->SetMarkerStyle(25);
tg_trkIPresZ_PtReco_forward_run4->SetLineStyle(2);
*/


}
