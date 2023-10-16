#include <Math/Vector4D.h>

Double_t GetQ(ROOT::Math::PtEtaPhiMVector &p1, ROOT::Math::PtEtaPhiMVector &p2){
   ROOT::Math::PtEtaPhiMVector Sum4V = p1+p2;
   Double_t q = Sum4V.M2() - 4*1.865*1.865;
   return (  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}

void check_D0_femto_pPb(){

TFile* f = new TFile("/eos/cms/store/group/phys_heavyions/caber/dataMBAll_pPb.root");
TFile* fout = new TFile("out.root","recreate");

TTree* t = (TTree*) f->Get("d0ana_newreduced/VertexCompositeNtuple");

Int_t candSize;
Float_t pT[1000];
Float_t y[1000];
Float_t eta[1000];
Float_t phi[1000];
Float_t mass[1000];
Float_t mva[1000];
Float_t pTD1[1000];
Float_t pTD2[1000];

t->SetBranchAddress("candSize",&candSize);
t->SetBranchAddress("pT",&pT);
t->SetBranchAddress("y",&y);
t->SetBranchAddress("eta",&eta);
t->SetBranchAddress("phi",&phi);
t->SetBranchAddress("mass",&mass);
t->SetBranchAddress("mva",&mva);
t->SetBranchAddress("pTD1",&pTD1);
t->SetBranchAddress("pTD2",&pTD2);

//Histograms
TH1F *h_NselD0perEvent = new TH1F("NselD0perEvent","Selected D0/D0bar per event",10,0.5,10.5);
h_NselD0perEvent->SetLineColor(2);
h_NselD0perEvent->SetLineWidth(2);
h_NselD0perEvent->SetMarkerColor(2);
h_NselD0perEvent->GetXaxis()->SetTitle("Number of D0/D0bar / event");
h_NselD0perEvent->GetYaxis()->SetTitle("Number of D0/D0bar / bin");
TH1F *h_MASSselD0 = new TH1F("MASSselD0","Mass of selected D0/D0bar",50,1.75,2.00);
h_MASSselD0->SetLineColor(2);
h_MASSselD0->SetLineWidth(2);
h_MASSselD0->SetMarkerColor(2);
h_MASSselD0->GetXaxis()->SetTitle("Mass [GeV]");
h_MASSselD0->GetYaxis()->SetTitle("Number of D0/D0bar / bin");
TH1F *h_QselD0 = new TH1F("QselD0","qinv of selected pairs of D0/D0bar",50,0.0,2.0);
h_QselD0->SetLineColor(2);
h_QselD0->SetLineWidth(2);
h_QselD0->SetMarkerColor(2);
h_QselD0->GetXaxis()->SetTitle("qinv [GeV]");
h_QselD0->GetYaxis()->SetTitle("Number of D0/D0bar pairs / bin");

Int_t nevents = t->GetEntries();
///Int_t nevents = 100000;
//loop over all events
for (Int_t i=0;i<nevents;i++) {

   t->GetEntry(i);

   Int_t aux_nD0s = 0;
   //loop over all D0 candidates
   for (Int_t ii=0;ii<candSize;ii++) {
      if(std::abs(y[ii])>1) continue;
      if(mva[ii]<0.56) continue;
      h_MASSselD0->Fill(mass[ii]);
      if(!(1.84<mass[ii] && mass[ii]<1.89))continue;
      aux_nD0s++; 
      ROOT::Math::PtEtaPhiMVector aux1_D0FourVector;
      aux1_D0FourVector.SetM(1.865);
      aux1_D0FourVector.SetEta(eta[ii]);
      aux1_D0FourVector.SetPhi(phi[ii]);
      aux1_D0FourVector.SetPt(pT[ii]);   
       
      //to build qinv      
      for (Int_t iii=ii+1;iii<candSize;iii++) {
         if(std::abs(y[iii])>1) continue;
         if(mva[iii]<0.56) continue;
	 if(!(1.84<mass[iii] && mass[iii]<1.89))continue;
	 if(pTD1[ii]==pTD1[iii] || pTD2[ii]==pTD2[iii] || pTD1[ii]==pTD2[iii] || pTD2[ii]==pTD1[iii])continue;
         ROOT::Math::PtEtaPhiMVector aux2_D0FourVector;
	 aux2_D0FourVector.SetM(1.865);
         aux2_D0FourVector.SetEta(eta[iii]);
         aux2_D0FourVector.SetPhi(phi[iii]);
         aux2_D0FourVector.SetPt(pT[iii]);
	 Double_t aux_q = GetQ(aux1_D0FourVector, aux2_D0FourVector);
	 h_QselD0->Fill(aux_q);
      }	      
   }
   h_NselD0perEvent->Fill(aux_nD0s);

}

fout->Write();

}
