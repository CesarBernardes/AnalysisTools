// -*- C++ -*-
//
// Package:    NTuplizerHBT_GenOnly/NTuplizerHBT_GenOnly
// Class:      NTuplizerHBT_GenOnly
// 
/**\class NTuplizerHBT_GenOnly NTuplizerHBT_GenOnly.cc NTuplizerHBT_GenOnly/NTuplizerHBT_GenOnly/plugins/NTuplizerHBT_GenOnly.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Cesar Bernardes
//         Created:  Sat, 16 Jul 2016 09:10:17 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//new
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TROOT.h"
#include "TFile.h"
#include "TVector3.h"
#include <TRandom1.h>
#include <TLorentzVector.h>
#include "TProfile.h"
#include "Math/Point3D.h"
#include "THnSparse.h"

#include <map>

//
// class declaration
//

using namespace edm;
using namespace reco;
using namespace std;

class NTuplizerHBT_GenOnly : public edm::EDAnalyzer {
   public:
      explicit NTuplizerHBT_GenOnly(const edm::ParameterSet&);
      ~NTuplizerHBT_GenOnly();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      // define some constants
      Double_t pi_mass  = 0.1396;
      Double_t cos_cut = 0.99996;
      Double_t dpt_cut = 0.04;

      //to obtain Ntrkoffline in GenP
      TFile *f_trk_corr;
      TH2F *reff2D;
      TH2F *rsec2D;
      TH2F *rfak2D;
      TH2F *rmul2D;

      //functions
      static bool etaMixSort(std::pair<Double_t,std::vector<TLorentzVector>> & a, std::pair<Double_t,std::vector<TLorentzVector>> & b);
      static bool etaMixSort_ch(std::pair<Double_t,std::vector<int>> & a, std::pair<Double_t,std::vector<int>> & b);
      static bool etaMixSort_ntrkoff(std::pair<Double_t,int> & a, std::pair<Double_t,int> & b); 
      bool splitcomb(TLorentzVector &vec1,TLorentzVector &vec2);
      Double_t GetQ(const TLorentzVector &p1, const TLorentzVector &p2);
      Double_t GetQlong( const TLorentzVector &p1, const TLorentzVector &p2 );
      Double_t GetQlongLCMS( const TLorentzVector &p1, const TLorentzVector &p2 );
      Double_t GetQout( const TLorentzVector &p1, const TLorentzVector &p2 );
      Double_t GetQside( const TLorentzVector &p1, const TLorentzVector &p2 ); 
      const TLorentzVector InvertPVector(TLorentzVector &vec);
      const TLorentzVector InvertXYVector(TLorentzVector &vec);
      Double_t ComputeEventWeight(std::vector<TLorentzVector> GoodTrackFourVector);
      Double_t Encode( int w );
      void initHistos(const edm::Service<TFileService> & fs);
      double GetNtrkoffline(edm::Handle<reco::GenParticleCollection> genP,std::vector<TLorentzVector> &GoodTrackFourVector_gp_trkoff);     
      void MixEvents(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix);
      void MixEvents_eta(int ntrkoff_min, int ntrkoff_max); 
      double getTrkCorrWeight(double pT, double eta);     

      //histograms - Control plots
      ///events
      TH1F *gp_multiplicity_;
      TH1F *gp_ntrkoffline_;
      ///tracks
      TH1F *gp_Pt_;
      TH1F *gp_Eta_;
      TH1F *gp_Phi_;
      TH1F *gp_pairSS_M_;
      TH1F *gp_pairOS_M_;
      TH1F *gp_pairSS_mix_M_;
      TH1F *gp_pairOS_mix_M_;
      ///extra
      TH2F *gp_mixsch_KtVsNch;
      TH2F *gp_sig_KtVsNch_ntrkoff;
      TH2F *gp_ptVsNch_ntrkoff;
      TH2F *gp_etaVsNch_ntrkoff;
      TH2F *gp_phiVsNch_ntrkoff;
      TH2F *gp_sig_pairSS_MVsNch_ntrkoff;
      TH2F *gp_pairOS_MVsNch_ntrkoff;


      //THnSparse -- for correlation function 
      ///3D hist - 1D BEC analysis
      THnSparseF *gp_hs_qschVsKtVsNch_ntrkoff; //signal (qinv,kT,Ntrkoffline) -- MC: no Coulomb Correction
      THnSparseF *gp_hs_qdchVsKtVsNch_ntrkoff; //oppSigRef (qinv,kT,Ntrkoffline)
      THnSparseF *gp_hs_qINVschVsKtVsNch_ntrkoff; //InvertedRef (qinv,kT,Ntrkoffline)
      THnSparseF *gp_hs_qROTschVsKtVsNch_ntrkoff; //RotatedRef (qinv,kT,Ntrkoffline)
      THnSparseF *gp_hs_qinv_mixschVsKtVsNch_ntrkoff; //MixingRef Same-sign (qinv,kT,Ntrkoffline)
      THnSparseF *gp_hs_qinv_mixdchVsKtVsNch_ntrkoff; //MixingRef Opp-sign (qinv,kT,Ntrkoffline)
      ///5D hist - 3D BEC analysis
      THnSparseF *gp_hs_qLLCMSqOqSschVsKtVsNch_ntrkoff; //signal (qLLCMS,qOut,qSide,kT,Ntrkoffline) -- MC: no Coulomb Correction
      THnSparseF *gp_hs_qLLCMSqOqSdchVsKtVsNch_ntrkoff; //oppSigRef (qLLCMS,qOut,qSide,kT,Ntrkoffline) 
      THnSparseF *gp_hs_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff; //InvertedRef (qLLCMS,qOut,qSide,kT,Ntrkoffline) 
      THnSparseF *gp_hs_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff; //RotatedRef (qLLCMS,qOut,qSide,kT,Ntrkoffline)  
      THnSparseF *gp_hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff; //MixingRef Same-sign (qLLCMS,qOut,qSide,kT,Ntrkoffline) 
      THnSparseF *gp_hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff; //MixingRef Opp-sign (qLLCMS,qOut,qSide,kT,Ntrkoffline) 

  



      //tags for collections
      edm::EDGetTokenT<GenParticleCollection> genPToken_;

      //flags
      int Multiplicity_; //create histograms in bins of multiplicity 
      int mix_procedure_; //which mixing procedure to use. If equal to 1 is random and 2 is eta mix 

      //event information
      std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector_gp_vec;
      std::vector<std::vector<int>> ev_GoodTrackCharge_gp_vec;
      std::vector<double> ev_ntrkoff_vec;
      std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_gp_vec;
      std::vector<std::pair<Double_t, std::vector<int>> > ev_GoodTrackCharge_etaMixWeight_gp_vec;
      std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vec; 
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
NTuplizerHBT_GenOnly::NTuplizerHBT_GenOnly(const edm::ParameterSet& iConfig) : 
  genPToken_(consumes<GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genPTag"))){
   //now do what ever initialization is needed

  Multiplicity_       = iConfig.getParameter<uint32_t>("Multiplicity");
  mix_procedure_      = iConfig.getParameter<uint32_t>("Mix_procedure");

}


NTuplizerHBT_GenOnly::~NTuplizerHBT_GenOnly()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
NTuplizerHBT_GenOnly::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   std::vector<TLorentzVector> GoodTrackFourVector_gp;
   std::vector<int> GoodTrackCharge_gp; 
   std::vector<TLorentzVector> GoodTrackFourVector_gp_trkoff; //for eta mix

   edm::Handle<GenParticleCollection> genP;
   iEvent.getByToken(genPToken_,genP);

   //before doing loop in genParticles
   //get Ntrkoffline from the GenParticle multiplicity
   double Ntrkoffline = GetNtrkoffline(genP,GoodTrackFourVector_gp_trkoff);   
   //ev_ntrkoff_vec.push_back(Ntrkoffline);
   //std::cout<<"Ntrkoffline  :  "<<Ntrkoffline<<std::endl;


   if(Multiplicity_ == 0){if(80<=Ntrkoffline){return;}}//MB
   else if(Multiplicity_ == 1){if(80>Ntrkoffline){return;}}//HM

   //ev_ntrkoff_vec.push_back(Ntrkoffline); 

   // loop through "gen-only" particles
   for(size_t i_genPart=0; i_genPart<genP->size(); i_genPart++){

      const GenParticle & genPart = (*genP)[i_genPart];

      if(genPart.status() != 1 || genPart.charge()==0) continue; //only charged primaries
      if(genPart.pt()<0.2)continue;
      if(fabs(genPart.eta())>2.4)continue;

      gp_Pt_->Fill(genPart.pt());
      gp_Eta_->Fill(genPart.eta());
      gp_Phi_->Fill(genPart.phi());

      /*std::cout<<"genPart.status()  :  "<< genPart.status() << ";  genPart.charge()  :  " << genPart.charge() 
                   <<";  genPart.pdgId()  :  "<< genPart.pdgId()<<";  genPart.pt()  :   "<< genPart.pt() 
                                     <<";  genPart.eta()  :   "<<genPart.eta()
                                                       <<";  genPart.mass()  :  "<<genPart.mass()<<std::endl;*/
      TLorentzVector gp_vector;
      gp_vector.SetXYZM(genPart.px(),genPart.py(),genPart.pz(),pi_mass); //same as in reco
      GoodTrackFourVector_gp.push_back(gp_vector);
      GoodTrackCharge_gp.push_back(genPart.charge());

      TLorentzVector aux_gp_4vec; aux_gp_4vec.SetXYZM(genPart.px(),genPart.py(),genPart.pz(),pi_mass); //same as in reco 
      Double_t aux_pt = aux_gp_4vec.Pt();
      Double_t aux_eta = aux_gp_4vec.Eta();
      Double_t aux_phi = aux_gp_4vec.Phi();
      gp_ptVsNch_ntrkoff->Fill(aux_pt,Ntrkoffline);
      gp_etaVsNch_ntrkoff->Fill(aux_eta,Ntrkoffline);
      gp_phiVsNch_ntrkoff->Fill(aux_phi,Ntrkoffline);

      for(size_t ii_genPart=i_genPart+1; ii_genPart<genP->size(); ii_genPart++){

         const reco::GenParticle & genPart_assoc = (*genP)[ii_genPart];

         if(genPart_assoc.status() != 1 || genPart_assoc.charge()==0)continue; //only charged primaries      
         if(genPart_assoc.pt()<0.2)continue;
         if(fabs(genPart_assoc.eta())>2.4)continue;

         TLorentzVector aux_gp_vector;
         aux_gp_vector.SetXYZM(genPart_assoc.px(),genPart_assoc.py(),genPart_assoc.pz(),pi_mass); //same as in reco

         if(splitcomb(gp_vector,aux_gp_vector)){continue;}
          
         TLorentzVector gp_pair_vector = gp_vector + aux_gp_vector;
         Double_t q_gp = GetQ(gp_vector,aux_gp_vector);
         Double_t qo_gp= GetQ(gp_vector,InvertPVector(aux_gp_vector));
         Double_t qs_gp= GetQ(gp_vector,InvertXYVector(aux_gp_vector));
         //Double_t qlong_gp = GetQlong(GoodTrackFourVector[itk1],GoodTrackFourVector[itk2]);
         Double_t qlongLCMS_gp = GetQlongLCMS(gp_vector,aux_gp_vector);
         Double_t qout_gp = GetQout(gp_vector,aux_gp_vector);
         Double_t qside_gp = GetQside(gp_vector,aux_gp_vector);
         Double_t qlongLCMS_inv_gp = GetQlongLCMS(gp_vector,InvertPVector(aux_gp_vector));
         Double_t qlongLCMS_rot_gp = GetQlongLCMS(gp_vector,InvertXYVector(aux_gp_vector));
         Double_t qout_inv_gp = GetQout(gp_vector,InvertPVector(aux_gp_vector));
         Double_t qout_rot_gp = GetQout(gp_vector,InvertXYVector(aux_gp_vector));
         Double_t qside_inv_gp = GetQside(gp_vector,InvertPVector(aux_gp_vector));
         Double_t qside_rot_gp = GetQside(gp_vector,InvertXYVector(aux_gp_vector));

         Double_t kt_gp=(gp_pair_vector.Pt())/2.;
         TLorentzVector gp_pair_vector_rot = gp_vector + InvertXYVector(aux_gp_vector); 
         Double_t kt_gp_rot=(gp_pair_vector_rot.Pt())/2.;
         TLorentzVector gp_pair_vector_inv = gp_vector + InvertPVector(aux_gp_vector);
         Double_t kt_gp_inv = (gp_pair_vector_inv.Pt())/2.;

         Double_t x5D_LCMS_gp[5]={qlongLCMS_gp,qout_gp,qside_gp,kt_gp,(Double_t)Ntrkoffline};
         Double_t x5D_LCMS_inv_gp[5]={qlongLCMS_inv_gp,qout_inv_gp,qside_inv_gp,kt_gp_inv,(Double_t)Ntrkoffline};
         Double_t x5D_LCMS_rot_gp[5]={qlongLCMS_rot_gp,qout_rot_gp,qside_rot_gp,kt_gp_rot,(Double_t)Ntrkoffline};

         Double_t x3D_gp[3]={q_gp,kt_gp,(Double_t)Ntrkoffline};
         Double_t x3D_rot_gp[3]={qs_gp,kt_gp_rot,(Double_t)Ntrkoffline};
         Double_t x3D_inv_gp[3]={qo_gp,kt_gp_inv,(Double_t)Ntrkoffline};

         if( (genPart.charge())*(genPart_assoc.charge())>0 ){
            gp_pairSS_M_->Fill(gp_pair_vector.M());

            gp_sig_KtVsNch_ntrkoff->Fill(kt_gp,Ntrkoffline); 
            gp_sig_pairSS_MVsNch_ntrkoff->Fill(gp_pair_vector.M(),Ntrkoffline);

            //for 3D analysis
            gp_hs_qschVsKtVsNch_ntrkoff->Fill(x3D_gp);
            gp_hs_qINVschVsKtVsNch_ntrkoff->Fill(x3D_inv_gp);
            gp_hs_qROTschVsKtVsNch_ntrkoff->Fill(x3D_rot_gp);
            //for 5D analysis
            gp_hs_qLLCMSqOqSschVsKtVsNch_ntrkoff->Fill(x5D_LCMS_gp);
            gp_hs_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff->Fill(x5D_LCMS_inv_gp);
            gp_hs_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff->Fill(x5D_LCMS_rot_gp);

         }else{
            gp_pairOS_M_->Fill(gp_pair_vector.M());

            gp_pairOS_MVsNch_ntrkoff->Fill(gp_pair_vector.M(),Ntrkoffline);

            //for 3D analysis
            gp_hs_qdchVsKtVsNch_ntrkoff->Fill(x3D_gp);
            //for 5D analysis
            gp_hs_qLLCMSqOqSdchVsKtVsNch_ntrkoff->Fill(x5D_LCMS_gp); 

         }


      }


   }//first loop in genParticles

   ev_GoodTrackFourVector_gp_vec.push_back(GoodTrackFourVector_gp);
   ev_GoodTrackCharge_gp_vec.push_back(GoodTrackCharge_gp);
   ev_ntrkoff_vec.push_back(Ntrkoffline);

   //build vector of pairs to be sorted in etaMixWeight 
   Double_t aux_etaMix_w = ComputeEventWeight(GoodTrackFourVector_gp_trkoff);
   std::pair<Double_t, std::vector<TLorentzVector>> aux_pair_GoodTrackFourVector_etaMixWeight = make_pair(aux_etaMix_w, GoodTrackFourVector_gp);
   ev_GoodTrackFourVector_etaMixWeight_gp_vec.push_back(aux_pair_GoodTrackFourVector_etaMixWeight);
   std::pair<Double_t, std::vector<int>> aux_pair_GoodTrackCharge_etaMixWeight = make_pair(aux_etaMix_w,GoodTrackCharge_gp);
   ev_GoodTrackCharge_etaMixWeight_gp_vec.push_back(aux_pair_GoodTrackCharge_etaMixWeight);
   std::pair<Double_t,int> aux_pair_ntrkoff_etaMixWeight = make_pair(aux_etaMix_w,Ntrkoffline);
   ev_ntrkoff_etaMixWeight_vec.push_back(aux_pair_ntrkoff_etaMixWeight);

}


// ------------ method called once each job just before starting event loop  ------------
void 
NTuplizerHBT_GenOnly::beginJob()
{

edm::Service<TFileService> fs;
initHistos(fs); // Initilization of Histograms

f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_Pythia8.root"); //trkoffline

TH2F * corr2Deff = (TH2F*) f_trk_corr->Get(Form("%s/heff","HITrackCorrections"));
TH2F * corr2Dsim = (TH2F*) f_trk_corr->Get(Form("%s/hsim","HITrackCorrections"));
TH2F * corr2Drec = (TH2F*) f_trk_corr->Get(Form("%s/hrec","HITrackCorrections"));
TH2F * corr2Dfak = (TH2F*) f_trk_corr->Get(Form("%s/hfak","HITrackCorrections"));
TH2F * corr2Dsec = (TH2F*) f_trk_corr->Get(Form("%s/hsec","HITrackCorrections"));
TH2F * corr2Dmul = (TH2F*) f_trk_corr->Get(Form("%s/hmul","HITrackCorrections"));

reff2D = (TH2F*) corr2Deff->Clone("reff");
rmul2D = (TH2F*) corr2Dmul->Clone("rmul");
rfak2D = (TH2F*) corr2Dfak->Clone("rfak");
rsec2D = (TH2F*) corr2Dsec->Clone("rsec");

reff2D->Divide(corr2Deff,corr2Dsim,1,1,"B");
rmul2D->Divide(corr2Dmul,corr2Dsim,1,1,"B");
rfak2D->Divide(corr2Dfak,corr2Drec,1,1,"B");
rsec2D->Divide(corr2Dsec,corr2Drec,1,1,"B");


}

// ------------ method called once each job just after ending the event loop  ------------
void 
NTuplizerHBT_GenOnly::endJob() 
{

std::cout << " This is called at the end of the each job " << std::endl;
std::cout << " Now is the time for mixing the events "<< std::endl;

if(mix_procedure_==1){ //random mix

//default is to mix 40 events...used 25 and 100 for more systematics

if(Multiplicity_ == 0){
   MixEvents(0,9,100);
   MixEvents(10,19,100);
   MixEvents(20,29,100);
   MixEvents(30,39,100);
   MixEvents(40,49,100);
   MixEvents(50,59,100);
   MixEvents(60,69,100);
   MixEvents(70,79,100);  
}
else if(Multiplicity_ == 1){
   MixEvents(80,84,100);
   MixEvents(85,89,100);
   MixEvents(90,94,100);
   MixEvents(95,99,100);
   MixEvents(100,104,100);  
   MixEvents(105,109,100);
   MixEvents(110,114,100);
   MixEvents(115,119,100);
   MixEvents(120,124,100);
   MixEvents(125,129,100);
   MixEvents(130,134,100);
   MixEvents(135,139,100);
   MixEvents(140,144,100);
   MixEvents(145,149,100);
   MixEvents(150,154,100);
   MixEvents(155,159,100);
   MixEvents(160,164,100);
   MixEvents(165,169,100);
   MixEvents(170,174,100);
   MixEvents(175,179,100);
   MixEvents(180,250,100);
}
else{}


}//end "if" condition for random mix


//call function for eta mix
////Mix two-by-two events in the sorted vector by etaMixWeights for a given 
////ntrkoffline range
/////MixEvents_eta(ntrkoff min(inclusive), ntrkoff max(inclusive))
if(mix_procedure_==2){

if(Multiplicity_ == 0){
   MixEvents_eta(0,9);
   MixEvents_eta(10,19);
   MixEvents_eta(20,29);
   MixEvents_eta(30,39);
   MixEvents_eta(40,49);
   MixEvents_eta(50,59);
   MixEvents_eta(60,69);
   MixEvents_eta(70,79);  
}
else if(Multiplicity_ == 1){
   MixEvents_eta(80,84);
   MixEvents_eta(85,89);
   MixEvents_eta(90,94);
   MixEvents_eta(95,99);
   MixEvents_eta(100,104);  
   MixEvents_eta(105,109);
   MixEvents_eta(110,114);
   MixEvents_eta(115,119);
   MixEvents_eta(120,124);
   MixEvents_eta(125,129);
   MixEvents_eta(130,134);
   MixEvents_eta(135,139);
   MixEvents_eta(140,144);
   MixEvents_eta(145,149);
   MixEvents_eta(150,154);
   MixEvents_eta(155,159);
   MixEvents_eta(160,164);
   MixEvents_eta(165,169);
   MixEvents_eta(170,174);
   MixEvents_eta(175,179);
   MixEvents_eta(180,250);
}
else{}



}



}


bool
NTuplizerHBT_GenOnly::etaMixSort(std::pair<Double_t,std::vector<TLorentzVector>> & a, std::pair<Double_t,std::vector<TLorentzVector>> & b){

   return a.first < b.first ? true : false ;

}

bool  
NTuplizerHBT_GenOnly::etaMixSort_ch(std::pair<Double_t,std::vector<int>> & a, std::pair<Double_t,std::vector<int>> & b){

   return a.first < b.first ? true : false ;

}

bool
NTuplizerHBT_GenOnly::etaMixSort_ntrkoff(std::pair<Double_t,int> & a, std::pair<Double_t,int> & b){

   return a.first < b.first ? true : false ;

}

bool
NTuplizerHBT_GenOnly::splitcomb(TLorentzVector &vec1,TLorentzVector &vec2){
   bool issplit=false;
   Double_t cosa = TMath::Abs(vec1.Px()*vec2.Px() + vec1.Py()*vec2.Py() + vec1.Pz()*vec2.Pz())/(vec1.P()*vec2.P());
   Double_t deltapt = TMath::Abs(vec1.Pt() - vec2.Pt());
   if ( (cosa >cos_cut) && (deltapt < dpt_cut) ) { issplit = true;}
   return issplit;
}

Double_t
NTuplizerHBT_GenOnly::GetQ(const TLorentzVector &p1, const TLorentzVector &p2){
   TLorentzVector Sum4V = p1+p2;
   Double_t q = Sum4V.Mag2() - 4*pi_mass*pi_mass;
   return (  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}


Double_t
NTuplizerHBT_GenOnly::GetQlong(const TLorentzVector &p1, const TLorentzVector &p2){
   TLorentzVector Diff4V = p1-p2;
   Double_t qlong = fabs(Diff4V.Pz());
   return qlong;
}

Double_t
NTuplizerHBT_GenOnly::GetQlongLCMS(const TLorentzVector &p1, const TLorentzVector &p2){
   Double_t num = 2*( (p1.Pz())*(p2.E()) - (p2.Pz())*(p1.E()) );
   Double_t den = TMath::Sqrt( (p1.E()+p2.E())*(p1.E()+p2.E()) - (p1.Pz()+p2.Pz())*(p1.Pz()+p2.Pz()) );
   Double_t qlongLCMS = 0.0;
   if(den != 0) qlongLCMS = fabs(num/den);
   return qlongLCMS;
}

Double_t
NTuplizerHBT_GenOnly::GetQout(const TLorentzVector &p1, const TLorentzVector &p2){

   TVector3 qT;
   qT.SetXYZ(p1.Px()-p2.Px(),p1.Py()-p2.Py(),0.0);
   TVector3 kT;
   kT.SetXYZ( (p1.Px()+p2.Px())/2.0 , (p1.Py()+p2.Py())/2.0 ,0.0);
   TVector3 qout;
   qout = qT.Dot(kT.Unit())*kT.Unit();
   
   Double_t absValue = qout.Mag();

   return absValue; 


}

Double_t
NTuplizerHBT_GenOnly::GetQside(const TLorentzVector &p1, const TLorentzVector &p2){

   TVector3 qT;
   qT.SetXYZ(p1.Px()-p2.Px(),p1.Py()-p2.Py(),0.0);
   TVector3 kT;
   kT.SetXYZ( (p1.Px()+p2.Px())/2.0 , (p1.Py()+p2.Py())/2.0 ,0.0);
   TVector3 qout;
   qout = qT.Dot(kT.Unit())*kT.Unit();
   TVector3 qsid;
   qsid = qT - qout;

   Double_t absValue = qsid.Mag();

   return absValue;

}



const TLorentzVector
NTuplizerHBT_GenOnly::InvertPVector( TLorentzVector &vec){
   TLorentzVector ovec = vec;
   ovec.SetPx(-vec.Px());
   ovec.SetPy(-vec.Py());
   ovec.SetPz(-vec.Pz());
   return ovec;
}

const TLorentzVector
NTuplizerHBT_GenOnly::InvertXYVector( TLorentzVector &vec){
   TLorentzVector ovec = vec;
   ovec.SetX(-vec.X());
   ovec.SetY(-vec.Y());
   return ovec;
}

Double_t
NTuplizerHBT_GenOnly::ComputeEventWeight(std::vector<TLorentzVector> GoodTrackFourVector){
   int nFwd(0), nBkw(0), nCtr(0) ;
   for( unsigned int i=0 ; i< GoodTrackFourVector.size() ; i++ ){
      Double_t eta = GoodTrackFourVector[i].Eta();
      if( eta < -0.8 ){ nBkw++ ;}
      else if( eta > 0.8 ){ nFwd++ ;}
      else{ nCtr++ ;}
   }

   Double_t ReturnValue = (100*Encode( nFwd ) + 10 * Encode( nCtr ) + Encode( nBkw )) ;
   return ReturnValue ;
}

Double_t
NTuplizerHBT_GenOnly::Encode( int w ){
   int Range[]={3,5,8,10,13,16,20,25,30,200,10000} ;//added a protection for high-multiplicity 
   int i(0), j(0) ;
   while( w >= Range[i++] ) j++ ;
   return (Long64_t) j ;
}


double
NTuplizerHBT_GenOnly::getTrkCorrWeight(double pT, double eta){

  double eff = reff2D->GetBinContent(
                  reff2D->GetXaxis()->FindBin(eta),
                  reff2D->GetYaxis()->FindBin(pT) );
  if(eff >= 0.9999 || eff <= 0.0001) eff = 1;

  double sec = rsec2D->GetBinContent(
              rsec2D->GetXaxis()->FindBin(eta),
              rsec2D->GetYaxis()->FindBin(pT));
  if( sec >= 0.9999 || sec <= 0.0001) sec = 0;

  double fak = rfak2D->GetBinContent(
              rfak2D->GetXaxis()->FindBin(eta),
              rfak2D->GetYaxis()->FindBin(pT));
  if( fak >= 0.9999 || fak <= 0.0001) fak = 0;

  double mul = rmul2D->GetBinContent(
              rmul2D->GetXaxis()->FindBin(eta),
              rmul2D->GetYaxis()->FindBin(pT));
  if( mul >= 0.9999 || mul <= 0.0001) mul = 0;


  return (1. - fak ) * ( 1. - sec ) / eff  / (1. + mul );
  //return (1. - fak ) / eff;
  //return 1. / eff;
}


double
NTuplizerHBT_GenOnly::GetNtrkoffline(edm::Handle<reco::GenParticleCollection> genP,std::vector<TLorentzVector> &GoodTrackFourVector_gp_trkoff){


double aux_ntrkoffline=0.0;
int genP_multiplicity=0;

for(size_t i_genPart=0; i_genPart<genP->size(); i_genPart++){

      const GenParticle & genPart = (*genP)[i_genPart];

      if(genPart.status() != 1 || genPart.charge()==0) continue; //only charged primaries
      if(genPart.pt()<0.4)continue; //used in trackoffline selection
      if(fabs(genPart.eta())>2.4)continue;

      genP_multiplicity++;

      double aux_tk_corr = getTrkCorrWeight(genPart.pt(), genPart.eta()); 
      if(aux_tk_corr != 0) aux_ntrkoffline = aux_ntrkoffline + 1.0*(1.0/aux_tk_corr);
      else aux_ntrkoffline++;

      TLorentzVector pvector;
      pvector.SetXYZM(genPart.px(),genPart.py(),genPart.pz(),pi_mass);
      GoodTrackFourVector_gp_trkoff.push_back(pvector);  
}


gp_multiplicity_->Fill(genP_multiplicity);
gp_ntrkoffline_->Fill(aux_ntrkoffline);

return aux_ntrkoffline;

}


void NTuplizerHBT_GenOnly::initHistos(const edm::Service<TFileService> & fs){


TH1D::SetDefaultSumw2();
TFileDirectory trkQA = fs->mkdir( "TrackHistograms" );

gp_Pt_ = trkQA.make<TH1F>("gp_Pt_", "GenParticles: p_{T} spectrum", 100, 0, 20);
gp_Pt_->GetXaxis()->SetTitle("p_{T} [GeV]");
double etaBins[13] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0,
        0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
gp_Eta_ = trkQA.make<TH1F>("gp_Eta_", "GenParticles: Pseudorapidity distribution",12,etaBins); //same used in trk corrections
gp_Eta_->GetXaxis()->SetTitle("#eta distribution of the tracks");
gp_Phi_ = trkQA.make<TH1F>("gp_Phi_", "GenParticles: Azimuthal distribution",100, -3.15, 3.15);
gp_Phi_->GetXaxis()->SetTitle("#phi distribution of the tracks");

gp_pairSS_M_ = trkQA.make<TH1F>("gp_pairSS_M_", "GenParticles: invariant mass same-sign tracks", 100, 0, 1);
gp_pairSS_M_->GetXaxis()->SetTitle("M [GeV]");
gp_pairOS_M_ = trkQA.make<TH1F>("gp_pairOS_M_", "GenParticles: invariant mass opposite-sign tracks", 100, 0, 1);
gp_pairOS_M_->GetXaxis()->SetTitle("M [GeV]");

gp_multiplicity_ = trkQA.make<TH1F>("gp_multiplicity_","Primary and charged Gen Particle multiplicity",250,-0.5,249.5);
gp_multiplicity_->GetXaxis()->SetTitle("Multiplicity");
gp_ntrkoffline_  = trkQA.make<TH1F>("gp_ntrkoffline_","Ntrkoffline",250,-0.5,249.5);
gp_ntrkoffline_->GetXaxis()->SetTitle("Multiplicity");


TFileDirectory mixedEvents = fs->mkdir( "MixedHistograms" );

gp_pairSS_mix_M_ = mixedEvents.make<TH1F>("gp_pairSS_mix_M_", "GenParticles: invariant mass same-sign tracks after mix", 100, 0, 1);
gp_pairSS_mix_M_->GetXaxis()->SetTitle("M [GeV]");
gp_pairOS_mix_M_ = mixedEvents.make<TH1F>("gp_pairOS_mix_M_", "GenParticles: invariant mass opposite-sign tracks after mix", 100, 0, 1);
gp_pairOS_mix_M_->GetXaxis()->SetTitle("M [GeV]");

gp_mixsch_KtVsNch= mixedEvents.make<TH2F>("gp_mixsch_KtVsNch","gp_mixsch_KtVsNch",200,0.,20.,250,-0.5,249.5);

//THnSparse
Int_t bins5D[5]=   {40, 40, 40, 10, 25  };
Double_t xmin5D[5]={0., 0., 0., 0., 0.  };
Double_t xmax5D[5]={2., 2., 2., 1., 250.};

Int_t bins3D[3]=   {200,20,250  };
Double_t xmin3D[3]={0. ,0.,-0.5 };
Double_t xmax3D[3]={2. ,1.,249.5};

gp_hs_qinv_mixschVsKtVsNch_ntrkoff = mixedEvents.make<THnSparseF>("gp_hs_qinv_mixschVsKtVsNch_ntrkoff","gp_hs_qinv_mixschVsKtVsNch_ntrkoff",3, bins3D,  xmin3D, xmax3D);
gp_hs_qinv_mixdchVsKtVsNch_ntrkoff = mixedEvents.make<THnSparseF>("gp_hs_qinv_mixdchVsKtVsNch_ntrkoff","gp_hs_qinv_mixdchVsKtVsNch_ntrkoff",3, bins3D,  xmin3D, xmax3D);
gp_hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff = mixedEvents.make<THnSparseF>("gp_hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff","gp_hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
gp_hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff = mixedEvents.make<THnSparseF>("gp_hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff","gp_hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
 


TFileDirectory qinvHist = fs->mkdir( "qinvHistograms" );

gp_sig_KtVsNch_ntrkoff= qinvHist.make<TH2F>("gp_sig_KtVsNch_ntrkoff","gp_sig_KtVsNch_ntrkoff",200,0.,20.,250,-0.5,249.5);
gp_ptVsNch_ntrkoff = qinvHist.make<TH2F>("gp_ptVsNch_ntrkoff","gp_ptVsNch_ntrkoff",200,0.,20.,250,-0.5,249.5);
gp_etaVsNch_ntrkoff = qinvHist.make<TH2F>("gp_etaVsNch_ntrkoff","gp_etaVsNch_ntrkoff",12,etaBins,250,-0.5,249.5); //same used in trk corrections
gp_phiVsNch_ntrkoff = qinvHist.make<TH2F>("gp_phiVsNch_ntrkoff","gp_phiVsNch_ntrkoff",100, -3.15, 3.15,250,-0.5,249.5);
gp_sig_pairSS_MVsNch_ntrkoff= qinvHist.make<TH2F>("gp_sig_pairSS_MVsNch_ntrkoff","gp_sig_pairSS_MVsNch_ntrkoff",200,0,2,250,-0.5,249.5);
gp_pairOS_MVsNch_ntrkoff= qinvHist.make<TH2F>("gp_pairOS_MVsNch_ntrkoff","gp_pairOS_MVsNch_ntrkoff",200,0,2,250,-0.5,249.5);


//THnSparse
gp_hs_qLLCMSqOqSschVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("gp_hs_qLLCMSqOqSschVsKtVsNch_ntrkoff","gp_hs_qLLCMSqOqSschVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
gp_hs_qLLCMSqOqSdchVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("gp_hs_qLLCMSqOqSdchVsKtVsNch_ntrkoff","gp_hs_qLLCMSqOqSdchVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
gp_hs_qschVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("gp_hs_qschVsKtVsNch_ntrkoff","gp_hs_qschVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
gp_hs_qINVschVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("gp_hs_qINVschVsKtVsNch_ntrkoff","gp_hs_qINVschVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
gp_hs_qROTschVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("gp_hs_qROTschVsKtVsNch_ntrkoff","gp_hs_qROTschVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
gp_hs_qdchVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("gp_hs_qdchVsKtVsNch_ntrkoff","gp_hs_qdchVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
gp_hs_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("gp_hs_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff","gp_hs_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
gp_hs_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("gp_hs_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff","gp_hs_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);



}

void 
NTuplizerHBT_GenOnly::MixEvents(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix){

int aux_n_evts = (int)ev_ntrkoff_vec.size();
//std::cout<<"Number of events for mixing  :  " << aux_n_evts  << std::endl;


for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {

   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_gp_vec= ev_GoodTrackFourVector_gp_vec[nevt1];
   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_gp_vec.size();
   std::vector<int> ev_GoodTrackChargeTemp_nevt1_gp_vec = ev_GoodTrackCharge_gp_vec[nevt1];

   if(ev_ntrkoff_vec[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff_vec[nevt1])continue;

   int takeAssociated = 0;
   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
      if(ev_ntrkoff_vec[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff_vec[nevt_assoc])continue;

      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_gp_vec= ev_GoodTrackFourVector_gp_vec[nevt_assoc];
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_gp_vec.size();
      std::vector<int> ev_GoodTrackChargeTemp_nevt_assoc_gp_vec=ev_GoodTrackCharge_gp_vec[nevt_assoc];

      takeAssociated++;
      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing.

      for(int imix=0; imix<nMixmult_nevt1; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
            if(splitcomb(ev_GoodTrackFourVectorTemp_nevt1_gp_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_gp_vec[iimix])){continue;}
            Double_t qmix_gp = GetQ(ev_GoodTrackFourVectorTemp_nevt1_gp_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_gp_vec[iimix]);
            //Double_t qmixlong = GetQlong(ev_GoodTrackFourVectorTemp_nevt1_gp_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_gp_vec[iimix]);
            Double_t qmixlongLCMS_gp = GetQlongLCMS(ev_GoodTrackFourVectorTemp_nevt1_gp_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_gp_vec[iimix]);
            Double_t qmixout_gp = GetQout(ev_GoodTrackFourVectorTemp_nevt1_gp_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_gp_vec[iimix]);
            Double_t qmixside_gp = GetQside(ev_GoodTrackFourVectorTemp_nevt1_gp_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_gp_vec[iimix]);


            TLorentzVector psum2_gp = ev_GoodTrackFourVectorTemp_nevt1_gp_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_gp_vec[iimix];
            Double_t kt_gp=(psum2_gp.Pt())/2.;

            Double_t xmix5D_LCMS_gp[5]={qmixlongLCMS_gp,qmixout_gp,qmixside_gp,kt_gp,(Double_t)ev_ntrkoff_vec[nevt1]};
            Double_t xmix3D_gp[3]={qmix_gp,kt_gp,(Double_t)ev_ntrkoff_vec[nevt1]};

           
            if(ev_GoodTrackChargeTemp_nevt1_gp_vec[imix]*ev_GoodTrackChargeTemp_nevt_assoc_gp_vec[iimix]>0){
               gp_pairSS_mix_M_->Fill(psum2_gp.M());

               gp_mixsch_KtVsNch->Fill(kt_gp,ev_ntrkoff_vec[nevt1]); 

               gp_hs_qinv_mixschVsKtVsNch_ntrkoff->Fill(xmix3D_gp);
               gp_hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff->Fill(xmix5D_LCMS_gp);

            }
            else{
               gp_pairOS_mix_M_->Fill(psum2_gp.M());

               gp_hs_qinv_mixdchVsKtVsNch_ntrkoff->Fill(xmix3D_gp);
               gp_hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff->Fill(xmix5D_LCMS_gp);

            }
         }
      }

   } 
   
}

}//end of "MixEvents" function


void 
NTuplizerHBT_GenOnly::MixEvents_eta(int ntrkoff_min, int ntrkoff_max){


///sort vectors of maps for each ntrkoffline range
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > aux_ev_GoodTrackFourVector_etaMixWeight_vec;
std::vector<std::pair<Double_t, std::vector<int>> > aux_ev_GoodTrackCharge_etaMixWeight_vec;
std::vector<std::pair<Double_t,int> > aux_ev_ntrkoff_etaMixWeight_vec;


int aux_n_evts = ev_GoodTrackFourVector_etaMixWeight_gp_vec.size();

for(int ievt=0; ievt<aux_n_evts; ievt++){

   if((ev_ntrkoff_etaMixWeight_vec[ievt]).second<ntrkoff_min || ntrkoff_max<(ev_ntrkoff_etaMixWeight_vec[ievt]).second)continue; //only events in a given range

   aux_ev_GoodTrackFourVector_etaMixWeight_vec.push_back(ev_GoodTrackFourVector_etaMixWeight_gp_vec[ievt]);
   aux_ev_GoodTrackCharge_etaMixWeight_vec.push_back(ev_GoodTrackCharge_etaMixWeight_gp_vec[ievt]);
   aux_ev_ntrkoff_etaMixWeight_vec.push_back(ev_ntrkoff_etaMixWeight_vec[ievt]); 

}

std::sort( aux_ev_GoodTrackFourVector_etaMixWeight_vec.begin(), aux_ev_GoodTrackFourVector_etaMixWeight_vec.end(), NTuplizerHBT_GenOnly::etaMixSort );
std::sort( aux_ev_GoodTrackCharge_etaMixWeight_vec.begin(), aux_ev_GoodTrackCharge_etaMixWeight_vec.end(), NTuplizerHBT_GenOnly::etaMixSort_ch );
std::sort( aux_ev_ntrkoff_etaMixWeight_vec.begin(), aux_ev_ntrkoff_etaMixWeight_vec.end(), NTuplizerHBT_GenOnly::etaMixSort_ntrkoff );


int aux_n_evts_inNtrkoffRange = aux_ev_GoodTrackFourVector_etaMixWeight_vec.size();

//auxiliar vectors for mixing
std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec;
int nMixmult_nevt1;
std::vector<int> ev_GoodTrackChargeTemp_nevt1_vec;

std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec;
int nMixmult_nevt_assoc;
std::vector<int> ev_GoodTrackChargeTemp_nevt_assoc_vec;

//loop in the vector of maps for mixing - only in a given range of ntrkoffline
for(int nevt=0; nevt+1<aux_n_evts_inNtrkoffRange; nevt+=2) {


   ev_GoodTrackFourVectorTemp_nevt1_vec= (aux_ev_GoodTrackFourVector_etaMixWeight_vec[nevt]).second;
   nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   ev_GoodTrackChargeTemp_nevt1_vec = (aux_ev_GoodTrackCharge_etaMixWeight_vec[nevt]).second;


   ev_GoodTrackFourVectorTemp_nevt_assoc_vec= (aux_ev_GoodTrackFourVector_etaMixWeight_vec[nevt+1]).second;
   nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
   ev_GoodTrackChargeTemp_nevt_assoc_vec= (aux_ev_GoodTrackCharge_etaMixWeight_vec[nevt+1]).second;

   
   for(int imix=0; imix<nMixmult_nevt1; imix++){
      for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
         if(splitcomb(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix])){continue;}
         Double_t qmix_gp = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
         //Double_t qmixlong = GetQlong(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
         Double_t qmixlongLCMS_gp = GetQlongLCMS(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
         Double_t qmixout_gp = GetQout(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
         Double_t qmixside_gp = GetQside(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]); 

         TLorentzVector psum2_gp = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
         Double_t kt_gp=(psum2_gp.Pt())/2.;

         Double_t xmix3D_gp[3]={qmix_gp,kt_gp,(Double_t)(aux_ev_ntrkoff_etaMixWeight_vec[nevt]).second};
         Double_t xmix5D_LCMS_gp[5]={qmixlongLCMS_gp,qmixout_gp,qmixside_gp,kt_gp,(Double_t)(aux_ev_ntrkoff_etaMixWeight_vec[nevt]).second};

         if(ev_GoodTrackChargeTemp_nevt1_vec[imix]*ev_GoodTrackChargeTemp_nevt_assoc_vec[iimix]>0){
            gp_pairSS_mix_M_->Fill(psum2_gp.M());
            gp_mixsch_KtVsNch->Fill(kt_gp,(aux_ev_ntrkoff_etaMixWeight_vec[nevt]).second);

            gp_hs_qinv_mixschVsKtVsNch_ntrkoff->Fill(xmix3D_gp);
            gp_hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff->Fill(xmix5D_LCMS_gp);

         }else{

            gp_pairOS_mix_M_->Fill(psum2_gp.M());

            gp_hs_qinv_mixdchVsKtVsNch_ntrkoff->Fill(xmix3D_gp);
            gp_hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff->Fill(xmix5D_LCMS_gp);  

         }
      }
   }



}//end first for loop

//clear vectors for next ntrkoff range
aux_ev_GoodTrackFourVector_etaMixWeight_vec.clear();
aux_ev_GoodTrackCharge_etaMixWeight_vec.clear();
aux_ev_ntrkoff_etaMixWeight_vec.clear();


}//end MixEvents_eta function






// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NTuplizerHBT_GenOnly::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NTuplizerHBT_GenOnly);
