// -*- C++ -*-
//
// Package:    NTuplizerHBT/NtuplizerHBT
// Class:      NtuplizerHBT
// 
/**\class NtuplizerHBT NtuplizerHBT.cc NTuplizerHBT/NtuplizerHBT/plugins/NtuplizerHBT.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Cesar Bernardes
//         Created:  Thu, 10 Dec 2015 15:53:07 GMT
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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TROOT.h"
#include "TFile.h"
#include "TVector3.h"
#include <TRandom1.h>
#include <TLorentzVector.h>
#include "TProfile.h"
#include "Math/Point3D.h"

#include <map>

#include "NtuplizerSelections.h"


//
// class declaration
//

class NtuplizerHBT : public edm::EDAnalyzer {
   public:
      explicit NtuplizerHBT(const edm::ParameterSet&);
      ~NtuplizerHBT();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

      //define some constants
      Double_t pi_mass  = 0.1396;
      Double_t cos_cut = 0.99996;
      Double_t dpt_cut = 0.04; 

      //functions
      void initialize();        
      static bool vtxSort( const reco::Vertex &  a, const reco::Vertex & b );
      static bool etaMixSort(std::pair<Double_t,std::vector<TLorentzVector>> & a, std::pair<Double_t,std::vector<TLorentzVector>> & b);
      static bool etaMixSort_ch(std::pair<Double_t,std::vector<int>> & a, std::pair<Double_t,std::vector<int>> & b);
      static bool etaMixSort_ntrkoff(std::pair<Double_t,int> & a, std::pair<Double_t,int> & b); 
      bool splitcomb(TLorentzVector &vec1,TLorentzVector &vec2);
      Double_t GetQ( const TLorentzVector &p1, const TLorentzVector &p2 );
      Double_t GetQlong( const TLorentzVector &p1, const TLorentzVector &p2 );
      Double_t GetQlongLCMS( const TLorentzVector &p1, const TLorentzVector &p2 );
      Double_t GetQout( const TLorentzVector &p1, const TLorentzVector &p2 );
      Double_t GetQside( const TLorentzVector &p1, const TLorentzVector &p2 ); 
      const TLorentzVector InvertPVector( TLorentzVector &vec);
      const TLorentzVector InvertXYVector( TLorentzVector &vec); 
      const Double_t CoulombWpm( const Double_t& q);
      const Double_t CoulombW( const Double_t& q);
      Double_t ComputeEventWeight(std::vector<TLorentzVector> GoodTrackFourVector); 
      Double_t Encode( int w );
      int ntrkoff(edm::Handle<std::vector<reco::Track> > generalTracks,double vtx_x,double vtx_y,double vtx_z,std::vector<TLorentzVector> &GoodTrackFourVector_trkoff);
      double getTrkCorrWeight(double pT, double eta);
      void fillTreeWithCommonVariables(math::XYZPoint BSPosition, math::XYZPoint vtx,
                                       double vtx_xError, double vtx_yError, double vtx_zError,
                                       std::vector<reco::TrackCollection::const_iterator> selected_tracks);
      void fillHistosWithCommonVariables(math::XYZPoint BSPosition, math::XYZPoint vtx,
                                       double vtx_xError, double vtx_yError, double vtx_zError,
                                       std::vector<reco::TrackCollection::const_iterator> selected_tracks); 
      void initHistos(const edm::Service<TFileService> & fs);
      void MixEvents(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix);
      void MixEvents_eta(int ntrkoff_min, int ntrkoff_max);

      //main tree
      TTree *my_tree;

      //histograms - Control plots
      ///event
      TH1F *Selections_;
      TH1F *VertexX_;
      TH1F *VertexY_;
      TH1F *VertexZ_;
      TH1F *NVertex_;
      TH1F *N_tk_HPTracks_; //high-purity
      TH1F *N_tk_AllTracks_;//no selection
      TH1F *N_tk_SelTracks_;//selected tracks
      TH1F *N_tk_Offline_;//selected track-offline 
      ///tracks
      TH1F *tk_NHits_;
      TH1F *tk_Pt_;
      TH1F *tk_Eta_;
      TH1F *tk_Phi_;
      TH1F *tk_DzOverSigmaDz_;
      TH1F *tk_DxyOverSigmaDxy_;
      TH1F *tk_SigmaPtOverPt_;
      TH1F *tk_pixelLayersWithMeasurement_;   
      TH1F *tk_NormalizedChi2_;
      TH1F *tk_Ndof_; 
      TH1F *tk_offlineSelCorr_;  
      TH2F *tk_NtrkOfflineCorr_NtrkOffline_;
      TH1F *tk_pairSS_M_;
      TH1F *tk_pairOS_M_;
      ///Extra
      TH1F *DeltazVertex_; //between two events after mixing
      TH2F *mixsch_KtVsNch_ntrkoff;
      TH2F *sig_KtVsNch_ntrkoff;
      TH2F *ptVsNch_ntrkoff;
      TH2F *etaVsNch_ntrkoff;
      TH2F *phiVsNch_ntrkoff;
      TH2F *sig_pairSS_MVsNch_ntrkoff;
      TH2F *pairOS_MVsNch_ntrkoff;       

      
      //THnSparse -- for correlation function 
      ///3D hist - 1D BEC analysis -- (qinv,kT,Ntrkoffline)
      THnSparseF *hs_qschVsKtVsNch_ntrkoff; //signal -- with coulomb correction (coulombCorr)
      THnSparseF *hs_qschncVsKtVsNch_ntrkoff; //signal -- without coulomb correction (NoCoulombCorr) 
      THnSparseF *hs_qdchVsKtVsNch_ntrkoff; //oppSigRef (coulombCorr)
      THnSparseF *hs_qdchncVsKtVsNch_ntrkoff; //oppSigRef (NoCoulombCorr)
      THnSparseF *hs_qINVschVsKtVsNch_ntrkoff; //InvertedRef (coulombCorr)
      THnSparseF *hs_qINVschncVsKtVsNch_ntrkoff; //InvertedRef (NoCoulombCorr)
      THnSparseF *hs_qROTschVsKtVsNch_ntrkoff; //RotatedRef (coulombCorr)
      THnSparseF *hs_qROTschncVsKtVsNch_ntrkoff; //RotatedRef (NoCoulombCorr) 
      THnSparseF *hs_qinv_mixschVsKtVsNch_ntrkoff; //MixingRef Same-sign (coulombCorr)
      THnSparseF *hs_qinv_mixdchVsKtVsNch_ntrkoff; //MixingRef Opp-sign (NoCoulombCorr)   
      ///5D hist - 3D BEC analysis -- (qLLCMS,qOut,qSide,kT,Ntrkoffline)
      THnSparseF *hs_qLLCMSqOqSschVsKtVsNch_ntrkoff; //signal -- with coulomb correction (coulombCorr)
      THnSparseF *hs_qLLCMSqOqSschncVsKtVsNch_ntrkoff; //signal -- without coulomb correction (NoCoulombCorr) 
      THnSparseF *hs_qLLCMSqOqSdchVsKtVsNch_ntrkoff;  //oppSigRef (coulombCorr)
      THnSparseF *hs_qLLCMSqOqSdchncVsKtVsNch_ntrkoff;  //oppSigRef (NoCoulombCorr)
      THnSparseF *hs_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff; //InvertedRef (coulombCorr)
      THnSparseF *hs_qLLCMSqOqS_INV_schncVsKtVsNch_ntrkoff; //InvertedRef (NoCoulombCorr)
      THnSparseF *hs_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff; //RotatedRef (coulombCorr)
      THnSparseF *hs_qLLCMSqOqS_ROT_schncVsKtVsNch_ntrkoff; //RotatedRef (NoCoulombCorr)
      THnSparseF *hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff; //MixingRef Same-sign (coulombCorr)
      THnSparseF *hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff; //MixingRef Opp-sign (NoCoulombCorr) 

 
      //tags for collections 
      edm::EDGetTokenT<reco::TrackCollection> trackToken_;
      edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
      edm::EDGetTokenT<TrackingParticleCollection> tpSrc_;
  
      //flags
      bool isRealData; //run on data or mc - chose in config file
      int selection;   //your track selection - chose in config file
      std::string Selection_; //your track selection - chose in config file
      bool fillTree_; //If true and running on data it will fill trees of the order of TB!!!!
      int Multiplicity_; //create histograms in bins of multiplicity according to the trigger.
      int mix_procedure_; //which mixing procedure to use. If equal to 1 is random and 2 is eta mix

      //for track corrections
      //TrackCorrector2D testload;
      TFile *f_trk_corr;
      TH2F *reff2D;
      TH2F *rsec2D;
      TH2F *rfak2D;
      TH2F *rmul2D;

      //event information
      int evNumber, runNumber, LumiSection;
      double bs_x, bs_y, bs_z; 
      int N_vtx;
      double vtx_x, vtx_y, vtx_z, vtx_xError, vtx_yError, vtx_zError;  //primary vertex
      std::vector<double> ev_vtx_z_vec;
      std::vector<std::vector<TLorentzVector>> ev_GoodTrackFourVector_vec;
      std::vector<std::vector<int>> ev_GoodTrackCharge_vec;
      std::vector<int> ev_ntrkoff_vec; 
      std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > ev_GoodTrackFourVector_etaMixWeight_vec;
      std::vector<std::pair<Double_t, std::vector<int>> > ev_GoodTrackCharge_etaMixWeight_vec;
      std::vector<std::pair<Double_t,int> > ev_ntrkoff_etaMixWeight_vec;
 
      //tracks after selection - common variables to all selections
      int N_tk;
      std::vector<int> *tk_charge;
      std::vector<int> *tk_algo;  
      std::vector<double> *tk_pt; std::vector<double> *tk_ptError;
      std::vector<double> *tk_p; std::vector<double> *tk_px; std::vector<double> *tk_py; std::vector<double> *tk_pz;
      std::vector<double> *tk_eta; std::vector<double> *tk_etaError;
      std::vector<double> *tk_phi; std::vector<double> *tk_phiError; 
      std::vector<double> *tk_dxy; std::vector<double> *tk_dxy_bs; std::vector<double> *tk_dxy_vtx;
      std::vector<double> *tk_dxyError; std::vector<double> *tk_dxyError_vtx;
      std::vector<double> *tk_dz; std::vector<double> *tk_dz_bs; std::vector<double> *tk_dz_vtx;       
      std::vector<double> *tk_dzError; std::vector<double> *tk_dzError_vtx;
      std::vector<double> *tk_chi2; std::vector<double> *tk_ndof; std::vector<double> *tk_normalizedChi2;
      std::vector<int> *tk_numberOfValidHits; std::vector<int> *tk_numberOfLostHits; 

      //tracks after selection - specific selection for HBT cuts
      std::vector<double> *tk_inner_r;
      
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
NtuplizerHBT::NtuplizerHBT(const edm::ParameterSet& iConfig) : 
  trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("generalTrackTag"))),
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("offlinePrimaryVerticesTag"))),
  tpSrc_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("tpSrc"))){
  //now do what ever initialization is needed
  isRealData          = iConfig.getParameter<bool>("RealData");
  Selection_          = iConfig.getParameter<std::string>("Selection");
  fillTree_           = iConfig.getParameter<bool>("FillTree");
  Multiplicity_       = iConfig.getParameter<uint32_t>("Multiplicity");
  mix_procedure_      = iConfig.getParameter<uint32_t>("Mix_procedure"); 


  //prepare for switch
  if(Selection_ == "Baseline") selection = Baseline;
  else if(Selection_ == "HBT_Padova") selection = HBT_Padova;
  else if(Selection_ == "Ridge") selection = Ridge;
  else if(Selection_ == "HBT_Sprace") selection = HBT_Sprace;
  else {
     cms::Exception ex("InvalidConfiguration");
     ex << "Unknown selection " << Selection_  
        << ". Please check NtuplizerHBT.cc for allowed values.";
     throw ex; 
  }


}


NtuplizerHBT::~NtuplizerHBT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
NtuplizerHBT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;   

   Selections_->Fill(selection);

   //initialize variables and clean vectors
   initialize();

   //event ID
   runNumber    = iEvent.id().run();
   evNumber     = iEvent.id().event();
   LumiSection  = iEvent.id().luminosityBlock();

   //beamspot position
   math::XYZPoint BSPosition(0.0,0.0,0.0);
   edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot",recoBeamSpotHandle);
   reco::BeamSpot bs = *recoBeamSpotHandle;
   BSPosition = bs.position();
   //std::cout<<" bs x : "<<BSPosition.x()<<" bs y : "<<BSPosition.y()<<" bs z : "<<BSPosition.z()<<std::endl;
   bs_x = (double)BSPosition.x();
   bs_y = (double)BSPosition.y();
   bs_z = (double)BSPosition.z();

   //vertex collection - chose it in the config. file
   edm::Handle<std::vector<reco::Vertex> > vertexCollection;
   iEvent.getByToken(vertexToken_,vertexCollection);
   std::vector<reco::Vertex> vtx_sorted = *vertexCollection;
   std::sort( vtx_sorted.begin(), vtx_sorted.end(), NtuplizerHBT::vtxSort );
   if(vtx_sorted.size() == 0) return;
   if(fabs(vtx_sorted.begin()->position().z()) > 15.0) return; //default

   vtx_x = (double)vtx_sorted.begin()->position().x(); VertexX_->Fill(vtx_x);
   vtx_y = (double)vtx_sorted.begin()->position().y(); VertexY_->Fill(vtx_y);
   vtx_z = (double)vtx_sorted.begin()->position().z(); VertexZ_->Fill(vtx_z);
   vtx_xError = (double)vtx_sorted.begin()->xError();
   vtx_yError = (double)vtx_sorted.begin()->yError();
   vtx_zError = (double)vtx_sorted.begin()->zError();
   N_vtx = vtx_sorted.size(); NVertex_->Fill(N_vtx);

   //track collection - chose it in the config. file
   edm::Handle<reco::TrackCollection> generalTracks;
   iEvent.getByToken(trackToken_,generalTracks);
   //std::cout<<"Size of General Tracks   :  "<< generalTracks->size() << std::endl;

   //vector for selected tracks - to fill Tree and Histograms
   std::vector<reco::TrackCollection::const_iterator> selected_tracks;
   //vectors to compute Qinv in the same event
   std::vector<TLorentzVector> GoodTrackFourVector;
   std::vector<int> GoodTrackCharge;
   std::vector<TLorentzVector> GoodTrackFourVector_trkoff; //for eta mix
   //define the selections and assign values to the tracks' variables 
   //the switch will chose the selection you asked in the config file
   int aux_N_tk_highPurity = 0; //Number of HighPurityTracks
   switch(selection){
      case Baseline : {
         math::XYZPoint vtx(vtx_x,vtx_y,vtx_z);
         for(reco::TrackCollection::const_iterator iter_tk = generalTracks->begin(); iter_tk != generalTracks->end(); iter_tk++){
            if(iter_tk->pt()<0.2)continue;
            if(fabs(iter_tk->eta())>2.4)continue;
            if(!iter_tk->quality(reco::TrackBase::highPurity))continue;
            aux_N_tk_highPurity++;
            selected_tracks.push_back(iter_tk);
            TLorentzVector pvector;
            pvector.SetXYZM(iter_tk->px(),iter_tk->py(),iter_tk->pz(),pi_mass);
            GoodTrackFourVector.push_back(pvector);
            GoodTrackCharge.push_back(iter_tk->charge()); 
         }
         fillTreeWithCommonVariables(BSPosition,vtx,vtx_xError,vtx_yError,vtx_zError,selected_tracks);
         fillHistosWithCommonVariables(BSPosition,vtx,vtx_xError,vtx_yError,vtx_zError,selected_tracks);
         break;
      }  
      case HBT_Padova : {
         math::XYZPoint vtx(vtx_x,vtx_y,vtx_z);
         for(reco::TrackCollection::const_iterator iter_tk = generalTracks->begin(); iter_tk != generalTracks->end(); iter_tk++){
            if(iter_tk->pt()<0.2)continue;
            if(fabs(iter_tk->eta())>2.4)continue;
            if(!iter_tk->quality(reco::TrackBase::highPurity))continue;
            aux_N_tk_highPurity++;
            if(iter_tk->normalizedChi2()>=5.0)continue;
            if(iter_tk->ndof()<=5)continue;
            double aux_tk_AbsDxy_vtx = (double)fabs(iter_tk->dxy(vtx));
            double aux_tk_inner_r = (double)sqrt(iter_tk->innerPosition().x()*iter_tk->innerPosition().x() + iter_tk->innerPosition().y()*iter_tk->innerPosition().y());
            if(aux_tk_AbsDxy_vtx >= 0.15 || aux_tk_inner_r >= 20.0) continue;
            selected_tracks.push_back(iter_tk);
            tk_inner_r->push_back(aux_tk_inner_r);
            TLorentzVector pvector;
            pvector.SetXYZM(iter_tk->px(),iter_tk->py(),iter_tk->pz(),pi_mass);
            GoodTrackFourVector.push_back(pvector);
            GoodTrackCharge.push_back(iter_tk->charge());
         }
         fillTreeWithCommonVariables(BSPosition,vtx,vtx_xError,vtx_yError,vtx_zError,selected_tracks);
         fillHistosWithCommonVariables(BSPosition,vtx,vtx_xError,vtx_yError,vtx_zError,selected_tracks);
         break;
      }
      case Ridge : {
         math::XYZPoint vtx(vtx_x,vtx_y,vtx_z);
         for(reco::TrackCollection::const_iterator iter_tk = generalTracks->begin(); iter_tk != generalTracks->end(); iter_tk++){
            double aux_tk_dz_vtx = (double)iter_tk->dz(vtx);
            double aux_tk_dzError_vtx  = (double)sqrt(iter_tk->dzError()*iter_tk->dzError()+vtx_zError*vtx_zError);
            double aux_tk_dxy_vtx = (double)iter_tk->dxy(vtx);
            double aux_tk_dxyError_vtx  = (double)sqrt(iter_tk->dxyError()*iter_tk->dxyError()+vtx_xError*vtx_yError); 
            if(iter_tk->pt()<0.2)continue;
            if(fabs(iter_tk->eta())>2.4)continue;
            if(!iter_tk->quality(reco::TrackBase::highPurity))continue;
            aux_N_tk_highPurity++;
            if(fabs(iter_tk->ptError())/iter_tk->pt()>0.1)continue;
            if(fabs(aux_tk_dz_vtx/aux_tk_dzError_vtx)>3)continue;
            if(fabs(aux_tk_dxy_vtx/aux_tk_dxyError_vtx)>3)continue; 
            selected_tracks.push_back(iter_tk);
            TLorentzVector pvector;
            pvector.SetXYZM(iter_tk->px(),iter_tk->py(),iter_tk->pz(),pi_mass);
            GoodTrackFourVector.push_back(pvector);
            GoodTrackCharge.push_back(iter_tk->charge()); 
         }
         fillTreeWithCommonVariables(BSPosition,vtx,vtx_xError,vtx_yError,vtx_zError,selected_tracks);
         fillHistosWithCommonVariables(BSPosition,vtx,vtx_xError,vtx_yError,vtx_zError,selected_tracks);
         break;
      }
      case HBT_Sprace : {
         math::XYZPoint vtx(vtx_x,vtx_y,vtx_z);
         for(reco::TrackCollection::const_iterator iter_tk = generalTracks->begin(); iter_tk != generalTracks->end(); iter_tk++){
            double aux_tk_dz_vtx = (double)iter_tk->dz(vtx);
            double aux_tk_dzError_vtx  = (double)sqrt(iter_tk->dzError()*iter_tk->dzError()+vtx_zError*vtx_zError);
            double aux_tk_dxy_vtx = (double)iter_tk->dxy(vtx);
            double aux_tk_dxyError_vtx  = (double)sqrt(iter_tk->dxyError()*iter_tk->dxyError()+vtx_xError*vtx_yError);
            const reco::HitPattern& hit_pattern = iter_tk->hitPattern(); 
            if(iter_tk->pt()<0.2)continue;
            if(fabs(iter_tk->eta())>2.4)continue;
            if(!iter_tk->quality(reco::TrackBase::highPurity))continue;
            aux_N_tk_highPurity++;
            if(fabs(iter_tk->ptError())/iter_tk->pt()>0.1)continue;
            if(fabs(aux_tk_dz_vtx/aux_tk_dzError_vtx)>3)continue;
            if(fabs(aux_tk_dxy_vtx/aux_tk_dxyError_vtx)>3)continue;
            if(hit_pattern.pixelLayersWithMeasurement()==0)continue; 
            selected_tracks.push_back(iter_tk);
            TLorentzVector pvector;
            pvector.SetXYZM(iter_tk->px(),iter_tk->py(),iter_tk->pz(),pi_mass);
            GoodTrackFourVector.push_back(pvector);
            GoodTrackCharge.push_back(iter_tk->charge());
         }
         fillTreeWithCommonVariables(BSPosition,vtx,vtx_xError,vtx_yError,vtx_zError,selected_tracks);
         fillHistosWithCommonVariables(BSPosition,vtx,vtx_xError,vtx_yError,vtx_zError,selected_tracks);
         break;      
      }
   }//end of switch
   //std::cout<<"Size of Selected Tracks   :  "<< selected_tracks.size() << std::endl;
   N_tk_HPTracks_->Fill(aux_N_tk_highPurity); 
   int aux_N_tk_generalTracks = generalTracks->size();
   N_tk_AllTracks_->Fill(aux_N_tk_generalTracks);
   int aux_N_tk_selected = selected_tracks.size();
   N_tk_SelTracks_->Fill(aux_N_tk_selected);
   int aux_N_tk_offline = ntrkoff(generalTracks,vtx_x,vtx_y,vtx_z,GoodTrackFourVector_trkoff);
   N_tk_Offline_->Fill(aux_N_tk_offline);

   if(Multiplicity_ == 0){if(80<=aux_N_tk_offline){return;}}
   ///if(Multiplicity_ == 0){if(250<aux_N_tk_offline){return;}} //for extended MinBias 
   else if(Multiplicity_ == 1){if(80>aux_N_tk_offline || 105<=aux_N_tk_offline){return;}}
   else if(Multiplicity_ == 2){if(105>aux_N_tk_offline || 130<=aux_N_tk_offline){return;}}
   else if(Multiplicity_ == 3){if(130>aux_N_tk_offline){return;}}
   else{} 

   if(GoodTrackFourVector.size()<2)return; //event not used for signal, then do not use for mixing

   for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
      Double_t aux_tk1_corr = (Double_t)getTrkCorrWeight(GoodTrackFourVector[itk1].Pt(),GoodTrackFourVector[itk1].Eta());
      Double_t aux_pt = GoodTrackFourVector[itk1].Pt();
      Double_t aux_eta = GoodTrackFourVector[itk1].Eta(); 
      Double_t aux_phi = GoodTrackFourVector[itk1].Phi(); 
      ptVsNch_ntrkoff->Fill(aux_pt,aux_N_tk_offline,aux_tk1_corr);
      etaVsNch_ntrkoff->Fill(aux_eta,aux_N_tk_offline,aux_tk1_corr);
      phiVsNch_ntrkoff->Fill(aux_phi,aux_N_tk_offline,aux_tk1_corr);
      for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
         if(splitcomb(GoodTrackFourVector[itk1],  GoodTrackFourVector[itk2])){continue;}
         Double_t q = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);
         Double_t qo= GetQ(GoodTrackFourVector[itk1],InvertPVector(GoodTrackFourVector[itk2]));
         Double_t qs= GetQ(GoodTrackFourVector[itk1],InvertXYVector(GoodTrackFourVector[itk2]));
         //Double_t qlong = GetQlong(GoodTrackFourVector[itk1],GoodTrackFourVector[itk2]);
         Double_t qlongLCMS = GetQlongLCMS(GoodTrackFourVector[itk1],GoodTrackFourVector[itk2]);
         Double_t qout = GetQout(GoodTrackFourVector[itk1],GoodTrackFourVector[itk2]);
         Double_t qside = GetQside(GoodTrackFourVector[itk1],GoodTrackFourVector[itk2]); 
         Double_t qlongLCMS_inv = GetQlongLCMS(GoodTrackFourVector[itk1],InvertPVector(GoodTrackFourVector[itk2]));
         Double_t qlongLCMS_rot = GetQlongLCMS(GoodTrackFourVector[itk1],InvertXYVector(GoodTrackFourVector[itk2]));
         Double_t qout_inv = GetQout(GoodTrackFourVector[itk1],InvertPVector(GoodTrackFourVector[itk2]));
         Double_t qout_rot = GetQout(GoodTrackFourVector[itk1],InvertXYVector(GoodTrackFourVector[itk2]));
         Double_t qside_inv = GetQside(GoodTrackFourVector[itk1],InvertPVector(GoodTrackFourVector[itk2]));
         Double_t qside_rot = GetQside(GoodTrackFourVector[itk1],InvertXYVector(GoodTrackFourVector[itk2])); 

         TLorentzVector psum2 = GoodTrackFourVector[itk1] + GoodTrackFourVector[itk2];
         Double_t kt=(psum2.Pt())/2.;
         TLorentzVector psum2_rot = GoodTrackFourVector[itk1] + InvertXYVector(GoodTrackFourVector[itk2]);
         Double_t kt_rot = (psum2_rot.Pt())/2.;
         TLorentzVector psum2_inv = GoodTrackFourVector[itk1] + InvertPVector(GoodTrackFourVector[itk2]);
         Double_t kt_inv = (psum2_inv.Pt())/2.;

         Double_t x5D_LCMS[5]={qlongLCMS,qout,qside,kt,(Double_t)aux_N_tk_offline};
         Double_t x5D_LCMS_inv[5]={qlongLCMS_inv,qout_inv,qside_inv,kt_inv,(Double_t)aux_N_tk_offline};
         Double_t x5D_LCMS_rot[5]={qlongLCMS_rot,qout_rot,qside_rot,kt_rot,(Double_t)aux_N_tk_offline};
         
         Double_t x3D[3]={q,kt,(Double_t)aux_N_tk_offline};
         Double_t x3D_rot[3]={qs,kt_rot,(Double_t)aux_N_tk_offline};
         Double_t x3D_inv[3]={qo,kt_inv,(Double_t)aux_N_tk_offline};

         Double_t aux_tk2_corr = (Double_t)getTrkCorrWeight(GoodTrackFourVector[itk2].Pt(),GoodTrackFourVector[itk2].Eta()); 
         Double_t aux_tk12_corr= aux_tk1_corr*aux_tk2_corr;
         //aux_tk12_corr=1.0; ///no corrections. Just to compare with corrections

         if(GoodTrackCharge[itk1]*GoodTrackCharge[itk2]>0){
            tk_pairSS_M_->Fill(psum2.M(),aux_tk12_corr);           

            sig_KtVsNch_ntrkoff->Fill(kt,aux_N_tk_offline,aux_tk12_corr);
            sig_pairSS_MVsNch_ntrkoff->Fill(psum2.M(),aux_N_tk_offline,aux_tk12_corr);

            //for 3D analysis
            hs_qschVsKtVsNch_ntrkoff->Fill(x3D,CoulombW(q)*aux_tk12_corr);
            hs_qschncVsKtVsNch_ntrkoff->Fill(x3D,aux_tk12_corr);
            hs_qINVschVsKtVsNch_ntrkoff->Fill(x3D_inv,CoulombW(q)*aux_tk12_corr);
            hs_qINVschncVsKtVsNch_ntrkoff->Fill(x3D_inv,aux_tk12_corr);
            hs_qROTschVsKtVsNch_ntrkoff->Fill(x3D_rot,CoulombW(q)*aux_tk12_corr);
            hs_qROTschncVsKtVsNch_ntrkoff->Fill(x3D_rot,aux_tk12_corr);
            //for 5D analysis
            hs_qLLCMSqOqSschVsKtVsNch_ntrkoff->Fill(x5D_LCMS,CoulombW(q)*aux_tk12_corr);
            hs_qLLCMSqOqSschncVsKtVsNch_ntrkoff->Fill(x5D_LCMS,aux_tk12_corr);
            hs_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff->Fill(x5D_LCMS_inv,CoulombW(q)*aux_tk12_corr);
            hs_qLLCMSqOqS_INV_schncVsKtVsNch_ntrkoff->Fill(x5D_LCMS_inv,aux_tk12_corr);
            hs_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff->Fill(x5D_LCMS_rot,CoulombW(q)*aux_tk12_corr);
            hs_qLLCMSqOqS_ROT_schncVsKtVsNch_ntrkoff->Fill(x5D_LCMS_rot,aux_tk12_corr);
 
         }else{
            //now, the same for opposite charge
            tk_pairOS_M_->Fill(psum2.M(),aux_tk12_corr);
            pairOS_MVsNch_ntrkoff->Fill(psum2.M(),aux_N_tk_offline,aux_tk12_corr);

            //for 3D analysis
            hs_qdchVsKtVsNch_ntrkoff->Fill(x3D,CoulombWpm(q)*aux_tk12_corr);
            hs_qdchncVsKtVsNch_ntrkoff->Fill(x3D,aux_tk12_corr); 
            //for 5D analysis
            hs_qLLCMSqOqSdchVsKtVsNch_ntrkoff->Fill(x5D_LCMS,CoulombWpm(q)*aux_tk12_corr); 
            hs_qLLCMSqOqSdchncVsKtVsNch_ntrkoff->Fill(x5D_LCMS,aux_tk12_corr);

         }
      }
   }
   ev_vtx_z_vec.push_back(vtx_z);
   ev_GoodTrackFourVector_vec.push_back(GoodTrackFourVector);
   ev_GoodTrackCharge_vec.push_back(GoodTrackCharge);
   ev_ntrkoff_vec.push_back(aux_N_tk_offline);

   //build vector of pairs ordered by the etaMixWeight
   Double_t aux_etaMix_w = ComputeEventWeight(GoodTrackFourVector_trkoff);
   //std::cout<<"Event Eta Mix weight Using trkoff selection  :  "<<aux_etaMix_w<<std::endl;
   std::pair<Double_t, std::vector<TLorentzVector>> aux_pair_GoodTrackFourVector_etaMixWeight = make_pair(aux_etaMix_w, GoodTrackFourVector);
   ev_GoodTrackFourVector_etaMixWeight_vec.push_back(aux_pair_GoodTrackFourVector_etaMixWeight);
   std::pair<Double_t, std::vector<int>> aux_pair_GoodTrackCharge_etaMixWeight = make_pair(aux_etaMix_w,GoodTrackCharge); 
   ev_GoodTrackCharge_etaMixWeight_vec.push_back(aux_pair_GoodTrackCharge_etaMixWeight);
   std::pair<Double_t,int> aux_pair_ntrkoff_etaMixWeight = make_pair(aux_etaMix_w,aux_N_tk_offline);
   ev_ntrkoff_etaMixWeight_vec.push_back(aux_pair_ntrkoff_etaMixWeight);



   if(fillTree_){my_tree->Fill();}
   else{}

}

void
NtuplizerHBT::initialize()
{

// event information
evNumber      = -9999;
runNumber     = -9999;
LumiSection   = -9999;
/// beamspot position
bs_x = -9999.;
bs_y = -9999.;
bs_z = -9999.;
/// vertex
////number of vertices
N_vtx = 0;
////position of primary vertex - highest multiplicity
vtx_x = -9999.;
vtx_y = -9999.;
vtx_z = -9999.;
vtx_xError = -9999.;
vtx_yError = -9999.;
vtx_zError = -9999.;
/// tracks after selection
N_tk = 0;
tk_charge->clear();
tk_algo->clear();
tk_pt->clear();
tk_ptError->clear();
tk_p->clear();
tk_px->clear();
tk_py->clear();
tk_pz->clear();
tk_eta->clear();
tk_etaError->clear();
tk_phi->clear();
tk_phiError->clear();
tk_dxy->clear();
tk_dxy_bs->clear();
tk_dxy_vtx->clear();
tk_dxyError->clear();
tk_dxyError_vtx->clear();
tk_dz->clear();
tk_dz_bs->clear();
tk_dz_vtx->clear();
tk_dzError->clear();
tk_dzError_vtx->clear();
tk_chi2->clear();
tk_ndof->clear();
tk_normalizedChi2->clear();
tk_numberOfValidHits->clear();
tk_numberOfLostHits->clear();

tk_inner_r->clear();

}


// ------------ method called once each job just before starting event loop  ------------
void 
NtuplizerHBT::beginJob()
{

//for track corrections
//f_trk_corr = new TFile("/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_ntrkOffline_Pythia8/results/trkCorr_forHBTAnalysis_tot_Pythia8.root");
//f_trk_corr = new TFile("/afs/cern.ch/work/c/caber/HBT_13TeVAnalysis/pp13TeV_withSkimmedSamples/CMSSW_7_4_15_patch1/src/TrackingCode/HIRun2015Ana/test/TrkCorrHBTAnalysis_Projects/crab_trkCorrHBTAnalysis_pp13TeV_ntrkOffline_EPOS/results/trkCorr_forHBTAnalysis_tot_EPOS.root");

//f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_EPOS.root"); //to correct Ntrkoffline
//f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_Pythia8.root"); //to correct Ntrkoffline
//f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_V2.root"); //new SPRACE cuts - default
//f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_EPOS.root"); //new SPRACE cuts - EPOS for systematics
f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia6.root"); //new SPRACE cuts - pythia6 for systematics
//f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Tight.root"); //new SPRACE cuts - tight track selection for systematics
//f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Medium.root"); //new SPRACE cuts - medium track selection for systematics
//f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Loose.root"); //new SPRACE cuts - loose track selection for systematics
//f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Zvtx3cm.root"); //new SPRACE cuts - Zvtx3cm for systematics
//f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_3Zvtx15cm.root"); //new SPRACE cuts - 3Zvtx15cm for systematics

//f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Loose2.root"); //new SPRACE cuts - loose2 track selection for systematics
//f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Medium2.root"); //new SPRACE cuts - medium2 track selection for systematics
//f_trk_corr = new TFile("trkCorr_forHBTAnalysis_tot_SpraceMinusChi2NdofCuts_Pythia8_Tight2.root"); //new SPRACE cuts - tight2 track selection for systematics


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


edm::Service<TFileService> fs;
my_tree = fs->make<TTree>("track_tree","track_tree");
initHistos(fs); // Initilization of Histograms

tk_charge= new std::vector<int>;
tk_algo= new std::vector<int>;
tk_pt= new std::vector<double>;
tk_ptError= new std::vector<double>;
tk_p= new std::vector<double>;
tk_px= new std::vector<double>;
tk_py= new std::vector<double>;
tk_pz= new std::vector<double>;
tk_eta= new std::vector<double>;
tk_etaError= new std::vector<double>;
tk_phi= new std::vector<double>;
tk_phiError= new std::vector<double>;
tk_dxy= new std::vector<double>;
tk_dxy_bs= new std::vector<double>;
tk_dxy_vtx= new std::vector<double>;
tk_dxyError= new std::vector<double>;
tk_dxyError_vtx= new std::vector<double>;
tk_dz= new std::vector<double>;
tk_dz_bs= new std::vector<double>;
tk_dz_vtx= new std::vector<double>;
tk_dzError= new std::vector<double>;
tk_dzError_vtx= new std::vector<double>;
tk_chi2= new std::vector<double>;
tk_ndof= new std::vector<double>;
tk_normalizedChi2= new std::vector<double>;
tk_numberOfValidHits= new std::vector<int>;
tk_numberOfLostHits= new std::vector<int>;
tk_inner_r= new std::vector<double>;


my_tree->Branch("evNumber",    &evNumber,    "evNumber/I");
my_tree->Branch("runNumber",   &runNumber,   "runNumber/I");
my_tree->Branch("LumiSection", &LumiSection, "LumiSection/I");

my_tree->Branch("bs_x",    &bs_x,    "bs_x/D");
my_tree->Branch("bs_y",    &bs_y,    "bs_y/D");
my_tree->Branch("bs_z",    &bs_z,    "bs_z/D");

my_tree->Branch("N_vtx",    &N_vtx,    "N_vtx/I");
my_tree->Branch("vtx_x",    &vtx_x,    "vtx_x/D");
my_tree->Branch("vtx_y",    &vtx_y,    "vtx_y/D");
my_tree->Branch("vtx_z",    &vtx_z,    "vtx_z/D");
my_tree->Branch("vtx_xError",    &vtx_xError,    "vtx_xError/D");
my_tree->Branch("vtx_yError",    &vtx_yError,    "vtx_yError/D");
my_tree->Branch("vtx_zError",    &vtx_zError,    "vtx_zError/D");

my_tree->Branch("selection",    &selection,    "selection/I");

my_tree->Branch("N_tk",    &N_tk,    "N_tk/I");
my_tree->Branch("tk_charge", "std::vector<int>", &      tk_charge );
my_tree->Branch("tk_algo", "std::vector<int>", &      tk_algo );
my_tree->Branch("tk_pt", "std::vector<double>", &     tk_pt );
my_tree->Branch("tk_ptError", "std::vector<double>", &     tk_ptError );
my_tree->Branch("tk_p", "std::vector<double>", &     tk_p );
my_tree->Branch("tk_px", "std::vector<double>", &     tk_px );
my_tree->Branch("tk_py", "std::vector<double>", &     tk_py );
my_tree->Branch("tk_pz", "std::vector<double>", &     tk_pz );
my_tree->Branch("tk_eta", "std::vector<double>", &     tk_eta );
my_tree->Branch("tk_etaError", "std::vector<double>", &     tk_etaError );
my_tree->Branch("tk_phi", "std::vector<double>", &     tk_phi );
my_tree->Branch("tk_phiError", "std::vector<double>", &     tk_phiError );
my_tree->Branch("tk_dxy", "std::vector<double>", &     tk_dxy );
my_tree->Branch("tk_dxy_bs", "std::vector<double>", &     tk_dxy_bs );
my_tree->Branch("tk_dxy_vtx", "std::vector<double>", &     tk_dxy_vtx );
my_tree->Branch("tk_dxyError", "std::vector<double>", &     tk_dxyError );
my_tree->Branch("tk_dxyError_vtx", "std::vector<double>", &     tk_dxyError_vtx );
my_tree->Branch("tk_dz", "std::vector<double>", &     tk_dz );
my_tree->Branch("tk_dz_bs", "std::vector<double>", &     tk_dz_bs );
my_tree->Branch("tk_dz_vtx", "std::vector<double>", &     tk_dz_vtx );
my_tree->Branch("tk_dzError", "std::vector<double>", &     tk_dzError );
my_tree->Branch("tk_dzError_vtx", "std::vector<double>", &     tk_dzError_vtx );
my_tree->Branch("tk_chi2", "std::vector<double>", &     tk_chi2 );
my_tree->Branch("tk_ndof", "std::vector<double>", &     tk_ndof );
my_tree->Branch("tk_normalizedChi2", "std::vector<double>", &     tk_normalizedChi2 );
my_tree->Branch("tk_numberOfValidHits", "std::vector<int>", &      tk_numberOfValidHits );
my_tree->Branch("tk_numberOfLostHits", "std::vector<int>", &      tk_numberOfLostHits );

my_tree->Branch("tk_inner_r", "std::vector<double>", &     tk_inner_r );

}

// ------------ method called once each job just after ending the event loop  ------------
void 
NtuplizerHBT::endJob(){

std::cout << " This is called at the end of the each job " << std::endl;
std::cout << " Now is the time for mixing the events "<< std::endl;
//int aux_n_evts = (int)ev_ntrkoff_vec.size();
//std::cout<<"Number of events for mixing  :  " << aux_n_evts  << std::endl;

///MixEvents(ntrkoff min(inclusive), ntrkoff max(inclusive), max events to mix with a given event)

if(mix_procedure_==1){

if(Multiplicity_ == 0){
   MixEvents(0,9,40);
   MixEvents(10,19,40);
   MixEvents(20,29,40);
   MixEvents(30,39,40);
   MixEvents(40,49,40);
   MixEvents(50,59,40);
   MixEvents(60,69,40);
   MixEvents(70,79,40);  
}
/*if(Multiplicity_ == 0){ ////for extended MinBias 
   MixEvents(0,4,40);
   MixEvents(5,9,40);
   MixEvents(10,12,40);
   MixEvents(13,15,40);
   MixEvents(16,18,40);
   MixEvents(19,21,40);
   MixEvents(22,24,40);
   MixEvents(25,27,40);
   MixEvents(28,30,40);
   MixEvents(31,33,40);
   MixEvents(34,36,40);
   MixEvents(37,39,40);   
   MixEvents(40,44,40);
   MixEvents(45,49,40);
   MixEvents(50,54,40);
   MixEvents(55,59,40);
   MixEvents(60,69,40);
   MixEvents(70,79,40);
   MixEvents(80,89,40);
   MixEvents(90,99,40);
   MixEvents(100,109,40);
   MixEvents(110,119,40);
   MixEvents(120,250,40);
}*/
else if(Multiplicity_ == 1){
   MixEvents(80,84,40);
   MixEvents(85,89,40);
   MixEvents(90,94,40);
   MixEvents(95,99,40);
   MixEvents(100,104,40);  
}
else if(Multiplicity_ == 2){
   MixEvents(105,109,40);
   MixEvents(110,114,40);
   MixEvents(115,119,40);
   MixEvents(120,124,40);
   MixEvents(125,129,40);
}
else if(Multiplicity_ == 3){
   MixEvents(130,134,40);
   MixEvents(135,139,40);
   MixEvents(140,144,40);
   MixEvents(145,149,40);
   MixEvents(150,154,40);
   MixEvents(155,159,40);
   MixEvents(160,164,40);
   MixEvents(165,169,40);
   MixEvents(170,174,40);
   MixEvents(175,179,40);
   MixEvents(180,250,40);
}
else{}

}//end "if" condition for random mix


//call function for eta mix
//Mix two-by-two events in the sorted vector by etaMixWeights for a given 
//ntrkoffline range
///MixEvents_eta(ntrkoff min(inclusive), ntrkoff max(inclusive))
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
/*if(Multiplicity_ == 0){ ////for extended MinBias 
   MixEvents_eta(0,4);
   MixEvents_eta(5,9);
   MixEvents_eta(10,12);
   MixEvents_eta(13,15);
   MixEvents_eta(16,18);
   MixEvents_eta(19,21);
   MixEvents_eta(22,24);
   MixEvents_eta(25,27);
   MixEvents_eta(28,30);
   MixEvents_eta(31,33);
   MixEvents_eta(34,36);
   MixEvents_eta(37,39);
   MixEvents_eta(40,44);
   MixEvents_eta(45,49);
   MixEvents_eta(50,54);
   MixEvents_eta(55,59);
   MixEvents_eta(60,69);
   MixEvents_eta(70,79);
   MixEvents_eta(80,89);
   MixEvents_eta(90,99);
   MixEvents_eta(100,109);
   MixEvents_eta(110,119);
   MixEvents_eta(120,250);
}*/  
else if(Multiplicity_ == 1){
   MixEvents_eta(80,84);
   MixEvents_eta(85,89);
   MixEvents_eta(90,94);
   MixEvents_eta(95,99);
   MixEvents_eta(100,104);
}  
else if(Multiplicity_ == 2){
   MixEvents_eta(105,109);
   MixEvents_eta(110,114);
   MixEvents_eta(115,119);
   MixEvents_eta(120,124);
   MixEvents_eta(125,129);
}  
else if(Multiplicity_ == 3){
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

}//end "if" condition for eta mix


}//endjob

bool 
NtuplizerHBT::vtxSort( const reco::Vertex &  a, const reco::Vertex & b ){
   if( a.tracksSize() != b.tracksSize() )
      return  a.tracksSize() > b.tracksSize() ? true : false ;
   else
      return  a.chi2() < b.chi2() ? true : false ;
}

bool
NtuplizerHBT::etaMixSort(std::pair<Double_t,std::vector<TLorentzVector>> & a, std::pair<Double_t,std::vector<TLorentzVector>> & b){

   return a.first < b.first ? true : false ;

}

bool
NtuplizerHBT::etaMixSort_ch(std::pair<Double_t,std::vector<int>> & a, std::pair<Double_t,std::vector<int>> & b){

   return a.first < b.first ? true : false ;

}

bool
NtuplizerHBT::etaMixSort_ntrkoff(std::pair<Double_t,int> & a, std::pair<Double_t,int> & b){

   return a.first < b.first ? true : false ;

}


bool 
NtuplizerHBT::splitcomb(TLorentzVector &vec1,TLorentzVector &vec2){
   bool issplit=false;
   Double_t cosa = TMath::Abs(vec1.Px()*vec2.Px() + vec1.Py()*vec2.Py() + vec1.Pz()*vec2.Pz())/(vec1.P()*vec2.P());
   Double_t deltapt = TMath::Abs(vec1.Pt() - vec2.Pt());
   if ( (cosa >cos_cut) && (deltapt < dpt_cut)) { issplit = true;}
   //std::cout << "cosa: " << cosa << " dpt: " << deltapt << " is split: " << issplit << std::endl;
   return issplit;
}

Double_t 
NtuplizerHBT::GetQ(const TLorentzVector &p1, const TLorentzVector &p2){
   TLorentzVector Sum4V = p1+p2;
   Double_t q = Sum4V.Mag2() - 4*pi_mass*pi_mass;
   //  std::cout<<(  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  ) <<std::endl;
   return (  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}

Double_t
NtuplizerHBT::GetQlong(const TLorentzVector &p1, const TLorentzVector &p2){
   TLorentzVector Diff4V = p1-p2;
   Double_t qlong = fabs(Diff4V.Pz());
   return qlong;
}

Double_t
NtuplizerHBT::GetQlongLCMS(const TLorentzVector &p1, const TLorentzVector &p2){
   Double_t num = 2*( (p1.Pz())*(p2.E()) - (p2.Pz())*(p1.E()) );
   Double_t den = TMath::Sqrt( (p1.E()+p2.E())*(p1.E()+p2.E()) - (p1.Pz()+p2.Pz())*(p1.Pz()+p2.Pz()) );
   Double_t qlongLCMS = 0.0;
   if(den != 0) qlongLCMS = fabs(num/den);
   return qlongLCMS;
}

Double_t
NtuplizerHBT::GetQout(const TLorentzVector &p1, const TLorentzVector &p2){

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
NtuplizerHBT::GetQside(const TLorentzVector &p1, const TLorentzVector &p2){

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
NtuplizerHBT::InvertPVector( TLorentzVector &vec){
   TLorentzVector ovec = vec;
   ovec.SetPx(-vec.Px());
   ovec.SetPy(-vec.Py());
   ovec.SetPz(-vec.Pz());
   return ovec;
}

const TLorentzVector 
NtuplizerHBT::InvertXYVector( TLorentzVector &vec){
  TLorentzVector ovec = vec;
  ovec.SetX(-vec.X());
  ovec.SetY(-vec.Y());
  return ovec;
}

//return the weight factor due to Coloumb repulsion [Gamow] same charge
const Double_t 
NtuplizerHBT::CoulombW(const Double_t& q){
   const Double_t alpha=1./137.;
   Double_t x=2*TMath::Pi()*(alpha*pi_mass/q);
   //return (TMath::Exp(x)-1.)/x; // OLD MATTIA's DEFINITION
  
   //Double_t ws=scf*((exp(arg)-1)/arg-1)+1; // PAOLO's DEFINITION
   Double_t weight = 1;//0.85; // TEMPORARY SET TO 0.85 * GAMOW FACTOR
   //Double_t weight = 1.15; //for syst. +15%
   //Double_t weight = 0.85; //for syst. -15%
   return weight*( (TMath::Exp(x)-1.)/x -1 ) + 1;
}

//return the weight factor due to Coloumb attraction [Gamow] opposite charge
const  Double_t 
NtuplizerHBT::CoulombWpm(const Double_t& q){
   const Double_t alpha=1./137.;
   Double_t x=2*TMath::Pi()*(alpha*pi_mass/q);
   // return (1.-TMath::Exp(-x))/x; // OLD MATTIA's DEFINITION
  
   // Double_t wd=scf*((1-exp(-arg))/arg-1)+1; // PAOLO's DEFINITION
   Double_t weight = 1;//0.85; // TEMPORARY SET TO 0.85 * GAMOW FACTOR
   //Double_t weight = 1.15; //for syst. +15%
   //Double_t weight = 0.85; //for syst. -15%
   return weight*( (1.-TMath::Exp(-x))/x -1 ) + 1;
}

Double_t 
NtuplizerHBT::ComputeEventWeight(std::vector<TLorentzVector> GoodTrackFourVector) // compute eta of event
{
        int nFwd(0), nBkw(0), nCtr(0) ;
        for( unsigned int i=0 ; i< GoodTrackFourVector.size() ; i++ )
    {
                Double_t eta = GoodTrackFourVector[i].Eta() ;
                if( eta < -0.8 )
        { nBkw++ ;}
                else if( eta > 0.8 )
        { nFwd++ ;}
                else
        { nCtr++ ;}
    }


        Double_t ReturnValue = (100*Encode( nFwd ) + 10 * Encode( nCtr ) + Encode( nBkw )) ;
        return ReturnValue ;
}

Double_t 
NtuplizerHBT::Encode( int w )
{
        //int Range[]={3,5,8,10,13,16,20,25,30,200} ;
        int Range[]={3,5,8,10,13,16,20,25,30,200,10000} ;//added a protection for high-multiplicity 
        int i(0), j(0) ;
        while( w >= Range[i++] ) j++ ;
        //henc->Fill(w, (Double_t) j) ;
        return (Long64_t) j ;
}


int
NtuplizerHBT::ntrkoff(edm::Handle<std::vector<reco::Track> > generalTracks,double vtx_x,double vtx_y,double vtx_z,std::vector<TLorentzVector> &GoodTrackFourVector_trkoff){
   math::XYZPoint vtx(vtx_x,vtx_y,vtx_z);
   int aux_ntrkoffline=0;
   double aux_ntrkoffline_corr=0.0;
   for(reco::TrackCollection::const_iterator iter_tk = generalTracks->begin(); iter_tk != generalTracks->end(); iter_tk++){
      double aux_tk_dz_vtx = (double)iter_tk->dz(vtx);
      double aux_tk_dzError_vtx  = (double)sqrt(iter_tk->dzError()*iter_tk->dzError()+vtx_zError*vtx_zError);
      double aux_tk_dxy_vtx = (double)iter_tk->dxy(vtx);
      double aux_tk_dxyError_vtx  = (double)sqrt(iter_tk->dxyError()*iter_tk->dxyError()+vtx_xError*vtx_yError);
      if(iter_tk->pt()<0.4)continue;
      if(fabs(iter_tk->eta())>2.4)continue;
      if(!iter_tk->quality(reco::TrackBase::highPurity))continue;
      if(fabs(iter_tk->ptError())/iter_tk->pt()>0.1)continue;
      if(fabs(aux_tk_dz_vtx/aux_tk_dzError_vtx)>3)continue;
      if(fabs(aux_tk_dxy_vtx/aux_tk_dxyError_vtx)>3)continue;
      aux_ntrkoffline++; 
      double aux_tk_corr = getTrkCorrWeight(iter_tk->pt(), iter_tk->eta());
      //std::cout<<"trk eta : "<< iter_tk->eta()<<" tk pt : "<<iter_tk->pt()<<"  weight : "
      //               <<aux_tk_corr<<std::endl; 
      tk_offlineSelCorr_->Fill(aux_tk_corr);
      aux_ntrkoffline_corr = aux_ntrkoffline_corr + 1.0*aux_tk_corr;

      TLorentzVector pvector;
      pvector.SetXYZM(iter_tk->px(),iter_tk->py(),iter_tk->pz(),pi_mass);
      GoodTrackFourVector_trkoff.push_back(pvector);  

   }
   tk_NtrkOfflineCorr_NtrkOffline_->Fill(aux_ntrkoffline,aux_ntrkoffline_corr); 
   return aux_ntrkoffline;
}

double
NtuplizerHBT::getTrkCorrWeight(double pT, double eta){

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


void
NtuplizerHBT::fillTreeWithCommonVariables(math::XYZPoint BSPosition, math::XYZPoint vtx, 
                                          double vtx_xError, double vtx_yError, double vtx_zError, 
                                          std::vector<reco::TrackCollection::const_iterator> selected_tracks){

N_tk = selected_tracks.size();

for(unsigned int itk = 0; itk<selected_tracks.size(); itk++){
   tk_charge->push_back(selected_tracks[itk]->charge());
   tk_algo->push_back(selected_tracks[itk]->algo());
   tk_pt->push_back(selected_tracks[itk]->pt());
   tk_ptError->push_back(selected_tracks[itk]->ptError());
   tk_p->push_back(selected_tracks[itk]->p());
   tk_px->push_back(selected_tracks[itk]->px());
   tk_py->push_back(selected_tracks[itk]->py());
   tk_pz->push_back(selected_tracks[itk]->pz());
   tk_eta->push_back(selected_tracks[itk]->eta());
   tk_etaError->push_back(selected_tracks[itk]->etaError());
   tk_phi->push_back(selected_tracks[itk]->phi());
   tk_phiError->push_back(selected_tracks[itk]->phiError());
   tk_dxy->push_back(selected_tracks[itk]->dxy());
   tk_dxy_bs->push_back(selected_tracks[itk]->dxy(BSPosition));
   tk_dxy_vtx->push_back(selected_tracks[itk]->dxy(vtx));
   tk_dxyError->push_back(selected_tracks[itk]->dxyError());
   double aux_tk_dxyError_vtx  = (double)sqrt(selected_tracks[itk]->dxyError()*selected_tracks[itk]->dxyError()+vtx_xError*vtx_yError);
   tk_dxyError_vtx->push_back(aux_tk_dxyError_vtx);
   tk_dz->push_back(selected_tracks[itk]->dz());
   tk_dz_bs->push_back(selected_tracks[itk]->dz(BSPosition));
   tk_dz_vtx->push_back(selected_tracks[itk]->dz(vtx));
   tk_dzError->push_back(selected_tracks[itk]->dzError());
   double aux_tk_dzError_vtx  = (double)sqrt(selected_tracks[itk]->dzError()*selected_tracks[itk]->dzError()+vtx_zError*vtx_zError);
   tk_dzError_vtx->push_back(aux_tk_dzError_vtx);
   tk_chi2->push_back(selected_tracks[itk]->chi2());
   tk_ndof->push_back(selected_tracks[itk]->ndof());
   tk_normalizedChi2->push_back(selected_tracks[itk]->normalizedChi2());       
   tk_numberOfValidHits->push_back(selected_tracks[itk]->numberOfValidHits());
   tk_numberOfLostHits->push_back(selected_tracks[itk]->numberOfLostHits()); 
   
   
}


}


void
NtuplizerHBT::fillHistosWithCommonVariables(math::XYZPoint BSPosition, math::XYZPoint vtx,
                                          double vtx_xError, double vtx_yError, double vtx_zError,
                                          std::vector<reco::TrackCollection::const_iterator> selected_tracks){

for(unsigned int itk = 0; itk<selected_tracks.size(); itk++){


   double aux_tk_corr = getTrkCorrWeight(selected_tracks[itk]->pt(), selected_tracks[itk]->eta());

   tk_NHits_->Fill(selected_tracks[itk]->numberOfValidHits(),aux_tk_corr);
   tk_Pt_->Fill(selected_tracks[itk]->pt(),aux_tk_corr);
   tk_Eta_->Fill(selected_tracks[itk]->eta(),aux_tk_corr);
   tk_Phi_->Fill(selected_tracks[itk]->phi(),aux_tk_corr);
   double aux_tk_dz_vtx = selected_tracks[itk]->dz(vtx);
   double aux_tk_dzError_vtx  = (double)sqrt(selected_tracks[itk]->dzError()*selected_tracks[itk]->dzError()+vtx_zError*vtx_zError);
   double aux_tk_dzOverDzError_vtx = aux_tk_dz_vtx/aux_tk_dzError_vtx;
   tk_DzOverSigmaDz_->Fill(aux_tk_dzOverDzError_vtx,aux_tk_corr);
   double aux_tk_dxy_vtx = selected_tracks[itk]->dxy(vtx);
   double aux_tk_dxyError_vtx  = (double)sqrt(selected_tracks[itk]->dxyError()*selected_tracks[itk]->dxyError()+vtx_xError*vtx_yError);
   double aux_tk_dxyOverDxyError_vtx = aux_tk_dxy_vtx/aux_tk_dxyError_vtx;
   tk_DxyOverSigmaDxy_->Fill(aux_tk_dxyOverDxyError_vtx,aux_tk_corr); 
   double aux_tk_ptError = selected_tracks[itk]->ptError();
   double aux_tk_pt = selected_tracks[itk]->pt();
   double aux_tk_ptErrorOverPt = aux_tk_ptError/aux_tk_pt;
   tk_SigmaPtOverPt_->Fill(aux_tk_ptErrorOverPt,aux_tk_corr); 

   tk_NormalizedChi2_->Fill(selected_tracks[itk]->normalizedChi2(),aux_tk_corr); 
   tk_Ndof_->Fill(selected_tracks[itk]->ndof(),aux_tk_corr);
   //std::cout<<"nhit   :   "<< selected_tracks[itk]->numberOfValidHits() << std::endl;
   //std::cout<<"ndof   :   "<< selected_tracks[itk]->ndof() << std::endl;
   //std::cout<<"======================================================"<<std::endl;

   const reco::HitPattern& hit_pattern = selected_tracks[itk]->hitPattern();
   int aux_tk_pixelLayersWithMeasurement = hit_pattern.pixelLayersWithMeasurement();   
   tk_pixelLayersWithMeasurement_->Fill(aux_tk_pixelLayersWithMeasurement,aux_tk_corr); 
  
}

}


void NtuplizerHBT::initHistos(const edm::Service<TFileService> & fs){


TH1D::SetDefaultSumw2();
TFileDirectory trkQA = fs->mkdir( "TrackHistograms" );

Selections_ = trkQA.make<TH1F>("Selections_","Selection Number",7,-0.5,6.5);
Selections_->GetXaxis()->SetTitle("Selection number");

tk_offlineSelCorr_ = trkQA.make<TH1F>("tk_offlineSelCorr_","Per track corrections for trackOffline selection",100,-2,4);
tk_offlineSelCorr_->GetXaxis()->SetTitle("Track Corrections");

tk_NtrkOfflineCorr_NtrkOffline_ = trkQA.make<TH2F>("tk_NtrkOfflineCorr_NtrkOffline_","NtrkOfflineCorrected Vs NtrkOfflineRaw",250,-0.5,249.5,900,0.,300);

VertexX_ = trkQA.make<TH1F>("VertexX_", "x position of the PV", 100, -1, 1);
VertexX_->GetXaxis()->SetTitle("x [cm]");
VertexY_ = trkQA.make<TH1F>("VertexY_", "y position of the PV", 100, -1, 1);
VertexY_->GetXaxis()->SetTitle("y [cm]");
VertexZ_ = trkQA.make<TH1F>("VertexZ_", "z position of the PV", 100, -30, 30);
VertexZ_->GetXaxis()->SetTitle("z [cm]");
NVertex_ = trkQA.make<TH1F>("NVertex_", "Number of vertices", 11, -0.5,10.5);
NVertex_->GetXaxis()->SetTitle("N_{vtx}");

tk_NHits_ = trkQA.make<TH1F>("tk_NHits_", "Tracks by number of valid hits", 101, -0.5, 100.5);
tk_NHits_->GetXaxis()->SetTitle("Number of valid hits");
tk_Pt_ = trkQA.make<TH1F>("tk_Pt_", "p_{T} spectrum", 100, 0, 20);
tk_Pt_->GetXaxis()->SetTitle("p_{T} [GeV]");
//tk_Eta_ = trkQA.make<TH1F>("tk_Eta_", "Pseudorapidity distribution", 50, -2.5, 2.5);
double etaBins[13] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0,
        0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
tk_Eta_ = trkQA.make<TH1F>("tk_Eta_", "Pseudorapidity distribution",12,etaBins); //same used in trk corrections
tk_Eta_->GetXaxis()->SetTitle("#eta distribution of the tracks");
tk_Phi_ = trkQA.make<TH1F>("tk_Phi_", "Azimuthal distribution",100, -3.15, 3.15);
tk_Phi_->GetXaxis()->SetTitle("#phi distribution of the tracks");

tk_DzOverSigmaDz_ = trkQA.make<TH1F>("tk_DzOverSigmaDz_", "Significance of z separation between Trk and PV",100, -10.0, 10.0);
tk_DzOverSigmaDz_->GetXaxis()->SetTitle("d_{z}/#sigma(d_{z})");
tk_DxyOverSigmaDxy_ = trkQA.make<TH1F>("tk_DxyOverSigmaDxy_", "Significance of the impact parameter with respect to PV",100, -10.0, 10.0);
tk_DxyOverSigmaDxy_->GetXaxis()->SetTitle("d_{xy}/#sigma(d_{xy})");
tk_SigmaPtOverPt_ = trkQA.make<TH1F>("tk_SigmaPtOverPt_", "Relative uncertainty of the p_{T} measurement",100, -0.5, 0.5);
tk_SigmaPtOverPt_->GetXaxis()->SetTitle("#sigma(p_{T})/p_{T}");

N_tk_HPTracks_ = trkQA.make<TH1F>("N_tk_HPTracks_","High-Purity tracks",250,-0.5,249.5);
N_tk_HPTracks_->GetXaxis()->SetTitle("Track Multiplicity");

N_tk_AllTracks_ = trkQA.make<TH1F>("N_tk_AllTracks_","All tracks",250,-0.5,249.5);
N_tk_AllTracks_->GetXaxis()->SetTitle("Track Multiplicity");

N_tk_SelTracks_ = trkQA.make<TH1F>("N_tk_SelTracks_","Selected tracks",250,-0.5,249.5);
N_tk_SelTracks_->GetXaxis()->SetTitle("Track Multiplicity");

N_tk_Offline_= trkQA.make<TH1F>("N_tk_Offline_","Track Offline Selection",250,-0.5,249.5);
N_tk_Offline_->GetXaxis()->SetTitle("Track Multiplicity");

tk_pixelLayersWithMeasurement_ = trkQA.make<TH1F>("tk_pixelLayersWithMeasurement_", "Number of Pixel Layers with Measurements", 6, -0.5, 5.5);
tk_pixelLayersWithMeasurement_->GetXaxis()->SetTitle("Number of Layers");

tk_NormalizedChi2_ = trkQA.make<TH1F>("tk_NormalizedChi2_", "Chi-squared divided by ndof",100, 0, 6);
tk_NormalizedChi2_->GetXaxis()->SetTitle("#chi^{2}/ndof");

tk_Ndof_ = trkQA.make<TH1F>("tk_Ndof_", "Number of Degrees of Freedom of the Fit",101, -0.5, 100.5);
tk_Ndof_->GetXaxis()->SetTitle("ndof");

tk_pairSS_M_ = trkQA.make<TH1F>("tk_pairSS_M_", "Invariant mass same-sign tracks", 100, 0, 1);
tk_pairSS_M_->GetXaxis()->SetTitle("M [GeV]");
tk_pairOS_M_ = trkQA.make<TH1F>("tk_pairOS_M_", "Invariant mass opposite-sign tracks", 100, 0, 1);
tk_pairOS_M_->GetXaxis()->SetTitle("M [GeV]");


TFileDirectory mixedEvents = fs->mkdir( "MixedHistograms" );

DeltazVertex_ = mixedEvents.make<TH1F>("DeltazVertex_","DeltazVertex_", 100, -2.5, 2.5);
mixsch_KtVsNch_ntrkoff= mixedEvents.make<TH2F>("mixsch_KtVsNch_ntrkoff","mixsch_KtVsNch_ntrkoff",200,0.,20.,250,-0.5,249.5);

//THnSparse
Int_t bins5D[5]=   {40, 40, 40, 10, 25  };
Double_t xmin5D[5]={0., 0., 0., 0., 0.  };
Double_t xmax5D[5]={2., 2., 2., 1., 250.};

Int_t bins3D[3]=   {200,20,250  };
Double_t xmin3D[3]={0. ,0.,-0.5 };
Double_t xmax3D[3]={2. ,1.,249.5};

hs_qinv_mixschVsKtVsNch_ntrkoff = mixedEvents.make<THnSparseF>("hs_qinv_mixschVsKtVsNch_ntrkoff","hs_qinv_mixschVsKtVsNch_ntrkoff",3, bins3D,  xmin3D, xmax3D);
hs_qinv_mixdchVsKtVsNch_ntrkoff = mixedEvents.make<THnSparseF>("hs_qinv_mixdchVsKtVsNch_ntrkoff","hs_qinv_mixdchVsKtVsNch_ntrkoff",3, bins3D,  xmin3D, xmax3D);
hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff = mixedEvents.make<THnSparseF>("hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff","hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff = mixedEvents.make<THnSparseF>("hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff","hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
 
TFileDirectory qinvHist = fs->mkdir( "qinvHistograms" );

sig_KtVsNch_ntrkoff= qinvHist.make<TH2F>("sig_KtVsNch_ntrkoff","sig_KtVsNch_ntrkoff",200,0.,20.,250,-0.5,249.5);
ptVsNch_ntrkoff = qinvHist.make<TH2F>("ptVsNch_ntrkoff","ptVsNch_ntrkoff",200,0.,20.,250,-0.5,249.5);
etaVsNch_ntrkoff = qinvHist.make<TH2F>("etaVsNch_ntrkoff","etaVsNch_ntrkoff",12,etaBins,250,-0.5,249.5); //same used in trk corrections
phiVsNch_ntrkoff = qinvHist.make<TH2F>("phiVsNch_ntrkoff","phiVsNch_ntrkoff",100, -3.15, 3.15,250,-0.5,249.5);
sig_pairSS_MVsNch_ntrkoff= qinvHist.make<TH2F>("sig_pairSS_MVsNch_ntrkoff","sig_pairSS_MVsNch_ntrkoff",200,0,2,250,-0.5,249.5);
pairOS_MVsNch_ntrkoff= qinvHist.make<TH2F>("pairOS_MVsNch_ntrkoff","pairOS_MVsNch_ntrkoff",200,0,2,250,-0.5,249.5);

//THnSparse
hs_qLLCMSqOqSschVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("hs_qLLCMSqOqSschVsKtVsNch_ntrkoff","hs_qLLCMSqOqSschVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
hs_qLLCMSqOqSschncVsKtVsNch_ntrkoff= qinvHist.make<THnSparseF>("hs_qLLCMSqOqSschncVsKtVsNch_ntrkoff","hs_qLLCMSqOqSschncVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
hs_qLLCMSqOqSdchVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("hs_qLLCMSqOqSdchVsKtVsNch_ntrkoff","hs_qLLCMSqOqSdchVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
hs_qLLCMSqOqSdchncVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("hs_qLLCMSqOqSdchncVsKtVsNch_ntrkoff","hs_qLLCMSqOqSdchncVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
hs_qschVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("hs_qschVsKtVsNch_ntrkoff","hs_qschVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
hs_qschncVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("hs_qschncVsKtVsNch_ntrkoff","hs_qschncVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
hs_qINVschVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("hs_qINVschVsKtVsNch_ntrkoff","hs_qINVschVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
hs_qINVschncVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("hs_qINVschncVsKtVsNch_ntrkoff","hs_qINVschncVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
hs_qROTschVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("hs_qROTschVsKtVsNch_ntrkoff","hs_qROTschVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
hs_qROTschncVsKtVsNch_ntrkoff= qinvHist.make<THnSparseF>("hs_qROTschncVsKtVsNch_ntrkoff","hs_qROTschncVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
hs_qdchVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("hs_qdchVsKtVsNch_ntrkoff","hs_qdchVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
hs_qdchncVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("hs_qdchncVsKtVsNch_ntrkoff","hs_qdchncVsKtVsNch_ntrkoff",3,bins3D,xmin3D,xmax3D);
hs_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("hs_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff","hs_qLLCMSqOqS_INV_schVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
hs_qLLCMSqOqS_INV_schncVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("hs_qLLCMSqOqS_INV_schncVsKtVsNch_ntrkoff","hs_qLLCMSqOqS_INV_schncVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
hs_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff = qinvHist.make<THnSparseF>("hs_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff","hs_qLLCMSqOqS_ROT_schVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);
hs_qLLCMSqOqS_ROT_schncVsKtVsNch_ntrkoff= qinvHist.make<THnSparseF>("hs_qLLCMSqOqS_ROT_schncVsKtVsNch_ntrkoff","hs_qLLCMSqOqS_ROT_schncVsKtVsNch_ntrkoff",5,bins5D,xmin5D,xmax5D);

}


void NtuplizerHBT::MixEvents(int ntrkoff_min, int ntrkoff_max, int nEvt_to_mix){

int aux_n_evts = (int)ev_ntrkoff_vec.size();

for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
  
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector_vec[nevt1];
   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   std::vector<int> ev_GoodTrackChargeTemp_nevt1_vec = ev_GoodTrackCharge_vec[nevt1];

   if(ev_ntrkoff_vec[nevt1]<ntrkoff_min || ntrkoff_max<ev_ntrkoff_vec[nevt1])continue;

   int takeAssociated = 0;
   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
      if(ev_ntrkoff_vec[nevt_assoc]<ntrkoff_min || ntrkoff_max<ev_ntrkoff_vec[nevt_assoc])continue;
 
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector_vec[nevt_assoc]; 
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
      std::vector<int> ev_GoodTrackChargeTemp_nevt_assoc_vec=ev_GoodTrackCharge_vec[nevt_assoc];
       
      takeAssociated++;

      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.

      DeltazVertex_->Fill((ev_vtx_z_vec)[nevt1] - (ev_vtx_z_vec)[nevt_assoc]);
      for(int imix=0; imix<nMixmult_nevt1; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
            if(splitcomb(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix])){continue;}
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            //Double_t qmixlong = GetQlong(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            Double_t qmixlongLCMS = GetQlongLCMS(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            Double_t qmixout = GetQout(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            Double_t qmixside = GetQside(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);

            TLorentzVector psum2 = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
            Double_t kt=(psum2.Pt())/2.; 

            Double_t xmix5D_LCMS[5]={qmixlongLCMS,qmixout,qmixside,kt,(Double_t)ev_ntrkoff_vec[nevt1]};
            Double_t xmix3D[3]={qmix,kt,(Double_t)ev_ntrkoff_vec[nevt1]};

            Double_t aux_tk1_corr = (Double_t)getTrkCorrWeight(ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Pt(),ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Eta());
            Double_t aux_tk2_corr = (Double_t)getTrkCorrWeight(ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Pt(),ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Eta());
            Double_t aux_tk12_corr= aux_tk1_corr*aux_tk2_corr;
            //aux_tk12_corr=1.; //just to compare with applying corrections 
            if(ev_GoodTrackChargeTemp_nevt1_vec[imix]*ev_GoodTrackChargeTemp_nevt_assoc_vec[iimix]>0){
 
               mixsch_KtVsNch_ntrkoff->Fill(kt,ev_ntrkoff_vec[nevt1],aux_tk12_corr);             

               hs_qinv_mixschVsKtVsNch_ntrkoff->Fill(xmix3D,aux_tk12_corr); 
               hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff->Fill(xmix5D_LCMS,aux_tk12_corr);

            }else{

               hs_qinv_mixdchVsKtVsNch_ntrkoff->Fill(xmix3D,aux_tk12_corr); 
               hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff->Fill(xmix5D_LCMS,aux_tk12_corr);

            }
         }
      }
   }

    
}//end first for loop 

}//end MixEvents function


void NtuplizerHBT::MixEvents_eta(int ntrkoff_min, int ntrkoff_max){


///sort vectors of maps for each ntrkoffline range
std::vector<std::pair<Double_t, std::vector<TLorentzVector>> > aux_ev_GoodTrackFourVector_etaMixWeight_vec;
std::vector<std::pair<Double_t, std::vector<int>> > aux_ev_GoodTrackCharge_etaMixWeight_vec;
std::vector<std::pair<Double_t,int> > aux_ev_ntrkoff_etaMixWeight_vec;


int aux_n_evts = ev_GoodTrackFourVector_etaMixWeight_vec.size();

for(int ievt=0; ievt<aux_n_evts; ievt++){

   if((ev_ntrkoff_etaMixWeight_vec[ievt]).second<ntrkoff_min || ntrkoff_max<(ev_ntrkoff_etaMixWeight_vec[ievt]).second)continue; //only events in a given range

   aux_ev_GoodTrackFourVector_etaMixWeight_vec.push_back(ev_GoodTrackFourVector_etaMixWeight_vec[ievt]);
   aux_ev_GoodTrackCharge_etaMixWeight_vec.push_back(ev_GoodTrackCharge_etaMixWeight_vec[ievt]);
   aux_ev_ntrkoff_etaMixWeight_vec.push_back(ev_ntrkoff_etaMixWeight_vec[ievt]); 

}

std::sort( aux_ev_GoodTrackFourVector_etaMixWeight_vec.begin(), aux_ev_GoodTrackFourVector_etaMixWeight_vec.end(), NtuplizerHBT::etaMixSort );
std::sort( aux_ev_GoodTrackCharge_etaMixWeight_vec.begin(), aux_ev_GoodTrackCharge_etaMixWeight_vec.end(), NtuplizerHBT::etaMixSort_ch );
std::sort( aux_ev_ntrkoff_etaMixWeight_vec.begin(), aux_ev_ntrkoff_etaMixWeight_vec.end(), NtuplizerHBT::etaMixSort_ntrkoff );


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
         Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
         //Double_t qmixlong = GetQlong(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
         Double_t qmixlongLCMS = GetQlongLCMS(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
         Double_t qmixout = GetQout(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
         Double_t qmixside = GetQside(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);        

         TLorentzVector psum2 = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
         Double_t kt=(psum2.Pt())/2.;

         Double_t xmix3D[3]={qmix,kt,(Double_t)(aux_ev_ntrkoff_etaMixWeight_vec[nevt]).second}; 
         Double_t xmix5D_LCMS[5]={qmixlongLCMS,qmixout,qmixside,kt,(Double_t)(aux_ev_ntrkoff_etaMixWeight_vec[nevt]).second};

         Double_t aux_tk1_corr = (Double_t)getTrkCorrWeight(ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Pt(),ev_GoodTrackFourVectorTemp_nevt1_vec[imix].Eta());
         Double_t aux_tk2_corr = (Double_t)getTrkCorrWeight(ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Pt(),ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix].Eta());
         Double_t aux_tk12_corr= aux_tk1_corr*aux_tk2_corr;
         if(ev_GoodTrackChargeTemp_nevt1_vec[imix]*ev_GoodTrackChargeTemp_nevt_assoc_vec[iimix]>0){

            mixsch_KtVsNch_ntrkoff->Fill(kt,(aux_ev_ntrkoff_etaMixWeight_vec[nevt]).second,aux_tk12_corr);

            hs_qinv_mixschVsKtVsNch_ntrkoff->Fill(xmix3D,aux_tk12_corr);
            hs_qLLCMSqOqSmixschVsKtVsNch_ntrkoff->Fill(xmix5D_LCMS,aux_tk12_corr);

         }else{
            
            hs_qinv_mixdchVsKtVsNch_ntrkoff->Fill(xmix3D,aux_tk12_corr);
            hs_qLLCMSqOqSmixdchVsKtVsNch_ntrkoff->Fill(xmix5D_LCMS,aux_tk12_corr);

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
NtuplizerHBT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NtuplizerHBT);
