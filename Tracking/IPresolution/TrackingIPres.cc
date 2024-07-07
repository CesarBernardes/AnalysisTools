// -*- C++ -*-
//
// Package:    TrackingIPres/TrackingIPres
// Class:      TrackingIPres
//
/**\class TrackingIPres TrackingIPres.cc TrackingIPres/TrackingIPres/plugins/TrackingIPres.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Cesar Bernardes
//         Created:  Mon, 10 Jun 2024 18:45:31 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TH2F.h"
#include "TH1F.h"
#include "TVector3.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

///using reco::TrackCollection;

class TrackingIPres : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TrackingIPres(const edm::ParameterSet&);
  ~TrackingIPres() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void initHistos( const edm::Service< TFileService > & fs ); 

  // ----------member data ---------------------------
/*  edm::EDGetTokenT< edm::View< pat::PackedCandidate > > tracksToken_;
  edm::EDGetTokenT< edm::View< pat::PackedGenParticle > > genpToken_;
  edm::EDGetTokenT< reco::VertexCollection > vertexToken_;*/
  edm::EDGetTokenT< edm::View< reco::Track > > tracksToken_;
  edm::EDGetTokenT< edm::View< reco::GenParticle > > genpToken_;
  edm::EDGetTokenT< reco::VertexCollection > vertexToken_;


  TH1F* hist_trkpt_;
  TH1F* hist_trketa_;
  TH1F* hist_trkphi_;
  TH1F* hist_trkdzsig_;
  TH1F* hist_trkdxysig_;
  TH1F* hist_trkptres_;

  TH1F* hist_dR_reco_gen_;
  TH1F* hist_dpT_reco_gen_; 

  TH2F* hist_2D_trkIPresXY_PtReco_barrel_;
  TH2F* hist_2D_trkIPresZ_PtReco_barrel_;
  TH2F* hist_2D_trkIPresXY_PtReco_forward_;
  TH2F* hist_2D_trkIPresZ_PtReco_forward_;
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
TrackingIPres::TrackingIPres(const edm::ParameterSet& iConfig)
    : tracksToken_( consumes< edm::View< reco::Track > >( iConfig.getParameter< edm::InputTag >( "tracks" ) ) ),
      genpToken_( consumes< edm::View< reco::GenParticle > >( iConfig.getParameter< edm::InputTag >( "genparticles" ) ) ),
      vertexToken_( consumes< reco::VertexCollection >( iConfig.getParameter< edm::InputTag >( "vertices" ) ) ){
//	    tracksToken_( consumes< edm::View< pat::PackedCandidate > >( iConfig.getParameter< edm::InputTag >( "tracks" ) ) ),
//      genpToken_( consumes< edm::View< pat::PackedGenParticle > >( iConfig.getParameter< edm::InputTag >( "genparticles" ) ) ),
//      vertexToken_( consumes< reco::VertexCollection >( iConfig.getParameter< edm::InputTag >( "vertices" ) ) ){
  //now do what ever initialization is needed
}

TrackingIPres::~TrackingIPres() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void TrackingIPres::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;

  //track collection
  auto trks = iEvent.getHandle( tracksToken_ );

  //gen particle collection
  auto genp = iEvent.getHandle( genpToken_ );

  //vtx collection
  auto pvs = iEvent.getHandle( vertexToken_ );

  //best vertex
  double bestvzError;
  math::XYZPoint bestvtx;
  math::Error<3>::type vtx_cov;
  if ( !pvs->empty() ) {
     const reco::Vertex& vtx = (*pvs)[0];
     bestvzError = vtx.zError();
     bestvtx = vtx.position();
     vtx_cov = vtx.covariance();      
  }else { 
     return; 
  } 

  if ( std::abs( bestvtx.z() ) > 15. ) return; 

  //loop over reco tracks
  //int trkIndx2 = -1;
  for (auto const& trk2 : *trks) {   

     //trkIndx2++;      

     ///if ( !trk2.hasTrackDetails() ) continue;
     ///auto iter_tk2 = trk2.pseudoTrack();

     //auto iter_tk2 = trk2.at(trkIndx2);

     /*double trk2_eta = iter_tk2.eta();
     double trk2_pt = iter_tk2.pt();
     double trk2_phi = iter_tk2.phi();
     ///if ( std::abs( trk2_eta ) >= 3.0 ) continue;

     double dzvtx2 = iter_tk2.dz( bestvtx );
     double dxyvtx2 = iter_tk2.dxy( bestvtx );
     double dzerror2 = std::hypot( iter_tk2.dzError(), bestvzError );
     double dxyerror2 = iter_tk2.dxyError( bestvtx, vtx_cov );

     if ( std::abs( iter_tk2.ptError() ) / trk2_pt >= 0.1 ) continue;
     if ( std::abs( dzvtx2 / dzerror2 ) >= 3 ) continue;
     if ( std::abs( dxyvtx2 / dxyerror2 ) >= 3 ) continue;*/

     double trk2_eta = trk2.eta();
     double trk2_pt = trk2.pt();
     double trk2_phi = trk2.phi();
     if ( std::abs( trk2_eta ) >= 4.0 ) continue;

     //Cesar: removed cuts to avoid deformation of distribution
     double dzvtx2 = trk2.dz( bestvtx );
     double dxyvtx2 = trk2.dxy( bestvtx );
     double dzerror2 = std::hypot( trk2.dzError(), bestvzError );
     double dxyerror2 = trk2.dxyError( bestvtx, vtx_cov );

     /*if ( trk2.quality(reco::TrackBase::qualityByName("highPurity")) != 1 ) continue;
     if ( std::abs( trk2.ptError() ) / trk2_pt >= 0.1 ) continue;
     if ( std::abs( dzvtx2 / dzerror2 ) >= 3 ) continue;
     if ( std::abs( dxyvtx2 / dxyerror2 ) >= 3 ) continue;*/

     hist_trkpt_->Fill(trk2_pt);
     hist_trketa_->Fill(trk2_eta);
     hist_trkphi_->Fill(trk2_phi);
     hist_trkdzsig_->Fill(dzvtx2 / dzerror2);
     hist_trkdxysig_->Fill(dxyvtx2 / dxyerror2);
     //hist_trkptres_->Fill(std::abs( iter_tk2.ptError() ) / trk2_pt);
     hist_trkptres_->Fill(std::abs( trk2.ptError() ) / trk2_pt);

     for (auto const& gp2 : *genp){

        if(gp2.status()!=1)continue; //only final state particles
        if(gp2.charge()==0)continue;
        if(std::abs(gp2.eta())>=4)continue;

        double gp2_pt = gp2.pt();
        double gp2_eta = gp2.eta();
        double gp2_phi = gp2.phi();

        TVector3 vec_gp2;
        vec_gp2.SetPtEtaPhi(gp2_pt, gp2_eta, gp2_phi);
        TVector3 vec_trk2;
        vec_trk2.SetPtEtaPhi(trk2_pt, trk2_eta, trk2_phi);
        double dR2 = vec_gp2.DeltaR(vec_trk2);
        double dpT2 = gp2_pt - trk2_pt;
        hist_dR_reco_gen_->Fill(dR2);
        if(dR2<0.03)hist_dpT_reco_gen_->Fill(dpT2);
        if(dR2<0.03 && std::abs(dpT2)<0.1){ 
           //Based on CMS-TRK-11-001		
		//std::cout<<"gen trk (vx,vy,vz) = "<<gp2.vx()<<" "<<gp2.vy()<<" "<<gp2.vz()<<std::endl;
		//std::cout<<"gen trk (px,py,pz) = "<<gp2.px()<<" "<<gp2.py()<<" "<<gp2.pz()<<std::endl;
		//std::cout<<"reco PV (x,y,z) = "<<bestvtx.x()<<" "<<bestvtx.y()<<" "<<bestvtx.z()<<std::endl;
           double dxy_gen_vtx = (1.0*(gp2.vx() -bestvtx.x())*gp2.py() + 1.0*(gp2.vy()-bestvtx.y())*gp2.px())/gp2.pt();
           double dz_gen_vtx = (gp2.vz()-bestvtx.z()) - ((gp2.vx()-bestvtx.x())*gp2.px() + (gp2.vy()-bestvtx.y())*gp2.py())/gp2.pt() * (gp2.pz()/gp2.pt());
	   double IPres_xy =dxyvtx2 - dxy_gen_vtx; //in centimeter
	   double IPres_z = dzvtx2 - dz_gen_vtx; //in centimeter

	   if(std::abs(trk2_eta)<1){ 
	      hist_2D_trkIPresXY_PtReco_barrel_->Fill(trk2_pt,IPres_xy);
              hist_2D_trkIPresZ_PtReco_barrel_->Fill(trk2_pt, IPres_z);
	   }else{
	      hist_2D_trkIPresXY_PtReco_forward_->Fill(trk2_pt,IPres_xy);
              hist_2D_trkIPresZ_PtReco_forward_->Fill(trk2_pt, IPres_z);	   
           } 
           break; 
        }
     }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void TrackingIPres::beginJob() {
  // please remove this method if not needed
  edm::Service< TFileService > fs;
  initHistos( fs );
}

// ------------ method called once each job just after ending the event loop  ------------
void TrackingIPres::endJob() {
  // please remove this method if not needed
}

void TrackingIPres::initHistos(const edm::Service< TFileService > & fs ) {

   TFileDirectory Inf = fs->mkdir( "Histograms" );

   hist_trkpt_ = Inf.make< TH1F >( "TrkPt", "", 200, 0.0, 20.0 );
   hist_trketa_ = Inf.make< TH1F >( "TrkEta", "", 100, -6.0, 6.0 );
   hist_trkphi_ = Inf.make< TH1F >( "TrkPhi", "", 100, -3.2, 3.2 );
   hist_trkdzsig_ = Inf.make< TH1F >( "TrkDzSig", "", 100, -5, 5 );
   hist_trkdxysig_ = Inf.make< TH1F >( "TrkDxySig", "", 100, -5, 5 );
   hist_trkptres_ = Inf.make< TH1F >( "TrkPtRes", "", 100, 0, 0.15 ); 

   ///for efficiency and fake rate - numerators and denominators
   hist_dR_reco_gen_ = Inf.make< TH1F >( "dR_Reco_Gen", "", 1000, 0.0, 1.0 );
   hist_dpT_reco_gen_ = Inf.make< TH1F >( "dpT_Reco_Gen", "", 1000, -1.5, 1.5 );

   hist_2D_trkIPresXY_PtReco_barrel_= Inf.make< TH2F >("trkIPresXY_PtReco_barrel","",200,0,20,400,-0.4,0.4);
   hist_2D_trkIPresZ_PtReco_barrel_= Inf.make< TH2F >("trkIPresZ_PtReco_barrel","",200,0,20,400,-0.4,0.4);
   hist_2D_trkIPresXY_PtReco_forward_= Inf.make< TH2F >("trkIPresXY_PtReco_forward","",200,0,20,400,-0.4,0.4);
   hist_2D_trkIPresZ_PtReco_forward_= Inf.make< TH2F >("trkIPresZ_PtReco_forward","",200,0,20,400,-0.4,0.4);

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TrackingIPres::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  ///edm::ParameterSetDescription desc;
  ///desc.setUnknown();
  ///descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //descriptions.addWithDefaultLabel(desc);
  edm::ParameterSetDescription desc;
  ///desc.add<edm::InputTag>("vertices", {"offlineSlimmedPrimaryVertices"});
  ///desc.add<edm::InputTag>("tracks", {"packedPFCandidates"});
  ///desc.add<edm::InputTag>("genparticles", {"packedGenParticles"});
  desc.add<edm::InputTag>("vertices", {"offlinePrimaryVertices"});
  desc.add<edm::InputTag>("tracks", {"generalTracks"});
  desc.add<edm::InputTag>("genparticles", {"genParticles"});
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackingIPres);
