
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//------ EXTRA HEADER FILES--------------------//
#include "math.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"

// for tracking information
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

// for vertex information 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//for centrality
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

//for beamspot information
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// triggers
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// ROOT
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TTree.h>
#include <TDirectory.h>
#include "Math/Point3D.h"

//
// class declaration
//
class AnalyzerFlow : public edm::EDAnalyzer {
   public:
      explicit AnalyzerFlow(const edm::ParameterSet&);
      ~AnalyzerFlow();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      
      // user routines (detailed description given with the method implementations)
      int SelectEvent(const edm::Event& iEvent);
      int SelectPrimaryVertex(const edm::Handle<reco::VertexCollection>& primVertex);
      void InitBranchVars();
      //bool splitcomb(TLorentzVector &vec1,TLorentzVector &vec2);
      //Double_t GetQ( const TLorentzVector &p1, const TLorentzVector &p2 );
      //const Double_t CoulombWpm( const Double_t& q);
      //const Double_t CoulombW( const Double_t& q);
      //void MixEvents(float hfsum_min, float hfsum_max, int nEvt_to_mix);

      // input tags
      edm::InputTag _inputTagPrimaryVertex;
      edm::InputTag _inputTagCentrality;
      edm::InputTag _inputTagTrack;

      // general flags and variables
      int _flagRECO;
      int _nevents;
      int _neventsSelected;
      
      // storage
      TFile* _file;
      TTree* _tree;
      
      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      // >>>>>>>>>>>>>>>> event variables >>>>>>>>>>>>>>>>>>>>>>>
      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      // (their description given when tree branches are created)
      // event
      int _evRunNumber;
      int _evEventNumber;

      // primary vertex
      int _Npv;
      int _pvNDOF;
      float _pvZ;
      float _pvRho;
      /*std::vector<double> ev_vtx_z_vec;
      std::vector< std::vector<TLorentzVector> > ev_GoodTrackFourVector_vec;
      std::vector< std::vector<int> > ev_GoodTrackCharge_vec;
      std::vector<float> ev_hfsum_vec;*/
      // centrality
      float _HFsumETPlus;
      float _HFsumETMinus;
      float _HFsumET;
      // tracks
      static const int _maxNtrk = 10000;
      int _Ntrk;
      float _trkPt[_maxNtrk];
      float _trkEta[_maxNtrk];
      float _trkPhi[_maxNtrk];
      float _trkPtRes[_maxNtrk];
      float _trkDzSig[_maxNtrk];
      float _trkDxySig[_maxNtrk];
      float _trkNpixLayers[_maxNtrk];
      ///qinv
      /*static const int _maxNpair = _maxNtrk*_maxNtrk;
      int _NSSpair;
      int _NOSpair;
      int _NpairMix;
      float _qinvSigSS[_maxNpair];
      float _coulombWSS[_maxNpair];
      float _qinvSigOS[_maxNpair];
      float _coulombWOS[_maxNpair];
      float _qinvMixSS_OS[_maxNpair];*/ //here we will mix with just one random event to have similar stats as in the Signal sample

};


//
// constructor
//
AnalyzerFlow::AnalyzerFlow(const edm::ParameterSet& iConfig)
{
  // for proper log files writing (immediate output)
  setbuf(stdout, NULL);
  
  // input tags
  _inputTagPrimaryVertex = edm::InputTag("hiSelectedVertex"); //'hiSelectedVertex' is generally used for PbPb collisions
  _inputTagCentrality = edm::InputTag("hiCentrality");
  _inputTagTrack = edm::InputTag("hiSelectedTracks");

  // read configuration parameters
  _flagRECO = 1;//iConfig.getParameter<int>("reco"); // if true, RECO level processed
  _nevents = 0; // number of processed events
  _neventsSelected = 0; // number of selected events
  edm::Service<TFileService> fs;
  _tree = fs->make<TTree>("Flow", "Flow"); //make output tree

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>> tree branches >>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //
  // event
  _tree->Branch("evRunNumber", &_evRunNumber, "evRunNumber/I"); // run number
  _tree->Branch("evEventNumber", &_evEventNumber, "evEventNumber/I"); // event number

  if(_flagRECO)
  {
    // primary vertex
    _tree->Branch("Npv", &_Npv, "Npv/I"); // total number of primary vertices
    _tree->Branch("pvNDOF", &_pvNDOF, "pvNDOF/I"); // number of degrees of freedom of the primary vertex
    _tree->Branch("pvZ", &_pvZ, "pvZ/F"); // z component of the primary vertex
    _tree->Branch("pvRho", &_pvRho, "pvRho/F"); // rho of the primary vertex (projection on transverse plane)
    // centrality
    _tree->Branch("HFsumETPlus", &_HFsumETPlus, "HFsumETPlus/F");
    _tree->Branch("HFsumETMinus", &_HFsumETMinus, "HFsumETMinus/F");
    _tree->Branch("HFsumET", &_HFsumET, "HFsumET/F");
    // tracks
    _tree->Branch("Ntrk", &_Ntrk, "Ntrk/I");
    _tree->Branch("trkPt", _trkPt, "trkPt[Ntrk]/F"); 
    _tree->Branch("trkEta", _trkEta, "trkEta[Ntrk]/F"); 
    _tree->Branch("trkPhi", _trkPhi, "trkPhi[Ntrk]/F");
    _tree->Branch("trkPtRes", _trkPtRes, "trkPtRes[Ntrk]/F");
    _tree->Branch("trkDzSig", _trkDzSig, "trkDzSig[Ntrk]/F");
    _tree->Branch("trkDxySig", _trkDxySig, "trkDxySig[Ntrk]/F");
    _tree->Branch("trkNpixLayers", _trkNpixLayers, "trkNpixLayers[Ntrk]/F");
    ///qinv
    /*_tree->Branch("NSSpair", &_NSSpair, "NSSpair/I");
    _tree->Branch("NOSpair", &_NOSpair, "NOSpair/I");  
    _tree->Branch("NpairMix", &_NpairMix, "NpairMix/I");
    _tree->Branch("qinvSigSS", _qinvSigSS, "qinvSigSS[NSSpair]/F");
    _tree->Branch("coulombWSS", _coulombWSS, "coulombWSS[NSSpair]/F");
    _tree->Branch("qinvSigOS", _qinvSigOS, "qinvSigOS[NOSpair]/F");
    _tree->Branch("coulombWOS", _coulombWOS, "coulombWOS[NOSpair]/F");
    _tree->Branch("qinvMixSS_OS", _qinvMixSS_OS, "qinvMixSS_OS[NpairMix]/F");
     */
  }

}


// destructor
AnalyzerFlow::~AnalyzerFlow()
{
}


//
// member functions
//

// initialise event variables with needed default (zero) values; called in the beginning of each event
void AnalyzerFlow::InitBranchVars()
{
  _evRunNumber = 0;
  _evEventNumber = 0;
  _Npv= 0;
  _pvNDOF = 0;
  _pvZ = 0;
  _pvRho = 0;

  _HFsumETPlus=-999.9;
  _HFsumETMinus=-999.9;
  _HFsumET=-999.9;

  _Ntrk = 0;
  //_NSSpair = 0;
  //_NOSpair = 0;

}

// Store event info (fill corresponding tree variables)
int AnalyzerFlow::SelectEvent(const edm::Event& iEvent)
{
  _evRunNumber = iEvent.id().run();
  _evEventNumber = iEvent.id().event();
  return 0;
}

// select primary vertex
int AnalyzerFlow::SelectPrimaryVertex(const edm::Handle<reco::VertexCollection>& primVertex)
{
  // if no primary vertices in the event, return false status
  if(primVertex->size() == 0)
    return false;
  // take the first primary vertex
  reco::VertexCollection::const_iterator pv = primVertex->begin();
  // fill z and rho (projection on transverse plane)
  _pvZ = pv->z();
  _pvRho = TMath::Sqrt(TMath::Power(pv->x(), 2.0) + TMath::Power(pv->y(), 2.0));
  // fill number of primary veritces
  _Npv = primVertex->size();
  // fill number of degrees of freedom
  _pvNDOF = pv->ndof();
  // return true status
  return true;
}

// ------------ method called for each event  ------------
void AnalyzerFlow::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;

  // event counting, printout after each 1K processed events
  _nevents++;
  //printf("*** EVENT %6d ***\n", _nevents);
  if( (_nevents % 1000) == 0)
  {
    //printf("*****************************************************************\n");
    printf("************* NEVENTS = %d K, selected = %d *************\n", _nevents / 1000, _neventsSelected);
    //printf("*****************************************************************\n");
  }
  
  // declare event contents
  edm::Handle<reco::VertexCollection> primVertex;
  edm::Handle<reco::Centrality> cent;
  edm::Handle<reco::TrackCollection> tracks;

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>> event selection >>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //
  // initialise event variables with default values
  InitBranchVars();
  // process generator level, if needed
  // process reco level, if needed
  if(_flagRECO)
  {
    // centrality 
    iEvent.getByLabel(_inputTagCentrality, cent);
    _HFsumETPlus = cent->EtHFtowerSumPlus(); //transverse energy sum in the positive side of HF detector
    _HFsumETMinus = cent->EtHFtowerSumMinus(); //transverse energy sum in the negative side of HF detector
    _HFsumET = cent->EtHFtowerSum(); //transverse energy sum in both sides of HF detector
    //std::cout<<"_HFsumETPlus : "<<_HFsumETPlus<<"; _HFsumETMinus : "<<_HFsumETMinus<<"; _HFsumET : "<<_HFsumET<<std::endl;

    // primary vertex
    iEvent.getByLabel(_inputTagPrimaryVertex, primVertex);
    reco::VertexCollection::const_iterator pv = primVertex->begin();
    // fill primary vertex
    SelectPrimaryVertex(primVertex);

    // tracks
    iEvent.getByLabel(_inputTagTrack, tracks);

    ///vectors to compute Qinv in the same event
    //std::vector<TLorentzVector> GoodTrackFourVector;
    //std::vector<int> GoodTrackCharge;

    double vtx_x = (double)primVertex->begin()->position().x(); 
    double vtx_y = (double)primVertex->begin()->position().y();
    double vtx_z = (double)primVertex->begin()->position().z();
    double vtx_xError = (double)primVertex->begin()->xError();
    double vtx_yError = (double)primVertex->begin()->yError();
    double vtx_zError = (double)primVertex->begin()->zError();
    math::XYZPoint vtx(vtx_x,vtx_y,vtx_z);

    for(reco::TrackCollection::const_iterator iter_tk = tracks->begin(); iter_tk != tracks->end(); iter_tk++){
       double aux_tk_dz_vtx = (double)iter_tk->dz(vtx);
       double aux_tk_dzError_vtx  = (double)sqrt(iter_tk->dzError()*iter_tk->dzError()+vtx_zError*vtx_zError);
       double aux_tk_dxy_vtx = (double)iter_tk->dxy(vtx);
       double aux_tk_dxyError_vtx  = (double)sqrt(iter_tk->dxyError()*iter_tk->dxyError()+vtx_xError*vtx_yError);
       const reco::HitPattern& hit_pattern = iter_tk->hitPattern();
       if(iter_tk->pt()<0.3)continue;
       if(fabs(iter_tk->eta())>2.4)continue;
       if(!iter_tk->quality(reco::TrackBase::highPurity))continue;
       if(fabs(iter_tk->ptError())/iter_tk->pt()>0.1)continue;
       if(fabs(aux_tk_dz_vtx/aux_tk_dzError_vtx)>3)continue;
       if(fabs(aux_tk_dxy_vtx/aux_tk_dxyError_vtx)>3)continue;
       if(hit_pattern.pixelLayersWithMeasurement()==0)continue;

       /*TLorentzVector pvector;
       pvector.SetXYZM(iter_tk->px(),iter_tk->py(),iter_tk->pz(),0.1396);
       GoodTrackFourVector.push_back(pvector);
       GoodTrackCharge.push_back(iter_tk->charge());*/

       _trkPt[_Ntrk] = iter_tk->pt();
       _trkEta[_Ntrk] = iter_tk->eta(); 
       _trkPhi[_Ntrk] = iter_tk->phi();
       _trkPtRes[_Ntrk] = fabs(iter_tk->ptError())/iter_tk->pt();
       _trkDzSig[_Ntrk] = aux_tk_dz_vtx/aux_tk_dzError_vtx;
       _trkDxySig[_Ntrk] = aux_tk_dxy_vtx/aux_tk_dxyError_vtx;
       _trkNpixLayers[_Ntrk] = hit_pattern.pixelLayersWithMeasurement();

       _Ntrk++;
    } 

    /*if(GoodTrackFourVector.size()<2)return;
    for(unsigned int itk1=0; itk1<GoodTrackFourVector.size();itk1++){
       for(unsigned int itk2=itk1+1; itk2<GoodTrackFourVector.size();itk2++){
          if(splitcomb(GoodTrackFourVector[itk1],  GoodTrackFourVector[itk2])){continue;}
          Double_t q = GetQ(GoodTrackFourVector[itk1], GoodTrackFourVector[itk2]);
	  if(q>1.0)continue;//to reduce size of trees
          if(GoodTrackCharge[itk1]*GoodTrackCharge[itk2]>0){ //same charge
             _qinvSigSS[_NSSpair] = q;
	     _coulombWSS[_NSSpair] = CoulombW(q);
	     _NSSpair++;
	  }else{ //oposite charge
	     _qinvSigOS[_NOSpair] = q;
             _coulombWOS[_NOSpair] = CoulombWpm(q);	     
	     _NOSpair++;
	  }      
       }	       
    }	     

    ev_vtx_z_vec.push_back(vtx_z);
    ev_GoodTrackFourVector_vec.push_back(GoodTrackFourVector);
    ev_GoodTrackCharge_vec.push_back(GoodTrackCharge);
    ev_hfsum_vec.push_back(_HFsumET);*/
  }
  // fill event info
  SelectEvent(iEvent);
  // all done: store event
  _tree->Fill();
  _neventsSelected++;
}

/*
bool AnalyzerFlow::splitcomb(TLorentzVector &vec1,TLorentzVector &vec2){
   bool issplit=false;
   Double_t cosa = TMath::Abs(vec1.Px()*vec2.Px() + vec1.Py()*vec2.Py() + vec1.Pz()*vec2.Pz())/(vec1.P()*vec2.P());
   Double_t deltapt = TMath::Abs(vec1.Pt() - vec2.Pt());
   if ( (cosa >0.99996) && (deltapt < 0.04) ) { issplit = true;}
   //std::cout << "cosa: " << cosa << " dpt: " << deltapt << " is split: " << issplit << std::endl;
   return issplit;
}

Double_t AnalyzerFlow::GetQ(const TLorentzVector &p1, const TLorentzVector &p2){
   TLorentzVector Sum4V = p1+p2;
   Double_t q = Sum4V.Mag2() - 4*0.1396*0.1396;
   //  std::cout<<(  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  ) <<std::endl;
   return (  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}

//return the weight factor due to Coloumb repulsion [Gamow] same charge
const Double_t AnalyzerFlow::CoulombW(const Double_t& q){
   const Double_t alpha=1./137.;
   Double_t x=2*TMath::Pi()*(alpha*0.1396/q);
   //return (TMath::Exp(x)-1.)/x; // OLD MATTIA's DEFINITION

   //Double_t ws=scf*((exp(arg)-1)/arg-1)+1; // PAOLO's DEFINITION
   Double_t weight = 1;//0.85; // TEMPORARY SET TO 0.85 * GAMOW FACTOR
   //Double_t weight = 1.15; //for syst. +15%
   //Double_t weight = 0.85; //for syst. -15%
   return weight*( (TMath::Exp(x)-1.)/x -1 ) + 1;
}

//return the weight factor due to Coloumb attraction [Gamow] opposite charge
const  Double_t AnalyzerFlow::CoulombWpm(const Double_t& q){
   const Double_t alpha=1./137.;
   Double_t x=2*TMath::Pi()*(alpha*0.1396/q);
   // return (1.-TMath::Exp(-x))/x; // OLD MATTIA's DEFINITION

   // Double_t wd=scf*((1-exp(-arg))/arg-1)+1; // PAOLO's DEFINITION
   Double_t weight = 1;//0.85; // TEMPORARY SET TO 0.85 * GAMOW FACTOR
   //Double_t weight = 1.15; //for syst. +15%
   //Double_t weight = 0.85; //for syst. -15%
   return weight*( (1.-TMath::Exp(-x))/x -1 ) + 1;
}

*/

// ------------ method called when starting to processes a run  ------------
void AnalyzerFlow::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

// below is some default stuff, was not modified

// ------------ method called once each job just before starting event loop  ------------
void AnalyzerFlow::beginJob() {;}

// ------------ method called once each job just after ending the event loop  ------------
void AnalyzerFlow::endJob() {;}

// ------------ method called when ending the processing of a run  ------------
void AnalyzerFlow::endRun(edm::Run const& run, edm::EventSetup const& setup) {;}

// ------------ method called when starting to processes a luminosity block  ------------
void AnalyzerFlow::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {;}

// ------------ method called when ending the processing of a luminosity block  ------------
void AnalyzerFlow::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {;}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void AnalyzerFlow::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzerFlow);
