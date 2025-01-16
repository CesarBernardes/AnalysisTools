
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

// for muons
///#include "DataFormats/MuonReco/interface/Muon.h"
///#include "DataFormats/MuonReco/interface/MuonFwd.h"
///#include "DataFormats/MuonReco/interface/MuonSelectors.h"
///#include "DataFormats/MuonReco/interface/MuonIsolation.h"

//for centrality
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

// for electrons uncomment when implemented
//#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

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

// Trigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include <cassert>

//
// class declaration
//
class AnalyzerHBT : public edm::EDAnalyzer {
   public:
      explicit AnalyzerHBT(const edm::ParameterSet&);
      ~AnalyzerHBT();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual int analyzeTrigger(const edm::Event&, const edm::EventSetup&, const std::string& triggerName);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      
      // user routines (detailed description given with the method implementations)
      int SelectEvent(const edm::Event& iEvent);
      ///int SelectMu(const edm::Handle<reco::TrackCollection>& muons, const reco::VertexCollection::const_iterator& pv);
// Note that muons are taken from the TrackCollection one can collect other necessary data from the input root file by looking
// at its structure in the TBrowser:
// i.e. from the TBrowser we see a folder called:  recoTracks_globalMuons__RECO.
// This means using edm::Handle<reco::TrackCollection> is needed to get hte data (see line just above) and globalMuons will be the 
// relative input tag used in the analyzer function.
   //   int SelectEl(const edm::Handle<reco::GsfElectronCollection>& electrons, const reco::VertexCollection::const_iterator& pv);
      int SelectPrimaryVertex(const edm::Handle<reco::VertexCollection>& primVertex);
      ///const reco::Candidate* GetFinalState(const reco::Candidate* particle, const int id);
      ///void FillFourMomentum(const reco::Candidate* particle, float* p);
      void InitBranchVars();
      bool splitcomb(TLorentzVector &vec1,TLorentzVector &vec2);
      Double_t GetQ( const TLorentzVector &p1, const TLorentzVector &p2 );
      const Double_t CoulombWpm( const Double_t& q);
      const Double_t CoulombW( const Double_t& q);
      ///void MixEvents(float hfsum_min, float hfsum_max, int nEvt_to_mix);

      //define some constants
      ///static const double pi_mass  = 0.1396;
      ///static const double cos_cut = 0.99996;
      ///static const double dpt_cut = 0.04;

      // input tags
      ///edm::InputTag _inputTagMuons;
      ///edm::InputTag _inputTagElectrons;
      ///edm::InputTag _inputTagBtags;
      edm::InputTag _inputTagPrimaryVertex;
      edm::InputTag _inputTagCentrality;
      edm::InputTag _inputTagTrack;

      //trigger
      std::string   processName_;
      std::string   triggerName_;
      std::string   datasetName_;
      edm::InputTag triggerResultsTag_;
      edm::InputTag triggerEventTag_;
       // additional class data memebers 
       // these are actually the containers where we will store
       // the trigger information
      edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
      edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
      HLTConfigProvider hltConfig_;
      int* trigflag;
      int HltEvtCnt;

      // general flags and variables
      int _flagMC;
      int _flagRECO;
      int _flagGEN;
      int _nevents;
      int _neventsSelected;
      ///int _signLeptonP;
      ///int _signLeptonM;
      
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
      // muons
      ///static const int _maxNmu = 10;
      ///int _Nmu;
      ///int _Nmu0;
      ///float _muPt[_maxNmu];
      ///float _muEta[_maxNmu];
      ///float _muPhi[_maxNmu];
      ///float _muC[_maxNmu];
      ///float _muIso03[_maxNmu];
      ///float _muIso04[_maxNmu];
      ///int _muHitsValid[_maxNmu];
      ///int _muHitsPixel[_maxNmu];
      ///float _muDistPV0[_maxNmu];
      ///float _muDistPVz[_maxNmu];
      ///float _muTrackChi2NDOF[_maxNmu];
      // electrons
      ///static const int _maxNel = 10;
      ///int _Nel;
      ///float _elPt[_maxNel];
      ///float _elEta[_maxNel];
      ///float _elPhi[_maxNel];
      ///float _elIso03[_maxNel];
      ///float _elIso04[_maxNel];
      ///int _elConvFlag[_maxNel];
      ///float _elConvDist[_maxNel];
      ///float _elConvDcot[_maxNel];
      ///float _elMissHits[_maxNel];
      ///float _elDistPV0[_maxNel];
      ///float _elDistPVz[_maxNel];
      // primary vertex
      int _Npv;
      int _pvNDOF;
      float _pvZ;
      float _pvRho;
      std::vector<double> ev_vtx_z_vec;
      std::vector< std::vector<TLorentzVector> > ev_GoodTrackFourVector_vec;
      std::vector< std::vector<int> > ev_GoodTrackCharge_vec;
      std::vector<float> ev_hfsum_vec;
      // centrality
      float _HFsumETPlus;
      float _HFsumETMinus;
      float _HFsumET;
      float _zdcSum;
      float _zdcSumMinus;
      float _zdcSumPlus;
      float _pixelMultiplicity;
      // tracks
      static const int _maxNtrk = 20000;
      int _Ntrk;
      float _trkPt[_maxNtrk];
      float _trkEta[_maxNtrk];
      float _trkPhi[_maxNtrk];
      float _trkPtRes[_maxNtrk];
      float _trkDzSig[_maxNtrk];
      float _trkDxySig[_maxNtrk];
      float _trkNpixLayers[_maxNtrk];
      float _trkNstripLayers[_maxNtrk];
      float _trkNvalidHits[_maxNtrk];
      float _trkNormChi2[_maxNtrk];

      ///qinv
      /////static const int _maxNpair = _maxNtrk*_maxNtrk;
      /////int _NSSpair;
      /////int _NOSpair;
      ///int _NpairMix;
      /////float _qinvSigSS[_maxNpair];
      /////float _coulombWSS[_maxNpair];
      /////float _qinvSigOS[_maxNpair];
      /////float _coulombWOS[_maxNpair];
      ///float _qinvMixSS_OS[_maxNpair]; //for mixing
      ///float _HFsumETMixP1[_maxNpair]; //for mixing (particle 1 in the pair) - to map correctly the qinv obtained in the mixing procedure
      ///float _HFsumETMixP2[_maxNpair]; //for mixing (particle 2 in the pair) 	

};

//
// constants (particle masses)
//
///double _massMu = 0.105658;
///double _massEl = 0.000511;

//
// constructor
//
AnalyzerHBT::AnalyzerHBT(const edm::ParameterSet& iConfig):
processName_(iConfig.getParameter<std::string>("processName")),
triggerName_(iConfig.getParameter<std::string>("triggerName")),
datasetName_(iConfig.getParameter<std::string>("datasetName")),
triggerResultsTag_(iConfig.getParameter<edm::InputTag>("triggerResults")),
triggerEventTag_(iConfig.getParameter<edm::InputTag>("triggerEvent"))
{
  // for proper log files writing (immediate output)
  setbuf(stdout, NULL);
  
  // input tags
  ///_inputTagMuons = edm::InputTag("globalMuons");
  //_inputTagElectrons = edm::InputTag("gsfElectrons"); //use this to Analyze electrons
  ///_inputTagPrimaryVertex = edm::InputTag("offlinePrimaryVertices"); //'hiSelectedVertex' is generally used for PbPb collisions
   _inputTagPrimaryVertex = edm::InputTag("hiSelectedVertex"); 
  _inputTagCentrality = edm::InputTag("hiCentrality");
///  _inputTagTrack = edm::InputTag("generalTracks");
  _inputTagTrack = edm::InputTag("hiSelectedTracks");

  // read configuration parameters
  _flagMC = 0;//iConfig.getParameter<int>("mc"); // true for MC, false for data
  _flagRECO = 1;//iConfig.getParameter<int>("reco"); // if true, RECO level processed
  _flagGEN = 0;//iConfig.getParameter<int>("gen"); // if true, generator level processed (works only for MC)
  _nevents = 0; // number of processed events
  _neventsSelected = 0; // number of selected events
  edm::Service<TFileService> fs;
  _tree = fs->make<TTree>("TreeMBUCC", "TreeMBUCC"); //make output tree

  //trigger
  HltEvtCnt = 0;
  const int kMaxTrigFlag = 10000;
  trigflag = new int[kMaxTrigFlag];  

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>> tree branches >>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //
  // event
  _tree->Branch("evRunNumber", &_evRunNumber, "evRunNumber/I"); // run number
  _tree->Branch("evEventNumber", &_evEventNumber, "evEventNumber/I"); // event number

  if(_flagRECO)
  {
    // muons
    ///_tree->Branch("Nmu", &_Nmu, "Nmu/I"); // number of muons 
    ///_tree->Branch("muPt", _muPt, "muPt[Nmu]/F"); // muon pT
    ///_tree->Branch("muEta", _muEta, "muEta[Nmu]/F"); // muon pseudorapidity
    ///_tree->Branch("muPhi", _muPhi, "muPhi[Nmu]/F"); // muon phi
    ///_tree->Branch("muC", _muC, "muC[Nmu]/F"); // muon phi
    ///_tree->Branch("muIso03", _muIso03, "muIso03[Nmu]/F"); // muon isolation, delta_R=0.3
    ///_tree->Branch("muIso04", _muIso04, "muIso04[Nmu]/F"); // muon isolation, delta_R=0.4
    ///_tree->Branch("muHitsValid", _muHitsValid, "muHitsValid[Nmu]/I"); // muon valid hits number
    ///_tree->Branch("muHitsPixel", _muHitsPixel, "muHitsPixel[Nmu]/I"); // muon pixel hits number
    ///_tree->Branch("muDistPV0", _muDistPV0, "muDistPV0[Nmu]/F"); // muon distance to the primary vertex (projection on transverse plane)
    ///_tree->Branch("muDistPVz", _muDistPVz, "muDistPVz[Nmu]/F"); // muon distance to the primary vertex (z projection)
    ///_tree->Branch("muTrackChi2NDOF", _muTrackChi2NDOF, "muTrackChi2NDOF[Nmu]/F"); // muon track number of degrees of freedom
    // primary vertex
    _tree->Branch("Npv", &_Npv, "Npv/I"); // total number of primary vertices
    _tree->Branch("pvNDOF", &_pvNDOF, "pvNDOF/I"); // number of degrees of freedom of the primary vertex
    _tree->Branch("pvZ", &_pvZ, "pvZ/F"); // z component of the primary vertex
    _tree->Branch("pvRho", &_pvRho, "pvRho/F"); // rho of the primary vertex (projection on transverse plane)
    // centrality
    _tree->Branch("HFsumETPlus", &_HFsumETPlus, "HFsumETPlus/F");
    _tree->Branch("HFsumETMinus", &_HFsumETMinus, "HFsumETMinus/F");
    _tree->Branch("HFsumET", &_HFsumET, "HFsumET/F");
    _tree->Branch("zdcSum", &_zdcSum, "zdcSum/F");
    _tree->Branch("zdcSumMinus", &_zdcSumMinus, "zdcSumMinus/F");
    _tree->Branch("zdcSumPlus", &_zdcSumPlus, "zdcSumPlus/F");
    _tree->Branch("pixelMultiplicity", &_pixelMultiplicity, "pixelMultiplicity/F");
    // tracks
    _tree->Branch("Ntrk", &_Ntrk, "Ntrk/I");
    _tree->Branch("trkPt", _trkPt, "trkPt[Ntrk]/F"); 
    _tree->Branch("trkEta", _trkEta, "trkEta[Ntrk]/F"); 
    _tree->Branch("trkPhi", _trkPhi, "trkPhi[Ntrk]/F");
    _tree->Branch("trkPtRes", _trkPtRes, "trkPtRes[Ntrk]/F");
    _tree->Branch("trkDzSig", _trkDzSig, "trkDzSig[Ntrk]/F");
    _tree->Branch("trkDxySig", _trkDxySig, "trkDxySig[Ntrk]/F");
    _tree->Branch("trkNpixLayers", _trkNpixLayers, "trkNpixLayers[Ntrk]/F");
    _tree->Branch("trkNstripLayers", _trkNstripLayers, "trkNstripLayers[Ntrk]/F");
    _tree->Branch("trkNvalidHits", _trkNvalidHits, "trkNvalidHits[Ntrk]/F");
    _tree->Branch("trkNormChi2",_trkNormChi2, "trkNormChi2[Ntrk]/F");

    ///qinv
    /////_tree->Branch("NSSpair", &_NSSpair, "NSSpair/I");
    /////_tree->Branch("NOSpair", &_NOSpair, "NOSpair/I");  
    ///_tree->Branch("NpairMix", &_NpairMix, "NpairMix/I");
    /////_tree->Branch("qinvSigSS", _qinvSigSS, "qinvSigSS[NSSpair]/F");
    /////_tree->Branch("coulombWSS", _coulombWSS, "coulombWSS[NSSpair]/F");
    /////_tree->Branch("qinvSigOS", _qinvSigOS, "qinvSigOS[NOSpair]/F");
    /////_tree->Branch("coulombWOS", _coulombWOS, "coulombWOS[NOSpair]/F");
    ///_tree->Branch("qinvMixSS_OS", _qinvMixSS_OS, "qinvMixSS_OS[NpairMix]/F");
    ///_tree->Branch("HFsumETMixP1", _HFsumETMixP1, "HFsumETMixP1[NpairMix]/F");
    ///_tree->Branch("HFsumETMixP2", _HFsumETMixP2, "HFsumETMixP2[NpairMix]/F");

  }

}


// destructor
AnalyzerHBT::~AnalyzerHBT()
{
}


//
// member functions
//

// initialise event variables with needed default (zero) values; called in the beginning of each event
void AnalyzerHBT::InitBranchVars()
{
  _evRunNumber = 0;
  _evEventNumber = 0;
  ///_Nmu = 0;
  _Npv= 0;
  _pvNDOF = 0;
  _pvZ = 0;
  _pvRho = 0;

  _HFsumETPlus=-9999.9;
  _HFsumETMinus=-9999.9;
  _HFsumET=-9999.9;
  _zdcSum=-9999.9;
  _zdcSumMinus=-9999.9;
  _zdcSumPlus=-9999.9;
  _pixelMultiplicity=-9999.9;

  _Ntrk = 0;
  //_NSSpair = 0;
  //_NOSpair = 0;

}

// Store event info (fill corresponding tree variables)
int AnalyzerHBT::SelectEvent(const edm::Event& iEvent)
{
  _evRunNumber = iEvent.id().run();
  _evEventNumber = iEvent.id().event();
  return 0;
}

// muon selection
///int AnalyzerHBT::SelectMu(const edm::Handle<reco::TrackCollection>& muons, const reco::VertexCollection::const_iterator& pv)
///{
///  using namespace std;
///  _Nmu = 0;
///  _Nmu0 = 0;
///  // loop over muons
///  for (reco::TrackCollection::const_iterator it = muons->begin(); it != muons->end(); it++)
///  {
///	_Nmu0++;
///    if(_Nmu == _maxNmu)
///    {
///      printf("Maximum number of muons %d reached, skipping the rest\n", _maxNmu);
///      return 0;
///    }
///    _muHitsValid[_Nmu] = 0;
///    _muHitsPixel[_Nmu] = 0;
///    const reco::HitPattern& p = it->hitPattern();
///    for (int i = 0; i < p.numberOfHits(); i++) 
///    {
///      uint32_t hit = p.getHitPattern(i);
///      if (p.validHitFilter(hit) && p.pixelHitFilter(hit))
///        _muHitsPixel[_Nmu]++;
///      if (p.validHitFilter(hit))
///        _muHitsValid[_Nmu]++;
///    }
///    // fill three momentum (pT, eta, phi)
///    _muPt[_Nmu] = it->pt();// * it->charge();
///    _muEta[_Nmu] = it->eta();
///    _muPhi[_Nmu] = it->phi();
///    _muC[_Nmu]=it->charge();
///    // fill chi2/ndof
///    if (it->ndof()) _muTrackChi2NDOF[_Nmu] = it->chi2() / it->ndof();
///    // fill distance to primary vertex
///    _muDistPV0[_Nmu] = TMath::Sqrt(TMath::Power(pv->x() - it->vx(), 2.0) + TMath::Power(pv->y() - it->vy(), 2.0));
///    _muDistPVz[_Nmu] = TMath::Abs(pv->z() - it->vz());
///    // store muon
///    _Nmu++;
///    // determine muon sign (in the end the event will be stored only there are opposite signed leptons)
///    if(it->charge() == +1)
///        _signLeptonP = 1;
///    if(it->charge() == -1)
///        _signLeptonM = 1;
///  }
///	cout<<"Muons before selection: "<<_Nmu0<<endl;
///  return 0;
///}

// select primary vertex
int AnalyzerHBT::SelectPrimaryVertex(const edm::Handle<reco::VertexCollection>& primVertex)
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

// fill 4-momentum (p) with provided particle pointer
///void AnalyzerHBT::FillFourMomentum(const reco::Candidate* particle, float* p)
///{
///  // if NULL pointer provided, initialise with default (zero)
///  if(particle == NULL)
///  {
///    p[0] = p[1] = p[2] = p[3] = 0.0;
///    return;
///  }
///  
///  p[0] = particle->px();
///  p[1] = particle->py();
///  p[2] = particle->pz();
///  p[3] = particle->mass();
///}

// select MC generator level information
// (analysis specific ttbar dileptonic decay)

// ------------ method called for each event  ------------
void AnalyzerHBT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;

   //Get event products: 
   // In the following, the code is trying to access the information 
   // from the ROOT files and point the containers (that we created), 
   // namely triggerResultsHandle_ aed triggerEVentHandle_, 
   // to the correct "address", given at configuration time 
   // and assigned to triggerResultsTag_ and triggerEventTag_
 
   // After that, a simple sanity check is done.
 
   iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
   if (!triggerResultsHandle_.isValid()) {
     //cout << "HLTEventAnalyzerAOD::analyze: Error in getting TriggerResults product from Event!" << endl;
     return;
   }
   iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
   if (!triggerEventHandle_.isValid()) {
     //cout << "HLTEventAnalyzerAOD::analyze: Error in getting TriggerEvent product from Event!" << endl;
     return;
   }

   // sanity check
   assert(triggerResultsHandle_->size()==hltConfig_.size());

   bool pass_trigger=false;
   Int_t accept;
   TString trigName;
   const vector<string> triggerNamesInDS = hltConfig_.datasetContent(datasetName_);
   //cout << "trigger Names size: "<<triggerNamesInDS.size()<<endl;
   for (unsigned i = 0; i < triggerNamesInDS.size(); i++) {
      accept=0;
      trigName = triggerNamesInDS[i];
      trigflag[i]=0;
      if (HltEvtCnt==0){
         if(trigName=="HLT_HIUCC015_v1" || trigName=="HLT_HIUCC015_v2" || trigName=="HLT_HIUCC015_v3" || trigName=="HLT_HIUCC015_v4" || trigName=="HLT_HIUCC010_v1" || trigName=="HLT_HIUCC010_v2" || trigName=="HLT_HIUCC010_v3" || trigName=="HLT_HIUCC010_v4" || trigName=="HLT_HIMinBiasHfOrBSC_v1" || trigName=="HLT_HIMinBiasHfOrBSC_v2" || trigName=="HLT_HIMinBiasHfOrBSC_v3" || trigName=="HLT_HIMinBiasHfOrBSC_v4") _tree->Branch(trigName,&trigflag[i],trigName+"/I");
      }else{
	if (!_tree->GetBranch(trigName) && (trigName=="HLT_HIUCC015_v1" || trigName=="HLT_HIUCC015_v2" || trigName=="HLT_HIUCC015_v3" || trigName=="HLT_HIUCC015_v4" || trigName=="HLT_HIUCC010_v1" || trigName=="HLT_HIUCC010_v2" || trigName=="HLT_HIUCC010_v3" || trigName=="HLT_HIUCC010_v4" || trigName=="HLT_HIMinBiasHfOrBSC_v1" || trigName=="HLT_HIMinBiasHfOrBSC_v2" || trigName=="HLT_HIMinBiasHfOrBSC_v3" || trigName=="HLT_HIMinBiasHfOrBSC_v4")) _tree->Branch(trigName,&trigflag[i],trigName+"/I");
     }
     accept=  analyzeTrigger(iEvent,iSetup,triggerNamesInDS[i]);
     if(accept && (trigName=="HLT_HIUCC015_v1" || trigName=="HLT_HIUCC015_v2" || trigName=="HLT_HIUCC015_v3" || trigName=="HLT_HIUCC015_v4" || trigName=="HLT_HIUCC010_v1" || trigName=="HLT_HIUCC010_v2" || trigName=="HLT_HIUCC010_v3" || trigName=="HLT_HIUCC010_v4" || trigName=="HLT_HIMinBiasHfOrBSC_v1" || trigName=="HLT_HIMinBiasHfOrBSC_v2" || trigName=="HLT_HIMinBiasHfOrBSC_v3" || trigName=="HLT_HIMinBiasHfOrBSC_v4")){
        trigflag[i] = 1;
	pass_trigger=true;
	//cout << trigName<<endl;
     }
     //if((trigName=="HLT_HIMinBiasHfOrBSC_v1" || trigName=="HLT_HIMinBiasHfOrBSC_v2" || trigName=="HLT_HIMinBiasHfOrBSC_v3" || trigName=="HLT_HIMinBiasHfOrBSC_v4") && accept){ pass_trigger=true; break; }
     //if((trigName=="HLT_HIMinBiasHfOrBSC_v1" || trigName=="HLT_HIMinBiasHfOrBSC_v2" || trigName=="HLT_HIMinBiasHfOrBSC_v3" || trigName=="HLT_HIMinBiasHfOrBSC_v4") && !accept){ pass_trigger=false; break; }
     //if(trigName=="HLT_HIMinBiasHfOrBSC_v1" && accept){ pass_trigger=true; break; }
     //else if(trigName=="HLT_HIMinBiasHfOrBSC_v1" && !accept){ pass_trigger=false; break; }
     //else {}
   }
   HltEvtCnt++;

if(pass_trigger){ 

  // event counting, printout after each 1K processed events
  _nevents++;
  //printf("*** EVENT %6d ***\n", _nevents);
  if( (_nevents % 1000) == 0)
  {
    //printf("*****************************************************************\n");
    printf("************* NEVENTS = %d K, selected = %d *************\n", _nevents / 1000, _neventsSelected);
    //printf("*****************************************************************\n");
  }
  //return;
  
  // declare event contents
  edm::Handle<reco::VertexCollection> primVertex;
  ///edm::Handle<reco::TrackCollection> muons;
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
    _HFsumETPlus = cent->EtHFtowerSumPlus();
    _HFsumETMinus = cent->EtHFtowerSumMinus();
    _HFsumET = cent->EtHFtowerSum();
    _zdcSum = cent->zdcSum();
    //std::cout<<"_HFsumETPlus : "<<_HFsumETPlus<<"; _HFsumETMinus : "<<_HFsumETMinus<<"; _HFsumET : "<<_HFsumET<<std::endl;
    _zdcSumMinus = cent->zdcSumMinus();
    _zdcSumPlus = cent->zdcSumPlus();
    _pixelMultiplicity = cent->multiplicityPixel();

    // primary vertex
    iEvent.getByLabel(_inputTagPrimaryVertex, primVertex);
    reco::VertexCollection::const_iterator pv = primVertex->begin();
    // muons
    ///iEvent.getByLabel(_inputTagMuons, muons);
    ///SelectMu(muons, pv);
    // fill primary vertex
    SelectPrimaryVertex(primVertex);

    if(7.*_HFsumET+_zdcSum>=36000) {return;} //PU rejection

    // tracks
    iEvent.getByLabel(_inputTagTrack, tracks);

    ///vectors to compute Qinv in the same event
    std::vector<TLorentzVector> GoodTrackFourVector;
    std::vector<int> GoodTrackCharge;

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
       //if(iter_tk->pt()<0.3)continue;
       //if(iter_tk->pt()<0.4)continue;
       if(fabs(iter_tk->eta())>2.4)continue;
       if(!iter_tk->quality(reco::TrackBase::highPurity))continue;
       //if(fabs(iter_tk->ptError())/iter_tk->pt()>0.05)continue;
       //if(fabs(iter_tk->ptError())/iter_tk->pt()>0.15)continue;
       //if(fabs(aux_tk_dz_vtx/aux_tk_dzError_vtx)>5)continue;
       //if(fabs(aux_tk_dxy_vtx/aux_tk_dxyError_vtx)>5)continue;
       //if(hit_pattern.pixelLayersWithMeasurement()==0)continue;
       //if(iter_tk->numberOfValidHits()<13)continue;
       //if(iter_tk->normalizedChi2()>0.15*iter_tk->numberOfValidHits())continue;

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
       _trkNstripLayers[_Ntrk] = hit_pattern.stripLayersWithMeasurement();
       _trkNvalidHits[_Ntrk] = hit_pattern.numberOfValidHits();
       _trkNormChi2[_Ntrk] = iter_tk->normalizedChi2();
 
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

}///end - if pass trigger

}


bool AnalyzerHBT::splitcomb(TLorentzVector &vec1,TLorentzVector &vec2){
   bool issplit=false;
   Double_t cosa = TMath::Abs(vec1.Px()*vec2.Px() + vec1.Py()*vec2.Py() + vec1.Pz()*vec2.Pz())/(vec1.P()*vec2.P());
   Double_t deltapt = TMath::Abs(vec1.Pt() - vec2.Pt());
   if ( (cosa >0.99996) && (deltapt < 0.04) ) { issplit = true;}
   //std::cout << "cosa: " << cosa << " dpt: " << deltapt << " is split: " << issplit << std::endl;
   return issplit;
}

Double_t AnalyzerHBT::GetQ(const TLorentzVector &p1, const TLorentzVector &p2){
   TLorentzVector Sum4V = p1+p2;
   Double_t q = Sum4V.Mag2() - 4*0.1396*0.1396;
   //  std::cout<<(  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  ) <<std::endl;
   return (  q>0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}

//return the weight factor due to Coloumb repulsion [Gamow] same charge
const Double_t AnalyzerHBT::CoulombW(const Double_t& q){
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
const  Double_t AnalyzerHBT::CoulombWpm(const Double_t& q){
   const Double_t alpha=1./137.;
   Double_t x=2*TMath::Pi()*(alpha*0.1396/q);
   // return (1.-TMath::Exp(-x))/x; // OLD MATTIA's DEFINITION

   // Double_t wd=scf*((1-exp(-arg))/arg-1)+1; // PAOLO's DEFINITION
   Double_t weight = 1;//0.85; // TEMPORARY SET TO 0.85 * GAMOW FACTOR
   //Double_t weight = 1.15; //for syst. +15%
   //Double_t weight = 0.85; //for syst. -15%
   return weight*( (1.-TMath::Exp(-x))/x -1 ) + 1;
}



// ------------ method called when starting to processes a run  ------------
void AnalyzerHBT::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{

   using namespace std;
   using namespace edm;	

   //If the hltConfig can be initialized, then the below is an example of
   //how to extract the config information for the trigger from the 
   //so-called provenance.

   // The trigger configuration can change from 
   // run to run (during the run is the same),
   // so it needs to be called here.

   ///   "init" return value indicates whether intitialisation has succeeded
   ///   "changed" parameter indicates whether the config has actually changed

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      // check if trigger name in (new) config
      if (triggerName_!="@") { // "@" means: analyze all triggers in config
	const unsigned int n(hltConfig_.size());
	const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
	//if (triggerIndex>=n) {
	//  cout << "HLTEventAnalyzerAOD::analyze:"
	//       << " TriggerName " << triggerName_ 
	//       << " not available in (new) config!" << endl;
	//  cout << "Available TriggerNames are: " << endl;
	//  hltConfig_.dump("Triggers");
	//}
      }
	//hltConfig_.dump("Datasets");//use to check the Dataset name to analyze the triggers
    }
  } else {
    //cout << "HLTEventAnalyzerAOD::analyze:"
	// << " config extraction failure with process name "
	// << processName_ << endl;
  }


}

// below is some default stuff, was not modified

// ------------ method called once each job just before starting event loop  ------------
void AnalyzerHBT::beginJob() {;}

// ------------ method called once each job just after ending the event loop  ------------
void AnalyzerHBT::endJob() {

///std::cout << " This is called at the end of the each job " << std::endl;
///std::cout << " Now is the time for mixing the events "<< std::endl;

///First two define hfSumET (GeV) range to mix and third defines max number of events to mix
///MixEvents(3380.0,5000.0,10); //0-5% 
///MixEvents(2770.0,3380.0,10); //5-10% 
///MixEvents(2280.0,2770.0,10); //10-15%
///MixEvents(1850.0,2280.0,10); //15-20%
///MixEvents(1500.0,1850.0,10); //20-25%
///MixEvents(1210.0,1500.0,10); //25-30%
///MixEvents(900.0,1210.0,10); //30-35%
///MixEvents(240.0,900.0,10); //35-60%
///MixEvents(0.0,240.0,10); //60-100%

}
/*
void AnalyzerHBT::MixEvents(float hfsum_min, float hfsum_max, int nEvt_to_mix){

_NpairMix=0;
int aux_n_evts = (int)ev_hfsum_vec.size();

for(int nevt1=0; nevt1<aux_n_evts; nevt1++) {
  
   std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt1_vec= ev_GoodTrackFourVector_vec[nevt1];
   int nMixmult_nevt1=ev_GoodTrackFourVectorTemp_nevt1_vec.size();
   std::vector<int> ev_GoodTrackChargeTemp_nevt1_vec = ev_GoodTrackCharge_vec[nevt1];

   if(ev_hfsum_vec[nevt1]<hfsum_min || hfsum_max<ev_hfsum_vec[nevt1])continue;

   int takeAssociated = 0;
   for(int nevt_assoc=nevt1+1; nevt_assoc<aux_n_evts; nevt_assoc++) {
      if(ev_hfsum_vec[nevt_assoc]<hfsum_min || hfsum_max<ev_hfsum_vec[nevt_assoc])continue;
      if(std::abs((ev_vtx_z_vec)[nevt1] - (ev_vtx_z_vec)[nevt_assoc])>2.0) continue; //use only events close in z
 
      std::vector<TLorentzVector> ev_GoodTrackFourVectorTemp_nevt_assoc_vec= ev_GoodTrackFourVector_vec[nevt_assoc]; 
      int nMixmult_nevt_assoc=ev_GoodTrackFourVectorTemp_nevt_assoc_vec.size();
      std::vector<int> ev_GoodTrackChargeTemp_nevt_assoc_vec=ev_GoodTrackCharge_vec[nevt_assoc];
      
      takeAssociated++;

      if(takeAssociated > nEvt_to_mix)break; //for "nEvt_to_mix" events mixing. If you remove this "break", you will mix all the events.

      //DeltazVertex_->Fill((ev_vtx_z_vec)[nevt1] - (ev_vtx_z_vec)[nevt_assoc]);
      for(int imix=0; imix<nMixmult_nevt1; imix++){
         for(int iimix=0; iimix<nMixmult_nevt_assoc; iimix++){
            if(splitcomb(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix])){continue;}
            Double_t qmix = GetQ(ev_GoodTrackFourVectorTemp_nevt1_vec[imix], ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
	    if(qmix>1.0)continue;//to reduce size of trees
            _qinvMixSS_OS[_NpairMix] = qmix;
	    _HFsumETMixP1[_NpairMix] = ev_hfsum_vec[nevt1];
	    _HFsumETMixP2[_NpairMix] = ev_hfsum_vec[nevt_assoc];
	    _NpairMix++;
            //Double_t qmixlong = GetQlong(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            //Double_t qmixlongLCMS = GetQlongLCMS(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            //Double_t qmixout = GetQout(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);
            //Double_t qmixside = GetQside(ev_GoodTrackFourVectorTemp_nevt1_vec[imix],ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix]);

            //TLorentzVector psum2 = ev_GoodTrackFourVectorTemp_nevt1_vec[imix]+ev_GoodTrackFourVectorTemp_nevt_assoc_vec[iimix];
            //Double_t kt=(psum2.Pt())/2.; 
         }
      }
   }
}//end first for loop 
_tree->Fill();
}//end MixEvents function
*/
// ------------ method called when ending the processing of a run  ------------
void AnalyzerHBT::endRun(edm::Run const& run, edm::EventSetup const& setup) {;}

// ------------ method called when starting to processes a luminosity block  ------------
void AnalyzerHBT::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {;}

// ------------ method called when ending the processing of a luminosity block  ------------
void AnalyzerHBT::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {;}

//---------------------------Actual trigger analysis-------------
int AnalyzerHBT::analyzeTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName) 
//-----------------------------------------------------------------
{

  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;
  Int_t trig_accept=0;
  //cout<<"Currently analyzing trigger "<<triggerName<<endl;

  //Check the current configuration to see how many total triggers there are
  const unsigned int n(hltConfig_.size());
  //Get the trigger index for the current trigger
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  //check that the trigger in the event and in the configuration agree
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));
  // abort on invalid trigger name
  if (triggerIndex>=n) {
    //cout << "HLTEventAnalyzerAOD::analyzeTrigger: path "
	// << triggerName << " - not found!" << endl;
    return 0;
  }
  
//  const std::pair<int,int> prescales(hltConfig_.prescaleValues(iEvent,iSetup,triggerName));
//  cout << "HLTEventAnalyzerAOD::analyzeTrigger: path "
//      << triggerName << " [" << triggerIndex << "] "
//      << "prescales L1T,HLT: " << prescales.first << "," << prescales.second
//      << endl;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // EXAMPLE: Find out if the trigger was active, accepted, or in error.
  // We could also find out whether the trigger was active (wasrun), 
  // if it accepted the event (accept) or if it gave an error (error).
  // Results from TriggerResults product
  //Uncomment the lines below
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(triggerResultsHandle_->accept(triggerIndex)){
//cout <<triggerName<<"\t";
  //cout << " Trigger path status:"
    //   << " WasRun=" << triggerResultsHandle_->wasrun(triggerIndex)
    //   << " Accept=" << triggerResultsHandle_->accept(triggerIndex)
    //   << " Error =" << triggerResultsHandle_->error(triggerIndex)
    //   << endl;
   trig_accept = 1;
}
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  const unsigned int m(hltConfig_.size(triggerIndex));
  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  assert (moduleIndex<m);
  
  return trig_accept;
}//--------------------------analyzeTrigger()


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void AnalyzerHBT::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzerHBT);
