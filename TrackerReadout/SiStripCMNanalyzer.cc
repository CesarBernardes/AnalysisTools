// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include "DataFormats/SiStripDigi/interface/SiStripProcessedRawDigi.h"
#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"

#include "CondFormats/SiStripObjects/interface/SiStripPedestals.h"
#include "CondFormats/DataRecord/interface/SiStripPedestalsRcd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/Provenance/interface/RunLumiEventNumber.h"

//ROOT inclusion
#include "TH1F.h"


//cesar-begin
#include "CondFormats/SiStripObjects/interface/SiStripFedCabling.h"
#include "CondFormats/DataRecord/interface/SiStripFedCablingRcd.h"
#include "TProfile.h" 
#include "TH2.h"
#include "TTree.h"
//cesar-end


//
// class decleration
//

class SiStripCMNanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns> {
  public:
    explicit SiStripCMNanalyzer(const edm::ParameterSet&);
    ~SiStripCMNanalyzer() override;
    static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

  private:
    void beginJob() override ;
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endJob() override ;
    void beginRun(edm::Run const& iEvent, edm::EventSetup const&) override;
    void endRun(edm::Run const& iEvent, edm::EventSetup const&) override;

    edm::ESGetToken<SiStripFedCabling, SiStripFedCablingRcd> cablingToken_; //cesar

    edm::EDGetTokenT<edm::DetSetVector<SiStripRawDigi> > srcAPVCM_;
    edm::Service<TFileService> fs_;

    TFileDirectory sdMisc_;

    //const SiStripFedCabling* cabling_; //cesar
 
    TTree *my_tree;
    std::vector<std::string> *vec_detId_plus_APV;
    std::vector<double> *vec_APVCMN;
    std::vector<int> *vec_fedID;
    std::vector<int> *vec_lumiS;

};

SiStripCMNanalyzer::SiStripCMNanalyzer(const edm::ParameterSet& conf){

  cablingToken_ = esConsumes(); //cesar

  usesResource(TFileService::kSharedResource);

  srcAPVCM_ = consumes<edm::DetSetVector<SiStripRawDigi>>(conf.getParameter<edm::InputTag>("srcAPVCM"));

  sdMisc_= fs_->mkdir("Miscellanea");

  my_tree = sdMisc_.make<TTree>("track_tree","track_tree");

  vec_detId_plus_APV = new std::vector<std::string>;
  vec_APVCMN = new std::vector<double>;
  vec_fedID = new std::vector<int>;
  vec_lumiS = new std::vector<int>;

  my_tree->Branch("vec_detId_plus_APV", "std::vector<std::string>", &    vec_detId_plus_APV );
  my_tree->Branch("vec_APVCMN", "std::vector<double>", &    vec_APVCMN );
  my_tree->Branch("vec_fedID", "std::vector<int>", &  vec_fedID );
  my_tree->Branch("vec_lumiS", "std::vector<int>", &  vec_lumiS );

}


SiStripCMNanalyzer::~SiStripCMNanalyzer()
{

}

void
SiStripCMNanalyzer::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {

  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcAPVCM",edm::InputTag("siStripZeroSuppression","APVCMVirginRaw"));
  descriptions.add("siStripCMNanalyzer", desc);

}

void
SiStripCMNanalyzer::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  using namespace edm;

  vec_detId_plus_APV->clear();
  vec_APVCMN->clear();
  vec_fedID->clear();
  vec_lumiS->clear();

  //edm::ESHandle<SiStripFedCabling> c;
  //es.get<SiStripFedCablingRcd>().get( c );
  //cabling_ = c.product();

  // read the cabling map from the EventSetup
  const auto& cabling_ = &es.getData(cablingToken_);

  unsigned int lumiSection = e.luminosityBlock();
  //std::cout<<"lumiSection  :  "<<lumiSection<<std::endl;


  edm::Handle<edm::DetSetVector<SiStripRawDigi> > moduleCM; 
  e.getByToken(srcAPVCM_,moduleCM); 

  for(const auto & set : *moduleCM){
    ///std::cout<<" det id  :  "<< set.id <<std::endl;
    uint32_t detID = set.id;
    uint16_t countAPVs=0; 
    for(const auto & itCM : set){
      ///std::cout<<"countAPVs "<<countAPVs<<std::endl;
      std::string aux_detID_plus_apv_st="DetId"+std::to_string(detID)+"_APV"+std::to_string(countAPVs); 
      vec_detId_plus_APV->push_back(aux_detID_plus_apv_st);
      vec_lumiS->push_back(lumiSection);
      uint16_t adc = itCM.adc();
      ///std::cout<<"adc  =  "<<adc<<std::endl;
      vec_APVCMN->push_back(2*adc-1024);
      countAPVs++;
      std::vector<uint16_t>::const_iterator itfed = cabling_->fedIds().begin();
      bool foundFEDId=false;   
      for ( ; itfed != cabling_->fedIds().end(); itfed++ ) {
        // get the cabling connections for this FED
        auto conns1 = cabling_->fedConnections(*itfed); 
        std::vector<FedChannelConnection>::const_iterator itconn = conns1.begin();
        for ( ; itconn != conns1.end(); itconn++ ) {
          if ( !itconn->detId() || itconn->detId() == sistrip::invalid32_ ) continue;
          if (itconn->detId() == detID){
            int fedID = (*itfed);
            vec_fedID->push_back(fedID); 
            foundFEDId=true; 
            break;
          }
        }
        if(foundFEDId)break; 
      }      
    }

  }

  my_tree->Fill(); 

}

void
SiStripCMNanalyzer::beginRun(edm::Run const& iEvent, edm::EventSetup const&)
{

}

void
SiStripCMNanalyzer::endRun(edm::Run const& iEvent, edm::EventSetup const&)
{

}


// ------------ method called once each job just before starting event loop  ------------
void SiStripCMNanalyzer::beginJob()
{ 

}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripCMNanalyzer::endJob() 
{
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripCMNanalyzer);

