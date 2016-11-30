#include <memory>
#include <iostream>
#include "TTree.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"

#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"

using namespace std;

class NtupleProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit NtupleProducer(const edm::ParameterSet&);
  ~NtupleProducer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetToken HGCalTBRecHitCollection_;
  edm::EDGetToken HGCalTBClusterCollection_;
  int CERN_8layers_config;
  double minEnergy;
  int sensorSize;
  double mipToMeV; //obatained from simulation

  std::vector<float> layerZPosition;
  std::vector<double> skirocADCToMip;

  TTree* tree;
  int _evtID;
  int _nhit;
  std::vector<int> _cellID;
  std::vector<double> _x;
  std::vector<double> _y;
  std::vector<double> _z;
  std::vector<double> _energy;
  
  float _thrustX0;
  float _thrustX;
  float _thrustY0;
  float _thrustY;
  
  int _nclusters;
  std::vector<double> _cluster_x;
  std::vector<double> _cluster_y;
  std::vector<double> _cluster_z;
  std::vector<double> _cluster_energy;
  std::vector<int> _cluster_size;
  
  std::string mapfile_ = "HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt";
  struct {
    HGCalElectronicsMap emap_;
  } essource_;  
};

NtupleProducer::NtupleProducer(const edm::ParameterSet& iConfig) :
  CERN_8layers_config( iConfig.getUntrackedParameter<int>("CERN_8layers_config",0) ),
  minEnergy( iConfig.getUntrackedParameter<double>("minEnergy",30) ),
  sensorSize( iConfig.getUntrackedParameter<int>("sensorSize",128) ),
  mipToMeV( iConfig.getUntrackedParameter<double>("MIPToMeV",52.81e-03) )
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
  HGCalTBClusterCollection_ = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBCLUSTERS"));

  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(mapfile_);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  }
  std::set<int> layers=essource_.emap_.layersInMap();
  std::cout << "list of layers : \t";
  for( std::set<int>::const_iterator it=layers.begin(); it!=layers.end(); ++it )
    std::cout << (*it) << ",";
  std::cout << std::endl;

  std::set<int> skirocs=essource_.emap_.skirocsInMap();
  std::cout << "list of skirocs : \t";
  for( std::set<int>::const_iterator it=skirocs.begin(); it!=skirocs.end(); ++it )
    std::cout << (*it) << ",";
  std::cout << std::endl;
  
  if( CERN_8layers_config==0 ){
    float sum=0.;
    sum+=0.0 ; layerZPosition.push_back( sum );
    sum+=5.35; layerZPosition.push_back( sum );
    sum+=5.17; layerZPosition.push_back( sum );
    sum+=3.92; layerZPosition.push_back( sum );
    sum+=4.08; layerZPosition.push_back( sum );
    sum+=1.15; layerZPosition.push_back( sum );
    sum+=4.11; layerZPosition.push_back( sum );
    sum+=2.14; layerZPosition.push_back( sum );
  }
  else if( CERN_8layers_config==1 ){
    float sum=0.;
    sum+=0.0 ; layerZPosition.push_back( sum );
    sum+=4.67; layerZPosition.push_back( sum );
    sum+=5.17; layerZPosition.push_back( sum );
    sum+=4.43; layerZPosition.push_back( sum );
    sum+=4.98; layerZPosition.push_back( sum );
    sum+=1.15; layerZPosition.push_back( sum );
    sum+=5.40; layerZPosition.push_back( sum );
    sum+=5.60; layerZPosition.push_back( sum );
  }

  double adctomip[]={16.95,16.6933,16.0208,17.0226,17.6833,17.1882,16.4708,15.9629,17.1542,16.7324,16.5,17.5457,15.3652,16.273,16.4111,15.2706};
  std::vector<double> vec; vec.insert( vec.begin(), adctomip, adctomip+16 );
  skirocADCToMip = iConfig.getUntrackedParameter< std::vector<double> >("skirocADCToMip",vec);

  if( skirocADCToMip.size() != skirocs.size() || layerZPosition.size()!=layers.size() ){
    std::cout << "problem in parameter initialisation : \n"
	      << "nlayers.size() = " << layers.size() << "=?="
	      << "layerZPosition.size() = " << layerZPosition.size() << "\n"
	      << "skirocADCToMip.size() = " << skirocADCToMip.size() << "=?="
	      << "skirocs.size() = " << skirocs.size() << "=?="
	      << "=======> throw" << std::endl;
    throw;
  }
  
  _evtID = 0;
  tree = fs->make<TTree>("HGC_Events", "HGCAL TB variables tree");
  tree->Branch("evtID",&_evtID);
  tree->Branch("nhit",&_nhit);
  tree->Branch("cellID","std::vector<int>",&_cellID);
  tree->Branch("x","std::vector<double>", &_x);
  tree->Branch("y","std::vector<double>", &_y);
  tree->Branch("z","std::vector<double>", &_z);
  tree->Branch("energy","std::vector<double>", &_energy);
  tree->Branch("thrustX0",&_thrustX0);
  tree->Branch("thrustX",&_thrustX);
  tree->Branch("thrustY0",&_thrustY0);
  tree->Branch("thrustY",&_thrustY);
  tree->Branch("nclusters",_nclusters);
  tree->Branch("cluster_x","std::vector<double>",&_cluster_x);
  tree->Branch("cluster_y","std::vector<double>",&_cluster_y);
  tree->Branch("cluster_z","std::vector<double>",&_cluster_z);
  tree->Branch("cluster_energy","std::vector<double>",&_cluster_energy);
  tree->Branch("cluster_size","std::vector<int>",&_cluster_size);
}

NtupleProducer::~NtupleProducer()
{
}

void
NtupleProducer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  _evtID++;

  _nhit=0;
  _nclusters=0;
  _thrustX0=0;
  _thrustX=0;
  _thrustY0=0;
  _thrustY=0;

  _cellID.clear();
  _x.clear();
  _y.clear();
  _z.clear();
  _energy.clear();

  _cluster_x.clear();
  _cluster_y.clear();
  _cluster_z.clear();
  _cluster_energy.clear();
  _cluster_size.clear();

  edm::Handle<HGCalTBRecHitCollection> Rechits;
  event.getByToken(HGCalTBRecHitCollection_, Rechits);

  edm::Handle<reco::HGCalTBClusterCollection> clusters;
  event.getByToken(HGCalTBClusterCollection_, clusters);
  
  HGCalTBCellVertices cellVertice;
  std::pair<double, double> CellCentreXY;
  for( auto hit : *Rechits ){
    if( hit.energy()<minEnergy ) continue;
    _nhit++;
    _cellID.push_back( hit.id() );
    CellCentreXY=cellVertice.GetCellCentreCoordinatesForPlots( hit.id().layer(), hit.id().sensorIU(), hit.id().sensorIV(), hit.id().iu(), hit.id().iv(), sensorSize);
    _x.push_back( CellCentreXY.first );
    _y.push_back( CellCentreXY.second );
    _z.push_back( layerZPosition.at( hit.id().layer()-1 ) );
    HGCalTBElectronicsId eid= essource_.emap_.detId2eid(hit.id());
    _energy.push_back( hit.energy()/skirocADCToMip[ eid.iskiroc()-1 ]*mipToMeV );
  }

  _nclusters=clusters->size();
  for( auto cluster : *clusters ){
    _cluster_x.push_back( cluster.x() );
    _cluster_y.push_back( cluster.y() );
    _cluster_z.push_back( cluster.z() );
    float en=0.;
    for( std::vector< std::pair<DetId, float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it ){
      HGCalTBElectronicsId eid= essource_.emap_.detId2eid( (*it).first );
      en += (*it).second/skirocADCToMip[ eid.iskiroc()-1 ]*mipToMeV;
    }
    _cluster_energy.push_back( en );
    _cluster_size.push_back( cluster.size() );
  }
  tree->Fill();
}

void
NtupleProducer::beginJob()
{
}

void
NtupleProducer::endJob()
{

}

void
NtupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NtupleProducer);
