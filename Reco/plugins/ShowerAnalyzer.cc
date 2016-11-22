#include <memory>
#include <iostream>
#include "TTree.h"
#include <sstream>
#include <cmath>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBCaloTrack.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"
#include "HGCal/Geometry/interface/HGCalTBSpillParameters.h"

#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"

#include "HGCal/Reco/interface/HGCalTBSortingHelper.h"
#include "HGCal/Reco/interface/HGCalTBCaloTrackingUtil.h"
#include "HGCal/Reco/interface/Distance.h"

using namespace std;

class ShowerAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit ShowerAnalyzer(const edm::ParameterSet&);
  ~ShowerAnalyzer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetToken HGCalTBRecHitCollection_;
  edm::EDGetToken HGCalTBClusterCollection_;
  edm::EDGetToken HGCalTBClusterCollection7_;
  edm::EDGetToken HGCalTBClusterCollection19_;
  int nlayers;
  int nskirocsperlayer;
  int sensorSize;
  int CERN_8layers_config;
  double energyShowerThreshold ;
  float maxTransverseProfile;

  std::vector<float> layerZPosition;
  std::vector<float> layerZX0;
  std::vector<double> skirocADCToMip;

  TTree* tree;
  int _evtID;
  float _theta;
  float _phi;
  int _nhit;
  float _x0;
  float _y0;
  float _ax;
  float _ay;
  float _energyInCluster;
  std::vector<int> _clustersizelayer;
  std::vector<double> _transverseprofile;
  std::vector<double> _energylayer;
  std::vector<double> _clusterenergylayer;
  std::vector<double> _meanx;
  std::vector<double> _meany;
  std::vector<double> _rmsx;
  std::vector<double> _rmsy;
  
    std::string mapfile_ = "HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt";
  struct {
    HGCalElectronicsMap emap_;
  } essource_;  
  const float PI = 3.1415927;
  SortByEnergy<reco::HGCalTBCluster,reco::HGCalTBCluster> energySorter;
};

ShowerAnalyzer::ShowerAnalyzer(const edm::ParameterSet& iConfig) :
  nlayers( iConfig.getUntrackedParameter<int>("NLayers",8) ),
  nskirocsperlayer( iConfig.getUntrackedParameter<int>("NSkirocsPerLayer",2) ),
  sensorSize( iConfig.getUntrackedParameter<int>("SensorSize",128) ),
  CERN_8layers_config( iConfig.getUntrackedParameter<int>("CERN_8layers_config",0) ),
  energyShowerThreshold( iConfig.getUntrackedParameter<double>("energyShowerThreshold",200) ),
  maxTransverseProfile( iConfig.getUntrackedParameter<double>("maxTransverseProfile",20) )
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
  HGCalTBClusterCollection_ = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBCLUSTERS"));
  // HGCalTBClusterCollection7_ = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBCLUSTERS7"));
  // HGCalTBClusterCollection19_ = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBCLUSTERS19"));

  _evtID = 0;
 
  if( CERN_8layers_config==0 ){
    float x0sum=0;
    x0sum+=6.268; layerZX0.push_back( x0sum );
    x0sum+=1.131; layerZX0.push_back( x0sum );
    x0sum+=1.131; layerZX0.push_back( x0sum );
    x0sum+=1.362; layerZX0.push_back( x0sum );
    x0sum+=0.574; layerZX0.push_back( x0sum );
    x0sum+=1.301; layerZX0.push_back( x0sum );
    x0sum+=0.574; layerZX0.push_back( x0sum );
    x0sum+=2.420; layerZX0.push_back( x0sum );   
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
    float x0sum=0;
    x0sum+=5.048; layerZX0.push_back( x0sum ); 
    x0sum+=3.412; layerZX0.push_back( x0sum ); 
    x0sum+=3.412; layerZX0.push_back( x0sum ); 
    x0sum+=2.866; layerZX0.push_back( x0sum ); 
    x0sum+=2.512; layerZX0.push_back( x0sum ); 
    x0sum+=1.625; layerZX0.push_back( x0sum ); 
    x0sum+=2.368; layerZX0.push_back( x0sum ); 
    x0sum+=6.021; layerZX0.push_back( x0sum );
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

  double adctomip[]={16.9426, 16.6226, 15.8083, 16.9452, 17.95, 16.6588, 16.6542, 15.4429, 17.3042, 16.5757, 16.35, 17.4657, 15.3609, 16.2676, 16.5, 15.1118};
  std::vector<double> vec; vec.insert( vec.begin(), adctomip, adctomip+16 );
  skirocADCToMip = iConfig.getUntrackedParameter< std::vector<double> >("skirocADCToMip",vec);

  if( skirocADCToMip.size() != (unsigned int)nlayers*nskirocsperlayer ){
    std::cout << "problem in parameter initialisation : \n"
	      << "nlayers = " << nlayers << "\n"
	      << "nskirocsperlayer = " << nskirocsperlayer << "\n"
	      << "skirocADCToMip.size() = " << skirocADCToMip.size() << " while it should be equal to " << nlayers*nskirocsperlayer << " (nlayers*nskirocsperlayer) \n"
	      << "=======> throw" << std::endl;
    throw;
  }
  
  tree = fs->make<TTree>("tree", "HGCAL TB variables tree");
  tree->Branch( "evtID",&_evtID ); 
  tree->Branch( "theta",&_theta );
  tree->Branch( "phi",&_phi );
  tree->Branch( "nhit",&_nhit );
  tree->Branch( "vx0",&_x0 );
  tree->Branch( "vy0",&_y0 );
  tree->Branch( "ax",&_ax );
  tree->Branch( "ay",&_ay );
  tree->Branch( "energyInCluster",&_energyInCluster );

  tree->Branch( "clustersizelayer","std::vector<int>",&_clustersizelayer);
  tree->Branch( "transverseprofile","std::vector<double>",&_transverseprofile);
  tree->Branch( "energylayer","std::vector<double>",&_energylayer);
  tree->Branch( "clusterenergylayer","std::vector<double>",&_clusterenergylayer);
  tree->Branch( "meanx","std::vector<double>",&_meanx);
  tree->Branch( "meany","std::vector<double>",&_meany);
  tree->Branch( "rmsx","std::vector<double>",&_rmsx);
  tree->Branch( "rmsy","std::vector<double>",&_rmsy);
}

ShowerAnalyzer::~ShowerAnalyzer()
{
}

void
ShowerAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  _evtID++;
  _transverseprofile.clear();
  _energylayer.clear();
  _clusterenergylayer.clear();
  _clustersizelayer.clear();
  _meanx.clear();
  _meany.clear();
  _rmsx.clear();
  _rmsy.clear();
  _energyInCluster=0;
  edm::Handle<HGCalTBRecHitCollection> Rechits;
  event.getByToken(HGCalTBRecHitCollection_, Rechits);

  edm::Handle<reco::HGCalTBClusterCollection> clusters;
  event.getByToken(HGCalTBClusterCollection_, clusters);
  
  // edm::Handle<reco::HGCalTBClusterCollection> clusters7;
  // event.getByToken(HGCalTBClusterCollection7_, clusters7);
  // edm::Handle<reco::HGCalTBClusterCollection> clusters19;
  // event.getByToken(HGCalTBClusterCollection19_, clusters19);
        
  for( int ir=0; ir<maxTransverseProfile; ir++ ){
    _transverseprofile.push_back(0);
  }

  //
  
  reco::HGCalTBClusterCollection shower;
  _nhit=0;
  for( int ilayer=0; ilayer<nlayers; ilayer++){
    _energylayer.push_back(0.0);
    _meanx.push_back(0.0);
    _meany.push_back(0.0);
    _rmsx.push_back(0.0);
    _rmsy.push_back(0.0);
    _clusterenergylayer.push_back(0.0);
    _clustersizelayer.push_back(0);
    std::vector<reco::HGCalTBCluster> temp;
    for( auto cluster : *clusters ){
      if( cluster.layer()!=ilayer ) continue;
      temp.push_back(cluster);
      _nhit+=cluster.size();
    }
    if( temp.empty() ) continue;
    std::sort( temp.begin(), temp.end(), energySorter.sort );
    shower.push_back( (*temp.begin()) );
    _clusterenergylayer.at(ilayer) = (*temp.begin()).energy();
    _clustersizelayer.at(ilayer) = (*temp.begin()).size();
    _energyInCluster+=(*temp.begin()).energy();
  }

  Distance<HGCalTBRecHit,reco::HGCalTBCaloTrack> dist;
  if( shower.size()>=4 ){
    std::vector<HGCalTBDetId> detIds;
    reco::WeightedLeastSquare<reco::HGCalTBClusterCollection> wls;
    std::vector<float> trackPar;
    std::vector<float> trackParError;
    wls.run( shower, trackPar, trackParError);
    float z0=layerZPosition.at(0);
    Vector momentum = Vector(-1., 0., trackPar[1]).Cross( Vector(0., -1., trackPar[3]) );
    Point vertex = Point( trackPar[0]+trackPar[1]*z0, trackPar[2]+trackPar[3]*z0, z0 );
    reco::HGCalTBCaloTrack track( wls.chi2( shower, trackPar), shower.size(), vertex, momentum,/* cov,*/ detIds);
    _theta = track.momentum().theta()*180/PI;
    _phi = track.momentum().phi()*180/PI;
    _x0 = track.vertex().x();
    _y0 = track.vertex().y();
    _ax = track.momentum().x();
    _ay = track.momentum().y();    
  }
  
  //  HGCalTBCellVertive cellVertice;
  for( auto hit : *Rechits ){
    _energylayer[ hit.id().layer()-1 ] += hit.energy();
    //_meanx[ (*it).id().layer()-1 ]+=;
    //_meany[ (*it).id().layer()-1 ]+=;
    //_rmsx[ (*it).id().layer()-1 ];
    //_rmsy[ (*it).id().layer()-1 ];
  }
  tree->Fill();
}

void
ShowerAnalyzer::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(mapfile_);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  }
}

void
ShowerAnalyzer::endJob()
{

}

void
ShowerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ShowerAnalyzer);
