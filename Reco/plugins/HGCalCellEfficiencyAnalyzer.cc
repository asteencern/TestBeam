#include <memory>
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <sstream>
#include <cmath>
#include <limits>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/Point3D.h"

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/DataFormats/interface/HGCalTBCaloTrack.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"

#include "HGCal/Reco/interface/HGCalTBCaloTrackingUtil.h"

#include "HGCal/Geometry/interface/HGCalTBDetIdFinder.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBSpillParameters.h"

using namespace std;

class HGCalCellEfficiencyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit HGCalCellEfficiencyAnalyzer(const edm::ParameterSet&);
  ~HGCalCellEfficiencyAnalyzer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetToken HGCalTBRecHitCollection_;
  edm::EDGetToken HGCalTBClusterCollection_;
  
  int _nlayers;
  int _nskirocsperlayer;
  int _nchannelsperskiroc;
  int _sensorSize;
  int _minTouchedLayers;
  double _maxChi2;
  double _minEnergyEfficiency;
  double _clusterMaxEnergyThresholdForTrack;
  double _maxDistanceToRecoTrack;
  std::vector<double> _layerZPositions;

  TTree* tree;
  std::vector<double> _chi2;
  std::vector<double> _theta;
  std::vector<double> _phi;
  std::vector<int> _tracknhit;
  std::vector<double> _vx0;
  std::vector<double> _vy0;
  std::vector<double> _px;
  std::vector<double> _py;
  std::vector<int> _efficiencyLayer;
  int _efficiency;
  int _runNumber;
  int _evtID;

  std::map<int, TH1F*> h_energy_map;
  std::map<int, TH1F*> h_distance_map;
  TH2F* h_efficiency_distance;

  const float PI = 3.1415927;

  std::string mapfile_ = "HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt";
  struct {
    HGCalElectronicsMap emap_;
  } essource_;

  std::string ghostMapFile_ = "HGCal/CondObjects/data/map_emptyDetId.txt";
  struct {
    HGCalElectronicsMap emap_;
  } essourceGhost_;

};


HGCalCellEfficiencyAnalyzer::HGCalCellEfficiencyAnalyzer(const edm::ParameterSet& iConfig) :
  _nlayers( iConfig.getUntrackedParameter<int>("NLayers",8) ),
  _nskirocsperlayer( iConfig.getUntrackedParameter<int>("NSkirocsPerLayer",2) ),
  _nchannelsperskiroc( iConfig.getUntrackedParameter<int>("NChannelsPerSkiroc",64) ),
  _sensorSize( iConfig.getUntrackedParameter<int>("SensorSize",128) ),
  _minTouchedLayers( iConfig.getUntrackedParameter<int>("minTouchedLayers",4) ),
  _maxChi2( iConfig.getUntrackedParameter<double>("maxChi2",9.48) ),
  _minEnergyEfficiency( iConfig.getUntrackedParameter<double>("minEnergyEfficiency",3.1e-5) ),
  _clusterMaxEnergyThresholdForTrack( iConfig.getUntrackedParameter<double>("clusterMaxEnergyThreshold",2.1e-4) ),//GeV, 4 MIP
  _maxDistanceToRecoTrack( iConfig.getUntrackedParameter<double>("maxDistanceToRecoTrack",2.0) )
{
  HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
  HGCalTBClusterCollection_ = consumes<reco::HGCalTBClusterCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBCLUSTERS"));
  float sum=0.;
  std::vector<double> vec;
  sum+=0.0 ; vec.push_back( sum );
  sum+=5.35; vec.push_back( sum );
  sum+=5.17; vec.push_back( sum );
  sum+=3.92; vec.push_back( sum );
  sum+=4.08; vec.push_back( sum );
  sum+=1.15; vec.push_back( sum );
  sum+=4.11; vec.push_back( sum );
  sum+=2.14; vec.push_back( sum );
  //cern config 1 (5X0->15X0) is default
  _layerZPositions = iConfig.getUntrackedParameter< std::vector<double> >("LayerZPositions",vec);
  std::cout << iConfig.dump() << std::endl;

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "HGCAL TB variables tree");
  tree->Branch( "run",&_runNumber );
  tree->Branch( "event",&_evtID );
  tree->Branch( "chi2","std::vector<double>",&_chi2 ); 
  tree->Branch( "theta","std::vector<double>",&_theta );
  tree->Branch( "phi","std::vector<double>",&_phi );
  tree->Branch( "tracknhit","std::vector<int>",&_tracknhit );
  tree->Branch( "vx0","std::vector<double>",&_vx0 );
  tree->Branch( "vy0","std::vector<double>",&_vy0 );
  tree->Branch( "px","std::vector<double>",&_px );
  tree->Branch( "py","std::vector<double>",&_py );
  tree->Branch( "efficiencyLayer","std::vector<int>",&_efficiencyLayer );

  std::ostringstream os( std::ostringstream::ate );
  for( int ilayer=0; ilayer<_nlayers; ilayer++ ){
    for( int iskiroc=0; iskiroc<_nskirocsperlayer; iskiroc++ ){
      for( int ichannel=0; ichannel< _nchannelsperskiroc; ichannel++ ){
	os.str("");
	os << "Layer" << ilayer << "_Skiroc" << iskiroc << "_Channel" << ichannel << "";
	int key=ilayer*1000+iskiroc*100+ichannel;
	TH1F* h = fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 5000, -2, 8 );
	h_energy_map[key]=h;
      }
    }
    os.str("");
    os << "Layer" << ilayer << "_distance";
    h_distance_map[ilayer]=fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 600, -0.5, 5.5 );;
  }
  h_efficiency_distance = fs->make<TH2F>("EfficiencyVsDistance", "", 1000, 0, 3*sqrt(3)/2*HGCAL_TB_CELL::FULL_CELL_SIDE, 10,-0.5,1.5);
}

HGCalCellEfficiencyAnalyzer::~HGCalCellEfficiencyAnalyzer()
{
}

void
HGCalCellEfficiencyAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  _tracknhit.clear();
  _chi2.clear();
  _theta.clear();
  _phi.clear();
  _vx0.clear();
  _vy0.clear();
  _px.clear();
  _py.clear();
  _efficiencyLayer.clear();
  _runNumber=event.run();
  _evtID=(event.id().luminosityBlock()-1)*EVENTSPERSPILL + event.id().event();
  edm::Handle<HGCalTBRecHitCollection> Rechits;
  event.getByToken(HGCalTBRecHitCollection_, Rechits);
  edm::Handle<reco::HGCalTBClusterCollection> clusters;
  event.getByToken(HGCalTBClusterCollection_, clusters);

  HGCalTBDetIdFinder finder(essource_.emap_,_sensorSize);

  HGCalTBCellVertices cellVertice;
  for( int ilayer=0; ilayer<_nlayers; ilayer++ ){
    _efficiencyLayer.push_back(-1);
    _tracknhit.push_back(-1);;
    _chi2.push_back(-1);
    _theta.push_back(-100);
    _phi.push_back(-100);
    _vx0.push_back(-100);
    _vy0.push_back(-100);
    _px.push_back(-100);
    _py.push_back(-100);
    _efficiency=0;
    reco::HGCalTBClusterCollection tmp;
    std::set<int> touchedLayers;
    std::vector<HGCalTBDetId> detIds;
    for( auto cluster : *clusters ){
      if( cluster.energy()>_clusterMaxEnergyThresholdForTrack || cluster.layer()==ilayer+1 || cluster.size()>2 )
	continue;
      tmp.push_back(cluster);
      touchedLayers.insert( cluster.layer() );
      for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it )
	detIds.push_back( HGCalTBDetId((*it).first) );
    }
    if( (int)touchedLayers.size()<_minTouchedLayers ) continue;
    reco::LeastSquare<reco::HGCalTBClusterCollection> ls;
    std::vector<float> trackPar;
    std::vector<float> trackParError;
    ls.run( tmp, trackPar, trackParError);
    float chi2 = ls.chi2( tmp, trackPar);
    int ndof = tmp.size();
    Vector momentum = Vector(-1., 0., trackPar[1]).Cross( Vector(0., -1., trackPar[3]) );
    Point vertex = Point( trackPar[0], trackPar[2], 0.0 );
    reco::HGCalTBCaloTrack track( chi2, ndof, vertex, momentum,/* cov,*/ detIds);
    reco::HGCalTBClusterCollection cleancol;
    reco::TrackCleaner cleaner;
    cleaner.clean( tmp, cleancol, track, _maxDistanceToRecoTrack );
    touchedLayers.clear();
    detIds.clear();
    for( auto cluster : cleancol ){
      touchedLayers.insert( cluster.layer() );
      for( std::vector< std::pair<DetId,float> >::const_iterator it=cluster.hitsAndFractions().begin(); it!=cluster.hitsAndFractions().end(); ++it )
	detIds.push_back( HGCalTBDetId((*it).first) );
    }
    if( (int)touchedLayers.size()<_minTouchedLayers ) continue;
    ls.run( cleancol, trackPar, trackParError);
    chi2 = ls.chi2( cleancol, trackPar);
    ndof = cleancol.size();
    momentum = Vector(-1., 0., trackPar[1]).Cross( Vector(0., -1., trackPar[3]) );
    vertex = Point( trackPar[0], trackPar[2], 0.0 );
    track = reco::HGCalTBCaloTrack( chi2, ndof, vertex, momentum,/* cov,*/ detIds);
    
    if( track.isNull() || track.normalisedChi2()>_maxChi2 )
      continue;

    math::XYZPoint p=track.expectedTrackProjection( _layerZPositions.at(ilayer) );
    std::pair<double,double> xy(p.x(),p.y());
    HGCalTBDetId id;
    HGCalTBRecHit hit;
    finder.run(id,xy,ilayer+1);
    if( essourceGhost_.emap_.existsDetId( id ) ){
      int iu,iv;iu=iv=0;
      //non flipped sensor:
      if( id.iu()==4 && id.iv()==-6 ){ iu=4; iv=-7; }
      if( id.iu()==7 && id.iv()==-4 ){ iu=7; iv=-3; }
      if( id.iu()==3 && id.iv()==4 ){ iu=4; iv=3; }
      if( id.iu()==-3 && id.iv()==7 ){ iu=-2; iv=6; }
      if( id.iu()==-4 && id.iv()==7 ){ iu=-4; iv=6; }
      if( id.iu()==-7 && id.iv()==3 ){ iu=-7; iv=4; }
      if( id.iu()==-4 && id.iv()==-3 ){ iu=-3; iv=-4; }
      if( id.iu()==2 && id.iv()==-6 ){ iu=3; iv=-7; }
      //flipped sensor:
      if( id.iu()==-2 && id.iv()==6 ){ iu=-3; iv=7; }
      if( id.iu()==4 && id.iv()==-7 ){ iu=4; iv=-6; }
      if( id.iu()==3 && id.iv()==-7 ){ iu=2; iv=-6; }
      if( id.iu()==-4 && id.iv()==6 ){ iu=-4; iv=7; }
      //6th layer:
      if( id.iu()==-3 && id.iv()==-4 ){ iu=-4; iv=-3; }
      HGCalTBDetId newid( id.layer(),id.sensorIU(),id.sensorIV(),iu,iv,id.cellType() );
      if( (*Rechits).find(newid)!=(*Rechits).end() )
	hit=(*(*Rechits).find(newid));
      else{
	//std::cout << "big problem find here : \t id = " << id << "\t newid = " << newid << std::endl;
     	//getchar();
     	continue;
      }
    }
    else if( (*Rechits).find(id)!=(*Rechits).end() )
      hit=(*(*Rechits).find(id));
    else{
      //std::cout << "even bigger problem at id = " << id << "\t" << id.cellType() << std::endl;
      //getchar();
      continue;
    }
    HGCalTBTopology top;
    int maxDist=1;
    std::set<HGCalTBDetId> detids=top.getNeighboringCellsDetID(id, _sensorSize, maxDist, essource_.emap_);
    for( std::set<HGCalTBDetId>::iterator it=detids.begin(); it!=detids.end(); ++it ){
      if( (*Rechits).find(*it)==(*Rechits).end() ) continue;
      HGCalTBRecHit ahit=(*(*Rechits).find(*it));
      if( ahit.energy()>hit.energy() ){
	hit=ahit;
	id=(*it);
      }
    }
    if( hit.energy()>_minEnergyEfficiency )
      _efficiency=1;
    std::pair<double,double> xycentre=cellVertice.GetCellCentreCoordinatesForPlots( ilayer, id.sensorIU(), id.sensorIV(), id.iu(), id.iv(), _sensorSize);
    math::XYZPoint xyzcentre(xycentre.first,xycentre.second,_layerZPositions.at(ilayer));
    
    h_efficiency_distance->Fill( std::sqrt((xyzcentre-p).mag2()),_efficiency );
    if( _efficiency==1 )
      h_distance_map[ilayer]->Fill( std::sqrt((xyzcentre-p).mag2()) );

    uint32_t EID = essource_.emap_.detId2eid( hit.id() );
    HGCalTBElectronicsId eid(EID);
    int key=ilayer*1000 + (eid.iskiroc() - 1)%2*100 + eid.ichan();
    h_energy_map[key]->Fill( hit.energy()/51.9e-6 );
    _tracknhit[ilayer] = track.getDetIds().size();
    _theta[ilayer] = track.momentum().theta()*180/PI;
    _phi[ilayer] = track.momentum().phi()*180/PI;
    _vx0[ilayer] = track.vertex().x();
    _vy0[ilayer] = track.vertex().y();
    _px[ilayer] = track.momentum().x();
    _py[ilayer] = track.momentum().y();
    _chi2[ilayer]=track.normalisedChi2();
    _efficiencyLayer[ilayer] = _efficiency;
    
    //if( _efficiency==0 ){
    //  std::cout << event.id().luminosityBlock() << "\t" << event.id().event() << "\t"<< "\t" << id << "\t" << hit << "\t" << _chi2 << "\t" << _efficiency << "\t" << _tracknhit << "\t (" << p.x() << ", " << p.y() << ", " << p.z() << ")\t" << (int)touchedLayers.size() << std::endl;
    //  //for( std::vector<HGCalTBDetId>::iterator it=track.getDetIds().begin(); it!=track.getDetIds().end(); ++it ){
    //  //	//xy=cellVertice.GetCellCentreCoordinatesForPlots( (*it).layer(), (*it).sensorIU(), (*it).sensorIV(), (*it).iu(), (*it).iv(), _sensorSize);
    //  //	//xyz = math::XYZPoint(xycentre.first,xycentre.second,_layerZPositions.at( (*it).layer()-1));
    //  //	std::cout << (*it) << std::endl;
    //  //}
    //}
  }
  tree->Fill();
}

void
HGCalCellEfficiencyAnalyzer::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(mapfile_);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  } 

  fip = edm::FileInPath(ghostMapFile_);
  if (!io.load(fip.fullPath(), essourceGhost_.emap_)) {
    throw cms::Exception("Unable to load electronics map of ghost channels");
  } 
}

void
HGCalCellEfficiencyAnalyzer::endJob()
{

}

void
HGCalCellEfficiencyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalCellEfficiencyAnalyzer);
