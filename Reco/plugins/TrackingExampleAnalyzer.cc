#include <memory>
#include <iostream>
#include "TH1F.h"
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

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/DataFormats/interface/HGCalTBCaloTrack.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"

#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"

#include "HGCal/Reco/interface/HGCalTBCaloTrackingUtil.h"

using namespace std;

class TrackingExampleAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit TrackingExampleAnalyzer(const edm::ParameterSet&);
  ~TrackingExampleAnalyzer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetToken HGCalTBRecHitCollection_;
  int nlayers;
  int nskirocsperlayer;
  int sensorSize;
  int minMip;
  int maxMip;
  int cmThreshold;
  int CERN_8layers_config;
  bool prepareTreeForDisplay;
  bool doTrackCleaning;
  double maxDistanceToRecoTrack;
  double maxChi2;

  std::vector<float> layerZPosition;

  TTree* tree;
  int _evtID;
  float _chi2;
  float _theta;
  float _phi;
  int _nhit;
  int _tracknhit;
  float _x0;
  float _y0;
  float _ax;
  float _ay;
  bool _trackSuccess;

  std::vector<double> _deltas;

  std::vector<double> _nhitlayer;
  std::vector<double> _energylayer;
  std::vector<double> _meanx;
  std::vector<double> _meany;
  std::vector<double> _rmsx;
  std::vector<double> _rmsy;
  std::vector<double> _commonMode;
  
  std::vector<double> _x;
  std::vector<double> _y;
  std::vector<double> _z;
  std::vector<double> _energy;
  
  std::map<int,TH1F*> h_delta_layer;

  std::map<int, TH1F*> h_mip_map;
  std::map<int, TH1F*> h_noise_map;

  const float PI = 3.1415927;

  std::string mapfile_ = "HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt";
  struct {
    HGCalElectronicsMap emap_;
  } essource_;

};


TrackingExampleAnalyzer::TrackingExampleAnalyzer(const edm::ParameterSet& iConfig) :
  nlayers( iConfig.getUntrackedParameter<int>("NLayers",8) ),
  nskirocsperlayer( iConfig.getUntrackedParameter<int>("NSkirocsPerLayer",2) ),
  sensorSize( iConfig.getUntrackedParameter<int>("SensorSize",128) ),
  minMip( iConfig.getUntrackedParameter<int>("minMip",10) ), //~5 sigma away from noise
  maxMip( iConfig.getUntrackedParameter<int>("maxMip",48) ), //~3 mip 
  cmThreshold( iConfig.getUntrackedParameter<int>("CMThreshold",30) ),
  CERN_8layers_config( iConfig.getUntrackedParameter<int>("CERN_8layers_config",0) ),
  prepareTreeForDisplay( iConfig.getUntrackedParameter<bool>("PrepareTreeForDisplay",false) ),
  doTrackCleaning( iConfig.getUntrackedParameter<bool>("doTrackCleaning",true) ),
  maxDistanceToRecoTrack( iConfig.getUntrackedParameter<double>("maxDistanceToRecoTrack",2.0*HGCAL_TB_CELL::FULL_CELL_SIDE) ),
  maxChi2( iConfig.getUntrackedParameter<double>("maxChi2",9.48) )
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
  
  std::cout << iConfig.dump() << std::endl;
  
  tree = fs->make<TTree>("tree", "HGCAL TB variables tree");
  tree->Branch( "evtID",&_evtID ); 
  tree->Branch( "chi2",&_chi2 ); 
  tree->Branch( "theta",&_theta );
  tree->Branch( "phi",&_phi );
  tree->Branch( "nhit",&_nhit );
  tree->Branch( "tracknhit",&_tracknhit );
  tree->Branch( "vx0",&_x0 );
  tree->Branch( "vy0",&_y0 );
  tree->Branch( "px",&_ax );
  tree->Branch( "py",&_ay );
  tree->Branch( "trackSuccess",&_trackSuccess );
  tree->Branch( "nhitlayer","std::vector<double>",&_nhitlayer);
  tree->Branch( "energylayer","std::vector<double>",&_energylayer);
  tree->Branch( "meanx","std::vector<double>",&_meanx);
  tree->Branch( "meany","std::vector<double>",&_meany);
  tree->Branch( "rmsx","std::vector<double>",&_rmsx);
  tree->Branch( "rmsy","std::vector<double>",&_rmsy);
  tree->Branch( "deltas","std::vector<double>",&_deltas);
  
  if( prepareTreeForDisplay ){
    tree->Branch( "x","std::vector<double>",&_x );
    tree->Branch( "y","std::vector<double>",&_y );
    tree->Branch( "z","std::vector<double>",&_z );
    tree->Branch( "energy","std::vector<double>",&_energy );
  }

  _evtID=0;


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

  std::ostringstream os( std::ostringstream::ate );
  for( int ilayer=0; ilayer<nlayers; ilayer++ ){
    os.str("");
    os << "Delta_Layer" << ilayer;
    TH1F* h = fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 1000, 0, 10 );
    h_delta_layer[ilayer]=h;
    if( !prepareTreeForDisplay )
      for( int iskiroc=0; iskiroc<nskirocsperlayer; iskiroc++ ){
	for( int ichannel=0; ichannel<64; ichannel++ ){
	  os.str("");
	  os << "Layer" << ilayer << "_Skiroc" << iskiroc << "_Channel" << ichannel << "_MIP";
	  int key=ilayer*1000+iskiroc*100+ichannel;
	  h = fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 8000, -20, 60 );
	  h_mip_map[key]=h;
	  os.str("");
	  os << "Layer" << ilayer << "_Skiroc" << iskiroc << "_Channel" << ichannel << "_Noise";
	  h = fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 10000, -50, 50 );
	  h_noise_map[key]=h;
	}
      }
  }
  
}

TrackingExampleAnalyzer::~TrackingExampleAnalyzer()
{
}

void
TrackingExampleAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  _nhitlayer.clear(); 
  _energylayer.clear();
  _x.clear();
  _y.clear();
  _z.clear();
  _energy.clear();
  _meanx.clear();
  _meany.clear();
  _rmsx.clear();
  _rmsy.clear();
  _deltas.clear();
  _chi2 = _theta = _phi = _nhit = _x0 = _y0 = _ax = _ay = 0;

  HGCalTBCellVertices cellVertice;
  HGCalTBRecHitCollection coltmp[nlayers];
  HGCalTBRecHitCollection col;

  edm::Handle<HGCalTBRecHitCollection> Rechits;
  event.getByToken(HGCalTBRecHitCollection_, Rechits);
  
  for( int ilayer=0; ilayer<nlayers; ilayer++){
    for( auto hit : *Rechits ){
      if( hit.id().layer()-1!=ilayer ) continue;

      if( essource_.emap_.existsDetId( hit.id() )==false )
	std::cout << "problem at key = " << hit.id() << std::endl;

      coltmp[ ilayer ].push_back( hit );
    }
    _nhitlayer.push_back(0); 
    _energylayer.push_back(0.);
    _meanx.push_back(0.);
    _meany.push_back(0.);
    _rmsx.push_back(0.);
    _rmsy.push_back(0.);   
    for( std::vector<HGCalTBRecHit>::iterator it=coltmp[ ilayer ].begin(); it!=coltmp[ ilayer ].end(); ++it ){
      if( (*it).id().cellType()==1 ||
	  (*it).id().cellType()==2 || 
	  (*it).id().cellType()==3 || 
	  //(*it).id().cellType()==4 || 
	  (*it).id().cellType()==5 )
	continue;

      if( (*it).energy() > minMip && (*it).energy() < maxMip ){
	std::pair<double, double> CellCentreXY;
	CellCentreXY=cellVertice.GetCellCentreCoordinatesForPlots( (*it).id().layer(), 
								   (*it).id().sensorIU(), 
								   (*it).id().sensorIV(), 
								   (*it).id().iu(), 
								   (*it).id().iv(), 
								   sensorSize);
	(*it).setPosition( math::XYZPoint( CellCentreXY.first, 
					   CellCentreXY.second, 
					   layerZPosition[ ilayer ]) 
			   );
	col.push_back( (*it) );
	//_nhitlayer[ ilayer ]++;
	_meanx[ ilayer ] += (*it).x()*(*it).energy() ;
	_meany[ ilayer ] += (*it).y()*(*it).energy() ;
	_rmsx[ ilayer ] += (*it).x()*(*it).x()*(*it).energy();
	_rmsy[ ilayer ] += (*it).y()*(*it).y()*(*it).energy();
      }
      else if( (*it).energy() < minMip && !prepareTreeForDisplay ){
	uint32_t EID = essource_.emap_.detId2eid( (*it).id() );
	HGCalTBElectronicsId eid(EID);
	int key=( (*it).id().layer()-1 )*1000 + (eid.iskiroc() - 1)%2*100 + eid.ichan();
	h_noise_map[key]->Fill( (*it).energy() );
      }
    }
    if( _energylayer[ ilayer ]>0 ){
      _meanx[ ilayer ] /= _energylayer[ ilayer ];
      _meany[ ilayer ] /= _energylayer[ ilayer ];
      if( _rmsx[ ilayer ]/_energylayer[ ilayer ] - _meanx[ ilayer ]*_meanx[ ilayer ] < std::numeric_limits<double>::epsilon() )
	_rmsx[ ilayer ] = 0.0;
      else 
	_rmsx[ ilayer ] = std::sqrt( _rmsx[ ilayer ]/_energylayer[ ilayer ] - _meanx[ ilayer ]*_meanx[ ilayer ] );

      if( _rmsy[ ilayer ]/_energylayer[ ilayer ] - _meany[ ilayer ]*_meany[ ilayer ] < std::numeric_limits<double>::epsilon() )
	_rmsy[ ilayer ] = 0.0;
      else _rmsy[ ilayer ] = std::sqrt( _rmsy[ ilayer ]/_energylayer[ ilayer ] - _meany[ ilayer ]*_meany[ ilayer ] );
    }
  }
  
  _nhit=col.size();
  _trackSuccess=false;
  if( col.size()>=4 ){
    std::vector<HGCalTBDetId> detIds;
    reco::WeightedLeastSquare<HGCalTBRecHitCollection> wls;
    std::vector<float> trackPar;
    std::vector<float> trackParError;
    wls.run( col, trackPar, trackParError);
    float chi2 = wls.chi2( col, trackPar);
    int ndof = col.size();
    Vector momentum = Vector(-1., 0., trackPar[1]).Cross( Vector(0., -1., trackPar[3]) );
    float z0=layerZPosition.at(0);
    Point vertex = Point( trackPar[0]+trackPar[1]*z0, trackPar[2]+trackPar[3]*z0, z0 );
    reco::HGCalTBCaloTrack track( chi2, ndof, vertex, momentum,/* cov,*/ detIds);

    if( doTrackCleaning==true ){
      reco::TrackCleaner cleaner;
      HGCalTBRecHitCollection cleancol;
      cleaner.clean( col, cleancol, track, maxDistanceToRecoTrack );
      wls.run( cleancol, trackPar, trackParError);
      for( std::vector<HGCalTBRecHit>::iterator it=cleancol.begin(); it!=cleancol.end(); ++it )  
      	detIds.push_back( (*it).id() );
      chi2 = wls.chi2( cleancol, trackPar);
      ndof = cleancol.size();
      momentum = Vector(-1., 0., trackPar[1]).Cross( Vector(0., -1., trackPar[3]) );
      z0=layerZPosition.at(0);
      vertex = Point( trackPar[0]+trackPar[1]*z0, trackPar[2]+trackPar[3]*z0, z0 );
      track=reco::HGCalTBCaloTrack( chi2, ndof, vertex, momentum,/* cov,*/ detIds);
      //runPCA(cleancol);
    }
    else{
      for( std::vector<HGCalTBRecHit>::iterator it=col.begin(); it!=col.end(); ++it )  
     	detIds.push_back( (*it).id() );
      track=reco::HGCalTBCaloTrack( chi2, ndof, vertex, momentum,/* cov,*/ detIds);
    }
  
    _trackSuccess = (track.getDetIds().size()>=4) ? true:false;
    _tracknhit = track.getDetIds().size();
    _chi2 = track.normalisedChi2();
    _theta = track.momentum().theta()*180/PI;
    _phi = track.momentum().phi()*180/PI;
    _x0 = track.vertex().x();
    _y0 = track.vertex().y();
    _ax = track.momentum().x();
    _ay = track.momentum().y();    
    
    for( std::vector<HGCalTBDetId>::iterator it=track.getDetIds().begin(); it!=track.getDetIds().end(); ++it ){
      _nhitlayer[ (*it).layer()-1 ]++;
      HGCalTBRecHit hit=(*col.find(*it));
      _energylayer[ (*it).layer()-1  ]+=hit.energy();
      if( prepareTreeForDisplay ){
	_x.push_back( hit.x() );
	_y.push_back( hit.y() );
	_z.push_back( hit.z() );
	_energy.push_back( hit.energy() );
      }
      else if( _chi2 < maxChi2 ){
	uint32_t EID = essource_.emap_.detId2eid( *it );
	HGCalTBElectronicsId eid(EID);
	int key=( (*it).layer()-1 )*1000 + (eid.iskiroc() - 1)%2*100 + eid.ichan();
	h_mip_map[key]->Fill( hit.energy() );
      }
    }
  }
  tree->Fill();
  _evtID++;
}

void
TrackingExampleAnalyzer::beginJob()
{
  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(mapfile_);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  }
}

void
TrackingExampleAnalyzer::endJob()
{

}

void
TrackingExampleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackingExampleAnalyzer);
