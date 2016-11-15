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

#include "HGCal/Reco/plugins/HGCalTBClustering.h"
#include "HGCal/Reco/interface/HGCalTBSortingHelper.h"
#include "HGCal/Reco/plugins/HGCalTBCaloTrackFitter.h"
#include "HGCal/Reco/plugins/HGCalTBCommonModeSubtraction.h"
#include "HGCal/Reco/plugins/Distance.h"

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
  int nlayers;
  int nskirocsperlayer;
  int sensorSize;
  int cmThreshold;
  double minEnergy;
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
  std::vector<double> _commonMode;
  std::vector<double> _specCommonMode;
  
  std::vector<int> _channelForCM;
  
  std::string mapfile_ = "HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt";
  struct {
    HGCalElectronicsMap emap_;
  } essource_;  
  const float PI = 3.1415927;
  HGCalTBCaloTrackFitter *algo_HGCalTBCaloTrackFitter;
  HGCalTBClustering *algo_HGCalTBClustering;
  SortByEnergy<reco::HGCalTBCluster,reco::HGCalTBCluster> energySorter;
  HGCalTBClusteringParameterSetting m_HGCalTBClusteringParameterSetting;
};

ShowerAnalyzer::ShowerAnalyzer(const edm::ParameterSet& iConfig) :
  nlayers( iConfig.getUntrackedParameter<int>("NLayers",8) ),
  nskirocsperlayer( iConfig.getUntrackedParameter<int>("NSkirocsPerLayer",2) ),
  sensorSize( iConfig.getUntrackedParameter<int>("SensorSize",128) ),
  cmThreshold( iConfig.getUntrackedParameter<int>("CMThreshold",30) ),
  minEnergy( iConfig.getUntrackedParameter<double>("minEnergy",10) ),
  CERN_8layers_config( iConfig.getUntrackedParameter<int>("CERN_8layers_config",0) ),
  energyShowerThreshold( iConfig.getUntrackedParameter<double>("energyShowerThreshold",200) ),
  maxTransverseProfile( iConfig.getUntrackedParameter<double>("maxTransverseProfile",20) )
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
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
  
  m_HGCalTBClusteringParameterSetting.maxTransverse=1;
  algo_HGCalTBClustering = new HGCalTBClustering();
  algo_HGCalTBClustering->SetHGCalTBClusteringParameterSetting(m_HGCalTBClusteringParameterSetting);
  
  algo_HGCalTBCaloTrackFitter = new HGCalTBCaloTrackFitter( CERN_8layers_config,sensorSize );

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
  tree->Branch( "commonMode","std::vector<double>",&_commonMode);
  tree->Branch( "specCommonMode","std::vector<double>",&_specCommonMode);
  
  int channel[]={1,2,9,11,13,14,21,27,31,37,41,48,49,50,55,57,62,63};
  int border[]={3,15,29,35,47,53,58,61};
  for( unsigned int i=0; i<sizeof(channel)/sizeof(int); i++ )
  _channelForCM.push_back(channel[i]);
  for( unsigned int i=0; i<sizeof(border)/sizeof(int); i++ )
    _channelForCM.push_back(border[i]);
}

ShowerAnalyzer::~ShowerAnalyzer()
{
  delete algo_HGCalTBClustering;
  delete algo_HGCalTBCaloTrackFitter;
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
  _commonMode.clear();
  _specCommonMode.clear();
  _energyInCluster=0;
  edm::Handle<HGCalTBRecHitCollection> Rechits;
  event.getByToken(HGCalTBRecHitCollection_, Rechits);

  HGCalTBRecHitCollection hitcol;
  HGCalTBCellVertices cellVertice;
  HGCalTBCommonModeSubtraction subtraction( cmThreshold );
  std::vector<reco::HGCalTBCluster> outClusterCol;
  
  //hard-coded : ok for 6" layers
  for( int ir=0; ir<maxTransverseProfile; ir++ ){
    _transverseprofile.push_back(0);
  }

  for( int ilayer=0; ilayer<nlayers; ilayer++){
    _energylayer.push_back(0.0);
    _clusterenergylayer.push_back(0.0);
    _clustersizelayer.push_back(0);
    _meanx.push_back(0.0);
    _meany.push_back(0.0);
    _rmsx.push_back(0.0);
    _rmsy.push_back(0.0);
    _commonMode.push_back(0.0);
    _specCommonMode.push_back(0.0);

    HGCalTBRecHitCollection hitcollayer;
    HGCalTBRecHitCollection coltmp;
    double cm=0;
    int count=1;
    for( auto hit : *Rechits ){
      if( hit.id().layer()-1!=ilayer ) continue;
      uint32_t EID = essource_.emap_.detId2eid( hit.id() );
      HGCalTBElectronicsId eid(EID);
      int chan=eid.ichan();
      if( std::find( _channelForCM.begin(),_channelForCM.end(),chan ) != _channelForCM.end() && hit.energy()<cmThreshold ){
	cm+=hit.energy();
       	count++;
	//if( ilayer==0 && eid.iskiroc()==1 )std::cout << chan << std::endl;
      }
      coltmp.push_back( hit );
    }
    _specCommonMode[ ilayer ] = cm/count;
    subtraction.Run( coltmp );
    _commonMode[ ilayer ] = subtraction.commonMode();
    for( std::vector<HGCalTBRecHit>::iterator it=coltmp.begin(); it!=coltmp.end(); ++it ){
      if( //(*it).id().cellType()==1 ||
	  (*it).id().cellType()==2 || 
	  (*it).id().cellType()==3 || 
	  //(*it).id().cellType()==4 || 
	  (*it).id().cellType()==5 )
	continue;
      if( (*it).energy() > minEnergy ){
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
	hitcollayer.push_back( (*it) );
	hitcol.push_back( (*it) );
      }
    }
    if( hitcollayer.size() > 0 ){
      std::vector<reco::HGCalTBCluster> clusters;
      algo_HGCalTBClustering->Run(hitcollayer,clusters);
      std::sort( clusters.begin(), clusters.end(), energySorter.sort );
      //outClusterCol.insert( outClusterCol.end(), clusters.begin(), clusters.end() );
      outClusterCol.push_back( (*clusters.begin()) );
      _energyInCluster+=(*clusters.begin()).energy();
      _clusterenergylayer[ ilayer ]=(*clusters.begin()).energy();
      _clustersizelayer[ ilayer ]=(*clusters.begin()).size();
    }
  }
  Distance<HGCalTBRecHit,reco::HGCalTBCaloTrack> dist;
  if( hitcol.size()>=4 ){
    reco::HGCalTBCaloTrack track;
    algo_HGCalTBCaloTrackFitter->Run( track,hitcol );
    if( !track.isNull() ){
      _theta = track.momentum().theta()*180/PI;
      _phi = track.momentum().phi()*180/PI;
      _x0 = track.vertex().x();
      _y0 = track.vertex().y();
      _ax = track.momentum().x();
      _ay = track.momentum().y();    
      for( std::vector<reco::HGCalTBCluster>::iterator it=outClusterCol.begin(); it!=outClusterCol.end(); ++it ){
	for( std::vector< std::pair<DetId,float> >::const_iterator jt=(*it).hitsAndFractions().begin(); jt!=(*it).hitsAndFractions().end(); ++jt ){
	  HGCalTBRecHit hit=(*hitcol.find(jt->first));
	  int d=(int)fabs( dist.distance( hit,track )/HGCAL_TB_CELL::FULL_CELL_SIDE );
	  if( d<maxTransverseProfile )
	    _transverseprofile.at(d)+=hit.energy();
	}
      }
    }
  }
  _nhit=hitcol.size();
  for( std::vector<HGCalTBRecHit>::iterator it=hitcol.begin(); it!=hitcol.end(); ++it ){
    _energylayer[ (*it).id().layer()-1 ] += (*it).energy();
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
