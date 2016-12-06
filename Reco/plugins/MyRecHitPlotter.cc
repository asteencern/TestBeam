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

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"

#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"

#include "HGCal/Reco/interface/RecHitCommonMode.h"

using namespace std;

struct layerInfo{
  int id;
  double energy;
  double emax;
  double noiseMean;
  int noiseCount;
  inline void Reset()
  { 
    energy=0.;
    emax=0.;
    noiseMean=0.;
    noiseCount=0;
  }
};

struct cellInfo{
  int key;
  double noise;
  double signal;
  inline void Reset()
  { 
    noise=-1000.;
    signal=-1000.;
  }
};

class MyRecHitPlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit MyRecHitPlotter(const edm::ParameterSet&);
  ~MyRecHitPlotter();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetToken HGCalTBRecHitCollection_;
  std::string mapfile_ = "HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt";
  struct {
    HGCalElectronicsMap emap_;
  } essource_;  
  
  int _evtID;
  int _nhit;
  double _energytot;
  double signalminenergy;

  std::set<int> layersInMap;

  TTree* tree;
  std::map<int,cellInfo> cellInfos;
  std::map<int,layerInfo> layerInfos;
  std::map<int,double> layerCM;

  std::map<int,double> layerEnergy;
  std::map<int,double> layerMaxEnergy;

  std::map<int,double> layerNoise;
  std::map<int,int> layerNhitNoise;

  //std::map<int, TH1F*> h_noise_map; //one per cell
  //std::map<int, TH1F*> h_signal_map; //one per cell
  //
  //std::map<int, TH2F*> h_noise_vs_elayer; //one per layer
  //std::map<int, TH2F*> h_noise_vs_emax; //one per layer
  //
  //std::map<int, TH2F*> h_noise_vs_eventID; //one per layer

  RecHitCommonMode *rhcm;

};


MyRecHitPlotter::MyRecHitPlotter(const edm::ParameterSet& iConfig) :
  signalminenergy( iConfig.getUntrackedParameter<double>("signalMinEnergy",10.) )
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
  
  std::cout << iConfig.dump() << std::endl;

  _evtID=0;

  HGCalCondObjectTextIO io(0);
  edm::FileInPath fip(mapfile_);
  if (!io.load(fip.fullPath(), essource_.emap_)) {
    throw cms::Exception("Unable to load electronics map");
  }

  layersInMap = essource_.emap_.layersInMap();
  int nskirocsperlayer = essource_.emap_.skirocsInMap().size()/layersInMap.size();
  tree = fs->make<TTree>("tree", "HGCAL TB variables tree");
  tree->Branch( "evtID",&_evtID ); 
  tree->Branch( "nhit",&_nhit ); 
  tree->Branch( "energyTot",&_energytot ); 
  std::ostringstream os( std::ostringstream::ate );
  for( auto layer : layersInMap ){
    layerInfo lay;
    lay.id=layer;
    layerInfos.insert( std::pair<int,layerInfo>(layer,lay) );
    layerCM.insert( std::pair<int,double>(layer,0.) );
    os.str("");
    os << "energy_layer" << layer;
    tree->Branch( os.str().c_str(),&layerInfos[layer].energy); 
    os.str("");
    os << "energyMax_layer" << layer;
    tree->Branch( os.str().c_str(),&layerInfos[layer].emax); 
    os.str("");
    os << "noise_layer" << layer;
    tree->Branch( os.str().c_str(),&layerInfos[layer].noiseMean); 
    os.str("");
    os << "nhit_noise_layer" << layer;
    tree->Branch( os.str().c_str(),&layerInfos[layer].noiseCount); 
    os.str("");
    os << "cm_layer" << layer;
    tree->Branch( os.str().c_str(),&layerCM[layer]); 
    for( int iskiroc=0; iskiroc<nskirocsperlayer; iskiroc++ ){
      for( int ichannel=0; ichannel<64; ichannel++ ){
	int key=layer*1000+iskiroc*100+ichannel;
	cellInfo cell;
	cell.key=key;
	cellInfos.insert( std::pair<int,cellInfo>(key,cell) );
	os.str("");
	os << "noise_layer" << layer << "_skiroc" << iskiroc << "_channel" << ichannel;
	tree->Branch( os.str().c_str(),&cellInfos[key].noise); 
	os.str("");
	os << "signal_layer" << layer << "_skiroc" << iskiroc << "_channel" << ichannel;
	tree->Branch( os.str().c_str(),&cellInfos[key].signal); 
      }
    }
  }
  
  //for( auto layer : layersInMap ){
  //  os.str("");
  //  os << "Noise_vs_Elayer" << layer;
  //  TH2F* h2 = fs->make<TH2F>(os.str().c_str(), os.str().c_str(), 4000,0,40000, 1000, -50, 50);
  //  h_noise_vs_elayer.insert( std::pair<int,TH2F*>( layer,h2 ) );
  //  os.str("");
  //  os << "Noise_vs_EMax_layer" << layer;
  //  h2 = fs->make<TH2F>(os.str().c_str(), os.str().c_str(), 1500,0,15000, 1000, -50, 50);
  //  h_noise_vs_emax.insert( std::pair<int,TH2F*>( layer,h2 ) );
  //  os.str("");
  //  os << "Noise_vs_EventID_layer" << layer;
  //  h2 = fs->make<TH2F>(os.str().c_str(), os.str().c_str(), 6000,0,6000, 1000, -50, 50);
  //  h_noise_vs_eventID.insert( std::pair<int,TH2F*>( layer,h2 ) );
  //  
  //  for( int iskiroc=0; iskiroc<nskirocsperlayer; iskiroc++ ){
  //    for( int ichannel=0; ichannel<64; ichannel++ ){
  //	int key=layer*1000+iskiroc*100+ichannel;
  //	os.str("");
  //	os << "Layer" << layer << "_Skiroc" << iskiroc << "_Channel" << ichannel << "_Noise";
  //	TH1F *h1 = fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 2000, -100, 100 );
  //	h_noise_map.insert( std::pair<int,TH1F*>(key,h1) );
  //	os.str("");
  //	os << "Layer" << layer << "_Skiroc" << iskiroc << "_Channel" << ichannel << "_Signal";
  //	h1 = fs->make<TH1F>(os.str().c_str(), os.str().c_str(), 15100, -100, 15000 );
  //	h_signal_map.insert( std::pair<int,TH1F*>(key,h1) );
  //    }
  //  }
  //}
  rhcm = new RecHitCommonMode( essource_.emap_ );
}

MyRecHitPlotter::~MyRecHitPlotter()
{
  delete rhcm;
}

void
MyRecHitPlotter::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  edm::Handle<HGCalTBRecHitCollection> Rechits;
  event.getByToken(HGCalTBRecHitCollection_, Rechits);
  
  for( std::map<int,cellInfo>::iterator  it=cellInfos.begin(); it!=cellInfos.end(); ++it )
    it->second.Reset();
  for( std::map<int,layerInfo>::iterator  it=layerInfos.begin(); it!=layerInfos.end(); ++it ){
    it->second.Reset();
    layerCM[it->first]=0.0;
  }
  _nhit=0;
  _energytot=0.;
  
  HGCalTBRecHitCollection tmp=(*Rechits);
  rhcm->evaluate( tmp,30 );
  for( auto hit : *Rechits ){
    hit.setEnergy( hit.energy()-rhcm->getMeanCommonModeNoise(hit.id()) );
    HGCalTBElectronicsId eid=HGCalTBElectronicsId( essource_.emap_.detId2eid(hit.id()) );
    int key=1000*hit.id().layer() + 100*((eid.iskiroc()-1)%2) + eid.ichan();
    if( hit.energy() > signalminenergy ){
      _nhit++;
      _energytot+=hit.energy();
      cellInfos[key].signal=hit.energy();
      //h_signal_map[key]->Fill( hit.energy() );
      if( layerEnergy.find(hit.id().layer())==layerEnergy.end() ){
	layerEnergy[hit.id().layer()]=hit.energy();
	layerMaxEnergy[hit.id().layer()]=hit.energy();
      }
      else{
	layerEnergy[hit.id().layer()]+=hit.energy();
	if( hit.energy()>layerMaxEnergy[hit.id().layer()] )
	  layerMaxEnergy[hit.id().layer()]=hit.energy();
      }	
    }
    else{
      cellInfos[key].noise=hit.energy();
      //h_noise_map[key]->Fill( hit.energy() );
      if( layerNhitNoise.find(hit.id().layer())==layerNhitNoise.end() ){
	layerNhitNoise[hit.id().layer()]=1;
	layerNoise[hit.id().layer()]=hit.energy();
      }
      else{
	layerNhitNoise[hit.id().layer()]+=1;
	layerNoise[hit.id().layer()]+=hit.energy();
      }	
    }
    //std::cout << key << "\t" << h_signal_map[key]->GetEntries() << "\t" << h_noise_map[key]->GetEntries() << std::endl;
  }
  
  for( auto layer : layersInMap ){
    layerCM[layer] = rhcm->getMeanCommonModeNoise(layer,0);
    layerNoise[layer]/=layerNhitNoise[layer];
    // h_noise_vs_elayer[layer]->Fill(layerEnergy[layer],layerNoise[layer]);
    // h_noise_vs_emax[layer]->Fill(layerMaxEnergy[layer],layerNoise[layer]);
    // h_noise_vs_eventID[layer]->Fill(_evtID,layerNoise[layer]);
    layerInfos[layer].energy=layerEnergy[layer];
    layerInfos[layer].emax=layerMaxEnergy[layer];
    layerInfos[layer].noiseMean=layerNoise[layer];
    layerInfos[layer].noiseCount=layerNhitNoise[layer];
  }
  tree->Fill();
  _evtID++;
  layerEnergy.clear();
  layerMaxEnergy.clear();
  layerNoise.clear();
  layerNhitNoise.clear();

}

void
MyRecHitPlotter::beginJob()
{
}

void
MyRecHitPlotter::endJob()
{
  // for( std::map<int,TH1F*>::iterator it = h_signal_map.begin(); it != h_signal_map.end(); ++it ){
  //   std::cout << it->first/1000 << "\t" << (it->first%1000)/100 << "\t" << it->first%100 << "\t" 
  // 	      << it->second->GetEntries() << "\t" << h_noise_map[it->first]->GetEntries() << std::endl;
  // }
  
}

void
MyRecHitPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyRecHitPlotter);
