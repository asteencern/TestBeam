#include <memory>
#include <iostream>
#include "TH1F.h"
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

#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"
#include "HGCal/Geometry/interface/HGCalTBSpillParameters.h"

#include "HGCal/Reco/plugins/HGCalTBClustering.h"
#include "HGCal/Reco/interface/HGCalTBSortingHelper.h"

using namespace std;

class ClusteringExampleAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit ClusteringExampleAnalyzer(const edm::ParameterSet&);
  ~ClusteringExampleAnalyzer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetToken HGCalTBRecHitCollection_;
  int sensorsize;
  int nlayers;
  double ALLCELLS_THRESHOLD ;
  std::vector< TH1F* > h_sum_layer;
  std::vector< TH1F* > h_layer_seven;
  std::vector< TH1F* > h_layer_nineteen;
  std::vector< TH1F* > h_cluster_layer;
  std::vector< TH1F* > h_clusterSize_layer;
  std::vector< TH1F* > h_7clusterSize_layer;
  std::vector< TH1F* > h_19clusterSize_layer;
  TH1F *h_sum_all, *h_cluster7E_all, *h_cluster19E_all, *h_clusterE_all, *h_cluster7N_all, *h_cluster19N_all, *h_clusterN_all ;
  
  int SPILL , EVENT , LAYER ;
  std::map<int, double> sensorEnergy;
  std::map<int, double> clusterEnergy;
  std::map<int, double> cluster7Energy;
  std::map<int, double> cluster19Energy;
  std::map<int, int> clusterNhit;
  std::map<int, int> cluster7Nhit;
  std::map<int, int> cluster19Nhit;
  
  HGCalTBClustering *algo_HGCalTBClustering;
  SortByEnergy<reco::HGCalTBCluster,reco::HGCalTBCluster> energySorter;
  HGCalTBClusteringParameterSetting m_HGCalTBClusteringParameterSetting;
};


ClusteringExampleAnalyzer::ClusteringExampleAnalyzer(const edm::ParameterSet& iConfig) :
  nlayers( iConfig.getUntrackedParameter<int>("NLayers",8) )
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
  sensorsize=128;
  ALLCELLS_THRESHOLD=100;
  SPILL = EVENT = LAYER = 0;

  TH1F* h;
  for(int layer = 0; layer < nlayers; layer++) {
    stringstream name, sevenname, nineteenname, Xname, Yname, X_Y_name, clusterName, clusterSizeName, cluster19SizeName, cluster7SizeName, RMSX_name, RMSY_name, fullname, innername, outername, halfname, mousebitesname, mergedname;
    name << "AllCells_Sum_Layer" << layer + 1;
    sevenname << "Cells7_Sum_Layer" << layer + 1;
    nineteenname << "Cells19_Sum_Layer" << layer + 1;
    clusterName << "ClusterEnergy_Layer" << layer + 1;
    clusterSizeName << "ClusterSize_Layer" << layer + 1;
    cluster19SizeName << "19ClusterSize_Layer" << layer + 1;
    cluster7SizeName << "7ClusterSize_Layer" << layer + 1;
    h = fs->make<TH1F>(name.str().c_str(), name.str().c_str(), 40010, -10, 40000);
    h_sum_layer.push_back(h);
    h = fs->make<TH1F>(sevenname.str().c_str(), sevenname.str().c_str(), 40010, -10, 40000); 
    h_layer_seven.push_back(h);
    h = fs->make<TH1F>(nineteenname.str().c_str(), nineteenname.str().c_str(), 40010, -10, 40000);
    h_layer_nineteen.push_back(h);
    h = fs->make<TH1F>(clusterName.str().c_str(), clusterName.str().c_str(), 40010, -10, 40000);
    h_cluster_layer.push_back(h);
    h = fs->make<TH1F>(clusterSizeName.str().c_str(), clusterSizeName.str().c_str(), 100, 0, 100);
    h_clusterSize_layer.push_back(h);
    h = fs->make<TH1F>(cluster19SizeName.str().c_str(), cluster19SizeName.str().c_str(), 100, 0, 100);
    h_19clusterSize_layer.push_back(h);
    h = fs->make<TH1F>(cluster7SizeName.str().c_str(), cluster7SizeName.str().c_str(), 100, 0, 100);
    h_7clusterSize_layer.push_back(h);
  }
  h_sum_all = fs->make<TH1F>("AllCells_Sum_AllLayers", "AllCells_Sum_AllLayers", 40010, -10, 200000);
  h_sum_all->Sumw2(); 
  h_clusterE_all = fs->make<TH1F>("Cluster_Sum_AllLayers", "Cluster_Sum_AllLayers", 40010, -10, 200000); 
  h_clusterE_all->Sumw2(); 
  h_cluster19E_all = fs->make<TH1F>("Cluster19_Sum_AllLayers", "Cluster19_Sum_AllLayers", 40010, -10, 200000); 
  h_cluster19E_all->Sumw2(); 
  h_cluster7E_all = fs->make<TH1F>("Cluster7_Sum_AllLayers", "Cluster7_Sum_AllLayers", 40010, -10, 200000); 
  h_cluster7E_all->Sumw2(); 
  h_clusterN_all = fs->make<TH1F>("Nhit_Sum_AllLayers", "Nhit_Sum_AllLayers", 510, -10, 500); 
  h_clusterN_all->Sumw2(); 
  h_cluster19N_all = fs->make<TH1F>("Nhit19_Sum_AllLayers", "Nhit19_Sum_AllLayers", 510, -10, 500); 
  h_cluster19N_all->Sumw2(); 
  h_cluster7N_all = fs->make<TH1F>("Nhit7_Sum_AllLayers", "Nhit7_Sum_AllLayers", 510, -10, 500); 
  h_cluster7N_all->Sumw2(); 
  
  algo_HGCalTBClustering = new HGCalTBClustering();

}

ClusteringExampleAnalyzer::~ClusteringExampleAnalyzer()
{
  delete algo_HGCalTBClustering;
}

void
ClusteringExampleAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
  if(((event.id()).event() - 1) % (EVENTSPERSPILL * nlayers) == 0 && (event.id()).event() != 1) {
    SPILL++;
  }
  LAYER = (((event.id()).event() - 1) / EVENTSPERSPILL) % nlayers;
  EVENT = ((event.id()).event() - 1) % EVENTSPERSPILL + EVENTSPERSPILL * SPILL;

  if( LAYER == 0 ){
    sensorEnergy[EVENT]=0.;
    clusterEnergy[EVENT]=0.;
    cluster7Energy[EVENT]=0.;
    cluster19Energy[EVENT]=0.;
    clusterNhit[EVENT]=0;
    cluster7Nhit[EVENT]=0;
    cluster19Nhit[EVENT]=0;
  }
  
  edm::Handle<HGCalTBRecHitCollection> Rechits;
  event.getByToken(HGCalTBRecHitCollection_, Rechits);

  HGCalTBRecHitCollection coltmp;
  float allcells_sum=0.;
  for( auto hit : *Rechits ){
    coltmp.push_back(hit);
    allcells_sum += hit.energyHigh();
  }
  
  if( allcells_sum > ALLCELLS_THRESHOLD ){
    h_sum_layer[LAYER]->Fill( allcells_sum  ) ;
    sensorEnergy[EVENT] += allcells_sum;
  }
  
  m_HGCalTBClusteringParameterSetting.maxTransverse=1;
  algo_HGCalTBClustering->SetHGCalTBClusteringParameterSetting(m_HGCalTBClusteringParameterSetting);
  std::vector<reco::HGCalTBCluster> outClusterCol;
  algo_HGCalTBClustering->Run(coltmp,outClusterCol);
  if( outClusterCol.size()>0 ){
    std::sort( outClusterCol.begin(), outClusterCol.end(), energySorter.sort );
    if( (*outClusterCol.begin()).energyHigh()>ALLCELLS_THRESHOLD ) {
      clusterEnergy[EVENT] += (*outClusterCol.begin()).energyHigh() ;      
      clusterNhit[EVENT] += (*outClusterCol.begin()).size() ;
      h_cluster_layer[LAYER]->Fill( (*outClusterCol.begin()).energyHigh() ) ;
      h_clusterSize_layer[LAYER]->Fill( (*outClusterCol.begin()).size() ) ;
    }
  }

  reco::HGCalTBCluster cluster7;
  algo_HGCalTBClustering->RunSimple(coltmp,cluster7);
  if( cluster7.size()>0 && cluster7.energyHigh()>ALLCELLS_THRESHOLD){
    cluster7Energy[EVENT] += cluster7.energyHigh() ;      
    cluster7Nhit[EVENT] += cluster7.size() ;
    h_layer_seven[LAYER]->Fill( cluster7.energyHigh() );
    h_7clusterSize_layer[LAYER]->Fill( cluster7.size() );
  }	

  reco::HGCalTBCluster cluster19;
  m_HGCalTBClusteringParameterSetting.maxTransverse=2;
  algo_HGCalTBClustering->SetHGCalTBClusteringParameterSetting(m_HGCalTBClusteringParameterSetting);
  algo_HGCalTBClustering->RunSimple(coltmp,cluster19);
  if( cluster19.size()>0 && cluster19.energyHigh()>ALLCELLS_THRESHOLD){
    cluster19Energy[EVENT] += cluster19.energyHigh() ;      
    cluster19Nhit[EVENT] += cluster19.size() ;
    h_layer_nineteen[LAYER]->Fill( cluster19.energyHigh() );
    h_19clusterSize_layer[LAYER]->Fill( cluster19.size() );
  }
}

void
ClusteringExampleAnalyzer::beginJob()
{
}

void
ClusteringExampleAnalyzer::endJob()
{

  for(int event = 0; event < (SPILL + 1) * EVENTSPERSPILL; event++) {
    if( sensorEnergy.find(event)!=sensorEnergy.end() ) h_sum_all->Fill( sensorEnergy[event] );
    if( clusterEnergy.find(event)!=clusterEnergy.end() ) h_clusterE_all->Fill( clusterEnergy[event] );
    if( cluster19Energy.find(event)!=cluster19Energy.end() ) h_cluster19E_all->Fill( cluster19Energy[event] );
    if( cluster7Energy.find(event)!=cluster7Energy.end() ) h_cluster7E_all->Fill( cluster7Energy[event] );
    if( clusterNhit.find(event)!=clusterNhit.end() ) h_clusterN_all->Fill( clusterNhit[event] );
    if( cluster19Nhit.find(event)!=cluster19Nhit.end() ) h_cluster19N_all->Fill( cluster19Nhit[event] );
    if( cluster7Nhit.find(event)!=cluster7Nhit.end() ) h_cluster7N_all->Fill( cluster7Nhit[event] );
  }
}

void
ClusteringExampleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ClusteringExampleAnalyzer);
