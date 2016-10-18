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
  double cmThreshold ;
  std::vector< TH1F* > h_sum_layer;
  std::vector< TH1F* > h_layer_seven;
  std::vector< TH1F* > h_layer_nineteen;
  std::vector< TH1F* > h_cluster_layer;
  std::vector< TH1F* > h_clusterSize_layer;
  std::vector< TH1F* > h_7clusterSize_layer;
  std::vector< TH1F* > h_19clusterSize_layer;
  TH1F *h_sum_all, *h_cluster7E_all, *h_cluster19E_all, *h_clusterE_all, *h_cluster7N_all, *h_cluster19N_all, *h_clusterN_all ;
  
  int EVENT ;
  
  HGCalTBClustering *algo_HGCalTBClustering;
  SortByEnergy<reco::HGCalTBCluster,reco::HGCalTBCluster> energySorter;
  HGCalTBClusteringParameterSetting m_HGCalTBClusteringParameterSetting;
};


ClusteringExampleAnalyzer::ClusteringExampleAnalyzer(const edm::ParameterSet& iConfig) :
  nlayers( iConfig.getUntrackedParameter<int>("NLayers",8) ),
  cmThreshold( iConfig.getUntrackedParameter<int>("CMThreshold",30) )
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
  sensorsize=128;
  ALLCELLS_THRESHOLD=100;
  EVENT = 0;

  TH1F* h;
  for(int layer = 0; layer < nlayers; layer++) {
    stringstream name, sevenname, nineteenname, Xname, Yname, X_Y_name, clusterName, clusterSizeName, cluster19SizeName, cluster7SizeName, RMSX_name, RMSY_name, fullname, innername, outername, halfname, mousebitesname, mergedname;
    name << "AllCells_Sum_Layer" << layer ;
    sevenname << "Cells7_Sum_Layer" << layer ;
    nineteenname << "Cells19_Sum_Layer" << layer ;
    clusterName << "ClusterEnergy_Layer" << layer ;
    clusterSizeName << "ClusterSize_Layer" << layer ;
    cluster19SizeName << "19ClusterSize_Layer" << layer ;
    cluster7SizeName << "7ClusterSize_Layer" << layer ;
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
  EVENT++;
  double sensorEnergy = 0;
  double clusterEnergy = 0;
  double cluster7Energy = 0;
  double cluster19Energy = 0;
  int clusterNhit = 0;
  int cluster7Nhit = 0;
  int cluster19Nhit = 0;
  
  edm::Handle<HGCalTBRecHitCollection> Rechits;
  event.getByToken(HGCalTBRecHitCollection_, Rechits);

  std::map<int,HGCalTBRecHitCollection> colMap;

  double CM=0;
  int count=0;
  for( auto hit : *Rechits ){
    if( hit.id().cellType()!=0 || hit.energyHigh()<cmThreshold ) continue;
    count++;
    CM+=hit.energyHigh();
  }
  CM/=count;

  for( auto hit : *Rechits ){
    if( hit.energyHigh()-CM < cmThreshold ) continue;
    sensorEnergy+=hit.energyHigh()-CM;
    if( colMap.find( hit.id().layer()-1 )!=colMap.end() )
      colMap[ hit.id().layer()-1 ].push_back( hit );
    else {
      HGCalTBRecHitCollection coltmp;
      coltmp.push_back( hit );
      colMap[ hit.id().layer()-1 ]=coltmp;
    }
  }  

  for( std::map<int,HGCalTBRecHitCollection>::iterator it=colMap.begin(); it!=colMap.end(); ++it){
    double en=0;
    for( auto hit : it->second )
      en+=hit.energyHigh()-CM;
    h_sum_layer[ it->first ]->Fill( en ); 
   
    m_HGCalTBClusteringParameterSetting.maxTransverse=1;
    algo_HGCalTBClustering->SetHGCalTBClusteringParameterSetting(m_HGCalTBClusteringParameterSetting);
    std::vector<reco::HGCalTBCluster> outClusterCol;
    algo_HGCalTBClustering->Run(it->second,outClusterCol);
    if( outClusterCol.size()>0 ){
      std::sort( outClusterCol.begin(), outClusterCol.end(), energySorter.sort );
      if( (*outClusterCol.begin()).energyHigh()>ALLCELLS_THRESHOLD ) {
	clusterEnergy+=(*outClusterCol.begin()).energyHigh()-CM*(*outClusterCol.begin()).size();
	clusterNhit+=(*outClusterCol.begin()).size();
	h_cluster_layer[ it->first]->Fill( (*outClusterCol.begin()).energyHigh()-CM*(*outClusterCol.begin()).size() );
	h_clusterSize_layer[ it->first ]->Fill( (*outClusterCol.begin()).size() );
      }
    }
    reco::HGCalTBCluster cluster7;
    algo_HGCalTBClustering->RunSimple(it->second,cluster7);
    if( cluster7.size()>0 && cluster7.energyHigh()>ALLCELLS_THRESHOLD){
      cluster7Energy+=cluster7.energyHigh()-CM*cluster7.size();
      cluster7Nhit+=cluster7.size();
      h_layer_seven[ it->first ]->Fill( cluster7.energyHigh()-CM*cluster7.size() );
      h_7clusterSize_layer[ it->first ]->Fill( cluster7.size() );
    }

    reco::HGCalTBCluster cluster19;
    m_HGCalTBClusteringParameterSetting.maxTransverse=2;
    algo_HGCalTBClustering->SetHGCalTBClusteringParameterSetting(m_HGCalTBClusteringParameterSetting);
    algo_HGCalTBClustering->RunSimple(it->second,cluster19);
    if( cluster19.size()>0 && cluster19.energyHigh()>ALLCELLS_THRESHOLD){
      cluster19Energy+=cluster19.energyHigh()-CM*cluster19.size();
      cluster19Nhit+=cluster19.size();
      h_layer_nineteen[ it->first ]->Fill( cluster19.energyHigh()-CM*cluster19.size() );
      h_19clusterSize_layer[ it->first ]->Fill( cluster19.size() );
    }

  }
  
  if( sensorEnergy > ALLCELLS_THRESHOLD ){
    h_sum_all->Fill( sensorEnergy );
    h_clusterE_all->Fill( clusterEnergy );
    h_cluster19E_all->Fill( cluster19Energy );
    h_cluster7E_all->Fill( cluster7Energy );
    h_clusterN_all->Fill( clusterNhit );
    h_cluster19N_all->Fill( cluster19Nhit );
    h_cluster7N_all->Fill( cluster7Nhit );
  }  
}

void
ClusteringExampleAnalyzer::beginJob()
{
}

void
ClusteringExampleAnalyzer::endJob()
{

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
