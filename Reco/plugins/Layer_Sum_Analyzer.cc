/* Need full layer, 7 cell cluster, and 19 cell cluster histograms for each layer
 * also need each for all layers summed
 * use ADC to MIP conversion of 1 MIP = 10 ADC Counts */

/**
   @Author: Ryan Quinn <ryan>
   7 July 2016
   quinn@physics.umn.edu
*/



// system include files
#include <memory>
#include <iostream>
#include "TH2Poly.h"
#include "TH1F.h"
#include "TF1.h"
#include <sstream>
#include <fstream>
#include <math.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHit.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/CondObjects/interface/HGCalCondObjectTextIO.h"
#include "HGCal/DataFormats/interface/HGCalTBElectronicsId.h"
#include "HGCal/DataFormats/interface/HGCalTBDataFrameContainers.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"
#include "HGCal/Geometry/interface/HGCalTBSpillParameters.h"

#include "HGCal/Reco/plugins/HGCalTBClustering.h"
#include "HGCal/DataFormats/interface/HGCalTBClusterCollection.h"
#include "HGCal/Reco/interface/HGCalTBSortingHelper.h"
#include "HGCal/Reco/plugins/HGCalTBCommonModeSubstraction.h"

using namespace std;

class Layer_Sum_Analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{

public:
  explicit Layer_Sum_Analyzer(const edm::ParameterSet&);
  ~Layer_Sum_Analyzer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  void analyze(const edm::Event& , const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetToken HGCalTBRecHitCollection_;
  int sensorsize;
  std::pair<double, double> CellCentreXY;
  HGCalTBCellVertices TheCell;
  double ADCtoMIP[MAXLAYERS];
  int CMTHRESHOLD ;// anything less than this value is added to the commonmode sum
  double ALLCELLS_THRESHOLD ;
  

  TH1F *h_sum_layer[MAXLAYERS], *h_layer_seven[MAXLAYERS], *h_layer_nineteen[MAXLAYERS], *h_cluster_layer[MAXLAYERS], *h_clusterSize_layer[MAXLAYERS], *h_7clusterSize_layer[MAXLAYERS], *h_19clusterSize_layer[MAXLAYERS]; 
  TH1F *h_sum_all, *h_cluster7E_all, *h_cluster19E_all, *h_clusterE_all, *h_cluster7N_all, *h_cluster19N_all, *h_clusterN_all ;
  TH1F *h_x_layer[MAXLAYERS], *h_y_layer[MAXLAYERS];
  TH1F *h_fullcommon[MAXLAYERS];
  TH1F *h_innercalibcommon[MAXLAYERS];
  TH1F *h_outercalibcommon[MAXLAYERS];
  TH1F *h_halfcommon[MAXLAYERS];
  TH1F *h_mousebitescommon[MAXLAYERS];
  TH1F *h_mergedcommon[MAXLAYERS];
  TH2F *h_x_y_layer[MAXLAYERS], *h_rmsx_layer[MAXLAYERS], *h_rmsy_layer[MAXLAYERS];
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


Layer_Sum_Analyzer::Layer_Sum_Analyzer(const edm::ParameterSet& iConfig)
{

  // initialization
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  HGCalTBRecHitCollection_ = consumes<HGCalTBRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCALTBRECHITS"));
  sensorsize=128;
  CMTHRESHOLD=30;// anything less than this value is added to the commonmode sum
  ALLCELLS_THRESHOLD=100;
  SPILL = EVENT = LAYER = 0;

  //booking the histos
  for(int layer = 0; layer < MAXLAYERS; layer++) {
    stringstream name, sevenname, nineteenname, Xname, Yname, X_Y_name, clusterName, clusterSizeName, cluster19SizeName, cluster7SizeName, RMSX_name, RMSY_name, fullname, innername, outername, halfname, mousebitesname, mergedname;

    name << "AllCells_Sum_Layer" << layer + 1;
    sevenname << "Cells7_Sum_Layer" << layer + 1;
    nineteenname << "Cells19_Sum_Layer" << layer + 1;
    Xname << "X_Layer" << layer + 1;
    Yname << "Y_Layer" << layer + 1;
    X_Y_name<<"X_Y_Layer"<< layer + 1;
    RMSX_name<<"rms_x_Layer"<< layer + 1;
    RMSY_name<<"rms_y_Layer"<< layer + 1;
    clusterName << "ClusterEnergy_Layer" << layer + 1;
    clusterSizeName << "ClusterSize_Layer" << layer + 1;
    cluster19SizeName << "19ClusterSize_Layer" << layer + 1;
    cluster7SizeName << "7ClusterSize_Layer" << layer + 1;
    fullname << "FullCellCommonMode" << layer + 1;
    innername << "InnerCalibCellCommonMode" << layer + 1;
    outername << "OuterCalibCellCommonMode" << layer + 1;
    halfname << "HalfCellCommonMode" << layer + 1;
    mousebitesname << "MouseBitesCellCommonMode" << layer + 1;
    mergedname << "MergedCellCommonMode" << layer + 1;

    h_sum_layer[layer] = fs->make<TH1F>(name.str().c_str(), name.str().c_str(), 40010, -10, 40000);
    h_layer_seven[layer] = fs->make<TH1F>(sevenname.str().c_str(), sevenname.str().c_str(), 40010, -10, 40000);
    h_layer_nineteen[layer] = fs->make<TH1F>(nineteenname.str().c_str(), nineteenname.str().c_str(), 40010, -10, 40000);
    h_x_layer[layer] = fs->make<TH1F>(Xname.str().c_str(), Xname.str().c_str(),2000,-10.,10. );
    h_y_layer[layer] = fs->make<TH1F>(Yname.str().c_str(), Yname.str().c_str(),2000,-10.,10. );
    h_x_y_layer[layer] = fs->make<TH2F>(X_Y_name.str().c_str(), X_Y_name.str().c_str(),2000,-10.,10.,2000,-10.,10. );
    h_rmsx_layer[layer] = fs->make<TH2F>(RMSX_name.str().c_str(), RMSX_name.str().c_str(),6000,0.,6000.,1000,0.,10. );
    h_rmsy_layer[layer] = fs->make<TH2F>(RMSY_name.str().c_str(), RMSY_name.str().c_str(),6000,0.,6000.,1000,0.,10. );
    h_cluster_layer[layer] = fs->make<TH1F>(clusterName.str().c_str(), clusterName.str().c_str(), 40010, -10, 40000);
    h_clusterSize_layer[layer] = fs->make<TH1F>(clusterSizeName.str().c_str(), clusterSizeName.str().c_str(), 100, 0, 100);
    h_19clusterSize_layer[layer] = fs->make<TH1F>(cluster19SizeName.str().c_str(), cluster19SizeName.str().c_str(), 100, 0, 100);
    h_7clusterSize_layer[layer] = fs->make<TH1F>(cluster7SizeName.str().c_str(), cluster7SizeName.str().c_str(), 100, 0, 100);

    h_fullcommon[layer] = fs->make<TH1F>(fullname.str().c_str(), fullname.str().c_str() , 200, -100, 100);
    h_innercalibcommon[layer] = fs->make<TH1F>(innername.str().c_str(), innername.str().c_str() , 200, -100, 100);
    h_outercalibcommon[layer] = fs->make<TH1F>(outername.str().c_str(), outername.str().c_str() , 200, -100, 100);
    h_halfcommon[layer] = fs->make<TH1F>(halfname.str().c_str(), halfname.str().c_str() , 200, -100, 100);
    h_mousebitescommon[layer] = fs->make<TH1F>(mousebitesname.str().c_str(), mousebitesname.str().c_str() , 200, -100, 100);
    h_mergedcommon[layer] = fs->make<TH1F>(mergedname.str().c_str(), mergedname.str().c_str() , 200, -100, 100);
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

}//constructor ends here


Layer_Sum_Analyzer::~Layer_Sum_Analyzer()
{
  delete algo_HGCalTBClustering;

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


// ------------ method called for each event  ------------
void
Layer_Sum_Analyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{


  if(((event.id()).event() - 1) % (EVENTSPERSPILL * MAXLAYERS) == 0 && (event.id()).event() != 1) {
    SPILL++;
  }
  LAYER = (((event.id()).event() - 1) / EVENTSPERSPILL) % MAXLAYERS;
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
  HGCalTBCommonModeSubstraction substraction( CMTHRESHOLD );
  substraction.Run( Rechits, coltmp );
  h_fullcommon[LAYER]->Fill( substraction.fullCommonMode() );
  h_innercalibcommon[LAYER]->Fill( substraction.innerCalibCommonMode() );
  h_outercalibcommon[LAYER]->Fill( substraction.outerCalibCommonMode() );
  h_halfcommon[LAYER]->Fill( substraction.halfCommonMode() );
  h_mousebitescommon[LAYER]->Fill( substraction.mouseBitesCommonMode() );
  h_mergedcommon[LAYER]->Fill( substraction.mergedCommonMode() );
  
  double x_tmp = 0., y_tmp = 0., x2_tmp = 0., y2_tmp = 0.;  
  float allcells_sum=0.;
  for(std::vector<HGCalTBRecHit>::iterator it=coltmp.begin(); it!=coltmp.end(); ++it) {
    HGCalTBRecHit hit=(*it);
    allcells_sum += hit.energyHigh();
    CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots( hit.id().layer(), 
							     hit.id().sensorIU(), 
							     hit.id().sensorIV(), 
							     hit.id().iu(), 
							     hit.id().iv(), 
							     sensorsize);
    x_tmp += CellCentreXY.first*hit.energyHigh() ;
    y_tmp += CellCentreXY.second*hit.energyHigh() ;
    x2_tmp += CellCentreXY.first*CellCentreXY.first*hit.energyHigh() ;
    y2_tmp += CellCentreXY.second*CellCentreXY.second*hit.energyHigh() ;
  }

  if( allcells_sum>0 ){
    x_tmp /= allcells_sum;
    y_tmp /= allcells_sum;
    x2_tmp = std::sqrt( x2_tmp/allcells_sum - x_tmp*x_tmp );
    y2_tmp = std::sqrt( y2_tmp/allcells_sum - y_tmp*y_tmp );
    h_x_layer[LAYER]->Fill(x_tmp);
    h_y_layer[LAYER]->Fill(y_tmp);
    h_x_y_layer[LAYER]->Fill(x_tmp,y_tmp);
    h_rmsx_layer[LAYER]->Fill( EVENT, x2_tmp );
    h_rmsy_layer[LAYER]->Fill( EVENT, y2_tmp );
    if( allcells_sum > ALLCELLS_THRESHOLD ){
      h_sum_layer[LAYER]->Fill( allcells_sum  ) ;
      sensorEnergy[EVENT] += allcells_sum;
    }
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
  //h_commonmode_all->Fill(commonmode);

}// analyze ends here


// ------------ method called once each job just before starting event loop  ------------
void
Layer_Sum_Analyzer::beginJob()
{
//  for(int iii = 0; iii < 16; iii++){
//    ADCtoMIP[iii] = 1;
//    ADCtoMIP[iii] = ADCtoMIP[iii] / 1.3; // Converting response to 120 GeV protons to MIPs
//  }
//  
}

// ------------ method called once each job just after ending the event loop  ------------
void
Layer_Sum_Analyzer::endJob()
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

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Layer_Sum_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Layer_Sum_Analyzer);
