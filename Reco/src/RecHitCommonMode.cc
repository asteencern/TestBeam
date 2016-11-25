#include "HGCal/Reco/interface/RecHitCommonMode.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "sstream"
#include "iostream"

RecHitCommonMode::RecHitCommonMode(HGCalElectronicsMap& emap,bool useHistoForFullCell) : _useHistoForFullCell( useHistoForFullCell )
{
  emap_ = emap;
  layerlist=emap.layersInMap();
  std::ostringstream os( std::ostringstream::ate );
  for( std::set<int>::const_iterator it=layerlist.begin(); it!=layerlist.end(); ++it){
    LayerCommonMode lcm;
    lcm.layerID=(*it);
    layerCMMap.insert( std::pair<int,LayerCommonMode>( (*it),lcm ) );
    if( _useHistoForFullCell ){
      os.str("");
      os << "FullCell_Layer_" << (*it);
      TH1F* h=new TH1F(os.str().c_str(),"",400,-200,200);
      histMap.insert( std::pair<int,TH1F*>( (*it),h ) );
    }
  }
  HGCalTBTopology topo;
  ratioForCalibPad=topo.Cell_Area(1)/topo.Cell_Area(0);
}


RecHitCommonMode::~RecHitCommonMode()
{
  layerCMMap.clear();
  for( std::map<int,TH1F*>::iterator it=histMap.begin(); it!=histMap.end(); ++it)
    delete it->second;
  histMap.clear();
}

void 
RecHitCommonMode::resetCM()
{
  for( std::map<int,LayerCommonMode>::iterator it=layerCMMap.begin(); it!=layerCMMap.end(); ++it){
    it->second.reset();
    if( _useHistoForFullCell ) histMap[it->first]->Reset();
  }
}

void
RecHitCommonMode::evaluate(HGCalTBRecHitCollection &hits, float maxEcut)
{
  resetCM();
  for( auto hit : hits ) {
    if(hit.energy() > maxEcut)continue;
    switch( hit.id().cellType() ){
    case 0:     
      if( _useHistoForFullCell )
	histMap[ hit.id().layer() ]->Fill(hit.energy());
      layerCMMap[ hit.id().layer() ].full+=hit.energy();
      layerCMMap[ hit.id().layer() ].fullCount++;
      break; 
    case 1: 
      //CM in inner calib pad is calculated using full CM and cells area 
      break; 
    case 2:
      layerCMMap[ hit.id().layer() ].half+=hit.energy();
      layerCMMap[ hit.id().layer() ].halfCount++;
      break; 
    case 3:
      layerCMMap[ hit.id().layer() ].mousebite+=hit.energy();
      layerCMMap[ hit.id().layer() ].mousebiteCount++;
      break; 
    case 4:
      // outer calib pads look like full pad
      if( _useHistoForFullCell )
	histMap[ hit.id().layer() ]->Fill(hit.energy());
      layerCMMap[ hit.id().layer() ].full+=hit.energy();
      layerCMMap[ hit.id().layer() ].fullCount++;
      break;
    case 5:
      layerCMMap[ hit.id().layer() ].merged+=hit.energy();
      layerCMMap[ hit.id().layer() ].mergedCount++;
      break;
    default:
      throw cms::Exception("InvalidCellType") << "rechit celltype if out of range";
      break;
    }               
  }
	
  for( std::map<int,LayerCommonMode>::iterator it=layerCMMap.begin(); it!=layerCMMap.end(); ++it){
    if( !_useHistoForFullCell || histMap[it->first]->GetEntries()==0 ) 
      it->second.fitResult=0.;
    else{
      int fitstatus = histMap[it->first]->Fit("gaus", "Q");
      if(fitstatus ==0) it->second.fitResult = histMap[it->first]->GetFunction("gaus")->GetParameter(1);
      else it->second.fitResult = 0;
    }

    if(it->second.fullCount > 0) {
      it->second.full/=it->second.fullCount;
      it->second.calib = it->second.full*ratioForCalibPad;
    }
    if(it->second.halfCount > 0) it->second.half/=it->second.halfCount;
    if(it->second.mergedCount > 0) it->second.merged/=it->second.mergedCount;
    if(it->second.mousebiteCount > 0) it->second.mousebite/=it->second.mousebiteCount;
  }
}

//does not work -> after this block, hits energy is back to origin energy
void 
RecHitCommonMode::subtract(HGCalTBRecHitCollection &hits)
{
  for( auto hit : hits ){
    if( hit.id().cellType()==0 && _useHistoForFullCell )
      hit.setEnergy( hit.energy()-getGaussCommonModeNoise(hit.id()));
    else 
      hit.setEnergy( hit.energy()-getMeanCommonModeNoise(hit.id()));
  }
}


float 
RecHitCommonMode::getGaussCommonModeNoise(HGCalTBDetId id)
{
  if( checkLayer( id.layer() )==false ){
    std::cout << "Problem : no layer " << id.layer() << " in ElectronicsMap => return 0 for CM" << std::endl;
    return 0;
  }
  float CMNoise(0);
  switch( id.cellType() ){
  case 0: CMNoise = layerCMMap[ id.layer() ].fitResult; break;
  case 1: CMNoise = layerCMMap[ id.layer() ].calib; break;
  case 2: CMNoise = layerCMMap[ id.layer() ].half; break;
  case 3: CMNoise = layerCMMap[ id.layer() ].mousebite; break;
  case 4: CMNoise = layerCMMap[ id.layer() ].fitResult; break;
  case 5: CMNoise = layerCMMap[ id.layer() ].merged; break;
  default: CMNoise = 0; break;
  };      
  return CMNoise;
}

float 
RecHitCommonMode::getMeanCommonModeNoise(HGCalTBDetId id)
{
  if( checkLayer( id.layer() )==false ){
    std::cout << "Problem : no layer " << id.layer() << " in ElectronicsMap => return 0 for CM" << std::endl;
    return 0;
  }
  float CMNoise(0);
  switch( id.cellType() ){
  case 0: CMNoise = layerCMMap[ id.layer() ].full; break;
  case 1: CMNoise = layerCMMap[ id.layer() ].calib; break;
  case 2: CMNoise = layerCMMap[ id.layer() ].half; break;
  case 3: CMNoise = layerCMMap[ id.layer() ].mousebite; break;
  case 4: CMNoise = layerCMMap[ id.layer() ].full; break;
  case 5: CMNoise = layerCMMap[ id.layer() ].merged; break;
  default: CMNoise = 0; break;
  };      
  return CMNoise;
}

float 
RecHitCommonMode::getGaussCommonModeNoise(int layer, int type)
{
  if( checkLayer( layer )==false ){
    std::cout << "Problem : no layer " << layer << " in ElectronicsMap => return 0 for CM" << std::endl;
    return 0;
  }
  float CMNoise(0);
  switch( type ){
  case 0: CMNoise = layerCMMap[ layer ].fitResult; break;
  case 1: CMNoise = layerCMMap[ layer ].calib; break;
  case 2: CMNoise = layerCMMap[ layer ].half; break;
  case 3: CMNoise = layerCMMap[ layer ].mousebite; break;
  case 4: CMNoise = layerCMMap[ layer ].fitResult; break;
  case 5: CMNoise = layerCMMap[ layer ].merged; break;
  default: CMNoise = 0; break;
  };      
  return CMNoise;
}

float
RecHitCommonMode::getMeanCommonModeNoise(int layer, int type)
{
  if( checkLayer( layer )==false ){
    std::cout << "Problem : no layer " << layer << " in ElectronicsMap => return 0 for CM" << std::endl;
    return 0;
  }
  float CMNoise(0);
  switch( type ){
  case 0: CMNoise = layerCMMap[ layer ].full; break;
  case 1: CMNoise = layerCMMap[ layer ].calib; break;
  case 2: CMNoise = layerCMMap[ layer ].half; break;
  case 3: CMNoise = layerCMMap[ layer ].mousebite; break;
  case 4: CMNoise = layerCMMap[ layer ].full; break;
  case 5: CMNoise = layerCMMap[ layer ].merged; break;
  default: CMNoise = 0; break;
  };      
  return CMNoise;
}
