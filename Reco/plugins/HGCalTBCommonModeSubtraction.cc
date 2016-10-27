#include "HGCal/Reco/plugins/HGCalTBCommonModeSubtraction.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"

HGCalTBCommonModeSubtraction::HGCalTBCommonModeSubtraction(double thr) : threshold(thr) 
{
  fullCommonMode_ = 0. ;
  innerCalibCommonMode_ = 0. ;
  halfCommonMode_ = 0. ;
  mouseBitesCommonMode_ = 0. ;
  outerCalibCommonMode_ = 0. ;
  mergedCommonMode_ = 0. ;
}

void HGCalTBCommonModeSubtraction::Run( HGCalTBRecHitCollection &col )
{
  HGCalTBTopology topo;
  double fullCellArea=topo.Cell_Area(0);

  int fullCount=0;
  int innerCalibCount=0;
  int outerCalibCount=0;
  int halfCount=0;
  int mouseBitesCount=0;
  int mergedCount=0;

  for( std::vector<HGCalTBRecHit>::iterator it=col.begin(); it!=col.end(); ++it ){
    if( (*it).energyHigh()*fullCellArea/topo.Cell_Area( (*it).id().cellType() ) > threshold )
      continue;    
    
    if( (*it).id().cellType() == 0 ){ fullCommonMode_+=(*it).energyHigh(); fullCount++; }
    else if( (*it).id().cellType() == 1 ){ innerCalibCommonMode_+=(*it).energyHigh(); innerCalibCount++; }
    else if( (*it).id().cellType() == 2 ){ halfCommonMode_+=(*it).energyHigh(); halfCount++; }
    else if( (*it).id().cellType() == 3 ){ mouseBitesCommonMode_+=(*it).energyHigh(); mouseBitesCount++; }
    else if( (*it).id().cellType() == 4 ){ outerCalibCommonMode_+=(*it).energyHigh(); outerCalibCount++; }
    else if( (*it).id().cellType() == 5 ){ mergedCommonMode_+=(*it).energyHigh(); mergedCount++; }
  }

  fullCommonMode_ = ( fullCount>0 ) ? fullCommonMode_/fullCount : 0;
  innerCalibCommonMode_ = ( innerCalibCount>0 ) ? innerCalibCommonMode_/innerCalibCount : 0;
  halfCommonMode_ = ( halfCount>0 ) ? halfCommonMode_/halfCount : 0;
  mouseBitesCommonMode_ = ( mouseBitesCount>0 ) ? mouseBitesCommonMode_/mouseBitesCount : 0;
  outerCalibCommonMode_ = ( outerCalibCount>0 ) ? outerCalibCommonMode_/outerCalibCount : 0;
  mergedCommonMode_ = ( mergedCount>0 ) ? mergedCommonMode_/mergedCount : 0;

  
  for( std::vector<HGCalTBRecHit>::iterator it=col.begin(); it!=col.end(); ++it ){
    if( (*it).id().cellType() == 0 ){ (*it).setEnergy((*it).energy() - fullCommonMode_)           ; }
    else if( (*it).id().cellType() == 1 ){ (*it).setEnergy((*it).energy() - innerCalibCommonMode_); }
    else if( (*it).id().cellType() == 2 ){ (*it).setEnergy((*it).energy() - halfCommonMode_)      ; }
    else if( (*it).id().cellType() == 3 ){ (*it).setEnergy((*it).energy() - mouseBitesCommonMode_); }
    else if( (*it).id().cellType() == 4 ){ (*it).setEnergy((*it).energy() - outerCalibCommonMode_); }
    else if( (*it).id().cellType() == 5 ){ (*it).setEnergy((*it).energy() - mergedCommonMode_)    ; }
  }
}
