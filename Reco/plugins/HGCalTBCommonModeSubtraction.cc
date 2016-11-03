#include "HGCal/Reco/plugins/HGCalTBCommonModeSubtraction.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"

HGCalTBCommonModeSubtraction::HGCalTBCommonModeSubtraction(double thr) : threshold(thr) 
{
  commonMode_ = 0. ;
}

void HGCalTBCommonModeSubtraction::Run( HGCalTBRecHitCollection &col )
{
  int count=0;

  for( std::vector<HGCalTBRecHit>::iterator it=col.begin(); it!=col.end(); ++it )
    if( (*it).energyHigh() < threshold && (*it).id().cellType() == 0 ){
      commonMode_+=(*it).energyHigh(); 
      count++;
    }

  commonMode_ = ( count>0 ) ? commonMode_/count : 0;

  HGCalTBTopology topo;
  for( std::vector<HGCalTBRecHit>::iterator it=col.begin(); it!=col.end(); ++it ){
    double ratio=topo.Cell_Area( (*it).id().cellType() )/topo.Cell_Area( 0 );
    (*it).setEnergy((*it).energy() - commonMode_*ratio);
  }
}
