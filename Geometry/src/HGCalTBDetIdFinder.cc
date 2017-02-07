#include "HGCal/Geometry/interface/HGCalTBDetIdFinder.h"
#include "HGCal/Geometry/interface/HGCalTBCellParameters.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "cmath"

HGCalTBDetIdFinder::HGCalTBDetIdFinder( const HGCalElectronicsMap& map, int size ) : emap(map),
										     sensorSize(size)
{;}

bool HGCalTBDetIdFinder::insideHexagon( std::pair<double,double> &xypoint, std::pair<double,double> &hexacentre, double hexalength)
{
  float x=(xypoint.first-hexacentre.first)/hexalength;
  float y=(xypoint.second-hexacentre.second)/hexalength;
  float l2 = x * x + y * y;
  if (l2 > 1.0f) return false;
  if (l2 < 0.75f) return true; // (sqrt(3)/2)^2 = 3/4
  // Check borders
  float py = y * 1.15470053838f; // 2/sqrt(3)
  if (py > 1.0f || py < -1.0f) return false;
  float px = 0.5f * py + x;
  if (px > 1.0f || px < -1.0f) return false;
  if (px - py > 1.0f || px - py < -1.0f) return false;
  return true;
}


void HGCalTBDetIdFinder::run(HGCalTBDetId &id, std::pair<double,double> &xy, int layer)
{
  // check if xy is inside the module
  double l=11*HGCAL_TB_CELL::FULL_CELL_SIDE;
  if( sensorSize==256 )
    l=15*HGCAL_TB_CELL::FULL_CELL_SIDE;
  std::pair<double,double> moduleCentre(0.,0.);
  if( insideHexagon(xy,moduleCentre,l)==false )
    return;
  double u=xy.first/(3*HGCAL_TB_CELL::FULL_CELL_SIDE)-xy.second/(std::sqrt(3)*HGCAL_TB_CELL::FULL_CELL_SIDE);
  double v=-2*xy.first/(3*HGCAL_TB_CELL::FULL_CELL_SIDE);
  int iu=std::floor(u);
  int iv=std::floor(v);
  int iU=0;
  int iV=0;
  HGCalTBCellVertices cells;
  std::pair<double,double> centre=cells.GetCellCentreCoordinates(layer,iU,iV,iu,iv,sensorSize);
  if( insideHexagon(xy,centre,HGCAL_TB_CELL::FULL_CELL_SIDE)==true ){
    buildDetId(id, xy, layer, iU, iV, iu, iv);
    return;
  }
  centre=cells.GetCellCentreCoordinates(layer,iU,iV,iu+1,iv,sensorSize);
  if( insideHexagon(xy,centre,HGCAL_TB_CELL::FULL_CELL_SIDE)==true ){
    buildDetId(id, xy, layer, iU, iV, iu+1, iv);
    return;
  }
  centre=cells.GetCellCentreCoordinates(layer,iU,iV,iu+1,iv+1,sensorSize);
  if( insideHexagon(xy,centre,HGCAL_TB_CELL::FULL_CELL_SIDE)==true ){
    buildDetId(id, xy, layer, iU, iV, iu+1, iv+1);
    return;
  }
  centre=cells.GetCellCentreCoordinates(layer,iU,iV,iu,iv+1,sensorSize);
  if( insideHexagon(xy,centre,HGCAL_TB_CELL::FULL_CELL_SIDE)==true ){
    buildDetId(id, xy, layer, iU, iV, iu, iv+1);
    return;
  }
}

void HGCalTBDetIdFinder::buildDetId( HGCalTBDetId &id, std::pair<double,double> xy, int layer, int iU, int iV, int iu, int iv )
{
for(unsigned int cellType=0; cellType<6; cellType++){
  id=HGCalTBDetId(layer, iU, iV, iu, iv, cellType);
  if( emap.existsDetId(id)==true ){
    if( cellType==1 ){
      HGCalTBCellVertices cells;
      std::pair<double,double> centre=cells.GetCellCentreCoordinates(layer,iU,iV,iu,iv,sensorSize);
      if( insideHexagon(xy,centre,HGCAL_TB_CELL::CALIB_PAD_SIDE)==false )
	id=HGCalTBDetId(layer, iU, iV, iu, iv, 4);
    }
    break;
  }
 }
}
