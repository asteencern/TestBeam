#ifndef HGCAL_GEOMETRY_HGCALTBDETIDFINDER_H
#define HGCAL_GEOMETRY_HGCALTBDETIDFINDER_H 1

#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"

class HGCalTBDetIdFinder
{
 public:
  HGCalTBDetIdFinder( const HGCalElectronicsMap& map, int size );
  ~HGCalTBDetIdFinder(){;}
  void run( HGCalTBDetId &id, std::pair<double,double> &xy, int layer );
 private:
  bool insideHexagon( std::pair<double,double> &xypoint, std::pair<double,double> &hexacentre, double hexalength );
  void buildDetId(HGCalTBDetId &id, std::pair<double,double> xy, int layer, int iU, int iV, int iu, int iv);
  HGCalElectronicsMap emap;
  int sensorSize;
};

#endif
