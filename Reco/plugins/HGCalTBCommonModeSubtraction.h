#ifndef RECO_HGCALTBCOMMONMODESUBSTRACTION
#define RECO_HGCALTBCOMMONMODESUBSTRACTION 1

#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "iostream"
#include "vector"

class HGCalTBCommonModeSubtraction
{
 public:
 HGCalTBCommonModeSubtraction(double thr);
  ~HGCalTBCommonModeSubtraction(){;}
  void Run( HGCalTBRecHitCollection &col );
  inline double commonMode(){return commonMode_;}
 private:
  double threshold;
  double commonMode_;
};

#endif
