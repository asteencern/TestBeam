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
  inline double fullCommonMode(){return fullCommonMode_;}
  inline double innerCalibCommonMode(){return innerCalibCommonMode_;}
  inline double outerCalibCommonMode(){return outerCalibCommonMode_;}
  inline double halfCommonMode(){return halfCommonMode_;}
  inline double mouseBitesCommonMode(){return mouseBitesCommonMode_;}
  inline double mergedCommonMode(){return mergedCommonMode_;}
 private:
  double threshold;
  double fullCommonMode_;
  double innerCalibCommonMode_;
  double outerCalibCommonMode_;
  double halfCommonMode_;
  double mouseBitesCommonMode_;
  double mergedCommonMode_;
};

#endif
