#ifndef HGCAL_RECO_RECHITCOMMONMODE_H
#define HGCAL_RECO_RECHITCOMMONMODE_H
#include <TH1F.h>
#include <TF1.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "HGCal/CondObjects/interface/HGCalElectronicsMap.h"


struct LayerCommonMode{
  int layerID;
  float full;
  float calib;
  float merged;
  float half;
  float mousebite;
  int fullCount;
  int mergedCount;
  int halfCount;
  int mousebiteCount;
  float fitResult;
  void reset()
  {
    full=0.;
    calib=0.;
    merged=0.;
    half=0.;
    mousebite=0.;
    fullCount=0;
    mergedCount=0;
    halfCount=0;
    mousebiteCount=0;
  }
};

class RecHitCommonMode
{

 public:
  RecHitCommonMode(HGCalElectronicsMap& emap,bool useHistoForFullCell=false);
  ~RecHitCommonMode();
  void evaluate(HGCalTBRecHitCollection &hits, float maxEcut = 1000.0);
  void subtract(HGCalTBRecHitCollection &hits);
  float getGaussCommonModeNoise(HGCalTBDetId id);
  float getGaussCommonModeNoise(int ilayer, int type);
  float getMeanCommonModeNoise(HGCalTBDetId id);
  float getMeanCommonModeNoise(int ilayer, int type);
        
 private:

  HGCalElectronicsMap emap_;
  std::set<int> layerlist;
  float ratioForCalibPad;
  bool _useHistoForFullCell;

  std::map<int,LayerCommonMode> layerCMMap;
  std::map<int,TH1F*> histMap;
  
  void resetCM();
  inline bool checkLayer(int layer){ return layerCMMap.find(layer)!=layerCMMap.end(); }
};

#endif
