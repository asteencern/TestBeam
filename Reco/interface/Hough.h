#ifndef HOUGH_HH
#define HOUGH_HH

#include <HGCal/DataFormats/interface/HGCalTBRecHitCollections.h>
#include <HGCal/DataFormats/interface/HGCalTBClusterCollection.h>
#include <HGCal/DataFormats/interface/HGCalTBTrack.h>
#include <HGCal/DataFormats/interface/HGCalTBDetId.h>
#include <HGCal/Reco/interface/HGCalTBCaloTrackingUtil.h>
#include "HGCal/CondObjects/interface/HGCalTBDetectorLayout.h"

#include <vector>
#include <set>
#include <iostream>

namespace reco
{

  static const int nThetas=100;
  static const int maxRho=40;//!!should be fine when no more than seven 6" modules per layer!!
  
  struct HoughParameters
  {
  HoughParameters() : minimumNBins(4),
      tolerance(1), 
      maxNumberOfLayerBetween2Hits(2),
      sensorSize(128),
      maxChi2(10)
    {;}
    unsigned int minimumNBins;
    int tolerance;
    int maxNumberOfLayerBetween2Hits;
    int sensorSize;
    float maxChi2;
  };


  //ARNAUD j'aime ca: c'est toujours une allocation dynamic et je suis pas sur que ca ne crash pas
  struct HoughSpace{
    std::set<HGCalTBDetId> table[nThetas][maxRho];
  };
  
  struct HoughObject{
    HoughObject(){ for(int i=0;i<nThetas;i++){rho_zx[i]=0; rho_zy[i]=0;} }
    HGCalTBDetId id;
    float rho_zx[nThetas],rho_zy[nThetas];
    bool operator==(const HoughObject& ho) const
    {
      return id == ho.id ; 
    }
    bool operator<(const HoughObject& ho) const
    {
      return id.rawId() < ho.id.rawId();
    }
  };

  enum projection{z_x,z_y};
  
  class Hough{
  public:
    Hough(){;}
    Hough(HoughParameters params);
    Hough(HGCalTBDetectorLayout layout);
    Hough(HoughParameters params, HGCalTBDetectorLayout layout);
    ~Hough(){;}
    void run(HGCalTBRecHitCollection hits, std::vector<reco::HGCalTBCaloTrack>& trackCol);
    void run(HGCalTBClusterCollection clusters, std::vector<reco::HGCalTBCaloTrack>& trackCol);
    void run(std::set<HGCalTBDetId> detids, std::vector<reco::HGCalTBCaloTrack>& trackCol);
    void setParameters(HoughParameters params){m_params=params;}
    void setDetectorLayout(HGCalTBDetectorLayout layout){m_layout=layout;}
  private:
    void createHoughObjects(std::set<HGCalTBDetId> detids, std::set<HoughObject>& hgObjects);
    void fillHoughSpace(std::set<HoughObject>& hgObjects, HoughSpace& hgSpace, projection proj);
    void findMaxHoughBin(HoughSpace& hgSpace, std::set<HGCalTBDetId> &maxBinDetIds, int theta, int rho);
    void removeObjectsFromHgSpace( std::vector<HGCalTBDetId> &detIds, HoughSpace& hgSpace);
    float distanceBetween2DetIds(HGCalTBDetId id0, HGCalTBDetId id1);
    void createTrack(std::set<HGCalTBDetId> &detids, reco::HGCalTBCaloTrack &track);
    
    HoughParameters m_params;
    HGCalTBDetectorLayout m_layout;
  };
}

#endif
