#ifndef RECO_HGCALTBCALOTRACKFITTER
#define RECO_HGCALTBCALOTRACKFITTER 1

#include "HGCal/DataFormats/interface/HGCalTBCaloTrack.h"
#include "HGCal/DataFormats/interface/HGCalTBRecHitCollections.h"
#include "HGCal/DataFormats/interface/HGCalTBDetId.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Error.h"
#include "iostream"
#include "vector"

class HGCalTBCaloTrackFitter
{
 public:
  HGCalTBCaloTrackFitter(std::map<int,float> &map, int sensorsize=128);
  HGCalTBCaloTrackFitter(int setup=0, int sensorsize=128);
  ~HGCalTBCaloTrackFitter(){;}
  void Run(reco::HGCalTBCaloTrack &track, HGCalTBRecHitCollection hitcol, bool Cleaning=false, double maxDistance=1.0 );
  void deltaDistances( std::vector<double> &distances, HGCalTBRecHitCollection hitcol, reco::HGCalTBCaloTrack &track);

 private:
  void cleaning( HGCalTBRecHitCollection hitcol, HGCalTBRecHitCollection &cleancol, reco::HGCalTBCaloTrack &track, double maxDistance );

  int sensorSize; //default=128, needed for cell centre coordinate
  std::map<int,float> layerZPosition;

  void LeastSquare(HGCalTBRecHitCollection hitcol);
  
  double chi2;
  int ndof;
  Point vertex;
  Vector momentum;
  //CovarianceMatrix cov;
};

#endif
