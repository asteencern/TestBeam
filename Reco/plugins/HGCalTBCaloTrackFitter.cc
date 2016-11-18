#include "HGCal/Reco/plugins/HGCalTBCaloTrackFitter.h"
#include "HGCal/Reco/plugins/Distance.h"
#include "HGCal/Geometry/interface/HGCalTBCellVertices.h"
#include "cmath"

HGCalTBCaloTrackFitter::HGCalTBCaloTrackFitter(std::map<int,float> &map, int sensorsize)
{  
  layerZPosition=map;
  sensorSize=sensorsize;
}

HGCalTBCaloTrackFitter::HGCalTBCaloTrackFitter(int setup, int sensorsize)
{
  sensorSize=sensorsize;
  if( setup==0 ){
    float sum=0.;
    layerZPosition[1]=0.0+sum;sum+=layerZPosition[1];
    layerZPosition[2]=5.35+sum;sum+=layerZPosition[2];
    layerZPosition[3]=5.17+sum;sum+=layerZPosition[3];
    layerZPosition[4]=3.92+sum;sum+=layerZPosition[4];
    layerZPosition[5]=4.08+sum;sum+=layerZPosition[5];
    layerZPosition[6]=1.15+sum;sum+=layerZPosition[6];
    layerZPosition[7]=4.11+sum;sum+=layerZPosition[7];
    layerZPosition[8]=2.14+sum;sum+=layerZPosition[8];
  }
  else if( setup==1 ){
    float sum=0.;
    layerZPosition[1]=0.0+sum;sum+=layerZPosition[1];
    layerZPosition[2]=4.67+sum;sum+=layerZPosition[2];
    layerZPosition[3]=5.17+sum;sum+=layerZPosition[3];
    layerZPosition[4]=4.43+sum;sum+=layerZPosition[4];
    layerZPosition[5]=4.98+sum;sum+=layerZPosition[5];
    layerZPosition[6]=1.15+sum;sum+=layerZPosition[6];
    layerZPosition[7]=5.40+sum;sum+=layerZPosition[7];
    layerZPosition[8]=5.60+sum;sum+=layerZPosition[8];
  }
  else HGCalTBCaloTrackFitter(0);
}

void HGCalTBCaloTrackFitter::Run( reco::HGCalTBCaloTrack &track, HGCalTBRecHitCollection hitcol, bool Cleaning, double maxDistance )
{
  std::vector<HGCalTBDetId> detIds;

  LeastSquare( hitcol );
  
  track=reco::HGCalTBCaloTrack( chi2, ndof, vertex,
				 momentum,/* cov,*/ detIds);

  if( Cleaning==true ){
    HGCalTBRecHitCollection cleancol;
    cleaning( hitcol, cleancol, track, maxDistance );
    LeastSquare( cleancol );
    for( std::vector<HGCalTBRecHit>::iterator it=cleancol.begin(); it!=cleancol.end(); ++it )  
      detIds.push_back( (*it).id() );
    track=reco::HGCalTBCaloTrack( chi2, ndof, vertex,
				  momentum,/* cov,*/ detIds);
    return;
  }
  else{
    for( std::vector<HGCalTBRecHit>::iterator it=hitcol.begin(); it!=hitcol.end(); ++it )  
      detIds.push_back( (*it).id() );
    track=reco::HGCalTBCaloTrack( chi2, ndof, vertex,
				  momentum,/* cov,*/ detIds);
    return;
  }
}

void HGCalTBCaloTrackFitter::LeastSquare( HGCalTBRecHitCollection hitcol )
{
  float _params[] = {0,0,0,0};
  //float _paramsError[] = {0,0,0,0}; //don't know yet if it has any interest to keep it
  
  HGCalTBCellVertices cellVertice;

  float xsum = 0.0;
  float ysum = 0.0;
  float zsum = 0.0;
  float zzsum = 0.0;
  float xzsum = 0.0;
  float yzsum = 0.0;
  float esum = 0.0;
  for( std::vector<HGCalTBRecHit>::iterator it=hitcol.begin(); it!=hitcol.end(); ++it ){  
    xsum = xsum + (*it).x()*(*it).energy();
    ysum = ysum + (*it).y()*(*it).energy();
    zsum = zsum + (*it).z()*(*it).energy();
    xzsum = xzsum + (*it).x()*(*it).z()*(*it).energy();
    yzsum = yzsum + (*it).y()*(*it).z()*(*it).energy();
    zzsum = zzsum + (*it).z()*(*it).z()*(*it).energy();
    esum = esum + (*it).energy();
  }

  _params[0] = (zzsum*xsum-xzsum*zsum)/(esum*zzsum-zsum*zsum);
  _params[2] = (zzsum*ysum-yzsum*zsum)/(esum*zzsum-zsum*zsum);

  _params[1] = (xzsum*esum-xsum*zsum)/(esum*zzsum-zsum*zsum);
  _params[3] = (yzsum*esum-ysum*zsum)/(esum*zzsum-zsum*zsum);
  
  //_paramsError[0] = std::sqrt( zzsum/(hitcol.size()*zzsum-zsum*zsum) );
  //_paramsError[2] = std::sqrt( zzsum/(hitcol.size()*zzsum-zsum*zsum) );
  //
  //_paramsError[1] = std::sqrt( hitcol.size()/(hitcol.size()*zzsum-zsum*zsum) );
  //_paramsError[3] = std::sqrt( hitcol.size()/(hitcol.size()*zzsum-zsum*zsum) );

  chi2 = 0;
  ndof = hitcol.size();
  for( std::vector<HGCalTBRecHit>::iterator it=hitcol.begin(); it!=hitcol.end(); ++it ){
    Point p_track(  _params[0]+_params[1]*(*it).z(), _params[2]+_params[3]*(*it).z(), (*it).z() );
    chi2 += (p_track-(*it).position()).Mag2()*12;
  }
  //  chi2 /= ndof;
  
  momentum = Vector(-1., 0., _params[1]).Cross( Vector(0., -1., _params[3]) );
  float z0=layerZPosition.begin()->second;
  vertex = Point( _params[0]+_params[1]*z0, _params[2]+_params[3]*z0, z0 );
  
}

void HGCalTBCaloTrackFitter::deltaDistances( std::vector<double> &distances, HGCalTBRecHitCollection hitcol, reco::HGCalTBCaloTrack &track )
{
  distances.clear();
  Distance<HGCalTBRecHit,reco::HGCalTBCaloTrack> dist;
  for( std::vector<HGCalTBDetId>::const_iterator it=track.getDetIds().begin(); it!=track.getDetIds().end(); ++it ){
    //HGCalTBRecHit hit=(*col.find(*it));
    //for( std::vector<HGCalTBRecHit>::iterator it=hitcol.begin(); it!=hitcol.end(); ++it )
    distances.push_back( dist.distance( (*hitcol.find(*it)),track ) );
  }
}

void HGCalTBCaloTrackFitter::cleaning( HGCalTBRecHitCollection hitcol, HGCalTBRecHitCollection &cleancol, reco::HGCalTBCaloTrack &track, double maxDistance )
{
  Distance<HGCalTBRecHit,reco::HGCalTBCaloTrack> dist;
  for( std::vector<HGCalTBRecHit>::iterator it=hitcol.begin(); it!=hitcol.end(); ++it )
    if( dist.distance( (*it),track )<maxDistance )
      cleancol.push_back( *it );
}
