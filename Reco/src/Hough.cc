#include <HGCal/Reco/interface/Hough.h>
#include <HGCal/Geometry/interface/HGCalTBCellVertices.h>
#include <HGCal/Geometry/interface/HGCalTBCellParameters.h>
#include <HGCal/Reco/interface/HGCalTBCaloTrackingUtil.h>
#include <cmath>
#include <algorithm>

#define HOUGH_DEBUG

namespace reco
{
  
  const float PI = 3.1415927;

  Hough::Hough(HoughParameters params)
  {
    m_params=params;
  }

  Hough::Hough(HGCalTBDetectorLayout layout)
  {
    m_layout=layout;
  }

  Hough::Hough(HoughParameters params, HGCalTBDetectorLayout layout)
  {
    m_params=params;
    m_layout=layout;
  }

  void Hough::run(HGCalTBRecHitCollection hits, std::vector<reco::HGCalTBCaloTrack>& trackCol)
  {
    std::set<HGCalTBDetId> detids;
    for( auto hit : hits )
      detids.insert(hit.id());
    run(detids, trackCol);
  }
  void Hough::run(HGCalTBClusterCollection clusters, std::vector<reco::HGCalTBCaloTrack>& trackCol)
  {
    std::set<HGCalTBDetId> detids;
    for( auto cluster : clusters )
      detids.insert(HGCalTBDetId(cluster.seed()));
    run(detids, trackCol);
  }
  void Hough::run(std::set<HGCalTBDetId> detids, std::vector<reco::HGCalTBCaloTrack>& trackCol)
  {
    std::set<HoughObject> hgObjects;
    createHoughObjects(detids,hgObjects);
    HoughSpace hgSpace;
    fillHoughSpace(hgObjects,hgSpace,z_x);
    int theta_x(-1),rho_x(-1),theta_y(-1),rho_y(-1);
    while(1){
      std::set<HGCalTBDetId> maxBinDetIds;
      findMaxHoughBin(hgSpace, maxBinDetIds, theta_x, rho_x);
      if( maxBinDetIds.size()>m_params.minimumNBins ){
#ifdef DEBUG 
	std::cout << "theta = " << theta_x << "\t rho = " << rho_x << "\t nhit = " << maxBinDetIds.size() << std::endl;
	for( std::set<HGCalTBDetId>::iterator it=maxBinDetIds.begin(); it!=maxBinDetIds.end(); ++it )
	  std::cout << (*it) << std::endl;
#endif
	std::set<HoughObject> hgObjectsZY;
	for(std::set<HoughObject>::iterator it=hgObjects.begin(); it!=hgObjects.end(); ++it)
	  if( std::find( maxBinDetIds.begin(), maxBinDetIds.end(), (*it).id )!=maxBinDetIds.end() )
	    hgObjectsZY.insert(*it);
	HoughSpace hgSpaceZY;
	fillHoughSpace(hgObjectsZY,hgSpaceZY,z_y);
	std::set<HGCalTBDetId> maxBinDetIdsZY;
	findMaxHoughBin(hgSpaceZY, maxBinDetIdsZY, theta_y, rho_y);
#ifdef DEBUG 
	std::cout << "theta = " << theta_y << "\t rho = " << rho_y << "\t nhit = " << maxBinDetIdsZY.size() << std::endl;
	for( std::set<HGCalTBDetId>::iterator it=maxBinDetIdsZY.begin(); it!=maxBinDetIdsZY.end(); ++it )
	  std::cout << (*it) << std::endl;
#endif
	if( maxBinDetIdsZY.size()>m_params.minimumNBins ){
	  reco::HGCalTBCaloTrack track;
	  createTrack(maxBinDetIdsZY,track);
	  if( !track.isNull() && track.normalisedChi2()<m_params.maxChi2 ){
#ifdef DEBUG 
	    std::cout << "One track haas been created : \t" << track << std::endl;
#endif
	    trackCol.push_back(track);
	    removeObjectsFromHgSpace(track.getDetIds(),hgSpace);
	  }
	}
	maxBinDetIds.clear();
      }
      else break;
    }
  }

  void Hough::createHoughObjects(std::set<HGCalTBDetId> detids, std::set<HoughObject>& hgObjects)
  {
    std::pair<double, double> CellCentreXY;
    HGCalTBCellVertices TheCell;
    for( std::set<HGCalTBDetId>::iterator it=detids.begin(); it!=detids.end(); ++it ){
      HoughObject ho;
      ho.id=(*it);
      CellCentreXY = TheCell.GetCellCentreCoordinatesForPlots((*it).layer(),(*it).sensorIU(),(*it).sensorIV(),(*it).iu(),(*it).iv(),m_params.sensorSize);
      float x = CellCentreXY.first ;
      float y = CellCentreXY.second;
      float z = (m_layout.at( (*it).layer()-1 )).z();  
      for( int itheta=0; itheta<nThetas; itheta++ ){
      	ho.rho_zx[itheta]=z*std::cos(itheta*PI/nThetas) + x*std::sin(itheta*PI/nThetas);
      	ho.rho_zy[itheta]=z*std::cos(itheta*PI/nThetas) + y*std::sin(itheta*PI/nThetas);
      }
      hgObjects.insert(ho);
    }
  }

  void Hough::fillHoughSpace(std::set<HoughObject>& hgObjects, HoughSpace& hgSpace, projection proj)
  {
    int rho=0;
    for( std::set<HoughObject>::iterator it=hgObjects.begin(); it!=hgObjects.end(); ++it ){
      for( int itheta=0; itheta<nThetas; itheta++ ){
    	if( proj==z_x ) rho=int( (*it).rho_zx[itheta] );
    	else rho=int( (*it).rho_zy[itheta] );
	if( rho<maxRho )
	  hgSpace.table[itheta][rho].insert( (*it).id );
	else
	  std::cout << "Find rho = " << (*it).rho_zy[itheta] << " higher than maxRho = " << maxRho << std::endl;
      }
    }
  }

  void Hough::findMaxHoughBin(HoughSpace& hgSpace, std::set<HGCalTBDetId> &maxBinDetIds, int theta, int rho)
  {
#ifdef DEBUG
    if( maxBinDetIds.size()!=0 )
      std::cout << "Problem in findMaxHoughBin where maxBinDetIds should be empty : maxBinDetIds.size() = " << maxBinDetIds.size() << std::endl;
#endif
    for( int itheta=0; itheta<nThetas; itheta++ )
      for( int irho=0; irho<maxRho; irho++ )
    	if( hgSpace.table[itheta][irho].size()>maxBinDetIds.size() ){
    	  maxBinDetIds=hgSpace.table[itheta][irho];
	  theta=itheta;
	  rho=irho;
	}
    for( int itheta=theta-m_params.tolerance; itheta<=theta+m_params.tolerance; itheta++ )
      for( int irho=rho-m_params.tolerance; irho<=rho+m_params.tolerance; irho++ )
	for( std::set<HGCalTBDetId>::iterator it=hgSpace.table[itheta][irho].begin(); it!=hgSpace.table[itheta][irho].end(); ++it )
	  maxBinDetIds.insert(*it);
  }

  void Hough::removeObjectsFromHgSpace( std::vector<HGCalTBDetId> &detIds, HoughSpace& hgSpace)
  {
    for( std::vector<HGCalTBDetId>::iterator it=detIds.begin(); it!=detIds.end(); ++it )
      for( int itheta=0; itheta<nThetas; itheta++ )
    	for( int irho=0; irho<maxRho; irho++ )
    	  if( hgSpace.table[itheta][irho].find(*it)!=hgSpace.table[itheta][irho].end() )
    	    hgSpace.table[itheta][irho].erase(*it);
  }

  float Hough::distanceBetween2DetIds(HGCalTBDetId id0, HGCalTBDetId id1)
  {
    HGCalTBCellVertices m_cell;
    std::pair<double,double> xy0= m_cell.GetCellCentreCoordinatesForPlots(id0.layer(), id0.sensorIU(), id0.sensorIV(), id0.iu(), id0.iv(), m_params.sensorSize );
    std::pair<double,double> xy1= m_cell.GetCellCentreCoordinatesForPlots(id1.layer(), id1.sensorIU(), id1.sensorIV(), id1.iu(), id1.iv(), m_params.sensorSize );
    float z0 = (m_layout.at( id0.layer()-1 )).z();
    float z1 = (m_layout.at( id1.layer()-1 )).z();
    math::XYZPoint p0(xy0.first,xy0.second,z0);
    math::XYZPoint p1(xy1.first,xy1.second,z1);
    return std::sqrt( (p0-p1).mag2() );
  }

  void Hough::createTrack(std::set<HGCalTBDetId> &detids, reco::HGCalTBCaloTrack &track)
  {
    std::vector<math::XYZPoint> vec;
    std::vector<HGCalTBDetId> trackIds;
    for(std::set<HGCalTBDetId>::iterator it=detids.begin(); it!=detids.end(); ++it){
      HGCalTBCellVertices m_cell;
      std::pair<double,double> xy= m_cell.GetCellCentreCoordinatesForPlots((*it).layer(), (*it).sensorIU(), (*it).sensorIV(), (*it).iu(), (*it).iv(), m_params.sensorSize );
      float z = (m_layout.at( (*it).layer()-1 )).z();
      math::XYZPoint xyz(xy.first,xy.second,z);
      bool isolatedHit=true;
      for(std::set<HGCalTBDetId>::iterator jt=detids.begin(); jt!=detids.end(); ++jt){
	if( it==jt ) continue;
	if( std::abs( (*it).layer()-(*jt).layer() )<m_params.maxNumberOfLayerBetween2Hits ){
	  isolatedHit=false;
	  break;
	}
      }
      if(!isolatedHit){
	vec.push_back(xyz);
	trackIds.push_back(*it);
      }
    }
    if( vec.size()<m_params.minimumNBins ) return;
    reco::LeastSquare< std::vector<math::XYZPoint> > ls;
    std::vector<float> trackPar;
    std::vector<float> trackParError;
    ls.run( vec, trackPar, trackParError);
    float chi2 = ls.chi2( vec, trackPar);
    int ndof = vec.size();
    math::XYZVector momentum = Vector(-1., 0., trackPar[1]).Cross( Vector(0., -1., trackPar[3]) );
    math::XYZPoint vertex = Point( trackPar[0], trackPar[2], 0.0 );
    track=reco::HGCalTBCaloTrack( chi2, ndof, vertex, momentum,trackIds);    
  }
}
