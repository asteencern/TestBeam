#include "HGCal/Reco/plugins/HGCalTBClustering.h"
#include "HGCal/Geometry/interface/HGCalTBTopology.h"
#include "HGCal/Geometry/interface/HGCalTBGeometryParameters.h"

#include "algorithm"

void HGCalTBClustering::Run(HGCalTBRecHitCollection hitcol, std::vector<reco::HGCalTBCluster> &outClusterColl)
{
  std::vector<HGCalTBDetId> temp;
  for( std::vector<HGCalTBRecHit>::iterator it=hitcol.begin(); it!=hitcol.end(); ++it){
    if( std::find(temp.begin(), temp.end(), (*it).id())!=temp.end() )
      continue;
    std::vector<HGCalTBDetId> clusterDetIDs;
    temp.push_back( (*it).id() );
    clusterDetIDs.push_back( (*it).id() );
    BuildCluster(hitcol, temp, clusterDetIDs);

    reco::HGCalTBCluster cluster;
    cluster.setLayer( clusterDetIDs.at(0).layer() );
    float energyHigh=0.;
    float energyLow=0.;
    float energy=0.;
    float x,y,z; 
    x = y = z = 0.0;
    for( std::vector<HGCalTBDetId>::iterator jt=clusterDetIDs.begin(); jt!=clusterDetIDs.end(); ++jt){
      energyHigh+=(*hitcol.find(*jt)).energyHigh();
      energyLow+=(*hitcol.find(*jt)).energyLow();
      energy+=(*hitcol.find(*jt)).energy();
      x += (*hitcol.find(*jt)).x()*(*hitcol.find(*jt)).energy();
      y += (*hitcol.find(*jt)).y()*(*hitcol.find(*jt)).energy();
      z += (*hitcol.find(*jt)).z()*(*hitcol.find(*jt)).energy();
    }
    cluster.setPosition( math::XYZPoint(x,y,z)/energy );
    cluster.setEnergyLow(energyLow);
    cluster.setEnergyHigh(energyHigh);
    cluster.setEnergy(energy);

    for( std::vector<HGCalTBDetId>::iterator jt=clusterDetIDs.begin(); jt!=clusterDetIDs.end(); ++jt)
      cluster.addHitAndFraction( (*jt), (*hitcol.find(*jt)).energy()/energy );

    outClusterColl.push_back(cluster);
  }
}


void HGCalTBClustering::BuildCluster(HGCalTBRecHitCollection hitcol,
				     std::vector<HGCalTBDetId> &temp,
				     std::vector<HGCalTBDetId> &clusterDetIDs)
{ 
  HGCalTBTopology top;
  HGCalTBDetId detID=clusterDetIDs.back();
  std::set<HGCalTBDetId> neighbors=top.getNeighboringCellsDetID( detID, settings.sensorSize , settings.maxTransverse );
  for( std::set<HGCalTBDetId>::const_iterator it=neighbors.begin(); it!=neighbors.end(); ++it){
    if( settings.useDistanceInsteadCellID ) continue;
    if( std::find(temp.begin(), temp.end(), (*it))!=temp.end() || hitcol.find(*it)==hitcol.end() )
      continue;
    temp.push_back( (*it) );
    clusterDetIDs.push_back( (*it) );
    BuildCluster(hitcol, temp, clusterDetIDs);
  }
}

void HGCalTBClustering::RunSimple(HGCalTBRecHitCollection hitcol, reco::HGCalTBCluster &cluster)
{
  if( hitcol.size()==0 ) return;
  HGCalTBTopology top;
  HGCalTBRecHit hitMax=(*hitcol.begin());
  for( std::vector<HGCalTBRecHit>::iterator it=hitcol.begin(); it!=hitcol.end(); ++it){
    if( (*it).energy() > hitMax.energy() )
      hitMax=(*it);
  }
  cluster.setLayer( hitMax.id().layer() );
  float energyHigh=hitMax.energyHigh();
  float energyLow=hitMax.energyLow();
  float energy=hitMax.energy();
  float x = hitMax.x()*hitMax.energy();
  float y = hitMax.y()*hitMax.energy();
  float z = hitMax.z()*hitMax.energy();

  std::set<HGCalTBDetId> neighbors=top.getNeighboringCellsDetID( hitMax.id(), settings.sensorSize , settings.maxTransverse );
  for( std::set<HGCalTBDetId>::iterator jt=neighbors.begin(); jt!=neighbors.end(); ++jt){
    if( hitcol.find(*jt) != hitcol.end() ){
      energyHigh+=(*hitcol.find(*jt)).energyHigh();
      energyLow+=(*hitcol.find(*jt)).energyLow();
      energy+=(*hitcol.find(*jt)).energy();
      x += (*hitcol.find(*jt)).x()*(*hitcol.find(*jt)).energy();
      y += (*hitcol.find(*jt)).y()*(*hitcol.find(*jt)).energy();
      z += (*hitcol.find(*jt)).z()*(*hitcol.find(*jt)).energy();
    }
  }
  cluster.setPosition( math::XYZPoint(x,y,z)/energy );
  cluster.setEnergyLow(energyLow);
  cluster.setEnergyHigh(energyHigh);
  cluster.setEnergy(energy);

  cluster.addHitAndFraction( hitMax.id(), hitMax.energy()/energy );
  for( std::set<HGCalTBDetId>::iterator jt=neighbors.begin(); jt!=neighbors.end(); ++jt)
    if( hitcol.find(*jt) != hitcol.end() )
      cluster.addHitAndFraction( (*jt), (*hitcol.find(*jt)).energy()/energy );
  
}



