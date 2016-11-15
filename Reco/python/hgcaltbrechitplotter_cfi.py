import FWCore.ParameterSet.Config as cms

hgcaltbrechitsplotter = cms.EDAnalyzer("RecHitPlotter",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )

hgcaltbrechitsplotter_highgain_new = cms.EDAnalyzer("RecHitPlotter_HighGain_New",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
               HGCALTBTRACKS  = cms.InputTag("hgcaltbtracks","","unpack" )
                              )

hgcaltbeventdisplay = cms.EDAnalyzer("EventDisplay",
                                     HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                                     Nlayers = cms.untracked.int32( 8 ),
                                     SensorSize = cms.untracked.int32( 128 ),
                                     minEnergy = cms.untracked.double( 50 ),
                                     CMThreshold = cms.untracked.int32( 30 )
                                     )

hgcaltbrechitsplotter_highgain_cluster = cms.EDAnalyzer("RecHitPlotter_HighGain_Cluster",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )

hgcaltbrechitsplotter_highgain_cluster_telescope = cms.EDAnalyzer("RecHitPlotter_HighGain_Cluster_Telescope",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )

hgcaltbsimrechitsplotter_highgain_new = cms.EDAnalyzer("SimRecHitPlotter_HighGain_New",
               HGCALTBRECHITS = cms.InputTag("hgcaltbsimrechits","","unpack" )
                              )

hgcaltbrechitsplotter_highgain_correlation = cms.EDAnalyzer("RecHitPlotter_HighGain_Correlation",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )
hgcaltbrechitsplotter_highgain_correlation_cm = cms.EDAnalyzer("RecHitPlotter_HighGain_Correlation_CM",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )

FourLayerRecHitPlotterMax = cms.EDAnalyzer("FourLayerRecHitPlotterMax",
               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" )
                              )

LayerSumAnalyzer = cms.EDAnalyzer("Layer_Sum_Analyzer",
                                  HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                                  layers_config = cms.int32(1),
                                  mapFile_CERN = cms.string('HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt'),
                                  mapFile_FNAL = cms.string('')
                              )

hgcalclusteringexample = cms.EDAnalyzer("ClusteringExampleAnalyzer",
                                        HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                                        Nlayers = cms.untracked.int32( 8 )
                                        )
