import FWCore.ParameterSet.Config as cms

hgcaltbtrackingexample = cms.EDAnalyzer("TrackingExampleAnalyzer",
                                        HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                                        Nlayers = cms.untracked.int32( 8 ),
                                        #NSkirocsPerLayer = cms.untracked.int32( 2 ),
                                        SensorSize = cms.untracked.int32( 128 ),
                                        minMip = cms.untracked.int32( 10 ),
                                        maxMip = cms.untracked.int32( 48 ),
                                        CMThreshold = cms.untracked.int32( 30 ),
                                        CERN_8layers_config = cms.untracked.int32( 0 ),
                                        doTrackCleaning = cms.untracked.bool(True),
                                        #maxDistanceToRecoTrack = cms.untracked.double(13),
                                        #maxChi2 = cms.untracked.double(9.48),
                                        PrepareTreeForDisplay = cms.untracked.bool( True )
                                        )

hgcaltbtrackanalyzer = cms.EDAnalyzer("TrackAnalyzer",
                                      HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                                      HGCALTBCALOTRACKS = cms.InputTag("hgcaltbcalotracks","","unpack" ),
                                      Nlayers = cms.untracked.int32( 8 ),
                                      NSkirocsPerLayer = cms.untracked.int32( 2 ),
                                      NChannelsPerSkiroc = cms.untracked.int32( 64 ),
                                      maxChi2 = cms.untracked.double(5.0),
                                      noiseEnergyThreshold = cms.untracked.double( 9.0 ),
                                      LayerZPositions = cms.untracked.vdouble(0.0, 4.67, 9.84, 14.27, 19.25, 20.4, 25.8, 31.4),#cern config 2
                                      #LayerZPositions = cms.untracked.vdouble(0.0, 5.35, 10.52, 14.44, 18.52, 19.67, 23.78, 25.92),#cern config 1
                                      SensorSize = cms.untracked.int32(128)
                                      )

hgcaltbcellefficiency = cms.EDAnalyzer("HGCalCellEfficiencyAnalyzer",
                                       HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                                       HGCALTBCLUSTERS = cms.InputTag("hgcaltbclusters","","unpack" ),
                                       Nlayers = cms.untracked.int32(8),
                                       NSkirocsPerLayer = cms.untracked.int32(2),
                                       NChannelsPerSkiroc = cms.untracked.int32(64),
                                       SensorSize = cms.untracked.int32(128),
                                       minTouchedLayers = cms.untracked.int32(4),
                                       maxChi2 = cms.untracked.double(9.48),
                                       minEnergyEfficiency = cms.untracked.double(3.1e-5), #0.6 mip
                                       clusterMaxEnergyThresholdForTrack = cms.untracked.double(2.1e-4), #~4 mio
                                       maxDistanceToRecoTrack = cms.untracked.double(2.0),
                                       LayerZPositions = cms.untracked.vdouble(0.0, 4.67, 9.84, 14.27, 19.25, 20.4, 25.8, 31.4)#cern config 2
                                       )
