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
                                        PrepareTreeForDisplay = cms.untracked.bool( True )
                                        )
