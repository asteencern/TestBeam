import FWCore.ParameterSet.Config as cms

hgcaltbshower = cms.EDAnalyzer("ShowerAnalyzer",
                               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                               Nlayers = cms.untracked.int32( 8 ),
                               NSkirocsPerLayer = cms.untracked.int32( 2 ),
                               SensorSize = cms.untracked.int32( 128 ),
                               CMThreshold = cms.untracked.int32( 30 ),
                               CERN_8layers_config = cms.untracked.int32( 0 ),
                               minEnergy = cms.untracked.double( 100.0 )
                               )
