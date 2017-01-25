import FWCore.ParameterSet.Config as cms

hgcaltbntuple = cms.EDAnalyzer("NtupleProducer",
                               HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                               HGCALTBCLUSTERS = cms.InputTag("hgcaltbclusters","","unpack" ),
                               CERN_8layers_config = cms.untracked.int32( 0 ),
                               minEnergy = cms.untracked.double( 3.1e-5 ), #3.1e-5 ~ 0.6 MIP
                               sensorSize = cms.untracked.int32( 128 )
                               )
