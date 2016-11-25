import FWCore.ParameterSet.Config as cms

hgcaltbntuple = cms.EDAnalyzer("NtupleProducer",
                                       HGCALTBRECHITS = cms.InputTag("hgcaltbrechits","","unpack" ),
                                       HGCALTBCLUSTERS = cms.InputTag("hgcaltbclusters","","unpack" ),
                                       CERN_8layers_config = cms.untracked.int32( 0 ),
                                       minEnergy = cms.untracked.double( 30.0 ),
                                       sensorSize = cms.untracked.int32( 128 ),
                                       skirocADCToMip= cms.untracked.vdouble(16.95,16.6933,16.0208,17.0226,17.6833,17.1882,16.4708,15.9629,17.1542,16.7324,16.5,17.5457,15.3652,16.273,16.4111,15.2706)
                                       )
