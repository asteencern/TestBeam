import FWCore.ParameterSet.Config as cms

hgcaltbclusters = cms.EDProducer("HGCalTBClusterProducer",
                                 ElectronicMapFile = cms.untracked.string('HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt'),
                                 OutputCollectionName = cms.string(''),
                                 rechitCollection = cms.InputTag('hgcaltbrechits',"","unpack"),
                                 LayerZPositions= cms.untracked.vdouble(0.0, 4.67, 9.84, 14.27, 19.25, 20.4, 25.8, 31.4),#cern config 2
                                 minEnergy = cms.untracked.double(3.1e-5), #3.1e-5=~0.6 MIP=~9 ADC count
                                 RemoveSpecialCells = cms.untracked.bool( True ),
                                 PositionWeightsOption = cms.untracked.string('logarithmic'),
                                 LogWeightParams = cms.untracked.vdouble(5,1)
                                 )
