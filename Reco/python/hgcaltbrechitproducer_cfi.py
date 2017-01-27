
import FWCore.ParameterSet.Config as cms

hgcaltbrechits = cms.EDProducer("HGCalTBRecHitProducer",
                                OutputCollectionName = cms.string(''),
                                digiCollection = cms.InputTag('hgcaltbdigis'),
                                pedestalLow = cms.string('CondObjects/data/Ped_LowGain_L8.txt'),
                                pedestalHigh = cms.string('CondObjects/data/Ped_HighGain_L8.txt'),
                                gainLow = cms.string(''),
                                gainHigh = cms.string(''),
                                LG2HG_CERN = cms.vdouble(10.2, 10.2, 10., 10., 9.8, 8.8, 9.7, 9.7, 9.7, 9.7, 9.8, 9.8, 9.8, 9.8, 9.2, 9.2),
                                LG2HG_FNAL = cms.vdouble(1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.),
                                adcSaturation = cms.int32(1800),
                                mapFile_CERN = cms.string('HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt'),
                                mapFile_FNAL = cms.string(''),
                                layers_config = cms.int32(-1),
                                CommonModeThreshold = cms.untracked.double(2), # in MIP
                                doCommonMode = cms.untracked.bool(True),
                                MPVToMIP = cms.untracked.double(0.94),
                                MIPToMeV = cms.untracked.double(52.81e-06),
                                ConvertEnergyToGeV = cms.untracked.bool(True),
                                skirocADCToMip = cms.untracked.vdouble(15.926 ,15.6252 ,14.8598 ,15.9285 ,16.873 ,15.6593 ,15.6549 ,14.5163 ,16.2659 ,15.5812 ,15.369 ,16.4178 ,14.4392 ,15.2915 ,15.51 ,14.2051)
                                )
