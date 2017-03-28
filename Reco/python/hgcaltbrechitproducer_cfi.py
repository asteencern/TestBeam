
import FWCore.ParameterSet.Config as cms

hgcaltbrechits = cms.EDProducer("HGCalTBRecHitProducer",
                                OutputCollectionName = cms.string(''),
                                digiCollection = cms.InputTag('hgcaltbdigis'),
                                pedestalLow = cms.string('CondObjects/data/Ped_LowGain_L8.txt'),
                                pedestalHigh = cms.string('CondObjects/data/Ped_HighGain_L8.txt'),
                                gainLow = cms.string(''),
                                gainHigh = cms.string(''),
                                LG2HG_CERN = cms.vdouble(10.2, 10.2, 10., 10., 9.8, 9.8, 9.7, 9.7, 9.7, 9.7, 9.8, 9.8, 9.8, 9.8, 9.2, 9.2),
                                LG2HG_FNAL = cms.vdouble(1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.),
                                adcSaturation = cms.int32(1800),
                                mapFile_CERN = cms.string('HGCal/CondObjects/data/map_CERN_8Layers_Sept2016.txt'),
                                mapFile_FNAL = cms.string(''),
                                layers_config = cms.int32(-1),
                                CommonModeThreshold = cms.untracked.double(2), # in MIP
                                doCommonMode = cms.untracked.bool(True),
                                MPVToMIP = cms.untracked.double(0.94),
                                MIPToMeV = cms.untracked.double(51.9e-06),
                                ConvertEnergyToGeV = cms.untracked.bool(True),
                                skirocADCToMip = cms.untracked.vdouble(17.0265, 16.7563, 16.0833, 17.125, 17.7913, 17.2875, 16.55, 15.9688, 17.2652, 16.85, 16.5979, 17.6152, 15.4814, 16.4, 16.5356, 15.3364)
                                )
