import FWCore.ParameterSet.Config as cms

hgcaldigis = cms.EDProducer("HGCalTBRawToDigi",
                            InputLabel=cms.InputTag("source"), ###
                            fedId=cms.untracked.int32(1000), ## FED ID is hard coded in the HGCalTBTextSource TO BE FIXED
                            
                            electronicsMap=cms.untracked.string("HGCal/CondObjects/data/map_FNAL.txt"), ### electronicsMap written by hand
                            
)
