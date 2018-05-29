import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import os,sys

options = VarParsing.VarParsing('standard') # avoid the options: maxEvents, files, secondaryFiles, output, secondaryOutput because they are already defined in 'standard'
#Change the data folder appropriately to where you wish to access the files from:


options.register('inputFiles',
                [''],
                 VarParsing.VarParsing.multiplicity.list,
                 VarParsing.VarParsing.varType.string,
                 'Paths to the input files.'
                )

options.register('processedFile',
                 '/eos/user/t/tquast/outputs/Testbeam/July2017/rechits/RECHITS_1303.root',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where pedestal histograms are stored')

options.register('ntupleFile',
                 '',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file where ntuples are stored')

options.register('electronicMap',
                 'map_CERN_Hexaboard_July_6Layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the electronic map file in HGCal/CondObjects/data/')

options.register('beamEnergy',
                250,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Beam energy.'
                )

options.register('beamParticlePDGID',
                211,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Beam particles PDG ID.'
                )

options.register('setupConfiguration',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'setupConfiguration.'
                )

options.register('layerPositionFile',
                 'layer_distances_CERN_Hexaboard_July_6Layers.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'File indicating the layer positions in mm.')

options.register('hgcalLayout',
                 'layerGeom_desymarch2018_configuration4.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Name of the hgcal layout file in HGCal/CondObjects/data/')

options.register('physicsListUsed',
                "",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Specify the used physics list to be passed forward to the run data object.'
                )

options.register('areaSpecification',
                "H2",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Area which was used (for DWC simulation).'
                )

options.register('NHexaBoards',
                3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of hexaboards for analysis.'
                )

options.register('noisyChannelsFile',
                 '/home/tquast/tb2017/pedestals/noisyChannels_1190.txt',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Channels which are noisy and excluded from the reconstruction')

options.register('MaskNoisyChannels',
                1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Ignore noisy channels in the reconstruction.'
                )

options.register('reportEvery',
                10000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Path to the file from which the DWCs are read.'
                )

options.register('stopAtEvent',
                10000,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Stop processing after this event.'
                )

#options.register('outputFile',
#                 '/home/tquast/tb2017/analysis/energyReco_simTest.root',
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.string,
#                 'Output folder where analysis output are stored')



options.parseArguments()
print options



electronicMap="HGCal/CondObjects/data/%s" % options.electronicMap
hgcalLayout="HGCal/CondObjects/data/%s" % options.hgcalLayout

################################
process = cms.Process("gensim")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.stopAtEvent)
)

################################
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.ntupleFile))
####################################
# Reduces the frequency of event count couts
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
####################################


process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(options.processedFile),                                  
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_*_HGCALTBRECHITS_*',
                                                                         'keep *_*_HGCalTBDATURATracks_*',
                                                                         'keep *_*_FullRunData_*')
)

process.source = cms.Source("HGCalTBGenSimSource",
                        fileNames=cms.untracked.vstring(["file:%s" % file for file in options.inputFiles]),
                        RechitOutputCollectionName = cms.string('HGCALTBRECHITS'), 
                        produceDATURATracksInsteadOfDWCs = cms.untracked.bool(True),
                        DWCOutputCollectionName = cms.string(''), 
                        DATURAOutputCollectionName = cms.string('HGCalTBDATURATracks'), 
                        RunDataOutputCollectionName = cms.string('FullRunData'), 
                        e_mapFile_CERN = cms.untracked.string(electronicMap),
                        layerPositionFile=cms.string(options.layerPositionFile),
                        MaskNoisyChannels=cms.untracked.bool(bool(options.MaskNoisyChannels)),
                        ChannelsToMaskFileName=cms.untracked.string(options.noisyChannelsFile),
                        beamEnergy=cms.untracked.double(options.beamEnergy),
                        beamParticlePDGID=cms.untracked.int32(options.beamParticlePDGID),                        
                        energyNoise=cms.untracked.double(0.0),  #indicated in MIPs
                        setupConfiguration=cms.untracked.uint32(options.setupConfiguration),
                        energyNoiseResolution=cms.untracked.double(1./6.), #indicated in MIPs
                        GeVToMip=cms.untracked.double(1./(84.9*pow(10.,-6))),   #apply an overall scaling of the recorded intensities in the cells
                        areaSpecification = cms.untracked.string(options.areaSpecification),
                        physicsListUsed = cms.untracked.string(options.physicsListUsed),
                        datura_resolutions = cms.untracked.vdouble(6*[0.0184])        #set to the expected resolutions according to the manual
                        )


process.rechitntupler = cms.EDAnalyzer("RecHitNtupler",
                                       InputCollection=cms.InputTag("source","HGCALTBRECHITS"),
                                       RUNDATA = cms.InputTag("source", "FullRunData" ),
                                       ElectronicMap=cms.untracked.string(electronicMap),
                                       layerPositionFile = cms.untracked.string(options.layerPositionFile),
                                       DetectorLayout=cms.untracked.string(hgcalLayout),
                                       SensorSize=cms.untracked.int32(128),
                                       EventPlotter=cms.untracked.bool(True),
                                       MipThreshold=cms.untracked.double(5.0),
                                       NoiseThreshold=cms.untracked.double(0.5)
)

process.trackimpactntupler = cms.EDAnalyzer("ImpactPointNtupler",
                                       DATURATelescopeData = cms.InputTag("source","HGCalTBDATURATracks" ),
                                       RUNDATA = cms.InputTag("source", "FullRunData" ),
                                       nLayers=cms.untracked.int32(options.NHexaBoards)
)


####################################

process.p = cms.Path( process.rechitntupler * process.trackimpactntupler )


process.end = cms.EndPath(process.output)
