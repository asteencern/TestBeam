import FWCore.ParameterSet.Config as cms

from HGCal.Reco.hgcaltbrechitproducer_cfi import *
from HGCal.Reco.hgcaltbclusterproducer_cfi import *
from HGCal.Reco.hgcaltbcalotrackproducer_cfi import *
from HGCal.Reco.hgcaltbntupleproducer_cfi import *

LocalRecoSeq  = cms.Sequence(hgcaltbrechits)
ClusterRecoSeq  = cms.Sequence(hgcaltbclusters)
ClusterRecoSeq  = cms.Sequence(hgcaltbcalotracks)
NtupleRecoSeq  = cms.Sequence(hgcaltbntuple)

