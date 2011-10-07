# Marissa Rodenburg
# 2 October 2011

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
process = cms.Process("HALO")

# Straight from Ronny's halo_reco.py
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentCosmics_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')

process.MessageLogger.cerr.threshold = 'DEBUG'
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.GlobalTag.globaltag = 'GR_R_42_V14::All'

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:eventsNotTaggedBy2011.root')
)

# Pick number of events to run over
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
)

#############################################################
# Add CSC skims requiring events with segments
#############################################################
#process.load("DPGAnalysis/Skims/CSCSkim_cfi")
#process.cscSkim.minimumSegments = 1
#process.cscSkim.minimumHitChambers = 1

#############################################################
# Stuff to get the full halo data filled
#############################################################
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load("RecoMuon.MuonIdentification.muonIdProducerSequence_cff")
process.muonsFromCosmics = process.muontiming.clone( MuonCollection="muonsFromCosmics")
process.MessageLogger.suppressWarning = cms.untracked.vstring("CSCHaloData")

#############################################################
# Add Ronny's EDAnalyzer to record halo rates by fill
#############################################################
#process.haloanalyzer = cms.EDAnalyzer("CMSEventAnalyzer",
#                                      OutputFile = cms.string("eventsNotTaggedBy2011_HALO.root"),
#                                      BeamHaloSummaryLabel = cms.InputTag("BeamHaloSummary"),
#)

#############################################################
# Require events have zero-bias trigger
#############################################################
#process.load('L1Trigger.Skimmer.l1Filter_cfi')
#process.l1Filter.algorithms = cms.vstring('L1_ZeroBias')

#process.load(HLTrigger.HLTfilters.hltLevel1GTSeed_cfi)
#process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(False)
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('L1_ZeroBias')


#############################################################
# Tell the process what filename to use to save the output
#############################################################
process.Out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('eventsNotTaggedBy2011_HALO.root'),
#                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('skim_step')),
#                               outputCommands = cms.untracked.vstring('drop *',
#                                                                      "keep *_hltTriggerSummaryAOD_*_*",
#                                                                      "keep *_cscSegments_*_*",
#                                                                      "keep *_BeamHaloSummary_*_*",
#                                                                      "keep *_CSCHaloData_*_*",
#                                                                      "keep *_GlobalHaloData_*_*",
#                                                                      "keep *_*muon*_*_*",
#                                                                      "keep *_beamhaloTracks_*_*",
#                                                                      "keep *_*met*_*_*",
#                                                                      )
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

#############################################################
# Put all processes on the path and schedule them
#############################################################
# Path and EndPath definitions
# Use process.reconstructionCosmics for CMSSW_4_3_X
process.reconstruction_step = cms.Path(process.muonsFromCosmics)    
#process.skim_step = cms.Path(process.cscSkim)
#process.filter_step = cms.Path(process.hltLevel1GTSeed)
#process.filter_step = cms.Path(process.l1Filter)
process.beamhalo_step = cms.Path(process.CSCHaloData*process.BeamHaloId)
#process.analyze_step = cms.Path(process.haloanalyzer)
process.endjob_step = cms.EndPath(process.Out*process.endOfProcess)


# Schedule definition
process.schedule = cms.Schedule(
    process.reconstruction_step,
#    process.skim_step,
#    process.filter_step,
    process.beamhalo_step,
#    process.analyze_step,
    process.endjob_step
)
