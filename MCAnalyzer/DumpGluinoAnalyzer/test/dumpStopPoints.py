import FWCore.ParameterSet.Config as cms

process = cms.Process('SCAN')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:../../../STEP2_RAW2DIGI_L1Reco_RECO_PU_10_1_W51.root')
)

process.analyze = cms.EDAnalyzer('DumpGluinoAnalyzer',
                                 stoppedFile = cms.string("stoppedPointDump.txt")
)

process.p = cms.Path(process.analyze)
process.schedule = cms.Schedule(process.p)

#process.rhStopDump = cms.EDAnalyzer (
#    "RHStopDump",
#    stoppedFile = cms.string("src/stoppedPoint.txt")
#    )
#
#process.rhStopDumpstep = cms.Path (process. rhStopDump)
#process.shadule = cms.Schedule(process.rhStopDumpstep)
