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
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/HSCPgluino_M-200_415_GEN-SIM.root')
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
