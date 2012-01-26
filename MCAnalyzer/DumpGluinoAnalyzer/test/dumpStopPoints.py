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
    'file:/store/user/rodenm/data/mcgluino/EXO_HSCP_Gluino700_Summer11ReReproduce_GEN_SIM_exotica-EXO_HSCP_Gluino700_Summer11ReReproduce_GEN_SIM-67d0e10d7777587419655e23d92ea291_USER_HSCPgluino_M_700_7TeV_pythia6_cff_py_GEN_SIM_100_1_wPE.root')
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
