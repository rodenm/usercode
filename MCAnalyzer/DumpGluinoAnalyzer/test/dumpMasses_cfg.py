STOPMASS = 400
NEUTRALINOMASS = 320

import FWCore.ParameterSet.Config as cms

process = cms.Process('SCAN')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'file:stoppedHSCP_stage2_GEN-HLT_stop300_175_10_1_Xtm.root')
    #'file:stoppedHSCP_stage2_GEN-HLT_stop300_199_10_1_pOD.root')
    #'file:stoppedHSCP_stage2_GEN-HLT_stop400_175_10_1_j5S.root')
    #'file:stoppedHSCP_stage2_GEN-HLT_stop400_199_10_1_qon.root')
    #'file:stoppedHSCP_stage2_GEN-HLT_stop400_225_10_1_ILO.root')
    #'file:stoppedHSCP_stage2_GEN-HLT_stop400_237_10_1_8v6.root')
    #'file:stoppedHSCP_stage2_GEN-HLT_gluino300_245_10_1_apg.root')
    #'file:stoppedHSCP_stage2_GEN-HLT_gluino600_500_10_1_IQU.root')
    'file:stoppedHSCP_stage2_GEN-HLT_stop'+str(STOPMASS)+'_'+str(NEUTRALINOMASS)+'.root')
    #'file:HSCPstop_M_400_7TeV_pythia6_cff_py_GEN_SIM_100_1_eth.root')


)

process.analyze = cms.EDAnalyzer(
    'DumpGluinoAnalyzer',
    stoppedFile = cms.string('stoppedPointDump_stop'+str(STOPMASS)+'_'+str(NEUTRALINOMASS)+'.txt')
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
