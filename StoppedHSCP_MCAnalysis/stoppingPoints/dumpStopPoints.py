import FWCore.ParameterSet.Config as cms

process = cms.Process('SCAN')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'file:/uscms/home/rodenm/nobackup/mc_gluino/CMSSW_4_2_3_patch3/src/stopPoints/stoppedPoint_gluino1100.txt')
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/hscpGluino_M1100_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1.root')
    #'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/gluino1100/res/hscpGluino_M1100_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_1_1_0Zc.root',
    #'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/gluino1100/res/hscpGluino_M1100_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_2_1_xlP.root',
    #'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/gluino1100/res/hscpGluino_M1100_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_3_1_sX5.root',
    #'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/gluino1100/res/hscpGluino_M1100_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_4_1_GZK.root',
    #'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/gluino1100/res/hscpGluino_M1100_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_5_1_V2j.root',
    #'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/gluino1100/res/hscpGluino_M1100_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_6_1_XkP.root',
    #'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/gluino1100/res/hscpGluino_M1100_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_7_1_KrV.root',
    #'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/gluino1100/res/hscpGluino_M1100_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_8_1_L00.root',
    #'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/gluino1100/res/hscpGluino_M1100_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_9_1_L39.root',
    #'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/gluino1100/res/hscpGluino_M1100_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_10_1_Rd0.root',)
)

process.analyze = cms.EDAnalyzer('DumpGluinoAnalyzer',
                                 stoppedFile = cms.string("stoppedPointDump_gluino1100_custom.txt")
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
