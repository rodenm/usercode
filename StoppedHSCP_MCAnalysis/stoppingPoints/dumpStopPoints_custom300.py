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
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_10_1_e2I.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_11_1_vmG.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_12_1_3N1.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_13_1_W3i.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_14_1_QSq.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_15_1_RF8.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_16_1_sMz.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_17_1_WZX.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_18_1_xfk.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_19_1_lcW.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_1_1_jw5.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_20_1_bx0.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_2_1_m7I.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_3_1_Scd.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_4_1_WCl.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_5_1_tqd.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_6_1_VfH.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_7_1_c48.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_8_1_RdK.root',
    'file:~rodenm/work/stop_gluino/CMSSW_4_1_5/src/production/gluino300_stage1_custom/res/hscpGluino_M300_7TeV-pythia6_Summer11-START311_V2-v1_GEN-SIM_stage1_9_1_yRi.root',
    )
)

process.analyze = cms.EDAnalyzer('DumpGluinoAnalyzer',
                                 stoppedFile = cms.string("stoppedPointDump_gluino300_custom.txt")
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
