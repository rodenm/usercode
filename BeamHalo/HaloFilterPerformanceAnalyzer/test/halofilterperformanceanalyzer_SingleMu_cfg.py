import FWCore.ParameterSet.Config as cms

process = cms.Process("PERF")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration/StandardSequences/MagneticField_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')

process.MessageLogger.cerr.threshold = 'WARNING'
process.MessageLogger.categories.append('HaloFilterPerformanceAnalyzer')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
#process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.load("Configuration/StandardSequences/Magneticfield_cff")
#process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
#process.load("Configuration/StandardSequences/RawToDigi_cff")
#process.load('Configuration.EventContent.EventContent_cff')
#process.load("RecoMET/Configuration/RecoMET_BeamHaloId_cff")
#process.load("Configuration/StandardSequences/ReconstructionCosmics_cff")
#process.load("RecoMuon/Configuration/RecoMuon_cff")
#process.load("RecoMuon/Configuration/RecoMuonCosmics_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )
process.GlobalTag.globaltag ='GR_P_V32::All'
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/store/data/Run2012B/SingleMu/RECO/PromptReco-v1/000/194/076/3AEEBC32-AC9E-E111-BE36-485B3962633D.root' # 194076
    )
)

process.load("BeamHalo.HaloFilterPerformanceAnalyzer.halofilterperformanceanalyzer_cfi")
process.haloFilterAnalyzer.outputFile = cms.string('haloFilterPerformance_SingleMu.root')

process.p = cms.Path(process.haloFilterAnalyzer)
                                            
