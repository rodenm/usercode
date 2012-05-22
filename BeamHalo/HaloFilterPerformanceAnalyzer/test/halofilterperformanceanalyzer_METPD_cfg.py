import FWCore.ParameterSet.Config as cms

process = cms.Process("PERF")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration/StandardSequences/MagneticField_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')

process.MessageLogger.cerr.threshold = 'WARNING'
process.MessageLogger.categories.append('HaloFilterPerformanceAnalyzer')
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#process.load("RecoMET/Configuration/RecoMET_BeamHaloId_cff")
#process.load("Configuration/StandardSequences/ReconstructionCosmics_cff")
#process.load("RecoMuon/Configuration/RecoMuon_cff")
#process.load("RecoMuon/Configuration/RecoMuonCosmics_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )
process.GlobalTag.globaltag ='GR_P_V32::All'
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/data/Run2012A/MinimumBias/RECO/PromptReco-v1/000/193/339/66B54DF6-4D97-E111-9B8F-BCAEC5364C42.root' # Run 193339?
    '/store/data/Run2012A/MET/RECO/PromptReco-v1/000/193/336/AAD29C1C-7797-E111-A9A7-001D09F24664.root' # Run 193336
    )
)

process.load("BeamHalo.HaloFilterPerformanceAnalyzer.halofilterperformanceanalyzer_cfi")
#process.haloFilterAnalyzer = cms.EDAnalyzer('HaloFilterPerformanceAnalyzer')
#process.haloFilterAnalyzer.outputFile = cms.string('haloFilterPerformance_METPD_pfmet.root')
process.haloFilterAnalyzer.outputFile = cms.string('haloFilterPerformance.root')
                                            

process.p = cms.Path(process.haloFilterAnalyzer)
