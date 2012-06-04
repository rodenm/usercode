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
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.GlobalTag.globaltag ='GR_P_V32::All'
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'/store/data/Run2012A/MinimumBias/RECO/PromptReco-v1/000/193/339/66B54DF6-4D97-E111-9B8F-BCAEC5364C42.root' # 193339
    #
    # If WBM is correct, the HLT_ZeroBias trigger actually fired for this run...
    #
    '/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/F6BA4E81-849E-E111-923D-5404A63886C0.root', # 194076
    '/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/F6AD6EC4-889E-E111-90BD-5404A63886B9.root',
    '/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/E08A29C8-949E-E111-8A7F-003048F118AA.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/C03AE3B2-8D9E-E111-B3D8-003048D2C020.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/C0269A52-879E-E111-A5BD-5404A63886D4.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/BE5B5947-919E-E111-8C81-003048F024C2.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/B6F6EA3F-929E-E111-8C80-BCAEC5329703.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/AC8CDDD7-8A9E-E111-B004-BCAEC518FF6B.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/A20FE080-849E-E111-83E4-E0CB4E553673.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/9E540330-8A9E-E111-A1EC-5404A640A639.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/8ABBA708-889E-E111-B1BA-001D09F24D8A.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/88C022D1-AA9E-E111-A323-0025901D5D7E.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/7CD4C281-849E-E111-89EA-003048D2BE08.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/5C19110B-889E-E111-8F28-001D09F28D54.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/5AE091F7-969E-E111-B519-E0CB4E4408E7.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/407E1855-939E-E111-AEE2-002481E0D524.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/3A26052F-859E-E111-81C8-001D09F2441B.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/34BAC215-A89E-E111-8F85-0025B32445E0.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/2CDFC714-AF9E-E111-ABAB-0030486733B4.root',
    #'/store/data/Run2012B/MinimumBias/RECO/PromptReco-v1/000/194/076/1287ABB5-8D9E-E111-A812-0030486780AC.root',
    )
)

#
# Require ZeroBias trigger to get an unbiased calculation of the halo rate
#
process.load('HLTrigger.HLTfilters.hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(False)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('L1_ZeroBias')

#
# Try filtering for ZeroBias using HLT instead (the L1 filter gave nothing)
process.zerobiasfilter = cms.EDFilter("HLTHighLevel",
                                      TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                      HLTPaths = cms.vstring("HLT_ZeroBias*"),
                                      eventSetupPathsKey = cms.string(""),
                                      andOr = cms.bool(True),
                                      throw = cms.bool(False)
                                      )

process.load("BeamHalo.HaloFilterPerformanceAnalyzer.halofilterperformanceanalyzer_cfi")
process.haloFilterAnalyzer.outputFile = cms.string('haloFilterPerformance.root')

process.p = cms.Path(process.zerobiasfilter*process.haloFilterAnalyzer)
