import FWCore.ParameterSet.Config as cms

haloFilterAnalyzer = cms.EDAnalyzer('HaloFilterPerformanceAnalyzer',
                                    printEventInfo = cms.bool(False),
                                    minimumMET = cms.double(50.0),     # GeV
                                    outputFile = cms.string('haloFilterPerformance.root'),
)
