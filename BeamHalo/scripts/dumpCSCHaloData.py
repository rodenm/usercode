import sys
from ROOT import *

# Import what we need from FWLite
from DataFormats.FWLite import Events, Handle

#######################################################
# Handle commandline inputs & file io
#######################################################
if len(sys.argv) == 5:
    event_run = int(sys.argv[1])
    event_lb = int(sys.argv[2])
    event_id = int(sys.argv[3])
    inname = sys.argv[4]
    print 'Opening root file: ' + inname
else:
    print 'Wrong number of arguments.'
    print 'usage: python dumpCSCHaloData.py <run> <lb> <event> <RECO file>'
    sys.exit(2)
    
######################################################
# Handles and labels for recovering RECO variables
######################################################
halosummaryH = Handle('reco::BeamHaloSummary')
halosummaryL = ('BeamHaloSummary','','RECO')

cschalodataH = Handle('reco::CSCHaloData')
cschalodataL = ('CSCHaloData','','RECO')

######################################################
# Loop over events
######################################################
events = Events(inname)
nEvents = events.size()

for event in events:

    # Search for the input event
    eventAux = event.eventAuxiliary()
    id = eventAux.id()
    run = id.run()
    lumi = id.luminosityBlock()
    eventnum = id.event()

    if not (event_run == run and event_lb == lumi and event_id == eventnum):
        continue

    # Correct event: print out CSCHaloData
    event.getByLabel(halosummaryL, halosummaryH)
    halosummary = halosummaryH.product()

    event.getByLabel(cschalodataL, cschalodataH)
    cscdata = cschalodataH.product()

    print ''
    print str(run) + ':' + str(lumi) + ':' + str(eventnum)
    print ' --------------------------- '
    print 'CSCLoose: ' + str(halosummary.CSCLooseHaloId())
    print 'CSCTight: ' + str(halosummary.CSCTightHaloId())
    print ' --------------------------- '
    print 'NumberOfHaloTriggers():    ' + str(cscdata.NumberOfHaloTriggers())
    print 'NumberOfHaloTracks():      ' + str(cscdata.NumberOfHaloTracks())
    print 'NOutOfTimeHits() > 10:     ' + str(cscdata.NOutOfTimeHits())
    print 'SegmentsInBothEndcaps():   ' + str(cscdata.GetSegmentsInBothEndcaps())
    print 'NTracksSmalldT():          ' + str(cscdata.NTracksSmalldT())
    print 'NFlatHaloSegments() > 2/3: ' + str(cscdata.NFlatHaloSegments())
    print ''
    
