#! /usr/bin/env python
#
# Marissa Rodenburg, 2 October 2011
#
# examineCSCHaloData.py
#
# Script reads in CSCHaloData objects and looks at which of the 2010 and
# 2011 handles are triggered by each event.
#
# Requires RECO-level data for input, also requires that the 2011 CSCHaloData
# has been written into the root file.


# Import everything from ROOT
from ROOT import *
from DataFormats.FWLite import Events, Handle
import sys, os

if len(sys.argv) == 4:
    inname = sys.argv[1]
    process = sys.argv[2]
    maxEvents = int(sys.argv[3])
    print 'Opening root file: ' + inname
else:
    print 'Wrong number of arguments.'
    print 'usage: python examineCSCHaloData.py <input_RECO_file.root> <process> <max events to process>'
    sys.exit(2)

events = Events(inname)
print "Begin analyzing " + str(events.size()) + " events."

strlist1 = inname.split('/')
strlist2 = strlist1[len(strlist1)-1].split('.')
basename = strlist2[0]

cscHaloDataRECOH = Handle('reco::CSCHaloData')
cscHaloDataRECOL = ('CSCHaloData','', 'RECO')

cscHaloDataHALOH = Handle('reco::CSCHaloData')
cscHaloDataHALOL = ('CSCHaloData','', process)

cscsegmentsH = Handle('edm::RangeMap<CSCDetId,edm::OwnVector<CSCSegment,edm::ClonePolicy<CSCSegment> >,edm::ClonePolicy<CSCSegment> >')
cscsegmentsL = ('cscSegments','','RECO')

halosummaryRECOH = Handle('reco::BeamHaloSummary')
halosummaryRECOL = ('BeamHaloSummary','','RECO')

halosummaryHALOH = Handle('reco::BeamHaloSummary')
halosummaryHALOL = ('BeamHaloSummary','', process)

i = 1
eventsForFireworks = []
for event in events:

    event.getByLabel(halosummaryHALOL, halosummaryHALOH)
    halosummaryHALO = halosummaryHALOH.product()

    event.getByLabel(halosummaryRECOL, halosummaryRECOH)
    halosummaryRECO = halosummaryRECOH.product()

    # Looking for events that are tagged as halo in 2010, but not in 2011
    if halosummaryRECO.CSCLooseHaloId() and not halosummaryHALO.CSCLooseHaloId():
        event.getByLabel(cscHaloDataRECOL, cscHaloDataRECOH)
        cscHaloDataRECO = cscHaloDataRECOH.product()

        event.getByLabel(cscHaloDataHALOL, cscHaloDataHALOH)
        cscHaloDataHALO = cscHaloDataHALOH.product()
        
        event.getByLabel(cscsegmentsL, cscsegmentsH)
        cscSegments = cscsegmentsH.product()
        
        eventAux = event.eventAuxiliary()
        id = eventAux.id()
        run = id.run()
        lumi = id.luminosityBlock()
        eventnum = id.event()
        eventInfo = str(run) + ':' + str(lumi) + ':' + str(eventnum)
        eventsForFireworks.append(eventInfo + '\n')
        
        triggers = cscHaloDataRECO.NumberOfHaloTriggers()
        tracks   = cscHaloDataRECO.NumberOfHaloTracks()
        oottriggers = cscHaloDataRECO.NumberOfOutOfTimeTriggers()
        
        print '---------------------------------'+ str(i) + ' ' + eventInfo + '-----------------------------'
        print 'Event has ' + str(cscSegments.size()) + ' CSCSegments.'
        print 'NumberOfHaloTriggers(): ' + str(triggers)
        print 'NumberOfHaloTracks(): ' + str(tracks)
        print 'NumberOfOutOfTimeTriggers(): ' + str(oottriggers)
        print ''
        print 'NumberOfHaloTriggers(): ' + str(cscHaloDataHALO.NumberOfHaloTriggers())
        print 'NumberOfHaloTracks(): ' + str(cscHaloDataHALO.NumberOfHaloTracks())
        print 'NOutOfTimeHits() > 10 && NFlatHaloSegments() > 2: ' + str(cscHaloDataHALO.NOutOfTimeHits()) + ', ' + str(cscHaloDataHALO.NFlatHaloSegments())
        print 'GetSegmentsInBothEndcaps(): ' + str(cscHaloDataHALO.GetSegmentsInBothEndcaps())
        print 'NTracksSmalldT(): ' + str(cscHaloDataHALO.NTracksSmalldT())
        print ''
        
        if maxEvents != -1 and i >= maxEvents: break
        i = i + 1

print 'Creating file: ' + basename + '_pickEvents.txt'
file = open(basename + '_pickEvents.txt', 'w')
file.writelines(eventsForFireworks)
file.close()
