# Marissa Rodenburg, 28 July 2011
#
# Takes stoppedPoint.txt files as input.
# Reads in points, plots them on x-y and r-z 2D histograms.
#

import sys, os
from ROOT import *

print sys.argv

if len(sys.argv) == 2:
    inname = sys.argv[1]
    print 'Opening stopped points file: ' + inname
else:
    print 'Wrong number of arguments.'
    print 'usage: python plot_stopPoints.py <stoppedPoints filename>'
    sys.exit(2)

input = open(inname, 'r')

parts = os.path.basename(inname).split('.')
if len(parts) > 2:
    print 'Bad input filename: too many \'.\''
    sys.exit(1)

maxx = 0.
maxy = 0.
maxz = 0.
maxr = 0.
count = 0
cavern_count = 0
stop_count = 0    # include only points inside detector
tracker_count = 0
for line in input:
    words = line.split(' ')
    if len(words) != 4:
        print "Bad line: " + line
        count = count + 1
        continue
    # Want x, y, z, and r in cm (given in mm)
    x = float(words[1])/10.0
    y = float(words[2])/10.0
    z = float(words[3])/10.0
    r = sqrt(x*x + y*y)

    if math.fabs(x) > maxx: maxx = x
    if math.fabs(y) > maxy: maxy = y
    if math.fabs(z) > maxz: maxz = z
    if r > maxr: maxr = r
    if math.fabs(r) < 131.0 and math.fabs(z) < 300:
        tracker_count = tracker_count + 1

    if r > 728.5 or fabs(z) > fabs(1080): cavern_count = cavern_count + 1
    else: stop_count = stop_count + 1
    count = count + 1
    #if count > 10: break

#output.Write()
input.close()
print 'Total: ' + str(count) + ' events.'
print 'Max (x, y, z, r) : (' + str(maxx) + ', ' + str(maxy) + ', ' + str(maxz) + ', ' + str(maxr) + ')'
print 'Tracker_count = ' + str(tracker_count)
print 'Cavern_count = ' + str(cavern_count)
print 'Detector_count = '  + str(stop_count)


