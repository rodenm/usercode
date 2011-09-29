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
    
output = TFile(parts[0] + '.root', 'RECREATE')

# Distance in cm
hxy     = TH2F('hxy',     'x-y distribution',                          1600,-1600,1600,1600,-1600,1600)
hrz     = TH2F('hrz',     'r-z distribution',                          3000,-3000,3000,800,0,1600)
hrz2    = TH2F('hrz2',    'r-z distribuion, preserving sign(y)',       3000,-3000,3000,1600,-1600,1600)
hr2d    = TH1F('hr2d',    '|r|',                                       800,0,1600)
hr3d    = TH1F('hr3d',    '|sqrt(x^2+y^2+z^2)|',                       1700,0,3400)
hsmallr = TH1F('hsmallr', 'z for points with r < 10.0cm',              3000,-3000,3000)
hre     = TH2F('hre',     'Energy by r',                               400,0,800,300,700,1000)
hrpt    = TH2F('hrpt',    'pT by r',                                   400,0,800,300,0,600)
hrde    = TH2F('hrde',    'Difference in E_gluino and E_rHadron by r', 400,0,800,500,0,1000)

count = 0
maxx = 0
maxy = 0
maxz = 0
maxr = 0
badr = 0
rBeforeECAL = 0
goodlines = 0

print '*** Processing only lines with gluino energies ***'

for line in input:
    words = line.split(' ')
    # Some stopped points don't have associated gluino energies because the
    # gluino's phi value with > pi/2 away from the stopped point phi.
    # For these points, plot nothing.
    if len(words) != 15:
        #print "Bad line: " + line
        count = count + 1
        continue
    else: goodlines = goodlines + 1
    x = float(words[1])/10.0
    y = float(words[2])/10.0
    z = float(words[3])/10.0
    r = sqrt(x*x + y*y)
    
    if math.fabs(x) > maxx: maxx = x
    if math.fabs(y) > maxy: maxy = y
    if math.fabs(z) > maxz: maxz = z
    if r > maxr: maxr = r
    if math.fabs(r) < 10.0:
        badr = badr + 1
        #print "r, z = " + str(r) + ", " + str(z)
        hsmallr.Fill(z)
    if math.fabs(r) < 131.0:
        rBeforeECAL = rBeforeECAL + 1
        
    hxy.Fill(x, y)
    hrz.Fill(z, r)
        
    if y > 0: hrz2.Fill(z, r)
    else: hrz2.Fill(z, r*-1.0)
    
    hr2d.Fill(r)
    hr3d.Fill(sqrt(x*x + y*y + z*z))

    hre.Fill(r, float(words[9]))
    hrpt.Fill(r, float(words[8]))

    # Plot difference in energy between gluino and r-hadron by r
    hrde.Fill(r, fabs(float(words[9]) - float(words[14])))
        
    count = count + 1
    #if count > 10: break

hxy.SetMarkerSize(0.3)
hrz.SetMarkerSize(0.3)
hrz2.SetMarkerSize(0.3)
hr3d.SetMarkerSize(0.3)
hre.SetMarkerSize(0.3)
hrpt.SetMarkerSize(0.3)
hrde.SetMarkerSize(0.3)
#hrz.SetMarkerStyle(7)

output.Write()
input.close()
print 'Total: ' + str(count) + ' lines.'
print 'Max (x, y, z, r) : (' + str(maxx) + ', ' + str(maxy) + ', ' + str(maxz) + ', ' + str(maxr) + ')'
print 'Points with r < 10.0 cm: ' + str(badr) + '(' + str(100.0*badr/goodlines) + '%)'
print 'Points with r < 131.0 cm (stopped before ECAL): ' + str(rBeforeECAL) + '(' + str(100.0*rBeforeECAL/goodlines) + '%)'
print 'Lines correctly written:  ' + str(goodlines) 
print 'Finished.'


