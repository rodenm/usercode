# Marissa Rodenburg, 28 July 2011
#
# Takes stoppedPoint.txt files as input.
# Reads in points, plots them on x-y and r-z 2D histograms.
#

import sys, os
from ROOT import *

print sys.argv

if len(sys.argv) == 3:
    inname = sys.argv[1]
    inname2 = sys.argv[2]
    print 'Opening stopped points files: ' + inname + ', ' + inname2
else:
    print 'Wrong number of arguments.'
    print 'usage: python plot_stopPoints.py <stoppedPoints filename>'
    sys.exit(2)

input1 = open(inname, 'r')
input2 = open(inname2, 'r')

parts = os.path.basename(inname).split('.')
if len(parts) > 2:
    print 'Bad input filename: too many \'.\''
    sys.exit(1)
    
output = TFile(parts[0] + '_overlay.root', 'RECREATE')


# Handle input 1
# Distance in cm
hxy = TH2F('hxy','',320,-1600,1600,320,-1600,1600)
hrz = TH2F('hrz','',600,-3000,3000,160,0,1600)
hrz2 = TH2F('hrz2','',600,-3000,3000,320,-1600,1600)
hr2d = TH1F('hr2d','',800,0,1600)
hr3d = TH1F('hr3d','',1700,0,3400)
hsmallr = TH1F('hsmallr','z for points with r < 10.0cm',3000,-3000,3000)

count = 0
maxx = 0
maxy = 0
maxz = 0
maxr = 0
badr = 0
for line in input1:
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
    if math.fabs(r) < 10.0:
        badr = badr + 1
        #print "r, z = " + str(r) + ", " + str(z)
        hsmallr.Fill(z)
        
    hxy.Fill(x, y)
    hrz.Fill(z, r)
    
    if y > 0: hrz2.Fill(z, r)
    else: hrz2.Fill(z, r*-1.0)
    
    hr2d.Fill(r)
    hr3d.Fill(sqrt(x*x + y*y + z*z))
    
    count = count + 1
    #if count > 10: break

hxy.SetMarkerSize(0.3)
hrz.SetMarkerSize(0.3)
hrz2.SetMarkerSize(0.3)
hr2d.SetMarkerSize(0.3)
hr3d.SetMarkerSize(0.3)
#hrz.SetMarkerStyle(7)


# Handle input 2
# Distance in cm
hxy_i2 = TH2F('hxy_i2','',320,-1600,1600,320,-1600,1600)
hrz_i2 = TH2F('hrz_i2','',600,-3000,3000,160,0,1600)
hrz2_i2 = TH2F('hrz2_i2','',600,-3000,3000,320,-1600,1600)
hr2d_i2 = TH1F('hr2d_i2','',800,0,1600)
hr3d_i2 = TH1F('hr3d_i2','',1700,0,3400)
hsmallr_i2 = TH1F('hsmallr_i2','z for points with r < 10.0cm',3000,-3000,3000)

count = 0
maxx = 0
maxy = 0
maxz = 0
maxr = 0
badr = 0
for line in input2:
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
    if math.fabs(r) < 10.0:
        badr = badr + 1
        #print "r, z = " + str(r) + ", " + str(z)
        hsmallr_i2.Fill(z)
        
    hxy_i2.Fill(x, y)
    hrz_i2.Fill(z, r)
    
    if y > 0: hrz2_i2.Fill(z, r)
    else: hrz2_i2.Fill(z, r*-1.0)
    
    hr2d_i2.Fill(r)
    hr3d_i2.Fill(sqrt(x*x + y*y + z*z))
    
    count = count + 1
    #if count > 10: break

hxy_i2.SetMarkerSize(0.3)
hrz_i2.SetMarkerSize(0.3)
hrz2_i2.SetMarkerSize(0.3)
hr2d_i2.SetMarkerSize(0.3)
hr3d_i2.SetMarkerSize(0.3)
#hrz.SetMarkerStyle(7)

#draw and save the 2d r histogram
hr2d_i2.SetLineColor(kRed)
c = TCanvas()
hr2d_i2.DrawNormalized()
hr2d.DrawNormalized('same')
c.SaveAs(parts[0] + 'overlay.pdf')


output.Write()
#input1.close()
#input2.close()
print 'Total: ' + str(count) + ' events.'
print 'Max (x, y, z, r) : (' + str(maxx) + ', ' + str(maxy) + ', ' + str(maxz) + ', ' + str(maxr) + ')'
print 'Points with r < 10.0 cm: ' + str(badr) + '(' + str(100.0*badr/count) + '%)'
print 'Finished.'


