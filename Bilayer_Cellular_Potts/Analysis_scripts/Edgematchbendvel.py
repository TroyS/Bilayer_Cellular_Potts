#This will be used to compare the behavior of the bending runs to the full runs

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import math
import statistics
import Edgedatfunctions as ed
import MSDfunctions as msd
import BCPM_Classes as clas


import plotly.graph_objects as go

from matplotlib.colors import TABLEAU_COLORS, same_color
from matplotlib.pyplot import cm
from matplotlib import cm


lege = []

tle = "Edge match ratio vs timestep"


##############################################################

pval = ["p0_4p25", "p0_4p75", "p0_5p0"]

pval2 =['$p_{0}=4.25$','$p_{0}=4.75$','$p_{0}=5$']


pvald = [4.25, 4.75, 5.0] 

pnum = [135, 150, 159]

nshape = len(pvald)

Bi = ["6.0", "14.0"]


nbicouple = len(Bi)

Birange = [3.0,7.0]
BirangeS = ["3.0", "7.0"]

#Birangemagnif = [5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.75, 9.0, 9.25, 9.5, 9.75]

#snum = ["seed50","seed51", "seed52", "seed53", "seed54", "seed55", "seed57", "seed58", "seed59", "seed60", "seed61", "seed62", "seed63", "seed64", "seed65", "seed67", "seed68", "seed69", "seed70", "seed71"]

snum = ["seed50", "seed51", "seed52", "seed53", "seed54"]
#snum = ["seed55", "seed57" ,"seed58" ,"seed59", "seed60"]


nseeds = len(snum)


vlen = len(pval)
blen = len(Birange)


Binum = [0]*(vlen*nseeds)
colorssc = [0]*(vlen*nseeds)
colorsm = [0]*(vlen*nseeds)

colorwhl = ["red", "blue", "green" ,"purple", "orange", "pink", "brown", "gray", "olive", "cyan", "gold", "navy"]

#for i in range(0,nseeds):
#    for j in range(0,vlen):
#        colorssc[i+nseeds*j] = 10*(vlen+1)
#        colorsm[i+nseeds*j] = colorwhl[vlen]
#        Binum[i+nseeds*j] = Birange[i]




#Now the CFC runs

CFCtitle = ["40cellA1000"]

CFCshape = ["p0_4p25", "p0_4p75", "p0_5"]

CFCBi = ["6.0", "14.0"]

#CFCseeds = ["seed2", "seed3", "seed4", "seed5", "seed6"]

CFCseeds = ["seed2", "seed3", "seed4", "seed5", "seed6"]


Bendseeds = ["seed2", "seed3", "seed4", "seed5", "seed6", "seed7", "seed8", "seed9", "seed10", "seed11"]


BendseedsR2 = ["seed2", "seed3", "seed4", "seed5"]

updownseeds = ["seed2updown", "seed3updown", "seed4updown", "seed5updown", "seed6updown"]

nCFC = len(CFCseeds)


nBend = len(Bendseeds)

nBendR2 = len(BendseedsR2)

#Let's also work on a scatter plot
edgemax = []


#pval = ["p0_4p5"]

nlines = 0

xvals = np.array(pvald)
yvals = np.array(Birange)
#yvalsmagnif = np.array(Birangemagnif)
[X,Y] = np.meshgrid(xvals, yvals)
levels = np.linspace(len(xvals),len(yvals))
zcont = np.zeros((len(yvals),len(xvals)))
zcontstd = np.zeros((len(yvals),len(xvals)))
#[X2,Y2] = np.meshgrid(xvals,yvalsmagnif)
#levels2 = np.linspace(len(xvals),len(yvalsmagnif))
#zcontmagnif = np.zeros((len(yvalsmagnif), len(xvals)) )


#First let's get the number of lines in each file

sampfile = pval[0]+snum[0]+"Bi_" + Bi[0] + "edgematch.dat"

nlines1 = ed.getlines(sampfile)

sampfile2 = CFCtitle[0] + CFCshape[0] + CFCseeds[0] + "Bi_" + CFCBi[0] + "edgematch.dat"

nlines2 = ed.getlines(sampfile2)

sampfile3 = CFCtitle[0] + CFCshape[0] + "CFC" + CFCseeds[0] + "Bi_" + CFCBi[0] +"edgematch.dat"  

nlines3 = ed.getlines(sampfile3)


sampfile31 = CFCtitle[0] + CFCshape[0] + "CFC" + CFCseeds[0] + "Bi_" + CFCBi[0] +"region.dat"

npoints, nxx, latti, ncells = msd.getinforegion(sampfile31)


rad3 = math.sqrt(3)
Lat = 1
if latti == 2:
    Lat = math.sqrt(2/rad3)

ny = npoints/nxx
leny = rad3*Lat/2


bendlen = [ [9200, 18400, 27800, 37000, 40000], [10200, 20600, 30800, 40000], [10800, 21800, 32600, 40000] ]


updownlen = [8400.0, 10200.0, 12600.0, 15400.0, 16800.0, 20500.0, 20600.0, 21000.0, 21100.0, 24200.0, 25200.0, 25300.0, 25700.0, 28400.0, 29500.0, 30300.0, 32600.0, 33800.0, 35400.0, 36900.0, 40000.0]


sampfile4 = CFCtitle[0] + CFCshape[0] + "CFC" + CFCseeds[0] + "Bi_" + CFCBi[0] +"Bendingregimes.dat"

nlines4 = ed.getlinesneigh(sampfile4)

#sampfile5 = CFCtitle[0] + CFCshape[0] + CFCseeds[0] +  "updown" + "Bi_" + CFCBi[0] +"regextensions.txt"

sampfile5 = CFCtitle[0] + CFCshape[0] + CFCseeds[0] +  "updown" + "Bi_" + CFCBi[0] +"Bendingregimes.dat"

nlines5 = ed.getlinesneigh(sampfile5)


CFClabel = ["first_extension", "second_extension", "third_extension"]

timeslabel = [["first_extension", "second_extension", "third_extension", "fourth_extension"], ["first_extension", "second_extension", "third_extension"], ["first_extension", "second_extension", "third_extension"]]


#cmapf = plt.get_cmap('rainbow',nbicouple)
cmapf = plt.get_cmap('rainbow',nseeds)
endlist1 = nlines1-math.floor(nlines1/10)


endlist2 = nlines2-math.floor(nlines2/10)

endrise = nlines1-math.floor(nlines1/100)
#############################################################
#Now we analyze the data
#First, the normal cases


edgeregular = clas.Edgematchdat

edgeregular.tsdat, edgeregular.edgedat, edgeregular.edgediff, edgeregular.name = ed.getalldata(nlines1, nshape, nbicouple,nseeds, pval, Bi,snum, "edgematch.dat")


#####################################################

#Then, the long aspect ratio normal cases

edgelongasp = clas.Edgematchdat

CFCadjust =  ["40cellA1000p0_4p25", "40cellA1000p0_4p75", "40cellA1000p0_5"]

CFCadjust2 =  ["40cellA1000p0_4p25CFC", "40cellA1000p0_4p75CFC", "40cellA1000p0_5CFC"]

CFCadjust3 =  ["40cellA1000p0_4p25bendregionII", "40cellA1000p0_4p75bendregionII", "40cellA1000p0_5bendregionII"]



edgelongasp.tsdat, edgelongasp.edge, edgelongasp.edgediff, edgelongasp.name = ed.getalldata(nlines1, len(CFCshape), len(CFCBi),nCFC, CFCadjust2,CFCBi,CFCseeds, "edgematch.dat")


#################
#Then, the bending runs


edgeCFC = clas.Edgematchdat

edgeCFC.tsdat, edgeCFC.edge, edgeCFC.edgediff, edgeCFC.name = ed.getalldata(nlines3, len(CFCshape), len(CFCBi),nBend, CFCadjust2,CFCBi,Bendseeds, "edgematch.dat")


#bendCFC = clas.bendinfoseeds 

#bendCFC.bspeed, bendCFC.bentime, bendCFC.name = ed.getallbendspeedCFC(nlines4,len(CFCshape), nBend, len(CFCBi), CFCadjust2,Bendseeds,CFCBi, "Bendingregimes.dat",leny,bendlen,nxx)   


#bendCFCavg = clas.bendinfoavg

#bendCFCavgbspeed, bendCFCavgspeedstd = ed.avgallbendspeedCFC(len(CFCshape),nBend,len(CFCBi), bendlen, bendCFC.bspeed)


#print(bendCFCavgbspeed)


#print(bendCFCavg.bspeed[0][0])

bendCFCR2 = clas.bendinfoseeds

bendCFCR2.bspeed, bendCFCR2.bentime, bendCFCR2.name = ed.getallbendspeedCFC(nlines4,len(CFCshape), nBendR2, len(CFCBi), CFCadjust3,BendseedsR2,CFCBi, "Bendingregimes.dat",leny,bendlen,nxx)

bendCFCR2avg = clas.bendinfoavg

bendCFCR2avg.bspeed, bendCFCR2avg.speedstd = ed.avgallbendspeedCFC(len(CFCshape),nBendR2,len(CFCBi), bendlen, bendCFCR2.bspeed)

#bendCFCR2avg.get_bspeed(bendCFCR2avgbspeed)
#bendCFCR2avg.get_speedstd(bendCFCR2avgspeedstd)


#print(bendCFCR2avg.bspeed[0][0])

################

bendupdown = clas.bendvalues

bendupdown.regIpnts,bendupdown.name = ed.getallexttime(nlines5,len(CFCshape),nCFC,len(CFCBi),CFCadjust,updownseeds,CFCBi, "Bendingregimes.dat")

#print(bendupdown.regIpnts)

figloop = 0




#############################################

#############################################################


bwidth = 0.25

multiply = 0

bcolors = ['b', 'g']


#for p0loop in range(0,nshape):
#    labels = np.arange(len(timeslabel[p0loop]))
#    labels3 = [x + bwidth for x in labels]

#    fig = plt.subplots(figsize=(6,6))

#    for b0loop in range(0,len(Bi)):



#        xparams = [ bendCFCavg.bspeed[p0loop,b0loop], bup2diffavg[p0loop,b0loop], bup3diffavg[p0loop,b0loop], bup4diffavg[p0loop,biloop] ]
#        error = [bup1diffstd[p0loop,b0loop], bup2diffstd[p0loop,b0loop], bup3diffstd[p0loop,b0loop], bup4diffstd[p0loop,biloop] ]
#        bendspeeds = bendCFCavg.bspeed[p0loop][b0loop]
#        error = bendCFCavg.speedstd[p0loop][b0loop]
#        labels3 = [x +(b0loop+1)*bwidth for x in labels]
#        print(labels3)
#        print("p0loop", p0loop, " and b0loop ", b0loop)
#        print(bendspeeds)
#        plt.plot(labels3,bendspeeds )
#        plt.bar(labels3,bendspeeds, yerr = error, capsize=3, width=bwidth,label = BirangeS[b0loop] )

#        plt.grid()
#        plt.tick_params(axis='both', labelsize=28)
#        plt.xticks([])

#        print(updownrisetavg[p0loop,b0loop])

#        plt.axhline(y=updownrisetavg[p0loop,b0loop], color = bcolors[b0loop] )
##
#    plt.legend(loc='lower left', bbox_to_anchor=(0.99,0), fontsize=22)
#    plt.title("Up extension time differences")


#figloop = figloop + 1

############################################################
bwidth = 0.25

multiply = 0

bcolors = ['b', 'g']


for p0loop in range(0,nshape):
    labels = np.arange(len(timeslabel[p0loop]))
    labels3 = [x + bwidth for x in labels]

    fig = plt.subplots(figsize=(6,6))

    for b0loop in range(0,len(Bi)):


#        xparams = [ bendCFCavg.bspeed[p0loop,b0loop], bup2diffavg[p0loop,b0loop], bup3diffavg[p0loop,b0loop], bup4diffavg[p0loop,biloop] ]
#        error = [bup1diffstd[p0loop,b0loop], bup2diffstd[p0loop,b0loop], bup3diffstd[p0loop,b0loop], bup4diffstd[p0loop,biloop] ]
        error = bendCFCR2avg.speedstd[p0loop][b0loop]
        labels3 = [x +(b0loop+1)*bwidth for x in labels]
        plt.plot(labels3,bendCFCR2avg.bspeed[p0loop][b0loop] )
        plt.bar(labels3,bendCFCR2avg.bspeed[p0loop][b0loop], yerr = error, capsize=3, width=bwidth,label = BirangeS[b0loop] )

        plt.grid()
        plt.tick_params(axis='both', labelsize=28)
        plt.xticks([])

#        print(updownrisetavg[p0loop,b0loop])

#        plt.axhline(y=updownrisetavg[p0loop,b0loop], color = bcolors[b0loop] )
#
    plt.legend(loc='lower left', bbox_to_anchor=(0.99,0), fontsize=22)
#    plt.title("Up extension time differences")


figloop = figloop + 1





############################################################


plt.show()

