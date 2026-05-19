#This will be used to plot the bilayer edge match ratio and neighbor changes, along with a couple of regime maps

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import math
import statistics
import Edgedatfunctions as ed
import MSDfunctions as msd

import plotly.graph_objects as go

from matplotlib.colors import TABLEAU_COLORS, same_color
from matplotlib.pyplot import cm
from matplotlib import cm
import time

lege = []

tle = "Edge match ratio vs timestep"


##############################################################

pval = ["p0_4", "p0_4p12", "p0_4p25", "p0_4p37","p0_4p5", "p0_4p62", "p0_4p75", "p0_5p0", "p0_5p5", "p0_6"]

pval2 =['$p_{0}=4$', '$p_{0}=4.12$', '$p_{0}=4.25$','$p_{0}=4.37$','$p_{0}=4.5$', '$p_{0}=4.62', '$p_{0}=4.75$', '$p_{0}=5$', '$p_{0}=5.5$', '$p_{0}=6.0$']


pvald = [4.0, 4.12, 4.25, 4.37, 4.5, 4.62, 4.75, 5.0, 5.5, 6.0] 



nshape = len(pvald)

Bi = ["1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "4.0", "4.5", "5.0", "5.5", "6.0", "6.5", "7.0", "7.5", "8.0", "8.5", "9.0", "9.5", "10.0", "10.5", "11.0", "11.5", "12.0", "12.5", "13.0", "13.5", "14.0", "14.5", "15.0", "15.5", "16.0", "16.5", "17.0", "17.5", "18.0", "18.5", "19.0", "19.5"]

nbicouple = len(Bi)

Birange = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75]

Birangemagnif = [5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75]

vlen = len(pval)
blen = len(Birange)

Binum = [0]*(vlen*blen)
colorssc = [0]*(vlen*blen)
colorsm = [0]*(vlen*blen)

colorwhl = ["red", "blue", "green" ,"purple", "orange", "pink", "brown", "gray", "olive", "cyan", "gold", "navy", "lightsalmon", "chocolate", "navajowhite", "cornsilk", "beige", "lime", "aquamarine", "slateblue", "fuchsia"]

for i in range(0,blen):
    for j in range(0,vlen):
        colorssc[i+blen*j] = 10*(vlen+1)
        colorsm[i+blen*j] = colorwhl[vlen]
        Binum[i+blen*j] = Birange[i]



#Let's also work on a scatter plot
edgemax = []

nlines = 0


snum = ["seed50","seed51","seed52", "seed53", "seed54", "seed55", "seed57", "seed58", "seed59", "seed60", "seed61", "seed62", "seed63", "seed64", "seed65", "seed67", "seed68", "seed69", "seed70", "seed71"]


nseeds = len(snum)


xvals = np.array(pvald)
yvals = np.array(Birange)
yvalsmagnif = np.array(Birangemagnif)
[X,Y] = np.meshgrid(xvals, yvals)
levels = np.linspace(len(xvals),len(yvals))
zcont = np.zeros((len(yvals),len(xvals)))
zcontstd = np.zeros((len(yvals),len(xvals)))
[X2,Y2] = np.meshgrid(xvals,yvalsmagnif)
levels2 = np.linspace(len(xvals),len(yvalsmagnif))
zcontmagnif = np.zeros((len(yvalsmagnif), len(xvals)) )
[X3,Y3] = np.meshgrid(xvals,yvals)
levels3 = np.linspace(len(xvals),len(yvals))
rdistcont = np.zeros((len(yvals),len(xvals)))
rdispcont = np.zeros((len(yvals),len(xvals)))
distoverdisp = np.zeros((len(yvals),len(xvals)))
velcont = np.zeros((len(yvals),len(xvals)))


#First let's get the number of lines in each file

sampfile = pval[0]+snum[0]+"Bi_" + Bi[0] + "edgematch.dat"
sampfile2 = pval[0]+snum[0]+"Bi_" + Bi[0] + "neighchangesfix.dat"
sampfile3 = "p0_4seed50Bi_1.0region.dat"
nlines = ed.getlines(sampfile)
nlines2 = ed.getlinesneigh(sampfile2)
npoints, nxx, latti, ncells = msd.getinforegion(sampfile3)

rad3 = math.sqrt(3)
Lat = 1
if latti == 2:
    Lat = math.sqrt(2/rad3)

ny = npoints/nxx
Xmin = Lat
Xmax = (nxx+0.5)*Lat
Ymin = rad3*Lat/2
Ymax = rad3*Lat*ny/2

skip2 = 10

tsdat = np.zeros((nlines,nshape,nseeds,nbicouple))
edgedat = np.zeros((nlines,nshape,nseeds,nbicouple))
edgedatalt = np.zeros((nlines,nshape,nbicouple,nseeds))
edgediff = np.zeros((nlines,nshape,nseeds,nbicouple))
#edgediff2 = np.zeros((nlines/skip2,nshape,nseeds,nbicouple))
edgediffmax = np.zeros((nshape,nseeds,nbicouple))
edgediffmaxavg = np.zeros((nshape,nbicouple))
edgestd = np.zeros((nlines,nshape,nseeds,nbicouple))
neichange1 = np.zeros((nlines2,nshape,nseeds,nbicouple))
neichange2 = np.zeros((nlines2,nshape,nseeds,nbicouple))
neichangesboth = np.zeros((nlines2,nshape,nseeds,nbicouple))
neichangesum = np.zeros((nshape,nseeds,nbicouple))
neichangeatrise = np.zeros((nshape,nseeds,nbicouple))
neichangescontour = np.zeros((nbicouple,nshape))
edgeneicorrel = np.zeros((nshape,nseeds,nbicouple))

edgerisetime = np.zeros((nshape,nseeds,nbicouple))
edgerisetimeseeds = np.zeros((nshape,nbicouple))
edgerisetimestd = np.zeros((nshape,nbicouple))
edgerisetimecontour = np.zeros((nbicouple,nshape))
edgesettime = np.zeros((nlines,nshape,nseeds,nbicouple))

Time = np.zeros(nlines2)

xcomrawl1 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))
xdistcell1 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))
xdispcell1 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))
xabsdispcell1 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))

ycomrawl1 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))
ydistcell1 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))
ydispcell1 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))
yabsdispcell1 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))


xcomrawl2 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))
xdistcell2 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))
xdispcell2 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))
xabsdispcell2 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))

ycomrawl2 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))
ydistcell2 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))
ydispcell2 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))
yabsdispcell2 = np.zeros((nlines2,ncells,nshape,nseeds,nbicouple))

xdistseedsl1 = np.zeros((nlines2,nshape,nseeds,nbicouple))
ydistseedsl1 = np.zeros((nlines2,nshape,nseeds,nbicouple))
xabsdispseedsl1 = np.zeros((nlines2,nshape,nseeds,nbicouple))
yabsdispseedsl1 = np.zeros((nlines2,nshape,nseeds,nbicouple))

xdistseedsl2 = np.zeros((nlines2,nshape,nseeds,nbicouple))
ydistseedsl2 = np.zeros((nlines2,nshape,nseeds,nbicouple))
xabsdispseedsl2 = np.zeros((nlines2,nshape,nseeds,nbicouple))
yabsdispseedsl2 = np.zeros((nlines2,nshape,nseeds,nbicouple))

rdistseeds = np.zeros((nlines2,nshape,nseeds,nbicouple))
rabsdispseeds = np.zeros((nlines2,nshape,nseeds,nbicouple))

rdistatrise = np.zeros((nshape,nseeds,nbicouple))
rdispatrise = np.zeros((nshape,nseeds,nbicouple))

rdist = np.zeros((nlines2,nshape,nbicouple))

ravgdistrise = np.zeros((nshape,nbicouple))

ravgdistatrisestd = np.zeros((nshape,nbicouple))


rabsdisp = np.zeros((nlines2,nshape,nbicouple))

ravgdisprise = np.zeros((nshape,nbicouple))

ravgdispatrisestd = np.zeros((nshape,nbicouple))

cmapf = plt.get_cmap('rainbow',nbicouple)
cmaptemp = plt.get_cmap('tab20',nseeds)
endlist = nlines-math.floor(nlines/10)

endrise = nlines-math.floor(nlines/100)

#############################################################
#Now we get the data

getdatatimestar = time.perf_counter()

for p0loop in range(0,nshape):

    edgeAvg = []
    seeds = []
    edgeStd = []

    for seedloop in range(0,nseeds):
        for biloop in range(0,nbicouple):
            fnme =pval[p0loop]+snum[seedloop]+"Bi_" + Bi[biloop]+"edgematch.dat"
            lege.append(fnme)
            nchgfil = pval[p0loop] + snum[seedloop]+"Bi_" + Bi[biloop]+"neighchangesfix.dat"

            dispfil = pval[p0loop] + snum[seedloop]+"Bi_" + Bi[biloop]+"COMviadist.dat"

            tPB = ed.getdata(fnme,nlines)

            tsdat[:,p0loop,seedloop,biloop] = tPB[:,0]
            edgedat[:,p0loop,seedloop,biloop] = tPB[:,1]
            edgedatalt[:,p0loop,biloop,seedloop] = tPB[:,1]
            edgediff[1:,p0loop,seedloop,biloop] = np.diff(edgedat[:,p0loop,seedloop,biloop])

            edgediffmax[p0loop,seedloop,biloop] = np.max(abs(edgediff[:,p0loop,seedloop,biloop]))

            edgerisetime[p0loop,seedloop,biloop] = ed.risetimeavg(edgedat[:,p0loop,seedloop,biloop], endlist,nlines)

            neichange1[:,p0loop,seedloop,biloop], neichange2[:,p0loop,seedloop,biloop] = ed.getneighdata(nchgfil, nlines2)
            neichangesboth[:,p0loop,seedloop,biloop] = neichange1[:,p0loop,seedloop,biloop] + neichange2[:,p0loop,seedloop,biloop]
            neichangesum[p0loop,seedloop,biloop] = np.sum(neichangesboth[:,p0loop,seedloop,biloop])

            Time,xcomrawl1[:,:,p0loop,seedloop,biloop], ycomrawl1[:,:,p0loop,seedloop,biloop], xcomrawl2[:,:,p0loop,seedloop,biloop], ycomrawl2[:,:,p0loop,seedloop,biloop]= msd.getCOMdata(dispfil,nlines2,ncells)

#            st1 = pval[p0loop]
#            st2 = Bi[biloop]
#            st3 = snum[seedloop]
#            st4 = st1 + " " + st2 + " " + st3
#            print(st4)

            xdistcell1[:,:,p0loop,seedloop,biloop], ydistcell1[:,:,p0loop,seedloop,biloop],xdistcell2[:,:,p0loop,seedloop,biloop], ydistcell2[:,:,p0loop,seedloop,biloop] = msd.getdistallcellsalldim(xcomrawl1[:,:,p0loop,seedloop,biloop],ycomrawl1[:,:,p0loop,seedloop,biloop],xcomrawl2[:,:,p0loop,seedloop,biloop],ycomrawl2[:,:,p0loop,seedloop,biloop], nlines2,ncells)


            xdispcell1[:,:,p0loop,seedloop,biloop], ydispcell1[:,:,p0loop,seedloop,biloop],xdispcell2[:,:,p0loop,seedloop,biloop], ydispcell2[:,:,p0loop,seedloop,biloop] = msd.getdispallcellsalldim(xcomrawl1[:,:,p0loop,seedloop,biloop],ycomrawl1[:,:,p0loop,seedloop,biloop],xcomrawl2[:,:,p0loop,seedloop,biloop],ycomrawl2[:,:,p0loop,seedloop,biloop], nlines2,ncells)


#            corrtemp = np.correlate(edgediff[::10,p0loop,seedloop,biloop], neichangesboth[:nlines2-1,p0loop,seedloop,biloop], mode='full')

#            maxval = np.max(abs(corrtemp))

#            edgeneicorrel[p0loop,seedloop,biloop] = np.argmax(corrtemp == maxval)

            riseind = int(math.ceil(edgerisetime[p0loop,seedloop,biloop]/10))

            neichangeatrise[p0loop,seedloop,biloop] = np.sum(neichangesboth[0:riseind,p0loop,seedloop,biloop])


    for biloop in range(0,nbicouple):

        edgediffmaxavg[p0loop,biloop] = np.max(edgediffmax[p0loop,:,biloop])


#print(edgeneicorrel)

#######################################################


print("The max diff  is ", edgediffmax[:,:,:])


#Now we analyze



for p0loop in range(0,nshape):
    for seedloop in range(0,nseeds):
        for biloop in range(0,nbicouple):


            xdistseedsl1[:,p0loop,seedloop,biloop] = msd.getabsdistnospike2(xcomrawl1[:,:,p0loop,seedloop,biloop], xdistcell1[:,:,p0loop,seedloop,biloop],nlines2,ncells,Xmin,Xmax)
            ydistseedsl1[:,p0loop,seedloop,biloop] = msd.getabsdistnospike2(ycomrawl1[:,:,p0loop,seedloop,biloop], ydistcell1[:,:,p0loop,seedloop,biloop],nlines2,ncells,Ymin,Ymax)
            xdistseedsl2[:,p0loop,seedloop,biloop] = msd.getabsdistnospike2(xcomrawl2[:,:,p0loop,seedloop,biloop], xdistcell2[:,:,p0loop,seedloop,biloop],nlines2,ncells,Xmin,Xmax)
            ydistseedsl2[:,p0loop,seedloop,biloop] = msd.getabsdistnospike2(ycomrawl2[:,:,p0loop,seedloop,biloop], ydistcell2[:,:,p0loop,seedloop,biloop],nlines2,ncells,Ymin,Ymax)


            if np.count_nonzero(xdistseedsl1[:,p0loop,seedloop,biloop]) == 0:
                print("No cells found in p0 ")
                print(pval[p0loop])
                print("seed")
                print(snum[seedloop])
                print("Bi")
                print(Bi[biloop])
                print("X lay 1")


            if np.count_nonzero(ydistseedsl1[:,p0loop,seedloop,biloop]) == 0:
                print("No cells found in p0 ")
                print(pval[p0loop])
                print("seed")
                print(snum[seedloop])
                print("Bi")
                print(Bi[biloop])
                print("Y lay 1")


            if np.count_nonzero(xdistseedsl2[:,p0loop,seedloop,biloop]) == 0:
                print("No cells found in p0 ")
                print(pval[p0loop])
                print("seed")
                print(snum[seedloop])
                print("Bi")
                print(Bi[biloop])
                print("X lay 2")


            if np.count_nonzero(ydistseedsl2[:,p0loop,seedloop,biloop]) == 0:
                print("No cells found in p0 ")
                print(pval[p0loop])
                print("seed")
                print(snum[seedloop])
                print("Bi")
                print(Bi[biloop])
                print("Y lay 2")


            xabsdispseedsl1[:,p0loop,seedloop,biloop] = msd.getabsdispnospike(xcomrawl1[:,:,p0loop,seedloop,biloop], xdispcell1[:,:,p0loop,seedloop,biloop],nlines2,ncells,Xmin,Xmax)
            yabsdispseedsl1[:,p0loop,seedloop,biloop] = msd.getabsdispnospike(ycomrawl1[:,:,p0loop,seedloop,biloop], ydispcell1[:,:,p0loop,seedloop,biloop],nlines2,ncells,Ymin,Ymax)
            xabsdispseedsl2[:,p0loop,seedloop,biloop] = msd.getabsdispnospike(xcomrawl2[:,:,p0loop,seedloop,biloop], xdispcell2[:,:,p0loop,seedloop,biloop],nlines2,ncells,Xmin,Xmax)
            yabsdispseedsl2[:,p0loop,seedloop,biloop] = msd.getabsdispnospike(ycomrawl2[:,:,p0loop,seedloop,biloop], ydispcell2[:,:,p0loop,seedloop,biloop],nlines2,ncells,Ymin,Ymax)



            rdistseeds[:,p0loop,seedloop,biloop] = ( (xdistseedsl1[:,p0loop,seedloop,biloop] + xdistseedsl2[:,p0loop,seedloop,biloop]) **2 )/4 + ( ( ydistseedsl1[:,p0loop,seedloop,biloop] + ydistseedsl2[:,p0loop,seedloop,biloop] ) **2 )/4
            
            rdistseeds[:,p0loop,seedloop,biloop] = np.sqrt(rdistseeds[:,p0loop,seedloop,biloop] )

            riseind = int(math.ceil(edgerisetime[p0loop,seedloop,biloop]/10))

            rdistatrise[p0loop,seedloop,biloop] = rdistseeds[riseind,p0loop,seedloop,biloop]


            rabsdispseeds[:,p0loop,seedloop,biloop] = ( ( xabsdispseedsl1[:,p0loop,seedloop,biloop] + xabsdispseedsl2[:,p0loop,seedloop,biloop] ) **2 )/4 + ( ( yabsdispseedsl1[:,p0loop,seedloop,biloop] + yabsdispseedsl2[:,p0loop,seedloop,biloop] ) **2 )/4

            rabsdispseeds[:,p0loop,seedloop,biloop] = np.sqrt(rabsdispseeds[:,p0loop,seedloop,biloop])


            rdispatrise[p0loop,seedloop,biloop] = rabsdispseeds[riseind,p0loop,seedloop,biloop]


getdatatimeend = time.perf_counter()
getdatatime = getdatatimestar - getdatatimeend
print(f"Get data time: {getdatatime:0.4f} seconds")



############################################
#Now we get the averages 

getavgtimestar = time.perf_counter()

#First the average over seeds

tsavg = tsdat[:,0,0,0]
EMseedsavg = np.zeros((nlines,nshape,nbicouple))
EMseedsstd = np.zeros((nlines,nshape,nbicouple))

NCseeds = np.zeros((nlines2,nshape,nbicouple))
NCseedsstd = np.zeros((nlines2,nshape,nbicouple))


#Then the average over the last couple of timesteps
Regimemapavg, Regimemapstd = ed.Avgoverendtime(edgedatalt,nshape,nbicouple,nseeds,endlist,nlines)

#Print out the std info

precision = 3

for p0loop in range(0,nshape):
    for biloop in range(0,nbicouple):
        zcont[biloop,p0loop] = Regimemapavg[p0loop,biloop]
        edgerisetimeseeds[p0loop,biloop] =np.mean(edgerisetime[p0loop,:,biloop])
        edgerisetimestd[p0loop,biloop] = np.std(edgerisetime[p0loop,:,biloop])
        edgerisetimecontour[biloop,p0loop] = edgerisetimeseeds[p0loop,biloop]
        rdistcont[biloop,p0loop] = np.mean(rdistatrise[p0loop,:,biloop])
        rdispcont[biloop,p0loop] = np.mean(rdispatrise[p0loop,:,biloop])
        ravgdistatrisestd[p0loop,biloop] = np.std(rdistatrise[p0loop,:,biloop])
        neichangescontour[biloop,p0loop] = np.sum(neichangeatrise[p0loop,:,biloop])/float(nseeds)
        velcont[biloop,p0loop] = rdistcont[biloop,p0loop]/edgerisetimeseeds[p0loop,biloop]
        distoverdisp[biloop,p0loop] = rdistcont[biloop,p0loop]/rdispcont[biloop,p0loop]
        if biloop > 19:
            zcontmagnif[biloop-20,p0loop] = Regimemapavg[p0loop,biloop]
            

for lines in range(0,nlines2):
    for p0 in range(0,nshape):
        for bi in range(0,nbicouple):
            rdist[lines,p0,bi] = np.mean(rdistseeds[lines,p0,:,bi])

            NCseeds[lines,p0,bi] = np.sum(neichangesboth[lines,p0,:,bi])
            NCseedsstd[lines,p0,bi] = np.std(neichangesboth[lines,p0,:,bi])


for lines in range(0,nlines):
    for p0 in range(0,nshape):
        for bi in range(0,nbicouple):
            EMseedsavg[lines,p0,bi] = np.mean(edgedat[lines,p0,:,bi])
            EMseedsstd[lines,p0,bi] = np.std(edgedat[lines,p0,:,bi])



getavgtimeend = time.perf_counter()
getavgtime = getavgtimestar-getavgtimeend
print(f"Average time : {getavgtime:0.4f} seconds")

print(np.max(rdistcont))
print(np.max(rdispcont))
print(np.max(velcont))


#########################################################
#Now we plot

figloop = 1

nbihalf = nbicouple//2

pchk = 0
#bchk = 5 #Bi = 1.75
bchk = 12 #Bi = 3.5 
bi4p5 = 16 #Bi = 4.5
p5p0 = 7
skip = 10
tsavgskp = tsavg[::skip]


################################################################

fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )
custom_yticks1=[0.0, 0.3, 0.6, 0.9]
custom_yrange=[0.0,1.0]
custom_yticksens=[0,4,8,12]


for seedloop in range(0,nseeds): 

    PBseedskp = edgedat[::skip,pchk,seedloop,bchk]
    ax1.plot(tsavgskp,PBseedskp, color = cmaptemp(seedloop), label = snum[seedloop] )
    ax1.grid()
    ax1.tick_params(axis='both', labelsize = 28)
    ax1.set_ylim(custom_yrange)

    ax2.plot(Time,neichangesboth[:,pchk,seedloop,bchk], color = cmaptemp(seedloop), label = snum[seedloop])
    ax2.grid()
    ax2.tick_params(axis='both', labelsize = 28)


figloop = figloop + 1


########################################################
#The edge match difference
fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )
#custom_yticks1=[0.0, 0.3, 0.6, 0.9]
custom_yrange=[0.0,1.0]
custom_yticksens=[0,4,8,12]


for seedloop in range(0,nseeds): 

    #PBseedskp = edgedat[::skip,pchk,seedloop,bchk]
    PBseedskp = edgediff[::skip,pchk,seedloop,bchk]
    ax1.plot(tsavgskp,PBseedskp, color = cmaptemp(seedloop), label = snum[seedloop] )
    ax1.grid()
    ax1.tick_params(axis='both', labelsize = 28)
#    ax1.set_ylim(custom_yrange)

    ax2.plot(Time,neichangesboth[:,pchk,seedloop,bchk], color = cmaptemp(seedloop), label = snum[seedloop])
    ax2.grid()
    ax2.tick_params(axis='both', labelsize = 28)


figloop = figloop + 1


#############################################################

fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )
custom_yticks1=[0.0, 0.3, 0.6, 0.9]
custom_yrange=[0.0,1.0]

for seedloop in range(0,nseeds): 

    PBseedskp = edgedat[::skip,p5p0,seedloop,bchk]
    ax1.plot(tsavgskp,PBseedskp, color = cmaptemp(seedloop), label = snum[seedloop] )
    ax1.grid()
    ax1.tick_params(axis='both', labelsize = 28)
    ax1.set_ylim(custom_yrange)

    ax2.plot(Time,neichangesboth[:,p5p0,seedloop,bchk], color = cmaptemp(seedloop), label = snum[seedloop])
    ax2.grid()
    ax2.tick_params(axis='both', labelsize = 28)


figloop = figloop + 1



###########################################################

fig, (ax1, ax2) = plt.subplots(2,1, figsize=(6,6) )

for seedloop in range(0,nseeds): 

    PBseedskp = edgedat[::skip,nshape-1,seedloop,bchk]
    ax1.plot(tsavgskp,PBseedskp, color = cmaptemp(seedloop), label = snum[seedloop] )
    ax1.grid()
    ax1.tick_params(axis='both', labelsize = 28)
    ax1.set_ylim(custom_yrange)

    ax2.plot(Time,neichangesboth[:,nshape-1,seedloop,bchk], color = cmaptemp(seedloop), label = snum[seedloop])
    ax2.grid()
    ax2.tick_params(axis='both', labelsize = 28)



figloop = figloop + 1

#########################################
fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )

#plt.title("p0_4 ensemble averages")

nbihalf = nbicouple//2

for bi in range(0,nbihalf):
    bi2 = 2*bi

    PBskp = EMseedsavg[::skip,0,bi2]
    ax1.plot(tsavgskp,PBskp, color = cmapf(bi2), label = Birange[bi2] )
    ax1.grid()

    ax1.tick_params(axis='both', labelsize = 28)
    ax1.legend(loc='lower left', bbox_to_anchor=(0.99,-1.35), fontsize=22)
    ax1.set_ylim(custom_yrange)

    ax2.plot(Time,NCseeds[:,0,bi2], color = cmapf(bi2), label = Birange[bi2])
    ax2.grid()
    ax2.tick_params(axis='both', labelsize = 28)
    ax2.set_yticks(custom_yticksens)


figloop = figloop + 1

###############################################################
fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )

#plt.title("p0_6 ensemble averages")

nbihalf = nbicouple//2

for bi in range(0,nbihalf):
    bi2 = 2*bi

    PBskp = EMseedsavg[::skip,nshape-1,bi2]
    ax1.plot(tsavgskp,PBskp, color = cmapf(bi2), label = Birange[bi2] )
    ax1.grid()

    ax1.tick_params(axis='both', labelsize = 28)
    ax1.legend(loc='lower left', bbox_to_anchor=(0.99,-1.35), fontsize=22)
    ax1.set_ylim(custom_yrange)

    ax2.plot(Time,NCseeds[:,nshape-1,bi2], color = cmapf(bi2), label = Birange[bi2])
    ax2.grid()

    ax2.tick_params(axis='both', labelsize = 28)



figloop = figloop + 1

##############################

plt.figure(figloop)
plt.grid()
plt.title("Ensemble average walking distance")

skip = 10

for seeds in range(0,nseeds):
        plt.plot(Time,rdistseeds[:,p5p0,seeds,bi4p5], color = cmapf(seeds), label = snum[seeds])

figloop = figloop + 1


#################################################


plt.figure(figloop)
plt.grid()
plt.title("X2 walking distance")

skip = 10

for seeds in range(0,nseeds):
        plt.plot(Time,xdistseedsl2[:,p5p0,seeds,bi4p5], color = cmapf(seeds), label = snum[seeds])

figloop = figloop + 1



################################################

figregime = figloop+2

contour1 = plt.figure(figregime)
bounds = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
contour = plt.contourf(X,Y, zcont, levels=bounds, cmap='rainbow')
cbar = plt.colorbar(contour, location='right')
cbar.ax.tick_params(labelsize=32)
custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

plt.yticks(custom_yticks)
plt.tick_params(axis='both',labelsize=26)
plt.grid()
figregime = figregime+1

###########################################################


plt.figure(figregime)

for b0loop in range(0,nbihalf):
    b02 = 2*b0loop
    plt.plot(pvald,Regimemapstd[:,b0loop], color=cmapf(b02), label = Birange[b02])

    plt.grid()
    plt.tick_params(axis='both', labelsize = 28)

plt.legend(loc='lower left', bbox_to_anchor=(0.99,-0.02), fontsize=20)


figregime = figregime+1

###########################################################

plt.figure(figregime)

for b0loop in range(0,nbihalf):
    b02 = 2*b0loop

    plt.plot(pvald,edgerisetimestd[:,b0loop], color=cmapf(b02), label = Birange[b02])

    plt.grid()
    plt.tick_params(axis='both', labelsize = 28)

plt.legend(loc='lower left', bbox_to_anchor=(0.99,-0.02), fontsize=20)


figregime = figregime+1



###############################################

plt.figure(figregime)

for b0loop in range(0,nbihalf):
    b02 = 2*b0loop
#    plt.scatter(pvald,Regimemapstd[:,b0loop], color=cmapf(b02), label = Birange[b02])
    plt.plot(pvald,ravgdistatrisestd[:,b0loop], color=cmapf(b02), label = Birange[b02])

    plt.grid()
    plt.tick_params(axis='both', labelsize = 28)

plt.legend(loc='lower left', bbox_to_anchor=(0.99,-0.02), fontsize=20)


figregime = figregime+1


####################################################

contour2 = plt.figure(figregime)
bounds2 = [0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0]
contour = plt.contourf(X2,Y2, zcontmagnif, levels=bounds2, cmap='rainbow')
cbar = plt.colorbar(contour)
cbar.ax.tick_params(labelsize=32)
custom_yticks=[5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75]

plt.yticks(custom_yticks)
plt.tick_params(axis='both',labelsize=26)
plt.grid()
#

figregime = figregime+1

#############################################################

contour3 = plt.figure(figregime)
bounds3 = [0,200, 400, 600, 800,1000,1200, 1400, 1600, 1800, 2000, 2200]
contour = plt.contourf(X3,Y3, edgerisetimecontour, levels=bounds3, cmap='jet')
cbar = plt.colorbar(contour, location='right')
cbar.ax.tick_params(labelsize=32)
custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

plt.yticks(custom_yticks)
plt.tick_params(axis='both',labelsize=26)
plt.grid()

figregime = figregime+1

##########################################################

contour4 = plt.figure(figregime)
bounds4 = [0,15,30,45,60,75,90,105,120,135,150]
contour = plt.contourf(X3,Y3, rdistcont, levels=bounds4, cmap='plasma')
cbar = plt.colorbar(contour, location='right')
cbar.ax.tick_params(labelsize=32)
custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

plt.yticks(custom_yticks)
plt.tick_params(axis='both',labelsize=26)
plt.grid()


figregime = figregime+1



########################################################


contour5 = plt.figure(figregime)
#bounds5 = [0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,235,250,265,280]
bounds5 = [0,4,8,12,16,20,24,28,32,36,40]
#bounds5 = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40]
contour = plt.contourf(X3,Y3, neichangescontour, levels=bounds5, cmap='terrain')
cbar = plt.colorbar(contour, location='right')
cbar.ax.tick_params(labelsize=32)
custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

plt.yticks(custom_yticks)
plt.tick_params(axis='both',labelsize=26)
plt.grid()


figregime = figregime+1


####################################################
contour6 = plt.figure(figregime)
#bounds6 = [0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
bounds6 = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8, 8.5, 9]
contour = plt.contourf(X3,Y3, rdispcont, levels=bounds6, cmap='plasma')
cbar = plt.colorbar(contour, location='right')
cbar.ax.tick_params(labelsize=32)
custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

plt.yticks(custom_yticks)
plt.tick_params(axis='both',labelsize=26)
plt.grid()


figregime = figregime+1



##################################################################

contour7 = plt.figure(figregime)
bounds7 = [0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16]
contour = plt.contourf(X3,Y3, velcont, levels=bounds7, cmap='plasma')
#cbar = plt.colorbar(contour, location='right')
cbar = plt.colorbar(contour, location='left')
cbar.ax.tick_params(labelsize=32)
custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

plt.yticks(custom_yticks)
plt.tick_params(axis='both',labelsize=26)
plt.grid()


figregime = figregime+1

################################################################
contour8 = plt.figure(figregime)
bounds8 = [0,2,4,6,8,10,12,14,16,18,20]
contour = plt.contourf(X3,Y3, distoverdisp, levels=bounds8, cmap='plasma')
#cbar = plt.colorbar(contour, location='right')
cbar = plt.colorbar(contour, location='left')
cbar.ax.tick_params(labelsize=32)
custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

plt.yticks(custom_yticks)
plt.tick_params(axis='both',labelsize=26)
plt.grid()


figregime = figregime+1






#################################################################
plt.show()

