#This will be used to plot the bilayer edge match ratio and neighbor changes, along with a couple of regime maps

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

cmapf = plt.get_cmap('rainbow',nbicouple)
cmaptemp = plt.get_cmap('tab20',nseeds)
endlist = nlines-math.floor(nlines/10)

endrise = nlines-math.floor(nlines/100)

#############################################################
#Now we get the data

getdatatimestar = time.perf_counter()


edgeinfo = clas.Edgematchdat

edgeinfo.tsdat, edgeinfo.edgedat, edgeinfo.edgediff, edgeinfo.name = ed.getalldata(nlines,nshape,nbicouple,nseeds,pval,Bi,snum,"edgematch.dat") 


edgeinfo.edgerise, edgeinfo.edgeriseavg, edgeinfo.edgeristd = ed.getallrisetime(edgeinfo.edgedat,endlist,nlines,nshape,nbicouple,nseeds)

enavgdat = clas.Edgedatavg

enavgdat.tsdat, enavgdat.edge, enavgdat.edgestd  = ed.edgeavgoverseed(edgeinfo,nlines,nshape,nbicouple)


###
#The neighbor change data
neighdat = clas.Neighdata

neighdat.time,neighdat.neichg1, neighdat.neichg2, neighdat.neichg12, neighdat.neichgsum,neighdat.neichgtotseed, neighdat.neichgtot,neighdat.neichgavg,neighdat.neichgstd = ed.getallneighdata(nlines2,nshape,nseeds,nbicouple,pval,snum,Bi,"neighchangesfix.dat")

###Let's look for seeds with double neighbor changes that lead to jumps

#Goodseeds = []

ckk = 0

for p0 in range(0,nshape):
    for bi in range(0,nbicouple):
        for se in range(0,nseeds):
            for lines in range(0,nlines2-1):
                if neighdat.neichg1[lines,p0,bi,se] > 3 and neighdat.neichg2[lines,p0,bi,se] > 3 and edgeinfo.edgedat[10*lines,p0,bi,se] > 0.85 and abs(edgeinfo.edgedat[-1,p0,bi,se] -edgeinfo.edgedat[10*lines,p0,bi,se]) > 0.05:
#                    Goodseeds.append(edgeinfo.name[p0][bi][se])
                    print("The neighbor changes at time ", 10*lines, " ",neighdat.neichg1[lines,p0,bi,se], " and ", neighdat.neichg2[lines,p0,bi,se], " , the final match ratio is ", edgeinfo.edgedat[-1,p0,bi,se], "and the current ratio is ", edgeinfo.edgedat[10*lines,p0,bi,se] )
                    print("For",edgeinfo.name[p0][bi][se]) 
                    print(" ")


#print(Goodseeds)

####
#The COM Data
comdata = clas.trajdata

comdata.timedat, comdata.xrawl1, comdata.yrawl1, comdata.xrawl2, comdata.yrawl2, comdata.delxl1, comdata.delyl1, comdata.delxl2, comdata.delyl2, comdata.name = msd.getallCOMdata(nlines2,nshape,nseeds,nbicouple,pval,snum,Bi,"COMviadist.dat",ncells)


comdata.xdist1,comdata.ydist1,comdata.xdist2,comdata.ydist2 = msd.getalldist(comdata.xrawl1,comdata.yrawl1,comdata.xrawl2,comdata.yrawl2,nlines2,ncells,nshape,nseeds,nbicouple)


distdata = clas.MSDanddistavg

diststats = clas.Datapoints

distdata.time, distdata.xdistseed, distdata.ydistseed, distdata.rdistseed,distdata.xdist, distdata.ydist, distdata.rdist, diststats.findistavg, diststats.findiststd = msd.getalldistnospikeinfo(comdata,nlines2,ncells,Xmin,Xmax,Ymin,Ymax,nshape,nseeds,nbicouple)

##
#The speed data
speedinfo = clas.Speeddata

speedinfo.distatriseseed, speedinfo.speedatriseseed, speedinfo.distatrise,speedinfo.speedatrise, speedinfo.distatrisestd, speedinfo.speedatrisestd = msd.getalldistatrise(distdata,edgeinfo,nshape,nbicouple,nseeds)


###
#Properties at rise time

riseprop = clas.Propsatrise

riseprop.neiatriseseed, riseprop.neiatrise, riseprop.neiatrisestd = ed.getallneiatrise(neighdat,edgeinfo,nshape,nbicouple,nseeds)

riseprop.distatriseseed,riseprop.speedatriseseed,riseprop.distatrise,riseprop.speedatrise,riseprop.distatrisestd,riseprop.speedatrisestd = msd.getalldistatrise(distdata,edgeinfo,nshape,nbicouple,nseeds)

#######################################################

#print("The max diff  is ", edgediffmax[:,:,:])

#Now we analyze

getdatatimeend = time.perf_counter()
getdatatime = getdatatimestar - getdatatimeend
print(f"Get data time: {getdatatime:0.4f} seconds")

############################################
#Now we get the averages 

getavgtimestar = time.perf_counter()

#First the average over seeds

#Then the average over the last couple of timesteps
edgeinfo.edgedatavg, edgeinfo.edgedatstd = ed.Avgoverendtime(edgeinfo.edgedat,nshape,nbicouple,nseeds,endlist,nlines)

#Print out the std info

precision = 3
    
getavgtimeend = time.perf_counter()
getavgtime = getavgtimestar-getavgtimeend
print(f"Average time : {getavgtime:0.4f} seconds")


#########################################################
#Now we plot

figloop = 0

nbihalf = nbicouple//2

pchk = 0
#bchk = 5 #Bi = 1.75
bchk = 12 #Bi = 3.5 
bi4p5 = 16 #Bi = 4.5
bi4 = 14 #Bi = 4
p04p5 = 4
schk = 2
p5p0 = 7
skip = 10
tsavgskp = edgeinfo.tsdat[::skip,0,0,0]
Ntime = neighdat.time[:,0,0,0]
usex = 0
usey = 1 


################################################################

#fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )
custom_yticks1=[0.0, 0.3, 0.6, 0.9]
custom_yrange=[0.0,1.0]
custom_yticksens=[0,4,8,12]

#####
#figloop = ed.plotneighchange(neighdat,p04p5,bi4,schk,nseeds,custom_yticks1,custom_yrange,custom_yticksens,snum,figloop)
#plt.show()
#####


figloop = ed.plotedgeandneighchange(edgeinfo,neighdat,skip,pchk,bchk,nseeds,custom_yticks1,custom_yrange,custom_yticksens,snum,figloop)


########################################################
#The edge match difference


figloop = ed.plotedgediffandnei(edgeinfo,neighdat,skip,pchk,bchk,nseeds,custom_yticks1,custom_yrange,custom_yticksens,snum,figloop)

#fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )
#custom_yticks1=[0.0, 0.3, 0.6, 0.9]
#custom_yrange=[0.0,1.0]
#custom_yticksens=[0,4,8,12]


#for seedloop in range(0,nseeds): 

    #PBseedskp = edgedat[::skip,pchk,seedloop,bchk]
#    PBseedskp = edgeinfo.edgediff[::skip,pchk,bchk,seedloop]
#    ax1.plot(tsavgskp,PBseedskp, color = cmaptemp(seedloop), label = snum[seedloop] )
#    ax1.grid()
#    ax1.tick_params(axis='both', labelsize = 28)
#    ax1.set_ylim(custom_yrange)

#    ax2.plot(Ntime,neighdat.neichg12[:,pchk,bchk,seedloop], color = cmaptemp(seedloop), label = snum[seedloop])
#    ax2.grid()
#    ax2.tick_params(axis='both', labelsize = 28)


#figloop = figloop + 1


#############################################################


figloop = ed.plotedgeandneighchange(edgeinfo,neighdat,skip,p5p0,bchk,nseeds,custom_yticks1,custom_yrange,custom_yticksens,snum,figloop)


#fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )
#custom_yticks1=[0.0, 0.3, 0.6, 0.9]
#custom_yrange=[0.0,1.0]

#for seedloop in range(0,nseeds): 

#    PBseedskp = edgeinfo.edgedat[::skip,p5p0,bchk,seedloop]
#    ax1.plot(tsavgskp,PBseedskp, color = cmaptemp(seedloop), label = snum[seedloop] )
#    ax1.grid()
#    ax1.tick_params(axis='both', labelsize = 28)
#    ax1.set_ylim(custom_yrange)

#    ax2.plot(Ntime,neighdat.neichg12[:,p5p0,bchk,seedloop], color = cmaptemp(seedloop), label = snum[seedloop])
#    ax2.grid()
#    ax2.tick_params(axis='both', labelsize = 28)


#figloop = figloop + 1



###########################################################

figloop = ed.plotedgeandneighchange(edgeinfo,neighdat,skip,nshape-1,bchk,nseeds,custom_yticks1,custom_yrange,custom_yticksens,snum,figloop)




#fig, (ax1, ax2) = plt.subplots(2,1, figsize=(6,6) )

#for seedloop in range(0,nseeds): 

#    PBseedskp = edgeinfo.edgedat[::skip,nshape-1,bchk,seedloop]
#    ax1.plot(tsavgskp,PBseedskp, color = cmaptemp(seedloop), label = snum[seedloop] )
#    ax1.grid()
#    ax1.tick_params(axis='both', labelsize = 28)
#    ax1.set_ylim(custom_yrange)

#    ax2.plot(Ntime,neighdat.neichg12[:,nshape-1,bchk,seedloop], color = cmaptemp(seedloop), label = snum[seedloop])
#    ax2.grid()
#    ax2.tick_params(axis='both', labelsize = 28)



#figloop = figloop + 1



#########################################



figloop =ed.plotensmblavg(enavgdat,neighdat,skip,0,nbicouple,custom_yticks1,custom_yrange,custom_yticksens,Birange,figloop) 




#fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )

#plt.title("p0_4 ensemble averages")

#nbihalf = nbicouple//2

#for bi in range(0,nbihalf):
#    bi2 = 2*bi

#    PBskp = enavgdat.edge[::skip,0,bi2]
#    ax1.plot(tsavgskp,PBskp, color = cmapf(bi2), label = Birange[bi2] )
#    ax1.grid()

#    ax1.tick_params(axis='both', labelsize = 28)
#    ax1.legend(loc='lower left', bbox_to_anchor=(0.99,-1.35), fontsize=22)
#    ax1.set_ylim(custom_yrange)

#    ax2.plot(Ntime,neighdat.neichgsum[:,0,bi2], color = cmapf(bi2), label = Birange[bi2])
#    ax2.grid()
#    ax2.tick_params(axis='both', labelsize = 28)
#    ax2.set_yticks(custom_yticksens)


#figloop = figloop + 1

###############################################################

figloop =ed.plotensmblavg(enavgdat,neighdat,skip,nshape-1,nbicouple,custom_yticks1,custom_yrange,custom_yticksens,Birange,figloop)



#fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )

#plt.title("p0_6 ensemble averages")

#nbihalf = nbicouple//2

#for bi in range(0,nbihalf):
#    bi2 = 2*bi

#    PBskp = enavgdat.edge[::skip,nshape-1,bi2]
#    ax1.plot(tsavgskp,PBskp, color = cmapf(bi2), label = Birange[bi2] )
#    ax1.grid()

#    ax1.tick_params(axis='both', labelsize = 28)
#    ax1.legend(loc='lower left', bbox_to_anchor=(0.99,-1.35), fontsize=22)
#    ax1.set_ylim(custom_yrange)

#    ax2.plot(Ntime,neighdat.neichgsum[:,nshape-1,bi2], color = cmapf(bi2), label = Birange[bi2])
#    ax2.grid()

#    ax2.tick_params(axis='both', labelsize = 28)



#figloop = figloop + 1


##############################


figloop+=1

figloop = msd.plotrdistseeds(distdata,p5p0,bi4p5,nseeds,figloop,"Ensemble average walking distance",snum)

#plt.figure(figloop)
#plt.grid()
#plt.title("Ensemble average walking distance")

#skip = 10

#for seeds in range(0,nseeds):
#        plt.plot(Time,rdistseeds[:,p5p0,seeds,bi4p5], color = cmapf(seeds), label = snum[seeds])

#figloop = figloop + 1


#################################################

figloop = msd.plotdistseed(distdata,p5p0,bi4p5,nseeds,figloop,"X2 walking distance",snum, usex)



#plt.figure(figloop)
#plt.grid()
#plt.title("X2 walking distance")

#skip = 10

#for seeds in range(0,nseeds):
#        plt.plot(Time,xdistseedsl2[:,p5p0,seeds,bi4p5], color = cmapf(seeds), label = snum[seeds])

#figloop = figloop + 1


################################################


custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

bounds = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]


figloop = ed.plotregimemap(edgeinfo.edgedatavg,xvals,yvals,bounds,custom_yticks,1,figloop)

#contour1 = plt.figure(figregime)
#bounds = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
#contour = plt.contourf(X,Y, zcont, levels=bounds, cmap='rainbow')
#cbar = plt.colorbar(contour, location='right')
#cbar.ax.tick_params(labelsize=32)
#custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

#plt.yticks(custom_yticks)
#plt.tick_params(axis='both',labelsize=26)
#plt.grid()
#figregime = figregime+1


#plt.show()

###########################################################


figloop = ed.plothalfthelines(edgeinfo.edgedatstd,pvald,nbicouple,Birange,figloop)


#plt.figure(figregime)

#for b0loop in range(0,nbihalf):
#    b02 = 2*b0loop
#    plt.plot(pvald,Regimemapstd[:,b0loop], color=cmapf(b02), label = Birange[b02])

#    plt.grid()
#    plt.tick_params(axis='both', labelsize = 28)

#plt.legend(loc='lower left', bbox_to_anchor=(0.99,-0.02), fontsize=20)


#figregime = figregime+1

###########################################################


figloop = ed.plothalfthelines(edgeinfo.edgeristd,pvald,nbicouple,Birange,figloop)


#plt.figure(figregime)

#for b0loop in range(0,nbihalf):
#    b02 = 2*b0loop

#    plt.plot(pvald,edgerisetimestd[:,b0loop], color=cmapf(b02), label = Birange[b02])

#    plt.grid()
#    plt.tick_params(axis='both', labelsize = 28)

#plt.legend(loc='lower left', bbox_to_anchor=(0.99,-0.02), fontsize=20)


#figregime = figregime+1



###############################################

figloop = ed.plothalfthelines(riseprop.distatrisestd,pvald,nbicouple,Birange,figloop)


#plt.figure(figregime)

#for b0loop in range(0,nbihalf):
#    b02 = 2*b0loop
#    plt.scatter(pvald,Regimemapstd[:,b0loop], color=cmapf(b02), label = Birange[b02])
#    plt.plot(pvald,ravgdistatrisestd[:,b0loop], color=cmapf(b02), label = Birange[b02])

#    plt.grid()
#    plt.tick_params(axis='both', labelsize = 28)

#plt.legend(loc='lower left', bbox_to_anchor=(0.99,-0.02), fontsize=20)


#figregime = figregime+1


####################################################

#bounds2 = [0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0]
#custom_yticks=[5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75]

#contour2 = plt.figure(figregime)
#bounds2 = [0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0]
#contour = plt.contourf(X2,Y2, zcontmagnif, levels=bounds2, cmap='rainbow')
#cbar = plt.colorbar(contour)
#cbar.ax.tick_params(labelsize=32)
#custom_yticks=[5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75]


#figloop = ed.plotregimemap(enavgdat.edge,X,Y,bounds,custom_yticks,1,figloop)



#plt.yticks(custom_yticks)
#plt.tick_params(axis='both',labelsize=26)
#plt.grid()
#

#figregime = figregime+1

#############################################################

bounds3 = [0,200, 400, 600, 800,1000,1200, 1400, 1600, 1800, 2000, 2200]
custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

figloop = ed.plotregimemap(edgeinfo.edgeriseavg,xvals,yvals,bounds3,custom_yticks,2,figloop)

#contour3 = plt.figure(figregime)
#bounds3 = [0,200, 400, 600, 800,1000,1200, 1400, 1600, 1800, 2000, 2200]
#contour = plt.contourf(X3,Y3, edgerisetimecontour, levels=bounds3, cmap='jet')
#cbar = plt.colorbar(contour, location='right')
#cbar.ax.tick_params(labelsize=32)
#custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]


#plt.yticks(custom_yticks)
#plt.tick_params(axis='both',labelsize=26)
#plt.grid()

#figregime = figregime+1

##########################################################

bounds4 = [0,15,30,45,60,75,90,105,120,135,150]
custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

figloop = ed.plotregimemap(diststats.findistavg,xvals,yvals,bounds4,custom_yticks,3,figloop)


#contour4 = plt.figure(figregime)
#bounds4 = [0,15,30,45,60,75,90,105,120,135,150]
#contour = plt.contourf(X3,Y3, rdistcont, levels=bounds4, cmap='plasma')
#cbar = plt.colorbar(contour, location='right')
#cbar.ax.tick_params(labelsize=32)
#custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

#plt.yticks(custom_yticks)
#plt.tick_params(axis='both',labelsize=26)
#plt.grid()


#figregime = figregime+1



########################################################


bounds5 = [0,4,8,12,16,20,24,28,32,36,40]
custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]


figloop = ed.plotregimemap(riseprop.neiatrise,xvals,yvals,bounds5,custom_yticks,1,figloop)


#contour5 = plt.figure(figregime)
#bounds5 = [0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,235,250,265,280]
#bounds5 = [0,4,8,12,16,20,24,28,32,36,40]
#bounds5 = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40]
#contour = plt.contourf(X3,Y3, neichangescontour, levels=bounds5, cmap='terrain')
#cbar = plt.colorbar(contour, location='right')
#cbar.ax.tick_params(labelsize=32)
#custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

#plt.yticks(custom_yticks)
#plt.tick_params(axis='both',labelsize=26)
#plt.grid()


#figregime = figregime+1

plt.show()


####################################################

bounds6 = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8, 8.5, 9]
custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

figloop = ed.plotregimemap(riseprop.distatrise,xvals,yvals,bounds5,custom_yticks,3,figloop)


#contour6 = plt.figure(figregime)
#bounds6 = [0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
#bounds6 = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8, 8.5, 9]
#contour = plt.contourf(X3,Y3, rdispcont, levels=bounds6, cmap='plasma')
#cbar = plt.colorbar(contour, location='right')
#cbar.ax.tick_params(labelsize=32)
#custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

#plt.yticks(custom_yticks)
#plt.tick_params(axis='both',labelsize=26)
#plt.grid()


#figregime = figregime+1



##################################################################


custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

bounds7 = [0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16]

figloop = ed.plotregimemap(riseprop.speedatrise,xvals,yvals,bounds7,custom_yticks,3,figloop)

#contour7 = plt.figure(figregime)
#bounds7 = [0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16]
#contour = plt.contourf(X3,Y3, velcont, levels=bounds7, cmap='plasma')
#cbar = plt.colorbar(contour, location='right')
#cbar = plt.colorbar(contour, location='left')
#cbar.ax.tick_params(labelsize=32)
#custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

#plt.yticks(custom_yticks)
#plt.tick_params(axis='both',labelsize=26)
#plt.grid()


#figregime = figregime+1

################################################################


###Maybe we'll put this in


#contour8 = plt.figure(figregime)
#bounds8 = [0,2,4,6,8,10,12,14,16,18,20]
#contour = plt.contourf(X3,Y3, distoverdisp, levels=bounds8, cmap='plasma')
#cbar = plt.colorbar(contour, location='right')
#cbar = plt.colorbar(contour, location='left')
#cbar.ax.tick_params(labelsize=32)
#custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

#plt.yticks(custom_yticks)
#plt.tick_params(axis='both',labelsize=26)
#plt.grid()


#figregime = figregime+1






#################################################################
plt.show()

