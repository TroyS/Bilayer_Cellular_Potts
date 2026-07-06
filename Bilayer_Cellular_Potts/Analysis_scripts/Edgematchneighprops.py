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


edgeinfo = clas.Edgematchdat()

edgeinfo.tsdat, edgeinfo.edgedat, edgeinfo.edgediff, edgeinfo.name = ed.getalldata(nlines,nshape,nbicouple,nseeds,pval,Bi,snum,"edgematch.dat") 


edgeinfo.edgerise, edgeinfo.edgeriseavg, edgeinfo.edgeristd = ed.getallrisetime(edgeinfo.edgedat,endlist,nlines,nshape,nbicouple,nseeds)


enavgdat = clas.Edgedatavg()

enavgdat.tsdat, enavgdat.edge, enavgdat.edgestd  = ed.edgeavgoverseed(edgeinfo,nlines,nshape,nbicouple)


#############
#Let's start looking at the frequencies of edge match values
#first, the bins:

Edgemins = clas.Edgefreq()

binlow = 0
binhigh = 1
nbins = 20

binwidth = (binhigh-binlow)/nbins

Edgeminvals = clas.Energymins()

Edgemins.edgebins, Edgemins.edgebinsnum, Edgemins.edgebarrseed, Edgemins.edgebarr, Edgemins.edgeseedsfreq, Edgemins.nummins, Edgemins.highmatch, Edgemins.modematch, Edgemins.fullprob, emv = ed.getallbins(edgeinfo.edgedat,binlow,binhigh,nbins, nshape, nbicouple, nseeds,nlines)

clas.Energymins.getmins(Edgeminvals,emv)

#print(Edgemins.fullprob[7])
#print(Edgemins.modematch[nshape-1][nbicouple//2])

###
#The neighbor change data
neighdat = clas.Neighdata()

neighdat.time,neighdat.neichg1, neighdat.neichg2, neighdat.neichg12, neighdat.neichgsum,neighdat.neichgtotseed, neighdat.neichgtot,neighdat.neichgavg,neighdat.neichgstd, neighdat.name = ed.getallneighdata(nlines2,nshape,nseeds,nbicouple,pval,snum,Bi,"neighchangesfix.dat")


#ed.findhalfneichanges(neighdat,nshape,nbicouple,nseeds,nlines2)

###Let's look for seeds with half neighbor changes

#Goodseeds = []

ckk = 0

####
#The COM Data
comdata = clas.trajdata()

comdata.timedat, comdata.xrawl1, comdata.yrawl1, comdata.xrawl2, comdata.yrawl2, comdata.delxl1, comdata.delyl1, comdata.delxl2, comdata.delyl2, comdata.name = msd.getallCOMdata(nlines2,nshape,nseeds,nbicouple,pval,snum,Bi,"COMviadist.dat",ncells)


#comdata.xdist1,comdata.ydist1,comdata.xdist2,comdata.ydist2 = msd.getalldist(comdata.xrawl1,comdata.yrawl1,comdata.xrawl2,comdata.yrawl2,nlines2,ncells,nshape,nseeds,nbicouple)


distdata = clas.MSDanddistavg()

diststats = clas.Datapoints()

#distdata.time, distdata.xdistseed, distdata.ydistseed, distdata.rdistseed,distdata.xdist, distdata.ydist, distdata.rdist, diststats.findistavg, diststats.findiststd = msd.getalldistnospikeinfo(comdata,nlines2,ncells,Xmin,Xmax,Ymin,Ymax,nshape,nseeds,nbicouple)


#distdata.time, distdata.xmsdseed, distdata.ymsdseed, distdata.xdistseed, distdata.ydistseed, distdata.xmsd, distdata.ymsd, distdata.rmsd, distdata.xdist, distdata.ydist, distdata.rdist, 



###
#Properties at rise time

riseprop = clas.Propsatrise()

afterriseprop = clas.Propsafterrise()

timerise1 = time.time()

riseprop.neiatriseseed, riseprop.neiatrise, riseprop.neiatrisestd, afterriseprop.neiafterriseseed, afterriseprop.neiafterrise, afterriseprop.neiafterrisestd = ed.getallneiatandafterrise(neighdat,edgeinfo,nshape,nbicouple,nseeds)

pchk2 = 0
bchk2 = 8
schk2 = 5

print("Rise time, ", riseprop.neiatriseseed[pchk2,bchk2,schk2])
rtimeseed = int(riseprop.neiatriseseed[pchk2,bchk2,schk2])
print("Energymins, ", Edgeminvals.engymin[pchk2][bchk2][schk2])
print("Neichanges, ", sum(neighdat.neichg12[:rtimeseed,pchk2,bchk2,schk2]))


Edgemins.neifreqseed, Edgemins.neifreq = ed.neichangeatjumpprob(edgeinfo.edgedat, neighdat.neichg12, Edgeminvals.engymin, nlines,nshape,nseeds,nbicouple)




#riseprop.distatriseseed,riseprop.speedatriseseed,riseprop.distatrise,riseprop.speedatrise,riseprop.distatrisestd,riseprop.speedatrisestd, speedalt = msd.getalldistatrise(distdata,edgeinfo,nshape,nbicouple,nseeds)


#risetraj = clas.trajdatarise

#risetraj.xmsd1, risetraj.ymsd1, risetraj.xmsd2, risetraj.ymsd2, risetraj.rmsd, risetraj.timedat = msd.getallMSDlagsatrise(comdata.xrawl1,comdata.yrawl1,comdata.xrawl2,comdata.yrawl2,comdata,ncells,nshape,nseeds,nbicouple,edgeinfo)


#risetraj.xmsd1seed, risetraj.ymsd1seed, risetraj.xmsd2seed, risetraj.ymsd2seed, risetraj.rmsdseed,riseprop.msdatrise,riseprop.msdatrisestd,riseprop.fitcoeff= msd.getallmsdatrise(comdata,ncells,Xmin,Xmax,Ymin,Ymax,nshape,nseeds,nbicouple,risetraj)


#print(riseprop.fitcoeff)

timerise2 = time.time()

print("Time to get rise", timerise2-timerise1)

#######
#Properties after rise
#matchprop = clas.Propsafterrise


#matchprop.distafterriseseed,matchprop.speedafterriseseed,riseprop.distafterrise,riseprop.speedafterrise,riseprop.distafterrisestd,riseprop.speedafterrisestd, speedalt = msd.getalldistafterrise(distdata,edgeinfo,nshape,nbicouple,nseeds)

#matchtraj = clas.trajafterrise

timeset1 = time.time()

#matchtraj.xmsd1, matchtraj.ymsd1, matchtraj.xmsd2, matchtraj.ymsd2, matchtraj.rmsd, matchtraj.timedat = msd.getallMSDlagsafterrise(comdata.xrawl1,comdata.yrawl1,comdata.xrawl2,comdata.yrawl2,comdata,ncells,nshape,nseeds,nbicouple,edgeinfo)

timeset2 = time.time()

#print("Time to set", timeset2-timeset1)

timemsd1 = time.time()


#matchtraj.xmsd1seed, matchtraj.ymsd1seed, matchtraj.xmsd2seed, matchtraj.ymsd2seed, matchtraj.rmsdseed,matchprop.msdafterrise,matchprop.msdafterrisestd, matchprop.fitcoeff= msd.getallmsdatrise(comdata,ncells,Xmin,Xmax,Ymin,Ymax,nshape,nseeds,nbicouple,matchtraj)

timemsd2 = time.time()

print("Time for the rest of the msd", timemsd2-timemsd1)


#riseprop.

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
bi4p75 = 17 #Bi = 4.75
bi4 = 14 #Bi = 4
p04p5 = 4
schk = 2
p4p75 = 6
p5p0 = 7
p5p5 = 8
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

######################

figloop = ed.plotsingleedgeandneighchange(edgeinfo,neighdat,skip,pchk2,bchk2,schk2,custom_yticks1,custom_yrange,custom_yticksens,snum,figloop)

###################

figloop = ed.plotedgeandneighchange(edgeinfo,neighdat,skip,pchk,bchk,nseeds,custom_yticks1,custom_yrange,custom_yticksens,snum,figloop)


########################################################
#The edge match difference


figloop = ed.plotedgediffandnei(edgeinfo,neighdat,skip,pchk,bchk,nseeds,custom_yticks1,custom_yrange,custom_yticksens,snum,figloop)


#############################################################


figloop = ed.plotedgeandneighchange(edgeinfo,neighdat,skip,p5p0,bchk,nseeds,custom_yticks1,custom_yrange,custom_yticksens,snum,figloop)



###########################################################

#p0=5, Bi=4.75



###################

figloop = ed.plotedgeandneighchange(edgeinfo,neighdat,skip,nshape-1,bchk,nseeds,custom_yticks1,custom_yrange,custom_yticksens,snum,figloop)


#########################################



figloop =ed.plotensmblavg(enavgdat,neighdat,skip,0,nbicouple,custom_yticks1,custom_yrange,custom_yticksens,Birange,figloop) 


###############################################################

figloop =ed.plotensmblavg(enavgdat,neighdat,skip,nshape-1,nbicouple,custom_yticks1,custom_yrange,custom_yticksens,Birange,figloop)


##############################


figloop+=1

#figloop = msd.plotrdistseeds(distdata,p5p5,bi4p5,nseeds,figloop,"Ensemble average walking distance",snum)


#################################################

#figloop = msd.plotdistseed(distdata,p5p5,bi4p5,nseeds,figloop,"X2 walking distance",snum, usex)


################################################

#figloop = msd.plotrdistseedsatrise(distdata,edgeinfo,p5p5,bi4p5,nseeds,figloop,"Ensemble average walking distance",snum)


##############################################


#figloop = msd.plotmsdratrise(risetraj,p5p5,bi4p5,nseeds,figloop,"MSD at rise",snum, usex)

################################################

#figloop = msd.plotmsdratrise(matchtraj,p5p5,bi4p5,nseeds,figloop,"MSD after matching",snum, usex)



####################################################

custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

bounds = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]


figloop = ed.plotregimemap(edgeinfo.edgedatavg,xvals,yvals,bounds,custom_yticks,1,figloop)

###########################################################


figloop = ed.plothalfthelines(edgeinfo.edgedatstd,pvald,nbicouple,Birange,figloop)


###########################################################


figloop = ed.plothalfthelines(edgeinfo.edgeristd,pvald,nbicouple,Birange,figloop)



###############################################

#figloop = ed.plothalfthelines(riseprop.distatrisestd,pvald,nbicouple,Birange,figloop)



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

##########################################################

bounds4 = [0,15,30,45,60,75,90,105,120,135,150]
custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

#figloop = ed.plotregimemap(riseprop.distatrise,xvals,yvals,bounds4,custom_yticks,3,figloop)


############################################################


########################################################


bounds5 = [0,4,8,12,16,20,24,28,32,36,40]
custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]


#print(riseprop.neiatrise)

figloop = ed.plotregimemap(riseprop.neiatrise,xvals,yvals,bounds5,custom_yticks,4,figloop)


####################################################

#bounds6 = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8, 8.5, 9]
#custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

#figloop = ed.plotregimemap(riseprop.distatrise,xvals,yvals,bounds6,custom_yticks,3,figloop)


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
#bounds7 = [0,0.15,0.3,0.45,0.6,0.75,0.9,1.15,1.3,1.45,1.6]

#print(riseprop.speedatrise)

#figloop = ed.plotregimemap(riseprop.speedatrise,xvals,yvals,bounds7,custom_yticks,3,figloop)


#figloop = ed.plotregimemap(speedalt,xvals,yvals,bounds7,custom_yticks,3,figloop)


################################################################

custom_yticks=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]

bounds77 = [0,0.015,0.03,0.045,0.06,0.075,0.09,0.105,0.13]
bounds777 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#bounds77 = [0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16]
#bounds77 = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4]

#bounds7 = [0,0.15,0.3,0.45,0.6,0.75,0.9,1.15,1.3,1.45,1.6]

#print(riseprop.speedatrise)

#figloop = ed.plotregimemap(riseprop.msdatrise,xvals,yvals,bounds77,custom_yticks,3,figloop)

#figloop = ed.plotregimemap(matchprop.msdafterrise,xvals,yvals,bounds77,custom_yticks,3,figloop)



##############################################################

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


################################################################

print("Figloop ", figloop, " starts the stats regime map")

bounds9 = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]

bounds10 = [0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.0,4.4,4.8]


figloop = ed.plotregimemap(Edgemins.fullprob,xvals,yvals,bounds9,custom_yticks,2,figloop)

figloop = ed.plotregimemap(Edgemins.modematch,xvals,yvals,bounds9,custom_yticks,1,figloop)

figloop = ed.plotregimemap(Edgemins.nummins,xvals,yvals,bounds10,custom_yticks,4,figloop)

figloop = ed.plotregimemap(Edgemins.highmatch,xvals,yvals,bounds9,custom_yticks,1,figloop)

figloop = ed.plotregimemap(Edgemins.neifreq,xvals,yvals,bounds9,custom_yticks,4,figloop)


figloop = ed.plotregimemap(afterriseprop.neiafterrise,xvals,yvals,bounds10,custom_yticks,4,figloop)

#################################################################
plt.show()

