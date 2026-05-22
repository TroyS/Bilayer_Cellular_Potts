#We analyze the COM and shape index data

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import math
import statistics

import Edgedatfunctions as ed
import MSDfunctions as msd
import plotly.graph_objects as go
import time

import BCPM_Classes as clas

from matplotlib.colors import TABLEAU_COLORS, same_color
from matplotlib.pyplot import cm
from matplotlib import cm

pval = ["p0_3p5", "p0_3p75", "p0_4", "p0_4p12", "p0_4p25", "p0_4p37","p0_4p5", "p0_4p62","p0_4p75", "p0_5p0", "p0_5p5", "p0_6"]

pval2 =['$p_{0}=3.5', '$p_{0}=3.75', '$p_{0}=4$', '$p_{0}=4.12$', '$p_{0}=4.25$','$p_{0}=4.37$','$p_{0}=4.5$', '$p_{0}=4.62$', '$p_{0}=4.75$', '$p_{0}=5$', '$p_{0}=5.5$', '$p_{0}=6.0$']


pvald = [3.5, 3.75, 4.0, 4.12, 4.25, 4.37, 4.5, 4.62, 4.75, 5.0, 5.5, 6.0]

nshape = len(pvald)

Bi = ["0.0"]

nbicouple = len(Bi)

Birange = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75]

vlen = len(pval)
blen = len(Birange)

Binum = [0]*(vlen*blen)
colorssc = [0]*(vlen*blen)
colorsm = [0]*(vlen*blen)

colorwhl = ["red", "blue", "green" ,"purple", "orange", "pink", "brown", "gray", "olive", "cyan", "gold", "navy", "black"]

for i in range(0,blen):
    for j in range(0,vlen):
        colorssc[i+blen*j] = 10*(vlen+1)
        colorsm[i+blen*j] = colorwhl[vlen]
        Binum[i+blen*j] = Birange[i]


nlines = 0

snum = ["seed50", "seed51", "seed52", "seed53", "seed54", "seed60", "seed61", "seed62", "seed63", "seed64"]

nseeds = len(snum)

method = ["viadist"]

#We get the data

sampfile1 = pval[0] + snum[0] + "Bi_" + Bi[0] + "region.dat"

npoints, nxx, latti, ncells = msd.getinforegion(sampfile1)

sampfile2 = pval[0]+snum[0]+"Bi_"+Bi[0]+"COM"+method[0]+".dat"

nlines = msd.getlinesCOM(sampfile2)

sampfile3 = pval[0] + snum[0] + "Bi_" + Bi[0] + "neighchanges.dat"
nlinesnei = ed.getlinesneigh(sampfile3)

sampfile4 = pval[0] + snum[0] + "Bi_" + Bi[0] + "shapeindx.dat"
nlinesshp = ed.getlinesneigh(sampfile4)


rad3 = math.sqrt(3)
Lat = 1
if latti == 2:
    Lat = math.sqrt(2/rad3)

ny = npoints/nxx
Xmin = Lat
Xmax = (nxx+0.5)*Lat
Ymin = rad3*Lat/2
Ymax = rad3*Lat*ny/2

sidelen = Lat/rad3

cmapf = plt.get_cmap('rainbow',ncells)
 
starttime = 0
endtime = math.floor(2*nlines/5)

endshape = nlines-math.ceil(nlines/10)

plen = len(pval)
bilen = len(Bi)


############################

comseeds = clas.trajdata

comseeds.timedat, comseeds.xrawl1, comseeds.yrawl1, comseeds.xrawl2, comseeds.yrawl2, comseeds.delxl1, comseeds.delyl1, comseeds.delxl2, comseeds.delyl2, comseeds.name = msd.getallCOMdata(nlines,plen,nseeds,bilen,pval,snum,Bi,"COMviadist.dat",ncells)


#To check, we'll print

pchk = 3
schk = 6
bchk = 0

figloop = 0



timemsd1 = time.time()

comseeds.xmsd1, comseeds.ymsd1, comseeds.xmsd2, comseeds.ymsd2 = msd.getallMSDlags(comseeds.xrawl1,comseeds.yrawl1,comseeds.xrawl2,comseeds.yrawl2,nlines,ncells,plen,nseeds,bilen)


timemsd2 = time.time()
print("Time to get MSDs ", timemsd2-timemsd1)

comseeds.xdist1,comseeds.ydist1,comseeds.xdist2,comseeds.ydist2 = msd.getalldist(comseeds.xrawl1,comseeds.yrawl1,comseeds.xrawl2,comseeds.yrawl2,nlines,ncells,plen,nseeds,bilen)


neiseeds = clas.Neighdata

neiseeds.time, neiseeds.neichg1, neiseeds.neichg2, neiseeds.neichg12,neiseeds.neichgsum,neiseeds.neichgtotseed, neiseeds.neichgtot,neiseeds.neichgavg,neiseeds.neichgstd = ed.getallneighdata(nlines,plen,nseeds,bilen,pval,snum,Bi,"neighchanges.dat")


Shapeseeds = clas.Shapedata

Shapeseeds.Pl1, Shapeseeds.Pl2,Shapeseeds.p0l1,Shapeseeds.p0l2, Shapeseeds.Pavg, Shapeseeds.Pstd, Shapeseeds.p0avg, Shapeseeds.p0std = msd.getallshapedata(nlines,ncells,plen,nseeds,bilen,pval,snum,Bi,"shapeindx.dat")


mobility = clas.MSDanddistavg

Dpoints = clas.Datapoints


timeavg1 = time.time()

mobility.time, mobility.xmsdseed, mobility.ymsdseed, mobility.xdistseed, mobility.ydistseed,mobility.xmsd, mobility.ymsd, mobility.rmsd, mobility.xdist, mobility.ydist, mobility.rdist, Dpoints.diffcoeffavg, Dpoints.diffcoeffstd, Dpoints.findistavg, Dpoints.findiststd = msd.getallmsdanddist(comseeds,nlines,ncells,Xmin,Xmax,Ymin,Ymax,plen,nseeds,bilen,starttime,endtime)


timeavg2 = time.time()
print("Time to get averages ", timeavg2-timeavg1)




###############################################


figloop = 0
plt.figure(figloop)
plt.grid()
plt.title("Neighbor changes")
plt.scatter(pvald,neiseeds.neichgavg[:,0])
plt.plot(pvald,neiseeds.neichgavg[:,0])
figloop = figloop + 1



#######################################

#To check, we'll print

pchk = 3
schk = 6

figloop = 0

usex = 0
usey = 1

figloop = msd.plotCOM(comseeds,ncells,pchk,0,schk,figloop,"Original X COMS",usex)


   


#########################################################

#############################################
#X walking distance

figloop = msd.plotdist(comseeds,ncells,pchk,0,schk,figloop,"|X COM|",usex)



##########################################
#Ensemble averaged X walking distance


figloop = msd.plotenavgdist(mobility,plen,0,figloop,"<|X COM |>", pval,usex)

##############################################


figloop = msd.plotenavgmsd(mobility,plen,0,figloop,"<|X^{2} |>", pval,usey)


###########################################


figloop = msd.plotCOM(comseeds,ncells,pchk,0,schk,figloop,"Original Y COMS",usey)

   


#######################################################
#Y walking distance

figloop = msd.plotdist(comseeds,ncells,pchk,0,schk,figloop,"|Y COM|",usey)

###############################################


figloop = msd.plotenavgdist(mobility,plen,0,figloop,"<|Y COM |>", pval,usey)



######################################################

#Mean squared displacement Y


figloop = msd.plotenavgmsd(mobility,plen,0,figloop,"<|Y^{2} |>", pval,usey)


################################################################

#Mean squared displacement R

figloop = msd.plotmsdr(mobility,plen,0,figloop,"<|R^{2} |>", pval,endtime)



########################################################


#Now, the diffusion coefficient and observed shape index on the same plot

figloop = msd.plotshapeanddiffcoeff(Shapeseeds,Dpoints,pvald,0,figloop)


#######################################################################

figloop = msd.plotneighanddistfin(neiseeds,Dpoints,pvald,0,figloop)



plt.show()

