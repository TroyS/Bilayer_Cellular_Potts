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


updownseeds = ["seed2updown", "seed3updown", "seed4updown", "seed5updown", "seed6updown"]

nCFC = len(CFCseeds)


nBend = len(Bendseeds)

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


#bendlen = [0, 0, 0]
#for nu in range(0,len(pnum)):
#    bendlen[nu] = pnum[nu]/float(math.pi) 

#print(bendlen)


bendlen = [ [9200, 18400, 27800, 37000, 40000], [10200, 20600, 30800, 40000], [10800, 21800, 32600, 40000] ]


updownlen = [8400.0, 10200.0, 12600.0, 15400.0, 16800.0, 20500.0, 20600.0, 21000.0, 21100.0, 24200.0, 25200.0, 25300.0, 25700.0, 28400.0, 29500.0, 30300.0, 32600.0, 33800.0, 35400.0, 36900.0, 40000.0]


sampfile4 = CFCtitle[0] + CFCshape[0] + "CFC" + CFCseeds[0] + "Bi_" + CFCBi[0] +"Bendingregimes.dat"

nlines4 = ed.getlinesneigh(sampfile4)

#sampfile5 = CFCtitle[0] + CFCshape[0] + CFCseeds[0] +  "updown" + "Bi_" + CFCBi[0] +"regextensions.txt"

sampfile5 = CFCtitle[0] + CFCshape[0] + CFCseeds[0] +  "updown" + "Bi_" + CFCBi[0] +"Bendingregimes.dat"

nlines5 = ed.getlinesneigh(sampfile5)


#sdat_obj = np.dtype([('name', 'U10'), ('tsdat', object), ('edgedat', object),])


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

edgeregular.tsdat, edgeregular.edgedat, edgeregular.edgediff = ed.getalldata(nlines1, nshape, nbicouple,nseeds, pval, Bi,snum, "edgematch.dat")


#tsdat1, edgedat1 = ed.getalldata(nlines1, nshape, nseeds, nbicouple, pval, snum, Bi, "edgematch.dat")



#####################################################

#Then, the long aspect ratio normal cases

edgelongasp = clas.Edgematchdat

CFCadjust =  ["40cellA1000p0_4p25", "40cellA1000p0_4p75", "40cellA1000p0_5"]

CFCadjust2 =  ["40cellA1000p0_4p25CFC", "40cellA1000p0_4p75CFC", "40cellA1000p0_5CFC"]

edgelongasp.tsdat, edgelongasp.edge, edgelongasp.edgediff = ed.getalldata(nlines1, len(CFCshape), len(CFCBi),nCFC, CFCadjust2,CFCBi,CFCseeds, "edgematch.dat")


#################
#Then, the bending runs


edgeCFC = clas.Edgematchdat

edgeCFC.tsdat, edgeCFC.edge, edgeCFC.edgediff = ed.getalldata(nlines3, len(CFCshape), len(CFCBi),nBend, CFCadjust2,CFCBi,Bendseeds, "edgematch.dat")


bendCFC = clas.bendinfoseeds 

bendCFC.bspeed, bendCFC.bentime, bendCFC.name = ed.getallbendspeedCFC(nlines4,len(CFCshape), nBend, len(CFCBi), CFCadjust2,Bendseeds,CFCBi, "Bendingregimes.dat",leny,bendlen,nxx)   

#print(bendCFC.bspeed[0][0][1][1])


bendCFCavg = clas.bendinfoavg

bendCFCavg.bspeed, bendCFCavg.speedstd = ed.avgallbendspeedCFC(len(CFCshape),nBend,len(CFCBi), bendlen, bendCFC.bspeed)


#print(bendCFCavg.bspeed)
#print(bendCFCavg.speedstd)

bendupdown = clas.bendvalues

bendupdown.regIpnts,bendupdown.name = ed.getallexttime(nlines5,len(CFCshape),nCFC,len(CFCBi),CFCadjust,updownseeds,CFCBi, "Bendingregimes.dat")

print(bendupdown.regIpnts)

figloop = 0

################
#Then the updown runs
#
#for p0loop in range(0,len(CFCshape)):
#    for seedloop in range(0,nCFC):
#        for biloop in range(0,len(CFCBi)):
#
#
#            fsharp = CFCtitle[0] + CFCshape[p0loop] + CFCseeds[seedloop] +  "updown" +  "Bi_" + CFCBi[biloop] +"edgematch.dat"
#
#            updownfile = CFCtitle[0] + CFCshape[p0loop] + CFCseeds[seedloop] +   "updown" + "Bi_" + CFCBi[biloop] + "regextensions.txt"
#
#            tPB4 = ed.getdata(fsharp,nlines3)
#
#            tsdat4[:,p0loop,seedloop,biloop] = tPB4[:,0]
#            edgedat4[:,p0loop,seedloop,biloop] = tPB4[:,1]
#
##            upswitch, downswitch = ed.getbendtimeupdown(updownfile)
#
#            updownriset[p0loop,seedloop,biloop] = ed.risetimeavg(edgedat4[:,p0loop,seedloop,biloop],endlist2,nlines2)
#
##            updownspeed = ed.getbendspeedupdown(bendfile,nlines4,leny,nxx)
#
#
##            updownspeed1[p0loop,seedloop,biloop],updownspeed2[p0loop,seedloop,biloop],updownspeed3[p0loop,seedloop,biloop], updownspeed4[p0loop,seedloop,biloop],updownspeed5[p0loop,seedloop,biloop], updownspeed6[p0loop,seedloop,biloop], updownspeed7[p0loop,seedloop,biloop], updownspeed8[p0loop,seedloop,biloop] = ed.storebendspeeds(updownspeed)
#
#
#
#
#############################################
#print("For p0 = 4.25 and Bi 14, we have")
#print(Regimemapavg[0,1])
#
#
#
#
##Print out the std info
#
#precision = 3
#
#
#
##for p0loop in range(0,nshape):
##    for biloop in range(0,nbicouple):
##        zcont[biloop,p0loop] = Regimemapavg[p0loop,biloop]
##        if biloop > 19:
##            zcontmagnif[biloop-20,p0loop] = Regimemapavg[p0loop,biloop]
#            
#
##        print("P0 is ", end=' ')
##        print(pvald[p0loop], end=' ')
##        print(" , Bi is ", end=' ')
##        print(Birange[biloop], end=' ')
##        print(" , and the standard deviation is ",end=' ' )
##        print(f"{Regimemapstd[p0loop,biloop]:.{precision}f}")
#
#
##########################################################
##Now we plot
#
##nbihalf = nbicouple//2
##for p0loop in range(0,nshape):
#
##    tle2 = tle + ", for " + pval2[p0loop]
#
##    plt.figure(p0loop)
##    plt.grid()
##    plt.xlabel("MC Timestep")
#    #plt.ylabel("Edge_match_ratio")
##    plt.ylabel("$\sum{e_{1}e_{2}} / \sum{Ne}$")
#    #plt.title(" Bilayer Edge match ratio vs Timestep " )
##    plt.title(tle2)
#    #plt.show()
#    #plt.close()
#
##    skip = 10
##    tsavgskp = tsavg[::skip]
##    PBratavgskp = PBratavg[::skip]
##    PBratstdskp = PBratstd[::skip]
#
##    for biloop in range(0,nbihalf):
##        biloop2 = 2*biloop
##        PBratavgskp = EMseedsavg[::skip,p0loop,biloop2]
#        #plt.plot(tsavg[i],PBrat[i], color = cmapf(t), label=inle)
#        #plt.plot(tsavg[:,t],PBratavg[:,t], color = cmapf(t), label = Bi[t])
##        plt.plot(tsavgskp,PBratavgskp, color = cmapf(biloop2), label = Birange[biloop2])
##        plt.loglog(tsavgskp[:,t],PBratavgskp[:,t], color = cmapf(t), label = Bi[t])
#
##        plt.errorbar(tsavgskp[:,t],PBratavgskp[:,t], yerr=PBratstdskp[:,t], fmt="-")
##    plt.tick_params(axis='both', labelsize=28)
##    plt.legend(loc='upper right', bbox_to_anchor=(-0.05,1.1), fontsize=16)
##    plt.tight_layout()
#
#
#
#################
##We plot individual seeds
#
#figloop = 0
#
#
#for p0loop in range(0,nshape):
#
#    tle2 = tle + ", for " + pval2[p0loop]
#
#    plt.figure(figloop+p0loop)
#    plt.grid()
##    plt.xlabel("MC Timestep")
#    #plt.ylabel("Edge_match_ratio")
##    plt.ylabel("$\sum{e_{1}e_{2}} / \sum{Ne}$")
#    #plt.title(" Bilayer Edge match ratio vs Timestep " )
#    plt.title(tle2)
#    #plt.show()
#    #plt.close()
#
#    skip = 10
##    tsavgskp = tsavg[::skip]
##    PBratavgskp = PBratavg[::skip]
##    PBratstdskp = PBratstd[::skip]
#
#    #for biloop in range(0,nbihalf):
#    #    biloop2 = 2*biloop
#    #bi[14] is 8.0, effective value 4
#    biind = 1
#    for seedloop in range(0,nseeds):
#        EMskip1 = edgedat1[::skip,p0loop,seedloop,biind]
#        tsskip1 = tsdat1[::skip,p0loop,seedloop,biind]
##        PBratavgskp = EMseedsavg[::skip,p0loop,biloop2]
#        #plt.plot(tsavg[i],PBrat[i], color = cmapf(t), label=inle)
#        #plt.plot(tsavg[:,t],PBratavg[:,t], color = cmapf(t), label = Bi[t])
#        plt.plot(tsskip1,EMskip1, color = cmapf(seedloop), label = snum[seedloop])
##        plt.plot(tsskip,EMskip, color = 'k', label = snum[seedloop])
##        plt.loglog(tsavgskp[:,t],PBratavgskp[:,t], color = cmapf(t), label = Bi[t])
#
##        plt.errorbar(tsavgskp[:,t],PBratavgskp[:,t], yerr=PBratstdskp[:,t], fmt="-")
#    plt.tick_params(axis='both', labelsize=28)
##    plt.xlim(0,1000)
##    plt.legend(loc='upper right', bbox_to_anchor=(-0.05,1.1), fontsize=16)
##    plt.tight_layout()
#
#figloop = figloop + p0loop + 1
#
###############################################################
##########################################################
#
#############################################################
#
#custom_yticks1=[0.0, 0.25, 0.5, 0.75, 1.0]
#custom_yrange=[0.0,1.0]
#
#
#for p0loop in range(0,nshape):
#
#    plt.figure(figloop+p0loop)
#    plt.grid()
#
#    skip = 10
#
#    biind = 0
#
#    #for bloop in range(0,len(BirangeS)):
#    for bloop in range(0,2):
#
##        EMavg1 = EMseedsavg[::skip,p0loop,bloop]
##        tsskip1 = tsdat1[::skip,p0loop,0,bloop]
#
##        lab1 = "Square, " + r"$\lambda_{B}$" + "  = " + BirangeS[bloop]
#
#        EMavg2 = EMlongavg[::skip,p0loop,bloop]
#        tsskip2 = tsdat2[::skip,p0loop,0,bloop]
#
#        lab2 = "Full" + r"$\lambda_{B}$" +  " = " + BirangeS[bloop]
#
#        EMavg3 = EMCFCavg[::skip,p0loop,bloop]
#        tsskip3 = tsdat3[::skip,p0loop,0,bloop]
#
#        lab3 = "Bending" + r"$\lambda_{B}$" + " = " + BirangeS[bloop]
#
#        EMavg4 = EMupdownavg[::skip,p0loop,bloop]
#        tsskip4 = tsdat4[::skip,p0loop,0,bloop]
#
#        lab4 = "Updown" + r"$\lambda_{B}$" + " = " + BirangeS[bloop]
#
#
#
#        plt.plot(tsskip2,EMavg2, color = cmapf(bloop), linestyle = 'solid', label = lab2 )
#        plt.plot(tsskip3,EMavg3, color = cmapf(bloop), linestyle = 'dashed', label = lab3)
#        plt.plot(tsskip4,EMavg4, color = cmapf(bloop), linestyle = 'dotted', label = lab4)
#
#
#        plt.tick_params(axis='both', labelsize=28)
#        plt.legend(loc='lower left', bbox_to_anchor=(0.4,0.0), fontsize=20, ncol = 2)
#        plt.ylim(custom_yrange)        
#
##    plt.xlim(0,1000)
##    plt.legend(loc='upper right', bbox_to_anchor=(-0.05,1.1), fontsize=16)
##    plt.tight_layout()
#
#        figloop = figloop + p0loop + 1
#
#############################################################


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
        error = bendCFCavg.speedstd[p0loop][b0loop]
        labels3 = [x +(b0loop+1)*bwidth for x in labels]
        plt.plot(labels3,bendCFCavg.bspeed[p0loop][b0loop] )
        plt.bar(labels3,bendCFCavg.bspeed[p0loop][b0loop], yerr = error, capsize=3, width=bwidth,label = Bi[b0loop] )

        plt.grid()
        plt.tick_params(axis='both', labelsize=28)
        plt.xticks([])

#        print(updownrisetavg[p0loop,b0loop])

#        plt.axhline(y=updownrisetavg[p0loop,b0loop], color = bcolors[b0loop] )
##
    plt.legend(loc='lower left', bbox_to_anchor=(0.99,0), fontsize=22)
#    plt.title("Up extension time differences")


figloop = figloop + 1


############################################################
##bwidth = 0.25
##
##multiply = 0
##
##
##for p0loop in range(0,nshape):
##    labels = np.arange(len(timeslabel))
##    labels3 = [x + bwidth for x in labels]
##
##    fig = plt.subplots(figsize=(6,6))
##
##    for b0loop in range(0,len(Bi)):
##
##        xparams = [ bdn1diffavg[p0loop,b0loop], bdn2diffavg[p0loop,b0loop], bdn3diffavg[p0loop,b0loop], bdn4diffavg[p0loop,biloop] ]
##        error = [bdn1diffstd[p0loop,b0loop], bdn2diffstd[p0loop,b0loop], bdn3diffstd[p0loop,b0loop], bdn4diffstd[p0loop,biloop] ]
##        labels3 = [x +(b0loop+1)*bwidth for x in labels]
##        plt.plot(labels3,xparams )
##        plt.bar(labels3,xparams, yerr = error, capsize=3, width=bwidth,label = Bi[b0loop] )
##
##        plt.tick_params(axis='both', labelsize=28)
##        plt.xticks([])
##
##
##    plt.axhline(updownrisetavg[p0,bi])
##
##    plt.legend(loc='lower left', bbox_to_anchor=(0.99,0), fontsize=22)
##    plt.title("Down extension times")
##
##
##figloop = figloop + 1
##
##
############################################################
#
#
##print("Actual times")
##print(btimes1avg)
##print(btimes2avg)
##print(btimes3avg)
#
##print("Time differences")
##print(btimediff1)
##print(btimediff2)
##print(btimediff3)
#
##Let's show the data with a scatter plot
#
#figloop = figloop+1
#
#ExtensionsCFC = [1, 2, 3]
#Extensionsupdown = [1, 2, 3, 4]
#
#markers = ['o', 's', '^']
#
#lambdab = r"$\lambda_{B}$"
#
#plt.figure(figloop)
#
#
#for p0loop in range(0,nshape):
#    for b0loop in range(0,len(Bi)):
#        labe = pval2[p0loop] + ", " + lambdab + " " +  Bi[b0loop]
#        bendCFC = [bspeed1avg[p0loop,b0loop], bspeed2avg[p0loop,b0loop], bspeed3avg[p0loop,b0loop] ]
##        bendCFC = [btimes1avg[p0loop,b0loop], btimes2avg[p0loop,b0loop], btimes3avg[p0loop,b0loop] ]
#
#        plt.scatter(ExtensionsCFC, bendCFC, marker = markers[p0loop], color=cmapf(b0loop), label = labe )
#        plt.plot(ExtensionsCFC, bendCFC, color=cmapf(b0loop))
#
#
#plt.legend(loc='upper right', fontsize=20, ncol = 2)
#plt.tick_params(axis='both', labelsize=28)
#
#
##bbox_to_anchor=(0.0,0.0)
#
#figloop = figloop+1
#
#
###########################################################
#
#bwidth = 0.25
#
#multiply = 0
#
#
#for p0loop in range(0,nshape):
#    labels = np.arange(len(CFClabel))
#    labels3 = [x + bwidth for x in labels]
#
#    fig = plt.subplots(figsize=(6,6))
#
#    for b0loop in range(0,len(Bi)):
#
#        xparams = [ bspeed1avg[p0loop,b0loop], bspeed2avg[p0loop,b0loop], bspeed3avg[p0loop,b0loop] ]
#        error = [bspeed1std[p0loop,b0loop], bspeed2std[p0loop,b0loop], bspeed3std[p0loop,b0loop] ]
#        labels3 = [x +(b0loop+1)*bwidth for x in labels]
#        plt.plot(labels3,xparams )
#        plt.bar(labels3,xparams, yerr = error, capsize=3, width=bwidth,label = Bi[b0loop] )
#        plt.tick_params(axis='both', labelsize=28)
#        plt.xticks([])
#    plt.legend(loc='lower left', bbox_to_anchor=(0.99,0), fontsize=22)
#
#figloop = figloop + 1
#
#
#########################################################



plt.show()

