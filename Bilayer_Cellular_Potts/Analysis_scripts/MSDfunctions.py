#This will contain functions for analyzing the MSD

import numpy as ny
import math
import statistics
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import plotly.graph_objects as go
from matplotlib.colors import TABLEAU_COLORS, same_color
from matplotlib.pyplot import cm
from matplotlib import cm

######################################################

#First we get the number of lines in a file

def getlinesCOM(fle):
    with open(fle,'r') as cfile:
        clines = cfile.readlines()
    nlines = len(clines)-1

    return nlines


#########################################################

#We can get data from the region file

def getinforegion(fle):
    with open(fle,'r') as rfile:
        rlines = rfile.readlines()
        rdata = rlines[1].split()
        ns = int(rdata[0])
        nx = int(rdata[1])
        lat = int(rdata[2])
        nc = int(rdata[3])

    return ns, nx, lat, nc



#########################################################
#We get the data from the files

def getposfromdat(fle,n):

    cellnum1 = ny.zeros((n), dtype='int')
    cellnum2 = ny.zeros((n), dtype='int')
    xlay1 = ny.zeros((n))
    ylay1 = ny.zeros((n))
    xlay2 = ny.zeros((n))
    ylay2 = ny.zeros((n))

    with open(fle,'r') as Cfile:
        Clines = Cfile.readlines()
        for lines in range(1,2*n+1):
            Cdata = Clines[lines].split();
            if lines <= n:
                cellnum1[lines-1] = int(Cdata[1])
                xlay1[lines-1] = float(Cdata[2])
                ylay1[lines-1] = float(Cdata[3])
            if lines > n:
                cellnum2[lines-n-1] = int(Cdata[1])
                xlay2[lines-n-1] = float(Cdata[2])
                ylay2[lines-n-1] = float(Cdata[3])


    return cellnum1, cellnum2, xlay1, ylay1, xlay2, ylay2 

##################################################################

#Then we get the COM data

def getCOMdata(fle,n,cnum):

    time = ny.zeros(n)
    xlay1 = ny.zeros((n,cnum))
    ylay1 = ny.zeros((n,cnum))
    xlay2 = ny.zeros((n,cnum))
    ylay2 = ny.zeros((n,cnum))

    with open(fle,'r') as Cfile:
        Clines = Cfile.readlines()
        for lines in range(1,n+1):
            Cdata = Clines[lines].split();
            time[lines-1] = float(Cdata[0])
            for terms in range(0,cnum):
                index1 = 2*terms+1
                index2 = 2*terms+2
                index3 = 2*terms+1+2*cnum
                index4 = 2*terms+1+2*cnum+1
                xlay1[lines-1,terms] = float(Cdata[index1]) 
                ylay1[lines-1,terms] = float(Cdata[index2])
                xlay2[lines-1,terms] = float(Cdata[index3])
                ylay2[lines-1,terms] = float(Cdata[index4])


    return time, xlay1, ylay1, xlay2, ylay2


######################################################
def getallCOMdata(n,np0,nse,nbi,namep0,namesee,namebi,name,cnum):



    time = ny.zeros((n,np0,nbi,nse))
    xlay1 = ny.zeros((n,cnum,np0,nbi,nse))
    ylay1 = ny.zeros((n,cnum,np0,nbi,nse))
    xlay2 = ny.zeros((n,cnum,np0,nbi,nse))
    ylay2 = ny.zeros((n,cnum,np0,nbi,nse))
    xdiffl1 = ny.zeros((n-1,cnum,np0,nbi,nse))
    ydiffl1 = ny.zeros((n-1,cnum,np0,nbi,nse))
    xdiffl2 = ny.zeros((n-1,cnum,np0,nbi,nse))
    ydiffl2 = ny.zeros((n-1,cnum,np0,nbi,nse))

    Namu=[[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]

    for p0 in range(0,np0):
        for bi in range(0,nbi):
            for se in range(0,nse):

                filenme = namep0[p0] + namesee[se] + "Bi_" + namebi[bi] + name

                Namu[p0][bi][se] = filenme

                with open(filenme,'r') as Cfile:
                    Clines = Cfile.readlines()
                    for lines in range(1,n+1):
                        Cdata = Clines[lines].split();
                        time[lines-1,p0,bi,se] = float(Cdata[0])
                        for terms in range(0,cnum):
                            index1 = 2*terms+1
                            index2 = 2*terms+2
                            index3 = 2*terms+1+2*cnum
                            index4 = 2*terms+1+2*cnum+1
                            xlay1[lines-1,terms,p0,bi,se] = float(Cdata[index1])
                            ylay1[lines-1,terms,p0,bi,se] = float(Cdata[index2])
                            xlay2[lines-1,terms,p0,bi,se] = float(Cdata[index3])
                            ylay2[lines-1,terms,p0,bi,se] = float(Cdata[index4])
                
                            
                for terms2 in range(0,cnum):
                    xdiffl1[:,terms2,p0,bi,se] = ny.diff(xlay1[:,terms2,p0,bi,se])
                    ydiffl1[:,terms2,p0,bi,se] = ny.diff(ylay1[:,terms2,p0,bi,se])
                    xdiffl2[:,terms2,p0,bi,se] = ny.diff(xlay2[:,terms2,p0,bi,se])
                    ydiffl2[:,terms2,p0,bi,se] = ny.diff(ylay2[:,terms2,p0,bi,se])



    return time, xlay1, ylay1, xlay2, ylay2, xdiffl1, ydiffl1, xdiffl2, ydiffl2, Namu


########################################################
def getallCOMatrise(n,np0,nse,nbi,namep0,namesee,namebi,name,cnum):



    time = ny.zeros((n,np0,nbi,nse))
#    xlay1 = ny.zeros((n,cnum,np0,nbi,nse))
#    ylay1 = ny.zeros((n,cnum,np0,nbi,nse))
#    xlay2 = ny.zeros((n,cnum,np0,nbi,nse))
#    ylay2 = ny.zeros((n,cnum,np0,nbi,nse))
#    xdiffl1 = ny.zeros((n-1,cnum,np0,nbi,nse))
#    ydiffl1 = ny.zeros((n-1,cnum,np0,nbi,nse))
#    xdiffl2 = ny.zeros((n-1,cnum,np0,nbi,nse))
#    ydiffl2 = ny.zeros((n-1,cnum,np0,nbi,nse))

    xlay1 = ny.zeros((n,cnum,np0,nbi,nse))
    ylay1 = ny.zeros((n,cnum,np0,nbi,nse))
    xlay2 = ny.zeros((n,cnum,np0,nbi,nse))
    ylay2 = ny.zeros((n,cnum,np0,nbi,nse))
    xdiffl1 = ny.zeros((n-1,cnum,np0,nbi,nse))
    ydiffl1 = ny.zeros((n-1,cnum,np0,nbi,nse))
    xdiffl2 = ny.zeros((n-1,cnum,np0,nbi,nse))
    ydiffl2 = ny.zeros((n-1,cnum,np0,nbi,nse))




    Namu=[[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]

    for p0 in range(0,np0):
        for bi in range(0,nbi):
            for se in range(0,nse):

                filenme = namep0[p0] + namesee[se] + "Bi_" + namebi[bi] + name

                Namu[p0][bi][se] = filenme

                with open(filenme,'r') as Cfile:
                    Clines = Cfile.readlines()
                    for lines in range(1,n+1):
                        Cdata = Clines[lines].split();
                        time[lines-1,p0,bi,se] = float(Cdata[0])
                        for terms in range(0,cnum):
                            index1 = 2*terms+1
                            index2 = 2*terms+2
                            index3 = 2*terms+1+2*cnum
                            index4 = 2*terms+1+2*cnum+1
                            xlay1[lines-1,terms,p0,bi,se] = float(Cdata[index1])
                            ylay1[lines-1,terms,p0,bi,se] = float(Cdata[index2])
                            xlay2[lines-1,terms,p0,bi,se] = float(Cdata[index3])
                            ylay2[lines-1,terms,p0,bi,se] = float(Cdata[index4])
                
                            
                for terms2 in range(0,cnum):
                    xdiffl1[:,terms2,p0,bi,se] = ny.diff(xlay1[:,terms2,p0,bi,se])
                    ydiffl1[:,terms2,p0,bi,se] = ny.diff(ylay1[:,terms2,p0,bi,se])
                    xdiffl2[:,terms2,p0,bi,se] = ny.diff(xlay2[:,terms2,p0,bi,se])
                    ydiffl2[:,terms2,p0,bi,se] = ny.diff(ylay2[:,terms2,p0,bi,se])



    return time, xlay1, ylay1, xlay2, ylay2, xdiffl1, ydiffl1, xdiffl2, ydiffl2, Namu




##########################################################

#We get the shape index data

def getshapedata(fle,n,cnum):

    Play1 = ny.zeros((n,cnum))
    p0lay1 = ny.zeros((n,cnum))
    Play2 = ny.zeros((n,cnum))
    p0lay2 = ny.zeros((n,cnum))



    with open(fle,'r') as Cfile:
        Clines = Cfile.readlines()
        for lines in range(1,n+1):
            Cdata = Clines[lines].split();
            for terms in range(0,cnum):
                index1 = 2*terms+1
                index2 = 2*terms+2
                index3 = 2*terms+1+2*cnum
                index4 = 2*terms+1+2*cnum+1
                Play1[lines-1,terms] = float(Cdata[index1]) 
                p0lay1[lines-1,terms] = float(Cdata[index2])
                Play2[lines-1,terms] = float(Cdata[index3])
                p0lay2[lines-1,terms] = float(Cdata[index4])
    
    return Play1, p0lay1, Play2, p0lay2

#######################################################

def getallshapedata(n,cnum,np0,nse,nbi,namep0,namesee,namebi,name):

    Play1 = ny.zeros((n,cnum,np0,nbi,nse))
    p0lay1 = ny.zeros((n,cnum,np0,nbi,nse))
    Play2 = ny.zeros((n,cnum,np0,nbi,nse))
    p0lay2 = ny.zeros((n,cnum,np0,nbi,nse))


    Play12 = ny.zeros((np0,nbi,nse))
    Pavg = ny.zeros((np0,nbi))
    Pstd = ny.zeros((np0,nbi))

    p0lay12 = ny.zeros((np0,nbi,nse))

    p0avg = ny.zeros((np0,nbi))
    p0std = ny.zeros((np0,nbi))

    for p0 in range(0,np0):
        for bi in range(0,nbi):
            for se in range(0,nse):

                filenme = namep0[p0] + namesee[se] + "Bi_" + namebi[bi] + name

                Pcells = []
                p0cells = []
                with open(filenme,'r') as Cfile:
                    Clines = Cfile.readlines()
                    for lines in range(1,n+1):
                        Cdata = Clines[lines].split();
                        for terms in range(0,cnum):
                            index1 = 2*terms+1
                            index2 = 2*terms+2
                            index3 = 2*terms+1+2*cnum
                            index4 = 2*terms+1+2*cnum+1
                            Play1[lines-1,terms,p0,bi,se] = float(Cdata[index1]) 
                            p0lay1[lines-1,terms,p0,bi,se] = float(Cdata[index2])
                            Play2[lines-1,terms,p0,bi,se] = float(Cdata[index3])
                            p0lay2[lines-1,terms,p0,bi,se] = float(Cdata[index4])
                      


                Play12[p0,bi,se] = ( ny.mean(Play1[-1,:,p0,bi,se]) + ny.mean(Play2[-1,:,p0,bi,se]) )/2
                p0lay12[p0,bi,se] = ( ny.mean(p0lay1[-1,:,p0,bi,se]) + ny.mean(p0lay2[-1,:,p0,bi,se]) )/2


            Pavg[p0,bi] = ny.mean(Play12[p0,bi,:])
            p0avg[p0,bi] = ny.mean(p0lay12[p0,bi,:])

            Pstd[p0,bi] = ny.std(Play12[p0,bi,:])
            p0std[p0,bi] = ny.std(p0lay12[p0,bi,:])


    return Play1, p0lay1, Play2, p0lay2, Pavg, Pstd, p0avg, p0std




#################################################################

def avgovercellsnospike(dat,dif,xmin, xmax, nli,ncells):


    L = xmax-xmin

    xavg = ny.zeros(nli)

    ckeep = 0
    for cells in range(0,ncells):
        Dif = dif[:,cells]
        if ny.max(Dif) < L/4 and ny.min(Dif) > -L/4:
            xavg+= dat[:,cells]
            ckeep +=1

    if ckeep > 0:
        xavg = xavg/float(ckeep)

    if ckeep == 0:
        print("No reasonable cells found without spikes")
        xavg[:] = ny.nan

    return xavg

#######################################################
def avgovercellsnospikelist(dat,dif,xmin, xmax,n,ncells):


    L = xmax-xmin

    xavg = ny.zeros(n) 

    ckeep = 0
    for cells in range(0,ncells):
        Dif = dif[:,cells]
        if ny.max(Dif) < L/4 and ny.min(Dif) > -L/4:
            xavg+= ny.array(dat[cells])
            ckeep +=1

    if ckeep > 0:
        xavg = list(xavg/float(ckeep))

    if ckeep == 0:
        print("No reasonable cells found without spikes")
        xavg = []

    return xavg



#####################################################

def avgovernonzeroarray(dat,nli,nsed):


    avgary = ny.zeros(nli)


    num = 0

    for seeds in range(0,nsed):
        if ny.any(dat[:,seeds]):
            avgary += dat[:,seeds]

            num+=1

    if num != 0:
        avgary = avgary/float(num)

    if num == 0:
        print("None of the seeds had values")


    return avgary

#########################################################
def avgoverexistingarray(dat,nli,nsed):


    avgary = ny.zeros(nli)


    num = 0

    for seeds in range(0,nsed):
            avgary += dat[:,seeds]
            num+=1

    if num != 0:
        avgary = avgary/float(num)

    if num == 0:
        print("There were no usable seeds")


    return avgary



#########################################################
#Here we do the intercept and the slope for linear regression
#Taken from https://www.geeksforgeeks.org/machine-learning/linear-regression-python-implementation/
#https://www.lancaster.ac.uk/staff/drummonn/PHYS281/demo-line-fitting/

#B1 is the slope

def linreg_coeff(X, Y, n1, n2):

    Xs = X[n1:n2]
    Ys = Y[n1:n2]
    n = ny.size(Xs)
    xmean = ny.mean(Xs)
    ymean = ny.mean(Ys)

    SS_xy = ny.sum(Ys*Xs) - n*ymean*xmean
    SS_xx = ny.sum(Xs*Xs) - n*xmean*xmean

    if SS_xx != 0:

        b1 = SS_xy / SS_xx
    else:
        print("SS_xx was zero", Xs[n1], " ", Xs[n2])
        b1 = 0

    b0 = ymean - b1*xmean

    return b0, b1



#################################################################
#We'll just do a linear polynomial fit

def getpolyfit(X, Y, n1, n2):
    
    b1 = 0
    b0 = 0

    if n2 > n1:

        Xs = X[n1:n2]
        Ys = Y[n1:n2]
        n = ny.size(Xs)

        b1, b0 = ny.polyfit(Xs,Ys,1)

    return b1, b0


###################################################################

def getregline(X, B0, B1, n1, n2):

    xtemp = X[n1:n2]


    Yreg = B0 + B1*xtemp


    return Yreg

###########################################################
#Here we get the final distance for seeds

def getdistavgoverseeds(dat2,nli,nsed):


    distlist = []

    num = 0

 #   check2 = ny.zeros(nli)

    if nsed > 0:
        for seeds in range(0,nsed):
            if len(dat2[seeds] > 0):
                check2 = ny.array(dat2[seeds])
                distlist.append(check2[-1])

                num+=1

    if num == 0:
        print("None of the seeds were good")


    return distlist



#############################################################
#Here we get diffusion coefficients

def getdiffcoeff(time,dat1,n1,n2):

    D0 = 0

    if len(dat1) > 0 and n2 > n1:
        check1 = ny.array(dat1)
        check2 = ny.array(time)
        B1,B0 = getpolyfit(check2, check1, n1, n2)
        D0 = B1

    if len(dat1) == 0:
        print("The seed was not good")


    return D0, B0


##########################################################
#Here we get the exponent of a power law curve

def getpowerexp(time,dat1,n1,n2):

    beta = 0

    if len(dat1) > 0:
        checktime = ny.array(time)
        checkdat = ny.array(dat1)
        n11 = n2
        if checktime[0] == 0:
            checktime = checktime[1:]
            checkdat = checkdat[1:]
            n11 = n2-1
        logdat = ny.log(checkdat)
        logtime = ny.log(checktime)
        #B0,B1 = linreg_coeff(logtime, logdat, n1, n11)
        B1,B0 = getpolyfit(logtime,logdat,n1,n11)
        beta = B1

    if len(dat1) == 0:
        print("The seed was not good")


    return beta


##############################################################

#We get diffusion coefficients for each seed
def getdiffcoefseed(time,dat1,nli,nsed,n1,n2):


    difflist = []

    num = 0

    check1 = ny.zeros(nli)

    for seeds in range(0,nsed):
        if len(dat1[seeds] > 0):
            check1 = ny.array(dat1[seeds])
            B1,B0 = getpolyfit(time, check1, n1, n2)
            D0 = B1
            difflist.append(D0)
            num+=1


    if num == 0:
        print("None of the seeds were good")


    return difflist




##############################################################

#Here we get diffusion coefficients and the final distance for seeds

def getdiffcoefandist(time,dat1,dat2,nli,nsed,n1,n2):


    difflist = []
    distlist = []

    num = 0

    check1 = ny.zeros(nli)
    check2 = ny.zeros(nli)

    for seeds in range(0,nsed):
        if len(dat1[seeds] > 0):
            check1 = ny.array(dat1[seeds])
            B1,B0 = getpolyfit(time, check1, n1, n2)
            D0 = B1
            difflist.append(D0)
            num+=1
        if len(dat2[seeds] > 0):
            check2 = ny.array(dat2[seeds])
            distlist.append(check2[-1])


    if num == 0:
        print("None of the seeds were good")


    return difflist, distlist



###########################################################
##############################################################
#We get the average msd and dist over seeds that don't jump

def getalldistnospikeinfo(comdat,n,ncells,xmin, xmax,ymin,ymax,np0,nse,nbi):

    distkeepx=[[ None for _ in range(nbi) ] for _ in range(np0)]
    distkeepy=[[ None for _ in range(nbi) ] for _ in range(np0)]
    distkeepr=[[ None for _ in range(nbi) ] for _ in range(np0)]
    Distr = ny.zeros((np0,nbi))
    Distrstd = ny.zeros((np0,nbi))

    Lx = xmax-xmin
    Ly = ymax-ymin

    ckeep = 0

    xdistl1seed = ny.zeros((n,np0,nbi,nse))
    ydistl1seed = ny.zeros((n,np0,nbi,nse))
    xdistl2seed = ny.zeros((n,np0,nbi,nse))
    ydistl2seed = ny.zeros((n,np0,nbi,nse))

    rdistseed = ny.zeros((n,np0,nbi,nse))

    xdistl12seed = ny.zeros((n,np0,nbi,nse))
    ydistl12seed = ny.zeros((n,np0,nbi,nse))


    xdist1seedtemp = ny.zeros((n,nse))
    ydist1seedtemp = ny.zeros((n,nse))

    xdistl1 = ny.zeros((n,np0,nbi))
    ydistl1 = ny.zeros((n,np0,nbi))
    xdistl2 = ny.zeros((n,np0,nbi))
    ydistl2 = ny.zeros((n,np0,nbi))
    xdistl12 = ny.zeros((n,np0,nbi))
    ydistl12 = ny.zeros((n,np0,nbi))
    rdist = ny.zeros((n,np0,nbi))

    Time = comdat.timedat

    for p0 in ny.arange(0,np0):
        for bi in ny.arange(0,nbi):


            xmsdtemp = []
            ymsdtemp = []

            xdisttemp = []
            ydisttemp = []
            rdisttemp = []

            for se in ny.arange(0,nse):
               
                xdistl1seed[:,p0,bi,se] = avgovercellsnospike(comdat.xdist1[:,:,p0,bi,se],comdat.delxl1[:,:,p0,bi,se],xmin,xmax,n,ncells)

                ydistl1seed[:,p0,bi,se] = avgovercellsnospike(comdat.ydist1[:,:,p0,bi,se],comdat.delyl1[:,:,p0,bi,se],ymin,ymax,n,ncells)
 
                xdistl2seed[:,p0,bi,se] = avgovercellsnospike(comdat.xdist2[:,:,p0,bi,se],comdat.delxl2[:,:,p0,bi,se],xmin,xmax,n,ncells)

                ydistl2seed[:,p0,bi,se] = avgovercellsnospike(comdat.ydist2[:,:,p0,bi,se],comdat.delyl2[:,:,p0,bi,se],ymin,ymax,n,ncells)


                xdistl12seed[:,p0,bi,se] = xdistl1seed[:,p0,bi,se] + xdistl2seed[:,p0,bi,se]

                ydistl12seed[:,p0,bi,se] = ydistl1seed[:,p0,bi,se] + ydistl2seed[:,p0,bi,se]


#                if not ny.isnan(xdistl1seed[:,p0,bi,se]).any() and if not ny.isnan(xdistl2seed[:,p0,bi,se]).any() and if not ny.isnan(ydistl1seed[:,p0,bi,se]).any() and if not ny.isnan(ydistl2seed[:,p0,bi,se]).any():

                rdistseed[:,p0,bi,se] = ((xdistl1seed[:,p0,bi,se]+xdistl2seed[:,p0,bi,se])**2 + (ydistl1seed[:,p0,bi,se] + ydistl2seed[:,p0,bi,se])**2 )/4
                rdistseed[:,p0,bi,se] = ny.sqrt(rdistseed[:,p0,bi,se])

                rdisttemp.append(rdistseed[:,p0,bi,se])

#                rdisttemp.append(rdistseed[:,p0,bi,se])

#                xdistl12seed[:,p0,bi,se] = (xdistl1seed[:,p0,bi,se]+xdistl2seed[:,p0,bi,se])/2

#                ydistl12seed[:,p0,bi,se] = (ydistl1seed[:,p0,bi,se]+ydistl2seed[:,p0,bi,se])/2
                
                if not ny.isnan(xdistl1seed[:,p0,bi,se]).any():
                    xdisttemp.append(xdistl1seed[:,p0,bi,se])
                if not ny.isnan(xdistl2seed[:,p0,bi,se]).any():
                    xdisttemp.append(xdistl2seed[:,p0,bi,se])
                if not ny.isnan(ydistl1seed[:,p0,bi,se]).any():
                    ydisttemp.append(ydistl1seed[:,p0,bi,se])
                if not ny.isnan(ydistl2seed[:,p0,bi,se]).any():
                    ydisttemp.append(ydistl2seed[:,p0,bi,se])


            distkeepx[p0][bi] = getdistavgoverseeds(xdisttemp,n,len(xdisttemp))

            distkeepy[p0][bi] = getdistavgoverseeds(ydisttemp,n,len(ydisttemp))

            distkeepr[p0][bi] = getdistavgoverseeds(rdisttemp,n,nse)

            distrlist = ny.array(distkeepr[p0][bi]) 
            Distr[p0,bi] = ny.mean(distrlist) 
#            diststd = ny.std(distrlist)
            Distrstd[p0,bi] = ny.std(distrlist)



#            for se in ny.arange(0,nse):
#                xdist1seedtemp[:,se] = 

            xdistl1[:,p0,bi] = avgoverexistingarray(xdistl1seed[:,p0,bi,:],n,nse)
            ydistl1[:,p0,bi] = avgoverexistingarray(ydistl1seed[:,p0,bi,:],n,nse)
            xdistl2[:,p0,bi] = avgoverexistingarray(xdistl2seed[:,p0,bi,:],n,nse)
            ydistl2[:,p0,bi] = avgoverexistingarray(ydistl2seed[:,p0,bi,:],n,nse)





            xdistl12[:,p0,bi] = (xdistl1[:,p0,bi]+xdistl2[:,p0,bi])/2

            ydistl12[:,p0,bi] = (ydistl1[:,p0,bi]+ydistl2[:,p0,bi])/2



            rdist[:,p0,bi] = ny.sqrt((xdistl12[:,p0,bi])**2 + (ydistl12[:,p0,bi] )**2 ) 



    return Time,xdistl12seed, ydistl12seed, rdistseed,xdistl12,ydistl12,rdist,Distr,Distrstd

#

#############################################################
#We get the average msd and dist over seeds that don't jump

def getallmsdanddist(comdat,n,ncells,xmin, xmax,ymin,ymax,np0,nse,nbi,starttime,endtime):


    diffkeepx=[[ None for _ in range(nbi) ] for _ in range(np0)]
    diffkeepy=[[ None for _ in range(nbi) ] for _ in range(np0)]
    diffkeepr=[[ None for _ in range(nbi) ] for _ in range(np0)]
    distkeepx=[[ None for _ in range(nbi) ] for _ in range(np0)]
    distkeepy=[[ None for _ in range(nbi) ] for _ in range(np0)]
    distkeepr=[[ None for _ in range(nbi) ] for _ in range(np0)]

    Diffr = ny.zeros((np0,nbi))
    Diffrstd = ny.zeros((np0,nbi))
    Distr = ny.zeros((np0,nbi))
    Distrstd = ny.zeros((np0,nbi))


    Lx = xmax-xmin
    Ly = ymax-ymin
    xdist = ny.zeros(n)
    ckeep = 0

    xmsdl1seed = ny.zeros((n,np0,nbi,nse))
    ymsdl1seed = ny.zeros((n,np0,nbi,nse))
    xmsdl2seed = ny.zeros((n,np0,nbi,nse))
    ymsdl2seed = ny.zeros((n,np0,nbi,nse))

    xmsdl12seed = ny.zeros((n,np0,nbi,nse))
    ymsdl12seed = ny.zeros((n,np0,nbi,nse))

    rmsdseed = ny.zeros((n,np0,nbi,nse))

    xdistl1seed = ny.zeros((n,np0,nbi,nse))
    ydistl1seed = ny.zeros((n,np0,nbi,nse))
    xdistl2seed = ny.zeros((n,np0,nbi,nse))
    ydistl2seed = ny.zeros((n,np0,nbi,nse))

    xdistl12seed = ny.zeros((n,np0,nbi,nse))
    ydistl12seed = ny.zeros((n,np0,nbi,nse))

    rdistseed = ny.zeros((n,np0,nbi,nse))

    xmsdl1 = ny.zeros((n,np0,nbi))
    ymsdl1 = ny.zeros((n,np0,nbi))
    xmsdl2 = ny.zeros((n,np0,nbi))
    ymsdl2 = ny.zeros((n,np0,nbi))


    xmsdl12 = ny.zeros((n,np0,nbi))
    ymsdl12 = ny.zeros((n,np0,nbi))
    rmsd = ny.zeros((n,np0,nbi))


    xdistl1 = ny.zeros((n,np0,nbi))
    ydistl1 = ny.zeros((n,np0,nbi))
    xdistl2 = ny.zeros((n,np0,nbi))
    ydistl2 = ny.zeros((n,np0,nbi))
    xdistl12 = ny.zeros((n,np0,nbi))
    ydistl12 = ny.zeros((n,np0,nbi))
    rdist = ny.zeros((n,np0,nbi))

    time = ny.zeros((n,np0,nbi))


    for p0 in ny.arange(0,np0):
        for bi in ny.arange(0,nbi):

            xmsdtemp = []
            ymsdtemp = []

            rmsdtemp = []

            xdisttemp = []
            ydisttemp = []

            rdisttemp = []

            time[:,p0,bi] = comdat.timedat[:,p0,bi,0]

            for se in ny.arange(0,nse):

    
                xmsdl1seed[:,p0,bi,se] = avgovercellsnospike(comdat.xmsd1[:,:,p0,bi,se],comdat.delxl1[:,:,p0,bi,se],xmin,xmax,n,ncells)
                   
                xdistl1seed[:,p0,bi,se] = avgovercellsnospike(comdat.xdist1[:,:,p0,bi,se],comdat.delxl1[:,:,p0,bi,se],xmin,xmax,n,ncells)

                ymsdl1seed[:,p0,bi,se] = avgovercellsnospike(comdat.ymsd1[:,:,p0,bi,se],comdat.delyl1[:,:,p0,bi,se],ymin,ymax,n,ncells)

                ydistl1seed[:,p0,bi,se] = avgovercellsnospike(comdat.ydist1[:,:,p0,bi,se],comdat.delyl1[:,:,p0,bi,se],ymin,ymax,n,ncells)


                xmsdl2seed[:,p0,bi,se] = avgovercellsnospike(comdat.xmsd2[:,:,p0,bi,se],comdat.delxl2[:,:,p0,bi,se],xmin,xmax,n,ncells)
 
                xdistl2seed[:,p0,bi,se] = avgovercellsnospike(comdat.xdist2[:,:,p0,bi,se],comdat.delxl2[:,:,p0,bi,se],xmin,xmax,n,ncells)

                ymsdl2seed[:,p0,bi,se] = avgovercellsnospike(comdat.ymsd2[:,:,p0,bi,se],comdat.delyl2[:,:,p0,bi,se],ymin,ymax,n,ncells)

                ydistl2seed[:,p0,bi,se] = avgovercellsnospike(comdat.ydist2[:,:,p0,bi,se],comdat.delyl2[:,:,p0,bi,se],ymin,ymax,n,ncells)


                rmsdseed[:,p0,bi,se] = (xmsdl1seed[:,p0,bi,se]+xmsdl2seed[:,p0,bi,se])/2 + (ymsdl1seed[:,p0,bi,se]+ymsdl2seed[:,p0,bi,se])/2 


                rdistseed[:,p0,bi,se] = ((xdistl1seed[:,p0,bi,se]+xdistl2seed[:,p0,bi,se])**2 + (ydistl1seed[:,p0,bi,se] + ydistl2seed[:,p0,bi,se])**2)/4
                rdistseed[:,p0,bi,se] = ny.sqrt(rdistseed[:,p0,bi,se])

                rmsdtemp.append(rmsdseed[:,p0,bi,se])

                rdisttemp.append(rdistseed[:,p0,bi,se])

#                xmsdl12seed[:,p0,bi,se] = (xmsdl1seed[:,p0,bi,se]+xmsdl2seed[:,p0,bi,se])/2
#                ymsdl12seed[:,p0,bi,se] = (ymsdl1seed[:,p0,bi,se]+ymsdl2seed[:,p0,bi,se])/2
#                xdistl12seed[:,p0,bi,se] = (xdistl1seed[:,p0,bi,se]+xdistl2seed[:,p0,bi,se])/2
#                ydistl12seed[:,p0,bi,se] = (ydistl1seed[:,p0,bi,se]+ydistl2seed[:,p0,bi,se])/2
  
                if not ny.isnan(xmsdl1seed[:,p0,bi,se]).any():
                    xmsdtemp.append(xmsdl1seed[:,p0,bi,se])
                if not ny.isnan(xmsdl2seed[:,p0,bi,se]).any():
                    xmsdtemp.append(xmsdl2seed[:,p0,bi,se])
                if not ny.isnan(ymsdl1seed[:,p0,bi,se]).any():
                    ymsdtemp.append(ymsdl1seed[:,p0,bi,se])
                if not ny.isnan(ymsdl2seed[:,p0,bi,se]).any():
                    ymsdtemp.append(ymsdl2seed[:,p0,bi,se])
                if not ny.isnan(xdistl1seed[:,p0,bi,se]).any():
                    xdisttemp.append(xdistl1seed[:,p0,bi,se])
                if not ny.isnan(xdistl2seed[:,p0,bi,se]).any():
                    xdisttemp.append(xdistl2seed[:,p0,bi,se])
                if not ny.isnan(ydistl1seed[:,p0,bi,se]).any():
                    ydisttemp.append(ydistl1seed[:,p0,bi,se])
                if not ny.isnan(xdistl1seed[:,p0,bi,se]).any():
                    ydisttemp.append(ydistl2seed[:,p0,bi,se])


#            diffkeepx[p0][bi], distkeepx[p0][bi] = getdiffcoefandist(comdat.timedat[:,p0,bi,se],xmsdtemp,xdisttemp,n,len(xmsdtemp),starttime,endtime)

#            diffkeepy[p0][bi], distkeepy[p0][bi] = getdiffcoefandist(comdat.timedat[:,p0,bi,se],ymsdtemp,ydisttemp,n,len(xmsdtemp),starttime,endtime)


            diffkeepr[p0][bi] = getdiffcoefseed(comdat.timedat[:,p0,bi,se],rmsdtemp,n,nse,starttime,endtime)

            distkeepr[p0][bi] = getdistavgoverseeds(rdisttemp,n,nse)


            Diffr[p0,bi] = ny.mean(ny.array(diffkeepr[p0][bi])/4)
            Diffrstd[p0,bi] = ny.std(ny.array(diffkeepr[p0][bi])/4)

            Distr[p0,bi] = ny.mean(ny.array(distkeepr[p0][bi]))
            Distrstd[p0,bi] = ny.std(ny.array(distkeepr[p0][bi]))


            xmsdl1[:,p0,bi] = avgoverexistingarray(xmsdl1seed[:,p0,bi,:],n,nse)
            ymsdl1[:,p0,bi] = avgoverexistingarray(ymsdl1seed[:,p0,bi,:],n,nse)
            xmsdl2[:,p0,bi] = avgoverexistingarray(xmsdl2seed[:,p0,bi,:],n,nse)
            ymsdl2[:,p0,bi] = avgoverexistingarray(ymsdl2seed[:,p0,bi,:],n,nse)

            xdistl1[:,p0,bi] = avgoverexistingarray(xdistl1seed[:,p0,bi,:],n,nse)
            ydistl1[:,p0,bi] = avgoverexistingarray(ydistl1seed[:,p0,bi,:],n,nse)
            xdistl2[:,p0,bi] = avgoverexistingarray(xdistl2seed[:,p0,bi,:],n,nse)
            ydistl2[:,p0,bi] = avgoverexistingarray(ydistl2seed[:,p0,bi,:],n,nse)


            xmsdl12[:,p0,bi] = (xmsdl1[:,p0,bi]+xmsdl2[:,p0,bi])/2

            ymsdl12[:,p0,bi] = (ymsdl1[:,p0,bi]+ymsdl2[:,p0,bi])/2

            xdistl12[:,p0,bi] = (xdistl1[:,p0,bi]+xdistl2[:,p0,bi])/2

            ydistl12[:,p0,bi] = (ydistl1[:,p0,bi]+ydistl2[:,p0,bi])/2



            rmsd[:,p0,bi] = xmsdl12[:,p0,bi] + ymsdl12[:,p0,bi]
            rdist[:,p0,bi] = ny.sqrt((xdistl12[:,p0,bi])**2 + (ydistl12[:,p0,bi] )**2 ) 



    return time, xmsdl12seed, ymsdl12seed, xdistl12seed, ydistl12seed, xmsdl12, ymsdl12,rmsd,xdistl12,ydistl12,rdist,Diffr,Diffrstd,Distr,Distrstd

##################################################################
#We get the distance traveled at rise time

def getalldistatrise(Dis,rtim,np0,nbi,nse):

    rdistatriseseed = ny.zeros((np0,nbi,nse))

    rdistatrise = ny.zeros((np0,nbi))
    rdistatrisestd = ny.zeros((np0,nbi))

    speatriseseed = ny.zeros((np0,nbi,nse))
    speatrisetemp = []
    speatrise = ny.zeros((np0,nbi))
    speatrisestd = ny.zeros((np0,nbi))
    speatrisealt = ny.zeros((np0,nbi))


    for p0 in range(0,np0):
        for bi in range(0,nbi):
            speatrisetemp = []
            for se in range(0,nse):
                Rtime = math.ceil(rtim.edgerise[p0,bi,se]/10)
                rdistatriseseed[p0,bi,se] = Dis.rdistseed[Rtime,p0,bi,se]
                if Rtime >0:
                    speatriseseed[p0,bi,se] = Dis.rdistseed[Rtime,p0,bi,se]/float(rtim.edgerise[p0,bi,se])
                    sped = Dis.rdistseed[Rtime,p0,bi,se]/float(rtim.edgerise[p0,bi,se])
                    speatrisetemp.append(sped)
                else:
                    speatriseseed[p0,bi,se] = 0


            edgeriseseedavg = ny.mean(rtim.edgerise[p0,bi,:])
            rdistatrise[p0,bi] = ny.mean(rdistatriseseed[p0,bi,:])
            rdistatrisestd[p0,bi] = ny.mean(rdistatriseseed[p0,bi,:])

            if len(speatrisetemp) > 0:
                speatrise[p0,bi] = ny.mean(ny.array(speatrisetemp))
                speatrisestd[p0,bi] = ny.std(speatriseseed[p0,bi,:])
            else:
                speatrise[p0,bi] = 10000
                speatrisestd[p0,bi] = 0
#            speatrise[p0,bi] = ny.mean(speatriseseed[p0,bi,:])
#            speatrisestd[p0,bi] = ny.std(speatriseseed[p0,bi,:])
            speatrisealt[p0,bi] = rdistatrise[p0,bi]/float(edgeriseseedavg)




    return rdistatriseseed,speatriseseed,rdistatrise,speatrise,rdistatrisestd,speatrisestd, speatrisealt


##############################################################
#We get the MSD at rise time

def getallmsdatrise(comdat,ncells,xmin,xmax,ymin,ymax,np0,nse,nbi,rismd):


    diffkeepx=[[ None for _ in range(nbi) ] for _ in range(np0)]
    diffkeepy=[[ None for _ in range(nbi) ] for _ in range(np0)]
    diffkeepr=[[ None for _ in range(nbi) ] for _ in range(np0)]
#    distkeepx=[[ None for _ in range(nbi) ] for _ in range(np0)]
#    distkeepy=[[ None for _ in range(nbi) ] for _ in range(np0)]
#    distkeepr=[[ None for _ in range(nbi) ] for _ in range(np0)]

    Diffr = ny.zeros((np0,nbi))
    Diffrstd = ny.zeros((np0,nbi))
    Coef = ny.zeros((np0,nbi))
    Coefstd = ny.zeros((np0,nbi))
#    Distr = ny.zeros((np0,nbi))
#    Distrstd = ny.zeros((np0,nbi))


    Lx = xmax-xmin
    Ly = ymax-ymin
    ckeep = 0

    xmsdl1seed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]

    ymsdl1seed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]

    xmsdl2seed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]

    ymsdl2seed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]


    xmsdl12seed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]

    ymsdl12seed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]


    rmsdseed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]


#    xdistl1seed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]

#    ydistl1seed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]

#    xdistl2seed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]

#    ydistl2seed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]


#    xdistl12seed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]

#    ydistl12seed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]


#    rdistseed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]


#    time = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]


    for p0 in ny.arange(0,np0):
        for bi in ny.arange(0,nbi):

            xmsdtemp = []
            ymsdtemp = []

            rmsdtemp = []

            xdisttemp = []
            ydisttemp = []

            rdisttemp = []

#            time[:,p0,bi] = comdat.timedat[:,p0,bi,0]

            dcotemp = []

            coeftemp = []

            for se in ny.arange(0,nse):

                if rismd.xmsd1[p0][bi][se][0] is not None:

                    n = len(rismd.xmsd1[p0][bi][se][0])
    
                    xmsdl1seed[p0][bi][se] = avgovercellsnospikelist(rismd.xmsd1[p0][bi][se],comdat.delxl1[:,:,p0,bi,se],xmin,xmax,n,ncells)
                  
                    ymsdl1seed[p0][bi][se] = avgovercellsnospikelist(rismd.ymsd1[p0][bi][se],comdat.delyl1[:,:,p0,bi,se],ymin,ymax,n,ncells)

                    xmsdl2seed[p0][bi][se] = avgovercellsnospikelist(rismd.xmsd2[p0][bi][se],comdat.delxl2[:,:,p0,bi,se],xmin,xmax,n,ncells)

                    ymsdl2seed[p0][bi][se] = avgovercellsnospikelist(rismd.ymsd2[p0][bi][se],comdat.delyl2[:,:,p0,bi,se],ymin,ymax,n,ncells)

                    rtemp = (ny.array(xmsdl1seed[p0][bi][se])+ny.array(xmsdl2seed[p0][bi][se]))/2 + (ny.array(ymsdl1seed[p0][bi][se]) + ny.array(ymsdl2seed[p0][bi][se]))/2

                    rmsdseed[p0][bi][se] = list(rtemp)


#                    xdistl1seed[:,p0,bi,se] = avgovercellsnospike(comdat.xdist1[:,:,p0,bi,se],comdat.delxl1[:,:,p0,bi,se],xmin,xmax,n,ncells)

#                    ydistl1seed[:,p0,bi,se] = avgovercellsnospike(comdat.ydist1[:,:,p0,bi,se],comdat.delyl1[:,:,p0,bi,se],ymin,ymax,n,ncells)


#                    xdistl2seed[:,p0,bi,se] = avgovercellsnospike(comdat.xdist2[:,:,p0,bi,se],comdat.delxl2[:,:,p0,bi,se],xmin,xmax,n,ncells)

#                    ydistl2seed[:,p0,bi,se] = avgovercellsnospike(comdat.ydist2[:,:,p0,bi,se],comdat.delyl2[:,:,p0,bi,se],ymin,ymax,n,ncells)


#                    rdistseed[:,p0,bi,se] = ((xdistl1seed[:,p0,bi,se]+xdistl2seed[:,p0,bi,se])**2 + (ydistl1seed[:,p0,bi,se] + ydistl2seed[:,p0,bi,se])**2)/4
#                    rdistseed[:,p0,bi,se] = ny.sqrt(rdistseed[:,p0,bi,se])

#                    rmsdtemp.append(rmsdseed[:,p0,bi,se])

#                    rdisttemp.append(rdistseed[:,p0,bi,se])

#                xmsdl12seed[:,p0,bi,se] = (xmsdl1seed[:,p0,bi,se]+xmsdl2seed[:,p0,bi,se])/2
#                ymsdl12seed[:,p0,bi,se] = (ymsdl1seed[:,p0,bi,se]+ymsdl2seed[:,p0,bi,se])/2
#                xdistl12seed[:,p0,bi,se] = (xdistl1seed[:,p0,bi,se]+xdistl2seed[:,p0,bi,se])/2
#                ydistl12seed[:,p0,bi,se] = (ydistl1seed[:,p0,bi,se]+ydistl2seed[:,p0,bi,se])/2
  
#                if not ny.isnan(xmsdl1seed[:,p0,bi,se]).any():
#                    xmsdtemp.append(xmsdl1seed[:,p0,bi,se])
#                if not ny.isnan(xmsdl2seed[:,p0,bi,se]).any():
#                    xmsdtemp.append(xmsdl2seed[:,p0,bi,se])
#                if not ny.isnan(ymsdl1seed[:,p0,bi,se]).any():
#                    ymsdtemp.append(ymsdl1seed[:,p0,bi,se])
#                if not ny.isnan(ymsdl2seed[:,p0,bi,se]).any():
#                    ymsdtemp.append(ymsdl2seed[:,p0,bi,se])
#                if not ny.isnan(xdistl1seed[:,p0,bi,se]).any():
#                    xdisttemp.append(xdistl1seed[:,p0,bi,se])
#                if not ny.isnan(xdistl2seed[:,p0,bi,se]).any():
#                    xdisttemp.append(xdistl2seed[:,p0,bi,se])
#                if not ny.isnan(ydistl1seed[:,p0,bi,se]).any():
#                    ydisttemp.append(ydistl1seed[:,p0,bi,se])
#                if not ny.isnan(xdistl1seed[:,p0,bi,se]).any():
#                    ydisttemp.append(ydistl2seed[:,p0,bi,se])


#            diffkeepx[p0][bi], distkeepx[p0][bi] = getdiffcoefandist(comdat.timedat[:,p0,bi,se],xmsdtemp,xdisttemp,n,len(xmsdtemp),starttime,endtime)

#            diffkeepy[p0][bi], distkeepy[p0][bi] = getdiffcoefandist(comdat.timedat[:,p0,bi,se],ymsdtemp,ydisttemp,n,len(xmsdtemp),starttime,endtime)


                    tuse = ny.array(rismd.timedat[p0][bi][se])

                    if n > 3 :

                        Dco, intercpt = getdiffcoeff(tuse,rtemp,0,n//2)

                        Bet = getpowerexp(tuse,rtemp,0,n//2)

                        dcotemp.append(Dco)

                        coeftemp.append(Bet)

            diffkeepr[p0][bi] = dcotemp

#            distkeepr[p0][bi] = getdistavgoverseeds(rdisttemp,n,nse)


            Diffr[p0,bi] = ny.mean(ny.array(diffkeepr[p0][bi])/4)
            Diffrstd[p0,bi] = ny.std(ny.array(diffkeepr[p0][bi])/4)

            Coef[p0,bi] = ny.mean(ny.array(coeftemp))

#            Distr[p0,bi] = ny.mean(ny.array(distkeepr[p0][bi]))
#            Distrstd[p0,bi] = ny.std(ny.array(distkeepr[p0][bi]))


#            xmsdl1[:,p0,bi] = avgoverexistingarray(xmsdl1seed[:,p0,bi,:],n,nse)
#            ymsdl1[:,p0,bi] = avgoverexistingarray(ymsdl1seed[:,p0,bi,:],n,nse)
#            xmsdl2[:,p0,bi] = avgoverexistingarray(xmsdl2seed[:,p0,bi,:],n,nse)
#            ymsdl2[:,p0,bi] = avgoverexistingarray(ymsdl2seed[:,p0,bi,:],n,nse)

#            xdistl1[:,p0,bi] = avgoverexistingarray(xdistl1seed[:,p0,bi,:],n,nse)
#            ydistl1[:,p0,bi] = avgoverexistingarray(ydistl1seed[:,p0,bi,:],n,nse)
#            xdistl2[:,p0,bi] = avgoverexistingarray(xdistl2seed[:,p0,bi,:],n,nse)
#            ydistl2[:,p0,bi] = avgoverexistingarray(ydistl2seed[:,p0,bi,:],n,nse)


#            xmsdl12[:,p0,bi] = (xmsdl1[:,p0,bi]+xmsdl2[:,p0,bi])/2

#            ymsdl12[:,p0,bi] = (ymsdl1[:,p0,bi]+ymsdl2[:,p0,bi])/2

#            xdistl12[:,p0,bi] = (xdistl1[:,p0,bi]+xdistl2[:,p0,bi])/2

#            ydistl12[:,p0,bi] = (ydistl1[:,p0,bi]+ydistl2[:,p0,bi])/2



#            rmsd[:,p0,bi] = xmsdl12[:,p0,bi] + ymsdl12[:,p0,bi]
#            rdist[:,p0,bi] = ny.sqrt((xdistl12[:,p0,bi])**2 + (ydistl12[:,p0,bi] )**2 ) 



    return xmsdl1seed, ymsdl1seed, xmsdl2seed, ymsdl2seed, rmsdseed, Diffr,Diffrstd, Coef



#################################################################

#We plot COM curves

def plotCOM(COM,ncell,P,B,S,floo,titl,xy):
    
    cmap = plt.get_cmap('rainbow',ncell)

    plt.figure(floo)
    plt.grid()
    plt.title(titl)

    if xy == 0:

        for cell in range(0,ncell):
            plt.plot(COM.timedat[:,P,B,S],COM.xrawl1[:,cell,P,B,S], color = cmap(cell), label = cell )
            plt.plot(COM.timedat[:,P,B,S],COM.xrawl2[:,cell,P,B,S], color = cmap(cell+ncell), label = cell )

    else:
        for cell in range(0,ncell):
            plt.plot(COM.timedat[:,P,B,S],COM.yrawl1[:,cell,P,B,S], color = cmap(cell), label = cell )
            plt.plot(COM.timedat[:,P,B,S],COM.yrawl2[:,cell,P,B,S], color = cmap(cell+ncell), label = cell )


    plt.tick_params(axis='both', labelsize=28)
    plt.legend(loc='upper left')

    floo+= 1


    return floo


################################################################
#We plot distance curves for the x direction

def plotdist(COM,ncell,P,B,S,floo,titl,xy):
 
    cmap = plt.get_cmap('rainbow',ncell)

    plt.figure(floo)
    plt.grid()
    plt.title(titl)

    if xy == 0:

        for cell in range(0,ncell):
            plt.plot(COM.timedat[:,P,B,S],COM.xdist1[:,cell,P,B,S], color = cmap(cell), label = cell )
            plt.plot(COM.timedat[:,P,B,S],COM.xdist2[:,cell,P,B,S], color = cmap(cell), label = cell )

    else:
        for cell in range(0,ncell):
            plt.plot(COM.timedat[:,P,B,S],COM.ydist1[:,cell,P,B,S], color = cmap(cell), label = cell )
            plt.plot(COM.timedat[:,P,B,S],COM.ydist2[:,cell,P,B,S], color = cmap(cell), label = cell )


    plt.tick_params(axis='both', labelsize=28)
    plt.legend(loc='lower left')

    floo += 1
    return floo

#################################################################
#We plot seed averaged walking distance

def plotdistseed(COM,P,B,Sran,floo,titl,slabel,xy):
 
    cmap = plt.get_cmap('rainbow',Sran)

    plt.figure(floo)
    plt.grid()
    plt.title(titl)


    if xy == 0:

        for seed in range(0,Sran):
            plt.plot(COM.time[:,P,B,seed],COM.xdistseed[:,P,B,seed], color = cmap(seed), label = slabel[seed] )


    else:
        for seed in range(0,Sran):
            plt.plot(COM.time[:,P,B,seed],COM.ydistseed[:,P,B,seed], color = cmap(seed), label = slabel[seed] )



    plt.tick_params(axis='both', labelsize=28)
    plt.legend(loc='lower left')

    floo += 1
    return floo




################################################################
#We plot ensemble averaged walking distance

def plotenavgdist(COM,Pran,B,floo,titl,plabel,xy):
 
    cmap = plt.get_cmap('rainbow',Pran)

    plt.figure(floo)
    plt.grid()
    plt.title(titl)


    if xy == 0:

        for p0 in range(0,Pran):
            plt.plot(COM.time[:,p0,B],COM.xdist[:,p0,B], color = cmap(p0), label = plabel[p0] )


    else:
        for p0 in range(0,Pran):
            plt.plot(COM.time[:,p0,B],COM.ydist[:,p0,B], color = cmap(p0), label = plabel[p0] )



    plt.tick_params(axis='both', labelsize=28)
    plt.legend(loc='lower left')

    floo += 1
    return floo

###############################################################
#We plot msd

def plotenavgmsd(COM,Pran,B,floo,titl,plabel,xy):
 
    cmap = plt.get_cmap('rainbow',Pran)

    plt.figure(floo)
    plt.grid()
    plt.title(titl)


    if xy == 0:

        for p0 in range(0,Pran):
            plt.plot(COM.time[:,p0,B],COM.xmsd[:,p0,B], color = cmap(p0), label = plabel[p0] )


    else:
        for p0 in range(0,Pran):
            plt.plot(COM.time[:,p0,B],COM.ymsd[:,p0,B], color = cmap(p0), label = plabel[p0] )


    plt.tick_params(axis='both', labelsize=28)
    plt.legend(loc='lower left')

    floo += 1
    return floo

############################################################



###########################################################

def plotrdistseeds(COM,P,B,Sran,floo,titl,slabel):

    cmap = plt.get_cmap('rainbow',Sran)

    plt.figure(floo)
    plt.grid()
    plt.title(titl)

    for seed in range(0,Sran):
        plt.plot(COM.time[:,P,B,seed],COM.rdistseed[:,P,B,seed], color = cmap(seed), label = slabel[seed] )
#        plt.axvline(x=etime*10, color='k', linestyle='--')


    plt.tick_params(axis='both', labelsize=28)
#    plt.legend(loc='lower left')

    floo += 1
    return floo


############################################################
def plotrdistseedsatrise(COM,Ris,P,B,Sran,floo,titl,slabel):

    cmap = plt.get_cmap('rainbow',Sran)

    plt.figure(floo)
    plt.grid()
    plt.title(titl)

    for seed in range(0,Sran):
        Rtime = math.ceil(Ris.edgerise[P,B,seed]/10)
        plt.plot(COM.time[:Rtime,P,B,seed],COM.rdistseed[:Rtime,P,B,seed], color = cmap(seed), label = slabel[seed] )
#        plt.axvline(x=etime*10, color='k', linestyle='--')


    plt.tick_params(axis='both', labelsize=28)
#    plt.legend(loc='lower left')

    floo += 1
    return floo





#############################################################

def plotdistr(COM,Pran,B,floo,titl,plabel):

    cmap = plt.get_cmap('rainbow',Pran)

    plt.figure(floo)
    plt.grid()
    plt.title(titl)

    for p0 in range(0,Pran):
        plt.plot(COM.time[:,p0,B],COM.rdist[:,p0,B], color = cmap(p0), label = plabel[p0] )
#        plt.axvline(x=etime*10, color='k', linestyle='--')


    plt.tick_params(axis='both', labelsize=28)
    plt.legend(loc='lower left')

    floo += 1
    return floo


############################################################
#We plot msdR

def plotmsdr(COM,Pran,B,floo,titl,plabel,etime):
 
    cmap = plt.get_cmap('rainbow',Pran)

    plt.figure(floo)
    plt.grid()
    plt.title(titl)


    for p0 in range(0,Pran):
        plt.plot(COM.time[:,p0,B],COM.rmsd[:,p0,B], color = cmap(p0), label = plabel[p0] )
#        plt.axvline(x=etime*10, color='k', linestyle='--')


    plt.tick_params(axis='both', labelsize=28)
    plt.legend(loc='upper left', fontsize = 20)

    floo += 1
    return floo


#############################################################
#We plot msdR at rise

def plotmsdratrise(Tra,P,B,Sran,floo,titl,plabel,etime):
 
    cmap = plt.get_cmap('rainbow',Sran)

    plt.figure(floo)
    plt.grid()
    plt.title(titl)


    for seed in range(0,Sran):
        plt.plot(Tra.timedat[P][B][seed],Tra.rmsdseed[P][B][seed], color = cmap(seed), label = plabel[seed] )
#        plt.axvline(x=etime*10, color='k', linestyle='--')


    plt.tick_params(axis='both', labelsize=28)
#    plt.legend(loc='upper left', fontsize = 20)

    floo += 1
    return floo



################################################################

#We plot the shape index and diffusion coefficient

def plotshapeanddiffcoeff(Shape,Diffc,pvalues,Bfix, floop):


    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    color1 = 'tab:blue'
    ax1.scatter(pvalues,Shape.p0avg[:,Bfix], color = color1, s=50)
    ax1.plot(pvalues,Shape.p0avg[:,Bfix], color = color1)
    ax1.errorbar(pvalues,Shape.p0avg[:,Bfix], yerr=Shape.p0std[:,Bfix], capsize=3, fmt = 'o', color = color1)
    ax1.tick_params(axis='both', labelsize = 28)
    ax1.grid(True, which='both', axis='both', linestyle='-', linewidth=0.5)


    color2 = 'tab:red'
    ax2.scatter(pvalues,Diffc.diffcoeffavg[:,Bfix], color=color2, s=50)
    ax2.plot(pvalues,Diffc.diffcoeffavg[:,Bfix], color=color2)
    ax2.errorbar(pvalues,Diffc.diffcoeffavg[:,Bfix], yerr=Diffc.diffcoeffstd[:,Bfix], capsize=3, fmt = 's', color=color2)
    ax2.tick_params(axis='both', labelsize = 28)
    ax2.grid(True, which='both', axis='both', linestyle='-', linewidth=0.5)

    fig.tight_layout()


    floop +=1
    return floop


##################################################################
#We plot neighbor changes and dist traveled

def plotneighanddistfin(Nei,Dist,pvalues,Bfix, floop):

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    color1 = 'tab:blue'
    ax1.scatter(pvalues,Nei.neichgavg[:,Bfix], color = color1, s=50)
    ax1.plot(pvalues,Nei.neichgavg[:,Bfix], color = color1)
    ax1.errorbar(pvalues,Nei.neichgavg[:,Bfix], yerr=Nei.neichgstd[:,Bfix], capsize=3, fmt = 'o', color = color1)
    ax1.tick_params(axis='both', labelsize = 28)
    ax1.grid(True, which='both', axis='both', linestyle='-', linewidth=0.5)


    color2 = 'tab:red'
    ax2.scatter(pvalues,Dist.findistavg[:,Bfix], color=color2, s=50)
    ax2.plot(pvalues,Dist.findistavg[:,Bfix], color=color2)
    ax2.errorbar(pvalues,Dist.findistavg[:,Bfix], yerr=Dist.findiststd[:,Bfix], capsize=3, fmt = 's', color=color2)
    ax2.tick_params(axis='both', labelsize = 28)
    ax2.grid(True, which='both', axis='both', linestyle='-', linewidth=0.5)

    fig.tight_layout()


    floop +=1
    return floop




###################################################################

#We need to smooth the data to account for periodic boundary conditions

def Smoothxcomdata1(Xcom,nlines,ncells,xmin,xmax):

    Lx = xmax - xmin

    Xnew = ny.zeros((nlines,ncells))
    crossing = ny.zeros((nlines,ncells))

    for lines in range (1,nlines):
        for cells in range (0,ncells):
            delX = Xcom[lines,cells]-Xcom[0,cells]
#            delX = Xcom[lines,cells]-Xcom[0,cells]
            if delX > 3*Lx/4:
                crossing[lines,cells] = 1
            elif delX < -3*Lx/4:
                crossing[lines,cells] = -1


    for cells in range(0,ncells):
        Xnew[0,cells] = Xcom[0,cells]

    for lines in range (1,nlines):
        for cells in range (0,ncells):
            if crossing[lines,cells]==1:
                Xnew[lines,cells] = Xcom[lines,cells] -Lx
            elif crossing[lines,cells]==-1:
                Xnew[lines,cells] = Xcom[lines,cells] + Lx
            else:
                Xnew[lines,cells] = Xcom[lines,cells]


    return Xnew

   
##############################################
def Smoothxcomdata2(Xcom,nlines,ncells,xmin,xmax):

    Lx = xmax - xmin
    perioddist = ny.zeros((nlines,ncells))  

    Xnew = ny.zeros((nlines,ncells))

    crossing = ny.zeros((nlines,ncells))

    for lines in range (1,nlines):
        for cells in range (0,ncells):
            delX = Xcom[lines,cells]-Xcom[lines-1,cells]
            if delX > Lx:
                crossing[lines,cells] = 1
                perioddist[lines,cells] = Xcom[lines,cells]-xmax
            elif delX < Lx:
                crossing[lines,cells] = -1
                perioddist[lines,cells] = Xcom[lines,cells]-xmin


    for cells in range(0,ncells):
        Xnew[0,cells] = Xcom[0,cells]

    for lines in range (1,nlines):
        for cells in range (0,ncells):
            if crossing[lines,cells]==1:
                Xnew[lines,cells] = xmin-abs(perioddist[lines,cells])
            elif crossing[lines,cells]==-1:
                Xnew[lines,cells] = xmax+abs(perioddist[lines,cells])
            else:
                Xnew[lines,cells] = Xcom[lines,cells]


    return Xnew

   
###########################################################

def Smoothxcomdata3(Xcom,nlines,ncells,xmin,xmax):

    Lx = xmax - xmin

    Xnew = ny.array(Xcom)

    for lines in range (1,nlines):
        for cells in range (0,ncells):
            delX = Xnew[lines,cells]-Xnew[lines-1,cells]
            if delX > Lx/2:
                crossdist = Xnew[lines,cells]-xmax
                Xnew[lines,cells]= xmin - crossdist
            elif delX < -Lx/2:
                crossdist = xmin - Xnew[lines,cells]
                Xnew[lines,cells] = xmax + crossdist



    return Xnew

###################################################


def Smoothxcomdata3doub(Xlay1, Xlay2, nlines,ncells,xmin,xmax):

    Lx = xmax - xmin

    Xnew1 = ny.array(Xlay1)
    Xnew2 = ny.array(Xlay2)


    for lines in range (1,nlines):
        for cells in range (0,ncells):
            delX1 = Xnew1[lines,cells]-Xnew1[lines-1,cells]
            delX2 = Xnew2[lines,cells]-Xnew2[lines-1,cells]
            if delX1 > Lx/2:
                crossdist = Xnew1[lines,cells]-xmax
                Xnew1[lines,cells]= xmin - crossdist
            if delX2 > Lx/2:
                crossdist = Xnew1[lines,cells]-xmax
                Xnew2[lines,cells]= xmin - crossdist
            if delX1 < -Lx/2:
                crossdist = xmin - Xnew1[lines,cells]
                Xnew1[lines,cells] = xmax + crossdist
            if delX2 < -Lx/2:
                crossdist = xmin - Xnew1[lines,cells]
                Xnew2[lines,cells] = xmax + crossdist


    return Xnew1, Xnew2

#########################################################

def Smoothxcomdata3all(Xlay1, Xlay2, Ylay1, Ylay2, nlines,ncells,xmin,xmax,ymin,ymax):

    Lx = xmax - xmin
    Ly = ymax - ymin

    Xnew1 = ny.array(Xlay1)
    Xnew2 = ny.array(Xlay2)
    Ynew1 = ny.array(Ylay1)
    Ynew2 = ny.array(Ylay2)

    for lines in range (1,nlines):
        for cells in range (0,ncells):
            delX1 = Xnew1[lines,cells]-Xnew1[lines-1,cells]
            delX2 = Xnew2[lines,cells]-Xnew2[lines-1,cells]
            delY1 = Ynew1[lines,cells]-Ynew1[lines-1,cells]
            delY2 = Ynew2[lines,cells]-Ynew2[lines-1,cells]

            if delX1 > Lx/2:
                crossdist = Xnew1[lines,cells]-xmax
                Xnew1[lines,cells]= xmin - crossdist
            if delX2 > Lx/2:
                crossdist = Xnew1[lines,cells]-xmax
                Xnew2[lines,cells]= xmin - crossdist
            if delY1 > Ly/2:
                crossdist = Ynew1[lines,cells]-ymax
                Ynew1[lines,cells]= ymin - crossdist
            if delY2 > Ly/2:
                crossdist = Ynew1[lines,cells]-ymax
                Ynew2[lines,cells]= ymin - crossdist
            if delX1 < -Lx/2:
                crossdist = xmin - Xnew1[lines,cells]
                Xnew1[lines,cells] = xmax + crossdist
            if delX2 < -Lx/2:
                crossdist = xmin - Xnew1[lines,cells]
                Xnew2[lines,cells] = xmax + crossdist
            if delY1 < -Ly/2:
                crossdist = ymin - Ynew1[lines,cells]
                Ynew1[lines,cells] = ymax + crossdist
            if delY2 < -Ly/2:
                crossdist = ymin - Ynew1[lines,cells]
                Ynew2[lines,cells] = ymax + crossdist



    return Xnew1, Xnew2, Ynew1, Ynew2




###################################################

def Smoothxcomdatazero(Xcom,nlines,ncells,xmin,xmax):

    Lx = xmax - xmin

    Xnew = ny.array(Xcom)

    for lines in range (1,nlines):
        for cells in range (0,ncells):
            delX = Xnew[lines,cells]-Xnew[lines-1,cells]
            if delX > Lx/2:
                Xnew[lines,cells] = Xnew[lines-1,cells]
            elif delX < -Lx/2:
                Xnew[lines,cells] = Xnew[lines-1,cells]



    return Xnew





################################################
#This will split the grid into 10 local regions and track the evolution of the edge matching in each region 

def getbimatchdata(fle,n,cnum):
    time = ny.zeros(n)
    xlay1 = ny.zeros((n,cnum))
    ylay1 = ny.zeros((n,cnum))
    rlay1 = ny.zeros((n,cnum))
    xlay2 = ny.zeros((n,cnum))
    ylay2 = ny.zeros((n,cnum))
    rlay2 = ny.zeros((n,cnum))
    with open(fle,'r') as Cfile:
        Clines = Cfile.readlines()
        for lines in range(1,n):
            Cdata = Clines[lines].split();
            time[lines] = float(Cdata[0])
            lim1 = cnum+1
            lindx1 = 1
            for clay1 in range(1,lim1):
                cindx1 = clay1-1
                xlay1[lines,cindx1] = float(Cdata[lindx1])
                lindx1 = lindx1+1
                ylay1[lines,cindx1] = float(Cdata[lindx1])
                rlay1[lines,cindx1] = xlay1[lines,cindx1]*xlay1[lines,cindx1]+ylay1[lines,cindx1]*ylay1[lines,cindx1]
                rlay1[lines,cindx1] = math.sqrt(rlay1[lines,cindx1])
                lindx1 = lindx1+1

            lim2 = 2*(cnum)+1
            lindx2 = lindx1
            for clay2 in range(lim1,lim2):
                cindx2 = clay2-lim1
                xlay2[lines,cindx2] = float(Cdata[lindx2])
                lindx2 = lindx2 + 1
                ylay2[lines,cindx2] = float(Cdata[lindx2])
                rlay2[lines,cindx2] = xlay2[lines,cindx2]*xlay2[lines,cindx2]+ylay2[lines,cindx2]*ylay2[lines,cindx2]
                rlay2[lines,cindx2] = math.sqrt(rlay2[lines,cindx2])
                lindx2 = lindx2+1


    return time, xlay1, ylay1, rlay1, xlay2, ylay2, rlay2



#################################################

#Now we get the MSD for an ensemble 

def getMSDensemble(X, n, cells):

    xmsd = ny.zeros(n)
    for times in range(0,n):
        for seeds in range(0,cells):
            dx = X[times,seeds] - X[0,seeds]
            xmsd[times] = xmsd[times] + dx*dx/cells

    return xmsd


################################################

#Now we get the MSD by averaging over time lags


def getMSDlags(X, n):

    xmsd = ny.zeros(n)
    for times in range(1,n):   #This is the loop over all lags
        tau = n-times
        for i in range(0,tau):   #This is the loop for each lag
            dX = (X[i+times]-X[i])**2

            xmsd[times] = xmsd[times] + dX/tau


    xmsd = xmsd - xmsd[0]
    return xmsd

######################################################

def getdist(X, n):

    xdist = ny.zeros(n)

    for times in range(1,n):
        dX = abs(X[times]-X[times-1])
        xdist[times] = xdist[times-1]+dX

    return xdist

#####################################################

#We try to get the MSD with a faster method

#def getMSDlagslist(X, n):

#    xmsd = ny.zeros(n)
#    for times in X:   #This is the loop over all lags
#        if times != 
#        tau = n-times
#        for i in range(1,tau):   #This is the loop for each lag
#            dX = (X[i+times]-X[i])*(X[i+times]-X[i])
#            if times>n-3:
#                print("The displacments are")
#                print(X[i+times])
#                print(X[i])
#                print("With magnitude")
#                print(dX)
#            xmsd[times] = xmsd[times] + dX/tau


#    xmsd = xmsd - xmsd[0]
#    return xmsd




#########################################################

def getMSDlagsallcells(X, n,ncells):

    xmsd = ny.zeros((n,ncells))
    for times in range(1,n-1):   #This is the loop over all lags
        tau = n-times
        for i in range(1,tau):   #This is the loop for each lag
            for cells in range(0,ncells):
                dX = (X[i+times,cells]-X[i,cells])*(X[i+times,cells]-X[i,cells])
#            if times>n-3:
#                print("The displacments are")
#                print(X[i+times])
#                print(X[i])
#                print("With magnitude")
#                print(dX)
                xmsd[times,cells] = xmsd[times,cells] + dX/tau

    for cells in range(0,ncells):
        xmsd[:,cells] = xmsd[:,cells] - xmsd[0,cells]

    return xmsd



#################################################
def getMSDlagsallcellsalldim(X1, Y1, X2, Y2,n,ncells):

    x1msd = ny.zeros((n,ncells))
    y1msd = ny.zeros((n,ncells))
    x2msd = ny.zeros((n,ncells))
    y2msd = ny.zeros((n,ncells))

    for times in range(1,n):   #This is the loop over all lags
        tau = n-times
        for i in range(1,tau):   #This is the loop for each lag
            for cells in range(0,ncells):
                dX1 = (X1[i+times,cells]-X1[i,cells])*(X1[i+times,cells]-X1[i,cells])
                dY1 = (Y1[i+times,cells]-Y1[i,cells])*(Y1[i+times,cells]-Y1[i,cells])
                dX2 = (X2[i+times,cells]-X2[i,cells])*(X2[i+times,cells]-X2[i,cells])
                dY2 = (Y2[i+times,cells]-Y2[i,cells])*(Y2[i+times,cells]-Y2[i,cells])

                x1msd[times,cells] = x1msd[times,cells] + dX1/tau
                y1msd[times,cells] = y1msd[times,cells] + dY1/tau
                x2msd[times,cells] = x2msd[times,cells] + dX2/tau
                y2msd[times,cells] = y2msd[times,cells] + dY2/tau


    for cells in range(0,ncells):
        x1msd[:,cells] = x1msd[:,cells] - x1msd[0,cells]
        y1msd[:,cells] = y1msd[:,cells] - y1msd[0,cells]
        x2msd[:,cells] = x2msd[:,cells] - x2msd[0,cells]
        y2msd[:,cells] = y2msd[:,cells] - y2msd[0,cells]



    return x1msd, y1msd, x2msd, y2msd

#####################################################
def getallMSDlags(X1, Y1, X2, Y2,n,ncells,np0,nse,nbi):

    x1msd = ny.zeros((n,ncells,np0,nbi,nse))
    y1msd = ny.zeros((n,ncells,np0,nbi,nse))
    x2msd = ny.zeros((n,ncells,np0,nbi,nse))
    y2msd = ny.zeros((n,ncells,np0,nbi,nse))


    for p0 in range(0,np0):
        for bi in range(0,nbi):
            for se in range(0,nse):

                for cells in ny.arange(0,ncells):

                    X1temp = X1[:,cells,p0,bi,se]
                    Y1temp = Y1[:,cells,p0,bi,se]
                    X2temp = X2[:,cells,p0,bi,se]
                    Y2temp = Y2[:,cells,p0,bi,se]

                    x1msd[:,cells,p0,bi,se] = getMSDlags(X1temp,n)
                    y1msd[:,cells,p0,bi,se] = getMSDlags(Y1temp,n)
                    x2msd[:,cells,p0,bi,se] = getMSDlags(X2temp,n)
                    y2msd[:,cells,p0,bi,se] = getMSDlags(Y2temp,n)



#                for times in ny.arange(1,n):   #This is the loop over all lags
#                    tau = n-times
#                    for i in ny.arange(0,tau):   #This is the loop for each lag
#                        for cells in ny.arange(0,ncells):
#                            dX1 = (X1[i+times,cells,p0,bi,se]-X1[i,cells,p0,bi,se])**2
#                            dY1 = (Y1[i+times,cells,p0,bi,se]-Y1[i,cells,p0,bi,se])**2
#                            dX2 = (X2[i+times,cells,p0,bi,se]-X2[i,cells,p0,bi,se])**2
#                            dY2 = (Y2[i+times,cells,p0,bi,se]-Y2[i,cells,p0,bi,se])**2

#                            x1msd[times,cells,p0,bi,se] += dX1/tau
#                            y1msd[times,cells,p0,bi,se] += dY1/tau
#                            x2msd[times,cells,p0,bi,se] += dX2/tau
#                            y2msd[times,cells,p0,bi,se] += dY2/tau


#                for cells in range(0,ncells):
#                    x1msd[:,cells,p0,bi,se] = x1msd[:,cells,p0,bi,se] - x1msd[0,cells,p0,bi,se]
#                    y1msd[:,cells,p0,bi,se] = y1msd[:,cells,p0,bi,se] - y1msd[0,cells,p0,bi,se]
#                    x2msd[:,cells,p0,bi,se] = x2msd[:,cells,p0,bi,se] - x2msd[0,cells,p0,bi,se]
#                    y2msd[:,cells,p0,bi,se] = y2msd[:,cells,p0,bi,se] - y2msd[0,cells,p0,bi,se]



    return x1msd, y1msd, x2msd, y2msd

##########################################
#Here we get the MSD at rise time
def getallMSDlagsatrise(X1, Y1, X2, Y2, CT, ncells,np0,nse,nbi,rtim):

    x1msd = [[[[None for _ in range(ncells) ] for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)] 
    y1msd = [[[[None for _ in range(ncells) ] for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)] 
    x2msd = [[[[None for _ in range(ncells) ] for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)] 
    y2msd = [[[[None for _ in range(ncells) ] for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)] 

    rmsd = [[[[None for _ in range(ncells) ] for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)] 

    time = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)] 


    for p0 in range(0,np0):
        for bi in range(0,nbi):
            for se in range(0,nse):
#                for cells in range(0,ncells):
#                    x1temp = list(X1[:,cells,p0,bi,se])
#                    for times in ny.arange(1,n).tolist():
#                        tau = n-times
#                        for dt in x1temp:
#                            dX1 = x1temp[dt+i

                Rtime = math.ceil(rtim.edgerise[p0,bi,se]/10)

                time[p0][bi][se] = CT.timedat[:Rtime,p0,bi,se]


                for cells in range(0,ncells):

                    X1temp = []
                    Y1temp = []
                    X2temp = []
                    Y2temp = []

                    if Rtime > 1:

                        X1temp = X1[:Rtime,cells,p0,bi,se]
                        Y1temp = Y1[:Rtime,cells,p0,bi,se]
                        X2temp = X2[:Rtime,cells,p0,bi,se]
                        Y2temp = Y2[:Rtime,cells,p0,bi,se]

                        x1msd[p0][bi][se][cells] = getMSDlags(X1temp,Rtime)
                        y1msd[p0][bi][se][cells] = getMSDlags(Y1temp,Rtime)
                        x2msd[p0][bi][se][cells] = getMSDlags(X2temp,Rtime)
                        y2msd[p0][bi][se][cells] = getMSDlags(Y2temp,Rtime)

                        rtemp = (ny.array(x1msd[p0][bi][se][cells]) + ny.array(x2msd[p0][bi][se][cells]))/2 + (ny.array(y1msd[p0][bi][se][cells]) + ny.array(y2msd[p0][bi][se][cells]))/2

                        rmsd[p0][bi][se][cells] = list(rtemp)


#                for times in ny.arange(1,Rtime):   #This is the loop over all lags
#                    tau = Rtime-times
#                    for i in ny.arange(0,tau):   #This is the loop for each lag
#                        for cells in ny.arange(0,ncells):
#                            dX1 = (X1[i+times,cells,p0,bi,se]-X1[i,cells,p0,bi,se])**2
#                            dY1 = (Y1[i+times,cells,p0,bi,se]-Y1[i,cells,p0,bi,se])**2
#                            dX2 = (X2[i+times,cells,p0,bi,se]-X2[i,cells,p0,bi,se])**2
#                            dY2 = (Y2[i+times,cells,p0,bi,se]-Y2[i,cells,p0,bi,se])**2

#                            x1msd[times,cells,p0,bi,se] += dX1/tau
#                            y1msd[times,cells,p0,bi,se] += dY1/tau
#                            x2msd[times,cells,p0,bi,se] += dX2/tau
#                            y2msd[times,cells,p0,bi,se] += dY2/tau


#                for cells in range(0,ncells):
#                    x1msd[:,cells,p0,bi,se] = x1msd[:,cells,p0,bi,se] - x1msd[0,cells,p0,bi,se]
#                    y1msd[:,cells,p0,bi,se] = y1msd[:,cells,p0,bi,se] - y1msd[0,cells,p0,bi,se]
#                    x2msd[:,cells,p0,bi,se] = x2msd[:,cells,p0,bi,se] - x2msd[0,cells,p0,bi,se]
#                    y2msd[:,cells,p0,bi,se] = y2msd[:,cells,p0,bi,se] - y2msd[0,cells,p0,bi,se]



    return x1msd, y1msd, x2msd, y2msd, rmsd, time

##########################################
#Here we get the MSD after rise time
def getallMSDlagsafterrise(X1, Y1, X2, Y2, CT, ncells,np0,nse,nbi,rtim):

    x1msd = [[[[None for _ in range(ncells) ] for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)] 
    y1msd = [[[[None for _ in range(ncells) ] for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)] 
    x2msd = [[[[None for _ in range(ncells) ] for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)] 
    y2msd = [[[[None for _ in range(ncells) ] for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)] 

    rmsd = [[[[None for _ in range(ncells) ] for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)] 

    time = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)] 


    for p0 in range(0,np0):
        for bi in range(0,nbi):
            for se in range(0,nse):
#                for cells in range(0,ncells):
#                    x1temp = list(X1[:,cells,p0,bi,se])
#                    for times in ny.arange(1,n).tolist():
#                        tau = n-times
#                        for dt in x1temp:
#                            dX1 = x1temp[dt+i

                Rtime = math.ceil(rtim.edgerise[p0,bi,se]/10)

                time[p0][bi][se] = CT.timedat[Rtime:,p0,bi,se]

                ftime = len(CT.timedat[Rtime:,p0,bi,se])


                for cells in range(0,ncells):

                    X1temp = []
                    Y1temp = []
                    X2temp = []
                    Y2temp = []

                    if Rtime > 1 and Rtime < ftime:

                        X1temp = X1[Rtime:,cells,p0,bi,se]
                        Y1temp = Y1[Rtime:,cells,p0,bi,se]
                        X2temp = X2[Rtime:,cells,p0,bi,se]
                        Y2temp = Y2[Rtime:,cells,p0,bi,se]

                        x1msd[p0][bi][se][cells] = getMSDlags(X1temp,ftime)
                        y1msd[p0][bi][se][cells] = getMSDlags(Y1temp,ftime)
                        x2msd[p0][bi][se][cells] = getMSDlags(X2temp,ftime)
                        y2msd[p0][bi][se][cells] = getMSDlags(Y2temp,ftime)

                        rtemp = (ny.array(x1msd[p0][bi][se][cells]) + ny.array(x2msd[p0][bi][se][cells]))/2 + (ny.array(y1msd[p0][bi][se][cells]) + ny.array(y2msd[p0][bi][se][cells]))/2

                        rmsd[p0][bi][se][cells] = list(rtemp)


#                for times in ny.arange(1,Rtime):   #This is the loop over all lags
#                    tau = Rtime-times
#                    for i in ny.arange(0,tau):   #This is the loop for each lag
#                        for cells in ny.arange(0,ncells):
#                            dX1 = (X1[i+times,cells,p0,bi,se]-X1[i,cells,p0,bi,se])**2
#                            dY1 = (Y1[i+times,cells,p0,bi,se]-Y1[i,cells,p0,bi,se])**2
#                            dX2 = (X2[i+times,cells,p0,bi,se]-X2[i,cells,p0,bi,se])**2
#                            dY2 = (Y2[i+times,cells,p0,bi,se]-Y2[i,cells,p0,bi,se])**2

#                            x1msd[times,cells,p0,bi,se] += dX1/tau
#                            y1msd[times,cells,p0,bi,se] += dY1/tau
#                            x2msd[times,cells,p0,bi,se] += dX2/tau
#                            y2msd[times,cells,p0,bi,se] += dY2/tau


#                for cells in range(0,ncells):
#                    x1msd[:,cells,p0,bi,se] = x1msd[:,cells,p0,bi,se] - x1msd[0,cells,p0,bi,se]
#                    y1msd[:,cells,p0,bi,se] = y1msd[:,cells,p0,bi,se] - y1msd[0,cells,p0,bi,se]
#                    x2msd[:,cells,p0,bi,se] = x2msd[:,cells,p0,bi,se] - x2msd[0,cells,p0,bi,se]
#                    y2msd[:,cells,p0,bi,se] = y2msd[:,cells,p0,bi,se] - y2msd[0,cells,p0,bi,se]



    return x1msd, y1msd, x2msd, y2msd, rmsd, time







########################################
def getalldist(X1, Y1, X2, Y2,n,ncells,np0,nse,nbi):

    x1dist = ny.zeros((n,ncells,np0,nbi,nse))
    y1dist = ny.zeros((n,ncells,np0,nbi,nse))
    x2dist = ny.zeros((n,ncells,np0,nbi,nse))
    y2dist = ny.zeros((n,ncells,np0,nbi,nse))


    for p0 in range(0,np0):
        for bi in range(0,nbi):
            for se in range(0,nse):

                for cells in range(0,ncells):

                    x1dist[:,cells,p0,bi,se] = getdist(X1[:,cells,p0,bi,se],n)
                    y1dist[:,cells,p0,bi,se] = getdist(Y1[:,cells,p0,bi,se],n)
                    x2dist[:,cells,p0,bi,se] = getdist(X2[:,cells,p0,bi,se],n)
                    y2dist[:,cells,p0,bi,se] = getdist(Y2[:,cells,p0,bi,se],n)


#                for times in range(1,n):   
#                        for cells in range(0,ncells):
#                            dX1 = abs(X1[times,cells,p0,bi,se]-X1[times-1,cells,p0,bi,se])
#                            dY1 = abs(Y1[times,cells,p0,bi,se]-Y1[times-1,cells,p0,bi,se])
#                            dX2 = abs(X2[times,cells,p0,bi,se]-X2[times-1,cells,p0,bi,se])
#                            dY2 = abs(Y2[times,cells,p0,bi,se]-Y2[times-1,cells,p0,bi,se])

#                            x1dist[times,cells,p0,bi,se] = x1dist[times-1,cells,p0,bi,se] + dX1
#                            y1dist[times,cells,p0,bi,se] = y1dist[times-1,cells,p0,bi,se] + dY1
#                            x2dist[times,cells,p0,bi,se] = x2dist[times-1,cells,p0,bi,se] + dX2
#                            y2dist[times,cells,p0,bi,se] = y2dist[times-1,cells,p0,bi,se] + dY2



    return x1dist, y1dist, x2dist, y2dist




#########################################

def getdistallcellsalldim(X1, Y1, X2, Y2,n,ncells):

    x1dist = ny.zeros((n,ncells))
    y1dist = ny.zeros((n,ncells))
    x2dist = ny.zeros((n,ncells))
    y2dist = ny.zeros((n,ncells))

    for times in range(1,n):
        for cells in range(0,ncells):
            dX1 = abs(X1[times,cells]-X1[times-1,cells])
            dY1 = abs(Y1[times,cells]-Y1[times-1,cells])
            dX2 = abs(X2[times,cells]-X2[times-1,cells])
            dY2 = abs(Y2[times,cells]-Y2[times-1,cells])

            x1dist[times,cells] = x1dist[times-1,cells] + dX1
            y1dist[times,cells] = y1dist[times-1,cells] + dY1
            x2dist[times,cells] = x2dist[times-1,cells] + dX2
            y2dist[times,cells] = y2dist[times-1,cells] + dY2




    return x1dist, y1dist, x2dist, y2dist


##################################################
def getdispallcellsalldim(X1, Y1, X2, Y2,n,ncells):

    x1disp = ny.zeros((n,ncells))
    y1disp = ny.zeros((n,ncells))
    x2disp = ny.zeros((n,ncells))
    y2disp = ny.zeros((n,ncells))

#    for times in range(1,n):
#        for cells in range(0,ncells):
#            dX1 = X1[times,cells]-X1[times-1,cells]
#            dY1 = Y1[times,cells]-Y1[times-1,cells]
#            dX2 = X2[times,cells]-X2[times-1,cells]
#            dY2 = Y2[times,cells]-Y2[times-1,cells]

            #x1disp[times,cells] = x1disp[times-1,cells] + dX1
            #y1disp[times,cells] = y1disp[times-1,cells] + dY1
            #x2disp[times,cells] = x2disp[times-1,cells] + dX2
            #y2disp[times,cells] = y2disp[times-1,cells] + dY2

    for cells in range(0,ncells):

        x1disp[:,cells] = X1[:,cells] - X1[0,cells]
        y1disp[:,cells] = Y1[:,cells] - Y1[0,cells]
        x2disp[:,cells] = X2[:,cells] - X2[0,cells]
        y2disp[:,cells] = Y2[:,cells] - Y2[0,cells]


    return x1disp, y1disp, x2disp, y2disp



###############################################
#We ignore the cells with spikes

def getMSDlagnospikes(X, Xmsd,n,ncells,xmin, xmax):

    Lx = xmax-xmin
    xmsdavg = ny.zeros(n)
    ckeep = 0
    for cells in range(0,ncells):
        diffX = ny.diff(X[:,cells])
        if ny.max(diffX) < Lx/4 and ny.min(diffX) > -Lx/4:
            xmsdavg = xmsdavg + Xmsd[:,cells]
            ckeep = ckeep + 1

    if ckeep > 0:
        xmsdavg = xmsdavg/ckeep

    return xmsdavg, ckeep

#####################################################

#Let's also get the distance walked

def getabsdisp(X, n):

    xdisp = ny.zeros(n)
    for times in range(1,n):   #This is the loop over all lags
        dX = abs(X[times]-X[times-1])
        xdisp[times] = xdisp[times-1] + dX


    return xdisp


#####################################################

def getabsdispallcells(X,n,ncells):
    
    xdisp = ny.zeros((n,ncells))
    dX = ny.zeros((n-1,ncells))
    for cells in range(0,ncells):
        dX[:,cells] = ny.diff(X[:,cells])
    dXabs = ny.abs(dX)
    for times in range(1,n):
        for cells in range(0,ncells):
            #dX = abs(X[times,cells]-X[times-1,cells])
            xdisp[times,cells] = xdisp[times-1,cells] + dXabs[times-1,cells]

    return xdisp


##########################################
def rmsdavgovercells(X1, Y1, X2, Y2, chkx1, chky1, chkx2, chky2,n,ncells):

    xmsd = ny.zeros((n))
    ymsd = ny.zeros((n))
    rmsd = ny.zeros((n))

    Nx = 0
    Ny = 0

    for cells in range(0,ncells):
        if chkx1[cells] == 1:
            xmsd = xmsd + X1[:,cells]
            Nx = Nx+1
        if chky1[cells] == 1:
            ymsd = ymsd + Y1[:,cells]
            Ny = Ny + 1
        if chkx2[cells] == 1:
            xmsd = xmsd + X2[:,cells]
            Nx = Nx+1
        if chky2[cells] == 1:
            ymsd = ymsd + Y2[:,cells]
            Ny = Ny+1

    if Nx > 0:
        xmsd = xmsd/float(Nx)
    if Nx == 0:
        print("No usable cells in X")
    if Ny > 0:
        ymsd = ymsd/float(Ny)
    if Ny == 0:
        print("No usable cells in Y")
    if Nx == 0 and Ny == 0:
        print("No usable cells at all")

    rmsd = xmsd + ymsd 
 
    return rmsd

#######################################################
def rdistavgovercells(X1, Y1, X2, Y2, chkx1, chky1, chkx2, chky2,n,ncells):

    xdist = ny.zeros((n))
    ydist = ny.zeros((n))
    rdist = ny.zeros((n))

    Nx = 0
    Ny = 0

    for cells in range(0,ncells):
        if chkx1[cells] == 1:
            xdist = xdist + X1[:,cells]
            Nx = Nx+1
        if chky1[cells] == 1:
            ydist = ydist + Y1[:,cells]
            Ny = Ny + 1
        if chkx2[cells] == 1:
            xdist = xdist + X2[:,cells]
            Nx = Nx+1
        if chky2[cells] == 1:
            ydist = ydist + Y2[:,cells]
            Ny = Ny+1

    if Nx > 0:
        xdist = xdist/float(Nx)
    if Ny > 0:
        ydist = ydist/float(Ny)

    rdist = xdist **2 + ydist **2
    rdist = ny.sqrt(rdist)
 
    return rdist


###############################################
def rmsdavgoverseeds(X1, Y1, X2, Y2, chkx1, chky1, chkx2, chky2,n,nsee):

    xmsd = ny.zeros((n))
    ymsd = ny.zeros((n))
    rmsd = ny.zeros((n))

    Nx = 0
    Ny = 0

    for seeds in range(0,nsee):
        if ny.sum(chkx1[seeds]) > 0:
            xmsd = xmsd + X1[:,seeds]
            Nx = Nx+1
        if ny.sum(chky1[seeds]) > 0:
            ymsd = ymsd + Y1[:,seeds]
            Ny = Ny + 1
        if ny.sum(chkx2[seeds]) > 0:
            xmsd = xmsd + X2[:,seeds]
            Nx = Nx+1
        if ny.sum(chky2[seeds]) > 0:
            ymsd = ymsd + Y2[:,seeds]
            Ny = Ny+1

    if Nx > 0:
        xmsd = xmsd/float(Nx)
    if Ny > 0:
        ymsd = ymsd/float(Ny)

    rmsd = xmsd + ymsd 
 
    return rmsd

#######################################################
def rdistavgoverseeds(X1, Y1, X2, Y2, chkx1, chky1, chkx2, chky2,n,nsee):

    xdist = ny.zeros((n))
    ydist = ny.zeros((n))
    rdist = ny.zeros((n))

    Nx = 0
    Ny = 0

    for seeds in range(0,nsee):
        if ny.sum(chkx1[seeds]) > 0:
            xdist = xdist + X1[:,seeds]
            Nx = Nx+1
        if ny.sum(chky1[cells]) > 0:
            ydist = ydist + Y1[:,cells]
            Ny = Ny + 1
        if ny.sum(chkx2[cells]) > 0:
            xdist = xdist + X2[:,cells]
            Nx = Nx+1
        if ny.sum(chky2[cells]) > 0:
            ydist = ydist + Y2[:,cells]
            Ny = Ny+1

    if Nx > 0:
        xdist = xdist/float(Nx)
    if Ny > 0:
        ydist = ydist/float(Ny)

    rdist = xdist ** 2 + ydist ** 2
    rdist = ny.sqrt(rdist)
 
    return rdist







####################################################
def getabsdistbulk(X1,Y1,X2,Y2,n,ncells,nshape,nseeds,nbi):
    
    xdist1 = ny.zeros((n,ncells,nshape,nseeds,nbi))

    ydist1 = ny.zeros((n,ncells,nshape,nseeds,nbi))

    xdist2 = ny.zeros((n,ncells,nshape,nseeds,nbi))

    ydist2 = ny.zeros((n,ncells,nshape,nseeds,nbi))



    for lines in range(1,n):
        for shape in range (0,nshape):
            for seed in range(0,nseeds):
                for bi in range(0,nbi):
                    for cells in range(0,ncells):
                        xdist1[lines,cells,shape,seed,bi] = xdist1[lines-1,cells,shape,seed,bi] + abs(X1[lines,cells,shape,seed,bi] - X1[lines-1,cells,shape,seed,bi])

                        ydist1[lines,cells,shape,seed,bi] = ydist1[lines-1,cells,shape,seed,bi] + abs(Y1[lines,cells,shape,seed,bi] - Y1[lines-1,cells,shape,seed,bi])

                        xdist2[lines,cells,shape,seed,bi] = xdist2[lines-1,cells,shape,seed,bi] + abs(X2[lines,cells,shape,seed,bi] - X2[lines-1,cells,shape,seed,bi])

                        ydist2[lines,cells,shape,seed,bi] = ydist2[lines-1,cells,shape,seed,bi] + abs(Y2[lines,cells,shape,seed,bi] - Y2[lines-1,cells,shape,seed,bi])



    return xdist1, ydist1, xdist2, ydist2 

###############################################
def getdistanddispbulk(X1,Y1,X2,Y2,n,ncells,nshape,nseeds,nbi):
    
    xdist1 = ny.zeros((n,ncells,nshape,nseeds,nbi))

    ydist1 = ny.zeros((n,ncells,nshape,nseeds,nbi))

    xdist2 = ny.zeros((n,ncells,nshape,nseeds,nbi))

    ydist2 = ny.zeros((n,ncells,nshape,nseeds,nbi))

    xdisp1 = ny.zeros((n,ncells,nshape,nseeds,nbi))

    ydisp1 = ny.zeros((n,ncells,nshape,nseeds,nbi))

    xdisp2 = ny.zeros((n,ncells,nshape,nseeds,nbi))

    ydisp2 = ny.zeros((n,ncells,nshape,nseeds,nbi))



    for lines in range(1,n):
        for shape in range (0,nshape):
            for seed in range(0,nseeds):
                for bi in range(0,nbi):
                    for cells in range(0,ncells):
                        xdist1[lines,cells,shape,seed,bi] = xdisp1[lines-1,cells,shape,seed,bi] + abs(X1[lines,cells,shape,seed,bi] - X1[lines-1,cells,shape,seed,bi])
                        ydist1[lines,cells,shape,seed,bi] = ydisp1[lines-1,cells,shape,seed,bi] + abs(Y1[lines,cells,shape,seed,bi] - Y1[lines-1,cells,shape,seed,bi])

                        xdist2[lines,cells,shape,seed,bi] = xdisp2[lines-1,cells,shape,seed,bi] + abs(X2[lines,cells,shape,seed,bi] - X2[lines-1,cells,shape,seed,bi])

                        ydist2[lines,cells,shape,seed,bi] = ydisp2[lines-1,cells,shape,seed,bi] + abs(Y2[lines,cells,shape,seed,bi] - Y2[lines-1,cells,shape,seed,bi])

                        xdisp1[lines,cells,shape,seed,bi] = xdisp1[lines-1,cells,shape,seed,bi] + (X1[lines,cells,shape,seed,bi] - X1[lines-1,cells,shape,seed,bi])
                        ydisp1[lines,cells,shape,seed,bi] = ydisp1[lines-1,cells,shape,seed,bi] + (Y1[lines,cells,shape,seed,bi] - Y1[lines-1,cells,shape,seed,bi])

                        xdisp2[lines,cells,shape,seed,bi] = xdisp2[lines-1,cells,shape,seed,bi] + (X2[lines,cells,shape,seed,bi] - X2[lines-1,cells,shape,seed,bi])

                        ydisp2[lines,cells,shape,seed,bi] = ydisp2[lines-1,cells,shape,seed,bi] + (Y2[lines,cells,shape,seed,bi] - Y2[lines-1,cells,shape,seed,bi])



    return xdist1, ydist1, xdist2, ydist2, xdisp1, ydisp1, xdisp2, ydisp2 


################################################

def getabsdistnospike(X,n,ncells,xmin, xmax):
    
    Lx = xmax-xmin
    xdist = ny.zeros(n)
    ckeep = 0
    toadd = ny.zeros(ncells)
    for cells in range(0,ncells):
        diffX = ny.diff(X[:,cells])
        if ny.max(diffX) < Lx/4 and ny.min(diffX) > -Lx/4:
            toadd[cells] = 1
            ckeep = ckeep + 1

    for li in range (1, n):
        for cells in range(0,ncells):
            if toadd[cells] == 1:
                xdist[li] = xdist[li] + abs(X[li,cells] - X[li-1,cells])

    if ckeep > 0:
        xdist = xdist/float(ckeep)


    return xdist



################################################

def getabsdistnospike2(X,Xd,n,ncells,xmin, xmax):
    
    Lx = xmax-xmin
    xdist = ny.zeros(n)
    ckeep = 0
    for cells in range(0,ncells):
        diffX = ny.diff(X[:,cells])
        if ny.max(diffX) < Lx/4 and ny.min(diffX) > -Lx/4:
            xdist = xdist + Xd[:,cells]
            ckeep = ckeep + 1

    if ckeep > 0:
        xdist = xdist/float(ckeep)


    if ckeep == 0:
        print("No reasonable cells found")

    return xdist



#################################################

def getabsdispnospike(X,Xd,n,ncells,xmin, xmax):
    
    Lx = xmax-xmin
    xdisp = ny.zeros(n)
    ckeep = 0
    for cells in range(0,ncells):
        diffX = ny.diff(X[:,cells])
        if ny.max(abs(diffX)) < Lx/4 and ny.min(diffX) > -Lx/4:
            xdisp = xdisp + abs(Xd[:,cells])
            ckeep = ckeep + 1

    if ckeep > 0:
        xdisp = xdisp/float(ckeep)

    if ckeep == 0:
        print("No reasonable cells found")

    return xdisp


#####################################################

def getmsdanddistnospike(X,Xsd,Xst,n,ncells,xmin, xmax):
    
    Lx = xmax-xmin
    xmsd = ny.zeros(n)
    xdist = ny.zeros(n)
    ckeep = 0
    keeplist = ny.zeros(ncells)
    for cells in range(0,ncells):
        diffX = ny.diff(X[:,cells])
        if ny.max(diffX) < Lx/4 and ny.min(diffX) > -Lx/4:
            xmsd = xmsd + Xsd[:,cells]
            xdist = xdist + Xst[:,cells]
            keeplist[cells] = 1
            ckeep = ckeep + 1


    if ckeep > 0:
        xmsd = xmsd/float(ckeep)
        xdist = xdist/float(ckeep)


    return xmsd, xdist, keeplist  



###################################################

def getMSDR(X,Y,n):

    rmsd = ny.zeros(n)
    for t in range(0,n):
        rmsd[t] = rmsd[t] + X[t] + Y[t]


    return rmsd


#######################################################

