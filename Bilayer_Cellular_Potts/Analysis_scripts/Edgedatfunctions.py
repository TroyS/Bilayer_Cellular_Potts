#This will contain functions for data analysis

import numpy as ny
import statistics
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import plotly.graph_objects as go
from matplotlib.colors import TABLEAU_COLORS, same_color
from matplotlib.pyplot import cm
from matplotlib import cm
import math

###################################################
#Here we get the number of lines in the edgematch file
#We subtract 3 lines to account for the header
def getlines(fle):
    with open(fle, 'r') as gfile:
        glines = gfile.readlines()
        nlines = len(glines)-3

    return nlines

#####################################################
#Here we get the number of lines in the neighbor file

def getlinesneigh(fle):
    with open(fle, 'r') as nfile:
        nlines = nfile.readlines()
        totlines = len(nlines)-1

    return totlines

####################################################

#Here we get the timestep and edgematch ratio data from a file
def getdata(fle,n):

    ts = ny.zeros(n)
    Pm = ny.zeros(n)
    Bim = ny.zeros(n)
    PBrat = ny.zeros(n)

    tsPB = ny.zeros((n,2))

    with open(fle, 'r') as dfile:
        dlines = dfile.readlines()

        for i1 in range(3, n+3):
            ddata = dlines[i1].split();

            ts[i1-3] = float(ddata[0])
            Pm[i1-3] = float(ddata[1])
            Bim[i1-3] = float(ddata[2])
            PBrat[i1-3] = float(ddata[3])

    tsPB[:,0] = ts
    tsPB[:,1] = PBrat

    return ny.array(tsPB)

#######################################
#We get all the data for a configuration at once

def getalldata(n, np0,nbi,nse, namep0,namebi,namesee,name):

    ts = ny.zeros((n,np0, nbi, nse))
    Pm = ny.zeros((n, np0, nbi, nse))
    Bim = ny.zeros((n, np0, nbi, nse))
    PBrat = ny.zeros((n, np0, nbi, nse))

    diff = ny.zeros((n,np0,nbi,nse))

    tsPB = ny.zeros((n,2))

    Namu=[[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]


    for p0 in ny.arange(0,np0):
        for bi in ny.arange(0, nbi):
            for se in ny.arange(0, nse):
                filenme = namep0[p0] + namesee[se] + "Bi_" + namebi[bi] + name

                Namu[p0][bi][se] = filenme

                with open(filenme, 'r') as dfile:
                    dlines = dfile.readlines()

                    for i1 in ny.arange(3, n+3):
                        ddata = dlines[i1].split();

                        ts[i1-3,p0,bi,se] = float(ddata[0])
                        Pm[i1-3,p0,bi,se] = float(ddata[1])
                        Bim[i1-3,p0,bi,se] = float(ddata[2])
                        PBrat[i1-3,p0,bi,se] = float(ddata[3])

                diff[1:,p0,bi,se] = ny.diff(PBrat[:,p0,bi,se])

    tsPB = ny.array(ts)
    PB = ny.array(PBrat)
    Dif = ny.array(diff)


    return ts, PB, Dif, Namu

###########################################
#We get the bending extention times and match ratios in each region

def getallextandmatch(nli,np0,nse,nbi,namep0,namesee,namebi,name):

    x, y, z = np0, nbi, nse

#    extvals = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]

    oldreg1 = 0
    oldreg2 = 0
    extnum = 0

#    exttimes = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]

    numex = [None for _ in range(np0) ]

    extvals = [None for _ in range(np0)]

    names = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]


    r1match = ny.zeros((np0,nbi,nse,nli-1))
    r2match = ny.zeros((np0,nbi,nse,nli-1))
    timematch = ny.zeros(nli-1)

    for p0 in range (0,np0):
        oldlen = 0
        for bi in range(0, nbi):
            valslist = []
            numext = []
            for se in range(0, nse):

                filenme = namep0[p0] + namesee[se] + "Bi_" + namebi[bi] + name

                timeslist = []

                with open(filenme, 'r') as bfile:
                    blines = bfile.readlines()

                    sdata = blines[2].split()

                    start = float(sdata[1]) 

                    timeold = 0
                    

                    for i1 in range(2, nli+1):
                        
                        bdata = blines[i1].split();
                        time = int(bdata[0])

                        chk = 0

                        r1match[p0,bi,se,i1-2] = float(bdata[3])
                        r2match[p0,bi,se,i1-2] = float(bdata[4])

                        if p0 == 0 and bi == 0 and se == 0:
                            timematch[i1-2] = float(bdata[0])

                        if i1 > 2:
                            num1 = float(bdata[1])
                            num2 = float(bdata[2])
        

                            if num1 > oldregI:

                                valslist.append(num1)

                                oldregI = num1

                                timeold = time
                

                        else:
                            oldregI = float(bdata[1]) 
                            valslist.append(oldregI)

                            
                        #benspeed.append(speedlist)
                    names[p0][bi][se] = filenme

        uniqvals = list(set(valslist))
        uniqvals.sort()
        extvals[p0] = uniqvals
        if len(uniqvals) > oldlen:
            oldlen = len(uniqvals)
        exlab = ny.linspace(1,oldlen-1,oldlen-1)
        numex[p0] = list(exlab)


    return extvals, numex, names, r1match, r2match, timematch


##########################################

#This is a function for binning data

def setupbins(X,lowbound,highbound,numbins):

    delx = (highbound-lowbound)/float(numbins)
    Bins = ny.zeros(numbins)

    for ele in X.tolist():
        check = 0
        for j in range(0,numbins):
            dx1 = lowbound + j*delx
            dx2 = lowbound + (j+1)*delx
            if check == 0:
                if ele < dx2 and ele >= dx1:
                    Bins[j]+= 1
                    check = 1

    return Bins

###############################################
#Then we get bins for every seed

def getallbins(X,lowbound,highbound,numbins,np0,nbi,nse,nli):


    delx = (highbound-lowbound)/float(numbins)
    Bins = ny.linspace(lowbound,highbound,numbins)
    Freq = ny.zeros((np0,nbi,nse,numbins))
    Freqseed = ny.zeros((np0,nbi,numbins))
    Localminseed = ny.zeros((np0,nbi,nse,numbins))
    Localmin = ny.zeros((np0,nbi,numbins))
    Numengymins = ny.zeros((np0,nbi,nse))
    Modeengymins = ny.zeros((np0,nbi))
    Avgnumengymins = ny.zeros((np0,nbi))
    Maxengyminseed = ny.zeros((np0,nbi,nse))
    Maxengymin = ny.zeros((np0,nbi))
    Maxengyprob = ny.zeros((np0,nbi))
    Engyminvals = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]


    for p0 in range(0,np0):
        for bi in range(0,nbi):
            maxindx = 0
            for se in range(0,nse):
                maxindxseed = 0
                Freq[p0,bi,se,:] = setupbins(X[:,p0,bi,se],lowbound,highbound,numbins)
                Freqseed[p0,bi,:] += Freq[p0,bi,se,:]


############
#Then we get the metastable states
                eva =[]
                for bins in range(0,numbins):
                    if Freq[p0,bi,se,bins] >= nli/float(10):
                        Numengymins[p0,bi,se]+=1
                        Localminseed[p0,bi,se,bins] = 1
                        Localmin[p0,bi,bins] = 1
                        eva.append(bins*delx)
                        if bins > maxindxseed:
                            maxindxseed = bins
                        if bins > maxindx:
                            maxindx = bins


                Maxengyminseed[p0,bi,se] = maxindxseed*delx
                Engyminvals[p0][bi][se] = eva



            Maxengymin[p0,bi] = maxindx*delx
            modeindx = 0
            maxval = 0

            for bins in range(0,numbins):
                if Freqseed[p0,bi,bins] >= (nli/float(10*nse)) and Freqseed[p0,bi,bins] > maxval:
                    maxval = Freqseed[p0,bi,bins]
                    modeindx = bins

            Maxengymin[p0,bi] = maxindx*delx
            Avgnumengymins[p0,bi] = ny.mean(Numengymins[p0,bi,:])
            Modeengymins[p0,bi] = modeindx*delx

##########
#Now we get the probabilities

    for p0 in range(0,np0):
        for bi in range(0,nbi):
            numreach = 0
            for se in range(0,nse):
                if Maxengyminseed[p0,bi,se] == Maxengymin[p0,bi]:
                    numreach +=1

            Maxengyprob[p0,bi] = numreach/float(nse) 



    return Bins, Freq, Freqseed, Localminseed, Localmin, Avgnumengymins, Maxengymin, Modeengymins, Maxengyprob, Engyminvals 


##############################################

#We average over seeds

def edgeavgoverseed(EM,n, np0,nbi):

    
    PBavg = ny.zeros((n, np0, nbi))
    PBstd = ny.zeros((n, np0, nbi))
    time = ny.zeros((n,np0,nbi))

    for lines in ny.arange(0,n):
        for p0 in ny.arange(0,np0):
            for bi in ny.arange(0, nbi):
                PBavg[lines,p0,bi] = ny.mean(EM.edgedat[lines,p0,bi,:])
                PBstd[lines,p0,bi] = ny.std(EM.edgedat[lines,p0,bi,:])
                time[lines,p0,bi] = EM.tsdat[lines,p0,bi,0]


    return time,PBavg,PBstd

######################################
def getbinsforalledge(dat,n, np0,nbi, nse, lowbound,highbound,interval):

    delx = highbound-lowbound/interval
    nelements = math.ceil((highbound-lowbound)/interval)
 
    Freq = ny.zeros((nelements,np0,nbi,nse)) 

    for p0 in ny.arange(0,np0):
        for bi in ny.arange(0, nbi):
            for se in ny.arange(0,nse):
                Freq[:,p0,bi,se] = setupbins(dat.edgeinfo[:,p0,bi,se],lowbound,highbound,interval,n)
                


    return Freq



#########################################

#We get the neighbor changes

def getneighdata(fle,n):

    neichg1 = ny.zeros(n)
    neichg2 = ny.zeros(n)

    with open(fle, 'r') as nfile:
        nlines = nfile.readlines()

        for i1 in range(1, n+1):
            ndata = nlines[i1].split();

            neichg1[i1-1] = float(ndata[1])
            neichg2[i1-1] = float(ndata[2])


    return neichg1, neichg2



##########################################

def getallneighdata(n,np0,nse,nbi,namep0,namesee,namebi,name):


    Time = ny.zeros((n,np0,nbi,nse))
    neichg1 = ny.zeros((n,np0,nbi,nse))
    neichg2 = ny.zeros((n,np0,nbi,nse))

    neichg12 = ny.zeros((n,np0,nbi,nse)) 

    neichgsum = ny.zeros((n,np0,nbi))

    neichgtotseed = ny.zeros((np0,nbi,nse))
    neichgtot = ny.zeros((np0,nbi))

    neichgavg = ny.zeros((np0,nbi))
    neichgstd = ny.zeros((np0,nbi))

    Namu=[[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]



    for p0 in ny.arange(0,np0):
        for bi in ny.arange(0,nbi):
            for se in ny.arange(0,nse):

                filenme = namep0[p0] + namesee[se] + "Bi_" + namebi[bi] + name

                Namu[p0][bi][se] = filenme

                with open(filenme, 'r') as nfile:
                    nlines = nfile.readlines()

                    for i1 in ny.arange(1, n+1):
                        ndata = nlines[i1].split();
                        Time[i1-1,p0,bi,se] = float(ndata[0])
                        neichg1[i1-1,p0,bi,se] = float(ndata[1])
                        neichg2[i1-1,p0,bi,se] = float(ndata[2])
                        neichg12[i1-1,p0,bi,se]=float(ndata[1])+float(ndata[2])
                        neichgsum[i1-1,p0,bi]+=float(ndata[1])+float(ndata[2])

                neichgtotseed[p0,bi,se] = ny.sum(neichg12[:,p0,bi,se])

            neichgtot[p0,bi] = ny.sum(neichgtotseed[p0,bi,:])
            neichgavg[p0,bi] = ny.mean(neichgtotseed[p0,bi,:])
            neichgstd[p0,bi] = ny.std(neichgtotseed[p0,bi,:])


    return Time,neichg1, neichg2, neichg12, neichgsum,neichgtotseed,neichgtot,neichgavg, neichgstd, Namu



#######################################################
#We check the probability of a neighbor change occuring when a system leaves a metastable state

def neichangeatjumpprob(EM,Nei,Mstab,n,np0,nse,nbi):

    neiprobseed = ny.zeros((np0,nbi,nse))
    neiprob = ny.zeros((np0,nbi))
    timerange = math.floor(n/25)
    endtime = n - timerange

    for p0 in ny.arange(0,np0):
        for bi in ny.arange(0,nbi):
            for se in ny.arange(0,nse):
                check = 0
                numnei = 0
                numjump = 0
                minindx = 0
                tempEM = EM[:,p0,bi,se]
                tempend = len(Mstab[p0][bi][se])
#                tempengy = Mstab[p0][bi][se]
#                tempnei = Nei[p0][bi][se]
                if tempend > 0:
                    for lines in ny.arange(0,n-endtime):
                        if minindx < tempend:
                            if EM[lines,p0,bi,se] < Mstab[p0][bi][se][minindx] and EM[lines+timerange-1,p0,bi,se] >= Mstab[p0][bi][se][minindx]:

                                numjump+=1
                                minindx+=1
                                lfix = lines//10
                                rangefix = timerange//10
                                neich = ny.sum(Nei[lfix:lfix+rangefix,p0,bi,se])

                                if neich > 0:
                                    numnei+=1
                                lines+=timerange
                else:
                    numjump = 1
                    neich = ny.sum(Nei[:1,p0,bi,se])
                    if neich > 0:
                        numnei+=1
                if numjump > 0:
                    neiprobseed[p0,bi,se] = numnei/float(numjump)


            neiprob[p0,bi] = ny.mean(neiprobseed[p0,bi,:])


    return neiprobseed, neiprob


########################################################


#We get the bending extention times

def getallexttime(nli,np0,nse,nbi,namep0,namesee,namebi,name):

    x, y, z = np0, nbi, nse

#    extvals = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]

    oldreg1 = 0
    oldreg2 = 0
    extnum = 0

#    exttimes = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]

    numex = [None for _ in range(np0) ]

    extvals = [None for _ in range(np0)]

    names = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]


    for p0 in range (0,np0):
        oldlen = 0
        for bi in range(0, nbi):
            valslist = []
            numext = []
            for se in range(0, nse):

                filenme = namep0[p0] + namesee[se] + "Bi_" + namebi[bi] + name

                timeslist = []

                with open(filenme, 'r') as bfile:
                    blines = bfile.readlines()


                    sdata = blines[2].split()

                    start = float(sdata[1]) 

                    timeold = 0
                    

                    for i1 in range(2, nli+1):
                        
                        bdata = blines[i1].split();
                        time = int(bdata[0])

                        chk = 0

            
                        if i1 > 2:
                            num1 = float(bdata[1])
                            num2 = float(bdata[2])
        

                            if num1 > oldregI:

                                valslist.append(num1)

                                oldregI = num1

                                timeold = time
                

                        else:
                            oldregI = float(bdata[1]) 
                            valslist.append(oldregI)

                            
                        #benspeed.append(speedlist)
                    names[p0][bi][se] = filenme

        uniqvals = list(set(valslist))
        uniqvals.sort()
        extvals[p0] = uniqvals
        if len(uniqvals) > oldlen:
            oldlen = len(uniqvals)
        exlab = ny.linspace(1,oldlen-1,oldlen-1)
        numex[p0] = list(exlab)


    return extvals, numex, names

#################################################
#We get the bending extention times for up-down runs

def getupdownexttime(nli,np0,nse,nbi,namep0,namesee,namebi,name):

    x, y, z = np0, nbi, nse

#    extvals = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]

    oldreg1 = 0
    oldreg2 = 0
    extnum = 0

#    exttimes = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]

    numexup = [None for _ in range(np0) ]

    numexdown = [None for _ in range(np0) ]

    extvals1 = [None  for _ in range(np0)]

    extvals2 = [None for _ in range(np0)]

    names = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]


    for p0 in range (0,np0):
        oldlen1 = 0
        oldlen2 = 0
        for bi in range(0, nbi):
            upvalslist = []
            downvalslist = []
            numext = []
            for se in range(0, nse):

                filenme = namep0[p0] + namesee[se] + "Bi_" + namebi[bi] + name

                timeslist = []

                with open(filenme, 'r') as bfile:
                    blines = bfile.readlines()


                    sdata = blines[2].split()

                    start = float(sdata[1]) 

                    timeold = 0
                    

                    for i1 in range(2, nli+1):
                        
                        bdata = blines[i1].split();
                        time = int(bdata[0])

                        chk = 0

            
                        if i1 > 2:
                            num1 = float(bdata[3])
                            num2 = float(bdata[4])


                            if num1 > oldregI:

                                upvalslist.append(num1)

                                oldregI = num1

                                timeold = time

                            if num2 > oldregII:

                                downvalslist.append(num2)

                                oldregII = num2

                                timeold = time
                

                        else:
                            oldregI = float(bdata[3])
                            oldregII = float(bdata[4])
                            upvalslist.append(oldregI)
                            downvalslist.append(oldregII)

                            
                        #benspeed.append(speedlist)
                    names[p0][bi][se] = filenme

        uniqvals1 = list(set(upvalslist))
        uniqvals1.sort()
        uniqvals2 = list(set(downvalslist))
        uniqvals2.sort()
        extvals1[p0] = uniqvals1
        extvals2[p0] = uniqvals2
        if len(uniqvals1) > oldlen1:
            oldlen1 = len(uniqvals1)
        if len(uniqvals2) > oldlen2:
            oldlen2 = len(uniqvals2)

        exlab1 = ny.linspace(1,oldlen1-1,oldlen1-1)
        numexup[p0] = list(exlab1)
        exlab2 = ny.linspace(1,oldlen2-1,oldlen2-1)
        numexdown[p0] = list(exlab2)


    return extvals1, extvals2, numexup, numexdown,names



################################################
#We get the zippering speed for the CFC runs

def getallbendspeedCFC(nli,np0, nse, nbi, namep0, namesee,namebi,name,dy,bleny,Nx):
    
    x, y, z = np0, nbi, nse

    benspeed = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]

    oldreg1 = 0
    oldreg2 = 0
    extnum = 0

    exttimes = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]

    names = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]


    for p0 in range (0,np0):
        for bi in range(0, nbi):
            for se in range(0, nse):

                filenme = namep0[p0] + namesee[se] + "Bi_" + namebi[bi] + name

                extvals = bleny[p0]

                nextend = len(extvals)-1

                timeold = 0

                speedlist = [0 for _ in range(nextend)]


                timeslist = [0 for _ in range(nextend)]

                with open(filenme, 'r') as bfile:
                    blines = bfile.readlines()


                    sdata = blines[2].split()

                    start = float(sdata[1]) 

                    

                    for i1 in range(2, nli+1):
                        
                        bdata = blines[i1].split();
                        time = int(bdata[0])

                        chk = 0

            
                        if i1 > 2:
                            num1 = float(bdata[3])
                            num2 = float(bdata[4])
        

                            if num1 > oldregI:


                                for j in range(0,nextend):
                                    if num1 - extvals[j] == 0:

                                        move = (num1-oldregI)*dy/float(Nx)
                                        speed = move/(time-timeold)
            
                                        for k in range(0,j-oldindx):
                                            speedlist[k+oldindx] = speed
                                            timeslist[k+oldindx] = time

                                        oldindx = j

                                oldregI = num1
                                timeold = time
                

                        else:
                            oldregI = float(bdata[3]) 
                            oldregII = float(bdata[4])
                            oldindx = 0

                            
                        #benspeed.append(speedlist)
                    benspeed[p0][bi][se] = speedlist
                    exttimes[p0][bi][se] = timeslist
                    names[p0][bi][se] = filenme

    return benspeed, exttimes, names

#############################################
#We get the zippering speed for the CFC runs

def getupdownbendspeedCFC(nli,np0, nse, nbi, namep0, namesee,namebi,name,dy,uplen,downlen,Nx):
    
    x, y, z = np0, nbi, nse

    upbenspeed = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]


    downbenspeed = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]


    extnum = 0

    nreachup = ny.zeros((np0,nbi,nse))
    nreachdown = ny.zeros((np0,nbi,nse))

    upexttimes = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]

    downexttimes = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]

    names = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]


    for p0 in range (0,np0):
        for bi in range(0, nbi):

            upvals = uplen[p0]

            downvals = downlen[p0]

            nextendup = len(upvals)-1

            nextendown = len(downvals)-1

            for se in range(0, nse):

                filenme = namep0[p0] + namesee[se] + "Bi_" + namebi[bi] + name
                upspeedlist = [0 for _ in range (nextendup)] 
                uptimeslist = [0 for _ in range (nextendup)]
                downspeedlist = [0 for _ in range (nextendown)]
                downtimeslist = [0 for _ in range (nextendown)]

                with open(filenme, 'r') as bfile:
                    blines = bfile.readlines()

                    sdata = blines[2].split()

                    start = float(sdata[1]) 
 

                    timeold1 = 0

                    timeold2 = 0

                    for i1 in range(2, nli+1):
                        
                        bdata = blines[i1].split();
                        time = int(bdata[0])

                        chk = 0

            
                        if i1 > 2:
                            num1 = float(bdata[3])
                            num2 = float(bdata[4])

                            if num1 > oldregI:


                                for j in range(0,nextendup+1):

                                    if num1 - upvals[j] == 0:

                                        move = (num1-oldregI)*dy/float(Nx)
                                        speed = move/(time-timeold1)
            
                                        for k in range(0,j-oldindx1): 
                                            upspeedlist[k+oldindx1] = speed
                                            uptimeslist[k+oldindx1] = time
                                            nreachup[p0,bi,se]+=1

                                        oldindx1 = j

                                oldregI = num1
                                timeold1 = time

                            if num2 > oldregII:


                                for j in range(0,nextendown+1):
                                    if num2 - downvals[j] == 0:

                                        move = (num2-oldregII)*dy/float(Nx)
                                        speed = move/(time-timeold2)
            
                                        for k in range(0,j-oldindx2):
                                            downspeedlist[k+oldindx2] = speed
                                            downtimeslist[k+oldindx2] = time
                                            nreachdown[p0,bi,se]+=1

                                        oldindx2 = j

                                oldregII = num2
                                timeold2 = time
                

                        else:
                            oldregI = float(bdata[3]) 
                            oldregII = float(bdata[4])
                            oldindx1 = 0
                            oldindx2 = 0

                            
                        #benspeed.append(speedlist)
                upbenspeed[p0][bi][se] = upspeedlist
                downbenspeed[p0][bi][se] = downspeedlist
                upexttimes[p0][bi][se] = uptimeslist
                downexttimes[p0][bi][se] = downtimeslist
                names[p0][bi][se] = filenme

    return upbenspeed, downbenspeed, upexttimes, downexttimes, nreachup, nreachdown, names


#################################################
def collectbendspeedsCFC(np0,nbi,nse,speed):

    colspeed = [[None for _ in range (nbi)] for _ in range (np0) ]


    for p0 in range (0,np0):
        for bi in range(0, nbi):
            temp = []
            for se in range(0, nse):                
                if len(speed[p0][bi][se]) > 0:
                    for ele in speed[p0][bi][se]:
                        if ele > 0:
                            temp.append(ele)
            colspeed[p0][bi] = temp


    return colspeed


#################################################
def combinebendspeedsCFC(np0,nbi,nse,speed1,speed2):

    colspeed = [[None for _ in range (nbi)] for _ in range (np0) ]

    spemax = 0


    for p0 in range (0,np0):
        for bi in range(0, nbi):
            temp1 = []
#            for se in range(0, nse):
            if speed1[p0][bi]:
                for ele in speed1[p0][bi]:
                    temp1.append(ele)
            if speed2[p0][bi]:
                for ele in speed2[p0][bi]:
                    temp1.append(ele)
            colspeed[p0][bi] = temp1
            tempmax = max(colspeed[p0][bi])
            if tempmax > spemax:
                spemax = tempmax


    return colspeed, spemax


##################################################

#Here we get the frequency of extension speeds
def getspeedfreq(lowbound,highbound,numbins,np0,nbi,nse,bleny,bspe):


    delx = (highbound-lowbound)/float(numbins)
    Bins = ny.linspace(lowbound,highbound,numbins)
    Freq = ny.zeros((np0,nbi,nse,numbins))
    Freqseed = ny.zeros((np0,nbi,numbins))

    Binlabel = [None for _ in range(np0)]


    for p0 in range(0,np0):
        for bi in range(0,nbi):
            maxindx = 0
#            numext = len(bleny[p0])-1
            for se in range(0,nse):
                maxindxseed = 0
                X = ny.array(bspe[p0][bi][se])
                Freq[p0,bi,se,:] = setupbins(X,lowbound,highbound,numbins)
                Freqseed[p0,bi,:] += Freq[p0,bi,se,:]

        Binlabel[p0] = Bins


    return Bins, Binlabel, Freq, Freqseed

################################################
#Here we get the frequency of extension speeds
def getspeedfreq2(lowbound,highbound,delx,np0,nbi,nse,bleny,bspe):


    numbins = math.ceil((highbound-lowbound)/float(delx))
    Bins = ny.linspace(lowbound,highbound,numbins)
    Freq = ny.zeros((np0,nbi,nse,numbins))
    Freqseed = ny.zeros((np0,nbi,numbins))

    Binlabel = [None for _ in range(np0)]


    for p0 in range(0,np0):
        for bi in range(0,nbi):
            maxindx = 0
#            numext = len(bleny[p0])-1
            for se in range(0,nse):
                maxindxseed = 0
                X = ny.array(bspe[p0][bi][se])
                Freq[p0,bi,se,:] = setupbins(X,lowbound,highbound,numbins)
                Freqseed[p0,bi,:] += Freq[p0,bi,se,:]

        Binlabel[p0] = Bins


    return Bins, Binlabel, Freq, Freqseed



################################
def getupdownbendextenprob(np0,nbi,nse,nex,nactup, nactdown,upspe,downspe):

#    probseed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]

#    condprobseed = [[[None for _ in range(nse)] for _ in range(nbi) ] for _ in range(np0)]



    probtot = [[None for _ in range(nbi) ] for _ in range(np0)]
    probnext = [[None for _ in range(nbi) ] for _ in range(np0)]
    probupaccelseed = ny.zeros((np0,nbi,nse))
    probdownaccelseed = ny.zeros((np0,nbi,nse))
    probnextaccelseed = ny.zeros((np0,nbi,nse))
    probupaccel = ny.zeros((np0,nbi))
    probdownaccel = ny.zeros((np0,nbi))
    probnextaccel = ny.zeros((np0,nbi))


    for p0 in range(0,np0):
        nposs = len(nex[p0])

        for bi in range(0,nbi):
            numr = ny.zeros(nposs)
            extpro = ny.zeros(nposs)
            nextpro = ny.zeros(nposs)
            
            for se in range(0,nse):
                numupspeedup = 0
                numdownspeedup = 0
                upmade = int(nactup[p0,bi,se])
                downmade = int(nactdown[p0,bi,se])
                for extup in range(0,upmade):
                    numr[extup]+=1
                    if extup > 0:
                        if upspe[p0][bi][se][extup] >upspe[p0][bi][se][extup-1]:
                            numupspeedup+=1
                for extdown in range(0,downmade):
                    numr[extdown]+=1
                    if extdown > 0:
                        if downspe[p0][bi][se][extdown] > downspe[p0][bi][se][extdown-1]:
                            numdownspeedup+=1
                if upmade > 1:
                    probupspeedup = numupspeedup/float(upmade-1)
                    probupaccelseed[p0,bi,se] = probupspeedup
                else:
                    probupaccelseed[p0,bi,se] = ny.nan

                if downmade > 1:
                    probdownspeedup = numdownspeedup/float(downmade-1)
                    probdownaccelseed[p0,bi,se] = probdownspeedup
                else:
                    probdownaccelseed[p0,bi,se] = ny.nan

                probnextaccelseed[p0,bi,se] = (probupspeedup + probdownspeedup)/float(2)

            



            for ext in range(0,nposs):
                extpro[ext] = numr[ext]/float(2*nse)
                if ext > 0:
                    nextpro[ext] = extpro[ext]/float(extpro[ext-1])
                else:
                    nextpro[ext] = extpro[ext]


            probtot[p0][bi] = extpro
            probnext[p0][bi] = nextpro
            probupaccel[p0,bi] = ny.nanmean(probupaccelseed[p0,bi,:])
            probdownaccel[p0,bi] = ny.nanmean(probdownaccelseed[p0,bi,:])
            probnextaccel[p0,bi] = ny.nanmean(probnextaccelseed[p0,bi,:])



    return probtot, probnext, probupaccelseed,probdownaccelseed

################################################
#We average over seeds

def avgallbendspeedCFC(np0, nse, nbi, bleny, bspe):



    x, y = np0, nbi

    avgspeed = [[None for _ in range(y) ] for _ in range(x)]

    avgspeedstd = [[None for _ in range(y) ] for _ in range(x)]

#    Retained = []

    for p0 in range (0,np0):
        for bi in range(0, nbi):

            numext = len(bleny[p0][bi])-1

            speed = [0]*numext
#            retain = []



            speedstd = [0]*numext

            keep = [0]*numext

            for se in range(0, nse):


                if len(bspe[p0][bi][se]) > 0:

                    for ext in range(0, len(bspe[p0][bi][se])):

                        speed[ext] = speed[ext] + bspe[p0][bi][se][ext]

#                       retain.append(bspe[indx+se])
                        keep[ext] = keep[ext] + 1

            for j in range(0, numext):

                elekeep = []

                for sed in range(0,nse):
                    if len(bspe[p0][bi][sed]) > j:
                        if bspe[p0][bi][sed][j] is not None:
                            elekeep.append(bspe[p0][bi][sed][j])
                if len(elekeep) > 0:
                    speedstd[j] = statistics.pstdev(elekeep)


                if keep[j] >0:

                    speed[j] = speed[j]/float(keep[j])




            avgspeed[p0][bi] = speed

            avgspeedstd[p0][bi] = speedstd


    return avgspeed, avgspeedstd

#################################################
#We average over the nonzero speeds in each seed

def avgallbendspeedalt(np0, nse, nbi, bleny, bspe):

    x, y = np0, nbi

    avgspeed = [[0 for _ in range(y) ] for _ in range(x)]

    avgspeedstd = [[None for _ in range(y) ] for _ in range(x)]


#    Retained = []

    for p0 in range (0,np0):
        for bi in range(0, nbi):

            numext = len(bleny[p0])-1

            speed = ny.zeros(numext)
#            retain = []



            speedstd = [0]*numext

            keep = ny.zeros(numext)


            for se in range(0, nse):
                for ext in range(0, numext):
                    if bspe[p0][bi][se][ext] != 0:

                        speed[ext] = speed[ext] + bspe[p0][bi][se][ext]

#                       retain.append(bspe[indx+se])
                            
                        keep[ext] = keep[ext] + 1


            for ext in range(0, numext):


                if keep[ext] >0:

                    speed[ext] = speed[ext]/float(keep[ext])


                elekeep = []

                for se in range(0,nse):

                    if bspe[p0][bi][se][ext] != 0:

                        elekeep.append(bspe[p0][bi][se][ext])


                if len(elekeep) > 0:

                    speedstd[ext] = statistics.pstdev(elekeep)


            avgspeed[p0][bi] = list(speed)

            avgspeedstd[p0][bi] = speedstd


    return avgspeed, avgspeedstd






####################################################
#We get the zippering speed for the up-down runs

def getbendspeedupdown(fle,ele,nli,dy,Nx):

    bentimesup = []
    bentimesdown = []


    with open(ele, 'r') as bfile:
        blines = bfile.readlines()

        bdata = blines[1].split();
        ex1 = int(bdata[0])
        ex2 = int(bdata[1])

        if ex1>0:
            l1data = blines[2].split();
            for i1 in range(0,ex1):
                bentimesup.append(l1data[i1+1])

        if ex2>0:
            l2data = blines[3].split();
            for i2 in range(0,ex2):
                bentimesdown.append(l2data[i2+1])

    benspeedup = []
    benspeeddown = []
    oldreg = 0

    timeold = 0

    with open(fle, 'r') as bfile:
        blines = bfile.readlines()

        for i1 in range(2, nli+1):
            bdata = blines[i1].split();
            time = int(bdata[0])

            
            if i1 > 2:
                num1 = float(bdata[1])


                if num1 > oldreg:
                    Leny = (num1-oldreg)*dy/Nx
                    speed = Leny/float(time-timeold)
                    benspeed.append(speed)
                    oldregI = num1
                    oldregII = num2
                    timeold = time

                if num1 > oldregI + oldregII and num2 != 0:  #multiple extensions
                    Leny = (num1-oldregI)*dy/Nx
                    speed = Leny/float(time-timeold)
                    benspeed.append(speed)
                    benspeed.append(speed)
                    oldregI = num1
                    oldregII = num2
                    timeold = time

                if num1 > oldregI + oldregII and num2 == 0:

                    Leny = (num1-oldregI)*dy/Nx
                    speed = Leny/float(time-timeold)
                    benspeed.append(speed)
                    benspeed.append(speed)
                    benspeed.append(speed)
                    oldregI = num1
                    oldregII = num2
                    timeold = time


            else:
                oldregI = float(bdata[1]) 



    return benspeed




#####################################################
#Here we set the bending times

def storebendspeeds(bsp):

    b1 = 0
    b2 = 0
    b3 = 0

    if len(bsp) > 0:
        b1 = bsp[0]

    if len(bsp) > 1:
        b2 = bsp[1]
    
    if len(bsp) > 2:
        b3 = bsp[2]


    return b1, b2, b3


###################################
#This will grab the nonzero values

def grabnonzeroele(ary, l):

    keep = []    

    for i in range(0,l):
        if ary[i] != 0:
            keep.append(ary[i])


    return keep



###############################################

#This will average over the non-zero values
def avgovernonzero(ary, l):

    Avg = 0
    nele = 0

    for i in range(0,l):
        if ary[i] != 0:
            Avg = Avg + ary[i]
            nele = nele+1
            

    if nele != 0:
        Avg = Avg/float(nele)

    return Avg


#############################################
#We get the bending times for the up-down runs 

def getbendtimeupdown(fle):
    
    bentimesup = []
    bentimesdown = []


    with open(fle, 'r') as bfile:
        blines = bfile.readlines()

        bdata = blines[1].split();
        ex1 = int(bdata[0])
        ex2 = int(bdata[1])

        if ex1>0:
            l1data = blines[2].split();
            for i1 in range(0,ex1):
                bentimesup.append(l1data[i1+1])

        if ex2>0:
            l2data = blines[3].split();
            for i2 in range(0,ex2):
                bentimesdown.append(l2data[i2+1])



    return bentimesup, bentimesdown


#########################################

#Here we get the averages over files

def getavgs(nli, nfi, ts, eg):

    tavg = ny.zeros((nli,nfi))
    PBavg = ny.zeros((nli,nfi))
    PBstd = ny.zeros((nli,nfi))

    PBratfit = ny.zeros((nli,nfi))

    for i in range(0,nli):
        for j in range(0,nfi):
            tavg[i][j] = ny.mean(ts[i,j,:])
            PBavg[i][j] = ny.mean(eg[i,j,:])
            PBstd[i][j] = ny.std(eg[i,j,:])

    return tavg, PBavg, PBstd


###########################################
#Here we'll get the average over the last couple of timesteps


def edgeendavg(nfi,ns, eg, nli, nind):

    PBend = ny.zeros((ns,nfi))
    PBestd = ny.zeros(nfi)

    for i in range(0,nfi):
        for j in range(0,ns):
            eee = eg[nind:nli,i,j]
            PBend[j][i] = ny.mean(eee)

        PBestd[i] = ny.std(PBend[:,i])


    return PBend, PBestd



#################################################
#Here we get the rise time

def risetimeavg(PB,tstart,nline):

#First we find the average of the last couple timesteps

    enavg = ny.mean(PB[tstart:nline])

#Then the times we reach 10% and 90% of that value

    t1chk = 0
    t2chk = 0
    t1 = 0
    t2 = nline
    for times in range(0,nline):
        if PB[times] >= 0.1*enavg and t1chk == 0:
            t1 = times
            t1chk = 1
        if PB[times] >= 0.9*enavg and t2chk == 0:
            t2 = times
            t2chk = 1


    trise = t2-t1


    return trise


##########################################

def getallrisetime(EM,tstart,nline,np0,nbi,nse):


    trise = ny.zeros((np0,nbi,nse))


    triseavg = ny.zeros((np0,nbi))

    trisestd = ny.zeros((np0,nbi))

    for p0 in range(0,np0):
        for bi in range(0,nbi):
            for se in range(0,nse):

#First we find the average of the last couple timesteps

                enavg = ny.mean(EM[tstart:,p0,bi,se])

#Then the times we reach 10% and 90% of that value

                t1chk = 0
                t2chk = 0
                t1 = 0
                t2 = nline
                for times in range(0,nline):
                    if EM[times,p0,bi,se] >= 0.1*enavg and t1chk == 0:
                        t1 = times
                        t1chk = 1
                    if EM[times,p0,bi,se] >= 0.9*enavg and t2chk == 0:
                        t2 = times
                        t2chk = 1


                trise[p0,bi,se] = t2-t1

            triseavg[p0,bi] = math.ceil(ny.mean(trise[p0,bi,:]))

            trisestd[p0,bi] = ny.std(trise[p0,bi,:])


    return trise,triseavg, trisestd

#########################################
#We get the number of neighbor changes at rise time
def getallneiatandafterrise(Nei,Tri,np0,nbi,nse):


    Natrise = ny.zeros((np0,nbi,nse))
    Nafrise = ny.zeros((np0,nbi,nse))


    Natriseavg = ny.zeros((np0,nbi))
    Natrisestd = ny.zeros((np0,nbi))
    
    Nafriseavg = ny.zeros((np0,nbi))
    Nafrisestd = ny.zeros((np0,nbi))

    for p0 in range(0,np0):
        for bi in range(0,nbi):
            for se in range(0,nse):

                rtime = int(math.ceil(Tri.edgerise[p0,bi,se]/10))
                Natrise[p0,bi,se] = ny.sum(Nei.neichg12[:rtime,p0,bi,se])
                Nafrise[p0,bi,se] = ny.sum(Nei.neichg12[rtime,p0,bi,se])



            Natriseavg[p0,bi] = ny.mean(Natrise[p0,bi,:])

            Natrisestd[p0,bi] = ny.std(Natrise[p0,bi,:])

            Nafriseavg[p0,bi] = ny.mean(Nafrise[p0,bi,:])

            Nafrisestd[p0,bi] = ny.std(Nafrise[p0,bi,:])

    return Natrise,Natriseavg, Natrisestd, Nafrise, Nafriseavg, Nafrisestd

############################################


#We get the shape index and actual perimeter

def getshapeperi(fle, nline, ncel):

    Perim = []
    TS = []

    with open(fle, 'r') as dfile:
        dlines = dfile.readlines()

        for i1 in range(1, nline+3):
            ddata = dlines[i1].split();

            if i1 == 1:
                A1 = float(ddata[2])
                A2 = float(ddata[3])

                P1 = float(ddata[5])
                P2 = float(ddata[6])

#            if i1 == n+3:
            if i1 > 2: 
                
                TS.append(float(ddata[0]))
                Perim.append(float(ddata[1]))
                Perim[i1-3] = Perim[i1-3]/float(ncel)


    return A1, A2, P1, P2, TS,Perim


########################################################
#We get the shape index and actual perimeter at the end

def getshapeperiend(fle, nline, ncel, tail):

    Perim = []
    TS = []

    with open(fle, 'r') as dfile:
        dlines = dfile.readlines()

        for i1 in range(1, nline+3):
            ddata = dlines[i1].split();

            if i1 == 1:
                A1 = float(ddata[2])
                A2 = float(ddata[3])

                P1 = float(ddata[5])
                P2 = float(ddata[6])

#            if i1 == n+3:
            if i1 > 2: 
                
                TS.append(float(ddata[0]))
                Perim.append(float(ddata[1]))
                Perim[i1-3] = Perim[i1-3]/float(ncel)
        Peri = sum(Perim[tail:nline])/len(Perim[tail:nline])


    return A1, A2, P1, P2, TS,Peri


#######################################

def Avgoverseeds(EM, lenp0, lenbi, nlines):
    
    Eavg = ny.zeros((nlines,lenp0,lenbi))
    Estd = ny.zeros((nlines,lenp0,lenbi))
    for lines in range(0,nlines):
        for p0 in range(0,lenp0):
            for bi in range(0,lenbi):
                Eavg[lines,p0,bi] = ny.mean(EM[lines,p0,:,bi])
                Estd[lines,p0,bi] = ny.std(EM[lines,p0,:,bi])
    


    return Eavg, Estd

#############################################
        
    
def Avgoverendtime(EM, lenp0, lenbi, nseeds, nlstart, nlend):


    Eavg = ny.zeros((lenp0,lenbi))
    Estd = ny.zeros((lenp0,lenbi))
    Eseeds = ny.zeros((lenp0,lenbi,nseeds))

    for p0 in range(0,lenp0):
        for bi in range(0,lenbi):
            for seeds in range(0,nseeds):
                Eseeds[p0,bi,seeds] = ny.mean(EM[nlstart:nlend,p0,bi,seeds])

            Eavg[p0,bi] = ny.mean(Eseeds[p0,bi,:])
            Estd[p0,bi] = ny.std(Eseeds[p0,bi,:])



    return Eavg, Estd


###############################################
def freqoverendtime(X,lowbound,highbound,numbins,np0,nbi,nse,nlstart,nlend):


    delx = (highbound-lowbound)/float(numbins)
    Bins = ny.linspace(lowbound,highbound,numbins)
    Freq = ny.zeros((np0,nbi,numbins))

    Freqlist = [[0 for _ in range (nbi) ] for _ in range(np0) ]

    for p0 in range(0,np0):
        for bi in range(0,nbi):
            freqtemp = []
            for se in range(0,nse):
                Mf = ny.mean(X[nlstart:nlend,p0,bi,se])

                freqtemp.append(Mf)

            Freqlist[p0][bi] = freqtemp


    return Freqlist, Bins
    



########################################

def findhalfneichanges(Nei,np0,nbi,nse,nli):


    for p0 in ny.arange(0,np0):
        for bi in ny.arange(0,nbi):
            check = 0
            for se in ny.arange(0,nse):
                for lines in ny.arange(0,nli):
                    if Nei.neichg1[lines,p0,bi,se] == 2 and check == 0:
                        print("Half neighbor change in ", Nei.name[p0][bi][se], " at ts = " , lines )
                        check = 1

########################################

def convertlisttoarry(X, np0,nbi,nse,nli):

    Nary = ny.zeros((nli,np0,nbi,nse))

    for p0 in range(0,np0):
        for bi in range(0,nbi):
            for se in range(0,nse):
                Nary[:,p0,bi,se] = ny.array(X[p0][bi][se])



    return Nary

########################################
def plotedgeslong(EM,TM,P,B,Nse,ytick,yrange,selabel,floop):

    plt.figure(floop)
    cyticks1=ytick
    cyrange=yrange

    cmap = plt.get_cmap('tab20',Nse)


    for seed in range(0,Nse):

        tsavgskp = TM

        PBseedskp = EM[P,B,seed,:]
        plt.plot(tsavgskp,PBseedskp, color = cmap(seed), label = selabel[seed] )
        plt.grid()
        plt.tick_params(axis='both', labelsize = 28)
        plt.ylim(cyrange)
        plt.legend()
#        ax1.plot(tsavgskp,PBseedskp, color = cmap(seed), label = selabel[seed] )
#        ax1.grid()
#        ax1.tick_params(axis='both', labelsize = 28)
#        ax1.set_ylim(cyrange)

#        ax2.plot(Nei.time[:,P,B,seed],Nei.neichg12[:,P,B,seed], color = cmap(seed), label = selabel[seed])
#        ax2.grid()
#        ax2.tick_params(axis='both', labelsize = 28)

    floop+=1
    return floop





#########################################

def plotedges(EM,skip,P,B,Nse,ytick,yrange,selabel,floop):

    plt.figure(floop)
    cyticks1=ytick
    cyrange=yrange

    cmap = plt.get_cmap('tab20',Nse)


    for seed in range(0,Nse):

        tsavgskp = EM.tsdat[::skip,P,B,seed]

        PBseedskp = EM.edgedat[::skip,P,B,seed]
        plt.plot(tsavgskp,PBseedskp, color = cmap(seed), label = selabel[seed] )
        plt.grid()
        plt.tick_params(axis='both', labelsize = 28)
        plt.ylim(cyrange)

#        ax1.plot(tsavgskp,PBseedskp, color = cmap(seed), label = selabel[seed] )
#        ax1.grid()
#        ax1.tick_params(axis='both', labelsize = 28)
#        ax1.set_ylim(cyrange)

#        ax2.plot(Nei.time[:,P,B,seed],Nei.neichg12[:,P,B,seed], color = cmap(seed), label = selabel[seed])
#        ax2.grid()
#        ax2.tick_params(axis='both', labelsize = 28)

    floop+=1
    return floop

######################################################
def plotsingleedge(EM,skip,P,B,S,ytick,yrange,selabel,floop):

    plt.figure(floop)
    cyticks1=ytick
    cyrange=yrange

#    cmap = plt.get_cmap('tab20',Nse)


    tsavgskp = EM.tsdat[::skip,P,B,S]

    PBseedskp = EM.edgedat[::skip,P,B,S]
#    plt.plot(tsavgskp,PBseedskp, color = cmap(seed), label = selabel[seed] )
    plt.plot(tsavgskp,PBseedskp, label = selabel[S] )
    plt.grid()
    plt.tick_params(axis='both', labelsize = 28)
    plt.ylim(cyrange)


    floop+=1
    return floop


########################################


def plotneighchange(Nei,P,B,S,Nse,ytick,yrange,ytickneigh,selabel,floop):

    cyticks1=ytick
    cyrange=yrange
    cyticksens=ytickneigh

    cmap = plt.get_cmap('tab20',Nse)

#    for seed in range(0,Nse):

    plt.plot(Nei.time[:,P,B,S],Nei.neichg12[:,P,B,S], color = cmap(S), label = selabel[S])
    plt.grid()
    plt.tick_params(axis='both', labelsize = 28)

    floop+=1
    return floop


##########################################

def plotedgeandneighchange(EM,Nei,skip,P,B,Nse,ytick,yrange,ytickneigh,selabel,floop):

    fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )
    cyticks1=ytick
    cyrange=yrange
    cyticksens=ytickneigh

    cmap = plt.get_cmap('tab20',Nse)


    for seed in range(0,Nse):

        tsavgskp = EM.tsdat[::skip,P,B,seed]

        PBseedskp = EM.edgedat[::skip,P,B,seed]
        ax1.plot(tsavgskp,PBseedskp, color = cmap(seed), label = selabel[seed] )
        ax1.grid()
        ax1.tick_params(axis='both', labelsize = 28)
        ax1.set_ylim(cyrange)

        ax2.plot(Nei.time[:,P,B,seed],Nei.neichg12[:,P,B,seed], color = cmap(seed), label = selabel[seed])
        ax2.grid()
        ax2.tick_params(axis='both', labelsize = 28)

    floop+=1
    return floop

####################################################

def plotsingleedgeandneighchange(EM,Nei,skip,P,B,S,ytick,yrange,ytickneigh,selabel,floop):

    fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )
    cyticks1=ytick
    cyrange=yrange
    cyticksens=ytickneigh

#    cmap = plt.get_cmap('tab20',Nse)


    tsavgskp = EM.tsdat[::skip,P,B,S]

    PBseedskp = EM.edgedat[::skip,P,B,S]
#    ax1.plot(tsavgskp,PBseedskp, color = cmap(seed), label = selabel[seed] )
    ax1.plot(tsavgskp,PBseedskp, label = selabel[S] )
    ax1.grid()
    ax1.tick_params(axis='both', labelsize = 28)
    ax1.set_ylim(cyrange)

#    ax2.plot(Nei.time[:,P,B,seed],Nei.neichg12[:,P,B,seed], color = cmap(seed), label = selabel[seed])
    ax2.plot(Nei.time[:,P,B,S],Nei.neichg12[:,P,B,S],label = selabel[S])
    ax2.grid()
    ax2.tick_params(axis='both', labelsize = 28)

    floop+=1
    return floop





########################################################

def plotedgediffandnei(EM,Nei,skip,P,B,Nse,ytick,yrange,ytickneigh,selabel,floop):

    cmap = plt.get_cmap('tab20',Nse)
    fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )
    cyrange=yrange
    cyticksens=ytick


    for seed in range(0,Nse):

        PBseedskp = EM.edgediff[::skip,P,B,seed]
        ax1.plot(EM.tsdat[::skip,P,B,seed],PBseedskp, color = cmap(seed), label = selabel[seed] )
        ax1.grid()
        ax1.tick_params(axis='both', labelsize = 28)

        ax2.plot(Nei.time[:,P,B,seed],Nei.neichg12[:,P,B,seed], color = cmap(seed), label = selabel[seed])
        ax2.grid()
        ax2.tick_params(axis='both', labelsize = 28)


    floop +=1
    return floop


############################################################
#Now we plot the ensemble avgs

def plotensmblavg(EM,Nei,skip,P,Bran,ytick,yrange,ytickneigh,bilabel,floop):

    cmap = plt.get_cmap('rainbow',Bran)
    fig, (ax1, ax2 ) = plt.subplots(2,1, figsize=(6,6) )
    cyrange=yrange
    cyticksens=ytick

    bihalf = Bran//2
    for bi in range(0,bihalf):
        bi2 = 2*bi

        PBseedskp = EM.edge[::skip,P,bi2]
        ax1.plot(EM.tsdat[::skip,P,bi2],PBseedskp, color = cmap(bi2), label = bilabel[bi2] )
        ax1.grid()
        ax1.tick_params(axis='both', labelsize = 28)
        ax1.set_ylim(cyrange)

        ax2.plot(Nei.time[:,P,bi2,0],Nei.neichgsum[:,P,bi2], color = cmap(bi2), label = bilabel[bi2])
        ax2.grid()
        ax2.tick_params(axis='both', labelsize = 28)
        ax2.set_yticks(cyticksens)

    floop +=1
    return floop

##################################################
#We plot contour maps

def plotregimemap(Reg,Xv,Yv,bound,ytick,mapindx,floop):

    cyticks1=ytick
    if mapindx == 1:
        cmap2 = 'rainbow'
    elif mapindx == 2:
        cmap2 = 'jet'
    elif mapindx == 3:
        cmap2 = 'plasma'
    else:
        cmap2 = 'terrain'


    [X,Y] = ny.meshgrid(Xv,Yv)

    regimefix = ny.zeros((len(Yv),len(Xv)))
    for x in range(0,len(Xv)):
        for y in range(0,len(Yv)):
            regimefix[y,x] = Reg[x,y]

    plt.figure(floop)
    contour1 = plt.figure(floop)
    contour = plt.contourf(X,Y,regimefix,levels=bound, cmap=cmap2)
    cbar = plt.colorbar(contour, location='right')
    cbar.ax.tick_params(labelsize=32)
    plt.yticks(ytick)
    plt.tick_params(axis='both',labelsize=26)
    plt.grid()

    floop+=1
    return floop

###################################################
#We plot half the lines

def plothalfthelines(Lin,xvals,Bran,bilabel,floop):

    cmap = plt.get_cmap('rainbow',Bran)
    plt.figure(floop)
    bihalf = Bran//2
    for bi in range(0,bihalf):
        bi2 = 2*bi
        plt.plot(xvals,Lin[:,bi2],color=cmap(bi2), label=bilabel[bi2])
        plt.grid()
        plt.tick_params(axis='both', labelsize = 28)


    plt.legend(loc='lower left', bbox_to_anchor=(0.99,-0.02),fontsize=20)

    floop +=1
    return floop

###############################################################

def plotbargraph(xpar,labe, widt, color, Labe, figl):

    bwidth = widt
    multiply = 0
    bcolors = color
    plt.figure(figl)
    fig = plt.subplot(figsize=(6,6))

    plt.plot(labe, xpar)
    plt.bar(labe, xpar, capsize=3, width = bwidth, label= Labe)




##################################################################

def plotbargraphallp0(xpar,labeup,widt, color, Labe, pove, bove, titl, figl, chk1,chk2):

    bwidth = widt
    multiply = 0
    bcolors = color


    for p0 in range(0,pove):
        labelsup = ny.arange(len(labeup[p0]))
        plt.figure(figl)
        plt.title(titl)
#        fig = plt.subplot(figsize=(6,6))

        for b0 in range(0,bove):

            labelsuse = [x + (b0+1)*bwidth for x in labelsup ]

            if chk1 == 0:
                plt.plot(labelsuse, xpar[p0][b0])
            plt.bar(labelsuse, xpar[p0][b0], capsize=3, width = bwidth, label= Labe[b0])


            plt.grid()
            plt.tick_params(axis='both', labelsize=28)
            if chk1 == 1:
                plt.xticks(ticks=labelsuse)
            if chk2 == 0:
                plt.xticks([])

        figl+=1
        plt.legend(loc='lower left', bbox_to_anchor=(0.99,0), fontsize=22)


    figl +=1

    return figl

############################################
def plothistgraphallp0(xpar,labeup,widt, color, Labe, pove, bove, titl, figl, chk1,chk2):

    bwidth = widt
    multiply = 0
    bcolors = color

    for p0 in range(0,pove):
#        labelsup = ny.arange(len(labeup[p0]))
        labelsup = ny.arange(len(labeup))
#        nbins = len(labeup[p0])
        nbins = len(labeup)
        plt.figure(figl)
        plt.title(titl)
#        fig = plt.subplot(figsize=(6,6))

        dats = [ ] 


        for b0 in range(0,bove):

            if len(xpar[p0][b0]) > 0:
                dats.append(xpar[p0][b0])

        labelsuse = [x + (b0+1)*bwidth for x in labelsup ]

        if chk1 == 0:
            plt.plot(labelsuse, xpar[p0][b0])
        #plt.bar(labelsuse, xpar[p0][b0], capsize=3, width = bwidth, label= Labe[b0])
        plt.hist(dats, bins=labeup, label = Labe,stacked=False )

        plt.grid()
#        plt.xticks(labeup[:-1],xpar[p0][b0])
        plt.tick_params(axis='both', labelsize=28)
        plt.tight_layout()

#            if chk1 == 1:
#                plt.xticks(ticks=labelsuse)
#            if chk2 == 0:
#                plt.xticks([])

        figl+=1
        plt.legend(loc='upper left', bbox_to_anchor=(0.88,1.01), fontsize=22)


    figl +=1

    return figl






#############################################
def plotbargraphupdown(xpar1,xpar2,labeup,labedown,widt, color, Labe, pove, bove, titl, figl):

    bwidth = widt
    multiply = 0
    bcolors = color


    for p0 in range(0,pove):
        labelsup = ny.arange(len(labeup[p0]))
        labelsdown = ny.arange(len(labedown[p0]))
        plt.figure(figl)
        plt.title(titl)
#        fig = plt.subplot(figsize=(6,6))

        for b0 in range(0,bove):

            labelsupuse = [x + (b0+3)*bwidth for x in labelsup ]
            labelsdownuse = [x + (b0+4)*bwidth for x in labelsdown]


            plt.plot(labelsupuse, xpar1[p0][b0])
            plt.bar(labelsupuse, xpar1[p0][b0], capsize=3, width = bwidth, label= Labe[b0])

            plt.plot(labelsdownuse, xpar2[p0][b0])
            plt.bar(labelsdownuse, xpar2[p0][b0], capsize=3, width = bwidth, label= Labe[b0])




            plt.grid()
            plt.tick_params(axis='both', labelsize=28)
            plt.xticks([])

            figl+=1
        plt.legend(loc='lower left', bbox_to_anchor=(0.99,0), fontsize=22)


    figl +=1

    return figl

##############################################

#We make a scatter plot

def plotscatter(yvals1,yvals2,yerr1,yerr2,B,bilabel,plabel,floop):

    cmap = plt.get_cmap('rainbow',B+1)
    plt.figure(floop)

#    for bi in range(0,Bran):
    Y1 = ny.zeros(len(plabel)) 
    Y2 = ny.zeros(len(plabel))
    Ye1 = ny.zeros(len(plabel))
    Ye2 = ny.zeros(len(plabel))
    for p0 in range(0,len(plabel)):
        Y1[p0] = yvals1[p0][B]
        Y2[p0] = yvals2[p0][B]
        Ye1[p0] = yerr1[p0][B]
        Ye2[p0] = yerr2[p0][B]
#        plt.scatter(plabel,Y1, yerr=Ye1, fmt='o', capsize=2, color=cmap(bi), label=bilabel[bi])
#        plt.scatter(plabel,Y2, yerr=Ye2, fmt='o', capsize=2,  marker = "s", color=cmap(bi), label=bilabel[bi])
        plt.errorbar(plabel,Y1, yerr=Ye1, fmt='o', capsize=2, color=cmap(B), label=bilabel[B])
        plt.errorbar(plabel,Y2, yerr=Ye2, fmt='o', capsize=2,  marker = "s", color=cmap(B+1), label=bilabel[B])

        plt.grid()
        plt.tick_params(axis='both', labelsize = 28)


    plt.legend(loc='lower left', bbox_to_anchor=(0.99,-0.02),fontsize=20)

    floop +=1
    return floop





