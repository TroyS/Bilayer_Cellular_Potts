#This will contain functions for data analysis

import numpy as ny
import statistics
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import plotly.graph_objects as go
from matplotlib.colors import TABLEAU_COLORS, same_color
from matplotlib.pyplot import cm
from matplotlib import cm


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

    for p0 in ny.arange(0,np0):
        for bi in ny.arange(0, nbi):
            for se in ny.arange(0, nse):
                filenme = namep0[p0] + namesee[se] + "Bi_" + namebi[bi] + name

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


    return ts, PB, Dif


############################################
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


    for p0 in ny.arange(0,np0):
        for bi in ny.arange(0,nbi):
            for se in ny.arange(0,nse):

                filenme = namep0[p0] + namesee[se] + "Bi_" + namebi[bi] + name


                with open(filenme, 'r') as nfile:
                    nlines = nfile.readlines()

                    for i1 in ny.arange(1, n+1):
                        ndata = nlines[i1].split();
                        Time[i1-1,p0,bi,se] = float(ndata[0])
                        neichg1[i1-1,p0,bi,se] = float(ndata[1])
                        neichg2[i1-1,p0,bi,se] = float(ndata[2])
                        neichg12[i1-1,p0,bi,se]=float(ndata[1])+float(ndata[2])
                        neichgsum[i1-1,p0,bi]+=(float(ndata[1])+float(ndata[2]))

                neichgtotseed[p0,bi,se] = ny.sum(neichg12[:,p0,bi,se])

            neichgtot[p0,bi] = ny.sum(neichgtotseed[p0,bi,:])
            neichgavg[p0,bi] = ny.mean(neichgtotseed[p0,bi,:])
            neichgstd[p0,bi] = ny.std(neichgtotseed[p0,bi,:])


    return Time,neichg1, neichg2, neichg12, neichgsum,neichgtotseed,neichgtot,neichgavg, neichgstd




########################################################


#We get the bending extention times

def getallexttime(nli,np0,nse,nbi,namep0,namesee,namebi,name):

    x, y, z = np0, nbi, nse

#    extvals = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]

    oldreg1 = 0
    oldreg2 = 0
    extnum = 0

#    exttimes = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]

    extvals = []

    names = [[[None for _ in range(z)] for _ in range(y) ] for _ in range(x)]


    for p0 in range (0,np0):
        for bi in range(0, nbi):
            for se in range(0, nse):

                filenme = namep0[p0] + namesee[se] + "Bi_" + namebi[bi] + name

                valslist = []

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

                                extvals.append(num1)

                                oldregI = num1

                                timeold = time
                

                        else:
                            oldregI = float(bdata[1]) 
                            extvals.append(oldregI)

                            
                        #benspeed.append(speedlist)
                    names[p0][bi][se] = filenme

    
    uniqvals = list(set(extvals))
    uniqvals.sort()

    return uniqvals, names

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

                speedlist = []

                timeslist = []

                with open(filenme, 'r') as bfile:
                    blines = bfile.readlines()


                    sdata = blines[2].split()

                    start = float(sdata[1]) 

                        
                    extvals = bleny[p0]

                    nextend = len(extvals)

                    timeold = 0

                    

                    for i1 in range(2, nli+1):
                        
                        bdata = blines[i1].split();
                        time = int(bdata[0])

                        chk = 0

            
                        if i1 > 2:
                            num1 = float(bdata[1])
                            num2 = float(bdata[2])
        

                            if num1 > oldregI:


                                for j in range(0,nextend):
                                    if num1 - extvals[j] == 0:

                                        move = (num1-oldregI)*dy/float(Nx)
                                        speed = move/(time-timeold)
            
                                        for k in range(0,j-oldindx):
                                            speedlist.append(speed)
                                            timeslist.append(time)

                                        oldindx = j

                                oldregI = num1
                                timeold = time
                

                        else:
                            oldregI = float(bdata[1]) 
                            oldregII = float(bdata[2])
                            oldindx = 0

                            
                        #benspeed.append(speedlist)
                    benspeed[p0][bi][se] = speedlist
                    exttimes[p0][bi][se] = timeslist
                    names[p0][bi][se] = filenme

    return benspeed, exttimes, names

################################################
#We average over seeds

def avgallbendspeedCFC(np0, nse, nbi, bleny, bspe):



    x, y = np0, nbi

    avgspeed = [[None for _ in range(y) ] for _ in range(x)]

    avgspeedstd = [[None for _ in range(y) ] for _ in range(x)]

#    Retained = []

    for p0 in range (0,np0):
        for bi in range(0, nbi):

            numext = len(bleny[p0])-1

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

            triseavg[p0,bi] = ny.mean(trise[p0,bi,:])

            trisestd[p0,bi] = ny.std(trise[p0,bi,:])

    return trise,triseavg, trisestd

#########################################

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
#We plot regimes maps

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

    contour1 = plt.figure(floop)
    contour = plt.contourf(Xv,Yv,Reg,levels=bound, cmap=cmap2)
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
    bihalf = Bran//2
    for bi in range(0,bihalf):
        bi2 = 2*bi
        plt.plot(xvals,Lin[:,bi],color=cmap(bi2), label=bilabel)
        plt.grid()
        plt.tick_params(axis='both', labelsize = 28)


    plt.legend(loc='lower left', bbox_to_anchor=(0.99,-0.02),fontsize=20)

    floop +=1
    return floop

###############################################################





