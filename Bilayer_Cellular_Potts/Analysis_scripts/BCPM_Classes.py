#This will hold classes used for analyzing data

class Edgematchdat:
#    def __init__(self, tdat, edat, ediff,edstd,erise,eriseavg,eristd,nam,edatavg,edatstd):
    def __init__(self):

        self.tsdat = 0
        self.edgedat = 0
        self.edgediff = 0
        self.edgestd = 0
        self.edgerise = 0
        self.edgeriseavg = 0
        self.edgeristd = 0
        self.edgedatavg = 0
        self.edgedatstd = 0
        self.name = "string"


#        self.tsdat = tdat
#        self.edgedat = edat
#        self.edgediff = ediff
#        self.edgestd = edstd
#        self.edgerise = erise
#        self.edgeriseavg = eriseavg
#        self.edgeristd = eristd
#        self.edgedatavg = edatavg
#        self.edgedatstd = edatstd
#        self.name = nam


class Edgedatavg:
    def __init__(self):
        
        self.tsdat = 0
        self.edge = 0
        self.edgestd = 0

        #self.tsdat = np.array(tdat)
        #self.edge = np.array(edat)
        #self.edgestd = np.array(edstd)


class Neighdata:
    def __init__(self):
        self.time = 0
        self.neichg1 = 0
        self.neichg2 = 0
        self.neichg12 = 0
        self.neichgtotseed = 0
        self.neichgtot = 0
        self.neichgsum = 0
        self.neichgavg = 0
        self.neichgstd = 0
        self.name = "string"


class Edgefreq:
    
    def __init__(self):

        self.edgebins = 0
        self.edgebinsnum = 0
        self.edgebarrseed = 0
        self.edgebarr = 0
        self.edgeseedsfreq = 0
        self.edgefinal = 0
        self.nummins = 0
        self.highmatch = 0
        self.modematch = 0
        self.fullprop = 0
        self.neifreqseed = 0
        self.neifreq = 0


    #def __init__(self, edgef, edgedf, edgeEsee, edgeE, edgesen, numin, edgemc, edgemd, fp, nfresee, nfre):
        
        #self.edgebins = edgef
        #self.edgebinsnum = edgedf
        #self.edgebarrseed = edgeEsee
        #self.edgebarr = edgeE
        #self.edgeseedsfreq = edgesen
        #self.nummins = numin
        #self.highmatch = edgemc
        #self.modematch = edgemd
        #self.fullprop = fp
        #self.neifreqseed = nfresee
        #self.neifreq = nfre

class Energymins:
    def __init__(self):
        self.engymin = []

    def getmins(self, engy):
        self.engymin = engy


##########################################################

class Shapedata:
    def __init__(self):
        self.Pl1 = 0
        self.Pl2 = 0
        self.p0l1 = 0
        self.p0l2 = 0
        self.Pavg = 0
        self.Pstd = 0
        self.p0avg = 0
        self.p0std = 0



class trajdata:
    def __init__(self):
        self.timedat = 0
        self.xrawl1 = 0
        self.yrawl1 = 0
        self.xrawl2 = 0
        self.yrawl2 = 0
        self.delxl1 = 0
        self.delyl1 = 0
        self.delxl2 = 0
        self.delxl2 = 0
        self.xmsd1 = 0
        self.ymsd1 = 0
        self.xmsd2 = 0
        self.ymsd2 = 0
        self.xdist1 = 0
        self.ydist1 = 0
        self.xdist2 = 0
        self.ydist2 = 0
        self.name = "string"

class trajdatarise:
    def __init__(self):
        self.timedat = 0
        self.xrawl1 = 0
        self.yrawl1 = 0
        self.xrawl2 = 0
        self.yrawl2 = 0
        self.xmsd1 = 0
        self.ymsd1 = 0
        self.xmsd2 = 0
        self.ymsd2 = 0
        self.xmsd1seed = 0
        self.ymsd1seed = 0
        self.xmsd2seed = 0
        self.xmsd2seed = 0
        self.rmsd = 0
        self.rmsdseed = 0
        self.xdist1 = 0
        self.ydist1 = 0
        self.xdist2 = 0
        self.ydist2 = 0
        self.name = "string"


class trajafterrise:
    def __init__(self):
        self.timedat = 0
        self.xrawl1 = 0
        self.yrawl1 = 0
        self.xrawl2 = 0
        self.yrawl2 = 0
        self.xmsd1 = 0
        self.ymsd1 = 0
        self.xmsd2 = 0
        self.ymsd2 = 0
        self.xmsd1seed = 0
        self.ymsd1seed = 0
        self.xmsd2seed = 0
        self.xmsd2seed = 0
        self.rmsd = 0
        self.rmsdseed = 0
        self.xdist1 = 0
        self.ydist1 = 0
        self.xdist2 = 0
        self.ydist2 = 0
        self.name = "string"



class MSDanddistavg:
    def __init__(self):
        self.time = 0
        self.xmsdseed = 0
        self.ymsdseed = 0
        self.rmsdseed = 0
        self.xdistseed = 0
        self.ydistseed = 0
        self.rdistseed = 0
        self.xmsd = 0
        self.ymsd = 0
        self.rmsd = 0
        self.xdist = 0
        self.ydist = 0
        self.rdist = 0

class Datapoints:
    def __init__(self):
        self.diffcoeff = 0
        self.diffcoeffavg = 0
        self.diffcoeffstd = 0
        self.findist = 0
        self.findistavg = 0
        self.findiststd = 0
    
class Speeddata:
    def __init__(self):
        self.speed = 0
        self.speedavg = 0
        self.name = "string"




class Propsatrise:
    def __init__(self):
        self.distatriseseed = 0
        self.distatrise = 0
        self.distatrisestd = 0
        self.speedatriseseed = 0
        self.speedatrise = 0
        self.speedatrisestd = 0
        self.speedatrisealt = 0
        self.neiatriseseed = 0
        self.neiatrise = 0
        self.neiatrisestd = 0
        self.msdatriseseed = 0
        self.msdatrise = 0
        self.msdatrisestd = 0
        self.fitcoeff = 0


class Propsafterrise:
    def __init__(self):
        self.distafterriseseed = 0
        self.distafterrise = 0
        self.distafterrisestd = 0
        self.speedafterriseseed = 0
        self.speedafterrise = 0
        self.speedafterrisestd = 0
        self.speedafterrisealt = 0
        self.neiafterriseseed = 0
        self.neiafterrise = 0
        self.neiafterrisestd = 0
        self.msdafterriseseed = 0
        self.msdafterrise = 0
        self.msdafterrisestd = 0
        self.fitcoeff = 0


##############################################################


class bendinfoseeds:

    def __init__(self):
        self.numex = 0
        self.benlist = 0
        self.benuplist = 0
        self.bendownlist = 0
        self.bentime = 0
        self.benuptime = 0
        self.benddowntime = 0
        self.bendiff = 0
        self.bspeed = 0
        self.bupspeed = 0
        self.bdownspeed = 0
        self.regImatch = 0
        self.regIImatch = 0
        self.tsshort = 0
        self.numexup = 0
        self.numexdown = 0
        self.name = "string"


class bendspeedavg:

    def __init__(self):
        self.bspeed = []

    def add_speed(self,spe):
        #self.labelsnum.append(lab)
        self.bspeed = spe

class bendspeedstd:

    def __init__(self):
        self.speedstd = []

    def add_speed(self,spestd):
        #self.labelsnum.append(lab)
        self.speedstd = spestd

class bendinfofreq:

    def __init__(self):
        self.numex = 0
        self.bins = 0
        self.binlabel = []
        self.benspeedupdowncollect = []
        self.benspeedfreqseed = 0
        self.benspeedfreq = 0
        self.extendfreq = 0
        self.bspeed = 0
        self.speedmax = 0
        self.bendupaccelprob = 0
        self.benddownaccelprob = 0
        self.extprob = []
        self.nextextprob = []


class bendinfoavg:

    def __init__(self):
        self.numex = 0
        self.bentime = 0
        self.bendiff = 0
        self.bspeed = 0
        self.speedstd = 0
#        self.speedstd = []

#    def get_bspeed(speed):
#        self.bspeed.append(speed)

#    def get_speedstd(self, std):
#        self.speedstd.append(std)

class bendlabels:

    def __init__(self):
        self.labelsnum = []
    
    def add_label(self,lab):
        #self.labelsnum.append(lab)
        self.labelsnum = lab


class benduplabels:

    def __init__(self):
        self.labelsnum = []

    def add_label(self,lab):
        #self.labelsnum.append(lab)
        self.labelsnum = lab

class benddownlabels:

    def __init__(self):
        self.labelsnum = []

    def add_label(self,lab):
        #self.labelsnum.append(lab)
        self.labelsnum = lab



########################################################        
