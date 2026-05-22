#This will hold classes used for analyzing data

class Edgematchdat:
    def _init_(self, tdat, edat, ediff,edstd,erise,eriseavg,eristd,nam,edatavg,edatstd):
        self.tsdat = tdat
        self.edgedat = edat
        self.edgediff = ediff
        self.edgestd = edstd
        self.edgerise = erise
        self.edgeriseavg = eriseavg
        self.edgeristd = eristd
        self.edgedatavg = edatavg
        self.edgedatstd = edatstd
        self.name = nam


class Edgedatavg:
    def _init_(self, tdat, edat, edstd):
        self.tsdat = np.array(tdat)
        self.edge = np.array(edat)
        self.edgestd = np.array(edstd)


class Neighdata:
    def _init_(self,tim,nei1,nei2,nei12,neisum,neisumtot,neised,neiboth,neibothstd):
        self.time = tim
        self.neichg1 = nei1
        self.neichg2 = nei2
        self.neichg12 = nei12
        self.neichgtotseed = neisum
        self.neichgtot = neisumtot
        self.neichgsum = neised
        self.neichgavg = neiboth
        self.neichgstd = neibothstd


##########################################################

class Shapedata:
    def _init_(self,P1,P2,p01,p02,P12,P12std,p012,p012std):
        self.Pl1 = P1
        self.Pl2 = P2
        self.p0l1 = p01
        self.p0l2 = p02
        self.Pavg = P12
        self.Pstd = P12std
        self.p0avg = p012
        self.p0std = p012std



class trajdata:
    def _init_(self,tdat,xrl1,yrl1,xrl2,yrl2,xdiffl1,ydiffl1,xdiffl2,ydiffl2,Xmsd1,Ymsd1,Xmsd2,Ymsd2,Xdist1,Ydist1,Xdist2,Ydist2,nei1,nei2,nam):
        self.timedat = tdat
        self.xrawl1 = xrl1
        self.yrawl1 = yrl1
        self.xrawl2 = xrl2
        self.yrawl2 = yrl2
        self.delxl1 = xdiffl1
        self.delyl1 = ydiffl1
        self.delxl2 = xdiffl2
        self.delxl2 = ydiffl2
        self.xmsd1 = Xmsd1
        self.ymsd1 = Ymsd1
        self.xmsd2 = Xmsd2
        self.ymsd2 = Ymsd2
        self.xdist1 = Xdist1
        self.ydist1 = Ydist1
        self.xdist2 = Xdist2
        self.ydist2 = Ydist2
        self.name = nam



class MSDanddistavg:
    def _init_(self, tim,Xmsdsee, Ymsdsee, Rmsdsee, Xdistsee, Ydistsee, Rdistsee,Xmsd, Ymsd, Rmsd, Xdist, Ydist, Rdist, Rrisedist):
        self.time = tim
        self.xmsdseed = Xmsdsee
        self.ymsdseed = Ymsdsee
        self.rmsdseed = Rmsdsee
        self.xdistseed = Xdistsee
        self.ydistseed = Ydistsee
        self.rdistseed = Rdistsee
        self.xmsd = Xmsd
        self.ymsd = Ymsd
        self.rmsd = Rmsd
        self.xdist = Xdist
        self.ydist = Ydist
        self.rdist = Rdist


class Datapoints:
    def _init_(self, dco, dcoavg, dcostd, dis, disavg, disstd):
        self.diffcoeff = dco
        self.diffcoeffavg = dcoavg
        self.diffcoeffstd = dcostd
        self.findist = dis
        self.findistavg = disavg
        self.findiststd = disstd
    

class Speeddata:
    def _init_(self, spe,speavg,nam):
        self.speed = spe
        self.speedavg = speavg
        self.name = nam


class Propsatrise:
    def _init_(self, ratrisseed,ratris, ratrisstd, satrisseed,satris,satrisstd):
        self.distatriseseed = ratrisseed
        self.distatrise = ratris
        self.distatrisestd = ratrisstd
        self.speedatriseseed = satrisseed
        self.speedatrise = satris
        self.speedatrisestd = satrisstd

##############################################################


class bendvalues:
    def _init_(self, vals, nam):
        self.regIpnts = vals
        self.name = nam


class bendinfoseeds:

    def _init_(self, nex, bendt, bendd, bspee, nam):
        self.numex = nex
        self.bentime = bendt
        self.bendiff = bendd
        self.bspeed = bspee
        self.name = nam



class bendinfoavg:

    def _init_(self, nex, bendt, bendd, bspee, spestd):
        self.numex = nex
        self.bentime = bendt
        self.bendiff = bendd
        self.bspeed = bspee
        self.speedstd = spestd


########################################################        
