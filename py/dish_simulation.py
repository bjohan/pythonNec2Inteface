import band
import compute
import necFileGenerator
import plothelper
import numpy as np

class Reflector:
    def __init__(self, x_bottom, dia, f_over_d, mesh_res = 0.01):
        self.dia = dia;
        self.f_over_d = f_over_d
        self.f = self.f_over_d*self.dia
        self.mesh_res = mesh_res
        self.x_bottom = x_bottom
    
    def getX(self, x, y):
        r = np.sqrt(x*x+y*y)
        return self.f - r*r/(4*self.f)
        #f - r*r/(4*f) = x

    def addWire(self, fg, x1, y1, z1, x2, y2, z2, dia):
        #c = 299792458
        #dx = x2-x1
        #dy = y2-y1
        #dz = z2-z1
        #l = np.sqrt(dx**2 + dy**2 + dz**2)
        #lseg = (c/self.fh)/10
        #nseg = int(l/lseg)
        #if minSegments > 0:
        #    if nseg < minSegments:
        #        nseg=minSegments
        #fg.wire(x1,y1,z1, x2,y2,z2, dia/2, segments = nseg)
        fg.wire(x1,y1,z1, x2,y2,z2, dia/2, segments = 1)
        #return nseg


    def addNecGeometry(self, fg):
        ng = int(self.dia/self.mesh_res+0.5)
        grid = np.linspace(-self.dia/2, self.dia/2, ng)
        for yi in range(ng-1):
            for zi in range(ng-1):
                y = grid[yi]
                z = grid[zi]
                x = self.getX(z,y)+self.x_bottom

                yp = grid[yi+1]
                xyp = self.getX(z,yp)+self.x_bottom

                zp = grid[zi+1]
                xzp = self.getX(zp,y)+self.x_bottom
                r = np.sqrt(y*y+z*z)
                if r < self.dia/2:# and r > 0.35:
                    self.addWire(fg, x, y, z, xyp, yp, z, 0.001)
                    self.addWire(fg, x, y, z, xzp, y, zp, 0.001)


class Logperiodic:
    def __init__(self, dia, frequencyLow, frequencyHigh, sigma, tau, zc=50):
        self.name = "Log periodic"
        self.dia = dia
        self.fl = frequencyLow
        self.fh = frequencyHigh
        self.s = sigma
        self.t = tau
        self.zc = zc

    def addWire(self, fg, x1, y1, z1, x2, y2, z2, dia, minSegments=1):
        c = 299792458
        dx = x2-x1
        dy = y2-y1
        dz = z2-z1
        l = np.sqrt(dx**2 + dy**2 + dz**2)
        lseg = (c/self.fh)/10
        nseg = int(l/lseg)
        if minSegments > 0:
            if nseg < minSegments:
                nseg=minSegments
        fg.wire(x1,y1,z1, x2,y2,z2, dia/2, segments = nseg)
        return nseg


    def addNecGeometry(self, fg):
        B = self.fh/self.fl
        print("Relative bandwidth: %f"%(B))

        if self.t < 0.8 or self.t > 0.98:
            raise ValueError("Tau must be in range 0.8 to 0.98 and not %f"%(self.t))
        sopt = 0.243*self.t - 0.051
        print("Sigma optimum %f"%(sopt))

        if self.s < 0.03 or self.s > sopt:
            raise ValueError("Sigma must be in range 0.03 to %f and not %f"%(sopt, self.s))

        cotAlph = 4.0*self.s*(1.0-self.t)**-1
        print("cotalph %f"%(cotAlph))
        
        Bar = 1.1+7.7*((1.0-self.t)**2)*cotAlph
        print("B_ar %f"%(Bar))
    
        Bs = B*Bar
        print("B_s %f"%(Bs))


        N = 1+np.log(Bs)/np.log(1/self.t)
        print("N %f"%(N))
        if (N-np.floor(N)) < 0.3:
            N = int(np.floor(N))
        else:
            N = int(np.ceil(N))
        print("N %d"%(N))

        c = 299792458
        print("c %f"%(c))

        l1 = c/(2*self.fl)
        print("l1 %f"%(l1))

        ls = [l1]
        for i in range(N-1):
            li = self.t*ls[-1]
            ls.append(li)
            print("l%d %f"%(i+2, li))


        ds = []
        for i in range(N-1):
            di = cotAlph*(ls[i]-ls[i+1])/2
            ds.append(di)
            print("d%d,%d %f"%(i+1,i+2, di))

        print("L", np.sum(ds))

        lTerm = (c/self.fl)/8
        print("L_term %f"%(lTerm))

        zcn = 120*(np.log(ls[-1]/self.dia)-2.25)
        sp = self.s/np.sqrt(self.t)

        print("Zcn %f"%(zcn))
        
        print("sp %f"%(sp))

        zcfeed = (self.zc**2)/(8*sp*zcn)+self.zc*np.sqrt((self.zc/(8*sp*zcn))**2+1)
        print("Zc feed %f"%(zcfeed))

        mu = 4*np.pi*1e-7
        eps = 1/(mu*c*c)
        zfs = mu*c
        print("mu %E eps %E Zfs %f dia %f"%(mu, eps, zfs, self.dia))
        D = 0.001098*1.5;
        ztl = np.arccosh(D/self.dia)*zfs/(np.pi*np.sqrt(1.00054))
        print("Boom transmission line impedance %f"%(ztl))

        self.addWire(fg, 0,0,0, -lTerm,0,0, self.dia)
        self.addWire(fg, 0,0,D, -lTerm,0,D, self.dia)
        self.addWire(fg, -lTerm,0, 0, -lTerm, 0, D, self.dia)


        x = 0;
        self.addWire(fg, x, 0, 0,  x,ls[0]/2, 0, self.dia)
        for l, d,n in zip(ls[1:], ds, range(len(ds))):
            self.addWire(fg, x, 0, 0, x+d, 0, 0, self.dia)
            x+=d;
            sign = 1
            if not n%2:
                sign = -1
            self.addWire(fg, x, 0, 0, x, sign*l/2, 0, self.dia)

        x = 0;
        self.addWire(fg, x, 0, D,  x,-ls[0]/2, D, self.dia)
        for l, d, n in zip(ls[1:], ds, range(len(ds))):
            self.addWire(fg, x, 0, D, x+d, 0, D, self.dia)
            x+=d;
            sign = 1
            if n %2:
                sign = -1
            self.addWire(fg, x, 0, D, x, sign*l/2, D, self.dia)

        nseg = self.addWire(fg, x, 0, 0, x, 0, D, self.dia, minSegments=5)
        print("Exciting tag %d"%(fg.tag-1))
        fg.excite(fg.tag-1, int(nseg/2)+1)
 

        r = Reflector(0.1+0.07+0*np.sum(ds)+0.38*1.2*0, 0.3, 0.38, 0.02)
        print(r.f)
        #r.addNecGeometry(fg)
        #dishOffs = 0.5
        #dishD = 2.4
        #dishRes = 0.04
        #nGrid = int(dishD/dishRes)
        #grid = np.linspace(0, dishD, nGrid)-dishD/2

        #for zi in range(nGrid-1):
        #    for yi in range(nGrid-1):
        #        self.addWire(fg, dishOffs, grid[yi], grid[zi], dishOffs, grid[yi+1], grid[zi], 0.001)
        #        self.addWire(fg, dishOffs, grid[yi], grid[zi], dishOffs, grid[yi], grid[zi+1], 0.001)
                
            
        print("Geom done")
        return fg;
    
f0 = 850
#f1 = 26500
f1 = 1300
band = band.Band(f0, f1, 'entire')
#band = band.Band(f0, 15000, "low")
l = Logperiodic(0.001, f0*1e6*0.8, f1*1e6*1.2, 0.09, 0.88)
cpt = compute.Compute()
cpt.setAntenna(l)
cpt.addBands([band])
result = cpt.compute(11)
azg = result.frequencies[0].getAzimuthGrid()
elg = result.frequencies[0].getElevationGrid()
dbg = result.frequencies[0].getMajorDbGrid()
ph = plothelper.PlotHelper()
ph.plotColorMap(azg, elg, dbg)
