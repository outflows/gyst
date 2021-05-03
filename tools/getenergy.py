import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
from scipy import interpolate
import tqdm
import nmmn.lsd, nmmn.misc
import os.path
import pandas as pd
import astropy.visualization
import glob
import matplotlib as mpl

mpl.rcParams['font.serif'] = 'helvetica'
mpl.rcParams['font.sans-serif'] = 'helvetica'

path_to_dumps = "/home/physics/rodrgust/harmpi/MAD52_beta80/dumps/"
path_to_run = "/home/physics/rodrgust/gyst/run/"

def mdot(a,b):
    """
    Computes a contraction of two tensors/vectors.  Assumes
    the following structure: tensor[m,n,i,j,k] OR vector[m,i,j,k],
    where i,j,k are spatial indices and m,n are variable indices.
    """
    if (a.ndim == 3 and b.ndim == 3) or (a.ndim == 4 and b.ndim == 4):
          c = (a*b).sum(0)
    elif a.ndim == 5 and b.ndim == 4:
          c = np.empty(np.maximum(a[:,0,:,:,:].shape,b.shape),dtype=b.dtype)
          for i in range(a.shape[0]):
                c[i,:,:,:] = (a[i,:,:,:,:]*b).sum(0)
    elif a.ndim == 4 and b.ndim == 5:
          c = np.empty(np.maximum(b[0,:,:,:,:].shape,a.shape),dtype=a.dtype)
          for i in range(b.shape[1]):
                c[i,:,:,:] = (a*b[:,i,:,:,:]).sum(0)
    elif a.ndim == 5 and b.ndim == 5:
          c = np.empty((a.shape[0],b.shape[1],a.shape[2],a.shape[3],max(a.shape[4],b.shape[4])),dtype=a.dtype)
          for i in range(c.shape[0]):
                for j in range(c.shape[1]):
                      c[i,j,:,:,:] = (a[i,:,:,:,:]*b[:,j,:,:,:]).sum(0)
    elif a.ndim == 5 and b.ndim == 6:
          c = np.empty((a.shape[0],b.shape[1],b.shape[2],max(a.shape[2],b.shape[3]),max(a.shape[3],b.shape[4]),max(a.shape[4],b.shape[5])),dtype=a.dtype)
          for mu in range(c.shape[0]):
              for k in range(c.shape[1]):
                  for l in range(c.shape[2]):
                      c[mu,k,l,:,:,:] = (a[mu,:,:,:,:]*b[:,k,l,:,:,:]).sum(0)
    else:
           raise Exception('mdot', 'wrong dimensions')
    return c

def myfloat(f,acc=1):
    """ acc=1 means np.float32, acc=2 means np.float64 """
    if acc==1:
        return( np.float32(f) )
    else:
        return( np.float64(f) )

def rd(dump):
    read_file(dump,type="dump")

def read_file(dump,type=None,savedump=True,saverdump=False,noround=False):
    if type is None:
        if dump.startswith("dump"):
            type = "dump"
            print("Reading a dump file %s ..." % dump)
        elif dump.startswith("gdump2"):
            type = "gdump2"
            print("Reading a gdump2 file %s ..." % dump)
        elif dump.startswith("gdump"):
            type = "gdump"
            print("Reading a gdump file %s ..." % dump)
        elif dump.startswith("rdump"):
            type = "rdump"
            print("Reading a rdump file %s ..." % dump)
        elif dump.startswith("fdump"):
            type = "fdump"
            print("Reading a fdump file %s ..." % dump)
        else:
            print("Couldn't guess dump type; assuming it is a data dump")
            type = "dump"
    #normal dump
    if os.path.isfile( path_to_dumps + dump ):
        headerline = read_header(path_to_dumps + dump, returnheaderline = True)
        gd = read_body(path_to_dumps + dump,nx=N1+2*N1G,ny=N2+2*N2G,nz=N3+2*N3G,noround=1)
        if noround:
            res = data_assign(         gd,type=type,nx=N1+2*N1G,ny=N2+2*N2G,nz=N3+2*N3G)
        else:
            res = data_assign(myfloat(gd),type=type,nx=N1+2*N1G,ny=N2+2*N2G,nz=N3+2*N3G)
        return res
    #MPI-type dump that is spread over many files
    else:
        flist = np.sort(glob.glob( path_to_dumps + dump + "_[0-9][0-9][0-9][0-9]" ))
        if len(flist) == 0:
            print( "Could not find %s or its MPI counterpart" % dump )
            return
        sys.stdout.write( "Reading %s (%d files)" % (dump, len(flist)) )
        sys.stdout.flush()
        ndots = 10
        dndot = len(flist)/ndots
        if dndot == 0: dndot = 1
        for i,fname in enumerate(flist):
            #print( "Reading file %d out of %d..." % (i,len(flist)) )
            #header for each file might be different, so read each
            header = read_header(fname,issilent=1)
            if header is None:
                print( "Error reading header of %s, aborting..." % fname )
                return
            lgd = read_body(fname,nx=N1+2*N1G,ny=N2+2*N2G,nz=N3+2*N3G)
            #this gives an array of dimensions (-1,N1,N2,N3)+potentially ghost cells
            if 0 == i:
                #create full array: of dimensions, (-1,nx,ny,nz)
                fgd = np.zeros( (lgd.shape[0], nx+2*N1G, ny+2*N2G, nz+2*N3G), dtype=np.float32)
            if not type == "rdump":
                #construct full indices: ti, tj, tk
                #fti,ftj,ftk = mgrid[0:nx,0:ny,0:nz]
                lti,ltj,ltk = lgd[0:3,:,:].view();
                lti = np.int64(lti)
                ltj = np.int64(ltj)
                ltk = np.int64(ltk)
                fgd[:,lti+N1G,ltj+N2G,ltk+N3G] = lgd[:,:,:,:]
            else:
                print(starti,startj,startk)
                fgd[:,starti:starti+N1+2*N1G,startj:startj+N2+2*N2G,startk:startk+N3+2*N3G] = lgd[:,:,:,:]
            del lgd
            if i%dndot == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
        res = data_assign(fgd,type=type,nx=nx+2*N1G,ny=ny+2*N2G,nz=nz+2*N3G)
        if savedump:
            #if the full dump file does not exist, create it
            dumpfullname = path_to_dumps + dump
            if (type == "dump" or type == "gdump") and not os.path.isfile(dumpfullname):
                sys.stdout.write("Saving full dump to %s..." % dumpfullname)
                sys.stdout.flush()
                header[1] = header[4] #N1 = nx
                header[2] = header[5] #N2 = ny
                header[3] = header[6] #N3 = nz
                fout = open( dumpfullname, "wb" )
                #join header items with " " (space) as a glue
                #see http://stackoverflow.com/questions/12377473/python-write-versus-writelines-and-concatenated-strings
                #write it out with a new line char at the end
                fout.write(" ".join(header) + "\n")
                fout.flush()
                os.fsync(fout.fileno())
                #reshape the dump content
                gd1 = fgd.transpose(1,2,3,0)
                gd1.tofile(fout)
                fout.close()
                print( " done!" )
                if res is not None:
                    return res
        return res

#read in a header
def read_header(dump,issilent=True,returnheaderline=False):
    global t,nx,ny,nz,N1,N2,N3,N1G,N2G,N3G,starti,startj,startk,_dx1,_dx2,_dx3,a,gam,Rin,Rout,hslope,R0,ti,tj,tk,x1,x2,x3,r,h,ph,gcov,gcon,gdet,drdx,gn3,gv3,guu,gdd,dxdxp, games, startx1, startx2, startx3, x10, x20, tf, NPR, DOKTOT, BL,dt
    global fractheta
    global fracphi
    global rbr
    global npow2
    global cpow2
    #read image
    fin = open( dump, "rb" )
    headerline = fin.readline()
    header = headerline.split()
    nheadertot = len(header)
    fin.close()
    if not dump.startswith(path_to_dumps + "rdump"):
        if not issilent: print( "dump header: len(header) = %d" % len(header) )
        nheader = 57
        n = 0
        t = myfloat(np.float64(header[n])); n+=1
        #per tile resolution
        N1 = int(header[n]); n+=1
        N2 = int(header[n]); n+=1
        N3 = int(header[n]); n+=1
        #total resolution
        nx = int(header[n]); n+=1
        ny = int(header[n]); n+=1
        nz = int(header[n]); n+=1
        #numbers of ghost cells
        N1G = int(header[n]); n+=1
        N2G = int(header[n]); n+=1
        N3G = int(header[n]); n+=1
        startx1 = myfloat(float(header[n])); n+=1
        startx2 = myfloat(float(header[n])); n+=1
        startx3 = myfloat(float(header[n])); n+=1
        _dx1=myfloat(float(header[n])); n+=1
        _dx2=myfloat(float(header[n])); n+=1
        _dx3=myfloat(float(header[n])); n+=1
        tf=myfloat(float(header[n])); n+=1
        nstep=myfloat(float(header[n])); n+=1
        a=myfloat(float(header[n])); n+=1
        gam=myfloat(float(header[n])); n+=1
        cour=myfloat(float(header[n])); n+=1
        DTd=myfloat(float(header[n])); n+=1
        DTl=myfloat(float(header[n])); n+=1
        DTi=myfloat(float(header[n])); n+=1
        DTr=myfloat(float(header[n])); n+=1
        DTr01=myfloat(float(header[n])); n+=1
        dump_cnt=myfloat(float(header[n])); n+=1
        image_cnt=myfloat(float(header[n])); n+=1
        rdump_cnt=myfloat(float(header[n])); n+=1
        rdump01_cnt=myfloat(float(header[n])); n+=1
        dt=myfloat(float(header[n])); n+=1
        lim=myfloat(float(header[n])); n+=1
        failed=myfloat(float(header[n])); n+=1
        Rin=myfloat(float(header[n])); n+=1
        Rout=myfloat(float(header[n])); n+=1
        hslope=myfloat(float(header[n])); n+=1
        R0=myfloat(float(header[n])); n+=1
        NPR=int(header[n]); n+=1
        DOKTOT=int(header[n]); n+=1
        DOCYLINDRIFYCOORDS=int(header[n]); n+=1
        fractheta = myfloat(header[n]); n+=1
        fracphi   = myfloat(header[n]); n+=1
        rbr       = myfloat(header[n]); n+=1
        npow2     = myfloat(header[n]); n+=1
        cpow2     = myfloat(header[n]); n+=1
        x10 = myfloat(header[n]); n+=1
        x20 = myfloat(header[n]); n+=1
        fracdisk = myfloat(header[n]); n+=1
        fracjet = myfloat(header[n]); n+=1
        r0disk = myfloat(header[n]); n+=1
        rdiskend = myfloat(header[n]); n+=1
        r0jet = myfloat(header[n]); n+=1
        rjetend = myfloat(header[n]); n+=1
        jetnu = myfloat(header[n]); n+=1
        rsjet = myfloat(header[n]); n+=1
        r0grid = myfloat(header[n]); n+=1
        BL = myfloat(header[n]); n+=1
    else:
        print("rdump header")
        nheader = 48
        n = 0
        #per tile resolution
        N1 = int(header[n]); n+=1
        N2 = int(header[n]); n+=1
        N3 = int(header[n]); n+=1
        #total resolution
        nx = int(header[n]); n+=1
        ny = int(header[n]); n+=1
        nz = int(header[n]); n+=1
        #numbers of ghost cells
        N1G = int(header[n]); n+=1
        N2G = int(header[n]); n+=1
        N3G = int(header[n]); n+=1
        #starting indices
        starti = int(header[n]); n+=1
        startj = int(header[n]); n+=1
        startk = int(header[n]); n+=1
        t = myfloat(header[n]); n+=1
        tf = myfloat(header[n]); n+=1
        nstep = int(header[n]); n+=1
        a = myfloat(header[n]); n+=1
        gam = myfloat(header[n]); n+=1
        game = myfloat(header[n]); n+=1
        game4 = myfloat(header[n]); n+=1
        game5 = myfloat(header[n]); n+=1
        cour = myfloat(header[n]); n+=1
        DTd = myfloat(header[n]); n+=1
        DTl = myfloat(header[n]); n+=1
        DTi = myfloat(header[n]); n+=1
        DTr = myfloat(header[n]); n+=1
        DTr01 = myfloat(header[n]); n+=1
        dump_cnt = myfloat(header[n]); n+=1
        image_cnt = myfloat(header[n]); n+=1
        rdump_cnt = myfloat(header[n]); n+=1
        rdump01_cnt=myfloat(float(header[n])); n+=1
        dt = myfloat(header[n]); n+=1
        lim = myfloat(header[n]); n+=1
        failed = myfloat(header[n]); n+=1
        Rin = myfloat(header[n]); n+=1
        Rout = myfloat(header[n]); n+=1
        hslope = myfloat(header[n]); n+=1
        R0 = myfloat(header[n]); n+=1
        fractheta = myfloat(header[n]); n+=1
        fracphi = myfloat(header[n]); n+=1
        rbr = myfloat(header[n]); n+=1
        npow2 = myfloat(header[n]); n+=1
        cpow2 = myfloat(header[n]); n+=1
        x10 = myfloat(header[n]); n+=1
        x20 = myfloat(header[n]); n+=1
        mrat = myfloat(header[n]); n+=1
        fel0 = myfloat(header[n]); n+=1
        felfloor = myfloat(header[n]); n+=1
        tdump = myfloat(header[n]); n+=1
        trdump = myfloat(header[n]); n+=1
        timage = myfloat(header[n]); n+=1
        tlog  = myfloat(header[n]); n+=1
    if n < len(header):
        nheader = 60
        global_fracdisk   = myfloat(header[n]); n+=1
        global_fracjet    = myfloat(header[n]); n+=1
        global_r0disk     = myfloat(header[n]); n+=1
        global_rdiskend   = myfloat(header[n]); n+=1
        global_r0jet      = myfloat(header[n]); n+=1
        global_rjetend    = myfloat(header[n]); n+=1
        global_jetnu      = myfloat(header[n]); n+=1
        global_rsjet      = myfloat(header[n]); n+=1
        global_r0grid     = myfloat(header[n]); n+=1
    if n != nheader or n != nheadertot:
        print("Wrong number of elements in header: nread = %d, nexpected = %d, nototal = %d: incorrect format?"
              % (n, nheader, nheadertot) )
        return headerline
    if returnheaderline:
        return headerline
    else:
        return header

def read_body(dump,nx=None,ny=None,nz=None,noround=False):
        fin = open( dump, "rb" )
        header = fin.readline()
        if dump.startswith(path_to_dumps + "rdump"):
            dtype = np.float64
            body = np.fromfile(fin,dtype=dtype,count=-1)
            gd = body.view().reshape((nx,ny,nz,-1), order='C')
            if noround:
                gd=gd.transpose(3,0,1,2)
            else:
                gd=myfloat(gd.transpose(3,0,1,2))
        elif dump.startswith(path_to_dumps + "gdump2"):
            dtype = np.float64
            body = np.fromfile(fin,dtype=dtype,count=-1)
            gd = body.view().reshape((nx,ny,nz,-1), order='C')
            if noround:
                gd=gd.transpose(3,0,1,2)
            else:
                gd=myfloat(gd.transpose(3,0,1,2))
        elif dump.startswith(path_to_dumps + "fdump"):
            dtype = np.int64
            body = np.fromfile(fin,dtype=dtype,count=-1)
            gd = body.view().reshape((-1,nz,ny,nx), order='F')
            gd=myfloat(gd.transpose(0,3,2,1))
        else:
            dtype = np.float32
            body = np.fromfile(fin,dtype=dtype,count=-1)
            gd = body.view().reshape((-1,nz,ny,nx), order='F')
            gd=myfloat(gd.transpose(0,3,2,1))
        return gd

def data_assign(gd,type=None,**kwargs):
    if type is None:
        print("Please specify data type")
        return
    if type == "gdump":
        gdump_assign(gd,**kwargs)
        return None
    elif type == "gdump2":
        gdump2_assign(gd,**kwargs)
        return None
    elif type == "dump":
        dump_assign(gd,**kwargs)
        return None
    elif type == "rdump":
        gd = rdump_assign(gd,**kwargs)
        return gd
    elif type == "fdump":
        gd = fdump_assign(gd,**kwargs)
        return gd
    else:
        print("Unknown data type: %s" % type)
        return gd

def gdump_assign(gd,**kwargs):
    global t,nx,ny,nz,N1,N2,N3,_dx1,_dx2,_dx3,a,gam,Rin,Rout,hslope,R0,ti,tj,tk,x1,x2,x3,r,h,ph,gcov,gcon,gdet,drdx,gn3,gv3,guu,gdd,dxdxp, games
    nx = kwargs.pop("nx",nx)
    ny = kwargs.pop("ny",ny)
    nz = kwargs.pop("nz",nz)
    ti,tj,tk,x1,x2,x3,r,h,ph = gd[0:9,:,:].view();  n = 9
    gv3 = gd[n:n+16].view().reshape((4,4,nx,ny,nz),order='F').transpose(1,0,2,3,4); n+=16
    gn3 = gd[n:n+16].view().reshape((4,4,nx,ny,nz),order='F').transpose(1,0,2,3,4); n+=16
    gcov = gv3
    gcon = gn3
    guu = gn3
    gdd = gv3
    gdet = gd[n]; n+=1
    drdx = gd[n:n+16].view().reshape((4,4,nx,ny,nz),order='F').transpose(1,0,2,3,4); n+=16
    dxdxp = drdx
    if n != gd.shape[0]:
        print("rd: WARNING: nread = %d < ntot = %d: incorrect format?" % (n, gd.shape[0]) )
        return 1
    return 0

def gdump2_assign(gd,**kwargs):
    global t,nx,ny,nz,N1,N2,N3,_dx1,_dx2,_dx3,a,gam,Rin,Rout,hslope,R0,ti,tj,tk,x1,x2,x3,gdet,games,rf1,hf1,phf1,rf2,hf2,phf2,rf3,hf3,phf3,rcorn,hcord,phcorn,re1,he1,phe1,re2,he2,phe2,re3,he3,phe3
    nx = kwargs.pop("nx",nx)
    ny = kwargs.pop("ny",ny)
    nz = kwargs.pop("nz",nz)
    ti,tj,tk,x1,x2,x3 = gd[0:6,:,:].view();  n = 6
    rf1,hf1,phf1,rf2,hf2,phf2,rf3,hf3,phf3 = gd[0:9,:,:].view();  n += 9
    rcorn,hcord,phcorn,rcent,hcent,phcen = gd[0:6,:,:].view();  n += 6
    re1,he1,phe1,re2,he2,phe2,re3,he3,phe3 = gd[0:9,:,:].view();  n += 9
    gdet = gd[n]; n+=1
    if n != gd.shape[0]:
        print("rd: WARNING: nread = %d < ntot = %d: incorrect format?" % (n, gd.shape[0]) )
        return 1
    return 0

#read in a dump file
def dump_assign(gd,**kwargs):
    global t,nx,ny,nz,_dx1,_dx2,_dx3,gam,hslope,a,R0,Rin,Rout,ti,tj,tk,x1,x2,x3,r,h,ph,rho,ug,vu,B,pg,cs2,Sden,U,gdetB,divb,uu,ud,bu,bd,v1m,v1p,v2m,v2p,gdet,bsq,gdet,alpha,rhor, ktot, pg, ju, jd, jsq, jcurr
    nx = kwargs.pop("nx",nx)
    ny = kwargs.pop("ny",ny)
    nz = kwargs.pop("nz",nz)
    ti,tj,tk,x1,x2,x3,r,h,ph,rho,ug = gd[0:11,:,:].view(); n = 11
    pg = (gam-1)*ug
    lrho=np.log10(rho)
    vu=np.zeros_like(gd[0:4])
    B=np.zeros_like(gd[0:4])
    vu[1:4] = gd[n:n+3]; n+=3
    B[1:4] = gd[n:n+3]; n+=3
    #if total entropy equation is evolved (on by default)
    if DOKTOT == 1:
      ktot = gd[n]; n+=1
    divb = gd[n]; n+=1
    uu = gd[n:n+4]; n+=4
    ud = gd[n:n+4]; n+=4
    bu = gd[n:n+4]; n+=4
    bd = gd[n:n+4]; n+=4
    bsq = mdot(bu,bd)
    v1m,v1p,v2m,v2p,v3m,v3p=gd[n:n+6]; n+=6
    gdet=gd[n]; n+=1
    ju = gd[n:n+4]; n+=4
    jd = gd[n:n+4]; n+=4
    jsq = mdot(ju,jd)
    jcurr = np.sqrt(jsq)
    jdotu = mdot(ju,ud)
    Jsq = jsq + jdotu*jdotu
    gJsq = gdet*Jsq
    rhor = 1+(1-a**2)**0.5
    if "guu" in globals():
        #lapse
        alpha = (-guu[0,0])**(-0.5)
    if n != gd.shape[0]:
        print("rd: WARNING: nread = %d < ntot = %d: incorrect format?" % (n, gd.shape[0]) )
        return 1
    return 0




# Constants, definitions and units

# many are related to grmonty, but we can set a black hole mass if we want to have a "feeling" for the physical values

global MP, ME, CL, GNEWT, KBOL, SIGMA_THOMSON, MSUN, LSUN, YEAR, MBH
global TPTE_DISK, TPTE_JET, THETAE_MAX
global M_unit, L_unit, T_unit, RHO_unit, U_unit, B_unit, Ne_unit

# all constants in cgs units
ME = 9.1093826e-28 # electron mass
MP = 1.67262171e-24 # proton mass
CL = 2.99792458e10 # speed of light
GNEWT = 6.6742e-8 # gravitational constant
KBOL = 1.3806505e-16 # Boltzmann constant
SIGMA_THOMSON = 0.665245873e-24 # Thomson cross-section
MSUN = 1.989e33 # solar mass
LSUN = 3.827e33 # solar luminosity
YEAR = 31536000 # seconds in a year

# temperature and beta-prescription (Moscibrodzka 2016)
TPTE_DISK = 20. # R_high
TPTE_JET = 1. # R_low
THETAE_MAX = 1000.
TP_OVER_TE = 100.0

# grmonty units and BH mass
MBH = 4.5e6 * MSUN # Sgr A*
#MBH = 6.2e9 * MSUN # M87
#MBH = 5.0e9 * MSUN
#MBH = 10.0 * MSUN
M_unit = 1.0e19
#M_unit = 2.0*10e11
L_unit = GNEWT * MBH / (CL * CL)
T_unit = L_unit / CL
RHO_unit = M_unit / (L_unit*L_unit*L_unit)
U_unit = RHO_unit * CL * CL
B_unit = CL * np.sqrt(4. * np.pi * RHO_unit)
Ne_unit = RHO_unit / (MP + ME)

N1, N2, N3 = 1024, 512, 1
rthr = 624
dumpmin = 0
dumpmax = 1000
imin = 0
imax = N1

disc = 0
outflows = 0
jet = 1

if (disc):
    print("DISC")
elif (outflows):
    print("OUTFLOWS")
elif (jet):
    print("JET")

magenergy_arr = []
mean_energy_arr = []
max_energy=0
for dumpno in range(dumpmin,dumpmax+1):
    dumpfile = "dump%03d" %(dumpno)
    rd(dumpfile)

    sum_sheets_arr = []

    for i in range(imin,imax):
        for j in range(0,N2):
            for k in range(0,N3):
                if (disc):
                    myfilename = path_to_run + "sheets_disc/dump%03d_%d_%d_%d_s" %(dumpno,i,j,k)
                elif (outflows):
                    myfilename = path_to_run + "sheets_outflows/dump%03d_%d_%d_%d_s" %(dumpno,i,j,k)
                elif (jet):
                    myfilename = path_to_run + "sheets_jet/dump%03d_%d_%d_%d_s" %(dumpno,i,j,k)
                if (os.path.isfile(myfilename)):
                    myfile = open(myfilename, "rb")
                    dtype = np.float64
                    body = np.fromfile(myfile,dtype=dtype,count=-1)
                    sheetrows = (body.size)/10
                    mydata = body.view().reshape((int(10),int(sheetrows)), order='F')

                    myBr       = mydata[0]
                    myBth      = mydata[1]
                    myBphi     = mydata[2]
                    myJ        = mydata[3]
                    mybeta     = mydata[4]
                    mysigma    = mydata[5]
                    mysigmaphi = mydata[6]
                    myi        = mydata[7]
                    myj        = mydata[8]
                    myk        = mydata[9]
                    #myVrec     = mydata[10]
                    myfile.close()

                    sum_sheet = 0
                    for ii in range(len(myBr)):
                        if ii == 0:
                            distance = 0
                            distanceold = 0
                        else:
                            r1 = r[int(myi[ii])][int(myj[ii])][0]
                            r2 = r[int(myi[ii-1])][int(myj[ii-1])][0]
                            th1 = h[int(myi[ii])][int(myj[ii])][0]
                            th2 = h[int(myi[ii-1])][int(myj[ii-1])][0]
                            distance = np.sqrt(r1*r1 + r2*r2 - 2*r1*r2*np.cos(th1 - th2)) + distanceold
                            B2_8pi = (myBr[ii]*myBr[ii] + myBth[ii]*myBth[ii] + myBphi[ii]*myBphi[ii])/(8*np.pi)
                            mysum = B2_8pi * (distance-distanceold)
                            sum_sheet = sum_sheet + mysum
                            distanceold = distance
                    sum_sheets_arr.append(sum_sheet)#/len(myBr-1))
                    mean_energy_arr.append(sum_sheet)#/len(myBr-1)
                    if (sum_sheet > max_energy):
                        max_energy = sum_sheet

    sum_snapshot = sum(sum_sheets_arr)
    magenergy_arr.append(sum_snapshot)

#if (disc):
#    energyfile = open("energy_disc_2D.dat", "w+")
#elif (jet):
#    energyfile = open("energy_jet_2D.dat", "w+")
mean_energy = np.mean(mean_energy_arr)
std_energy = np.std(mean_energy_arr)
print("mean energy [code units]",mean_energy)
print("mean energy [erg] =",mean_energy*B_unit*B_unit*L_unit*L_unit*L_unit)
print("std error energy [code units]",std_energy)
print("std error energy [erg] =",std_energy*B_unit*B_unit*L_unit*L_unit*L_unit)
print("max energy [code units] =",max_energy)
print("max energy [erg] =",max_energy*B_unit*B_unit*L_unit*L_unit*L_unit)


magenergy_arr = [element for element in magenergy_arr if str(element) != 'nan']
magenergy_arr = [element for element in magenergy_arr if str(element) != 'inf']
magenergy_arr = np.array(magenergy_arr)

if (disc):
    magenergy_file = open(path_to_run + "magenergy_2D_disc.dat","w+")
elif (outflows):
    magenergy_file = open(path_to_run + "magenergy_2D_outflows.dat","w+")
elif (jet):
    magenergy_file = open(path_to_run + "magenergy_2D_jet.dat","w+")
for line in range(len(magenergy_arr)):
    magenergy_file.write(str(magenergy_arr[line]))
    magenergy_file.write("\n")
magenergy_file.close()
