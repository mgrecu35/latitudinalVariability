fname='2A.GPM.DPR.V8-20180723.20180118-S125407-E142641.022105.V06A.HDF5'
fname='2A.GPM.DPR.V8-20180723.20180215-S134018-E151250.022541.V06A.HDF5'

fname='2A-CS-Tasmania.GPM.DPR.V8-20180723.20200619-S194636-E194813.035844.V06A.HDF5'
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import combAlg as cmb
#from julia.api import Julia
#jl=Julia()#compiled_modules=False)
from numpy import *
cmb.mainfortpy()
cmb.initp2()
#if(dmm<0.9) epst=epst*(0.9/dmm)**0.25
#if(dmm>1.0) epst=epst*(1.0/dmm)**0.25     
dnw=cmb.tablep2.dmj*0.
for i in range(289):
    if cmb.tablep2.dmj[i]<0.8:
        dmm=cmb.tablep2.dmj[i]
        dnw[i]=1*log((0.8/dmm)**0.25)
    if cmb.tablep2.dmj[i]>0.8:
        dmm=cmb.tablep2.dmj[i]
        dnw[i]=0.5*log((0.8/dmm)**0.25)

dnw=dnw*2.0
zku1=cmb.tablep2.zkusj[:289]+10*dnw[:289]
att1=cmb.tablep2.attkuj[:289]*10**dnw[:289]
rJ=cmb.tablep2.rj[:289]*10**dnw[:289]

attKuCoeff=array([ 0.06972763, -3.16804927])
attKuCoeff=polyfit(zku1[50:],log10(att1[50:]),1)
rrateCoeff=polyfit(zku1[50:],log10(rJ[50:]),1)
#stop
import glob
fst=[]
#for i in range(30):
#    fs=glob.glob("/gpmdata/2018/08/%2.2i/radar/2A.GPM.DPR.V8*HDF5"%(i+1))
#    fst.extend(sorted(fs))
fst=glob.glob("../monthly/SEAsia/2A*DPR*HDF5")
fs=glob.glob("../monthly/SEAsiaCS/2A*DPR*HDF5")
fst=sorted(fst)
fst.extend(sorted(fs))
fstKu=glob.glob("../ORO/DPR-CS/2A-CS-CON*DPR*HDF5")
#fstKu=glob.glob("../ORO/DPR-CS/2A-CS-CON*Ku*HDF5")
print(len(fst))
#stop
icount=0
iplot=1
zKuL=[]
zKaL=[]
hL=[]
r1L=[]
ibbL=[]
rd1L=[]
#for f in sorted(fs)[0:1]:
iplot=0
ifile=0
dFL=[]
DPR_RDmCoeff=[ 0.01608098, -0.82884294]
fst=fst[:]
import pickle
#(175-binSf[i1,j1])*0.125, (binSf[i1,j1]-binClut[i1,j1])*0.125,\
#(175-zeroDeg[i1,j1])*0.125,binSf[i1,j1],zeroDeg[i1,j1],binClut[i1,j1],\
#piaSRT[i1,j1],dm[i1,j1,binClut[i1,j1]],# sfcRain[i1,j1],stormT[i1,j1],\
#sfcType[i1,j1],relFlag[i1,j1], piaSRTKu[i1,j1],relFlagKu[i1,j1]]

z1L,zKuLt,zKaLt,dmL,nwL,rateL,addInfoL=pickle.load(open("../monthly/cvProfs2018_20.pklz","rb"))
addInfoL=array(addInfoL)
zmL,kgainL,xL,kmeans,ic=pickle.load(open("../monthly/kFilterCvClasses20.pklz","rb"))

z1L=array(z1L)
z1L[z1L<0]=0
labels=kmeans.predict(z1L[:,:-1])

#btop=175-int(addInfoL[:,9]/125)
#sfcRain=addInfop[:,8]
#bsfc,bzd,bcf=addInfoL[:,3:6]

piaLs=[[]for i in range(16)]
NwDm=[]
rLs=[[]for i in range(16)]
dnL=[[] for i in range(16)]
import numpy as np

dmNw=np.array([[0.106177606,6.300601],[0.22007722,5.819639],[0.34555984,5.3547096],[0.45945945,4.9378757],[0.6216216,4.4248495],[0.7895753,4.0240483],[0.95173746,3.7194388],[1.1081082,3.5911825],[1.2857143,3.4949899],[1.438224,3.3667336],[1.6969112,3.1583166],[1.8880309,3.0781562],[2.0752895,2.9819639],[2.2818532,2.9178357],[2.4401546,2.8537073],[2.5019305,2.8216434], [3.5,2.6]])

NwInt=0.5*(np.interp(cmb.tablep2.dmj[:289],dmNw[:,0],dmNw[:,1])-np.log10(8000))

zkuNwDm=cmb.tablep2.zkusj[:289]+10*NwInt
#
#plt.plot(cmb.tablep2.dmj[:289],cmb.tablep2.zkusj[:289]+10*NwInt)
#stop
piaLs=[]
zetaHL=[]
epsL=[]
rsfcL=[]
zdnw=[]
for i in range(z1L.shape[0]):
    dnCoeff=array([ 0.00827692, -0.19116723])
    dnCoeff=array([ 0.01608098, -0.82884294])
    dnCoeff=array([-0.016238408, -0.01593891])
    zKu1=zKuLt[i][1]
    zKa1=zKaLt[i][1]
    bsfc,bzd,bcf=addInfoL[i,3:6]
    bsfc=int(bsfc)
    bzd=int(bzd)
    bcf=int(bcf)
    zmax=zKu1[50:min(bzd,bcf-4)].max()
    if zmax<35:
        dn=0.3
    else:
        if zmax<40:
            dn=0.3+0.2*(zmax-35)/5.
        else:
            dn=0.5
    dn=0.1
    dn=0.1+(ic[labels[i]]-8)/8*0.2
    dn=0.0
    reliabFlag=0
    dsrtPIA=-99.9
    pType=2
    dr=0.125
    eps=0.95
    bst=175-int(addInfoL[i,9]/125)
    imu=3
    binBB=0
    binBBT=0
    eps=1.0
    if (addInfoL[i,-1]==1 or addInfoL[i,-2]==1) and (labels[i] in ic[0:]):
        dn1d,dm1d,rrate1d0,zkuc,zkasim,\
            epst,piaku0,piaka = cmb.prof1d(bst,bzd,bcf,bsfc,\
                                           binBB,binBBT,zKu1,zKa1,pType,dr,eps,imu,dnCoeff,dn,\
                                           dsrtPIA,reliabFlag,addInfoL[i,-2],addInfoL[i,-1])
        
        
        #stop
        piaLs.append([piaku0,addInfoL[i,-2]])
        rsfcL.append([rrate1d0[bcf],addInfoL[i,8]])
        #stop


if 1==2:
    zeta=10**(attKuCoeff[0]*zKu1[bzd:bcf]+attKuCoeff[1])
    q=0.2*log(10)
    beta=0.69
    dr=0.125
    piaHB=-10/beta*log10(1-q*beta*zeta.sum()*dr)
    piaSRT=addInfoL[i,-2]
    eps=(1-10**(-0.1*piaSRT*beta))/(q*beta*zeta.sum()*dr)
    if eps>5:
        eps=5
    if eps<0.01:
        eps=0.01
    epsL.append(eps)
    epsD=10**(0.0+random.randn(100)*0.35)
    piaHBw=0.0
    rsfcw=0.0
    wprob=0.0
    zsfcm=0
    for eps1 in epsD:
        dn=0.999*eps1**(1/(1-beta))
        if eps1*q*beta*zeta.sum()*dr<0.985:
            piaHB=-10/beta*log10(1-eps1*q*beta*zeta.sum()*dr)
            zsfc=zKu1[bcf]+piaHB
            rsfc=10**(rrateCoeff[0]*zsfc+rrateCoeff[1])*dn**(1-0.616)
            prob1=exp(-0.5*(piaHB-addInfoL[i,-2])**2/4)
            rsfcw+=prob1*rsfc
            piaHBw+=prob1*piaHB
            zsfcm+=zsfc*prob1
            wprob+=prob1
    zsfcm/=wprob
    rsfcw/=wprob
    rsfcE=10**(rrateCoeff[0]*zsfcm+rrateCoeff[1])
    dnEst=1/(1-0.616)*log10(rsfcw/rsfcE)
    zdnw.append([zsfcm,dnEst])
    zetaHL.append(q*beta*zeta.sum()*0.125)
