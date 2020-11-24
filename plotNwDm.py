import matplotlib.pyplot as plt
import pickle
import numpy as np

[zKuL_SP,zKaL_SP,nodesL_SP,retr_SP]=pickle.load(open("SP_Retrievals.pklz","rb"))

hist_dm_SP=np.histogram(retr_SP[:,3],bins=0.25+0.125*np.arange(23))
hist_Nw_SP=np.histogram(retr_SP[:,4],bins=-2+0.125*np.arange(33))


[zKuL_T,zKaL_T,nodesL_T,retr_T]=pickle.load(open("Tropics_Retrievals.pklz","rb"))

hist_dm_T=np.histogram(retr_T[:,3],bins=0.25+0.125*np.arange(23))
hist_Nw_T=np.histogram(retr_T[:,4],bins=-2+0.125*np.arange(33))

plt.step(hist_dm_SP[1][:-1],100*hist_dm_SP[0]/sum(hist_dm_SP[0]))
plt.step(hist_dm_SP[1][:-1],100*hist_dm_T[0]/sum(hist_dm_T[0]))
plt.xlabel('Dm (mm)')
plt.ylabel('Percentage')
plt.legend(['Southern Pacific','Tropics'])
plt.title('Normalized PDF')

plt.figure()
plt.step(hist_Nw_SP[1][:-1],100*hist_Nw_SP[0]/sum(hist_Nw_SP[0]))
plt.step(hist_Nw_SP[1][:-1],100*hist_Nw_T[0]/sum(hist_Nw_T[0]))
plt.xlabel('log10(Nw/Nw_ref)')
plt.ylabel('Percentage')
plt.title('Normalized PDF')
plt.legend(['Southern Pacific','Tropics'])

zbbL=[]
z1L=[]
z2L=[]
piaLs=[]
# top, zero degree, clutter free, brightband top, bright band peak and surface bins
# rrate_cmb,dprsfcRate,pType,dm1d,dn1d,dmDPR,piaku,piaka,piaSRT,epsDPR,reliabFlag,binBB
#         0,    1     ,  2  , 3  , 4  ,  5,    6,    7,    8,    9,       10,      11    
for zku1,zka1,nodes,retr1 in zip(zKuL_SP,zKaL_SP,nodesL_SP,retr_SP):
    zbb=zku1[nodes[-2]]-zka1[nodes[-2]]
    if zku1[nodes[-2]]>10 and zka1[nodes[-2]]>10:
        zbbL.append([zbb,zku1[nodes[-2]]])
        if nodes[2]-nodes[1]>8:
            nodes2=nodes[1]+8
            z1L.append(zku1[nodes2-60:nodes2])
            z2L.append(zka1[nodes2-60:nodes2])
            piaLs.append([retr1[6],retr1[7], retr1[8]])

from sklearn.cluster import KMeans
z1L=np.array(z1L)
z2L=np.array(z2L)
piaLs=np.array(piaLs)
z1L[z1L<0]=0
z2L[z2L<0]=0
kmeans = KMeans(n_clusters=16, random_state=0).fit(z1L[:,:])
zS=[]
for i in range(16):
    a=np.nonzero(kmeans.labels_==i)
    zS.append(np.mean((10.0**(0.1*z1L[a[0],-8:]*0.71)).sum(axis=1)))

ind=np.argsort(zS)
for i in range(16):
    if i%4==0:
        plt.figure(figsize=(6,6))
    plt.subplot(2,2,i%4+1)
    a=np.nonzero(kmeans.labels_==ind[i])
    plt.plot(z1L[a[0],:].mean(axis=0)[::-1],np.arange(60)*0.125-1.0)
    plt.plot(z2L[a[0],:].mean(axis=0)[::-1],np.arange(60)*0.125-1.0)
    plt.title("PIA(ku)=(%6.2f,%6.2f)"%(piaLs[a[0],0].mean(axis=0),piaLs[a[0],-1].mean(axis=0)))
    plt.tight_layout()

import combAlg as cmb
cmb.mainfortpy()
cmb.initp2()

dmNw=np.array([[0.30136546,	2.7706423],[0.34136546,	2.7706423],[0.4437751,	2.6972477],[0.49497992,	2.6972477],[0.5803213,	2.5504587],[0.6997992,	2.3302753],[0.7851406,	2.0856268],[0.9558233,	1.7920489],[1.0753012,	1.6697248],[1.2630522,	1.5229357],[1.4166666,	1.351682],[1.6385542,	1.2293578],[1.809237,	1.1070336],[2.0993977,	0.9357798],[2.3554218,	0.8379205],[2.5943775,	0.7155963],[2.816265,	0.6666667],[2.8333333,	0.6911315],[3.0722892,	0.56880736],[3.498996,	0.5198777]])
dmNw=np.array([[0.106177606,	6.300601],[0.22007722,	5.819639],[0.34555984,	5.3547096],[0.45945945,	4.9378757],[0.6216216,	4.4248495],[0.7895753,	4.0240483],[0.95173746,	3.7194388],[1.1081082,	3.5911825],[1.2857143,	3.4949899],[1.438224,	3.3667336],[1.6969112,	3.1583166],[1.8880309,	3.0781562],[2.0752895,	2.9819639],[2.2818532,	2.9178357],[2.4401546,	2.8537073],[2.5019305,	2.8216434], [3.5,2.6]])

NwInt=np.interp(cmb.tablep2.dmj[:289],dmNw[:,0],dmNw[:,1])-np.log10(8000)
plt.plot(cmb.tablep2.dmj[:289],cmb.tablep2.zkusj[:289]+10*NwInt)
