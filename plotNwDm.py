import matplotlib.pyplot as plt
import pickle
import numpy as np

[zKuL_SP,zKaL_SP,nodesL_SP,retr_SP]=pickle.load(open("SP_Retrievals.pklz","rb"))

hist_dm_SP=np.histogram(retr_SP[:,3],bins=0.25+0.125*np.arange(23))
hist_Nw_SP=np.histogram(retr_SP[:,4],bins=-2+0.125*np.arange(33))


[zKuL_T,zKaL_T,nodesL_T,retr_T]=pickle.load(open("Tropics_Retrievals.pklz","rb"))

hist_dm_T=np.histogram(retr_T[:,3],bins=0.25+0.125*np.arange(23))
hist_Nw_T=np.histogram(retr_T[:,4],bins=-2+0.125*np.arange(33))

#plt.step(hist_dm_SP[1][:-1],100*hist_dm_SP[0]/sum(hist_dm_SP[0]))
#plt.step(hist_dm_SP[1][:-1],100*hist_dm_T[0]/sum(hist_dm_T[0]))
#plt.xlabel('Dm (mm)')
#plt.ylabel('Percentage')
#plt.legend(['Southern Pacific','Tropics'])
#plt.title('Normalized PDF')

#plt.figure()
#plt.step(hist_Nw_SP[1][:-1],100*hist_Nw_SP[0]/sum(hist_Nw_SP[0]))
#plt.step(hist_Nw_SP[1][:-1],100*hist_Nw_T[0]/sum(hist_Nw_T[0]))
#plt.xlabel('log10(Nw/Nw_ref)')
#plt.ylabel('Percentage')
#plt.title('Normalized PDF')
#plt.legend(['Southern Pacific','Tropics'])

zbbL=[]
z1L=[]
z2L=[]
piaLs=[]
zku_subL=[]
zka_subL=[]
pia_subL=[]
nodes_subL=[]
# top, zero degree, clutter free, brightband top, bright band peak and surface bins
# rrate_cmb,dprsfcRate,pType,dm1d,dn1d,dmDPR,piaku,piaka,piaSRT,epsDPR,reliabFlag,binBB
#         0,    1     ,  2  , 3  , 4  ,  5,    6,    7,    8,    9,       10,      11
retrsL=[]
for zku1,zka1,nodes,retr1 in zip(zKuL_SP,zKaL_SP,nodesL_SP,retr_SP):
    zbb=zku1[nodes[-2]]-zka1[nodes[-2]]
    #if zku1[nodes[-2]]>-1000 and zka1[nodes[-2]]>-1000:
    #    zbbL.append([zbb,zku1[nodes[-2]]])
    if nodes[2]-nodes[1]>6:
        nodes2=nodes[1]+6
        z1L.append(zku1[nodes2-60:nodes2])
        z2L.append(zka1[nodes2-60:nodes2])
        zku_subL.append(zku1)
        zka_subL.append(zku1)
        nodes_subL.append(nodes)
        pia_subL.append(retr1[8])
        piaLs.append([retr1[6],retr1[7], retr1[8]])
        retrsL.append(retr1)

from plotNwDm_inc import cluster

import combAlg as cmb
cmb.mainfortpy()
cmb.initp2()
NwInt,ind,kmeans=cluster(z1L,z2L,piaLs,cmb)
zkuNwDm=cmb.tablep2.zkusj[:289]+10*NwInt
cmb.tablep2.zkudn[0:289]=zkuNwDm
z1L=np.array(z1L)
dmLs=[[]for i in range(16)]
dmLt=[]
dnCoeff=np.array([-0.01257341, -0.00933038])
pType=1
dr=0.125
eps=1
imu=3
dmLt2=[]
piaLs2=[]
for i in ind:
    a=np.nonzero(kmeans.labels_==i)
    for j1 in a[0]:
        i1 = cmb.bisection2(zkuNwDm,z1L[j1,-1])
        dmLs[i].append(cmb.tablep2.dmj[i1-1])
        dmLt.append(cmb.tablep2.dmj[i1-1])
        # top, zero degree, clutter free, brightband top, bright band peak and surface bins
        # rrate_cmb,dprsfcRate,pType,dm1d,dn1d,dmDPR,piaku,piaka,piaSRT,epsDPR,reliabFlag,binBB
        #         0,    1     ,  2  , 3  , 4  ,  5,    6,    7,    8,    9,       10,      11
        nodes1=nodes_subL[j1]
        bst,bzd,bcf,binBBT,binBB,bsfc=nodes1
        zKu1=zku_subL[j1]
        zKa1=zka_subL[j1]
        reliabFlag=1
        dn=0.0
        piaSRT=pia_subL[j1]
        pType=1
        if retrsL[j1][-2]==1:
            dn1d,dm1d,rrate1d,zkuc,zkasim,\
                epst,piaku,piaka = cmb.prof1d(bst,bzd,bcf,bsfc,binBB,binBBT,zKu1,\
                                              zKa1,pType,dr,eps,imu,dnCoeff,dn,\
                                              piaSRT*5,reliabFlag,piaSRT,-1)
            dmLt2.append(dm1d[bcf])
            print(retrsL[j1][2],piaku,piaSRT)
            
            piaLs2.append([piaku,piaSRT,rrate1d[bcf],retrsL[j1][1]])
