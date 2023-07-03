import math
import random
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate

def Integ1(maxt,maxd,maxdif,ddf,gaussd,PMn):
    Gsn = np.zeros((maxd,maxt))
    for it in range(1,maxt):
        for ir in range(1,maxd):
         
            nsum1=0.0
            for idf in range(1,maxdif):
                if idf == 1 or idf == maxdif:
                    fc=0.50
                else :
                    fc =1.0
                nsum1 += fc*gaussd[ir][idf]*PMn[idf][it]*ddf
            Gsn[ir][it]= nsum1             
    return (Gsn)

def Integ2(maxt,maxd,maxdif,ddr,PM0,gaussd,Gs,Gsn):
     nsum1=np.zeros((maxdif,maxt)) 
     PMn=np.zeros((maxdif,maxt)) 
     tol = 0.00001
     for it in range(1,maxt):
        for idf in range(1,maxdif):
            nsum1[idf][it]=0.0
            for ir in range(1,maxd):
                if ir == 1 or ir == maxd:
                    fc=0.50
                else :
                    fc =1.0
                if Gsn[ir][it]> tol : 
                    gratio=Gs[ir][it]/Gsn[ir][it]
                else :
                    gratio = 1.0
                    
                nsum1[idf][it] += fc*gratio*gaussd[ir][idf]*ddr*(2.0*math.pi*ir*ddr)

     for it in range(1,maxt):
         for idf in range(1,maxdif):
             PMn[idf][it] = nsum1[idf][it]*PM0[idf][it]
            # print("PM",PMn[idf][it], nsum1[idf][it],PM0[idf][it])
     return (PMn)

def normalizeP(maxt,maxdif,PMn):
    nsum1 = np.zeros(maxt)
    for it in range(maxt):
        nsum1[it]=1.0
        for idf in range(1,maxdif):
            nsum1[it] += PMn[idf][it]

    for it in range(maxt):                  
        for idf in range(1,maxdif):
            PMn[idf][it]=PMn[idf][it]/nsum1[it]
    return (PMn) 
def normalizeG(maxt,maxd,Gsn):
    nsum1 = np.zeros(maxt)
    for it in range(maxt):
        nsum1[it]=1.0
        for ir in range(1,maxd):
            nsum1[it] += Gsn[ir][it] 
    for it in range(maxt):
        for ir in range(1,maxd):
            Gsn[ir][it] = Gsn[ir][it]/nsum1[it]
    return (Gsn)


def diffP(maxt,maxdif,PMn,PM0):
        deltad=0.0
        for idf in range(1,maxdif):
            for it in range(maxt):
                deltad += (PMn[idf][it]-PM0[idf][it])**2
        return(deltad)


if __name__ == '__main__':

    import numpy as np

    timecalc = 9
    pindT, timeT, xT, yT = np.genfromtxt("./20170202c059_RPE1_CAG_H2B_Halotag_TMR80pM_nonlam_non_starve_ctrl_info.txt", unpack=True)
    
    np.array(pindT)
    for x in pindT:
        pindT = pindT.astype(int)
        timeT = timeT.astype(int)
    lastnucIndex=pindT[-1]
    firstnucIndex=pindT[0]
    nucx = [[ip for ip in range(1,lastnucIndex)] for itt in  range(len(timeT))]
    nucy = [[ip for ip in range(1,lastnucIndex)] for itt in  range(len(timeT))]
    nuct = [[ip for ip in range(1,lastnucIndex)] for itt in  range(len(timeT))]

   
    for ix in range(0,len(pindT)):
        for ind in range(firstnucIndex,lastnucIndex):
            if(ind==pindT[ix]):
                nucx[ind][timeT[ix]]=xT[ix]
                nucy[ind][timeT[ix]]=yT[ix]
                nuct[ind][timeT[ix]]=timeT[ix]
              #  print(pindT[ix],nuct[ind][timeT[ix]])
                
    ji=0
    ix=0
    intime=np.empty(lastnucIndex,dtype=int)
    outtime=np.empty(lastnucIndex,dtype=int)
    while pindT[ix] < lastnucIndex :
          if pindT[ix] !=ji:
             intime[pindT[ix]]=timeT[ix]
             ji = pindT[ix]
          else :    
             outtime[pindT[ix]]=timeT[ix]
             ix=ix+1
                  
    delta = np.zeros_like(intime)            
    for ix in range(1,lastnucIndex):
        delta[ix]=(outtime[ix]-intime[ix]).astype(int)

    maxt=np.amax(delta)
             
    # Calculate MSD
    r2=np.zeros(maxt)
    nrm2=np.ones(maxt)
    
    for ind in range(firstnucIndex,lastnucIndex-1):
       for k1 in range(intime[ind],outtime[ind]-1):
           for k2 in range(k1,outtime[ind]):
               ddx = nucx[ind][k1]-nucx[ind][k2]
               ddy = nucy[ind][k1]-nucy[ind][k2]
               r2[k2-k1] += ddx*ddx+ddy*ddy 
               nrm2[k2-k1]+=  1.0

    for i in range(maxt): 
           r2[i]=r2[i]/nrm2[i]
    # plot MSD       
    #plt.plot(r2) 
    #plt.show()


   # Calculate Self part van Hove Correlation
   # 
    rmax=10.0
    maxd=500
    ddr=rmax/maxd
    Gs=np.zeros((maxd,maxt))
    print("shape",np.shape(Gs))
    np.shape(nrm2)
    for ind in range(firstnucIndex,lastnucIndex-1):
       for k1 in range(intime[ind],outtime[ind]-1):
           for k2 in range(k1,outtime[ind]):
               ddx = nucx[ind][k1]-nucx[ind][k2]
               ddy = nucy[ind][k1]-nucy[ind][k2]
               rd= math.sqrt(ddx*ddx+ddy*ddy) 
               k0 = int(rd/ddr)
               if k0 < maxd :
                  Gs[k0][k2-k1] +=  1.0
    for k0 in range(1,maxd):
       for k1 in range(maxt):
           Gs[k0][k2-k1] = Gs[k0][k2-k1]/(2.0*i*ddr*math.pi)

    Gs[0][:]=0.0

            
     
    # Normalize Gs
    nrm7=np.ones((maxt))
    for i in range(maxt):
        for j in range(maxd):
            nrm7[i] += Gs[j][i]  
    for i in range(maxt):
        for j in range(maxd):
            Gs[j][i]=Gs[j][i]/nrm7[i]
#    for i in range(maxd):
#        xplt[i]= i*ddr  
#        yplt[i]=Gs[i][timecalc]



    # Gaussian 
    maxdif=500
    difmax=25.0
    ddf=difmax/maxdif
    gaussd=np.zeros((maxd,maxdif))
    #  Note array starts from 1 not 0
    for i in range(1,maxdif):
        ms = i*ddf
        for j in range(1,maxd):
            r2 = j*j*ddr*ddr
            gaussd[i][j]=(1.0/(ms*math.pi))*math.exp(-r2/ms)
   
    

    # trial P(M) 
    a1=1.15
    a=1.0
    a0=0.59
    PM=np.zeros((maxdif,maxt))
    for it in range(1,maxt):
        for im in range(maxdif):
            d = im*ddf
            arg=a*(d-a0)*(d-a0)
            PM[im][it] = a1*math.exp(-arg)
    xplt= np.zeros(maxd)         
    yplt= np.zeros(maxd)        
           
   # Richardson-Lucy loop
    

    deltad = 1.0
    tol = 0.0000001
    Gsn=np.zeros((maxd,maxt))
    PMn=np.zeros((maxdif,maxt))
    PM0=np.zeros((maxdif,maxt))
    PMn=PM
    Gsn = Gs
    while deltad > tol:
        PM0 = PMn
        Gsn = Integ1(maxt,maxd,maxdif,ddf,gaussd,PMn)
        Gsn = normalizeG(maxt,maxd,Gsn)
        PMn = Integ2(maxt,maxd,maxdif,ddr,PM0,gaussd,Gs,Gsn)
      
        PMn = normalizeP(maxt,maxdif,PMn)
        deltad = diffP(maxt,maxdif,PMn,PM0) 
        print("deltad",deltad)

    for i in range(maxdif):
           xplt[i]= i*ddf  
           yplt[i]=PMn[i][20]

    plt.plot(xplt,yplt)
    plt.show()      
        


           
          
         
  








