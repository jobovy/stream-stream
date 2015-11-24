# stream2_util.py: some utility functions
import os, os.path
import csv
import numpy
import copy
from galpy.potential import LogarithmicHaloPotential
from galpy.util import bovy_coords, multi
from galpy.actionAngle import actionAngleIsochroneApprox
R0= 8.
V0= 220.
lp= LogarithmicHaloPotential(normalize=1.,q=0.9)
def rectangular_to_cylindrical(xv):
    R,phi,Z= bovy_coords.rect_to_cyl(xv[:,0],xv[:,1],xv[:,2])
    vR,vT,vZ= bovy_coords.rect_to_cyl_vec(xv[:,3],xv[:,4],xv[:,5],
                                          R,phi,Z,cyl=True)
    out= numpy.empty_like(xv)
    # Preferred galpy arrangement of cylindrical coordinates
    out[:,0]= R
    out[:,1]= vR
    out[:,2]= vT
    out[:,3]= Z
    out[:,4]= vZ
    out[:,5]= phi
    return out

def calc_aA_sim(RvR,filename,snap_gc):
    # Calculate the action angle variables for a simulation and store
    if not os.path.exists(filename):
        aAI= actionAngleIsochroneApprox(pot=lp,b=0.8)
        nbatch= 20
        multiOut= multi.parallel_map(\
            lambda x: aAI.actionsFreqsAngles(RvR[x*nbatch:(x+1)*nbatch,0]/R0,
                                             RvR[x*nbatch:(x+1)*nbatch,1]/V0,
                                             RvR[x*nbatch:(x+1)*nbatch,2]/V0,
                                             RvR[x*nbatch:(x+1)*nbatch,3]/R0,
                                             RvR[x*nbatch:(x+1)*nbatch,4]/V0,
                                             RvR[x*nbatch:(x+1)*nbatch,5]),
            range(len(snap_gc)//nbatch),
            numcores=25)
        acfs= numpy.reshape(numpy.swapaxes(numpy.array(multiOut),0,1),
                            (9,numpy.prod(numpy.array(multiOut).shape)//9))
        # Write to file
        csvfile= open(filename,'w')
        writer= csv.writer(csvfile,delimiter=',')
        for jj in range(len(acfs[0])):
            writer.writerow([acfs[0][jj],acfs[1][jj],acfs[2][jj],
                                     acfs[3][jj],acfs[4][jj],acfs[5][jj],
                                     acfs[6][jj],acfs[7][jj],acfs[8][jj]])
            csvfile.flush()
        csvfile.close()
    else:
        acfs= numpy.loadtxt(filename,delimiter=',').T
    return acfs

def calc_apar(acfs,angle=None,freq=False,debrisThreshold=6.):
    # Calculate the parallel angle offset, 
    # of angle if set (otherwise of the entire simulation), 
    # angle is a frequency if freq
    thetar= acfs[6]
    thetap= acfs[7]
    thetaz= acfs[8]
    if not angle is None:
        if not freq:
            angle[0]= (numpy.pi+(angle[0]-numpy.median(thetar))) % (2.*numpy.pi)
            angle[1]= (numpy.pi+(angle[1]-numpy.median(thetap))) % (2.*numpy.pi)
            angle[2]= (numpy.pi+(angle[2]-numpy.median(thetaz))) % (2.*numpy.pi)
    thetap= (numpy.pi+(thetap-numpy.median(thetap))) % (2.*numpy.pi)
    debrisIndx= numpy.fabs(thetap-numpy.pi) > (debrisThreshold*numpy.median(numpy.fabs(thetap-numpy.median(thetap))))
    if angle is None:
        thetar= (numpy.pi+(thetar-numpy.median(thetar))) % (2.*numpy.pi)
        thetaz= (numpy.pi+(thetaz-numpy.median(thetaz))) % (2.*numpy.pi)
        #center around 0 (instead of pi)
        thetar-= numpy.pi
        thetap-= numpy.pi
        thetaz-= numpy.pi
    elif freq:
        thetar= angle[0]
        thetap= angle[1]
        thetaz= angle[2]
    else:
        thetar= angle[0]-numpy.pi
        thetap= angle[1]-numpy.pi
        thetaz= angle[2]-numpy.pi
    #Frequencies
    Or= acfs[3]
    Op= acfs[4]
    Oz= acfs[5]
    dOr= Or[debrisIndx]-numpy.median(Or)
    dOp= Op[debrisIndx]-numpy.median(Op)
    dOz= Oz[debrisIndx]-numpy.median(Oz)
    dO= numpy.vstack((dOr,dOp,dOz))
    dO4dir= copy.copy(dO)
    dO4dir[:,dO4dir[0] < 0.]*= -1.
    dOdir= numpy.median(dO4dir,axis=1)
    dOdir/= numpy.sqrt(numpy.sum(dOdir**2.))
    #parallel angle
    dangle= numpy.vstack((thetar,thetap,thetaz))
    return numpy.dot(dangle.T,dOdir)
