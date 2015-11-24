# calc_deltav.py: calculate kicks on the GC stream in a variety of ways
import sys
import os, os.path
import pickle
import numpy
from galpy.orbit import Orbit
from galpy.util import save_pickles, multi, bovy_conversion, bovy_coords
from galpy.actionAngle import actionAngleIsochroneApprox
from galpy.snapshot import nemo_util
from galpy.potential import MovingObjectPotential
from stream2_util import rectangular_to_cylindrical, R0, V0, lp, \
    calc_apar, calc_aA_sim
def vxta_deltav(snap_gc,apar,acfs):
    """
    NAME:
       vxta_deltav
    PURPOSE:
       construct the (v,x,t,apar) for the points at which to compute deltav for the GC stream
    INPUT:
       snap_gc - the snapshot of the GC stream
       apar - parallel angle along the stream
       acfs - the action-angle coordinates for the GC stream
    OUTPUT:
       (v,x,t,apar)
    HISTORY:
       2015-11-25 - Written - Bovy (UofT)
    """
    # First we get the action-angle coordinates for the cluster orbit 
    # (approximating the cluster stream as an orbit)
    close_to_impact_indx= numpy.fabs(apar+2.4) < 0.2
    stream_orbit_RvR= rectangular_to_cylindrical(\
        numpy.median(snap_gc[close_to_impact_indx,1:,-1],axis=0)[:,numpy.newaxis].T)[0,:]
    stream_orbit= Orbit([stream_orbit_RvR[0]/R0,stream_orbit_RvR[1]/V0,stream_orbit_RvR[2]/V0,
                         stream_orbit_RvR[3]/R0,stream_orbit_RvR[4]/V0,stream_orbit_RvR[5]])
    aAI= actionAngleIsochroneApprox(pot=lp,b=0.8)
    stream_orbit_aA= aAI.actionsFreqsAngles(stream_orbit())
    impact_apar= calc_apar(acfs,list(stream_orbit_aA[6:]))
    impact_Opar= calc_apar(acfs,list(stream_orbit_aA[3:6]),freq=True)
    # Now we integrate this cluster orbit back and forth to build position and velocity of the track
    times= numpy.linspace(0.,0.12/bovy_conversion.time_in_Gyr(V0,R0),101)
    stream_orbit.integrate(times,lp)
    stream_orbit_back= stream_orbit.flip()
    stream_orbit_back.integrate(times,lp)
    x_gc= numpy.zeros((len(times)*2-1,3))
    v_gc= numpy.zeros((len(times)*2-1,3))
    t_gc= numpy.zeros((len(times)*2-1))
    apar_gc= numpy.zeros((len(times)*2-1))
    x_gc[len(times)-1:,0]= stream_orbit.x(times)
    x_gc[len(times)-1:,1]= stream_orbit.y(times)
    x_gc[len(times)-1:,2]= stream_orbit.z(times)
    x_gc[:len(times)-1,0]= stream_orbit_back.x(times[1:][::-1])
    x_gc[:len(times)-1,1]= stream_orbit_back.y(times[1:][::-1])
    x_gc[:len(times)-1,2]= stream_orbit_back.z(times[1:][::-1])
    v_gc[len(times)-1:,0]= stream_orbit.vx(times)
    v_gc[len(times)-1:,1]= stream_orbit.vy(times)
    v_gc[len(times)-1:,2]= stream_orbit.vz(times)
    v_gc[:len(times)-1,0]= -stream_orbit_back.vx(times[1:][::-1])
    v_gc[:len(times)-1,1]= -stream_orbit_back.vy(times[1:][::-1])
    v_gc[:len(times)-1,2]= -stream_orbit_back.vz(times[1:][::-1])
    t_gc[len(times)-1:]= times
    t_gc[:len(times)-1]= -times[1:][::-1]
    apar_gc= impact_apar+impact_Opar*t_gc
    return (v_gc,x_gc,t_gc,apar_gc)

def impulse_deltav_plummerint(v,x,galpot,GM,x0,v0,
                              rs=0.1,
                              tmax=0.25*0.9777922212082034/\
                                  bovy_conversion.time_in_Gyr(V0,R0)):
    """
    NAME:
       impulse_deltav_plummerint
    PURPOSE:
       calculate the velocity kick using direct plummer integration, potentially for multiple particles
    INPUT:
       v - [nstar,3] velocity of the star to be kicked at the central time
       x - [nstar,3] position of the star to be kicked at the central time
       GM - mass of the Plummer sphere doing the kicking
       x0 - [nplum,3] position of the Plummer sphere doing the kicking *at the start*
       v0 - [nplum,3] velocity of the Plummer sphere doing the kicking *at the start*
       rs - scale radius of the Plummer sphere(s)
       tmax - time to integrate the interaction over
    OUTPUT:
       deltav [nstar,3]
    HISTORY:
       2015-11-23 - Written - Bovy (UofT)
    """
    # For each dm particle, setup Plummer orbit
    plumpot= []
    times= numpy.linspace(0.,tmax,1001)
    halftimes= numpy.linspace(0.,tmax/2.,1001)
    for ii in range(len(x0)):
        tx0, tv0= x0[ii], v0[ii]
        R, phi, z= bovy_coords.rect_to_cyl(tx0[0],tx0[1],tx0[2])
        vR, vT, vz= bovy_coords.rect_to_cyl_vec(tv0[0],tv0[1],tv0[2],R,phi,z,cyl=True)
        oplum= Orbit(vxvv=[R,vR,vT,z,vz,phi])
        oplum.integrate(times,galpot,method='symplec4_c')
        plumpot.append(MovingObjectPotential(orbit=oplum,GM=GM,
                                             softening_model='plummer',softening_length=rs))
    plumpot.append(galpot) # Need to add this to!
    # Now integrate each (v,x) first backwards in galpot, then forwards in galpot+plumpot, then backwards in galpot
    deltav = numpy.zeros((len(v),3))
    R, phi, z= bovy_coords.rect_to_cyl(x[:,0],x[:,1],x[:,2])
    vR, vT, vz= bovy_coords.rect_to_cyl_vec(v[:,0],v[:,1],v[:,2],
                                            R,phi,z,cyl=True)
    for ii in range(len(v)):
        ostar= Orbit(vxvv=[R[ii],-vR[ii],-vT[ii],z[ii],-vz[ii],phi[ii]])
        ostar.integrate(halftimes,galpot,method='leapfrog_c')
        oboth= ostar(halftimes[-1]).flip()
        oboth.integrate(times,plumpot,method='leapfrog_c')
        ogalpot = oboth(times[-1]).flip()
        ogalpot.integrate(halftimes,galpot,method='leapfrog_c')
        deltav[ii,0] = -ogalpot.vx(halftimes[-1])-v[ii,0]
        deltav[ii,1] = -ogalpot.vy(halftimes[-1])-v[ii,1]
        deltav[ii,2] = -ogalpot.vz(halftimes[-1])-v[ii,2]
    return deltav

def calc_fullplummer_deltav(time,savefilename,ndm,v_gc,x_gc,rs=0.2):
    """
    NAME:
       calc_fullplummer_deltav
    PURPOSE:
       Calculate the velocity kicks using full Plummer integration using a set of DM particles; each DM particle is considered separately
    INPUT:
       time - (string) time string identifying the DM simulation (0.125, 0.25, 0.375, or 0.50)
       savefilename - name of the file to save the deltav to
       ndm - number of DM particles to compute the kick for
    OUTPUT:
       deltav [nstar,3,ndm]
    HISTORY:
       2015-11-24 - Written - Bovy (UofT)
    """
    if not os.path.exists(savefilename):
        # Load DM snapshot
        filename= os.path.join(os.getenv('DATADIR'),'bovy','stream-stream',
                               'dm_evol_%s_untilimpact.dat' % time)
        snap_dm= nemo_util.read(filename,swapyz=True)
        numpy.random.seed(1)
        permIndx= numpy.random.permutation(len(snap_dm))
        x0= snap_dm[permIndx,1:4,-1]/R0
        v0= snap_dm[permIndx,4:,-1]/V0
        multiOut= multi.parallel_map(\
            lambda x: impulse_deltav_plummerint(v_gc,x_gc,lp,
                                                10.**-4./bovy_conversion.mass_in_1010msol(V0,R0), # pretend all 10^6 Msolar, adjust later
                                                numpy.reshape(x0[x],(1,3)),
                                                numpy.reshape(v0[x],(1,3)),
                                                rs=rs,
                                                tmax=0.25*0.9777922212082034/\
                                                    bovy_conversion.time_in_Gyr(V0,R0)),
            range(ndm),numcores=25)
        deltav= numpy.swapaxes(numpy.array(multiOut),1,2).T
        save_pickles(savefilename,deltav)
    else:
        with open(savefilename,'rb') as savefile:
            deltav= pickle.load(savefile)
    return deltav

if __name__ == '__main__':
    # Load GC snapshot
    filename= os.path.join(os.getenv('DATADIR'),'bovy','stream-stream',
                       'gc_evol_unp_atimpact.dat')
    snap_gc= nemo_util.read(filename,swapyz=True)
    # Load GC action-angle
    aa_filename= os.path.join(os.getenv('DATADIR'),'bovy','stream-stream',
                              'gc_evol_unp_aa.dat')
    acfs= calc_aA_sim(None,aa_filename,snap_gc)
    apar= calc_apar(acfs)
    v_gc, x_gc, t_gc, a_gc= vxta_deltav(snap_gc,apar,acfs)
    if sys.argv[1].lower() == 'fullplummer':
        # Inputs: fullplummer 0.50 SAVEFILENAME 100
        calc_fullplummer_deltav(sys.argv[2],sys.argv[3],int(sys.argv[4]),
                                v_gc,x_gc)
