Notes on how to run the simulations using NEMO
===============================================

We simulate the King profile that represents a globular cluster
starting a time 10.25 time units ago, in gyrfalcON units. This
corresponds to 10.25*0.9777922212082034 ~ 10.02 Gyr. At this time, the
GC's phase-space coordinates are::

     X, Y, Z = 29.936262232141431, -0.93933856619906753, -0.93918238832906931
     vX, vY, vZ=  14.340618161448306, 105.52412254413353, 105.47140831277834

We then create the King cluster::

   mkking out=gc.nemo nbody=10000 W0=5. mass=20000 r_c=0.013 WD_units=t

and we shift it to the correct position and velocicy given above::

    snapshift gc.nemo gc_shifted.nemo rshift=29.936262232141431,-0.93933856619906753,-0.93918238832906931 vshift=14.340618161448306,105.52412254413353,105.47140831277834

Unperturbed stream
--------------------

We simulate the unperturbed stream evolution as::

   gyrfalcON in=gc_shifted.nemo out=$DATADIR/bovy/stream-stream/gc_evol_unp.nemo tstop=10.750 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > gc_evol_unp.log 2>&1

Perturbed stream
-----------------

We simulate the perturbed stream in a few different steps:

   * We simulate the GC separately until 0.125 time units before the
     impact

   * We simulate the DM halo separately until 0.125 time units before
     the impact; the DM halo is started at different times to create
     different length tails

   * We then combine the two snapshots and integrate them forward for
     0.25 time units

   * Then we separate out the GC again and integrate it forward for
     another 0.375 time units to end up at the same time as the
     unperturbed stream above.

Unperturbed evolution up to impact
+++++++++++++++++++++++++++++++++++

Do::

	gyrfalcON in=gc_shifted.nemo out=$DATADIR/bovy/stream-stream/gc_evol_untilimpact.nemo tstop=10.10.125 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > gc_evol_untilimpact.log 2>&1