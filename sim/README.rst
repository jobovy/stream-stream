Notes on how to run the simulations using NEMO
===============================================

We simulate the King profile that represents a globular cluster
starting a time 10.25 time units ago, in gyrfalcON units. This
corresponds to 10.25*0.9777922212082034 ~ 10.02 Gyr. At this time, the
GC's phase-space coordinates are::

     X, Y, Z = 29.577809361172694, -2.4079802843036058, -2.4053154109323742
     vX, vY, vZ=  36.904368548683991, 104.25395168642832, 103.9042939586286

We then create the King cluster::

   mkking out=gc.nemo nbody=100000 W0=5. mass=100000 r_c=0.013 WD_units=t

and we shift it to the correct position and velocity given above
(remember that we need to flip y and z for using the logarithmic
potential in NEMO)::

    snapshift gc.nemo gc_shifted.nemo rshift=29.577809361172694,-2.4053154109323742,-2.4079802843036058 vshift=36.904368548683991,103.9042939586286,104.25395168642832

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

Unperturbed evolution up to impact
+++++++++++++++++++++++++++++++++++

Do::

	gyrfalcON in=gc_shifted.nemo out=$DATADIR/bovy/stream-stream/gc_evol_untilimpact.nemo tstop=10.125 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > gc_evol_untilimpact.log 2>&1
	snaptrim in=$DATADIR/bovy/stream-stream/gc_evol_untilimpact.nemo out=-  times=10.125 | snaptime in=- out=gc_beforeimpact.nemo
DM perturber
+++++++++++++

The DM perturber is set on an orbit that interacts with the GC at
10.25 time units of the GC evolution. We model the DM perturber as a
Plummer sphere with mass 10^8 and rs = 0.625 kpc. We generate the DM
halo as::

     mkplum out=dm.nemo nbody=100000 r_s=0.625 mass=100000000 WD_units=t

If the DM halo starts out 0.50 time units ago, we need to shift it
as::

	snapshift dm.nemo dm_shifted_0.50.nemo rshift=7.9948807918667892,0.9626499711978439,-11.344426189362455 vshift=-100.4814119301404,-151.74007155179646,-84.037056563181252

If the DM halo starts out 0.375 time units ago, we need to shift it
as::

	snapshift dm.nemo dm_shifted_0.375.nemo rshift=-8.5662119198502609,-6.4724441270571687,0.50120970587644287 vshift=-50.169060756613071,112.96727941016003,214.43730878238947

If the DM halo starts out 0.25 time units ago, we need to shift it
as::

	snapshift dm.nemo dm_shifted_0.25.nemo rshift=4.8965792014085778,10.033551072654888,8.9181416272959186 vshift=149.50458227428786,21.235006395774576,-97.714578287010966

If the DM halo starts out 0.125 time units ago, we need to shift it
as::

	snapshift dm.nemo dm_shifted_0.125.nemo rshift=5.2570183967301833,-4.5786210693408176,-7.8848620802774407 vshift=-195.34540991555195,-155.93944808745755,-51.644977186288536

We integrate these snapshots until 0.125 time units before the
impact (we can skip the 0.125 one, because it's there already)::

	gyrfalcON in=dm_shifted_0.50.nemo out=$DATADIR/bovy/stream-stream/dm_evol_0.50_untilimpact.nemo tstop=0.375 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > dm_evol_0.50_untilimpact.log 2>&1
	gyrfalcON in=dm_shifted_0.375.nemo out=$DATADIR/bovy/stream-stream/dm_evol_0.375_untilimpact.nemo tstop=0.25 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > dm_evol_0.375_untilimpact.log 2>&1
	gyrfalcON in=dm_shifted_0.25.nemo out=$DATADIR/bovy/stream-stream/dm_evol_0.25_untilimpact.nemo tstop=0.125 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > dm_evol_0.25_untilimpact.log 2>&1

and we get the final outputs with ``snaptrim``::

    	snaptrim in=$DATADIR/bovy/stream-stream/dm_evol_0.50_untilimpact.nemo  out=- times=0.375 | snaptime in=- out=dm_0.50_beforeimpact.nemo
    	snaptrim in=$DATADIR/bovy/stream-stream/dm_evol_0.375_untilimpact.nemo out=- times=0.25 | snaptime in=- out=dm_0.375_beforeimpact.nemo
    	snaptrim in=$DATADIR/bovy/stream-stream/dm_evol_0.25_untilimpact.nemo out=-  times=0.125 | snaptime in=- out=dm_0.25_beforeimpact.nemo

These are added to the GC simulation with ``snapstack``::

      snapstack in1=gc_beforeimpact.nemo in2=dm_0.50_beforeimpact.nemo out=gcdm_0.50_beforeimpact.nemo
      snapstack in1=gc_beforeimpact.nemo in2=dm_0.375_beforeimpact.nemo out=gcdm_0.375_beforeimpact.nemo
      snapstack in1=gc_beforeimpact.nemo in2=dm_0.25_beforeimpact.nemo out=gcdm_0.25_beforeimpact.nemo
      snapstack in1=gc_beforeimpact.nemo in2=dm_shifted_0.125.nemo out=gcdm_0.125_beforeimpact.nemo

Now we can run these forward::

	gyrfalcON in=gcdm_0.50_beforeimpact.nemo out=$DATADIR/bovy/stream-stream/gcdm_evol_0.50_impact.nemo tstop=0.250 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > gcdm_evol_0.50_impact.log 2>&1
	gyrfalcON in=gcdm_0.375_beforeimpact.nemo out=$DATADIR/bovy/stream-stream/gcdm_evol_0.375_impact.nemo tstop=0.250 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > gcdm_evol_0.375_impact.log 2>&1
	gyrfalcON in=gcdm_0.25_beforeimpact.nemo out=$DATADIR/bovy/stream-stream/gcdm_evol_0.25_impact.nemo tstop=0.250 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > gcdm_evol_0.25_impact.log 2>&1
	gyrfalcON in=gcdm_0.125_beforeimpact.nemo out=$DATADIR/bovy/stream-stream/gcdm_evol_0.125_impact.nemo tstop=0.250 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > gcdm_evol_0.125_impact.log 2>&1

The final output files are then generated as::

    snaptrim in=$DATADIR/bovy/stream-stream/gcdm_evol_0.50_impact.nemo out=- times=0.25 | s2a - $DATADIR/bovy/stream-stream/gcdm_evol_0.50_afterimpact.dat
    snaptrim in=$DATADIR/bovy/stream-stream/gcdm_evol_0.375_impact.nemo out=- times=0.25 | s2a - $DATADIR/bovy/stream-stream/gcdm_evol_0.375_afterimpact.dat
    snaptrim in=$DATADIR/bovy/stream-stream/gcdm_evol_0.25_impact.nemo out=- times=0.25 | s2a - $DATADIR/bovy/stream-stream/gcdm_evol_0.25_afterimpact.dat
    snaptrim in=$DATADIR/bovy/stream-stream/gcdm_evol_0.125_impact.nemo out=- times=0.25 | s2a - $DATADIR/bovy/stream-stream/gcdm_evol_0.125_afterimpact.dat