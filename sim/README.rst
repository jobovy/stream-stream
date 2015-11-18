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

   * Then we separate out the GC again and integrate it forward for
     another 0.375 time units to end up at the same time as the
     unperturbed stream above.

Unperturbed evolution up to impact
+++++++++++++++++++++++++++++++++++

Do::

	gyrfalcON in=gc_shifted.nemo out=$DATADIR/bovy/stream-stream/gc_evol_untilimpact.nemo tstop=10.125 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > gc_evol_untilimpact.log 2>&1

DM perturber
+++++++++++++

The DM perturber is set on an orbit that interacts with the GC at
10.25 time units of the GC evolution. We model the DM perturber as a
Plummer sphere with mass 10^8 and rs = 0.625 kpc. We generate the DM
halo as::

     mkplum out=dm.nemo nbody=100000 r_s=0.625 mass=100000000 WD_units=t

If the DM halo starts out 10.25 time units ago, we need to shift it
as::

	snapshift dm.nemo dm_shifted_10.25.nemo rshift=-3.5959531011355614,-0.3076301806847701,-12.676070584635923 vshift=-122.87209051890167,-159.27470248909279,70.700287494516076

If the DM halo starts out 5.125 time units ago, we need to shift it
as::

	snapshift dm.nemo dm_shifted_5.125.nemo rshift=-7.2904594277232393,-8.0993202639480195,-2.6997962555439798 vshift=-112.18960536729145,0.28635224520599717,206.96653603845965

If the DM halo starts out 1.125 time units ago, we need to shift it
as::

	snapshift dm.nemo dm_shifted_1.125.nemo rshift=-7.4115651013300745,-5.8569155874230869,-4.9401086177590141 vshift=-57.329247868892914,-129.90969599626123,206.23950245614526

f the DM halo starts out 0.50 time units ago, we need to shift it
as::

	snapshift dm.nemo dm_shifted_0.50.nemo rshift=7.9948807918667892,0.9626499711978439,-11.344426189362455 vshift=-100.4814119301404,-151.74007155179646,-84.037056563181252

If the DM halo starts out 0.375 time units ago, we need to shift it
as::

	snapshift dm.nemo dm_shifted_0.375.nemo rshift=-8.5662119198502609,-6.4724441270571687,0.50120970587644287 vshift=-50.169060756613071,112.96727941016003,214.43730878238947

If the DM halo starts out 0.25 time units ago, we need to shift it
as::

	snapshift dm.nemo dm_shifted_0.25.nemo rshift=4.8965792014085778,10.033551072654888,8.9181416272959186 vshift=149.50458227428786,21.235006395774576,-97.714578287010966

We integrate these snapshots until 0.125 time units before the
impact::

	gyrfalcON in=dm_shifted_10.25.nemo out=$DATADIR/bovy/stream-stream/dm_evol_10.25_untilimpact.nemo tstop=10.125 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > dm_evol_10.25_untilimpact.log 2>&1
	gyrfalcON in=dm_shifted_5.125.nemo out=$DATADIR/bovy/stream-stream/dm_evol_5.125_untilimpact.nemo tstop=5. eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > dm_evol_5.125_untilimpact.log 2>&1
	gyrfalcON in=dm_shifted_1.125.nemo out=$DATADIR/bovy/stream-stream/dm_evol_1.125_untilimpact.nemo tstop=1. eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > dm_evol_1.125_untilimpact.log 2>&1
	gyrfalcON in=dm_shifted_0.50.nemo out=$DATADIR/bovy/stream-stream/dm_evol_0.50_untilimpact.nemo tstop=0.375 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > dm_evol_0.50_untilimpact.log 2>&1
	gyrfalcON in=dm_shifted_0.375.nemo out=$DATADIR/bovy/stream-stream/dm_evol_0.375_untilimpact.nemo tstop=0.25 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > dm_evol_0.375_untilimpact.log 2>&1
	gyrfalcON in=dm_shifted_0.25.nemo out=$DATADIR/bovy/stream-stream/dm_evol_0.25_untilimpact.nemo tstop=0.125 eps=0.0015 step=0.125 kmax=6 Nlev=10 fac=0.01 accname=LogPot accpars=0,48400.,0.,1.0,0.9 > dm_evol_0.25_untilimpact.log 2>&1
