
Commands to pre-compute computationally-expensive parts of the notebook
------------------------------------------------------------------------

The kicks using orbit-integration of Plummer spheres are computed using::

    python calc_deltav.py fullplummer 0.50 $DATADIR/bovy/stream-stream/deltav_0.50_fullplummer.sav 300
