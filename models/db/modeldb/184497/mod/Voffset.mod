COMMENT

    Voffset.mod

    Shifts voltage by the steady state level just before the onset of a current step.  This is useful when fitting passive/subthreshold parameters.

    Implemented by Christina Weaver, 2007 (christina.weaver@mssm.edu)

ENDCOMMENT

NEURON {
  SUFFIX offst
  RANGE Voff, Vshift, Vsum
}

PARAMETER { 
	Vraise = -71.3	(mV)	: how much to shift voltage overall
	on = 300	(ms)	: onset time of current step
	W = 50		(ms)	: length of window to average over
	we = 10		(ms)	: time before step onset to end averaging window
	didAvg = 0	(1)	: flag to minimize computation
}

UNITS { 
	(mV) = (millivolt) 
}

ASSIGNED {
  v (millivolt)
  : t (ms)
	Vsum		(mV)
	npts		(1)
	Voff		(mV)
	Vshift 		(mV)
}

INITIAL {
	Vsum = 0
	Vshift = 0
	npts = 0
	Voff = v - Voff + Vraise
        didAvg = 0
: printf("Set npts = %g, Voff %g\n",npts,Voff)
}

BREAKPOINT {
    if( t >= on-we-W && t <= on-we ) {
        Vsum = Vsum + v
        npts = npts + 1
: printf("time %g\tv %g\tVsum %g\tnpts %g\n",t,v,Vsum,npts)
    } 
    if( didAvg == 0 && t > on-we ) {		: past averaging window
        didAvg = 1
        Voff = Vsum / npts
:printf("***TIME %g\tFound offset %g = %g / %g\n",t,Voff,Vsum,npts)
    }
    Vshift = v - Voff + Vraise
:     if( t > on-we ) {
: printf("SH time %g\tVshift %g = v %g\t- Voff %g\t+ Vraise %g\n",t,Vshift,v,Voff,Vraise)
:    } 
}
