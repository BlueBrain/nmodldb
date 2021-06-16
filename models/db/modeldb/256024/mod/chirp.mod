TITLE chirp current

COMMENT
-----------------------------------------------------------------------------

    chirp current injection
    ==================================================

 IMPLEMENTATION

  This mechanism is implemented as a nonspecific current defined as a
  point process, mimicking a current-clamp stimulation protocol, injecting
  a chirp waveform with a constant frequency or linearly or exponentially changing frequencies.
  
	f = t-t1	: time relative to start of chirp

	'exponential'
	chirp = amp*sin(f*exp(f*beta)*Finit*pi)
	'linear'
	chirp = amp*sin(pi*beta*f^2)
	'constant'
	chirp = amp*sin(2*pi*Finit*f)

	amp, beta, linear, Finit, t1, dur

 PARAMETERS

  This mechanism takes the following parameters:

	amp	=	0.0 (nA)	: amplitude of injected current. positive values of i depolarize the cell
	t1 = 	0.0 (ms)	: starting time of the stimulation.
	dur =	0.0 (ms)	: duration of chirp current
	Finit =	0.0 (Hz)	: initial frequency of the chirp current.
	beta =	0.0 (Hz/s)	: rate of change of chirp frequency
	ctype =	1 (1)		: chirp type. 0 = constant, 1 = linear (default), 2 = exponential
	
-----------------------------------------------------------------------------
ENDCOMMENT


:INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS chirp
    RANGE amp, t1, dur, Finit, beta, ctype
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp) 
    (mV) = (millivolt)
}

PARAMETER {
	amp	=	1.0	(nA)	: amplitude of injected current
	t1 = 	1000 (ms)	: starting time of the stimulation.
	dur =	20000 (ms)	: duration of chirp current
	Finit =	0.05 (Hz)	: initial frequency of the chirp current.
	beta =	0.24 (Hz/s)	: rate of change of chirp frequency
	ctype =	2	(1)		: chirp type. 0 = constant, 1 = linear, 2 = exponential
}

ASSIGNED {
	i	(nA)			: fluctuating current
}


BREAKPOINT {
	LOCAL tc,  pi, uf, mstos
	:tc is time from start of chirp (s)

	mstos = 1000 (ms/s)	: conversion of ms to s
	uf = 1 (s)			: unit conversion for exponential chirp. beta is in Hz for exponential chirp

	pi=3.14159265358979323846
	:pi=PI
	
    if ((t < t1) || (t > t1+dur)) {  
    	i = 0
    } else {
    	tc=(t-t1)/mstos
    	if (ctype == 0) {
    		i = -amp*sin(2*pi*Finit*tc)
    	} else if (ctype == 1) {
    		i = -amp*sin(pi*beta*(tc+Finit/beta)^2)
    	} else if (ctype == 2) {
    		i = -amp*sin(2*pi*Finit*tc*exp(tc*beta*uf))
    	}
    }
}


