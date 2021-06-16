
COMMENT
26 Ago 2002 Modification of original channel to allow variable time step and to correct an initialization error.
    Done by Michael Hines(michael.hines@yale.e) and Ruggero Scorcioni(rscorcio@gmu.edu) at EU Advance Course in Computational Neuroscience. Obidos, Portugal


na.mod

Sodium channel, Hodgkin-Huxley style kinetics.  

Kinetics were fit to data from Huguenard et al. (1988) and Hamill et
al. (1991)

qi is not well constrained by the data, since there are no points
between -80 and -55.  So this was fixed at 5 while the thi1,thi2,Rg,Rd
were optimized using a simplex least square proc

voltage dependencies are shifted approximately from the best
fit to give higher threshold

Author: Zach Mainen, Salk Institute, 1994, zach@salk.edu

tadj, the temperature adjustment was removed from instantaneous conductance term
in BREAKPOINT

steady-state inactivation was changed to more usual form: a/(a+b) and
inactivation time constant was significantly reduced to reflect recent data
from Kole, .... Stuart '08 Nat Neurosci.

Corey Acker, July 2008, Neuroscience, UConn Health Center

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX na
	USEION na READ ena WRITE ina
	RANGE m, h, gna, gbar
	GLOBAL tha, thi1, thi2, qa, qi, qinf, thinf
	RANGE minf, hinf, mtau, htau
	GLOBAL Ra, Rb, Rd, Rg
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
	gbar = 1000   	(pS/um2)	: 0.12 mho/cm2
:	vshift = -10	(mV)		: voltage shift (affects all)
	vshift = 0	(mV)		: voltage shift (affects all)
								
	tha  = -35	(mV)		: v 1/2 for act		(-42)
	qa   = 9	(mV)		: act slope		
	Ra   = 0.182	(/ms)		: open (v)		
	Rb   = 0.124	(/ms)		: close (v)		

:	thi1  = -50	(mV)		: v 1/2 for inact 	
      thi1 = -65 (mV)
:	thi2  = -75	(mV)		: v 1/2 for inact 	
      thi2 = -65 (mV)
	qi   = 6	(mV)	        : inact tau slope
:	thinf  = -65	(mV)		: inact inf slope	
	thinf  = -65	(mV)		: inact inf slope	
	qinf  = 6.2	(mV)		: inact inf slope
:	Rg   = 0.0091	(/ms)		: inact (v)	
	Rg   = 0.02	(/ms)		: inact (v)	
	Rd   = 0.024	(/ms)		: inact recov (v) 

	temp = 23	(degC)		: original temp 
	q10  = 2.3			: temperature sensitivity

	v 		(mV)
	dt		(ms)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	gna		(pS/um2)
	ena		(mV)
	minf 		hinf
	mtau (ms)	htau (ms)
	tadj
}
 

STATE { m h }

INITIAL { 
	trates(v+vshift)
	m = minf
	h = hinf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
:        gna = tadj*gbar*m*m*m*h : originally included tadj
        gna = gbar*m*m*m*h
	ina = (1e-4) * gna * (v - ena)
} 

LOCAL mexp, hexp 

DERIVATIVE states {   :Computes state variables m, h, and n 
        trates(v+vshift)      :             at the current v and dt.
        m' =  (minf-m)/mtau
        h' =  (hinf-h)/htau
}

PROCEDURE trates(v) {  
                      
        
        TABLE minf,  hinf, mtau, htau
	DEPEND  celsius, temp, Ra, Rb, Rd, Rg, tha, thi1, thi2, qa, qi, qinf
	
	FROM vmin TO vmax WITH 199

	rates(v): not consistently executed from here if usetable == 1

:        tinc = -dt * tadj

:        mexp = 1 - exp(tinc/mtau)
:        hexp = 1 - exp(tinc/htau)
}


PROCEDURE rates(vm) {  
        LOCAL  a, b

	a = trap0(vm,tha,Ra,qa)
	b = trap0(-vm,-tha,Rb,qa)

        tadj = q10^((celsius - temp)/10)

	mtau = 1/tadj/(a+b)
	minf = a/(a+b)

		:"h" inactivation 

	a = trap0(-vm,-thi1,Rd,qi)
	b = trap0(vm,thi2,Rg,qi)
	htau = 1/tadj/(a+b)
:	hinf = 1/(1+exp((vm-thinf)/qinf))
      hinf = a/(a+b)
}


FUNCTION trap0(v,th,a,q) {
	if (fabs(v/th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	




