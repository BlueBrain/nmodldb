COMMENT

Ih current	 - hyperpolarization-activated nonspecific Na and K channel
		 - contributes to the resting membrane potential
		 - controls the afterhyperpolarization

Conductance density is modified by Vlad as through database

Reference:
Saraga, F. (2003). Active dendrites and spike propagation in multi-compartment models of oriens-lacunosum/moleculare hippocampal interneurons. The Journal of physiology, 552(3), 673-689.

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Ih
        USEION h READ eh WRITE ih VALENCE 1
        RANGE gkhbar,ih,g
        GLOBAL rinf
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        p = 5 (degC)
        dt (ms)
	gkhbar = 0.000005 (mho/cm2)
        eh = -32.9 (mV)

	t1 = 1 (ms)
	t2 = 1
	t3 = -0.116 (1/mV)
	t4 = -17.9
	t5 = 1
	t6 = 0.09 (1/mV)
	t7 = -1.84
	t8 = 100 (ms)
}
 
STATE {
        r
}
 
ASSIGNED {
        ih (mA/cm2)
	rinf 
	tau_r	(ms)
	g (S/cm2)
}
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        ih = gkhbar*r*(v - eh)
	g = gkhbar*r
}
 
INITIAL {
	rates(v, t1, t2, t3, t4, t5, t6, t7, t8)
	r = rinf
}

DERIVATIVE state { :Computes state variable h at current v and dt.
	rates(v, t1, t2, t3, t4, t5, t6, t7, t8)
	r' = (rinf - r)/tau_r
}
UNITSOFF

PROCEDURE rates(v(mV), t1 (ms), t2, t3(1/mV), t4, t5, t6(1/mV), t7, t8(ms)) {
	rinf = 1/(1 + exp((v+84)/10.2))
	tau_r = t1/(t2*exp(t3*v+t4)+t5*exp(t6*v+t7)) + t8 + 1e-8
}
 
UNITSON
