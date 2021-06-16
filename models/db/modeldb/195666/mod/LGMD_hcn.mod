TITLE Ih channel for LGMD with cyclic nucleotide enhancement


UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

: gmax and g are range variables (i.e., can change in different compartments
: while e is global
NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _h
    SUFFIX hcn

    NONSPECIFIC_CURRENT i
	USEION cn READ cni VALENCE 1

    : these variables will be accessed as compartment.rangevar_h
    : note: to make the channel constant available add the following
    : to the next line: vhalf, s1, s2, tau_max
    RANGE gmax, tau, g, taumax, i

    : this will be accessed as e_h, taumax_h, vhalf_h, s1_h, s2_h 
    GLOBAL e, taumin, vhalf, s1, s2
}

PARAMETER {
    : all parameterss adjustable in hoc files
    gmax= 0.001 (S/cm2)
    e = -37 (mV)
    vhalf = -78.2 (mV)
    s1 = -12.5 (mV)
    s2 = 13 (mV)
    taumax = 1350 (ms)
    taumin = 10 (ms)
    
   	cnvm=15		(mV)
	lcp=2		(1)
	kD=6e-4		(mM)
	taucn=10	(ms)

}

ASSIGNED { 
    v (mV)
    cni (mM)
    
    i (mA/cm2)
   	vs	(mV)
    ninf
    tau (ms)
    g (S/cm2)
    
}

STATE {
    n
   	ov	(mV)
	ovs	(mV)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    n = 1/(1+(exp((vhalf-ov+ovs)/s1)))
    g  = gmax*n
    i  = gmax*n*(v-e)
}

INITIAL {
    settables(v)
    lci(cni)
    n = 1/(1+(exp((vhalf-v+vs)/s1)))
   	ov = v
	ovs = vs
}

DERIVATIVE states { 
    settables(v)
    lci(cni)
    
    ov' = (v-ov)/tau
	ovs' = (vs-ovs)/taucn
	
    :n' = (ninf - n)/tau
}


PROCEDURE settables(v (mV)) {
    :LOCAL vhalf, s1, s2, taumax
    :local variables take units of right hand side, see below

    TABLE tau DEPEND vhalf, s1, s2, taumax
          FROM -150 TO 50 WITH 750

    : steady-state activation of Ih in mV
    :vhalf = -77.8 (mV)
    :s1 = 13.8 (mV)
    ninf = 1/(1+(exp((vhalf-v)/s1)))

    : steady-state Ih time constant
    : slope in mV and time constant in ms
    :s2 = 19.7 (mV)
    :taumax = 1071.1 (ms)
    :tau = 2*taumax/( exp((v-vhalf)/s2) + exp((vhalf-v)/s2) )
    tau = 4*taumax/(1+exp((vhalf-v)/s2))*ninf+taumin

}


PROCEDURE lci(cni (mM)) { :callable from hoc
	TABLE vs DEPEND lcp, kD, cnvm
          FROM 0 TO 0.01 WITH 500
    
	vs = cnvm-cnvm/(1+(cni/kD)^lcp)
    
}
