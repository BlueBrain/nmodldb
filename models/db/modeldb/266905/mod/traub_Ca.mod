TITLE Ca channel Traub

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (mS) = (millisiemens)
}

NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _gCa
    
    SUFFIX gCa
    USEION ca READ eca WRITE ica
    RANGE gmax, g, i
    GLOBAL rss
	: GLOBAL s_inf, tau_s, r_inf, tau_r
}

PARAMETER {
    gmax=0.004 (mho/cm2)  : value at soma
    rss=0	(1)
}

STATE {
    s
    r
}

ASSIGNED {
    v (mV)
    eca (mV)
    ica		(mA/cm2)
    s_inf	(1)
    r_inf	(1)
    tau_s	(ms)
    tau_r	(ms)
    g		(S/cm2)
    i		(mA/cm2)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    if (rss==1) {
    	g = gmax*s^2
    } else {
	    g = gmax*s^2*r
	}
    i = g*(v-eca)
    ica = i
}

INITIAL {
    rates(v)
    s = s_inf
    r = r_inf
}


FUNCTION alpha_s(v(mV)) {
    alpha_s = 1.6/(1 + exp( -0.072*(v-65) ) )
}

FUNCTION beta_s(v(mV)) {
    beta_s = 0.02*(v - 51.1)/(exp( (v-51.1)/5 ) - 1)
}

FUNCTION alpha_r(v(mV)) {
    if ( v<=0 ){
		alpha_r = 0.005
    } else {
		alpha_r = exp(-v/20)/200
    }
    
}

FUNCTION beta_r(v(mV)) {
    if ( v<=0 ) {
		beta_r = 0
    } else {
		beta_r = 0.005 - exp(-v/20)/200
    }
}

DERIVATIVE states {  
    rates(v)
    s' = (s_inf - s)/tau_s
    r' = (r_inf - r)/tau_r 
}

PROCEDURE rates(v (mV)) { :callable from hoc
    LOCAL alpha, beta
    TABLE s_inf, tau_s, r_inf, tau_r  : DEPEND vhalfn, vhalfl, tlmax, tnmax
    FROM -20 TO 130 WITH 750
    
    alpha = alpha_s(v)
    beta = beta_s(v)
    tau_s = 1/(alpha + beta)     
    s_inf = alpha*tau_s
    
   
    alpha = alpha_r(v)
    beta = beta_r(v)
    tau_r = 1/(alpha + beta)     
    r_inf = alpha*tau_r
   
}

