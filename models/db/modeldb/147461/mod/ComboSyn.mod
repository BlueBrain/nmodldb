TITLE  NMDA receptor with Ca influx + AMPA receptor

COMMENT
Written by Eduardo Conde-Sousa, FCUP
econdesousa@gmail.com
and 
Paulo Aguiar, FCUP
pauloaguiar@fc.up.pt
in Dec 2011
ENDCOMMENT

NEURON {
		POINT_PROCESS ComboSyn	
		USEION ca WRITE ica	
		USEION mg READ mgo VALENCE 2
		RANGE tau1, tau2, e_AMPA, i_AMPA
		RANGE g_AMPA
		RANGE tau_rise, tau_decay
		RANGE i_NMDA, g_NMDA, e_NMDA, mg, i2, ica, ca_ratio
		RANGE rr
		NONSPECIFIC_CURRENT i
}
   
UNITS {
		(nA) = (nanoamp)
		(mV) = (millivolt)
		(molar) = (1/liter)
		(mM) = (millimolar)
}    
    
PARAMETER {
    tau1      = 0.1   (ms) <1e-9,1e9>
    tau2      = 5     (ms) <1e-9,1e9>
    e_AMPA    = 0	    (mV)
  	tau_rise  = 5.0   (ms)  <1e-9,1e9>  
    tau_decay = 70.0  (ms)  <1e-9,1e9>  
    e_NMDA    = 0.0   (mV)  : synapse reversal potential
    mgo		    = 1.0   (mM)  : external magnesium concentration
    ca_ratio  = 0.1   (1)   : ratio of Ca current to total current, Burnashev/Sakmann J.Phys.1995 485 403-418)
    rr        = 0.839486356
}
    
    
ASSIGNED {
		v		(mV)
		i		(nA)
		i2		(nA)
		g_AMPA  (umho)
		g_NMDA  (umho)
		factor	(1)
		ica		(nA)
		i_NMDA	(nA)
		i_AMPA	(nA)
  	factor_AMPA
  	factor_NMDA
}

STATE {
		A_AMPA (uS)
		B_AMPA (uS)
		A_NMDA (uS)
		B_NMDA (uS)
}


INITIAL{
		LOCAL tp_NMDA, tp_AMPA
		if (tau_rise/tau_decay > .9999) {
				tau_rise = .9999*tau_decay
		}
		A_NMDA = 0
		B_NMDA = 0
		tp_NMDA = (tau_rise*tau_decay)/(tau_decay-tau_rise)*log(tau_decay/tau_rise)
		factor_NMDA = -exp(-tp_NMDA/tau_rise)+exp(-tp_NMDA/tau_decay)
		factor_NMDA = 1/factor_NMDA
		
		if (tau1/tau2 > .9999) {
				tau1 = .9999*tau2
		}
		A_AMPA = 0
		B_AMPA = 0
		tp_AMPA = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
		factor_AMPA = -exp(-tp_AMPA/tau1) + exp(-tp_AMPA/tau2)
		factor_AMPA = 1/factor_AMPA
}


BREAKPOINT {
		SOLVE state METHOD cnexp
		g_NMDA = B_NMDA-A_NMDA
		i2 = g_NMDA*mgblock(v)*(v-e_NMDA)
		ica = ca_ratio*i2
		i_NMDA = (1-ca_ratio)*i2
		g_AMPA = B_AMPA - A_AMPA
		i_AMPA = g_AMPA*(v - e_AMPA)
  	i=i_NMDA+i_AMPA
}


DERIVATIVE state{
		A_NMDA' = -A_NMDA/tau_rise
		B_NMDA' = -B_NMDA/tau_decay
		A_AMPA' = -A_AMPA/tau1
		B_AMPA' = -B_AMPA/tau2
		
}

FUNCTION mgblock(v(mV)) {
		: from Jahr & Stevens 1990
		mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mgo / 3.57 (mM)))
}

NET_RECEIVE (weight (uS)){
		A_NMDA = A_NMDA + weight*rr*factor_NMDA
		B_NMDA = B_NMDA + weight*rr*factor_NMDA
		A_AMPA = A_AMPA + weight*(1-rr)*factor_AMPA
		B_AMPA = B_AMPA + weight*(1-rr)*factor_AMPA
}
