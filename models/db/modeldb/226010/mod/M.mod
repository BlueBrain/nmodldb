TITLE  a low-voltage activated, non-inactivating, potassium current (M-current)

COMMENT
written for NEURON by Antonios Dougalis, 23 Feb 2015, London, UK
based on voltage clamp data from Dougalis et al., 2017 J Comput Neurosci 
ENDCOMMENT

UNITS {
        (S) = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX kM
        USEION k READ ek WRITE ik
        RANGE gkMbar,ik,ek
        RANGE minf,tau_m
		RANGE vhalfAct,slopeAct
		RANGE vhalfTact, slopeTact
}
 
PARAMETER {
        v   (mV)
        dt  (ms)
		gkMbar = 0.001 (S/cm2)
        ek  = -73.0  (mV)
		vhalfAct = -35 (mV)
		slopeAct = 8.5
		vhalfTact = -27.9 (mV)
		slopeTact = -6.9
}
 
STATE {
        m 
}
 
ASSIGNED {
        ik (mA/cm2)
        il (mA/cm2)
        minf
		tau_m
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        ik = gkMbar*m*(v - ek)             
}
 
UNITSOFF

INITIAL {
        m = minf
}

DERIVATIVE states { 
        LOCAL minf,tau_m
        minf = 1/(1 + exp(-(v - vhalfAct)/slopeAct))
        tau_m = 50 + 48/(1 + exp(-(v - vhalfTact)/slopeTact))
		m' = (minf-m)/tau_m
}
 
FUNCTION boltz(v,y,z) {
                boltz = 1/(1 + exp(-(v - y)/z))
}
 
UNITSON

