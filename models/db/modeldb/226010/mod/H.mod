TITLE  Hyperpolarisation-activated cation current (H-current)

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
        SUFFIX H
        NONSPECIFIC_CURRENT ih 
        RANGE ghbar,ih,eh
        RANGE minf
		RANGE tau_m
		RANGE vhalfAct,slopeAct
		RANGE vhalfTact,slopeTact
}
 
PARAMETER {
        v   (mV)
        dt  (ms)
		ghbar  = 0.00008 (S/cm2)
        eh = -40 (mV) 
		vhalfAct = -114.7 (mV)
		slopeAct = -12.8
        vhalfTact = -112.7 (mV)
		slopeTact =  6.7
}
 
STATE {
        m
}
 
ASSIGNED {
        ih (mA/cm2)
        minf
	    tau_m 
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        ih = ghbar*m*(v - eh)      
}
 
UNITSOFF

INITIAL {
        m = minf
}

DERIVATIVE states { 
        LOCAL minf,tau_m
        minf = 1/(1 + exp(-(v - vhalfAct)/slopeAct))
		tau_m = 200 + 2000/(1 + exp(-(v - vhalfTact)/slopeTact))
        m' = (minf-m)/tau_m
}
 
UNITSON

