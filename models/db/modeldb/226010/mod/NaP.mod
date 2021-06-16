TITLE  Persistent sodium current (NaP-current)

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
        SUFFIX NaP
        USEION na READ ena WRITE ina
        RANGE gnaPbar,ena,ina
        RANGE minf
		RANGE tau_m
		RANGE vhalfAct,slopeAct
}
 
PARAMETER {
        v   (mV)
        dt  (ms)
		gnaPbar  = 0.00002  (S/cm2)
        ena = 50 (mV)
        tau_m = 0.1 (ms)
        vhalfAct = -57 (mV)
		slopeAct = 3.5
}
 
STATE {
        m
}
 
ASSIGNED {
        ina (mA/cm2)
        minf
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        ina = gnaPbar*m*(v - ena)      
}
 
UNITSOFF

INITIAL {
        m = minf
}

DERIVATIVE states { 
        LOCAL minf
        minf = 1/(1 + exp(-(v - vhalfAct)/slopeAct))
        m' = (minf-m)/tau_m
}
 

 
UNITSON

