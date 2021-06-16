TITLE  A-type transient potassium current (A-current)

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
        SUFFIX kA
        USEION k READ ek WRITE ik 
        RANGE gkAbar,ek,ik, ikA
        RANGE minf,hinf
		RANGE tau_m, tau_h
		RANGE vhalfAct,slopeAct,vhalfInact,slopeInact
		RANGE vhalfTact,slopeTact,vhalfTInact,slopeTInact
}
 
PARAMETER {
        v   (mV)
        dt  (ms)
		gkAbar  = 0.006 (S/cm2)
        ek = -73 (mV) 
		vhalfAct = -57.5 (mV)
		vhalfInact = -93.0 (mV)
		slopeAct = 7.7
		slopeInact = -6.1
		vhalfTact = -68.8 (mV)
		vhalfTinact = -24.6 (mV)
		slopeTact = -5.0
		slopeTinact = -8.6
}
 
STATE {
        m h
}
 
ASSIGNED {
        ik (mA/cm2)
		ikA (mA/cm2)
        minf
		hinf
	    tau_m 
		tau_h
		
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        ikA = gkAbar*m*h*(v - ek)
ik=ikA		
}
 
UNITSOFF

INITIAL {
        m = minf
		h = hinf
}

DERIVATIVE states { 
        LOCAL minf,hinf,tau_m, tau_h
        minf = 1/(1 + exp(-(v - vhalfAct)/slopeAct))
		hinf = 1/(1 + exp(-(v - vhalfInact)/slopeInact))
        tau_m = 0.95 + 4.45/(1 + exp(-(v - vhalfTact)/slopeTact))
        tau_h = 128 + 77/(1 + exp(-(v - vhalfTinact)/slopeTinact))
		m' = (minf-m)/tau_m
        h' = (hinf-h)/tau_h
}
 
UNITSON

