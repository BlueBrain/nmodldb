TITLE  Transient sodium current (NaT-current)

COMMENT
written for NEURON by Antonios Dougalis, 23 Feb 2016, Ulm
based on voltage clamp data from Dougalis et al., 2017 J Compu Neurosci 
ENDCOMMENT

UNITS {
        (S) = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX NaT
        USEION na READ ena WRITE ina
        RANGE gnaTbar,inaT,ina,ena
        RANGE minf,hinf
		RANGE tau_m,tau_h
		RANGE vhalfAct,slopeAct,vhalfInact,slopeInact
		RANGE vhalfTact,slopeTact,vhalfTInact,slopeTInact
		}
 
PARAMETER {
        v   (mV)
        dt  (ms)
		gnaTbar  = 0.02 (S/cm2)
        ena = 50 (mV)
		vhalfAct = -44.0 (mV)
		vhalfInact = -62.0 (mV)
		slopeAct = 4.5
		slopeInact = -6.5
		vhalfTact = -27.9 (mV)
		slopeTact = -6.9
		vhalfTInact = -14.5 (mV)
		slopeTInact = -9.5  
        	
}
 
STATE {
        m h 
}
 
ASSIGNED {
        ina (mA/cm2)
        inaT (mA/cm2)
        minf
		hinf 
		tau_m
		tau_h
		}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        inaT = gnaTbar*m*m*m*h*(v - ena)
        ina=inaT
		
}
 
UNITSOFF

INITIAL {
        m = minf
        h = hinf
        
}

DERIVATIVE states { 
        LOCAL minf,hinf,tau_m,tau_h
        minf = 1/(1 + exp(-(v - vhalfAct)/slopeAct))
		hinf = 1/(1 + exp(-(v - vhalfInact)/slopeInact))
        tau_m = 0.19* (1/(1 + exp(-(v - vhalfTact)/slopeTact))) 
		tau_h = 2.35* (1/(1 + exp(-(v - vhalfTInact)/slopeTInact))) 
		m' = (minf-m)/tau_m
        h' = (hinf-h)/tau_h
        
}
 
UNITSON
