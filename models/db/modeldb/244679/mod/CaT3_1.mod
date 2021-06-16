TITLE T-type VGCC
COMMENT
T-type VGCC. Starting from the formulation for T-type VGCCs in Anwar et al., 2012, we multiplied the
activation curve by a sigmoid function to account for the observed activation at hyperpolarised and
intermediate states. The activation time was decreased by 70% while the inactivation time was
doubled to fit the rise and decay of the experimental V m trace at the hyperpolarised state.

Current Model Reference: Karima Ait Ouares , Luiza Filipis , Alexandra Tzilivaki , Panayiota Poirazi , Marco Canepari
(2018) Two distinct sets of Ca 2+ and K + channels
are activated at different membrane potential by the climbing fibre synaptic potential in Purkinje neuron dendrites.
PubMed link:

Contact: Filipis Luiza (luiza.filipis@univ-grenoble-alpes.fr)
ENDCOMMENT


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
        SUFFIX CaT3_1
        USEION ca READ cai, cao WRITE ica VALENCE 2
        RANGE g, pcabar, minf, taum, hinf, tauh
	RANGE ica, m ,h, t1,t2

    }

UNITS {
        (molar) = (1/liter)
        (mV) =  (millivolt)
        (mA) =  (milliamp)
        (mM) =  (millimolar)

}

CONSTANT {
	F = 9.6485e4 (coulombs)
	R = 8.3145 (joule/kelvin)
}

PARAMETER {
        v               (mV)
        t1=0.6
        t2=55

        celsius (degC)
        eca (mV)
	pcabar  = 2.5e-4 (cm/s)
        cai  (mM)           : adjusted for eca=120 mV
	cao  (mM)
	
	v0_m_inf = -55 (mV)
	v0_h_inf = -72 (mV)
	k_m_inf = -5 (mV)
	k_h_inf = 7  (mV)
	
	C_tau_m = 1
	A_tau_m = 1.0
	v0_tau_m1 = -40 (mV)
	v0_tau_m2 = -102 (mV)
	k_tau_m1 = 9 (mV)
	k_tau_m2 = -18 (mV)
	
	C_tau_h = 15
	A_tau_h = 1.0
	v0_tau_h1 = -32 (mV)
	k_tau_h1 = 7 (mV)
	
    }
    

STATE {
        m h
}

ASSIGNED {
        ica     (mA/cm2)
	g        (coulombs/cm3) 
        minf
        taum   (ms)
        hinf
        tauh   (ms)
	T (kelvin)
	E (volt)
	zeta
}

BREAKPOINT {
	SOLVE castate METHOD cnexp 

        ica = (1e3) *pcabar*m*m *h * g
}

DERIVATIVE castate {
        evaluate_fct(v)

        m' = (minf - m) / taum
        h' = (hinf - h) / tauh
}

FUNCTION ghk( v (mV), ci (mM), co (mM), z )  (coulombs/cm3) {
    E = (1e-3) * v
      zeta = (z*F*E)/(R*T)


    if ( fabs(1-exp(-zeta)) < 1e-6 ) {
        ghk = (1e-6) * (z*F) * (ci - co*exp(-zeta)) * (1 + zeta/2)
    } else {
        ghk = (1e-6) * (z*zeta*F) * (ci - co*exp(-zeta)) / (1-exp(-zeta))
    }
}


UNITSOFF
INITIAL {
	
	T = kelvinfkt (celsius)

        evaluate_fct(v)
        m = minf
        h = hinf
}

PROCEDURE evaluate_fct(v(mV)) { 

        
        minf = ((1/(1+exp(-t1*(v+t2))))/(1 + exp(-(v+51)/6)))
        hinf = 1.0 / ( 1 + exp((v - v0_h_inf)/k_h_inf) )
        if (v<=-90) {
	taum = 1
	} else {
	taum = 0.2 / (exp((v +40)/ 9) + exp(-(v +108)/18))+0.3
	}
	tauh = 2*( C_tau_h + A_tau_h / exp((v - v0_tau_h1)/k_tau_h1) ) 
	g = ghk(v, cai, cao, 2)
}

FUNCTION kelvinfkt( t (degC) )  (kelvin) {
    kelvinfkt = 273.19 + t
}

UNITSON
