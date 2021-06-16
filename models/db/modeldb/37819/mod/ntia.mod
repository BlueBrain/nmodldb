: $Id: ntia.mod,v 1.3 2002/11/08 22:05:20 billl Exp $
TITLE rapidly inactivating potassium current
:
:   K+ current responsible for blocking rebound low threshold spikes (LTS)
:   LOCAL GABAERGIC INTERNEURONS IN THE THALAMUS
:   Differential equations 
:
:   Model of Huguenard & McCormick, J Neurophysiol 68: 1373-1383, 1992.
:   The kinetics is described by standard equations (NOT GHK)
:   using a m4h format, according to the voltage-clamp data
:   of Huguenard, Coulter & Prince, J Neurophysiol.
:   66: 1304-1315, 1991.
:
:    - Kinetics adapted to fit the A-channel of interneuron
:    - Q10 changed to 5 and 3
:    - Time constant tau_m and tau_h from experimental data (from TC)
:    - shift parameter for fitting interneuron data, according to the
:    - voltage-clamp data from premature rat by Pape et al. J.
:    - Physiol. 1994. 
:
:   ACTIVATION FUNCTIONS FROM EXPERIMENTS (NO CORRECTION)
:
:   Reversal potential taken from Nernst Equation
:
:   Written by Jun Zhu, University of Wisconsin, August 19, 1994, at MBL, Woods Hole, MA
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iao
	USEION k3 WRITE ik3 VALENCE 1
	RANGE gabar, i
	GLOBAL m_inf, tau_m, h_inf, tau_h, n_inf, tau_n, shm, shh, shn, hx, nx
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
}

PARAMETER {
	v		(mV)
	celsius	= 36	(degC)
	erev	= -95 	(mV)
	gabar	= 0.0	(mho/cm2)
	shm	= 15 	(mV)
	shh	= 15 	(mV)
	shn	= 15 	(mV)
	hx	= 0 	(mV)
	nx	= 0 	(mV)
}

STATE {
	m h n
}

ASSIGNED {
	ik3	(mA/cm2)
	i	(mA/cm2)
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
	n_inf
        tau_n	(ms)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = gabar * (m*m*m*m*h * (v-erev) * 0.6 + m*m*m*m*n * (v-erev) * 0.4)
        ik3 = i
}

DERIVATIVE states {
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h 
	n' = (h_inf - n) / tau_n 
}

UNITSOFF
INITIAL {
	evaluate_fct(v)
	m = m_inf
	h = h_inf
	n = n_inf
:
:   Activation functions and kinetics were obtained from
:   Huguenard & McCormick, and were at 35.5 deg.

}

PROCEDURE evaluate_fct(v(mV)) { 
:
:   Time constants were obtained from Huguenard & McCormick
:   not sure about 7.4 and 5.0
:

	m_inf = 1.0 / ( 1 + exp(-(v+shm+34)/7.4) )
	h_inf = 1.0 / ( 1 + exp((v+shh+78)/5.0) )
	n_inf = 1.0 / ( 1 + exp((v+shn+78)/5.0) )

	tau_m = ( 1.0 / ( exp((v+shm+35.8)/19.7) + exp(-(v+shm+79.7)/12.7) ) + 0.37 )  

        if (v < -80-shh) {
		tau_h = ( 1.0 / ( exp((v+shh+46)/5) + exp(-(v+shh+238)/37.5) ) )
	} else {	
	        tau_h = 70 + hx
	}
        if (v < -73-shn) {
                tau_n = ( 1.0 / ( exp((v+shn+46)/5) + exp(-(v+shn+238)/37.5) ) ) 
        } else {
                tau_n = 60 + nx
        }
}
UNITSON




