
TITLE KATP channel 
 
: Author: Chitaranjan Mahapatra (chitaranjan@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: Mahapatra C, Brain KL, Manchanda R, A biophysically constrained computational model of the action potential 
: of mouse urinary bladder smooth muscle. PLOS One (2018) 
 
NEURON {
	SUFFIX KATP
	USEION atp READ atpi VALENCE -1
	USEION k READ ek WRITE ik
	RANGE gbar, ik, tauo,oinf
	GLOBAL qdeltat,at,athf,ap, q10
}
 
UNITS {
	(mA) = (milliamp) 
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}
PARAMETER { 
    qdeltat = 1
    gbar = 0.001(mho/cm2)
	atpi = 0.0003 (mM)
	athf = 0.006 (mM)
	ek = -21 (mV)
	at = 10
	ap =3
	celsius	= 37	(degC)
	        dt (ms)
			q10=2.3
} 

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	oinf 
    tauo (ms) 
} 
 
STATE{
	o 
} 

INITIAL{ 
    rate(atpi)
    o = oinf 
} 
 
BREAKPOINT { 
    SOLVE states METHOD cnexp 
	ik = gbar * o * (v - ek)
	
} 

DERIVATIVE states { 
	rate(atpi) 
	o' = (o-oinf ) / tauo
} 

		

PROCEDURE rate(atpi(mM)) {
     LOCAL a,qt
        qt=q10^((celsius-37)/10)
    oinf = (1/ (1+ ((atpi/athf)^ap)))

   if (atpi < 0.005) {
       tauo = 1 - (186.67 * atpi)/qt
    } else {
        tauo = 2/qt
    }
    tauo = at *  (tauo / qdeltat)
} 

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	       trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	