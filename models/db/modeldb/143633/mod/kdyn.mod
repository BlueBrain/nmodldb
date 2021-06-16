: simple first-order model of potassium dynamics
:
: FROM:
: Durstewitz D, Seamans JK, Sejnowski TJ (2000) Dopamine-mediated
: stabilization of delay-period activity in a network model of
: prefrontal cortex. J Neurophysiol 83:1733-50


NEURON {
	SUFFIX kdyn
	USEION k READ ik WRITE ko, ki 
	RANGE ko, ki, tk, tk0, dep, tkfac
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F    = (faraday) (coul)
}

PARAMETER {
	tk0    = 1342.4  (ms)           : decay time constant
	tkfac = 0.025  : factor for multiplying tauk equation
	koinf = 3 	(mM)      : equilibrium k+ concentration
	dep   = 290e-3 (micron)     : depth of shell for k+ diffusion
	KAF   = 1 ()		  : K accumulation factor
}

ASSIGNED {
	ik     (mA/cm2)
	tk	   (ms)
	ki     (mM)	  :
}

STATE { ko (mM) }

BREAKPOINT { 
	SOLVE states METHOD cnexp
}

DERIVATIVE states {
    evaluate_tau()
        ko'= (1e4)*KAF*ik/(F*dep) + (koinf-ko)/tk
}

INITIAL {
	tk = tk0
	ko = koinf
	ki = 106
}

UNITSOFF
PROCEDURE evaluate_tau() { LOCAL tauk 
	tauk = (tkfac*((ko-koinf)^2)) + (tk0/1000) 	: Removed (-) piece to prevent (-) Taus. Based on regression of Cordingley & Sonjen, 1978 Fig. 3b
												: Changed from 0.025*... to 0.015*...
	tk = tauk*1000 : convert to ms from sec
	ki = 106 
}
UNITSON
