TITLE Ca-dependent potassium current, K(C)-current
 
COMMENT
  from Table 3 of "A branching dendritic model of a rodent CA3 pyramidal neurone." Traub RD et al. J Physiol. (1994) 
  implemented by Nikita Vladimirov <nikita.vladimirov@gmail.com>
ENDCOMMENT

NEURON {
        SUFFIX Kc
		USEION k READ ek WRITE ik
		USEION ca READ cai
        RANGE  gbar, g, i
		GLOBAL Vm
} 
 
UNITS {
		(S)  = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
		(molar) = (1/liter)			: moles do not appear in units
		(mM)	= (millimolar)
}

PARAMETER { 
		gbar = 0   (S/cm2) 
		Vm   = -65 (mV) : resting potential
}

ASSIGNED {
		v   (mV)
		ek  (mV)
		cai	(mM)
		ik  (mA/cm2)
		i   (mA/cm2)
		g   (S/cm2)
		minf
		mtau (ms) 
}

STATE { m }

BREAKPOINT {
		SOLVE states METHOD cnexp
		if( cai /( 250 (mM) ) < 1 ) {
			g = gbar * m * cai/( 250 (mM) )
		} else  {
			g = gbar * m
		}
		i = g * (v - ek)
		ik = i
}

INITIAL {
		rates(v)
		m = minf
}

DERIVATIVE states {
        rates(v)
        m' = (minf - m) / mtau
}

PROCEDURE rates( v(mV) ) {
		LOCAL  alpham, betam
        TABLE minf, mtau FROM -100 TO 50 WITH 200
		UNITSOFF
		if( v - Vm <= 50 ) {
			alpham = exp( (v-Vm-10) / 11 - (v - Vm - 6.5 ) /27 ) / 18.975
			betam = 2 * exp( - (v - Vm - 6.5) / 27) - alpham
		} else {
			alpham = 2 * exp( - (v - Vm - 6.5) / 27)
			betam  = 0
		}
		minf   = alpham / ( alpham + betam )
		mtau   = 1 / ( alpham + betam )
		UNITSON
}