TITLE Delayed rectifier K current, for soma-dendrites compartments
 
COMMENT
  from Table 3 of "A branching dendritic model of a rodent CA3 pyramidal neurone." Traub RD et al. J Physiol. (1994) 
  implemented by Nikita Vladimirov <nikita.vladimirov@gmail.com>
ENDCOMMENT

NEURON {
        SUFFIX Kdrsd
		USEION k READ ek WRITE ik
        RANGE  gbar, g, i
		GLOBAL Vm
} 
 
UNITS {
		(S)  = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
}

PARAMETER { 
		gbar = 0   (S/cm2)
		Vm   = -65 (mV) : resting potential
}

ASSIGNED {
		v   (mV)
		ek  (mV)
		ik  (mA/cm2)
		i   (mA/cm2)
		g   (S/cm2)
		ninf
		ntau  (ms)
}

STATE { n }

BREAKPOINT {
		SOLVE states METHOD cnexp
		g = gbar * n^2
		i = g * (v - ek)
		ik = i
}

INITIAL {
		rates(v)
		n = ninf
}

DERIVATIVE states {
        rates(v)
        n' = (ninf-n)/ntau
}

PROCEDURE rates(v(mV)) {
		LOCAL  alphan, betan, small
        TABLE ninf, ntau FROM -100 TO 50 WITH 200
		UNITSOFF
			small = (35.1 - (v - Vm) )/5 
			if ( fabs(small) > 1e-6 ) {
				alphan = 0.016 * ( 35.1 - (v - Vm) ) / ( exp( (35.1 - (v - Vm) )/5 ) - 1)
			} else {
				alphan = 0.016 * 5 / (1 + small/2)
			}
			betan  = 0.25 * exp( (20 - (v - Vm) ) / 40 )
			ninf   = alphan / (alphan + betan)
			ntau   = 1 / (alphan + betan)
		UNITSON
}