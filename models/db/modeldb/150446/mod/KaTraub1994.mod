TITLE Transient K current (A-current)
 
COMMENT
  from Table 3 of "A branching dendritic model of a rodent CA3 pyramidal neurone." Traub RD et al. J Physiol. (1994) 
  implemented by Nikita Vladimirov <nikita.vladimirov@gmail.com>
ENDCOMMENT

NEURON {
        SUFFIX Ka
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
		minf
		hinf
		mtau (ms) 
		htau (ms)
}

STATE { m h }

BREAKPOINT {
		SOLVE states METHOD cnexp
		g = gbar * m * h
		i = g * (v - ek)
		ik = i
}

INITIAL {
		rates(v)
		m = minf
		h = hinf
}

DERIVATIVE states {
        rates(v)
        m' = (minf - m) / mtau
        h' = (hinf - h) / htau
}

PROCEDURE rates(v(mV)) {
		LOCAL  alpham, betam, alphah, betah, small
        TABLE minf, mtau, hinf, htau FROM -100 TO 50 WITH 200
		UNITSOFF
			small = (13.1 - (v - Vm) )/10 
			if ( fabs(small) > 1e-6 ) {
				alpham =  0.02 * (13.1 - (v - Vm) ) / ( exp( (13.1 - (v - Vm) )/10 ) - 1 )
			} else {
				alpham =  0.02 * 10 / ( 1 + small/2)
			}
			small = ( (v - Vm) - 40.1)/10
			if ( fabs(small) > 1e-6 ) {
				betam  =  0.0175 * ( (v - Vm) - 40.1) / ( exp( ( (v - Vm) - 40.1)/10 ) - 1 ) 
			} else {
				betam  =  0.0175 * 10 / ( 1 + small/2 )
			}
			minf   = alpham / ( alpham + betam )
			mtau   = 1 / ( alpham + betam )
			alphah = 0.0016 * exp( (-13 - (v - Vm) ) / 18 )
			betah  = 0.05 / ( 1 + exp( (10.1 - (v - Vm) )/5 ) )
			hinf   = alphah / ( alphah + betah )
			htau   = 1 / ( alphah + betah )
		UNITSON
}