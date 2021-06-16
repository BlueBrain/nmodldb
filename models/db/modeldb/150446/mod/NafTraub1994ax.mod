TITLE Fast sodium current, for axonal compartments
 
COMMENT
  from Table 4 of "A branching dendritic model of a rodent CA3 pyramidal neurone." Traub RD et al. J Physiol. (1994) 
  implemented by Nikita Vladimirov <nikita.vladimirov@gmail.com>
ENDCOMMENT

NEURON {
        SUFFIX NaFax
		USEION na READ ena WRITE ina
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
		ena (mV)
		ina (mA/cm2)
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
		g = gbar * m^3 * h
		i = g * (v - ena)
		ina = i
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
			small = (17.2 - (v - Vm) )/4
			if ( fabs(small) > 1e-6 ) {
				alpham =  0.8 * (17.2 - (v - Vm) ) / ( exp( (17.2 - (v - Vm) )/4 ) - 1 )
			} else {
				alpham = 0.8 * 4 / (1 + small/2)
			}
			small = ( (v - Vm) - 42.2)/5 
			if ( fabs(small) > 1e-6 ) {
				betam  =  0.7 * ( (v - Vm) - 42.2) / ( exp( ( (v - Vm) - 42.2)/5 ) - 1 ) 			
			} else {
				betam  =  0.7 * 5 / ( 1 + small/2 )
			}
			minf   = alpham / ( alpham + betam )
			mtau   = 1 / ( alpham + betam )
			alphah = 0.32 * exp( (42 - (v - Vm) ) / 18 )
			betah  = 10 / ( 1 + exp( (42 - (v - Vm) )/5 ) )
			hinf   = alphah / ( alphah + betah )
			htau   = 1 / ( alphah + betah )
		UNITSON
}