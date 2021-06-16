TITLE Potassium delayed rectifer

NEURON {
	SUFFIX kdr
  	USEION k READ ek WRITE ik
  	RANGE gbar, g, i, jikdr
}

UNITS {
  	(S) = (siemens)
  	(mV) = (millivolt)
  	(mA) = (milliamp)	
}

PARAMETER {
    gbar = 1e-2 (S/cm2)
    nimid = 0
    nislope = 25
    ntmid = 27
    ntslope = 15
    eK=  -95 (mV)
    jikdr 	(mA/cm2)
}

ASSIGNED {
  	v	(mV)
  	ek	(mV)
  	ik 	(mA/cm2)
  	i 	(mA/cm2)
  	g	(S/cm2)
}

STATE {n}

BREAKPOINT {
  	SOLVE states METHOD cnexp
  	g = gbar*n*n*n*n
  	i = g*(v-eK)
  	ik = i
	jikdr = i
}
  
INITIAL {
  	n = ninf(v)
}

DERIVATIVE states {
 	n'= (ninf(v)-n)/ntau(v)
}

FUNCTION ninf (Vm (mV)) () {

  	UNITSOFF
    	ninf = 1/(1+exp(-(Vm+nimid)/nislope))
  	UNITSON
}

FUNCTION ntau (Vm (mV)) (ms) {

  	UNITSOFF
        if( v < -10.0 ) {
                ntau = 0.25 + 4.35 * exp( ( Vm + ntmid ) / ntslope )
        }else{
                ntau = 0.25 + 4.35 * exp( ( -Vm - ntmid ) / ntslope )
        }
  	UNITSON
}

