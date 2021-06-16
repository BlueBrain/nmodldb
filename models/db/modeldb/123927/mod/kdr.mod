NEURON {
  SUFFIX kdrG
  USEION k READ ek WRITE ik
  RANGE gbar, g, i
}

UNITS {
  (S) = (siemens)
  (mV) = (millivolt)
  (mA) = (milliamp)	
}

PARAMETER {
  gbar = 0.006 (S/cm2)
	eK=  -95 (mV)
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
}
  
INITIAL {
  n = ninf(v)
}

DERIVATIVE states {
 n'= (ninf(v)-n)/ntau(v)
}

FUNCTION ninf (Vm (mV)) () {

  UNITSOFF
    ninf = 1/(1+exp(-(Vm+35)/10))
  UNITSON

}

FUNCTION ntau (Vm (mV)) (ms) {

  UNITSOFF
    ntau = 0.1 + (0.5 / (1+exp((Vm+27)/15)))
  UNITSON

}

