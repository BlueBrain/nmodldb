: K-leakage current
NEURON {
  SUFFIX kleck
  USEION k READ ek WRITE ik
  RANGE gbar, g, i
}

UNITS {
  (S) = (siemens)
  (mV) = (millivolt)
  (mA) = (milliamp)	
}

PARAMETER {
  gbar = 0.001 (S/cm2) 
  eK = -95 (mV)	
}

ASSIGNED {
  v	(mV)
  ek	(mV)
  ik 	(mA/cm2)
  i 	(mA/cm2)
  g	(S/cm2)
  
}

BREAKPOINT {
  i = gbar*(v-eK)
  ik = i
}

