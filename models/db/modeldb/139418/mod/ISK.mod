TITLE Calcium Dependant Potassium Channel

UNITS {
      (molar) = (1/liter)
      (mV) = (millivolts)
      (mA) = (milliamp)
      (mM) = (millimolar)
	(S) = (siemens)
}

NEURON {
       SUFFIX sk
       USEION ca READ cai
       USEION k READ ek WRITE ik
       RANGE i, gbar, theta, ztau
       GLOBAL zinf
}

PARAMETER {
	  cai	(mM)
	  ek	(mV)
	  v	(mV)
	  gbar = 0.0175	(S/cm2)
	  ztau = 1	(ms)
	  i   (mA/cm2)
}

ASSIGNED {
	 ik	 (mA/cm2)
	 zinf	 
}

STATE {z}

BREAKPOINT {
	   SOLVE state METHOD cnexp
	   i = gbar*(z^2)*(v-ek)
	   ik = i
}

DERIVATIVE state {
	   rate(cai) 
	   z' = (zinf-z)/ztau
}

INITIAL {
	rate(cai)
	z = zinf
}

FUNCTION alp(ca (mM)) (/ms) {
	 alp = 1/(1+((0.003/ca)^2))
}

PROCEDURE rate(ca (mM)) {
	  zinf = alp(ca)
}
