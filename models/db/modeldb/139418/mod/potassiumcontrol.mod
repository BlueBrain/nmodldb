TITLE  Glial Buffering, and Potassium Diffusion Mechanism

UNITS {
      (molar) = (1/liter)
      (mV) = (millivolts)
      (mA) = (milliamp)
      (mM) = (millimolar)
      FARADAY = (faraday) (coulombs)
}

NEURON {
       SUFFIX potassium
       USEION k READ ik, ko  WRITE ko
       RANGE fhspace, txfer, bmax, k1N, k1,  kbath
}

PARAMETER {
	  kbath = 9			  (mM)
	  fhspace =300			  (angstrom) : thickness of F-H space
	  txfer = 50			  (ms) : tau for F-H space
	  bmax = 50			  (mM)
	  k1N =1.1			  (none)
	  k1 = 0.0008			  (none)
}

ASSIGNED {
	 ik	 (mA/cm2)
	 g	 (mM)
	 k2	 (none)
}

STATE {
        ko (mM)
	b  (mM)
 }

BREAKPOINT {
	   SOLVE states METHOD cnexp
}

DERIVATIVE states {
	   rates (ko)
	   b' = k1*(bmax - b) - k2*ko*b 
	   glial (b)
	   ko' = (10*ik)/(0.15*FARADAY) + (kbath - ko)/txfer + g
}

PROCEDURE rates (ko) {
	   k2 = k1/( 1+exp ((ko-kbath)/(-1.15)) )
}

PROCEDURE glial (b) {
       	   g = (k1*(bmax-b)/k1N) - k2*ko*b 
}
