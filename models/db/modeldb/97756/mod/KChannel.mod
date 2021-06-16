: Two state kinetic scheme for Potassium channel
: Contains four kinetic parameters and one max conductance parameter.
NEURON {
      SUFFIX KCHANNEL
      USEION k READ ek WRITE ik
      RANGE g, gbar,a12,a21,z12,z21
}
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
      gbar = 0     (pS/um2)
      a12 = 0.01   (/ms)
      a21 = 0.02   (/ms)
      z12 = 0.01   (/mV)
      z21 = 0.02   (/mV)
}

ASSIGNED {
      v    (mV)
      ek   (mV)
      g    (pS/um2)
      ik   (mA/cm2)
      k12  (/ms)
      k21  (/ms)
}

STATE { c o }

BREAKPOINT {
      SOLVE states METHOD sparse
      g = gbar*o
      ik = (1e-4)*g*(v - ek)
}

INITIAL { SOLVE states STEADYSTATE sparse}

KINETIC states {   		
        rates(v)
	~c <-> o (k12,k21)
	CONSERVE c+o=1
}

PROCEDURE rates(v(millivolt)) {

      k12 = a12*exp(z12*v)
      k21 = a21*exp(-z21*v)
}
