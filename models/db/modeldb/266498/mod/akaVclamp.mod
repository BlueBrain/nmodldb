: Ka is the Early Transient Outward K current in Schild 1994 

NEURON {
       SUFFIX aka
       USEION k READ ek WRITE ik
       RANGE gbar, ek, ik

}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
}

PARAMETER {
    gbar =0.000990297 (S/cm2): (S/cm2)
	Q10ka=1.93 :All gating variables have the same constant
	
	V0p5p=-28.0 (mV):As defined by Schild 1994, zinf=1.0/(1.0+exp((V0p5z-V)/S0p5z)
	S0p5p=28.0 (mV)
	
	V0p5q=-58.0 (mV)
	S0p5q=-7.0 (mV)
	
	
	A_taup=5.0	(ms)	:As defined by Schild 1994, tauz=A_tauz*exp(-B^2(V-Vpz)^2)+C
	B_taup=0.022	(/mV)
	C_taup=2.5	(ms)
	Vpp=-65.0		(mV)
	
	A_tauq=100.0	(ms)
	B_tauq=0.035	(/mV)
	C_tauq=10.5	(ms)
	Vpq=-30.0	(mV)

}
COMMENT
	The above Q10 constants were given in Schild 1994 with no indication of how they
	were implemented. It was decided, based on the value of the Q10 constants given,
	that the most likely answer was that the tau of each gating variable was divided by
	the Q10. This is reflected below where tau_x=tau_x*(1/Q10x). Note that Schild only
	provides a single constant, not any type of equation for q10. The equations are orginally
	given for 22C, and this constant changes the equation to 37C.
ENDCOMMENT

ASSIGNED {
	 v	(mV) : NEURON provides this
	 ik	(mA/cm2)
	 celsius  (degC)
	 g	(S/cm2)
	 tau_q	(ms)
	 tau_p	(ms)
	 pinf
	 qinf
     ek	(mV)
         
}

STATE { p q } 


BREAKPOINT {
	   SOLVE states METHOD cnexp
	   g = gbar * p^3 * q
	   ik = g * (v-ek)
}

INITIAL {
	rates(v) : set tau_m, tau_h, hinf, minf
	: assume that equilibrium has been reached
	

    p = pinf
	q = qinf
}

DERIVATIVE states {
	   rates(v)
	   p' = (pinf - p)/tau_p
	   q' = (qinf - q)/tau_q
}



FUNCTION rates(Vm (mV)) (/ms) {
	 tau_p = A_taup*exp(-(B_taup)^2*(Vm-Vpp)^2)+C_taup
         pinf = 1.0/(1.0+exp((Vm-V0p5p)/(-S0p5p)))

	 tau_q = A_tauq*exp(-(B_tauq)^2*(Vm-Vpq)^2)+C_tauq
         qinf = 1.0/(1.0+exp((Vm-V0p5q)/(-S0p5q)))
	
	if (celsius >= 37) {
	tau_p=tau_p*(1/Q10ka)
	tau_q=tau_q*(1/Q10ka)
	}
}