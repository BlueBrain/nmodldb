: Kd is the delayed rectifier current in Schild 1994 

NEURON {
       SUFFIX akd
       USEION k READ ek WRITE ik
       RANGE gbar, ek, ik

}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
}

PARAMETER {
    gbar =0.001061033 (S/cm2)
	Q10kdn=1.40

	
	V0p5n=-14.62 (mV):As defined by Schild 1994, zinf=1.0/(1.0+exp((V0p5z-V)/S0p5z)
	S0p5n=18.38 (mV)
	
	A_alphan=.001265 (/ms-mV) :From Schild 1994, alphan=A_alphan*(Vm+B_alphan)/(1.0-exp((Vm+B_alphan)/C_alphan)
	B_alphan=14.273 (mV)
	C_alphan=-10.0	(mV)
	
	A_betan=0.125 (/ms)	:From Schild 1994, betan=A_betan*exp((Vm+B_betan)/C_betan)
	B_betan=55.0 (mV)
	C_betan=-2.5 (mV)
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
	 celsius (degC)
	 g	(S/cm2)
	 tau_n	(ms)
	 ninf
	 alphan (/ms)
	 betan (/ms)
     ek	(mV)
         
}

STATE { n }

BREAKPOINT {
	   SOLVE states METHOD cnexp
	   g = gbar * n
	   ik = g * (v-ek)
}

INITIAL {
	rates(v) : set tau_m, tau_h, hinf, minf
	: assume that equilibrium has been reached
	
    n = ninf
}

DERIVATIVE states {
	   rates(v)
	   n' = (ninf - n)/tau_n
}

FUNCTION alpha(Vm (mV)) (/ms) {
	alphan=(A_alphan*(Vm+B_alphan))/(1.0-exp((Vm+B_alphan)/C_alphan))
}

FUNCTION beta(Vm (mV)) (/ms) {
	betan=A_betan*exp((Vm+B_betan)/C_betan)
}

FUNCTION rates(Vm (mV)) (/ms) {
	alpha(Vm)
	beta(Vm)
	
	tau_n = 1/(alphan+betan)+1.0
	ninf = 1.0/(1.0+exp((Vm-V0p5n)/(-S0p5n)))
	
	if (celsius >= 37) {
	tau_n=tau_n*(1/Q10kdn)
	}
}
