: Naf is the fast, TTX-sensitive current in Schild 1994 

NEURON {
       SUFFIX anaf
       USEION na READ ena WRITE ina
       RANGE gbar, ena, ina

}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
}

PARAMETER {
    gbar =0.001043349  (S/cm2): (S/cm2)
	Q10nafm=2.30
	Q10nafh=1.50
	
	V0p5m=-41.35 (mV):As defined by Schild 1994, zinf=1.0/(1.0+exp((V0p5z-V)/S0p5z)
	S0p5m=4.75 (mV)
	
	V0p5h=-62.00 (mV)
	S0p5h=-4.50 	(mV)
	
	V0p5j=-40.00 (mV)
	S0p5j=-1.50 	(mV)
	
	A_taum=0.75	(ms)	:As defined by Schild 1994, tauz=A_tauz*exp(-B^2(V-Vpz)^2)+C
	B_taum=0.0635	(/mV)
	C_taum=0.12	(ms)
	Vpm=-40.35	(mV)
	
	A_tauh=6.5	(ms)
	B_tauh=0.0295	(/mV)
	C_tauh=0.55	(ms)
	Vph=-75.00	(mV)
	
	A_tauj=25	(ms)
	B_tauj=4.50	(mV)
	C_tauj=0.01	(ms)
	Vpj=-20.00	(mV)

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
	 celsius (degC)
	 ina	(mA/cm2)
	 g	(S/cm2)
	 tau_h	(ms)
	 tau_m	(ms)
	 tau_j	(ms)
	 minf
	 hinf
	 jinf
     ena	(mV)
         
}

STATE { m h l } 

COMMENT
	:#####NOTE##### A compiling error was encountered,
	and the problem seemed to be that the variable j was
	used by the compilor for some other process. The error
	was "the variable j0 was declared in another state" or
	something to that effect. To test this theory, I tried 
	to change the gating variable to l in this script.  This
	allowed the file to compile with no errors.  So for now,
	note that l in this file refers to the reactivation
	variable j in Schild 1994.  The other variables 
	(e.g jinf, tau_j) used to define j (now l) were not changed.
ENDCOMMENT

BREAKPOINT {
	   SOLVE states METHOD cnexp
	   g = gbar * m^3 * h *l
	   ina = g * (v-ena)
}

INITIAL {
	rates(v) : set tau_m, tau_h, hinf, minf
	: assume that equilibrium has been reached
	

    m = minf
	h = hinf
	l = jinf
}

DERIVATIVE states {
	   rates(v)
	   m' = (minf - m)/tau_m
	   h' = (hinf - h)/tau_h
	   l' = (jinf - l)/tau_j
}



FUNCTION rates(Vm (mV)) (/ms) {
	 tau_m = A_taum*exp(-(B_taum)^2*(Vm-Vpm)^2)+C_taum
         minf = 1.0/(1.0+exp((Vm-V0p5m)/(-S0p5m)))

	 tau_h = A_tauh*exp(-(B_tauh)^2*(Vm-Vph)^2)+C_tauh
         hinf = 1.0/(1.0+exp((Vm-V0p5h)/(-S0p5h)))

	 tau_j = (A_tauj/(1.0+exp((Vm+Vpj)/B_tauj)))+C_tauj
         jinf = 1.0/(1.0+exp((Vm-V0p5j)/(-S0p5j)))
		 
	if (celsius >= 37) {
		tau_m=tau_m*(1/Q10nafm)
		tau_h=tau_h*(1/Q10nafh)
	}

}