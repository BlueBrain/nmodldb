: Nas is the slower, TTX-insensitive current in Schild 1994 

NEURON {
       SUFFIX anas
       USEION na READ ena WRITE ina
       RANGE gbar, ena, ina

}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
}

PARAMETER {
    gbar =0.000901878 (S/cm2): (S/cm2)
	Q10nasm=2.30
	Q10nash=1.50
	
	V0p5m=-20.35 (mV):As defined by Schild 1994, zinf=1.0/(1.0+exp((V0p5z-V)/S0p5z)
	S0p5m=4.45 (mV)
	
	V0p5h=-18.00 (mV)
	S0p5h=-4.50 (mV)
	
	
	A_taum=1.50	(ms)	:As defined by Schild 1994, tauz=A_tauz*exp(-B^2(V-Vpz)^2)+C
	B_taum=0.0595	(/mV)
	C_taum=0.15	(ms)
	Vpm=-20.35		(mV)
	
	A_tauh=4.95	(ms)
	B_tauh=0.0335	(/mV)
	C_tauh=0.75	(ms)
	Vph=-20.00	(mV)

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
	 ina	(mA/cm2)
	 celsius (degC)
	 g	(S/cm2)
	 tau_h	(ms)
	 tau_m	(ms)
	 minf
	 hinf
     ena	(mV)
         
}

STATE { m h } 


BREAKPOINT {
	   SOLVE states METHOD cnexp
	   g = gbar * m^3 * h
	   ina = g * (v-ena)
}

INITIAL {
	rates(v) : set tau_m, tau_h, hinf, minf
	: assume that equilibrium has been reached
	

    m = minf
	h = hinf
}

DERIVATIVE states {
	   rates(v)
	   m' = (minf - m)/tau_m
	   h' = (hinf - h)/tau_h
}



FUNCTION rates(Vm (mV)) (/ms) {
	 tau_m = A_taum*exp(-(B_taum)^2*(Vm-Vpm)^2)+C_taum
         minf = 1.0/(1.0+exp((Vm-V0p5m)/(-S0p5m)))

	 tau_h = A_tauh*exp(-(B_tauh)^2*(Vm-Vph)^2)+C_tauh
         hinf = 1.0/(1.0+exp((Vm-V0p5h)/(-S0p5h)))
		 
	if (celsius >= 37) {
	tau_m=tau_m*(1/Q10nasm)
	tau_h=tau_h*(1/Q10nash)
	}
}