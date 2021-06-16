: Kds is the slowly inactivating delay current in Schild 1994 (Part of Id transient current)

NEURON {
       SUFFIX akds
       USEION k READ ek WRITE ik
       RANGE gbar, ek, ik

}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
}

PARAMETER {
    gbar =0.000406729  (S/cm2): (S/cm2)
	Q10kds=1.93 :All gating variables have the same constant
	
	V0p5x=-39.59 (mV):As defined by Schild 1994, zinf=1.0/(1.0+exp((V0p5z-V)/S0p5z)
	S0p5x=14.68(mV)
	
	V0p5y=-48.0 (mV)
	S0p5y=-7.0 (mV)
	
	
	A_taux=5.0	(ms)	:As defined by Schild 1994, tauz=A_tauz*exp(-B^2(V-Vpz)^2)+C
	B_taux=0.022	(/mV)
	C_taux=2.5	(ms)
	Vpx=-65.0		(mV)
	
	tau_y22=7500 (ms) :This is tau_y at 22 degC

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
	 tau_x	(ms)
	 tau_y	(ms)
	 xinf
	 yinf
     ek	(mV)
         
}

STATE { x y1 } 


BREAKPOINT {
	   SOLVE states METHOD cnexp
	   g = gbar * x^3 * y1
	   ik = g * (v-ek)
}

INITIAL {
	rates(v) : set tau_x, yinf, xinf
	: assume that equilibrium has been reached
	

    x = xinf
	y1 = yinf
}

DERIVATIVE states {
	   rates(v)
	   x' = (xinf - x)/tau_x
	   y1' = (yinf - y1)/tau_y
}



FUNCTION rates(Vm (mV)) (/ms) {
	tau_x = A_taux*exp(-(B_taux)^2*(Vm-Vpx)^2)+C_taux
        xinf = 1.0/(1.0+exp((Vm-V0p5x)/(-S0p5x)))
		 
    tau_y=tau_y22
		yinf = 1.0/(1.0+exp((Vm-V0p5y)/(-S0p5y)))
	
	if (celsius >= 37) {
	tau_x=tau_x*(1/Q10kds)
	tau_y=tau_y*(1/Q10kds)
	}
}