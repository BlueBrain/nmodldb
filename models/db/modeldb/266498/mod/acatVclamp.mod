: Cat is the Low threshold, transient Ca current in Schild 1994 

NEURON {
       SUFFIX acat
       USEION ca READ cai, cao WRITE ica
       RANGE gbar, ecat, ica

}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
	  F = (faraday) (coulombs)
	  (molar) = (1/liter)
	  (mM)    = (millimolar)
}

PARAMETER {
    gbar =9.72614E-05 (S/cm2): (S/cm2)
	Q10catd=1.90
	Q10catf=2.20
	
	V0p5d=-54.00 (mV):As defined by Schild 1994, zinf=1.0/(1.0+exp((V0p5z-V)/S0p5z)
	S0p5d=5.75 (mV)
	
	V0p5f=-68.00 (mV)
	S0p5f=-6 (mV)
	
	
	A_taud=22.0	(ms)	:As defined by Schild 1994, tauz=A_tauz*exp(-B^2(V-Vpz)^2)+C
	B_taud=0.052	(/mV)
	C_taud=2.5	(ms)
	Vpd=-68.0		(mV)
	
	A_tauf=103.0	(ms)
	B_tauf=0.050	(/mV)
	C_tauf=12.5	(ms)
	Vpf=-58.0	(mV)
	
	R=8.314 (joule/degC): Gas Constant
	z=2 : Charge of Ca ion
	ecaoffset=78.7 (mV)
	

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
	 ica	(mA/cm2)
	 celsius (degC)
	 g	(S/cm2)
	 tau_f	(ms)
	 tau_d	(ms)
	 dinf
	 finf
     ecat	(mV)
	 cai (mM)
	 cao (mM)
         
}

STATE { d f } 


BREAKPOINT {
	   SOLVE states METHOD cnexp
	   g = gbar * d * f
	   ica = g * (v-ecat)
}

INITIAL {
	rates(v) : set tau_m, tau_h, hinf, minf
	: assume that equilibrium has been reached
	

    d = dinf
	f = finf
}

DERIVATIVE states {
	   rates(v)
	   d' = (dinf - d)/tau_d
	   f' = (finf - f)/tau_f
}



FUNCTION rates(Vm (mV)) (/ms) {
	 tau_d = A_taud*exp(-(B_taud)^2*(Vm-Vpd)^2)+C_taud
         dinf = 1.0/(1.0+exp((Vm-V0p5d)/(-S0p5d)))

	 tau_f = A_tauf*exp(-(B_tauf)^2*(Vm-Vpf)^2)+C_tauf
         finf = 1.0/(1.0+exp((Vm-V0p5f)/(-S0p5f)))
		 
	ecat=(1000)*(R*(celsius+273.15)/z/F*log(cao/cai))-ecaoffset : Equation for eca given in Schild 1994
		 
	if (celsius >= 37) {
		tau_d=tau_d*(1/Q10catd)
		tau_f=tau_f*(1/Q10catf)
	}
}