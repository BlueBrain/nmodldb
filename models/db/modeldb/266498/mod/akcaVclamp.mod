: KCa is the calcium-activated potassium current in Schild 1994 

NEURON {
       SUFFIX akca
       USEION k READ ek WRITE ik
	   USEION ca READ cai
       RANGE gbar, ek, ik

}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
	  (molar) = (/liter)
	  (mM) = (millimolar)
}

PARAMETER {
    gbar =0.00022989 (S/cm2)
	Q10kcac=2.30

	
	A_alphac=750.0 (/ms-mM) :From Schild 1994, alphac=A_alphac*cai*((Vm+B_alphan)/C_alphan)
	B_alphac=-10.0 (mV)
	C_alphac=12.0	(mV)
	
	A_betac=0.05 (/ms)	:From Schild 1994, betac=A_betan*exp((Vm+B_betan)/C_betan)
	B_betac=-10.0 (mV)
	C_betac=-60.0 (mV)
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
	 cai (mM)
	 celsius (degC)
	 g	(S/cm2)
	 tau_c	(ms)
	 cinf
	 alphac (/ms)
	 betac (/ms)
     ek	(mV)
         
}

STATE { c }

BREAKPOINT {
	   SOLVE states METHOD cnexp
	   g = gbar * c
	   ik = g * (v-ek)
}

INITIAL {
	rates(v) : set tau_c, cinf
	: assume that equilibrium has been reached
	
    c = cinf
}

DERIVATIVE states {
	   rates(v)
	   c' = (cinf - c)/tau_c
}

FUNCTION alpha(Vm (mV)) (/ms) {
	alphac=A_alphac*cai*exp((Vm+B_alphac)/C_alphac)
}

FUNCTION beta(Vm (mV)) (/ms) {
	betac=A_betac*exp((Vm+B_betac)/C_betac)
}

FUNCTION rates(Vm (mV)) (/ms) {
	alpha(Vm)
	beta(Vm)
	
	tau_c = 4.5/(alphac+betac)
	cinf = alphac/(alphac+betac)
	
	if (celsius >= 37) {
	tau_c=tau_c*(1/Q10kcac)
	}
}