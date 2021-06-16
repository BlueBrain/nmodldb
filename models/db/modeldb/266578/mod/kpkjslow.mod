: HH Slow TEA-insensitive Purkinje potassium current
: FORREST MD (2014) Two Compartment Model of the Cerebellar Purkinje Neuron

NEURON {
	SUFFIX kpkjslow
	USEION k READ ek WRITE ik
	RANGE gkbar
	GLOBAL ninf, ntau
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER {
	v		(mV)
	gkbar = .004	(mho/cm2)
	
	nivh = -16.5	(mV)
	nik = 18.4
	
	ek
Q10 = 3 (1) 
  Q10TEMP = 22 (degC) 

}

ASSIGNED {
	ik
	ninf
	ntau		(ms)
celsius (degC) 
  qt (1) 

}

STATE {
	n
}

INITIAL {
	rates(v)
	n = ninf
: qt = Q10^((celsius-Q10TEMP)/10) 
qt = 1
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkbar * n^4 * (v - ek)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n) / ntau
}

PROCEDURE rates(Vm (mV)) {
	LOCAL v
	v = Vm + 11	: Account for Junction Potential
	ninf = 1/(1+exp(-(v-nivh)/nik))
	ntau = (1000 * ntau_func(v)) / qt
}

FUNCTION ntau_func(v (mV)) {
	ntau_func = .000796 + 1/(exp((v+73.2)/11.7)+exp((v-306.7)/-74.2))
}