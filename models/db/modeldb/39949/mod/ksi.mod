COMMENT
Described by Gruber et al. 2003, which they based on Nisenbaum et al. 1996

Gruber, A.J., Solla, S.A., Surmeier, D.J., and Houk, J.C.
Modulation of striatal single units by expected reward: 
a spiny neuron model displaying dopamine-induced bistability.
J. Neurophysiol. 90:1095-1114, 2003.

Nisenbaum, E.S., Wilson, C.J., Foehring, R.C., and Surmeier, D.J.
Isolation and characterization of a persistent potassium current in neostriatal neurons.
J. Neurophysiol. 76:1180-1194, 1996.

Unlike the formulation used by Gruber et al., 
which assumed instantaneous activation, 
this implementation assumes a constant activation time constant 
that is relatively fast compared to the time scale of the model 
(100-1000 ms).
ENDCOMMENT

NEURON {
	SUFFIX ksi
	USEION k READ ek WRITE ik
	RANGE gbar, g, i
	GLOBAL ninf, ntau
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER {
	gbar = 0.45	(mS/cm2)	<0,1e9>
	ek = -90	(mV)
	vh = -13.5	(mV)	: half activation
	ve = 11.8	(mV)	: slope
	ntauconst = 0.1	(ms)	: n activates much faster than 100-1000 ms
}

ASSIGNED {
	v	(mV)
	g	(mho/cm2)
	i	(uA/cm2)	: for consistency with their usage of uA/cm2
	ik	(mA/cm2)
	ninf	(1)
	ntau	(ms)
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = (0.001)*gbar*n
	ik = g*(v - ek)     
	i = (1000)*ik
}

INITIAL {
	rates(v)
	n = ninf
}

DERIVATIVE states { 
	rates(v)
	n' = (ninf-n)/ntau
}


: rates() computes rate and other constants at present v
: call once from hoc to initialize inf at resting v
PROCEDURE rates(v(mV)) {
UNITSOFF
	: "n" potassium activation
	ntau = ntauconst
	ninf = 1/(1 + exp(-(v - vh)/ve))
}
UNITSON
