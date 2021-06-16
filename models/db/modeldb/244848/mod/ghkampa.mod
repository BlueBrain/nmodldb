COMMENT

The kinetics part is obtained from Exp2Syn of NEURON.

Two state kinetic scheme synapse described by rise time taur, and
decay time constant taud. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/taur and 1/taud is

 A = a*exp(-t/taur) and
 G = a*taud/(taud-taur)*(-exp(-t/taur) + exp(-t/taud))
	where taur < taud

If taud-taur -> 0 then we have a alphasynapse.
and if taur -> 0 then we have just single exponential decay.

The factor is evaluated in the initial block such that an event of
weight 1 generates a peak conductance of 1.

Because the solution is a sum of exponentials, the coupled equations
can be solved as a pair of independent equations by the more efficient
cnexp method.

GHK based ionic currents for AMPA current added by Rishikesh Narayanan.

ENDCOMMENT

NEURON {
	POINT_PROCESS ghkampa
	USEION na WRITE ina
	USEION k WRITE ik
	
	RANGE taur, taud
	RANGE iampa
	RANGE P, Pmax
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	taur=2 (ms) <1e-9,1e9>
	taud = 10 (ms) <1e-9,1e9>
	nai = 18	(mM)	: Set for a reversal pot of +55mV
	nao = 140	(mM)
	ki = 140	(mM)	: Set for a reversal pot of -90mV
	ko = 5		(mM)
	celsius		(degC)
	Pmax=1e-6   (cm/s)	
}

ASSIGNED {
	ina     (nA)
	ik      (nA)
	v (mV)
	P (cm/s)
	factor
	iampa	(nA)

	Area (cm2)
}

STATE {
	A (cm/s)
	B (cm/s)
}

INITIAL {
	LOCAL tp
	if (taur/taud > .9999) {
		taur = .9999*taud
	}
	A = 0
	B = 0
	tp = (taur*taud)/(taud - taur) * log(taud/taur)
	factor = -exp(-tp/taur) + exp(-tp/taud)
	factor = 1/factor
	Area=1
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	P=B-A

: Area is just for unit conversion of ghk output

	ina = P*ghk(v, nai, nao,1)*Area	
	ik = P*ghk(v, ki, ko,1)*Area
	iampa = ik + ina
}

DERIVATIVE state {
	A' = -A/taur
	B' = -B/taud
}

FUNCTION ghk(v(mV), ci(mM), co(mM),z) (0.001 coul/cm3) {
	LOCAL arg, eci, eco
	arg = (0.001)*z*FARADAY*v/(R*(celsius+273.15))
	eco = co*efun(arg)
	eci = ci*efun(-arg)
	ghk = (0.001)*z*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

NET_RECEIVE(weight (uS)) { 	: No use to weight, can be used instead of Pmax,
							: if you want NetCon access to the synaptic
							: conductance.
	state_discontinuity(A, A + Pmax*factor)
	state_discontinuity(B, B + Pmax*factor)
}
