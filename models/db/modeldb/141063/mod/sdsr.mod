TITLE Isd and Isr currents of the Huber-Braun Model
: Braun et al. Int J Bifurcation and Chaos 8(5):881-889 (1998)
: Slow and subthreshol-activated Na+ and K+ currents responsible for oscillation
: non-Inactivating INap

: ******
: sr (slow potassium) conductance was modified adding a saturating function
: now asr resembles calcium concentration and ik is a Hill function of asr,
: therefore asr can be greater than 1 without making the open probability of
: KCa greater than 1
: ****

: Written by Patricio Orio, Jul 2006
:

NEURON {
	SUFFIX sdsr
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	RANGE gsd, gsr, isd
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gsd  = 0.0004   (mho/cm2)
	gsr  = 0.0004   (mho/cm2)
	V0sd = -40      (mV)
	zsd = 0.1		(/mV)
    eta = 12        (cm2/mA)
    k = 0.17        (1)
    tsd = 10        (ms)
    tsr = 24        (ms)
	n = 2
	Kd = 0.4	
}

STATE {
	asd
	asr
}

ASSIGNED {
	celsius	(degC)
	ina     (mA/cm2)
	ik      (mA/cm2)
    v       (mV)
    rho     (1)
    ena     (mV)
    ek      (mV)
    isd     (mA/cm2)
}

INITIAL {
    rho = 1.3^((celsius - 25 (degC))/10(degC))
    asd = 1/(1+exp(-zsd*(v - V0sd)))
    asr = (-eta * asd * rho * gsd * (v-ena))/k
    if (asr < 0) {asr = 0}
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    rho = 1.3^((celsius - 25 (degC))/10(degC))
	isd = rho * gsd * asd * (v - ena)
	ina = isd
	ik  = rho * gsr * (v - ek) * asr^n /(Kd^n + asr^n)
}

DERIVATIVE states {
    LOCAL phi, asdinf
    phi = 3^((celsius - 25 (degC))/ 10 (degC))
    asdinf = 1/(1+exp(-zsd*(v - V0sd)))
    asd' = phi * (asdinf - asd) / tsd
    asr' = phi * (-eta * isd - k*asr)/tsr
}

