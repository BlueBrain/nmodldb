TITLE ERG CURRENT

COMMENT

ERG CURRENT
Coded up by Michael David Forrest using equations and parameters from:
Canavier CC, Oprisan SA, Callaway  JC, Ji H, Shepard PD (2007) Computational model predicts a role for ERG current in repolarizing plateau potentials in dopamine neurons: implications for modulation of neuronal activity. J Neurophysiol 98(5):3006-3022.

ENDCOMMENT


NEURON {
	SUFFIX erg
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT i
	RANGE gbar, g, ik, i, igate, nc, ca, cva, cka, cb, cvb, ckb, vth, delay, vhalf, vslope, vhalfhat, vslopehat
	GLOBAL ninf, tau
	GLOBAL gateCurrent, gunit
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(pA) = (picoamp)
	(S)  = (siemens)
	(mS) = (millisiemens)
	(nS) = (nanosiemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(molar) = (1/liter)
	(mM) = (millimolar)		
}

CONSTANT {
	e0 = 1.60217646e-19 (coulombs)
	q10 = 2.7
	zn = 1.9196 (1)		: valence of n-gate
}

PARAMETER {

	ca = 0.22 (1/ms)
	cva = 16 (mV)
	cka = -26.5 (mV)
	cb = 0.22 (1/ms)           : I MOVED ALL THESE OUT OF CONSTANT BLOCK, INTO PARAMETER BLOCK, SO THEY CAN BE MODULATED
	cvb = 16 (mV)
	ckb = 26.5 (mV)        


	gateCurrent = 0 (1)	: gating currents ON = 1 OFF = 0
	
	gbar = 0.005 (S/cm2)   <0,1e9>
	gunit = 16 (pS)		: unitary conductance 
        vth = -10 (mV)
        delay = 3 (ms)
      vhalf= -35 (mV)
      vslope = 5 (mV)
      vhalfhat= -70 (mV)
      vslopehat = -20 (mV)
}

ASSIGNED {
	celsius (degC)
	v (mV)
	
	ik (mA/cm2)
        neo (mA/cm2)
	igate (mA/cm2)
	i (mA/cm2)
 
	ek (mV)
	g (S/cm2)
	nc (1/cm2)
	qt (1)

	ninf (1)
	tau (ms)
	alpha (1/ms)
	beta (1/ms)
        gategate
        gatefkt
       mousegate 

        hinfhat (1)
	tauhat (ms)
	alphahat (1/ms)
	betahat (1/ms) 

}

STATE { n h }

INITIAL {
	nc = (1e12) * gbar / gunit
:	qt = q10^((celsius-22 (degC))/10 (degC))
        qt = 1
	rateConst(v)
	n = ninf
        h = hinfhat
        neo = ik
}

BREAKPOINT {
	SOLVE state METHOD cnexp
      g = gbar * n * h
	ik = g * (v - ek) 


:	igate = nc * (1e6) * e0 * 4 * zn * ngateFlip()
:	if (gateCurrent != 0) { 
:		i = igate
:	}


}

DERIVATIVE state {
	rateConst(v)
:	n' = alpha * (1-n) - beta * n
        n' = (ninf-n)/tau 
        h' = (hinfhat-h)/tauhat 
}

PROCEDURE rateConst(v (mV)) {
	alpha = qt * alphaFkt(v)
	beta = qt * betaFkt(v)
:	ninf = alpha / (alpha + beta) 
        ninf = ninfFkt(v)
	tau = 1 / (alpha + beta)

: /////////////

	alphahat = qt * alphaFkthat(v)
	betahat = qt * betaFkthat(v)
:	hinfhat = alphahat / (alphahat + betahat) 
        hinfhat = hinfFkthat(v)
	tauhat = 1 / (alphahat + betahat)

}

FUNCTION alphaFkt(v (mV)) (1/ms) {
	alphaFkt = 0.00225 * exp(0.12 * v)
}

FUNCTION betaFkt(v (mV)) (1/ms) {
	betaFkt = 0.00004 * exp(-0.05 * v)
}

FUNCTION ninfFkt(v (mV)) (1/ms) {
	ninfFkt = 1 / ( 1 + exp(-(v - vhalf)/vslope)   )
}

: ///////////////////////////


FUNCTION alphaFkthat(v (mV)) (1/ms) {
	alphaFkthat =  0.1 * exp(0.02 * v)
}

FUNCTION betaFkthat(v (mV)) (1/ms) {
	betaFkthat = 0.003 * exp(-0.03 * v)
}

FUNCTION hinfFkthat(v (mV)) (1/ms) {
	hinfFkthat =  1 / ( 1 + exp(-(v - vhalfhat)/vslopehat)   )
}




:FUNCTION ngateFlip() (1/ms) {
:	ngateFlip = (ninf-n)/tau 
:}




