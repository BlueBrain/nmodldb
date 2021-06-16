TITLE CaL channel
: Calcium Channel with Goldman- Hodgkin-Katz permeability


NEURON {
	SUFFIX 	caL
	USEION 	ca READ cai, cao WRITE ica
	RANGE 	pbar, p, i
	GLOBAL	inf_n, tau_n, inf_h, tau_h, inf_c, tau_c
}

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)


	:FARADAY = 96520 (coul)
	:R = 8.3134 (joule/degC)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
   : inf_n parameters
	pbar	= 0.241e-3	(cm/s)	: Maximum Permeability
	gates_n       = 2 
  	vhalf_n       = -30.1 (mV)
  	slope_n       = -4.77  (mV)
  	needAdj = 1
  	vhalfA_n      = 0 (mV)                : adjusted for ngate power; Set in ngate_adjust()
  	slopeA_n      = 0 (mV)
  	v5_adj        = 0 (mV)                : for return values in ngate_adjust
  	slp_adj       = 0 (mV)

  : tau_n parameters
  	tauA_n        = 1.5  (ms)
  	tauDv_n       = 0     (mV)    : Delta to vhalf_n
  	tauG_n        = 0.5          : Left-right bias. range is (0,1)
  	tauF_n        = 0             : Up-Down bias. range is ~ -3.5(cup-shape), -3(flat), 0(from k), 1(sharper)
  	tau0_n        = 0   (ms)    : minimum tau


  : inf_h parameters
	vhalf_h	= 	-59.5	(mV)
	slope_h	=	10	(mV)
:	vhat_h	=	-59.5	(mV)
:	shat_h	=	10	(mV)
	tauA_h	=	200	(ms)
	tauG_h	=	0.5		: Left-right bias.  range (0,1)
	tau0_h	=	0	(ms)
	tauF_h	=	0
	tauDv_h	=	0

	hill_c	=	2
	K_c		=	.001	(mM)
	tauA_c	=	10	(ms)
	tau0_c	=	10	(ms)

}

ASSIGNED { 
	celsius		(degC) : 32
	v		(mV)
	i		(mA/cm2)
	ica		(mA/cm2)
	cai		(mM)
	cao		(mM)
	p		(cm/s)
	inf_n
	tau_n		(ms)
	inf_h
	tau_h		(ms)
	inf_c
	tau_c		(ms)
}

STATE {
	n	: activation
	h	: inactivation
	c	: calcium dependent inactivation
}		

BREAKPOINT {
	SOLVE states METHOD cnexp
	p 	= 0
	p 	= pbar * n^2 * h * c
	i 	= p * ghk(v, cai, cao)
	ica 	= i
}

INITIAL {
	rates(v, cai)
	n = inf_n
	h = inf_h
	c = inf_c
}

DERIVATIVE states {
	rates(v, cai)
	n' = ( inf_n - n) / tau_n
	h' = ( inf_h - h) / tau_h
	c' = ( inf_c - c) / tau_c
}

PROCEDURE rates (v (mV), cai ( mM)) {
  if( needAdj > 0 ){
    needAdj = 0
    ngate_adjust( gates_n, vhalf_n, slope_n )
    vhalfA_n = v5_adj
    slopeA_n = slp_adj
  }
  	inf_n = Boltzmann( v, vhalfA_n, slopeA_n )
  	tau_n = BorgMod_tau( v, vhalfA_n, slopeA_n, tau0_n, tauA_n, tauG_n, tauF_n, tauDv_n )
	inf_h = 1 / (1 + exp((v - vhalf_h)/slope_h))
	tau_h = tau0_h + tauA_h*4*sqrt(tauG_h * (1 - tauG_h))/(exp(-(v - vhalf_h)/slope_h*(1 - tauG_h)) + exp((v - vhalf_h)/slope_h*tauG_h))
	inf_c = 1/(1 + (cai/K_c)^hill_c)
	tau_c = tau0_c + tauA_c/(1 + (cai/K_c)^hill_c)
	
}
FUNCTION Boltzmann( v (mV), v5 (mV), s (mV) ){
  Boltzmann = 1 / (1 + exp( (v - v5) / s ))
}

FUNCTION BorgMod_tau( v (mV), v5 (mV), s (mV), tau0 (ms), tauA (ms), tauG, tauF, tauDv (mV) ) (ms) {
  LOCAL kc, kr, Dv, wr, kf

:  kr = 1000
:  wr = 1000
  Dv = (v - ( v5 + tauDv ))
:  kc =  kr * 10^tauF / s *1(mV)
  kf =  10^tauF

  BorgMod_tau = tau0 + tauA * 4 * sqrt( tauG * (1-tauG))
                / ( exp( - Dv *(1-tauG)*kf/s ) + exp( Dv *tauG*kf/s ))
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
	LOCAL z, eci, eco
	z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
	eco = co*efun(z)
	eci = ci*efun(-z)
	:high cao charge moves inward
	:negative potential charge moves inward
	ghk = (.001)*2*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

: Boltzmann's inverse
FUNCTION Boltz_m1( x, v5 (mV), s (mV) ) (mV) {
  Boltz_m1 = s * log( 1/x - 1 ) + v5
}

: Find parameters for a Boltzmann eq that when taken to the ngate power matches one with a single power
: return result in v5_adj and slp_adj
: We solve for exact match on two points
PROCEDURE ngate_adjust( ng, vh (mV), slp (mV) ) {
  LOCAL x1, x2, v1, v2
  x1 = 0.3
  x2 = 0.7
  v1 = Boltz_m1( x1, vh, slp )
  v2 = Boltz_m1( x2, vh, slp )
  slp_adj = (v2 - v1)/( log( (1/x2)^(1/ng) - 1 ) - log( (1/x1)^(1/ng) - 1 ) )
  v5_adj = v1 - slp_adj * log( 1 / x1^(1/ng) - 1 )
}
