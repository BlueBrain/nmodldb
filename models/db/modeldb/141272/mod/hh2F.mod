COMMENT
Original Hodgkin and Huxley model (J.Physiol. (Lond.) 117:500-544 (1952))
with stochastic conductances, using uncoupled (2-state) activation particles and
Diffusion approximation (Fox) algorithm.

Membrane voltage is in absolute mV and has been reversed in polarity
from the original HH convention and shifted to reflect a resting potential
of -65 mV.
ENDCOMMENT
 
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}
 
NEURON {
	SUFFIX hh2F
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gkbar, gl, el, gna, gk
}
 
PARAMETER {
	gnabar = .12 (S/cm2)	<0,1e9>
	gkbar = .036 (S/cm2)	<0,1e9>
	gl = .0003 (S/cm2)	<0,1e9>
	el = -54.3 (mV)
	NNa = 500
	NK = 160
}
 
STATE {
	m h n
}
 
ASSIGNED {
	v (mV)
	celsius (degC)
	ena (mV)
	ek (mV)
	dt	(ms)
	gna (S/cm2)
	gk (S/cm2)
	ina (mA/cm2)
	ik (mA/cm2)
	il (mA/cm2)
	alpha_m	(/ms)
	alpha_h	(/ms)
	alpha_n	(/ms)
	beta_m	(/ms)
	beta_h	(/ms)
	beta_n	(/ms)
	SDn
	SDm
	SDh
}
 
BREAKPOINT {
	SOLVE states METHOD cnexp
	gna = gnabar*m*m*m*h
	ina = gna*(v - ena)
	gk = gkbar*n*n*n*n
	ik = gk*(v - ek)      
	il = gl*(v - el)
}
 
INITIAL {
	rates(v)
	m=alpha_m/(alpha_m + beta_m)
	h=alpha_h/(alpha_h + beta_h)
	n=alpha_n/(alpha_n + beta_n)
}

DERIVATIVE states {  
	rates(v)
	m' = (1-m)*alpha_m - m*beta_m + normrand(0,SDm)
	h' = (1-h)*alpha_h - h*beta_h + normrand(0,SDh)
	n' = (1-n)*alpha_n - n*beta_n + normrand(0,SDn)
}
 
LOCAL q10

PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.

	UNITSOFF
	q10 = 3^((celsius - 6.3)/10)
	alpha_n = q10*0.01*(v+55)/(1-exp(-(v+55)/10))
	beta_n = q10*0.125*exp(-(v+65)/80)
	alpha_m = q10*0.1*(v+40)/(1-exp(-(v+40)/10))
	beta_m = q10*4*exp(-(v+65)/18)
	alpha_h = q10*0.07*exp(-(v+65)/20) 
	beta_h = q10/(1+exp(-(v+35)/10))
	SDn = sqrt(fabs(alpha_n*(1-n)+beta_n*n)/(dt*NK*4))
	SDm = sqrt(fabs(alpha_m*(1-m)+beta_m*m)/(dt*NNa*3))
	SDh = sqrt(fabs(alpha_h*(1-h)+beta_h*h)/(dt*NNa))
	UNITSON
}
 
