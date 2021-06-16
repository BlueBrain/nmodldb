COMMENT
Original Hodgkin and Huxley model (J.Physiol. (Lond.) 117:500-544 (1952))
with stochastic conductances, using uncoupled (2-state) activation particles and
Markov Chain modeling.

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
	SUFFIX hh2CW
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gkbar, gl, el, NNa, NK, Nn, Nm, Nh, next_evm, next_evh, next_evn
}
 
PARAMETER {
	gnabar = .12 (S/cm2)	<0,1e9>
	gkbar = .036 (S/cm2)	<0,1e9>
	gl = .0003 (S/cm2)	<0,1e9>
	el = -54.3 (mV)
	NNa = 500
	NK = 160
}
 
ASSIGNED {
	v (mV)
	celsius (degC)
	ena (mV)
	ek (mV)

	ina (mA/cm2)
	ik (mA/cm2)
	il (mA/cm2)
	alpha_m	(/ms)
	alpha_h	(/ms)
	alpha_n	(/ms)
	beta_m	(/ms)
	beta_h	(/ms)
	beta_n	(/ms)
	m
	h
	n
	Nn
	Nm
	Nh
	next_evm	(ms)
	next_evh	(ms)
	next_evn	(ms)
	prev_ev		(ms)
}
 
BREAKPOINT {
	states()
	m = Nm / (NNa*3)
	h = Nh / NNa
	n = Nn / (NK*4)
	ina = gnabar*m*m*m*h*(v - ena)
	ik = gkbar*n*n*n*n*(v - ek)      
	il = gl*(v - el)
}
 
INITIAL {
	rates(v)
	m=alpha_m/(alpha_m + beta_m)
	h=alpha_h/(alpha_h + beta_h)
	n=alpha_n/(alpha_n + beta_n)
	Nm=floor(NNa*3*m+0.5)
	Nh=floor(NNa*h+0.5)
	Nn=floor(NK*4*n+0.5)
	next_evm = - log(scop_random())/(Nm*beta_m+(NNa*3-Nm)*alpha_m)
	next_evh = - log(scop_random())/(Nh*beta_h+(NNa-Nh)*alpha_h)
	next_evn = - log(scop_random())/(Nn*beta_n+(NK*4-Nn)*alpha_n)
}

PROCEDURE states() {  
	rates(v)
	while (t>= next_evm){
		transm()
	}
	while (t>= next_evh){
		transh()
	}
	while (t>= next_evn){
		transn()
	}
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
	UNITSON
}
 
PROCEDURE transm() {
	if (scop_random() >= Nm*beta_m/((NNa*3-Nm)*alpha_m + Nm*beta_m)) {Nm = Nm + 1} else {Nm = Nm - 1}
	prev_ev = next_evm
	next_evm = prev_ev - log(scop_random())/(Nm*beta_m+(NNa*3-Nm)*alpha_m)
}

PROCEDURE transh() {
	if (scop_random() >= Nh*beta_h/((NNa-Nh)*alpha_h + Nh*beta_h)) {Nh = Nh + 1} else {Nh = Nh - 1}
	prev_ev = next_evh
	next_evh = prev_ev - log(scop_random())/(Nh*beta_h+(NNa-Nh)*alpha_h)
}

PROCEDURE transn() {
	if (scop_random() >= Nn*beta_n/((NK*4-Nn)*alpha_n + Nn*beta_n)) {Nn = Nn + 1} else {Nn = Nn - 1}
	prev_ev = next_evn
	next_evn = prev_ev - log(scop_random())/(Nn*beta_n+(NK*4-Nn)*alpha_n)
}