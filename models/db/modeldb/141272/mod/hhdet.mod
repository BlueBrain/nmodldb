TITLE hh.mod   squid sodium, potassium, and leak channels
 
COMMENT
 This is the original Hodgkin-Huxley treatment for the set of sodium, 
  potassium, and leakage channels found in the squid giant axon membrane.
  ("A quantitative description of membrane current and its application 
  conduction and excitation in nerve" J.Physiol. (Lond.) 117:500-544 (1952).)
 Membrane voltage is in absolute mV and has been reversed in polarity
  from the original HH convention and shifted to reflect a resting potential
  of -65 mV.
 Remember to set celsius=6.3 (or whatever) in your HOC file.
 See squid.hoc for an example of a simulation using this model.
 SW Jaslove  6 March, 1992
ENDCOMMENT
 
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}
 
NEURON {
	SUFFIX hhdet
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
}
 
STATE {
	m h n
}
 
ASSIGNED {
	v (mV)
	celsius (degC)
	ena (mV)
	ek (mV)

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
	m' = (1-m)*alpha_m - m*beta_m
	h' = (1-h)*alpha_h - h*beta_h
	n' = (1-n)*alpha_n - n*beta_n
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
}
 
