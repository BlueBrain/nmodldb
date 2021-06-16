COMMENT
Original Hodgkin and Huxley model (J.Physiol. (Lond.) 117:500-544 (1952))
with stochastic conductances, using coupled activation particles (5-state K 
channels, 8-state Na channels) and Diffusion approximation (Fox) algorithm
with steady-state values of variables in the stochastic terms

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
	SUFFIX hh58F1ss
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gkbar, gl, el, NNa, NK, sumN, sumMH
}
 
PARAMETER {
	gnabar = .12 (S/cm2)	<0,1e9>
	gkbar = .036 (S/cm2)	<0,1e9>
	gl = .0003 (S/cm2)	<0,1e9>
	el = -54.3 (mV)
	NNa = 5000
	NK = 1600 
}
 
ASSIGNED {
	v (mV)
	celsius (degC)
	ena (mV)
	ek (mV)
	dt (ms)
	ina (mA/cm2)
	ik (mA/cm2)
	il (mA/cm2)
	am	(/ms)
	ah	(/ms)
	an	(/ms)
	bm	(/ms)
	bh	(/ms)
	bn	(/ms)
	stsum
	M
	N
	H
	R[14]	(/ms)
	sumN
	sumMH
	mh0
	n0
	
}
 
STATE {	
	mh1
	mh2
	mh3
	mh4
	mh5
	mh6
	mh7
	n1
	n2
	n3
	n4
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar*mh7*(v - ena)
	ik = gkbar*n4*(v - ek)
	il = gl*(v - el)
	sumN = n0 + n1 + n2 + n3 + n4
	sumMH = mh0 + mh1 + mh2 + mh3 + mh4 + mh5 + mh6 + mh7
}
 
INITIAL {
	rates(v)	
	
	M=am/bm
	H=ah/bh
	N=an/bn
	stsum=(1+H)*(1+M)^3
	mh0=1/stsum
	mh1=3*M/stsum
	mh2=3*M^2/stsum
	mh3=M^3/stsum
	mh4=H/stsum
	mh5=H*3*M/stsum
	mh6=H*3*M^2/stsum
	mh7=H*M^3/stsum
	
	stsum=(1+N)^4
	n0=1/stsum
	n1=4*N/stsum
	n2=6*N^2/stsum
	n3=4*N^3/stsum
	n4=N^4/stsum
	rates(v)
}

DERIVATIVE states {  
	rates(v)
	mh1' = (-2*am-bm-ah)*mh1 + 3*am*mh0 + 2*bm*mh2 + bh*mh5 -R[0]+R[1]+R[4]	
	mh2' = (-am-2*bm-ah)*mh2 + 2*am*mh1 + 3*bm*mh3 + bh*mh6 -R[1]+R[2]+R[5]
	mh3' = (-3*bm-ah)*mh3 + am*mh2 + bh*mh7 -R[2]+R[6]
	mh4' = (-3*am-bh)*mh4 + bm*mh5 + ah*mh0 + R[7]-R[3]
  mh5' = (-2*am-bm-bh)*mh5 + 3*am*mh4 + 2*bm*mh6 + ah*mh1 -R[7]+R[8]-R[4]
  mh6' = (-am-2*bm-bh)*mh6 + 2*am*mh5 + 3*bm*mh7 + ah*mh2 -R[8]+R[9]-R[5]
  mh7' = (-3*bm-bh)*mh7 + am*mh6 + ah*mh3 -R[9]-R[6]
  mh0 = 1-mh1-mh2-mh3-mh4-mh5-mh6-mh7
	
	n1' = (-3*an-bn)*n1 + 4*an*n0 + 2*bn*n2 - R[10] + R[11]
	n2' = (-2*an-2*bn)*n2 + 3*an*n1 + 3*bn*n3 -R[11] + R[12]
	n3' = (-an-3*bn)*n3 + 2*an*n2 + 4*bn*n4 -R[12] + R[13]
	n4' = -4*bn*n4 + an*n3 -R[13]
	n0 = 1-n1-n2-n3-n4		
}
 
LOCAL q10

PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
	LOCAL q10
	UNITSOFF
	q10 = 3^((celsius - 6.3)/10)
	am = q10*0.1*(v+40)/(1-exp(-(v+40)/10))
	bm = q10*4*exp(-(v+65)/18)
	ah = q10*0.07*exp(-(v+65)/20) 
	bh = q10/(1+exp(-(v+35)/10))
	an = q10*0.01*(v+55)/(1-exp(-(v+55)/10))
	bn = q10*0.125*exp(-(v+65)/80)
		
	FROM ii=0 TO 9 {R[ii]=normrand(0,1/sqrt(NNa*dt*(ah+bh)*(am+bm^3)))}
	FROM ii=10 TO 13 {R[ii]=normrand(0,1/sqrt(NK*dt*(an+bn)^4))}
	R[0] = R[0]*sqrt(6*am*bh*bm^3)
	R[1] = R[1]*sqrt(12*am^2*bh*bm^2)
	R[2] = R[2]*sqrt(6*am^3*bh*bm)
	R[3] = R[3]*sqrt(2*ah*bh*bm^3)
	R[4] = R[4]*sqrt(6*ah*bh*am*bm^2)
	R[5] = R[5]*sqrt(6*ah*bh*am^2*bm)
	R[6] = R[6]*sqrt(2*ah*bh*am^3)
	R[7] = R[7]*sqrt(6*am*ah*bm^3)
	R[8] = R[8]*sqrt(12*am^2*ah*bm^2)
	R[9] = R[9]*sqrt(6*am^3*ah*bm)
	R[10] = R[10]*sqrt(8*an*bn^4) 
	R[11] = R[11]*sqrt(24*an^2*bn^3)
	R[12] = R[12]*sqrt(24*an^3*bn^2)
	R[13] = R[13]*sqrt(8*an^4*bn)
	UNITSON 
}
