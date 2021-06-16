COMMENT
Original Hodgkin and Huxley model (J.Physiol. (Lond.) 117:500-544 (1952))
with stochastic conductances, using coupled activation particles (5-state K 
channels, 8-state Na channels) and Markov Chain modeling (Chow & White algorithm)


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
	SUFFIX hh58CW
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gkbar, gl, el, NNa, NK, next_evK, next_evNa, Nast, Kst
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

	ina (mA/cm2)
	ik (mA/cm2)
	il (mA/cm2)
	alpha_m	(/ms)
	alpha_h	(/ms)
	alpha_n	(/ms)
	beta_m	(/ms)
	beta_h	(/ms)
	beta_n	(/ms)
	Nast[8]
	Kst[5]
	Nart[20]	(/ms)
	Krt[8]		(/ms)
	sumrtNa		(/ms)
	sumrtK		(/ms)
	cumsumNa[20](/ms)
	cumsumK[8]	(/ms)
	next_evNa	(ms)
	next_evK	(ms)
	prev_ev		(ms)
	ev			(/ms)
	M
	N
	H
	
}
 
STATE {mock}

BREAKPOINT {
	SOLVE mula METHOD cnexp
	ina = gnabar*Nast[7]*(v - ena)/NNa
	ik = gkbar*Kst[4]*(v - ek)/NK      
	il = gl*(v - el)
}
 
INITIAL {
	LOCAL stsum, q10
	q10 = 3^((celsius - 6.3)/10)
	alpha_m = q10*0.1*(v+40)/(1-exp(-(v+40)/10))
	beta_m = q10*4*exp(-(v+65)/18)
	alpha_h = q10*0.07*exp(-(v+65)/20) 
	beta_h = q10/(1+exp(-(v+35)/10))
	alpha_n = q10*0.01*(v+55)/(1-exp(-(v+55)/10))
	beta_n = q10*0.125*exp(-(v+65)/80)
	
	M=alpha_m/beta_m
	H=alpha_h/beta_h
	N=alpha_n/beta_n
	stsum=(1+H)*(1+M)^3
	Nast[0]=floor(NNa/stsum+0.5)
	Nast[1]=floor(NNa*3*M/stsum+0.5)
	Nast[2]=floor(NNa*3*M^2/stsum+0.5)
	Nast[3]=floor(NNa*M^3/stsum+0.5)
	Nast[4]=floor(NNa*H/stsum+0.5)
	Nast[5]=floor(NNa*H*3*M/stsum+0.5)
	Nast[6]=floor(NNa*H*3*M^2/stsum+0.5)
	Nast[7]=floor(NNa*H*M^3/stsum+0.5)
	ratesNa(v)
	next_evNa = - log(scop_random())/sumrtNa
	
	stsum=(1+N)^4
	Kst[0]=floor(NK/stsum+0.5)
	Kst[1]=floor(NK*4*N/stsum+0.5)
	Kst[2]=floor(NK*6*N^2/stsum+0.5)
	Kst[3]=floor(NK*4*N^3/stsum+0.5)
	Kst[4]=floor(NK*N^4/stsum+0.5)
	ratesK(v)
	next_evK = - log(scop_random())/sumrtK
}

DERIVATIVE mula {  
	while (t>= next_evNa){
		transNa()
	}
	while (t>= next_evK){
		transK()
	}
	mock'=0
}
 
LOCAL q10

PROCEDURE ratesNa(v(mV)) {  :Computes rate and other constants at current v.
	UNITSOFF
	q10 = 3^((celsius - 6.3)/10)
	alpha_m = q10*0.1*(v+40)/(1-exp(-(v+40)/10))
	beta_m = q10*4*exp(-(v+65)/18)
	alpha_h = q10*0.07*exp(-(v+65)/20) 
	beta_h = q10/(1+exp(-(v+35)/10))
	Nart[0]=3*alpha_m*Nast[0]
	Nart[1]=beta_m*Nast[1]
	Nart[2]=2*alpha_m*Nast[1]
	Nart[3]=2*beta_m*Nast[2]
	Nart[4]=alpha_m*Nast[2]
	Nart[5]=3*beta_m*Nast[3]
	Nart[6]=alpha_h*Nast[0]
	Nart[7]=beta_h*Nast[4]
	Nart[8]=alpha_h*Nast[1]
	Nart[9]=beta_h*Nast[5]
	Nart[10]=alpha_h*Nast[2]
	Nart[11]=beta_h*Nast[6]
	Nart[12]=alpha_h*Nast[3]
	Nart[13]=beta_h*Nast[7]
	Nart[14]=3*alpha_m*Nast[4]
	Nart[15]=beta_m*Nast[5]
	Nart[16]=2*alpha_m*Nast[5]
	Nart[17]=2*beta_m*Nast[6]
	Nart[18]=alpha_m*Nast[6]
	Nart[19]=3*beta_m*Nast[7]
	sumrtNa=0
	FROM ii=0 TO 19 {
		sumrtNa = sumrtNa + Nart[ii]
		cumsumNa[ii] = sumrtNa
	}
	FROM ii=0 TO 19 {cumsumNa[ii] = cumsumNa[ii] / sumrtNa}
	UNITSON
}

PROCEDURE ratesK(v(mV)) {
	UNITSOFF
	q10 = 3^((celsius - 6.3)/10)
	alpha_n = q10*0.01*(v+55)/(1-exp(-(v+55)/10))
	beta_n = q10*0.125*exp(-(v+65)/80)

	Krt[0]=4*alpha_n*Kst[0]
	Krt[1]=beta_n*Kst[1]
	Krt[2]=3*alpha_n*Kst[1]
	Krt[3]=2*beta_n*Kst[2]
	Krt[4]=2*alpha_n*Kst[2]
	Krt[5]=3*beta_n*Kst[3]
	Krt[6]=alpha_n*Kst[3]
	Krt[7]=4*beta_n*Kst[4]
	sumrtK=0
	FROM ii=0 TO 7 {
		sumrtK = sumrtK + Krt[ii]
		cumsumK[ii] = sumrtK
	}
	FROM ii=0 TO 7 {cumsumK[ii] = cumsumK[ii] / sumrtK}

	UNITSON
}

PROCEDURE transK() {
	ratesK(v)
	ev = scop_random()*1(/ms)
	if (ev <= cumsumK[0]) {
		Kst[0]=Kst[0]-1
		Kst[1]=Kst[1]+1
	}
	if (cumsumK[0] < ev && ev <= cumsumK[1]) {
		Kst[0]=Kst[0]+1
		Kst[1]=Kst[1]-1
	}	
	if (cumsumK[1] < ev && ev <= cumsumK[2]) {
		Kst[1]=Kst[1]-1
		Kst[2]=Kst[2]+1
	}
	if (cumsumK[2] < ev && ev <= cumsumK[3]) {
		Kst[1]=Kst[1]+1
		Kst[2]=Kst[2]-1
	}	
	if (cumsumK[3] < ev && ev <= cumsumK[4]) {
		Kst[2]=Kst[2]-1
		Kst[3]=Kst[3]+1
	}
	if (cumsumK[4] < ev && ev <= cumsumK[5]) {
		Kst[2]=Kst[2]+1
		Kst[3]=Kst[3]-1
	}
	if (cumsumK[5] < ev && ev <= cumsumK[6]) {
		Kst[3]=Kst[3]-1
		Kst[4]=Kst[4]+1
	}
	if (cumsumK[6] < ev && ev <= cumsumK[7]) {
		Kst[3]=Kst[3]+1
		Kst[4]=Kst[4]-1
	}
	prev_ev = next_evK
	next_evK = prev_ev - log(scop_random())/sumrtK
}

 
PROCEDURE transNa() {
	ratesNa(v)
	ev = scop_random()*1(/ms)
	if (ev <= cumsumNa[0]) {
		Nast[0]=Nast[0]-1
		Nast[1]=Nast[1]+1
	}
	if (cumsumNa[0] < ev && ev <= cumsumNa[1]) {
		Nast[0]=Nast[0]+1
		Nast[1]=Nast[1]-1
	}	
	if (cumsumNa[1] < ev && ev <= cumsumNa[2]) {
		Nast[1]=Nast[1]-1
		Nast[2]=Nast[2]+1
	}
	if (cumsumNa[2] < ev && ev <= cumsumNa[3]) {
		Nast[1]=Nast[1]+1
		Nast[2]=Nast[2]-1
	}	
	if (cumsumNa[3] < ev && ev <= cumsumNa[4]) {
		Nast[2]=Nast[2]-1
		Nast[3]=Nast[3]+1
	}
	if (cumsumNa[4] < ev && ev <= cumsumNa[5]) {
		Nast[2]=Nast[2]+1
		Nast[3]=Nast[3]-1
	}
	if (cumsumNa[5] < ev && ev <= cumsumNa[6]) {
		Nast[0]=Nast[0]-1
		Nast[4]=Nast[4]+1
	}
	if (cumsumNa[6] < ev && ev <= cumsumNa[7]) {
		Nast[0]=Nast[0]+1
		Nast[4]=Nast[4]-1
	}
	if (cumsumNa[7] < ev && ev <= cumsumNa[8]) {
		Nast[1]=Nast[1]-1
		Nast[5]=Nast[5]+1
	}
	if (cumsumNa[8] < ev && ev <= cumsumNa[9]) {
		Nast[1]=Nast[1]+1
		Nast[5]=Nast[5]-1
	}
	if (cumsumNa[9] < ev && ev <= cumsumNa[10]) {
		Nast[2]=Nast[2]-1
		Nast[6]=Nast[6]+1
	}
	if (cumsumNa[10] < ev && ev <= cumsumNa[11]) {
		Nast[2]=Nast[2]+1
		Nast[6]=Nast[6]-1
	}				
	if (cumsumNa[11] < ev && ev <= cumsumNa[12]) {
		Nast[3]=Nast[3]-1
		Nast[7]=Nast[7]+1
	}
	if (cumsumNa[12] < ev && ev <= cumsumNa[13]) {
		Nast[3]=Nast[3]+1
		Nast[7]=Nast[7]-1
	}
	if (cumsumNa[13] < ev && ev <= cumsumNa[14]) {
		Nast[4]=Nast[4]-1
		Nast[5]=Nast[5]+1
	}
	if (cumsumNa[14] < ev && ev <= cumsumNa[15]) {
		Nast[4]=Nast[4]+1
		Nast[5]=Nast[5]-1
	}
	if (cumsumNa[15] < ev && ev <= cumsumNa[16]) {
		Nast[5]=Nast[5]-1
		Nast[6]=Nast[6]+1
	}
	if (cumsumNa[16] < ev && ev <= cumsumNa[17]) {
		Nast[5]=Nast[5]+1
		Nast[6]=Nast[6]-1
	}
	if (cumsumNa[17] < ev && ev <= cumsumNa[18]) {
		Nast[6]=Nast[6]-1
		Nast[7]=Nast[7]+1
	}
	if (cumsumNa[18] < ev && ev <= cumsumNa[19]) {
		Nast[6]=Nast[6]+1
		Nast[7]=Nast[7]-1
	}
	prev_ev = next_evNa
	next_evNa = prev_ev - log(scop_random())/sumrtNa
}
