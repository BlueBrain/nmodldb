
COMMENT

kd21.mod (k slow) in the paper

Voltage gated k+ channels in layer5 neocortical pyramidal neurons from young rats: subtypes and gradients

Alon Korngreen and Bert Sakmann


ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kdf
	USEION k READ ek WRITE ik
	RANGE n,h, gk, gkbar
	RANGE ninf, ntau,hinf,htau,vshift
	RANGE ikd
	GLOBAL q10, temp, vmin, vmax,tadj,N
	GLOBAL ntau_slow1,ntau_slow2

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gkbar = 2.5   	(pS/um2)	: 0.03 mho/cm2
	v 		(mV)
								
	dt		(ms)
	celsius		(degC)
	temp = 23	(degC)		: original temp 	
	q10  = 2.3			: temperature sensitivity

	vmin = -150	(mV)
	vmax = 100	(mV)



	ntau_C1=1.25	(ms)
	ntau_C2=171 	:1.15	(ms)
	ntau_C3=-0.026	(1/mV)
	ntau_C4=13	(ms)
	
	ninf_a=-24		:-14	(mV)
	ninf_b=14	(mV)

	htau_C1=360	(ms)
	htau_C2=810	:1010(ms)
	htau_C3=-75	(mV)
	htau_C4=48	(mV)
	htau_C5=22	:24(ms/mV)
	htau_C6=-55	(mV)

	hinf_a=-54	(mV)
	hinf_b=-11	(mV)

	vshift=0	(mV)
	ntau_slow1=1
	ntau_slow2=0
	N=3
	ek		(mV)	
} 


ASSIGNED {

	ik 		(mA/cm2)
	gk		(pS/um2)
	ninf
	ntau (ms)	
	htau (ms)
	hinf
	tadj
	ikd		(mA/cm2)
	
}
STATE { n h}

INITIAL { 
	trates(v-vshift)
	n = ninf
	h = hinf
	tadj=q10^((celsius - temp)/10)
		
}
BREAKPOINT {
      SOLVE states METHOD cnexp
	gk = gkbar*h*( tadj)*n^N
	ik = (1e-4) * gk * (v - ek)
	ikd=ik
} 

DERIVATIVE states {   :Computes state variable n 

   	trates(v-vshift)      :             at the current v and dt.
	ntau=ntau*ntau_slow1+ntau_slow2 
	
      n'= (1- exp((-dt*tadj)/(ntau)))*(ninf-n)/dt
	h' = (1- exp((-dt*tadj)/(htau)))*(hinf-h)/dt
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
      TABLE ninf, ntau,hinf,htau
	DEPEND dt, celsius, temp, ntau_C1, ntau_C2, ntau_C3,ntau_C4,ninf_a,ninf_b,htau_C1, htau_C2 ,htau_C3,htau_C4,htau_C5,htau_C6   , hinf_a , hinf_b 

	FROM vmin TO vmax WITH 199

	rates(v): not consistently executed from here if usetable_hh == 1
}

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL ntau_C2_50,ntau_C3_50 
	ntau_C2_50=ntau_C2
	ntau_C3_50=ntau_C3
	if (v>-50){
		ntau_C2_50=ntau_C4
		ntau_C3_50=-ntau_C3
	}
	ntau=ntau_C1+ntau_C2_50*exp(-v*ntau_C3_50)
	ntau=ntau :/2
	htau=htau_C1+(htau_C2+htau_C5*(v-htau_C6))*exp(-((v-htau_C3)/htau_C4)^2)

	ninf=1/(1+exp((ninf_a-v)/ninf_b))

	hinf=1/(1+exp((hinf_a-v)/hinf_b))
}

