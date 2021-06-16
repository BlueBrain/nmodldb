
COMMENT

kv.mod

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kv
	USEION k READ ek WRITE ik
	RANGE n, h, gk, gbar
	RANGE ninf, ntau, hinf, htau
	GLOBAL q10, temp, vmin, vmax, v05, za, v05h,zh, Tscale
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gbar = 0   	(pS/um2)	
	v 		(mV)
								
	v05 = -15 (mV)
	za = 11 (mV)

	v05h = -18 (mV)
	zh = -11 (mV)
	
	dt		(ms)
	celsius		(degC)
	temp = 33	(degC)		: original temp 	
	q10  = 2.3				: temperature sensitivity

	Tscale = 10	(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
} 


ASSIGNED {
	ik 		(mA/cm2)
	ek		(mV)
	gk   		(pS/um2)
	ninf
	hinf
	ntau (ms)	
	htau (ms)
}
 

STATE { n h }

INITIAL { 
	rates(v)
	n = ninf
	h = hinf	
}

BREAKPOINT {
      SOLVE states METHOD cnexp
	gk=gbar*n*h
	ik = (1e-4)*gk*(v - ek)
} 


DERIVATIVE states {    
        rates(v)      
        n' = (ninf - n)/ntau
        h' =  (hinf - h)/htau
}

PROCEDURE rates(vm) {  	  
	  LOCAL qt

          qt=q10^((celsius-temp)/Tscale)
	  ninf = 1/(1 + exp(-(vm - v05)/za))
          ntau = (0.31+1492/((vm+23.5)^2+565)) :/qt

	  hinf = 1/(1 + exp(-(vm - v05h)/zh))
	  htau=300
}
