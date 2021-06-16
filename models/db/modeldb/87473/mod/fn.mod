COMMENT
	Sodium transient current for Av-Ron and Vidal, 1999
	Implemented by C Weaver, 2003

	Includes Fitzhugh-Nagumo / Rinzel type simplification
	of m to minf, and to replace 'h' with '1-n'


	equations:

	I_Na = gbar_Na * minf(V) ^3 * (1 - n) * (V - E_Na)
	I_K = gbar_K * n ^ 4 * (V - E_K)

	dn/dt = ( ninf(V) - n ) / ntau(V)

	minf(V) = 1 / ( 1 + exp(-2a*(V-Vha)))
	ninf(V) = 1 / ( 1 + exp(-2b*(V-Vhb)))
	ntau(V) = 1 / ( c * exp( b*(V-Vhb)) + c * exp (-b*(V-Vhb)) )

ENDCOMMENT

NEURON {
	SUFFIX fn
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	RANGE gnabar, gkbar, gkmodbar
	RANGE fastNashift, vha, vhb
	RANGE am, an, lamb
	RANGE minf, ninf, ntau
	RANGE totna, totk
}

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
}


PARAMETER { 
	fastNashift = 0	(mV)	: -3.5 (mV)
	gnabar = 0.0	(mho/cm2)
	gkbar = 0.0	(mho/cm2)
	am = 0.055	(/mV)
	vhm = -33	(mV)
	an = 0.055	(/mV)
	vhn = -40	(mV)
	lamb = 0.1		(1)
	tea = 0
	washfac = -0.005 (/ms)
}

ASSIGNED { 
	v	(mV)
	ena	(mV)  
	ek	(mV)  
	ina 	(mA/cm2) 
	ik 	(mA/cm2) 
	totna 	(mA/cm2) 
	totk 		(mA/cm2) 
	minf	(1)
	ninf	(1)
	ntau	(ms)

}

STATE { 
        n 
	gkmodbar (mho/cm2)
}


INITIAL {
	rates(v)
	if (tea > 0) {
		gkmodbar = gkbar
	}
	n = ninf
: printf( "fn ik=%g, ina =%g\n", ik, ina)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	totna = gnabar * minf^3 * (1-n) * ( v - ena ) 
	totk = gkbar * n^4 * ( v - ek )
	ina = gnabar * minf^3 * (1-n) * ( v - ena ) 
	ik = gkbar * n^4 * ( v - ek )
	if (tea > 0) {
		: to simulate TEA wash, a K-DR channel blocker
		ik = gkmodbar * n^4 * ( v - ek )
	}
}

DERIVATIVE states {
	rates(v)
	n' = (ninf-n)/ntau
	
	gkmodbar' = washfac * gkmodbar
}


UNITSOFF 

PROCEDURE rates(V (mV)) {
	minf = 1 / ( 1 + exp( -2 * am * ( V - fastNashift - vhm )) )
	ninf = 1 / ( 1 + exp( -2 * an * ( V - vhn )) )
	ntau = 1 / lamb / (exp( an*(V-vhn)) + exp (-an*(V-vhn)) )
}

UNITSON
