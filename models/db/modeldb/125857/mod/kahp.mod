COMMENT

	Maciej T. Lazarewicz, mlazarew@seas.upenn.edu

ENDCOMMENT

NEURON {

	SUFFIX kahp
	USEION k WRITE ik
	USEION ca READ cai
	RANGE gkbar, ik, qinf, tauq
}
	
UNITS {

    (mollar) = (1/liter)
	(mM)     = (millimollar)
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {

    gkbar = 0.8 (mS/cm2)
    ek   = -75 (mV)
}
    
ASSIGNED {

    v    (mV)
    ik   (mA/cm2)
    cai  (mM)
    qinf (1)
    tauq (ms)
}

STATE { q }

INITIAL { 
    
    rates(v)
    q  = qinf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	
	ik = (1e-3) * gkbar * q * (v-ek)
}


DERIVATIVE states { 

    rates(v)
    q' = (qinf-q)/tauq
        
}

PROCEDURE rates(v(mV)) { LOCAL a,b
    
    a = 0.01(/ms) * min(cai/500(mM),1)
    b = 1(/ms)/1000
    
    qinf = a/(a+b)
    tauq = 1.0/(a+b)
}

COMMENT

	Maciej T. Lazarewicz, mlazarew@seas.upenn.edu

ENDCOMMENT



:-------------------------------------------------------------------
FUNCTION fun1(v(mV),V0(mV),A(/ms),B(mV))(/ms) {

	 fun1 = A*exp((v-V0)/B)
}

FUNCTION fun2(v(mV),V0(mV),A(/ms),B(mV))(/ms) {

	 fun2 = A/(exp((v-V0)/B)+1)
}

FUNCTION fun3(v(mV),V0(mV),A(/ms),B(mV))(/ms) {

    if(fabs((v-V0)/B)<1e-6) {
    :if(v==V0) {
        fun3 = A*B/1(mV) * (1- 0.5 * (v-V0)/B)
    } else {
        fun3 = A/1(mV)*(v-V0)/(exp((v-V0)/B)-1)
    }
}

FUNCTION min(x,y) { if (x<=y){ min = x }else{ min = y } }