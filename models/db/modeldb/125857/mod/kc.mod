COMMENT

	Maciej T. Lazarewicz, mlazarew@seas.upenn.edu

ENDCOMMENT



NEURON {

	SUFFIX kc
	USEION k WRITE ik
	USEION ca READ cai
	RANGE gkbar, ik
}
	
UNITS {    

    (mollar) = (1/liter)
	(mM)     = (millimollar)
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {

    gkbar = 15    (mS/cm2)
    ek = -75    (mV)
}
    
ASSIGNED { 

    ik   (mA/cm2)    
    v    (mV)
    cai  (mM)
    cinf (1)
    tauc (ms)
}

STATE { c }

INITIAL { 
    
    rates(v)
    c  = cinf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	
	ik = (1e-3) * gkbar * min(cai/250(mM),1) * c * (v-ek)
}


DERIVATIVE states { 

    rates(v)
    c' = (cinf-c)/tauc
}


PROCEDURE rates(v(mV)) { LOCAL a, b

    if (v<=-10) {
    
        a = 2(/ms) / 37.95 * ( exp( ( v + 50 ) / 11(mV) - ( v + 53.5 ) / 27(mV) ) )
        b = 2(/ms) * exp( ( - v - 53.5 ) / 27(mV) ) - a
    
    }else{
    
        a =  2(/ms) * exp( ( - v - 53.5 ) / 27(mV) )
        b = 0(/ms)
    
    }
    
    cinf = a/(a+b)
    tauc = 1.0/(a+b)    
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