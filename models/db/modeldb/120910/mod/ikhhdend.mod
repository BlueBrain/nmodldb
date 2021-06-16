TITLE ikhhdend.mod   squid potassium channels
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX ikhhdend
        USEION k READ ek WRITE ik
        RANGE gkbar
        GLOBAL ninf, nexp
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius (degC)
        dt (ms)
	gkbar = 0.001(mho/cm2)        
}
 
STATE {
        n c
}
 
ASSIGNED {
	ik (mA/cm2)
        ek  (mV)       
        ninf nexp 
}
 
BREAKPOINT {
        SOLVE states
        ik = gkbar*n*n*(v - ek)      
}
 
UNITSOFF
 
INITIAL {
     rates(v)
     n = ninf     
}
PROCEDURE states() {  :Computes state variables n 
        rates(v)      :             at the current v and dt.

        n = n + nexp*(ninf-n)
}
 
PROCEDURE rates(v) {:Computes rate and o
         : ther constants at current v.
         : Call once from HOC to 
         : initialize inf at resting v.
     LOCAL  q10, tinc, alpha, beta, sum
     TABLE ninf, nexp  
	DEPEND dt, celsius

FROM -100 TO 100 WITH 200
        q10 = 2.3^((celsius - 20)/10)
        tinc = -dt * q10

                :"n" potassium activation system
        alpha = .024*vtrap(-(v-17),8) 
        beta = 0.2*exp(-(v+48)/35)
        sum = alpha + beta
        ninf = alpha/sum
        nexp = 1 - exp(tinc*sum)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
}
}
 
UNITSON

 



