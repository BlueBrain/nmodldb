TITLE kdr.mod
 
COMMENT
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
NEURON {
        SUFFIX kdr
        USEION k READ ek WRITE ik
        RANGE gskbar,   gsk
        GLOBAL ninf,  ntau
}
 
PARAMETER {
        gskbar = .003 (S/cm2)	<0,1e9>
}
 
STATE {
         n
}
 
ASSIGNED {
        v (mV)
        :celsius (degC)
        ek (mV)
	gsk (S/cm2)
        ik (mA/cm2)
         ninf
	 ntau (ms)
}
 
LOCAL  nexp        
 

BREAKPOINT {
        SOLVE states METHOD cnexp
        gsk = gskbar*n*n*n*n
	ik = gsk*(v - ek)      
}
 
 
INITIAL {
	rates(v)
	n = ninf
}

DERIVATIVE states {  
        rates(v)
         n' = (ninf-n)/ntau
}
 
LOCAL q10


PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
       : TABLE  ninf, ntau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
	q10=1        
	:q10 = 3^((celsius - 6.3)/10)
                :"n" potassium activation system
        alpha = -.028*vtrap((v+70-35),-6) 
        beta = .1056*exp((v+70-10)/40)
	sum = alpha + beta
        ntau = 1/(q10*sum)
        ninf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
