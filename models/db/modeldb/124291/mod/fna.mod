TITLE fna.mod   
 
COMMENT
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
NEURON {
        SUFFIX fna
        USEION na READ ena WRITE ina
        RANGE gnabar,  gna, alpha, beta, alphah, betah, sum
        GLOBAL minf, hinf, mtau, htau
}
 
PARAMETER {
        gnabar = .12 (S/cm2)	<0,1e9>
        gl = .0003 (S/cm2)	<0,1e9>
        :celsius (degC)
}
 
STATE {
        m h 
}
 
ASSIGNED {
        ena (mV)
        v (mV)
        :celsius (degC)

	alpha
	beta
	alphah
	betah
	sum
	gna (S/cm2)
        ina (mA/cm2)
        minf hinf 
	mtau (ms) htau (ms) 
}
 
LOCAL mexp, hexp
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
	ina = gna*(v - ena)
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
}
 
LOCAL q10

PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        :LOCAL  alphah, betah
        :TABLE minf, mtau, hinf, htau,  DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        q10=1
	:q10 = 3^((celsius - 6.3)/10)
                :"m" sodium activation system
        alpha = -.3 * vtrap((v+70-25),-5)
        beta =  .3 * vtrap((v+70-53),5)
        sum = alpha + beta
	mtau = 1/(q10*sum)
        minf = alpha/sum

                :"h" sodium inactivation system
        alphah = .23 /exp((v+70-3)/20)
        betah = 3.33 / (1+exp((v+70-55.5)/(-10)) )
        sum = alphah + betah
	htau = 1/(q10*sum)
        hinf = alphah/sum
       
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*( 1-  x/y/2  )
        }else{
                vtrap = x/( exp(x/y) - 1 )
        }
}
 
UNITSON
