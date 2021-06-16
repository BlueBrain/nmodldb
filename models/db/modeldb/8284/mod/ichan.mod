TITLE ichan.mod  
 
COMMENT
Chen K, Aradi I, Thon N, Eghbal-Ahmadi M, Baram TZ, Soltesz I: Persistently modified
h-channels after complex febrile seizures convert the seizure-induced enhancement of
inhibition to hyperexcitability. Nature Medicine, 7(3) pp. 331-337, 2001.
(modeling by Ildiko Aradi, iaradi@uci.edu)
sodium and potassium channel kinetics for Control and HT cases 
ENDCOMMENT
 
UNITS {
        (mA) =(milliamp)
        (mV) =(millivolt)
        (uF) = (microfarad)
	(molar) = (1/liter)
	(nA) = (nanoamp)
	(mM) = (millimolar)
	(um) = (micron)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
? interface 
NEURON { 
SUFFIX ichan 
USEION nat READ enat WRITE inat VALENCE 1
USEION kf READ ekf WRITE ikf  VALENCE 1
NONSPECIFIC_CURRENT il 
RANGE  gnat, gkf
RANGE gnatbar, gkfbar
RANGE gl, el
RANGE minf, mtau, hinf, htau, nfinf, nftau
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
PARAMETER {
        v (mV) 
        celsius = 6.3 (degC)
        dt (ms) 
        enat  (mV)
	gnatbar (mho/cm2)   
        ekf  (mV)
	gkfbar (mho/cm2)
	gl (mho/cm2)    
 	el (mV)
}
 
STATE {
	m h nf
}
 
ASSIGNED {
         
        gnat (mho/cm2) 
        gkf (mho/cm2)


        inat (mA/cm2)
        ikf (mA/cm2)

	il (mA/cm2)

	minf hinf nfinf
 	mtau (ms) htau (ms) nftau (ms)
	mexp hexp nfexp
} 

? currents
BREAKPOINT {
	SOLVE states
        gnat = gnatbar*m*m*m*h  
        inat = gnat*(v - enat)
        gkf = gkfbar*nf*nf*nf*nf
        ikf = gkf*(v-ekf)
	il = gl*(v-el)
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	
	m = minf
	h = hinf
      nf = nfinf
	
	VERBATIM
	return 0;
	ENDVERBATIM
}

? states
PROCEDURE states() {	:Computes state variables m, h, and n 
        trates(v)	:      at the current v and dt.
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
        nf = nf + nfexp*(nfinf-nf)
        VERBATIM
        return 0;
        ENDVERBATIM
}
 
LOCAL q10

? rates
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
       q10 = 3^((celsius - 6.3)/10)
                :"m" sodium activation system
	alpha = -0.3*vtrap((v+65-25),-5)
	beta = 0.3*vtrap((v+65-53),5)
	sum = alpha+beta        
	mtau = 1/sum      minf = alpha/sum
                :"h" sodium inactivation system
	alpha = 0.23/exp((v+65-3)/20)
	beta = 3.33/(1+exp((v+65-55.5)/-10))
	sum = alpha+beta
	htau = 1/sum 
        hinf = alpha/sum 
             :"nf" fKDR activation system
        alpha = -0.07*vtrap((v+65-47),-6)
	beta = 0.264/exp((v+65-22)/40)
	sum = alpha+beta        
	nftau = 1/sum      nfinf = alpha/sum
	
}
 
PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
        TABLE minf, mexp, hinf, hexp, nfinf, nfexp, mtau, htau, nftau
	DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
	rates(v)	: not consistently executed from here if usetable_hh == 1
		: so don't expect the tau values to be tracking along with
		: the inf values in hoc

	       tinc = -dt * q10
        mexp = 1 - exp(tinc/mtau)
        hexp = 1 - exp(tinc/htau)
	nfexp = 1 - exp(tinc/nftau)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON

