: kht.mod codes a high-threshold K+ channel.
: Default parameters of a H-H equation are fitted to our experimental data
: by using our channel generator. 
:
: Takaki Watanabe
: wtakaki@m.u-tokyo.ac.jp

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX kht
        USEION k READ ek WRITE ik
        RANGE gkhtbar, gkht, ik
        GLOBAL ninf, pinf, ntau, ptau
		GLOBAL aa1,bb1,cc1,dd1,ee1,ff1,gg1,hh1,ii1,jj1,kk1,ll1,mm1,nn1,oo1,pp1,qq1,rr1,ss1,tt1,uu1
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 20 (degC)  
        dt (ms)
        ek = -90 (mV)
        gkhtbar = 0.01592 (mho/cm2) <0,1e9>
		aa1= 15 <0,100>
		bb1= 5 <0,100>
		cc1= 1 <0,10>
		dd1= 23 <0,100>
		ee1= 6 <0,100>
		ff1= 100 <0,1e3>
		gg1= 11 <0,100>
		hh1= 60 <0,100>
		ii1= 24 <0,100>
		jj1= 21 <0,100>
		kk1= 60 <0,100>
		ll1= 23 <0,100>
		mm1= 0.5 <0,100>
		nn1= 100 <0,1e3>
		oo1= 4 <0,100>
		pp1= 60 <0,100>
		qq1= 21 <0,100>
		rr1= 5 <0,100>
		ss1= 10 <0,100>
		tt1= 18 <0,100>
		uu1= 0.5 <0,100>
		nf = 0.85 <0,1> :proportion of n vs p kinetics
}

STATE {
        n p
}

ASSIGNED {
    ik (mA/cm) 
    gkht (mho/cm2)
    pinf ninf
    ptau (ms) ntau (ms)
    }

LOCAL nexp, pexp

BREAKPOINT {
	SOLVE states
    
	gkht = gkhtbar*(nf*(n^2) + (1-nf)*p)
    ik = gkht*(v - ek)

}

UNITSOFF

INITIAL {
    trates(v)
    p = pinf
    n = ninf
}

PROCEDURE states() {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	n = n + nexp*(ninf-n)
	p = p + pexp*(pinf-p)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
    
    ninf =   (1 + exp(-(v + aa1) / bb1))^-cc1
    pinf =  1 / (1 + exp(-(v + dd1) / ee1))

	ntau = (ff1 / (gg1*exp((v+hh1) / ii1) + jj1*exp(-(v+kk1) / ll1))) + mm1
    ptau = (nn1 / (oo1*exp((v+pp1) / qq1) + rr1*exp(-(v+ss1) / tt1))) + uu1
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
	TABLE ninf, nexp, pinf, pexp
	DEPEND dt, celsius FROM -150 TO 150 WITH 300
	
    q10 = 3^((celsius - 20)/10) 
    rates(v)    
	tinc = -dt * q10
	nexp = 1 - exp(tinc/ntau)
	pexp = 1 - exp(tinc/ptau)
	}


UNITSON
