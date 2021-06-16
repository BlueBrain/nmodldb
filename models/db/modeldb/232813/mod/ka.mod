: ka.mod codes an A-type voltage-gated K+ channel.
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
        SUFFIX ka
        USEION k READ ek WRITE ik
        RANGE gkabar, gka, ik
        GLOBAL ainf, binf, cinf, atau, btau, ctau
		GLOBAL aa3,bb3,cc3,dd3,ee3,ff3,gg3,hh3,ii3,jj3,kk3,ll3,mm3,nn3,oo3,pp3,qq3,rr3,ss3,tt3,uu3,vv3,ww3,xx3,yy3,zz3,aaa3,bbb3,ccc3
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 20 (degC)  
        dt (ms)
        ek = -90 (mV)
        gkabar = 0.00477 (mho/cm2) <0,1e9>
		aa3=31 <0,1e3>
		bb3=6 <0,1e3>
		cc3=0.25 <0,1e3>
		dd3=66 <0,1e3>
		ee3=7 <0,1e3>
		ff3=0.5 <0,1e3>
		gg3=66 <0,1e3>
		hh3=7 <0,1e3>
		ii3=0.5 <0,1e3>
		jj3=100 <0,1e3>
		kk3=7 <0,1e3>
		ll3=60 <0,1e3>
		mm3=14 <0,1e3>
		nn3=29 <0,1e3>
		oo3=60 <0,1e3>
		pp3=24 <0,1e3>
		qq3=0.1 <0,1e3>
		rr3=1000 <0,1e4>
		ss3=14 <0,1e3>
		tt3=60 <0,1e3>
		uu3=27 <0,1e3>
		vv3=29 <0,1e3>
		ww3=60 <0,1e3>
		xx3=24 <0,1e3>
		yy3=1 <0,1e3>
		zz3=90 <0,1e3>
		aaa3=66 <0,1e3>
		bbb3=17	<0,1e3>
		ccc3=10 <0,1e3>
}

STATE {
        a b c
}

ASSIGNED {
    ik (mA/cm2) 
    gka (mho/cm2)
    ainf binf cinf
    atau (ms) btau (ms) ctau (ms)
    }

LOCAL aexp, bexp, cexp

BREAKPOINT {
	SOLVE states
    
	gka = gkabar*(a^4)*b*c
    ik = gka*(v - ek)

}

UNITSOFF

INITIAL {
    trates(v)
    a = ainf
    b = binf
    c = cinf
}

PROCEDURE states() {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	a = a + aexp*(ainf-a)
	b = b + bexp*(binf-b)
	c = c + cexp*(cinf-c)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

    ainf = (1 / (1 + exp(-1*(v + aa3) / bb3)))^cc3
    binf = 1 / (1 + exp((v + dd3) / ee3))^ff3
    cinf = 1 / (1 + exp((v + gg3) / hh3))^ii3

    atau =  (jj3 / (kk3*exp((v+ll3) / mm3) + nn3*exp(-(v+oo3) / pp3))) + qq3
    btau =  (rr3 / (ss3*exp((v+tt3) / uu3) + vv3*exp(-(v+ww3) / xx3))) + yy3
    ctau = (zz3 / (1 + exp((-aaa3-v) / bbb3))) + ccc3
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
	TABLE ainf, aexp, binf, bexp, cinf, cexp
	DEPEND dt, celsius FROM -150 TO 150 WITH 300
    
	q10 = 3^((celsius - 20)/10) 
    rates(v)    
	tinc = -dt * q10
	aexp = 1 - exp(tinc/atau)
	bexp = 1 - exp(tinc/btau)
	cexp = 1 - exp(tinc/ctau)
	}


UNITSON
