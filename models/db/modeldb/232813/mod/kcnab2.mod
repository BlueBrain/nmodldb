: kcnab2.mod codes low-threshold K+ channel Kv1.1 coexpressed with Kvbeta2b/kcnab2b in zebrafish.
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
        SUFFIX kcnab2
        USEION k READ ek WRITE ik
        RANGE gkcnab2bar, gkcnab2, ik
        GLOBAL winf4, zinf4, wtau4, ztau4
		GLOBAL aa4,bb4,cc4,dd4,ee4,ff4,gg4,hh4,ii4,jj4,kk4,ll4,mm4,nn4,oo4,pp4,qq4,rr4,ss4
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 20 (degC)  
        dt (ms)
        ek = -90 (mV)
        gkcnab2bar = 0.01592 (mho/cm2) <0,1e9>
		aa4= 35 <0,1e3>
		bb4= 9 <0,1e3>
		cc4= 0.25 <0,1e3>
		dd4= 71 <0,1e3>
		ee4= 150 <0,1e3>
		ff4= 175 <0,1e3>
		gg4= 6 <0,1e3>
		hh4= 85 <0,1e3>
		ii4= 20 <0,1e3>
		jj4= 80 <0,1e3>
		kk4= 62 <0,1e3>
		ll4= 6 <0,1e3>
		mm4= 1 <0,1e3>
		nn4= 1000 <0,1e4>
		oo4= 60 <0,1e3>
		pp4= 20 <0,1e3>
		qq4= 60 <0,1e3>
		rr4= 8 <0,1e3>
		ss4= 50 <0,1e3>
        zss4 = 0.9   <0,1e3>   : steady state inactivation of glt
}

STATE {
        w4 z4
}

ASSIGNED {
    ik (mA/cm2) 
    gkcnab2 (mho/cm2)
    winf4 zinf4
    wtau4 (ms) ztau4 (ms)
    }

LOCAL wexp4, zexp4

BREAKPOINT {
	SOLVE states
    
	gkcnab2 = gkcnab2bar*(w4^4)*z4
    ik = gkcnab2*(v - ek)

}

UNITSOFF

INITIAL {
    trates(v)
    w4 = winf4
    z4 = zinf4
}

PROCEDURE states() {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	w4 = w4 + wexp4*(winf4-w4)
	z4 = z4 + zexp4*(zinf4-z4)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
  
    winf4 = (1 / (1 + exp(-(v + aa4) / bb4)))^cc4
    zinf4 = zss4 + ((1-zss4) / (1 + exp((v + dd4) / ee4)))

    wtau4 =  (ff4 / (gg4*exp((v+hh4) / ii4) + jj4*exp(-(v+kk4) / ll4))) + mm4
    ztau4 =  (nn4 / (exp((v+oo4) / pp4) + exp(-(v+qq4) / rr4))) + ss4
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
	TABLE winf4, wexp4, zinf4, zexp4
	DEPEND dt, celsius FROM -150 TO 150 WITH 300
	
    q10 = 3^((celsius - 20)/10) 
    rates(v)    
	tinc = -dt * q10
	wexp4 = 1 - exp(tinc/wtau4)
	zexp4 = 1 - exp(tinc/ztau4)
	}


UNITSON
