: kcna.mod codes low-threshold K+ channel Kv1.1 encoded in zebrafish kcna1a.
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
        SUFFIX kcna
        USEION k READ ek WRITE ik
        RANGE gkcnabar, gkcna, ik
        GLOBAL winf2, zinf2, wtau2, ztau2
		GLOBAL aa2,bb2,cc2,dd2,ee2,ff2,gg2,hh2,ii2,jj2,kk2,ll2,mm2,nn2,oo2,pp2,qq2,rr2,ss2
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 20 (degC)  
        dt (ms)
        ek = -90 (mV)
        gkcnabar = 0.01592 (mho/cm2) <0,1e9>
		aa2= 35 <0,1e3>
		bb2= 7 <0,1e3>
		cc2= 0.25 <0,1e3>
		dd2= 71 <0,1e3>
		ee2= 150 <0,1e3>
		ff2= 120 <0,1e3>
		gg2= 6 <0,1e3>
		hh2= 80 <0,1e3>
		ii2= 24 <0,1e3>
		jj2= 70 <0,1e3>
		kk2= 59 <0,1e3>
		ll2= 7 <0,1e3>
		mm2= 1.3 <0,1e3>
		nn2= 1000 <0,1e4>
		oo2= 60 <0,1e3>
		pp2= 20 <0,1e3>
		qq2= 60 <0,1e3>
		rr2= 8 <0,1e3>
		ss2= 50 <0,1e3>
        zss2 = 0.9   <0,1e3>   : steady state inactivation of glt
}

STATE {
        w2 z2
}

ASSIGNED {
    ik (mA/cm2) 
    gkcna (mho/cm2)
    winf2 zinf2
    wtau2 (ms) ztau2 (ms)
    }

LOCAL wexp2, zexp2

BREAKPOINT {
	SOLVE states
    
	gkcna = gkcnabar*(w2^4)*z2
    ik = gkcna*(v - ek)

}

UNITSOFF

INITIAL {
    trates(v)
    w2 = winf2
    z2 = zinf2
}

PROCEDURE states() {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	w2 = w2 + wexp2*(winf2-w2)
	z2 = z2 + zexp2*(zinf2-z2)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
  
    winf2 = (1 / (1 + exp(-(v + aa2) / bb2)))^cc2
    zinf2 = zss2 + ((1-zss2) / (1 + exp((v + dd2) / ee2)))

    wtau2 =  (ff2 / (gg2*exp((v+hh2) / ii2) + jj2*exp(-(v+kk2) / ll2))) + mm2
    ztau2 =  (nn2 / (exp((v+oo2) / pp2) + exp(-(v+qq2) / rr2))) + ss2
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
	TABLE winf2, wexp2, zinf2, zexp2
	DEPEND dt, celsius FROM -150 TO 150 WITH 300
	
    q10 = 3^((celsius - 20)/10) 
    rates(v)    
	tinc = -dt * q10
	wexp2 = 1 - exp(tinc/wtau2)
	zexp2 = 1 - exp(tinc/ztau2)
	}


UNITSON
