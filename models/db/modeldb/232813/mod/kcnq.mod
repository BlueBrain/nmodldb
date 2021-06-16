: kcnq.mod codes low-threshold K+ channel Kv7.4/KCNQ4 in zebrafish.
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
        SUFFIX kcnq
        USEION k READ ek WRITE ik
        RANGE gkcnqbar, gkcnq, ik
        GLOBAL winf1, zinf1, wtau1, ztau1
		GLOBAL aa0,bb0,cc0,dd0,ee0,ff0,gg0,hh0,ii0,jj0,kk0,ll0,mm0,nn0,oo0,pp0,qq0,rr0,ss0
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 20 (degC)  : model is defined on measurements made at room temp in Xenopus oocytes system
        dt (ms)
        ek = -90 (mV)
        gkcnqbar = 0.01592 (mho/cm2) <0,1e9>
		aa0= 29 <0,100>
		bb0= 10 <0,100>
		cc0= 0.25 <0,10>
		dd0= 71 <0,100>
		ee0= 150 <0,1e3>
		ff0= 320 <0,1e4>
		gg0= 4 <0,100>
		hh0= 0 <0,1e2>
		ii0= 70 <0,1e3>
		jj0= 10 <0,100>
		kk0= 63 <0,100>
		ll0= 6 <0,100>
		mm0= 10 <0,1e3>
		nn0= 1000 <0,1e4>
		oo0= 60 <0,100>
		pp0= 20 <0,1e3>
		qq0= 60 <0,100>
		rr0= 8 <0,100>
		ss0= 50 <0,100>
        zss0 = 0.9   <0,10>   : steady state inactivation of glt
}

STATE {
        w1 z1
}

ASSIGNED {
    ik (mA/cm2) 
    gkcnq (mho/cm2)
    winf1 zinf1
    wtau1 (ms) ztau1 (ms)
    }

LOCAL wexp1, zexp1

BREAKPOINT {
	SOLVE states
    
	gkcnq = gkcnqbar*(w1^4)*z1
    ik = gkcnq*(v - ek)

}

UNITSOFF

INITIAL {
    trates(v)
    w1 = winf1
    z1 = zinf1
}

PROCEDURE states() {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	w1 = w1 + wexp1*(winf1-w1)
	z1 = z1 + zexp1*(zinf1-z1)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

    winf1 = (1 / (1 + exp(-(v + aa0) / bb0)))^cc0
    zinf1 = zss0 + ((1-zss0) / (1 + exp((v + dd0) / ee0)))

    wtau1 =  (ff0 / (gg0*exp((v+hh0) / ii0) + jj0*exp(-(v+kk0) / ll0))) + mm0
    ztau1 =  (nn0 / (exp((v+oo0) / pp0) + exp(-(v+qq0) / rr0))) + ss0
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
	TABLE winf1, wexp1, zinf1, zexp1
	DEPEND dt, celsius FROM -150 TO 150 WITH 300
	
	q10 = 3^((celsius - 20)/10) 
    rates(v)    
	tinc = -dt * q10
	wexp1 = 1 - exp(tinc/wtau1)
	zexp1 = 1 - exp(tinc/ztau1)
	}

UNITSON
