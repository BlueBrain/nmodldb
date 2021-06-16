: na.mod codes a voltage-gated Na+ channel.
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
        SUFFIX na
        USEION na READ ena WRITE ina
        RANGE gnabar, gna, ina
        GLOBAL hinf, minf, htau, mtau
		GLOBAL aa4,bb4,cc4,dd4,ee4,ff4,gg4,hh4,ii4,jj4,kk4,ll4,mm4,nn4,oo4,pp4,qq4,rr4,ss4,tt4
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 20 (degC)  
        dt (ms)
        ena (mV)
        gnabar =  0.07958 (mho/cm2) <0,1e9>
		aa4= 39 <0,1e3>  
		bb4= 8 <0,1e3>
		cc4= 63 <0,1e3>  
		dd4= 7 <0,1e3>
		ee4= 70 <0,1e3>
		ff4= 5 <0,1e3>   
		gg4= 110 <0,1e3>
		hh4= 16 <0,1e3>
		ii4= 200 <0,1e3>
		jj4= 50 <0,1e3>
		kk4= 10 <0,1e3>
		ll4= 0.01 <0,1e3>
		mm4= 110 <0,1e3>
		nn4= 3 <0,1e3>
		oo4= 110 <0,1e3>
		pp4= 20 <0,1e3>
		qq4= 10 <0,1e3>
		rr4= 80 <0,1e3>
		ss4= 13 <0,1e3>
        tt4 = 0.1   <0,1e3>  
}

STATE {
        m h
}

ASSIGNED {
    ina (mA/cm2)
    gna (mho/cm2)
    minf hinf
    mtau (ms) htau (ms)
    }

LOCAL mexp, hexp

BREAKPOINT {
	SOLVE states
    gna = gnabar*(m^3)*h
    ina = gna*(v - ena)
}

UNITSOFF

INITIAL {
    trates(v)
    m = minf
    h = hinf
}

PROCEDURE states() {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	m = m + mexp*(minf-m)
	h = h + hexp*(hinf-h)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

    minf = 1 / (1+exp(-(v + aa4) / bb4))
    hinf = 1 / (1+exp((v + cc4) / dd4))

    mtau =  (ee4 / (ff4*exp((v+gg4) / hh4) + ii4*exp(-(v+jj4) / kk4))) + ll4
    htau =  (mm4 / (nn4*exp((v+oo4) / pp4) + qq4*exp(-(v+rr4) / ss4))) + tt4
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
	TABLE minf, mexp, hinf, hexp
	DEPEND dt, celsius FROM -150 TO 150 WITH 300

	q10 = 3^((celsius - 20)/10) 
    rates(v)    
	tinc = -dt * q10
	mexp = 1 - exp(tinc/mtau)
	hexp = 1 - exp(tinc/htau)
	}


UNITSON
