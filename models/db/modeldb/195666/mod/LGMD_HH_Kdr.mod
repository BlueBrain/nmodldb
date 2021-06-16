TITLE HH_Kdr channel for LGMD
: Altered by RBD

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    THREADSAFE
    SUFFIX HH_Kdr
    USEION k READ ek WRITE ik
    : Range variables can differ across neuron (non constants) and there values are accessible
    RANGE g, gmax, t1, i
    
    GLOBAL vhalf
}

PARAMETER {
    gmax= 0.008 (mho/cm2)
    
    vhalf=-39
	vn2=-60
	t1=35
	t2=0.15
	tns=-9
	zn=9
}

ASSIGNED { 
    v (mV)
    ek (mV)
    
    ik (mA/cm2)
	ninf (/ms)
    ntau (/ms)
    g (mS/cm2)
    i (mA/cm2)
}

STATE {
    n
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gmax*n^4
    i = g*(v-ek)
    ik = i
}

INITIAL {
    settables(v)
    n = ninf
}

DERIVATIVE states {  
    settables(v)    
    n' = (ninf - n)/ntau
}

UNITSOFF

PROCEDURE settables(v (mV)) {
    TABLE ninf, ntau DEPEND vhalf, t1
          FROM -100 TO 50 WITH 1500

	ninf = 1/(1 + exp((vhalf-v)/zn))
	ntau = t1/(1+exp((vn2-v)/tns))*ninf+t2

}

UNITSON


