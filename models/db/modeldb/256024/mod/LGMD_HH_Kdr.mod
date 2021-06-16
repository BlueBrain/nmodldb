TITLE HH_Kdr channel for LGMD
: Altered by RBD

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}

NEURON {
    THREADSAFE
    SUFFIX HH_Kdr
    USEION k READ ek WRITE ik
    : Range variables can differ across neuron (non constants) and their values are accessible
    RANGE g, gmax, ntau, t1
    
    GLOBAL vhalf
}

PARAMETER {
    gmax= 0.008 (mho/cm2)
    
    vhalf=-39	(mV)
	vn2=-60		(mV)
	t1=35		(1)
	t2=0.15		(ms)
	tns=-9		(mV)
	zn=9		(mV)
}

ASSIGNED { 
    v (mV)
    ek (mV)
    
    ik (mA/cm2)
	ninf (1)
    ntau (ms)
    g (S/cm2)
}

STATE {
    n
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gmax*n^4
    ik = g*(v-ek)
}

INITIAL {
    settables(v)
    n = ninf
}

DERIVATIVE states {  
    settables(v)
    n' = (ninf - n)/(t1*ntau+t2)
}

:UNITSOFF

PROCEDURE settables(v (mV)) {
    TABLE ninf, ntau DEPEND vhalf, zn, vn2, tns
          FROM -80 TO 50 WITH 650

	ninf = 1/(1 + exp((vhalf-v)/zn))
	ntau = (1 (ms))/(1+exp((vn2-v)/tns))*ninf

}

:UNITSON


