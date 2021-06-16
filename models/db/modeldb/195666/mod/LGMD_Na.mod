TITLE Fast Na channel for LGMD
: Altered by RBD 7/2016

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
	(S) = (siemens)
}

NEURON {
    THREADSAFE
    SUFFIX Na
    USEION na READ ena WRITE ina
    GLOBAL dam, dbm
    RANGE gmax, g
}

PARAMETER {
: all values can be adjusted in hoc files
    gmax= 0.045 (mho/cm2)

	kam=0.075
	Aam=26
	dam=-42
	kbm=0.07
	Abm=7.4
	dbm=-8

	kah=0.15
	Aah=2.1
	dah=-60
	kbh=0.14
	Abh=8.5
	dbh=-30
}

ASSIGNED { 
    v (mV)
    ena (mV)
    
    ina (mA/cm2)
    am (/ms)
    bm (/ms)
    ah (/ms)
    bh (/ms)
	g (S/cm2)
}

STATE {
    m h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g  = gmax*m^3*h
    ina  = g*(v-ena)
}

INITIAL {
    settables(v)
    m = am/(am+bm)
    h = ah/(ah+bh)
}

DERIVATIVE states {  
    settables(v)      
    m' = ((am*(1-m)) - (bm*m))
    h' = ((ah*(1-h)) - (bh*h))
}

UNITSOFF

PROCEDURE settables(v (mV)) {
    TABLE am, bm, ah, bh DEPEND dam, dbm, dah, dbh
          FROM -100 TO 50 WITH 750

    if (v > (dbm-0.01) && v < (dbm+0.01)) {
    	am = Aam
	} else {
		am = kam*-Aam*(dbm-v) / (1 - exp((dbm-v)*kam))
	}
	bm = Abm*exp( (dam-v)*kbm)
	
	ah = Aah*exp( (dah-v)*kah)
	bh = Abh/(exp( (dbh-v)* kbh) + 1)

}

UNITSON


