TITLE HH_Na4 channel for LGMD
: Altered by RBD

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    THREADSAFE
    SUFFIX HH_NaF
    USEION na READ ena WRITE ina
    RANGE gmax,g
}

PARAMETER {
    gmax= 0.045 (mho/cm2)

	kam=0.25
	Aam=1.3
	dam=-40
	kbm=0.1
	Abm=9
	dbm=-60

	kah=0.08
	Aah=12.0
	dah=-90
	kbh=0.15
	Abh=6.0
	dbh=-40
}

ASSIGNED { 
    v (mV)
    ena (mV)
    ina (mA/cm2)
    am (/ms)
    bm (/ms)
    ah (/ms)
    bh (/ms)
    g  (mho/cm2)
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
    TABLE am, bm, ah, bh
          FROM -100 TO 100 WITH 2000

    if (v > (dbm-.1) && v < (dbm+0.1)) {
    	am = Aam
	} else {
		am = kam*-Aam*(dbm-v) / (1 - exp((dbm-v)*kam))
	}
	bm = Abm*exp( (dam-v)*kbm)
	
	ah = Aah*exp( (dah-v)*kah)
	bh = Abh/(exp( (dbh-v)* kbh) + 1)

}

UNITSON


