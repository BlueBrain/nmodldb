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

	kam=0.075	(/mV)
	Aam=26		(/ms)
	dam=-10		(mV)
	kbm=0.07	(/mV)
	Abm=7.4		(/ms)
	dbm=-42		(mV)

	kah=0.15	(/mV)
	Aah=2.1		(/ms)
	dah=-60		(mV)
	kbh=0.14	(/mV)
	Abh=8.5		(/ms)
	dbh=-30		(mV)
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

:UNITSOFF

PROCEDURE settables(v (mV)) {
    TABLE am, bm, ah, bh DEPEND dam, dbm, dah, dbh
          FROM -100 TO 50 WITH 750

    if (v > (dam-0.01) && v < (dam+0.01)) {
    	am = Aam
	} else {
		am = kam*-Aam*(dam-v) / (1 - exp((dam-v)*kam))
	}
	bm = Abm*exp( (dbm-v)*kbm)
	
	ah = Aah*exp( (dah-v)*kah)
	bh = Abh/(exp( (dbh-v)* kbh) + 1)

}

:UNITSON


