TITLE cal :High-threshold Ca current

NEURON { 
	SUFFIX cal
	USEION ca READ eca WRITE ica
	RANGE  gcal, ica, qca
}

UNITS { 
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
        FARADAY	= 96485.309 (coul/mole)
	PI	= (pi) (1) 
}

INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER { 
	gcal = 0.0 	(mho/cm2)
}

ASSIGNED { 
	ica 	(mA/cm2) 
	v	(mV)
        eca 	(mV)
	diam	(um)
}
 
STATE { m c qca }

BREAKPOINT { 
	SOLVE castate METHOD sparse
	ica = gcal*m*m*(v-eca) 
}
 
INITIAL { 
	m=m_inf(v)
	c=1-m
	ica = gcal*m*m*(v-eca)
	qca=0
}

LOCAL a1, a2

KINETIC castate {
	a1 = a_m(v)
	a2 = a_c(v)
	~ c <-> m	(a1, a2)
	CONSERVE m + c = 1
        
	COMPARTMENT diam*diam*PI/4 { qca }
	~ qca << (-ica*diam*PI*(1e4)/FARADAY )
}

FUNCTION a_m(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	a_m = 1.6/(1+exp(-0.072*(v-65)))
}

FUNCTION a_c(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	a_c = 0.02*(v-51.1)/(exp((v-51.1)/5)-1)
}

FUNCTION m_inf(v(mV)) {
        m_inf = a_m(v)/(a_m(v)+a_c(v))
}
