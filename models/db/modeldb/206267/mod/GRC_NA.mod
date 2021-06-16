TITLE Cerebellum Granule Cell Model

COMMENT
basato sul modello di Raman a 13 stati. genera corrente di sodio transiente, persistente e risorgente
with Long-Term Inactivation States L3,L4,L5,L6. and Vshift = -10mV
ENDCOMMENT

NEURON {
	SUFFIX GRC_NA
	USEION na READ ena WRITE ina
	RANGE gnabar, ina, g
	RANGE alfa, beta, gamma, delta, epsilon, teta, Con, Coff, Oon, Ooff, Lon, Loff
	RANGE Aalfa, Valfa, Abeta, Vbeta, Ateta, Vteta, Agamma, Adelta, Aepsilon, ACon, ACoff, AOon, AOoff, ALon, ALoff, Vshift
	RANGE n1, n2, n3, n4, c, d, V_threshold, use_threshold
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	Vshift = -10    (mV)
	celsius = 20  	(degC)
	ena = 87.39		(mV)
	gnabar = 0.013	(mho/cm2)
	Aalfa = 353.91 ( /ms)
	Valfa = 13.99 ( /mV) 
	Abeta = 1.272  ( /ms)
	Vbeta = 13.99 ( /mV)
	Agamma = 150 ( /ms)
	Adelta = 40  ( /ms)
	Aepsilon = 1.75 ( /ms)
	Ateta = 0.0201 ( /ms)
	Vteta = 25
	ACon = 0.5    ( /ms)
	ACoff = 0.5     ( /ms)
	AOon = 7.5     ( /ms)
	AOoff = 0.0005   ( /ms)
	ALon = 0   ( /ms)
	ALoff = 1000    ( /ms)
	n1 = 5.422
	n2 = 3.279
	n3 = 1.83
	n4 = 0.738
	c = 20
	d = 0.075
	V_threshold = -65 (mV)
	use_threshold = 1	
	
}

ASSIGNED {
	ina  (mA/cm2)
	g   (mho/cm2)
	
	gamma
	delta
	epsilon
	Con
	Coff
	Oon
	Ooff
	Lon
	Loff
	a
	b
	Q10
	
}

STATE {
	C1
	C2
	C3
	C4
	C5
	O
	OB
	I1
	I2
	I3
	I4
	I5
	I6
	L3
	L4
	L5
	L6
}


INITIAL {
	C1=1
	C2=0
	C3=0
	C4=0
	C5=0
	O=0
	OB=0
	I1=0
	I2=0
	I3=0
	I4=0
	I5=0
	I6=0
	L3=0
	L4=0
	L5=0
	L6=0
	Q10 =3^((celsius-20(degC))/10 (degC))
	gamma = Q10 * Agamma
	delta = Q10 * Adelta
	epsilon = Q10 * Aepsilon
	Con = Q10 * ACon
	Coff = Q10 * ACoff
	Oon = Q10 * AOon
	Ooff = Q10 * AOoff
	Lon = Q10 * ALon
	Loff = Q10 * ALoff
	a = (Oon/Con)^0.25
	b = (Ooff/Coff)^0.25
	
    }
    
	
BREAKPOINT {	    
    if ( use_threshold ) {
	if (v < V_threshold) {
	    delta = 1e10
	    gamma = 1e-10
	    :printf("%f\t",v)
	} else {
	    delta = Q10 * Adelta
	    gamma = Q10 * Agamma
	}
    }
    
    SOLVE kstates METHOD sparse
    g = gnabar * O	      	: (mho/cm2)
    ina = g * (v - ena)  	: (mA/cm2)
}


FUNCTION alfa(v(mV))(/ms){ 
	alfa = Q10*Aalfa*exp((v-Vshift)/Valfa) 
}

FUNCTION beta(v(mV))(/ms){ 
	beta = Q10*Abeta*exp((-v+Vshift)/Vbeta) 
}

FUNCTION teta(v(mV))(/ms){ 
	teta = Q10*Ateta*exp(-v/Vteta) 
}
 

KINETIC kstates {
	: 1 riga
	~ C1 <-> C2 (n1*alfa(v),n4*beta(v))
	~ C2 <-> C3 (n2*alfa(v),n3*beta(v))
	~ C3 <-> C4 (n3*alfa(v),n2*beta(v))
	~ C4 <-> C5 (n4*alfa(v),n1*beta(v))
	~ C5 <-> O  (gamma,delta)
	~  O <-> OB (epsilon,teta(v))
	
	: 2 riga
	~ I1 <-> I2	(n1*alfa(v)*a,n4*beta(v)*b)
	~ I2 <-> I3	(n2*alfa(v)*a,n3*beta(v)*b)
	~ I3 <-> I4	(n3*alfa(v)*a,n2*beta(v)*b)
	~ I4 <-> I5 (n4*alfa(v)*a,n1*beta(v)*b)
	~ I5 <-> I6 (gamma,delta)
	
	: 3 riga
	~ L3 <-> L4 (n3*alfa(v)*c,n2*alfa(v)*d)
	~ L4 <-> L5 (n4*alfa(v)*c,n1*alfa(v)*d)
	~ L5 <-> L6 (gamma,delta)
	
	: connette 1 riga con 2 riga
	~ C1 <-> I1 (Con,Coff)
	~ C2 <-> I2 (Con*a,Coff*b)
	~ C3 <-> I3 (Con*a^2,Coff*b^2)
	~ C4 <-> I4 (Con*a^3,Coff*b^3)
	~ C5 <-> I5 (Con*a^4,Coff*b^4)
	~  O <-> I6 (Oon,Ooff)
	
	: connette 1 riga con 3 riga
	~ C3 <-> L3 (Lon,Loff)
	~ C4 <-> L4 (Lon*c,Loff*d)
	~ C5 <-> L5 (Lon*c^2,Loff*d^2)
	~  O <-> L6 (Lon*c^2,Loff*d^2)
	
	CONSERVE C1+C2+C3+C4+C5+O+OB+I1+I2+I3+I4+I5+I6+L3+L4+L5+L6=1
}

