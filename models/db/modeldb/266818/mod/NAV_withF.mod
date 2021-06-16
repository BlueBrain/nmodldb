TITLE Mouse Ventricular Myocyte Model

COMMENT
16-state Nav Markov model (+ FHF2)
ENDCOMMENT

NEURON {
	SUFFIX NAV_withF
	USEION na READ ena WRITE ina
	RANGE gnabar, ina, g, Oapp, Capp, Iapp
	RANGE alfa, beta, gamma, delta, Con, Coff, Oon, Ooff
	RANGE Aalfa, Valfa, Abeta, Vbeta, Agamma, Adelta, ACon, ACoff, AOon, AOoff, Vshift
	RANGE n1, n2, n3, n4, n5, n6, s1, s2, s3, s4, s5, s6, V_threshold, use_threshold
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	Vshift = -54    (mV)
	celsius = 25  	(degC)
	ena = 87.39		(mV)
	gnabar = 0.036	(mho/cm2)
	Aalfa = 2.44375 ( /ms)
	Valfa = 9 ( /mV) 
	Abeta = 0.01325  ( /ms)
	Vbeta = 9 ( /mV)
	Agamma = 150 ( /ms)
	Adelta = 40  ( /ms)
	ACon = 0.004    ( /ms)
	ACoff = 0.05     ( /ms)
	AOon = 0.85     ( /ms)
	AOoff = 0.0005   ( /ms)
	n1 = 100
	n2 = 20
	n3 = 20
	n4 = 3
	n5 = 1.5
	n6 = 0.75
	s1 = 0.5
	s2 = 0.5
	s3 = 0.5
	s4 = 1.5
	s5 = 1.5
	s6 = 1.5
	V_threshold = -85 (mV)
	use_threshold = 0	
	
}

ASSIGNED {
	ina  (mA/cm2)
	g   (mho/cm2)
	Oapp
	Iapp
	Capp
	gamma
	delta
	Con
	Coff
	Oon
	Ooff
	a
	b
	Q10gate
	Q10cond
	
}

STATE {
	C1
	C2
	C3
	C4
	C5
	C6
	C7
	O
	I1
	I2
	I3
	I4
	I5
	I6
	I7
	I8

}

INITIAL {
	: C1=1
	: C2=0
	: C3=0
	: C4=0
	: C5=0
	: C6=0
	: C7=0
	: O=0
	: I1=0
	: I2=0
	: I3=0
	: I4=0
	: I5=0
	: I6=0
	: I7=0
	: I8=0
	Q10gate =2.3^((celsius-32.76(degC))/10 (degC))
	Q10cond =1.4^((celsius-25 (degC))/10 (degC))
	gamma = Q10gate * Agamma
	delta = Q10gate * Adelta
	Con = Q10gate * ACon
	Coff = Q10gate * ACoff
	Oon = Q10gate * AOon
	Ooff = Q10gate * AOoff
	a = (Oon/Con)^0.16667
	b = (Ooff/Coff)^0.16667
	SOLVE kstates STEADYSTATE sparse	
	
    }
    
BREAKPOINT {	    
    if ( use_threshold ) {
	if (v < V_threshold) {
	    delta = 1e10
	    gamma = 1e-10
	    :printf("%f\t",v)
	} else {
	    delta = Q10gate * Adelta
	    gamma = Q10gate * Agamma
	}
    }
 
    SOLVE kstates METHOD sparse
    g = Q10cond * gnabar * O	      	: (mho/cm2)
    ina = g * (v - ena)  	: (mA/cm2)
    Oapp = O*(ACon+ACoff)/ACoff
    Capp = (C1+C2+C3+C4+C5+C6+C7)*(ACon+ACoff)/ACoff
    Iapp = 1-Oapp-Capp
}

FUNCTION alfa(v(mV))(/ms){ 
	alfa = Q10gate*Aalfa*exp((v-Vshift)/Valfa) 
}

FUNCTION beta(v(mV))(/ms){ 
	beta = Q10gate*Abeta*exp((-v+Vshift)/Vbeta) 
}

KINETIC kstates {
	: 1st row
	~ C1 <-> C2 (n1*alfa(v),n6*beta(v))
	~ C2 <-> C3 (n2*alfa(v),n5*beta(v))
	~ C3 <-> C4 (n3*alfa(v),n4*beta(v))
	~ C4 <-> C5 (n4*alfa(v),n3*beta(v))
	~ C5 <-> C6 (n5*alfa(v),n2*beta(v))	
	~ C6 <-> C7 (n6*alfa(v),n1*beta(v))	
	~ C7 <-> O  (gamma,delta)
	
	: 2nd row
	~ I1 <-> I2	(n1*alfa(v)*a^s1,n6*beta(v)*b^s1)
	~ I2 <-> I3	(n2*alfa(v)*a^s2,n5*beta(v)*b^s2)
	~ I3 <-> I4	(n3*alfa(v)*a^s3,n4*beta(v)*b^s3)
	~ I4 <-> I5 (n4*alfa(v)*a^s4,n3*beta(v)*b^s4)
	~ I5 <-> I6 (n5*alfa(v)*a^s5,n2*beta(v)*b^s5)	
	~ I6 <-> I7 (n6*alfa(v)*a^s6,n1*beta(v)*b^s6)	
	~ I7 <-> I8 (gamma,delta)
	
		
	: connect 1st row with 2nd row
	~ C1 <-> I1 (Con,Coff)
	~ C2 <-> I2 (Con*a^s1,Coff*b^s1)
	~ C3 <-> I3 (Con*a^(s1+s2),Coff*b^(s1+s2))
	~ C4 <-> I4 (Con*a^(s1+s2+s3),Coff*b^(s1+s2+s3))
	~ C5 <-> I5 (Con*a^(s1+s2+s3+s4),Coff*b^(s1+s2+s3+s4))
	~ C6 <-> I6 (Con*a^(s1+s2+s3+s4+s5),Coff*b^(s1+s2+s3+s4+s5))	
	~ C7 <-> I7 (Con*a^(s1+s2+s3+s4+s5+s6),Coff*b^(s1+s2+s3+s4+s5+s6))	
	~  O <-> I8 (Oon,Ooff)
	
	
		
	CONSERVE C1+C2+C3+C4+C5+C6+C7+O+I1+I2+I3+I4+I5+I6+I7+I8=1

	
}



