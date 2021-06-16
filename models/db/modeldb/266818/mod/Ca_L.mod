TITLE Mouse Ventricular Myocyte Model

COMMENT
12-state L-type Cav Markov model based on Xiangming Lin's recordings of 129SvPas ventricular cardiomyocytes
ENDCOMMENT

NEURON {
	SUFFIX Ca_L
	USEION ca READ eca WRITE ica
	RANGE gcabar, ica, g
	RANGE alfa, beta, gamma, delta, Con, Coff, Oon, Ooff
	RANGE Aalfa, Valfa, Abeta, Vbeta, Agamma, Adelta, ACon, ACoff, AOon, AOoff, Vshift
	RANGE n1, n2
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	Vshift = -12    (mV)
	celsius = 25  	(degC)
	eca = 80		(mV)
	gcabar = 0.0003	(mho/cm2)
	Aalfa = 11.7336 ( /ms)
	Valfa = 50 ( /mV) 
	Abeta = 0.0324  ( /ms)
	Vbeta = 5.5 ( /mV)
	Agamma = 150 ( /ms)
	Adelta = 40  ( /ms)
	ACon = 0.001    ( /ms)
	ACoff = 10     ( /ms)
	AOon = 0.2     ( /ms)
	AOoff = 0.001   ( /ms)
	n1 = 32.532
	n2 = 0.123
	
	
}

ASSIGNED {
	ica  (mA/cm2)
	g   (mho/cm2)
	
	gamma
	delta
	Con
	Coff
	Oon
	Ooff
	a
	b
	Q10
	
}

STATE {
	C1
	C2
	C3
	O
	I1
	I2
	I3
	I4
}


INITIAL {
	C1=1
	C2=0
	C3=0
	O=0
	I1=0
	I2=0
	I3=0
	I4=0
	Q10 =3^((celsius-32.76(degC))/10 (degC))
	gamma = Q10 * Agamma
	delta = Q10 * Adelta
	Con = Q10 * ACon
	Coff = Q10 * ACoff
	Oon = Q10 * AOon
	Ooff = Q10 * AOoff
	a = (Oon/Con)^0.5
	b = (Ooff/Coff)^0.5
	
    }
    
	
BREAKPOINT {	    
   
    SOLVE kstates METHOD sparse
    g = gcabar * O	      	: (mho/cm2)
    ica = g * (v - eca)  	: (mA/cm2)
}


FUNCTION alfa(v(mV))(/ms){ 
	alfa = Q10*Aalfa*exp((v-Vshift)/Valfa) 
}

FUNCTION beta(v(mV))(/ms){ 
	beta = Q10*Abeta*exp((-v+Vshift)/Vbeta) 
}



KINETIC kstates {
	: 1 riga
	~ C1 <-> C2 (n1*alfa(v),n2*beta(v))
	~ C2 <-> C3 (n2*alfa(v),n1*beta(v))
	~ C3 <-> O  (gamma,delta)
	
	: 2 riga
	~ I1 <-> I2	(n1*alfa(v)*a,n2*beta(v)*b)
	~ I2 <-> I3	(n2*alfa(v)*a,n1*beta(v)*b)
	~ I3 <-> I4 (gamma,delta)
	
	: connette 1 riga con 2 riga
	~ C1 <-> I1 (Con,Coff)
	~ C2 <-> I2 (Con*a,Coff*b)
	~ C3 <-> I3 (Con*a^2,Coff*b^2)
	~  O <-> I4 (Oon,Ooff)
	
	
	CONSERVE C1+C2+C3+O+I1+I2+I3+I4=1
}

