TITLE sodium calcium exchange for bladder small DRG neuron soma model

: Adapted from Forrest et al., 2014
: and Courtemanche et al.,1998 Am J Physiol  275:H301
: revised from Luo and Rudy 

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

NEURON {
	SUFFIX ncxsoma
	USEION ca READ cao, cai WRITE ica
	USEION na READ nao, nai WRITE ina
	RANGE ImaxNax, ica, ina
    RANGE KnNacx, KcNacx, itotalncx
}

UNITS {
	(mA)	= (milliamp)
    (mV)	= (millivolt)
	(molar) = (1/liter)
	(mM)	= (millimolar)
	F		= (faraday) (coulombs)
	R		= (k-mole)	(joule/degC)
}

PARAMETER {
	ImaxNax  =  1.1e-5	 (mA/cm2) <0,1e6>
	KnNacx   =  87.5     (mM)   <0,1e6>
	KcNacx   =  1.38     (mM)   <0,1e6>
}

ASSIGNED {
	ica		  (mA/cm2)
	ina		  (mA/cm2)
	itotalncx (mA/cm2)
        
    celsius	  (degC)
	v			(mV)
    cao 		(mM)
    cai 		(mM)
	nao 		(mM)
	nai 		(mM)
}
 
BREAKPOINT { LOCAL Kqa, KB, k
    k	= (1e3)*R*(celsius + 273.14)/(F)
	Kqa = exp(0.35*v/k)
	KB	= exp( - 0.65*v/k)
			
	itotalncx = ImaxNax*(Kqa*nai^3*cao-KB*nao^3*cai)/((KnNacx^3 + nao^3)*(KcNacx + cao)*(1 + 0.1*KB))
	ina =  3*itotalncx
	ica = -2*itotalncx
}






