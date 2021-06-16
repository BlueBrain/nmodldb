COMMENT
%W%                    %G%
ENDCOMMENT

COMMENT
The non specific current eliminates the electrical effect of the calcium current. To obtain the 
electrogenic Na/Ca exchanger, the non specific current should be multiplied by 3/2.
ENDCOMMENT
INDEPENDENT {v FROM -100 TO 50 WITH 50 (mV)}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
        (mV) =  (millivolt)
	(um) = (micron)
	(mol) = (1)
	PI = (pi) (1)
	FARADAY = (faraday) (coulomb)
}

NEURON {
	SUFFIX nacaexch
	USEION ca READ cao, cai WRITE ica
	NONSPECIFIC_CURRENT i
        RANGE i,iexch
	GLOBAL k, Nai, Nao, E1, E2, c1, c2, e1, e2,s
}


PARAMETER {
	cao    		(mM)
	cai		(mM)

	k  = 9.355e-7   (mA/mM4-cm2)
	Nai = 9.5710389      (mM)
	Nao = 152.43715     (mM)
	e1 = 0.3516 
	e2 = 0.6825
	E1
	E2
	c1
	c2
	s = 1
}

INITIAL {
	c1 = Nai*Nai*Nai
	c2 = Nao*Nao*Nao
	E1 = e1/26.73
	E2 = e2/26.73
}

ASSIGNED {
        iexch   (mA/cm2)
	ica	(mA/cm2)
	i	(mA/cm2)
}

UNITSOFF


BREAKPOINT {
	iexch = -k*(cao*c1*exp(E1*(v)) - cai*c2*exp(-E2*(v)))
	ica = iexch
	i = -s*iexch
}

UNITSON

