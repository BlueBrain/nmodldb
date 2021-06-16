TITLE NMDA
 
UNITS {
        (pA) = (picoamp)
        (molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX nmda
	USEION ca WRITE ica
	USEION na  READ nai WRITE ina
	USEION k  WRITE ik
	RANGE  ica,ina,ik, inmda,Pbar,mg,km,nai,base,pr,q
        GLOBAL pinf
 
}

UNITS {
	:FARADAY = 96520 (coul)
	:R = 8.3134 (joule/degC)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
        v (mV)
	celsius= 35  	(degC)
	Pbar=2.0e-6	(cm/s)	: Maximum Permeability
	cai = 70e-6	(mM)
	cao = 2		(mM)
	acao = 0.3 
        anai = 0.75
        anao = 0.75
        pr = 0.0225
	nao = 145	(mM)
	ki =  140	(mM)
	ko = 2.5	(mM)
        aki = 0.75
        ako = 0.75
        dt (ms)
        gnmdabar = .0011 (mho/cm2) 
        q=9 (mV)
        km=50.7 (mM)
        mg = 1.2 (mM)
        enmda = 0 (mV)
        base = -0.001 (mA/cm2)
 
}
STATE {
         p <1e-8>
}


ASSIGNED { 
	   nai          (mM)
           ica		(mA/cm2)
           ina		(mA/cm2)
           ik		(mA/cm2)
           inmda 	(mA/cm2)
            pinf
}


BREAKPOINT {
     SOLVE states METHOD cnexp
	ina = Pbar*p*gold(v, anai*nai, anao*nao,1)
	ica = 10.6*Pbar*p*gold(v, cai, acao*cao,2)
	ik = Pbar*p*gold(v, aki*ki, ako*ko,1)
        inmda = ica + ik + ina
}


FUNCTION gold(v(mV), ci(mM), co(mM),z) (0.001 coul/cm3) {
	LOCAL arg, eci, eco
	arg = (0.001)*z*FARADAY*v/(R*(celsius+273.15))
	eco = co*efun(arg)
	eci = ci*efun(-arg)
	gold = (0.001)*z*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

UNITSOFF
DERIVATIVE states {
	pinf = pr + (1.0 - pr)/(1 + (mg/km)*exp(-v/q))
	p' = pinf - p 
	}
UNITSON
