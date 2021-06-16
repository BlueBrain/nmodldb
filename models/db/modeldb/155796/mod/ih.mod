TITLE I-h channel
: Modelled from Rugiero et al 2002, J Physiol 538:2, p447
: Written by Jordan Chambers (jordandchambers@gmail.com)

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 		(mV)
        eih = -38	(mV)        
	celsius 	(degC)
	gbar = 1e-4 	(mho/cm2)
        vhalf = -72   	(mV)
        vslope1 = 8.2  	(mV)
        vslope2 = 11.9 	(mV)
        tmc1 = 537     	(ms)
        tmc2 = 56     	(ms)
	q10 = 4.5	:guessed from other Ih studies at similar temperatures
	jtmc = 1 (1)
	jiih (mA/cm2)
}


NEURON {
	SUFFIX ih
	NONSPECIFIC_CURRENT i
	GLOBAL gbar, vhalf, vslope1, vslope2, tmc1, tmc2, linf, taul
	RANGE jiih
	THREADSAFE linf, taul
}

STATE {
        l
}

ASSIGNED {
	i (mA/cm2)
        linf      
        taul
        ghd
}

INITIAL {
	rate(v)
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ghd = gbar*l
	i = ghd*(v-eih)
	jiih = i
}

DERIVATIVE states {
        rate(v)
        l' =  (linf - l)/taul
}

PROCEDURE rate(v (mV)) {
        LOCAL a,qt
        qt=q10^((celsius-34)/10)

	linf = 1/(1 + exp((v - vhalf)/vslope1))
	taul = jtmc*(tmc1 + (qt*tmc2/(exp((v-vhalf)/vslope2) + exp(-(v-vhalf)/vslope2))))
}

