TITLE sjg_na.mod  The low threshold conductance of cochlear nucleus neurons

COMMENT

NEURON implementation of Jason Rothman's measurements of VCN conductances.

Original implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.

File split implementaiton, April 1, 2004.

Contact: pmanis@med.unc.edu

fit to MNTB data by Sarah J. Griffin, MRC Toxicology Unit, Leicester

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX sjg_na
        USEION na WRITE ina
        RANGE gnabar, gna, ina, sjgena
        GLOBAL hinf, minf, htau, mtau, pinf, ptau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        dt (ms)
        sjgena (mV)
        gnabar =  0.07958 (mho/cm2) <0,1e9>
    	fi = 0.8448 <0,1>	: fast inactivating
    	si = 0.0352 <0,1>	: slow inactivation
}

STATE {
        m h p
}

ASSIGNED {
    ina (mA/cm2) 
    gna (mho/cm2)
    minf hinf pinf
    mtau (ms) htau (ms) ptau (ms)
    }

LOCAL mexp, hexp, pexp

BREAKPOINT {
	SOLVE states
    
    gna = gnabar*((m^3* 0.005) + (m^3*h*0.97) + (m^3*p*0.025))
    ina = gna*(v - sjgena)

}

UNITSOFF

INITIAL {
    trates(v)
    m = minf
    h = hinf
    p = pinf
}

PROCEDURE states() {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	m = m + mexp*(minf-m)
	h = h + hexp*(hinf-h)
	p = p + pexp*(pinf-p)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
		      
		      : Fixed to match jamie's data
    minf = 1 / (1+exp(-(v + 38) / 7))
    hinf = 1/(1+exp((v + 55.4) / 6.3))
    pinf = 1/(1+exp((v + 55.4) / 6.3))


    mtau =  (5 / (5*exp((v+60) / 18) + 36*exp(-(v+60) / 25))) + 0.01
    htau =  (50 / (7*exp((v+60) / 11) + 10*exp(-(v+60) / 25))) + 0.5
    ptau =  (50 / (7*exp((v+60) / 11) + 10*exp(-(v+60) / 25))) + 50

}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	TABLE minf, mexp, hinf, hexp:, pinf, pexp
	DEPEND dt, celsius FROM -150 TO 150 WITH 300

    rates(v)    : not consistently executed from here if usetable_hh == 1
        : so don't expect the tau values to be tracking along with
        : the inf values in hoc

	mexp = 1 - exp(-dt/mtau)
	hexp = 1 - exp(-dt/htau)
	pexp = 1 - exp(-dt/ptau)
	}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON
