: simple first-order model of calcium dynamics
: from Durstewitz & Gabriel (2006), Cerebral Cortex

NEURON {
        SUFFIX cadyn
        USEION ca READ cai, ica WRITE cai 
        RANGE CAF, tc, cai
}

UNITS {
        (mM) = (milli/liter)
        (mA) = (milliamp)
        F    = (faraday) (coul)
}

PARAMETER {
        tc= 70 (ms)           : decay time constant
        cainf= 50e-6   (mM)      : (50e-6)equilibrium ca2+ concentration
        dep= 2e-4 (micron)     : depth of shell for ca2+ diffusion
}

ASSIGNED {
        ica     (mA/cm2)
        diam    (micron)
        A       (/coul/cm)
        CAF     ()
}

INITIAL {
        A =(1e4)/(F*dep)
	cai=cainf
}

STATE { cai (mM) }

BREAKPOINT { 
        SOLVE states METHOD derivimplicit
}

DERIVATIVE states {
         cai'= -A*CAF*ica - (cai-cainf)/tc
}
