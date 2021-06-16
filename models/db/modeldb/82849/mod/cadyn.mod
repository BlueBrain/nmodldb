: simple first-order model of calcium dynamics

NEURON {
        SUFFIX cadyn
        USEION ca READ cai,ica WRITE cai 
        RANGE CAF, tca, cai
}

UNITS {
        (mM) = (milli/liter)
        (mA) = (milliamp)
        F    = (faraday) (coul)
}

PARAMETER {
        tca= 70 (ms)           : decay time constant
        cainf= 50e-6   (mM)      : (50e-6)equilibrium ca2+ concentration
        dep= 2e-4 (micron)     : depth of shell for ca2+ diffusion

}

ASSIGNED {
        ica     (mA/cm2)
        diam    (micron)
        A       (/coul/cm)
        CAF     ()
}

STATE { cai (mM) }

BREAKPOINT { 
        SOLVE states METHOD cnexp
}

DERIVATIVE states {
         cai'= -A*CAF*ica - (cai-cainf)/tca
}

INITIAL {
        A =(1e4)/(F*dep)
}




