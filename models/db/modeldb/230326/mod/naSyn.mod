TITLE naSyn.mod

COMMENT

Syanpse with Na+ current

ENDCOMMENT

NEURON {
	POINT_PROCESS naSyn
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
	RANGE onset, tau_onset, tau_offset, coeff_onset, coeff_offset, gmax, erev, ina, ik, q, gna_max, gk_max
}


UNITS { 
    (nA) = (nanoamp) 
    (mV) = (millivolt) 
    (uS) = (microsiemens)
}


PARAMETER { 
    onset=0 (ms) 
    tau_onset=2.837         (ms) 
    tau_offset=121.286     (ms)
    
    coeff_onset=2.26
    coeff_offset=1
    
    gmax=0.0005    (uS)
    erev=10.        (mV)

    v (mV) 
}

ASSIGNED {
    
    q       (1)
    
    ena     (mV)
    ek      (mV)

    gna_max (uS)
    gk_max  (uS)

    gna     (uS)
    gk      (uS)
  
    ina     (nA)
    ik      (nA)

}

INITIAL {

gk_max=gmax/(1+((ek-erev)/(erev-ena)))
gna_max=gmax-gk_max


}

BREAKPOINT { 

    q=syn_func(t - onset, tau_onset, tau_offset, coeff_onset, coeff_offset)
    
    gna = gna_max * q
    gk = gk_max * q
 
    ina = gna*(v - ena)
    ik = gk*(v - ek)
 
}

FUNCTION syn_func(x, tau_onset, tau_offset, coeff_onset, coeff_offset) { 
    LOCAL vv
    if (x < 0 || x > 500) { 
        vv = 0
    } 
    else { 
        vv = coeff_offset*exp(-1*x/tau_offset)-coeff_onset*exp(-1*x/tau_onset)
    }
    if (vv > 0) {
        syn_func = vv
    }
    else {
        syn_func = 0
    }
    
     
}

