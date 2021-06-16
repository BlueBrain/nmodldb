TITLE mitralSyn.mod

COMMENT

VNO sensory neuron to AOB mitral cell synapse

ENDCOMMENT

NEURON {
	POINT_PROCESS mitralSyn
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
	RANGE onset, tau_onset, tau_offset1, tau_offset2, coeff_onset, coeff_offset1, coeff_offset2, gmax, erev, ina, ik
}


UNITS { 
    (nA) = (nanoamp) 
    (mV) = (millivolt) 
    (uS) = (microsiemens)
}


PARAMETER { 
    onset=0 (ms) 
    tau_onset=2.837         (ms) 
    tau_offset1=121.286     (ms)
    tau_offset2=1077.702    (ms)
    
    coeff_onset=2.26
    coeff_offset1=1
    coeff_offset2=0.265
    
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

    q=syn_func(t - onset, tau_onset, tau_offset1, tau_offset2, coeff_onset, coeff_offset1, coeff_offset2)
    
    gna = gna_max * q
    gk = gk_max * q
 
    ina = gna*(v - ena)
    ik = gk*(v - ek)
 
}

FUNCTION syn_func(x, tau_onset, tau_offset1, tau_offset2, coeff_onset, coeff_offset1, coeff_offset2) { 
    if (x < 0 || x > 10000) { 
        syn_func = 0
    } 
    else { 
        syn_func = coeff_offset1*exp(-1*x/tau_offset1)+coeff_offset2*exp(-1*x/tau_offset2)-coeff_onset*exp(-1*x/tau_onset)
    } 
}

