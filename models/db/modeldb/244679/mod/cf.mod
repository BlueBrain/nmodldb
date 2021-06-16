: synaptic current of climbing fiber
COMMENT

The model of the current associated with a CF-EPSP, designed to mimic the shape of the current
reported by Llano et al., 1991.
The parameters  were tuned to obtain the match of experimental V m and Ca 2+ transients.

Current Model Reference: Karima Ait Ouares , Luiza Filipis , Alexandra Tzilivaki , Panayiota Poirazi , Marco Canepari (2018) Two distinct sets of Ca 2+ and K + channels 
are activated at different membrane potential by the climbing fibre synaptic potential in Purkinje neuron dendrites. 

PubMed link: 

Contact: Filipis Luiza (luiza.filipis@univ-grenoble-alpes.fr)
ENDCOMMENT
NEURON {
    SUFFIX climbing
    RANGE del,tauO,tauC, icin, delf,tauOf,tauCf,fastfact
    NONSPECIFIC_CURRENT  i
    RANGE i, e, ef, g
}

PARAMETER {
    g = 15e-3 (siemens/cm2)  < 0, 1e9 >
    e = 0    (millivolt)
    del = 1 (ms)
    tauO = 0.4 (ms)
    tauC = 14 (ms)  
        delf = 1 (ms)
    tauOf = 0.4 (ms)
    tauCf = 14 (ms)
    fastfact=0 
    icin=0   (nanoamp/cm2)
    ef = 0    (millivolt)

}
UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (nA) = (nanoamp)
    (pA) = (picoamp)
    (S)  = (siemens)
}
ASSIGNED {
    i   (nanoamp/cm2)
    v   (millivolt) 

}



BREAKPOINT {
    at_time(del)
    if (t < del) {
          i = icin
    } else {
        
          i =icin+(1-fastfact)*g*(1-exp(-(t-del)/tauO))*exp(-(t-del)/tauC)*(e)+fastfact*g*(1-exp(-(t-delf)/tauOf))*exp(-(t-delf)/tauCf)*(ef)
          
           
    }
}