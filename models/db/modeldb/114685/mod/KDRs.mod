TITLE potassium delayed rectifier Kv2.1 membrane channel for GPi neuron model

COMMENT
 Potassium Kv2.1 membrane channel for GPi exhibiting *slow* deactivating
 components.  Based on derived kinetics from Baranauskas 1999.  Their
 primary experiments were performed at 35degC.

 Baranauskas 1999 -- Say in paper: Vh=-18mV, Vc=7mV
    Based on my own fitting of their data: Vh=-21.98mV, Vc=8.245mV 
 Hernandez 1999 -- tau at -40mV is 28.80ms
    Scaled the Gillies2006 curve to match this point

 Disregard Q10 given the proximity of the experiments to 36degC

ENDCOMMENT


UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX KDRs
    USEION k READ ki,ek WRITE ik
    RANGE gk
}

PARAMETER {
    v             (mV)
    dt            (ms)
    gk = 0.0030   (mho/cm2)  : Baranauskas 1999
    ek
    ki
    celsius
}

STATE {
    pslowd
    pslowi
}

ASSIGNED { 
    ik (mA/cm2)
    pinf_sd
    ptau_sd (ms)
    pinf_si
    ptau_si (ms)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = (gk)*pslowd*pslowi*(v-ek)
}

UNITSOFF

INITIAL {
    settables(v)
    pslowd = pinf_sd
    pslowi = pinf_si
}

DERIVATIVE states {  
    settables(v)
    pslowd' = (pinf_sd-pslowd)/ptau_sd
    pslowi' = (pinf_si-pslowi)/ptau_si
}

PROCEDURE settables(v) {
	
    LOCAL rate_k
    TABLE pinf_sd, ptau_sd, pinf_si, ptau_si DEPEND celsius FROM -100 TO 100 WITH 400

    rate_k =  4^((36-23)/10)
    pinf_sd = 1/(1+exp(-(v+21.98)/8.245))                        : slow deactivating component (parameters from Baranauskus1999)
    :pinf_sd = 1/(1+exp((17.5-v)/10))                              : slow deactivating component (from Mohapatra/Held modelDB)
    ptau_sd = rate_k*5/(exp((v+30)/14) + exp(-(v+30)/20)) + 5     : slow deactivating component (from Mohapatra/Held modelDB)

    pinf_si = 0.74/(1 + exp((v+25)/12)) + 0.26                    : slow inactivating component (from Mohapatra/Held modelDB)
    ptau_si = rate_k*500/(exp((v-50)/20) + exp((50-v)/20))+800    : slow inactivating component (from Mohapatra/Held modelDB)

}

UNITSON

