TITLE potassium delayed rectifier Kv3.1/3.2 membrane channel for GPi model neuron

COMMENT

 Potassium Kv3.1/3.2 membrane channel exhibiting *fast* deactivating
 components.  Based on derived kinetics from Baranauskas 1999 and Hernandez1999,
 in which the experiments were performed at 35degC and 22degC, respectively.

 Baranauskas 1999 -- Say in paper: Vh=-13mV, Vc=6mV
    Based on my own fitting of their data: Vh=-16.2mV, Vc=8.636mV
 Hernandez 1999 -- tau at -40mV is 2.27ms; fit tau to figure 11D
    Q10 = 1.7 --> rate_k = exp(log(Q10)*((1/295)-(1/309))/((1/292)-(1/302))) = 2.05

ENDCOMMENT

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX KDRf
    USEION k READ ki,ek WRITE ik
    RANGE gk
    GLOBAL rate_k,gmax_k
}

PARAMETER {
    v             (mV)
    dt            (ms)
    gk = 0.0038   (mho/cm2)  : Baranauskas 1999
    ek
    ki
    celsius
}

STATE {
    pfast 
}

ASSIGNED { 
    ik (mA/cm2)
    pinf
    ptau (ms)
    rate_k
    gmax_k
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = (gk*gmax_k)*pfast*(v-ek)
}

UNITSOFF

INITIAL {
    rate_k = 2.05
    gmax_k = 2.05
    settables(v)
    pfast = pinf
}

DERIVATIVE states {  
    settables(v)
    pfast' = (pinf-pfast)/ptau
}

PROCEDURE settables(v) {
	
    TABLE pinf, ptau DEPEND celsius FROM -100 TO 100 WITH 400
    
    pinf = 1/(1+exp(-(v+16.2)/8.6))                              : fast component (parameters from Baranauskas1999)
    ptau = 6.7/(exp(-(v+21.7)/21.2)+exp((v-11.7)/21.2))/rate_k   : fast component (fit to Hernandez1999-Fig11D)

}

UNITSON