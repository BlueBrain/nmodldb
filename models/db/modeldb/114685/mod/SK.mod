TITLE small conductance calcium activated potassium channels for GPi neuron model

COMMENT

 Small-conductance Ca2+ activated K+ current (SK) with SK2 subunits are found
 to a high degree in rat EP neurons (Stocker/Pedarzani2000).  Kinetics were 
 based on SK2 channel activity from Hirschberg (1998), which was recorded at 
 room temperature (22-24degC).

 Used the steady-state function from Gillies2006, but changed the time constant
 to reflect the averaged open-time distributions; note the taus did not depend
 on voltage.

 Q10=1.5 --> rate_k=exp(log(Q10)*((1/296)-(1/309))/((1/292)-(1/302)))=1.66

ENDCOMMENT

NEURON {
    SUFFIX sKCa
    USEION ca READ cai
    USEION k READ ki,ek WRITE ik
    RANGE  gk,isKCa
    GLOBAL sKCatau,rate_k,gmax_k
}

UNITS {
    (mM) = (milli/liter)
    (mA) = (milliamp)
    F = (faraday) (coulombs)	: Faradays constant 
}

PARAMETER {
    v (mV)
    dt (ms)
    gk = 0.0001 (mho/cm2)
    isKCa = 0.0 (mA/cm2)
    sKCatau = 6.1 (ms)
    ek 
    ki
    cai
    celsius	
}

ASSIGNED {
    ica (mA/cm2)
    ik (mA/cm2)
    winf 
    wtau (ms)
    rate_k
    gmax_k
}

STATE {
    w
}

BREAKPOINT {
    SOLVE integrate METHOD cnexp
    ik = (gk*gmax_k)*w*(v-ek)
    isKCa = ik
}

UNITSOFF

INITIAL {
    rate_k = 1.66
    gmax_k = 1.66
    setinf(cai)
    w = winf
}

DERIVATIVE integrate {
    setinf(cai)
    w' = (winf - w)/wtau
}

PROCEDURE setinf(cai) {
    LOCAL wcai
    : these equations are for uM calcium concentrations
    wcai = cai*1000
    winf = 0.81/(1+exp((llog(wcai)+0.3)/-0.46))
    wtau = sKCatau/rate_k
}

FUNCTION llog(x) {  :returns log of x, but error checks first
    if (x>1e-11) {
        llog = log(x)
    }else{
        llog=0
    }
}

UNITSON
