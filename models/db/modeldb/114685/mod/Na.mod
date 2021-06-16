TITLE sodium membrane channels for GPi model neuron

COMMENT

 Sodium from a CA1/3 pyramidal neuron, Traub et al 1991.  
 Based on Sah (1988) data, which were at 22-24degC, and Gillies2006
 The Q10 measured from Sah (for the peak conductance)
     temp = [17.5,20,22.5,26.5], log(INa) = [0.65,0.92,1.03,1.22], INa = [1.9155,2.5093,2.8011,3.3872]
     slope = 0.1581 nA/degC
     Q10 = gmaxQ10 = 1.581 nA/10degC
     rate_k = exp(ln(Q10)*((1/296)-(1/309))/((1/292)-(1/302)))=1.78

ENDCOMMENT

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX Na
    USEION na READ nai,ena WRITE ina
    RANGE gna, m, h
    GLOBAL rest,rate_k
}

PARAMETER {
    v   (mV)
    dt  (ms)
    gna   = 1e-7    (mho/cm2)
    rest  = -60.0   (mV)
    ena             (mV)
    nai
    celsius	
}

STATE {
    m h  
}

ASSIGNED { 
    ina     (mA/cm2)
    alpham  (/ms)
    betam   (/ms)
    alphah  (/ms)
    betah   (/ms)
    rate_k
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina  = (gna)*m*m*h*(v-ena)
}

UNITSOFF

INITIAL {
    settables(v)
    m = alpham/(alpham+betam)
    h = alphah/(alphah+betah)
    rate_k = 1.78               : based on calculated Q10 measurement
}

DERIVATIVE states {
    settables(v)      :Computes state variables at the current v and dt.
    m' = alpham * (1-m) - betam * m
    h' = alphah * (1-h) - betah * h
}

PROCEDURE settables(v) {

    LOCAL vadj
    TABLE alpham, betam, alphah, betah DEPEND rest,celsius FROM -100 TO 100 WITH 400

    vadj = v - rest

    :"m" sodium activation system
    alpham = 0.32*vtrap((13.1-vadj),4)
    betam =  0.28*vtrap((vadj-40.1),5)  : NOTE used Traub1991 and not Gillies2006

    :"h" sodium inactivation system
    alphah = 0.128*exp((17-vadj)/18)
    betah = 4/(exp((40-vadj)/5)+1)

}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate equations
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{
        vtrap = x/(exp(x/y) - 1)
    }
}
 
UNITSON