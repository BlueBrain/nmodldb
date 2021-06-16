TITLE INa for horizontal cell
: Voltage gated fast sodium current (INa) for horizontal cells
: 
: Based on parameters of Aoyama et al. (2000)


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX HzINa
    USEION na READ ena WRITE ina
    RANGE gbar
    RANGE m_inf, h_inf
    RANGE tau_m, tau_h
    RANGE m_exp, h_exp
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    gbar = 2.2642e-05 (mho/cm2) : 2.4 nS total
}

STATE {
    m h
}

ASSIGNED {
    v       (mV)
    ena     (mV)
    ina     (mA/cm2)
    celsius (degC)
    dt      (ms)

    m_inf
    h_inf
    tau_m
    tau_h
    m_exp
    h_exp
    tadj
}

BREAKPOINT {
    SOLVE states
    ina = gbar * m*m*m*h * (v - ena)
}

PROCEDURE states() {
    : exact when v held constant
    evaluate_fct(v)
    m = m + m_exp * (m_inf - m)
    h = h + h_exp * (h_inf - h)
}

UNITSOFF

INITIAL {
    m = 0.026
    h = 0.922
    tadj = 3.0 ^ ((celsius-25)/10)  : correction for physio temp
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a, b, v2
    a = (200 * (38-v)) / ( exp( (38-v)/25 ) - 1)
    b = 2000 * exp( -(55+v) / 18 )
    tau_m = 1 / (a + b) / tadj
    m_inf = a / (a + b)

    a = 1000 * exp( -(80+v) / 8 )
    b = 800 / ( exp((80-v)/75) + 1 )
    tau_h = 1 / (a + b) / tadj
    h_inf = a / (a + b)

    m_exp = 1 - exp(-dt/tau_m)
    h_exp = 1 - exp(-dt/tau_h)
}

UNITSON

