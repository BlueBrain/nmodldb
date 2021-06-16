NEURON {
    SUFFIX dIN_ca
    USEION ca  WRITE ica
    RANGE perm, ica, T, Sin, Sout
    RANGE alpha_A, alpha_B, alpha_C, alpha_D, alpha_E
    RANGE beta_A1, beta_B1, beta_C1, beta_D1, beta_E1
    RANGE beta_A2, beta_B2, beta_C2, beta_D2, beta_E2
}

UNITS {
    (nA) = (nanoamp)
    (uA) = (microamp)
    (mA) = (milliamp)
    (A) = (amp)
    (mV) = (millivolt)
    (mS) = (millisiemens)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (kHz) = (kilohertz)
    (mM) = (millimolar)
    (um) = (micrometer)
    (S) = (siemens)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/kelvin)
}

PARAMETER {
    T = 300.0 (kelvin): temperature
    z = 2 : ionic valence
    Sin = 100e-9 (moles/cm2)
    Sout = 10e-6 (moles/cm2)
    perm = 1.425e-10 : not sure about units here
    
    alpha_A = 4.05
    alpha_B = 0.0
    alpha_C = 1.0
    alpha_D = -15.32
    alpha_E = -13.57
    
    beta_A1 =  0.98859:1.24
    beta_B1 = 0.093
    beta_C1 = -1.0
    beta_D1 = 10.63
    beta_E1 = 1.0

    beta_A2 = 1.28
    beta_B2 = 0.0
    beta_C2 = 1.0
    beta_D2 = 5.39
    beta_E2 = 12.11
}

ASSIGNED {
    v (mV)
    ica (mA/cm2)
    cainf
    catau (ms)
}

STATE {
    ca
}

INITIAL {
    rates()
    ca = cainf
}


BREAKPOINT {
    SOLVE states METHOD cnexp
    ica = pow(ca,2)*ghk()
}

DERIVATIVE states {
    rates()
    ca' = (cainf - ca) / catau
}

UNITSOFF
PROCEDURE rates() {
    LOCAL alpha, beta
    if (v < -25.0) {
      :alpha = alphabeta(4.05,0.0,1.0,-15.32,-13.57)
      :beta = alphabeta(0.98859,0.093,-1.0,10.63,1.0)
      alpha = (alpha_A+alpha_B*v)/(alpha_C+exp((alpha_D+v)/alpha_E))
      beta = (beta_A1+beta_B1*v)/(beta_C1+exp((beta_D1+v)/beta_E1))
    }else {
      :alpha = alphabeta(4.05,0.0,1.0,-15.32,-13.57)
      :beta = alphabeta(1.28,0.0,1.0,5.39,12.11)
      alpha = (alpha_A+alpha_B*v)/(alpha_C+exp((alpha_D+v)/alpha_E))
      beta = (beta_A2+beta_B2*v)/(beta_C2+exp((beta_D2+v)/beta_E2))
    }
     catau = 1.0/(alpha + beta)
     cainf = alpha * catau
}

FUNCTION alphabeta(A,B,C,D,E){
    alphabeta = (A + B*v)/(C + exp((v+D)/E))
}

: Goldman-Hodgkin-Katz equation
FUNCTION ghk(){
	LOCAL cv, e
	: Careful here, didn't check the units, correcting by (1e3) multiplication (to mA)
	cv = (1e-3)*v*(z*FARADAY)/(R*T)
	e = exp(-cv)
	ghk = (1e3)*perm*cv*z*FARADAY*(Sin - Sout*e)/(1-e)
}
UNITSON


