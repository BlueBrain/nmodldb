
NEURON {
  SUFFIX nas_rs
  USEION na READ ena WRITE ina
  RANGE gna, g
}

UNITS {
  (S) = (siemens)
  (mV) = (millivolt)
  (mA) = (milliamp)
}

PARAMETER {
  gna = 0.1125 (S/cm2) :0.1125
  theta_m = -22.8 (mV): -22.8 control -22.8
  sigma_m = 11.8 (mV) :11.8
  theta_h = -62.9 (mV) : 
  sigma_h = -10 (mV) :-10.7
  theta_t_h = -60 (mV)
  sigma_t_h = -12 (mV)
  :taum = 0.1 (ms) : for stability with dt>0.01 ms 0.001
  q10 =3
}

ASSIGNED {
  v (mV)
  ena (mV)
  ina (mA/cm2)
  g (S/cm2)
  qt
}

STATE {
  m
  h
}

BREAKPOINT {
  SOLVE states METHOD cnexp
  g = gna * h * m^3
  ina = g * (v-ena)
}

INITIAL {
  m = minfi(v)
  h = hinfi(v)
}

DERIVATIVE states {
  m' = (minfi(v)-m)/taum(v)
  h' = (hinfi(v)-h)/tauh(v)
}

FUNCTION hinfi(v (mV)) {
  UNITSOFF
  hinfi=1/(1 + exp(-(v-theta_h)/sigma_h))
  UNITSON
}

FUNCTION tauh(v (mV)) (ms) {
  UNITSOFF
  qt = q10^((celsius - 24)/10)

  tauh = (.31 + 14 / ( 1 + exp(-(v-theta_t_h)/sigma_t_h))) :0.5 0.71
  tauh = tauh/qt
  
  UNITSON
}

FUNCTION taum(v (mV)) (ms) {
  UNITSOFF
  qt = q10^((celsius - 24)/10)

  taum = (((0.022 + 3.6 / (1 + exp ((v+27.9)/7.6))) * (0.009 + 1.9 / (1 + exp (-(v-1.3)/12.7)))))
  taum = taum/qt
  :printf("%f\n", qt)
  UNITSON
}

FUNCTION minfi(v (mV)) {
  UNITSOFF
  minfi=1/(1 + exp(-(v-theta_m)/sigma_m))
  UNITSON
}
