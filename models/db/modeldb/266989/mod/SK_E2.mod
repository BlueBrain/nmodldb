: SK-type calcium-activated potassium current
: Reference : Kohler et al. 1996

: 2019: From ModelDB, accession no. 139653

NEURON {
       SUFFIX SK_E2
       USEION k READ ek WRITE ik
	USEION cal1 READ cal1i VALENCE 2
       RANGE gSK_E2bar, gSK_E2, ik, zTau
}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      (mM) = (milli/liter)
}

PARAMETER {
          v            (mV)
          gSK_E2bar = .000001 (mho/cm2)
          zTau = 1              (ms)
          ek           (mV)
          cal1i          (mM)
}

ASSIGNED {
         zInf
         ik            (mA/cm2)
         gSK_E2	       (S/cm2)
}

STATE {
      z   FROM 0 TO 1
}

BREAKPOINT {
           SOLVE states METHOD cnexp
           gSK_E2  = gSK_E2bar * z
           ik   =  gSK_E2 * (v - ek)
}

DERIVATIVE states {
        rates(cal1i)
        z' = (zInf - z) / zTau
}

PROCEDURE rates(ca(mM)) {
          if(ca < 1e-7){
	              ca = ca + 1e-07
          }
          zInf = 1/(1 + (0.00043 / ca)^4.8)
}

INITIAL {
        rates(cal1i)
        z = zInf
}
