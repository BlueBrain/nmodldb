TITLE Potasium dr type current for RD Traub, J Neurophysiol 89:909-921, 2003

COMMENT
       Implemented by Aniruddha Yadav 2007 (aniruddha.yadav@mssm.edu)

ENDCOMMENT

UNITS {
        (mV) = (millivolt)
}

NEURON {
         SUFFIX kdr
         RANGE gbar
         USEION k READ ek WRITE ik
         RANGE Vkd, ik, vrev
         RANGE a1, a2, b1, b2
}

PARAMETER {
           gbar=0.0    (mho/cm2)
           Vkd = 10.0  (mV)
           a1 =0.25    (ms)
           b1 = 4.35   (ms)
           a2 = 0.25   (ms)
           b2=  4.35   (ms)
           ek        (mV)
           vrev = 29.5 (mV)
}

ASSIGNED {
           ik     (mA/cm2)
           minf   (1)
           mtau   (ms)
           v      (mV)
}

STATE {
     m
}

INITIAL {
         rates(v)
         m=minf
}

BREAKPOINT {
             SOLVE states METHOD cnexp
             ik = gbar * m * m * m * m * (v - ek )
}

INITIAL {
          rates(v)
          m = minf
}
            
DERIVATIVE states {
        rates(v)
        m' = (minf - m ) / mtau
}

UNITSOFF

PROCEDURE rates(V (mV)) {

        minf  = 1 / ( 1 + exp( ( -V - vrev) / 10 ) )
        if( V < -Vkd ) {
                mtau = a1 + b1 * exp( ( V + Vkd ) / 10 )
        }else{
                mtau = a2 + b2 * exp( ( -V - Vkd ) / 10 )
        }
}

UNITSON
