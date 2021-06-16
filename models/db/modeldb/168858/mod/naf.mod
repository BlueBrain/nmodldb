TITLE Transient inactivating Sodium current (naf), J Neurophysiol 89:909-921, 2003

COMMENT
       Implemented by Aniruddha Yadav 2007 (aniruddha.yadav@mssm.edu)

ENDCOMMENT

UNITS {
        (mV) = (millivolt)
}

NEURON {
         SUFFIX naf
         RANGE gbar
         USEION na READ ena WRITE ina
         RANGE Vna, ina, vrev
         RANGE a1, a2, b1, b2, c1, c2
}

PARAMETER {
           gbar=0.0    (mho/cm2)
           Vna = 10.0  (mV)
           a1 =0.025    (ms)
           b1 = 0.14   (ms)
           a2 = 0.02   (ms)
           b2=  0.145   (ms)
           c1= 0.15    (ms)
           c2=0.15     (ms)
           ena        (mV)
           vrev = -3.5 (mV)
}

ASSIGNED {
           ina     (mA/cm2)
           minf   (1)
           mtau   (ms)
           v      (mV)
           hinf   (1)
           htau   (ms)
}

STATE {
     m h
}

INITIAL {
         rates(v)
         m=minf
         h=hinf
}

BREAKPOINT {
             SOLVE states METHOD cnexp
             ina = gbar * m * m * m * h * (v - ena )
}

            
DERIVATIVE states {
        rates(v)
        m' = (minf - m ) / mtau
        h' = (hinf - h ) / htau
}

UNITSOFF

PROCEDURE rates(V (mV)) {
        
         minf  = 1 / ( 1 + exp( ( - ( V + vrev) - 38 ) / 10 ) )
        if( ( V + vrev ) < -30.0 ) {
                mtau = a1 + b1 * exp( ( ( V + vrev ) + 30 ) / 10 )
        } else{
                mtau = a2 + b2 * exp( ( - ( V + vrev ) - 30 ) / 10 )
        }

        : hinf, and htau are shifted 3.5 mV comparing to the paper

        hinf  = 1 / ( 1 + exp( ( ( V + vrev  ) + 62.9 ) / 10.7 ) )
        htau = c1 + c2 / ( 1 + exp( ( ( V + vrev  ) + 37 ) / 15 ) )
}


UNITSON
