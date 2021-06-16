TITLE kml.mod  

COMMENT

Created by Christian Roessert

Implementation of a simple high-threshold potassium current (Prescott et al., 2008) using the
Morris-Lecar formalism (Morris and Lecar, 1981)

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX kml
        USEION k READ ek WRITE ik
        RANGE gbar, g, bn, gn, tn, ik
        GLOBAL ninf, ntau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        ek = -90 (mV)
        gbar = 0 (S/cm2) <0,1e9>
        bn = -20 (mV)
        gn = 10 (mV) 
        tn = 3 (ms)
}

STATE {
        n
}

ASSIGNED {
    ik (mA/cm2) 
    g (S/cm2)
    ninf
    ntau (ms)
    }

BREAKPOINT {
	SOLVE states METHOD cnexp 

	g = gbar*n
    ik = g*(v - ek)

}

INITIAL {
    rates(v)
    n = ninf
}


DERIVATIVE states {  
        rates(v)
        n' =  (ninf-n)/ntau
}


UNITSOFF

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
    
    : ninf = (1/2) * (1 + tanh( (v-bn) / gn ) )
    ninf = (1 / (1 + exp((bn - v) / (gn/2))))
    : ntau = tn / cosh( (v-bn) / (2*gn) )
    ntau = (2*tn) / ( exp((v-bn)/(2*gn)) + exp((-v+bn)/(2*gn)) )
}

UNITSON

: v=[-80:1:0]; bn=-20; gn=10; tn=2;
: ninf = (1/2) * (1 + tanh( (v-bn) ./ gn ) )
: ninf2 = (1 ./ (1 + exp((bn - v) ./ (gn/2))))
: ntau = tn ./ cosh( (v-bn) ./ (2*gn) )
: ntau2 = (2*tn) ./ ( exp((v-bn)./(2*gn)) + exp((-v+bn)./(2*gn)) )
: plot(v, ninf, v, ninf2, 'r--')
: plot(v, ntau, v, ntau2, 'r--')
