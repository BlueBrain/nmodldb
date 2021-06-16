COMMENT
current noise injection, changes between + and - imax every dt
Paul Bush 1995
Note this file apparently has to be modified for UNIX implementation
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
:     for unix  POINT_PROCESS NOISE
        SUFFIX NOISE
        NONSPECIFIC_CURRENT i
        RANGE imax
:       GLOBAL seed
}

ASSIGNED {
: for unix      rn
}

UNITS {
        (nA) = (nanoamp)
}

PARAMETER {
        imax=0          (umho)
:                               seed=1
}

INITIAL {
:                               srandom(seed)
: for unix      rn = (1/(2^31))*2
}

ASSIGNED { i (nA) }

BREAKPOINT {

SOLVE dum       :has to be in a proc otherwise error

}

PROCEDURE dum() {

        i = (scop_random()-0.5)*imax

        : for unix      i = ( random()*rn-1 )*imax
        :VERBATIM
        :       printf("%f      %f\n",i,t);
        :ENDVERBATIM
        VERBATIM
                return 0;
        ENDVERBATIM

}
