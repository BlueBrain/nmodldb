COMMENT

THIS IS THE MORE EFFICIENT VERSION of EPSPPlas.
Since short-term plasticity, and the synaptic dynamics is the same
for all the presynaptic terminals of a neuron, we can calculate these
parameters presynaptically, and then pass them to each synapse model.

For computational efficiency it is obviously better to process the release
parameters in the presynaptic mechanism, rather than do the same calculations
in each synaptic mechanism.  The problem is the synaptic delay.
To deal with this problem I have created a history (histR, histG) of the
synaptic conductances, that are accessed by the synaptic mechanisms. The vector
index accessed corresponds to the delay.
In the vector the first index (0) is always the current time step.
Thus if dt = 0.1 to implement a delay of 1ms you should setpointer:
setpointer syn.R_1, IN[0].soma.histR_ExIAF[10]
for a zero ms delay
setpointer syn.R_1, IN[0].soma.histR_ExIAF[0]



ENDCOMMENT


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
   SUFFIX EtoIPlasSom
        RANGE C, lastrelease, lastspike, releaseat, Delay
        GLOBAL Cdur, Deadtime, terror
        GLOBAL Alpha_1, Beta_1
        RANGE  ampa, R0_1, R1_1, Rinf_1, Rtau_1
        GLOBAL Alpha_2, Beta_2
   GLOBAL DurNMDA, tauNMDA
   :SHORT-TERM PLASTICITY
   GLOBAL U, trec, tfac
   RANGE  R, u, RG

}

UNITS {
        (nA) = (nanoamp)
        (mV) = (millivolt)
        (umho) = (micromho)
        (mM) = (milli/liter)
}

STATE {
        nmda                            : fraction of open NMDA channels
        R_2
}

PARAMETER {

   terror
        Cdur    = 1       (ms)          : transmitter duration (rising phase)
        Deadtime = 1  (ms)              : mimimum time between release events
        Delay = 1     (ms)

        Alpha_1 = 1.5   (/ms mM)        : AMPA forward (binding) rate
        Beta_1 = 0.75   (/ms)           : AMPA backward (unbinding) rate
        Alpha_2 = 0.25 (/ms mM)         : NMDA forward (binding) rate
        Beta_2 = 0.025 (/ms)            : NMDA backward (unbinding) rate
        DurNMDA = 0.4           : used in sigmoid R_G, small values long time-peak
        tauNMDA      =50            : in ms
        tfac = 500     (ms)             : this value should be close to 0 for no facilitation
        trec = 125              (ms)       : recovery from depression time constant
        U = 0.2                                            : percent of transmitter released on first pulse
}


ASSIGNED {

        dt               (ms)
        v                (mV)           : postsynaptic voltage
        i                (nA)           : current = g*(v - Erev)
        C                (mM)           : transmitter concentration

        ampa
        R0_1                            : open channels at start of release
        R1_1                            : open channels at end of release
        Rinf_1                  : steady state channels open
        Rtau_1  (ms)            : time constant of channel binding
        RG                              : reflects binding of Transm -> G protein

        lastrelease     (ms)            : time of last spike
        lastspike       (ms)
        releaseat

        R                               : Releasable pool
        u                               : for running value of U


}

INITIAL {
   terror = dt/10
        C = 0
        ampa = 0
        lastrelease = -9e4
        lastspike   = -9e4
        releaseat   = -9e4

        R = 1
        u = U

}

BREAKPOINT {
    SOLVE release
}

PROCEDURE release() { LOCAL q
    :will crash if user hasn't set pre with the connect statement
        :FIND OUT THERE WAS A SPIKE
        q = (t - lastspike)             : time since last release ended

        if (q > Deadtime) {             : ready for another release?
                if (v > 0) {            : spike occured?
                        lastspike = t
                        releaseat = t + Delay
                }
        }

        : CALCULATE RELEASE PARAMETERS
        q = (t - lastrelease -Cdur)                     : time since last spike with delay
        if (q > Deadtime) {          : start release

                if (t > releaseat - terror && t < releaseat + terror) {
                        lastrelease = t
                u = u*(exptable(-q/tfac)) + U*(1-u*exptable(-q/tfac))
                        R = R*(1 - u)*exptable(-q/trec) + 1 - exptable(-q/trec)
                        C = R*u                         : start new release, turn on
                        Rinf_1 = C*Alpha_1 / (C*Alpha_1 + Beta_1)
                        Rtau_1 = 1 / ((Alpha_1 * C) + Beta_1)
                        R0_1 = ampa
                        R_2=(1-R_2)*0.5+R_2
                }

        } else if (q < 0) {                     : still releasing?
                : do nothing
        } else if (C > 0) {                     : in dead time after release, turn off
                C = 0.
                R1_1 = ampa
        }
        if (C > 0) {                            : transmitter being released?
           ampa = Rinf_1 + (R0_1 - Rinf_1) * exptable (- (t - lastrelease) / Rtau_1)
        } else {                                        : no release occuring
           ampa = R1_1 * exptable (- Beta_1 * (t - (lastrelease + Cdur)))
        }

    SOLVE G_protein METHOD cnexp

        VERBATIM
        return 0;
        ENDVERBATIM
}

DERIVATIVE G_protein {                                                  : ready for anotherelease?

    R_2'=-(R_2/tauNMDA)

        if (R_2<0.01) {
                RG = 0.0
        } else {
                RG = 1/( 1 + exptable(-((R_2)-DurNMDA)/0.05) )          : binding of T -> G
        }
        nmda' = Alpha_2 * RG * (1-nmda) - Beta_2 * nmda

    :printf("D-------->%f\n",t)

}

FUNCTION exptable(x) {
        TABLE  FROM -10 TO 10 WITH 2000

        if ((x > -10) && (x < 10)) {
                exptable = exp(x)
        } else {
                exptable = 0.
        }
}

