TITLE GABA_A synapse with short-term plasticity

NEURON {
    POINT_PROCESS tmGabaA
    RANGE tau1, tau2, e, i, q
    RANGE tau, tauR, tauF, U, u0
    RANGE base, f_gaba
    POINTER pka
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
}

PARAMETER {
    tau1= 0.5 (ms)
    tau2 = 7.5 (ms)  : tau2 > tau1
    e = -60 (mV)
    tau = 3 (ms)
    tauR = 500 (ms)  : tauR > tau
    tauF = 0 (ms)    : tauF >= 0
    U = 0.1 (1) <0, 1>
    u0 = 0 (1) <0, 1>
    q = 2
    base   = 0.0      : set in simulation file    
	f_gaba = 0.0      : set in simulation file
}

ASSIGNED {
    v (mV)
    i (nA)
    g (uS)
    factor
    x
    pka (1)
}

STATE {
    A (uS)
    B (uS)
}

INITIAL {
    LOCAL tp
    A = 0
    B = 0
    tp = (tau1*tau2)/(tau2-tau1) * log(tau2/tau1)
    factor = -exp(-tp/tau1) + exp(-tp/tau2)
    factor = 1/factor
    tau1 = tau1/q
    tau2 = tau2/q
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    g = B - A
    i = modulation(f_gaba)*g*(v - e)
}

DERIVATIVE state {
    A' = -A/tau1
    B' = -B/tau2
}

NET_RECEIVE(weight (uS), y, z, u, tsyn (ms)) {
    INITIAL {
        y = 0
        z = 0
        u = u0
        tsyn = t
    }
    z = z*exp(-(t-tsyn)/tauR)
    z = z + (y*(exp(-(t-tsyn)/tau) - exp(-(t-tsyn)/tauR)) / (tau/tauR - 1) )
    y = y*exp(-(t-tsyn)/tau)
    x = 1-y-z
    if (tauF > 0) {
        u = u*exp(-(t-tsyn)/tauF)
        u = u + U*(1-u)
    } else {
        u = U
    }
    A = A + weight*factor*x*u
    B = B + weight*factor*x*u
    y = y + x*u
    tsyn = t
}

FUNCTION modulation(modFact) {
    
    : returns modulation factor
    
    modulation = 1 + modFact * (pka - base)
    
}

COMMENT

Implementation of GABA_A synapse model with short-term facilitation
and depression based on modified tmgsyn.mod [1] by Tsodyks et al [2].
Choice of time constants follows [3].  NEURON implementation by Alexander
Kozlov <akozlov@kth.se>.

[1] tmgsyn.mod, ModelDB (https://senselab.med.yale.edu/ModelDB/),
accession number 3815.

[2] Tsodyks M, Uziel A, Markram H (2000) Synchrony generation in recurrent
networks with frequency-dependent synapses. J Neurosci. 20(1):RC50.

[3] Wolf JA, Moyer JT, Lazarewicz MT, Contreras D, Benoit-Marand M,
O'Donnell P, Finkel LH (2005) NMDA/AMPA ratio impacts state transitions
and entrainment to oscillations in a computational model of the nucleus
accumbens medium spiny projection neuron. J Neurosci 25(40):9080-95.
ENDCOMMENT
