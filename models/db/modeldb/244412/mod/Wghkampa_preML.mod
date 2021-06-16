COMMENT

The kinetics part is obtained from Exp2Syn of NEURON.

Two state kinetic scheme synapse described by rise time taur, and
decay time constant taud. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/taur and 1/taud is

A = a*exp(-t/taur) and
G = a*taud/(taud-taur)*(-exp(-t/taur) + exp(-t/taud))
where taur < taud

If taud-taur -> 0 then we have a alphasynapse.
and if taur -> 0 then we have just single exponential decay.

The factor is evaluated in the initial block such that an event of
weight 1 generates a peak conductance of 1.

Because the solution is a sum of exponentials, the coupled equations
can be solved as a pair of independent equations by the more efficient
cnexp method.

Added by Rishikesh Narayanan:

1. GHK based ionic currents for AMPA current 
2. Weights, and their update, according Shouval et al., PNAS, 2002.

Details may be found in:

Narayanan R, Johnston D. The h current is a candidate mechanism for 
regulating the sliding modification threshold in a BCM-like synaptic 
learning rule.  J Neurophysiol. 2010 Aug;104(2):1020-33.

ENDCOMMENT

NEURON {
    POINT_PROCESS Wghkampa_preML
    USEION na WRITE ina
    USEION k WRITE ik
    USEION ca READ cai	: Weight update requires cai 
    USEION glut READ gluti WRITE iglut,gluti VALENCE 0
    
    RANGE taur, taud
    RANGE iampa,winit, iglut
    RANGE P, Pmax, lr
    
    :Presynaptic
    RANGE g_factor, I, U_SE_factor, glut_factor
    RANGE U_SE, U_SE_init, tau_takeover, tau_in
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (mM) = (millimolar)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

PARAMETER {
    taur=2 		(ms) <1e-9,1e9>
    taud = 10 	(ms) <1e-9,1e9>
    nai = 18	(mM)	: Set for a reversal pot of +55mV
    nao = 140	(mM)
    ki = 140	(mM)	: Set for a reversal pot of -90mV
    ko = 5		(mM)
    cai			(mM)
    celsius		(degC)
    Pmax=1e-6   (cm/s)	
    alpha1=0.35	:Parameters for the Omega function.
    beta1=80
    alpha2=0.55
    beta2=80
    : winit=1		(1)
    
    :Presynaptic
    tau_rec = 0.8 (ms)
    tau_in = 3 (ms)
    U_SE_init = 0.1 (1)
    U_SE_factor = 0 (1)
    tau_takeover = 0.1 (ms)
    
    : Glutamate
    glut_factor = 40 (1) : conversion from I to glutamate
    g_factor = 1 (1)
}


ASSIGNED {
    ina     (nA)
    ik      (nA)
    v (mV)
    P (cm/s)
    factor
    iampa	(nA)
    lr
    Area (cm2)
    U_SE (1)
    :Presinaptici
    I (1)
    
    : Glutamate
    iglut (nA)
}


STATE {
    A (cm/s)
    B (cm/s)
    : w (1)
    x  (1)
    y_rel  (1)
    z  (1)    
    gluti (mM)
}

INITIAL {
    LOCAL tp
    if (taur/taud > .9999) {
	taur = .9999*taud
    }
    A = 0
    B = 0
    tp = (taur*taud)/(taud - taur) * log(taud/taur)
    factor = -exp(-tp/taur) + exp(-tp/taud)
    factor = 1/factor
    Area=1
    U_SE = U_SE_init
    : w=winit
    
    :Presynaptic
    x = 1 (1)
    y_rel = 0 (1)
    z = 0 (1)
    I = 0 (1)

    gluti = 0 (mM)
}


BREAKPOINT {
    SOLVE state METHOD cnexp
    P=B-A
    
    : Area is just for unit conversion of ghk output
    
    : ina = P*w*ghk(v, nai, nao,1)*Area * g_factor	
    : ik = P*w*ghk(v, ki, ko,1)*Area * g_factor
    ina = P*ghk(v, nai, nao,1)*Area * g_factor	
    ik = P*ghk(v, ki, ko,1)*Area * g_factor
    iampa = ik + ina
    : printf("bp%g\t",gluti)
}

DERIVATIVE state {
    lr=eta(cai)
    : w' = lr*(Omega(cai)-w)
    A' = -A/taur
    B' = -B/taud
    
    U_SE = U_SE_init * (1 + U_SE_factor)
    x' = z / tau_rec - U_SE * x * I
    y_rel' = -y_rel / tau_in + U_SE * x * I
    z = 1 - x - y_rel
    gluti' = -gluti / tau_takeover + y_rel * glut_factor :U_SE * x * I
    :printf("x = %g,\t gluti = %g,\t I = %g\n",x,gluti, I)
}

FUNCTION ghk(v(mV), ci(mM), co(mM),z) (0.001 coul/cm3) {
    LOCAL arg, eci, eco
    arg = (0.001)*z*FARADAY*v/(R*(celsius+273.15))
    eco = co*efun(arg)
    eci = ci*efun(-arg)
    ghk = (0.001)*z*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
    if (fabs(z) < 1e-4) {
	efun = 1 - z/2
    }else{
	efun = z/(exp(z) - 1)
    }
}

FUNCTION eta(ci (mM)) { : when ci is 0, inv has to be 3 hours.
    LOCAL inv, P1, P2, P3, P4
    P1=100	
    P2=P1*1e-4	: There was a slip in the paper, which says P2=P1/1e-4
    P4=1e3
    P3=3		: Cube, directly multiplying, see below.
    
    ci=(ci-1e-4)*1e3 	: The function takes uM, and we get mM.
    
    inv=P4 + P1/(P2+ci*ci*ci) :As P3 is 3, set ci^P3 as ci*ci*ci.
    eta=1/inv
}	

FUNCTION Omega(ci (mM)) {
    ci=(ci-1e-4)*1e3	: The function takes uM, and we get mM.
    Omega=0.25+1/(1+exp(-(ci-alpha2)*beta2))-0.25/(1+exp(-(ci-alpha1)*beta1))
}

NET_RECEIVE(weight (uS)) { 	
    : Presynaptic terminal STP
    if( flag==0 ) {  
	I = 10
	net_send(0.1,2)
	}	
    if( flag==2 ) { 
	I = 0
	: gluti = 0
    }

    : No use to weight, can be used instead of Pmax,
    : if you want NetCon access to the synaptic
    : conductance.
    : printf("%g\t",y)
    A = A + Pmax * factor * gluti
    B = B + Pmax * factor * gluti
    
}   

    
