COMMENT
Two state kinetic scheme synapse described by rise time taur,
and decay time constant taud. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/taur and 1/taud is
A = a*exp(-t/taur) and
G = a*taud/(taud-taur)*(-exp(-t/taur) + exp(-t/taud))
where taur < taud

If taud-taur -> 0 then we have a alphasynapse.
and if taur -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
    POINT_PROCESS ghknmda
    USEION na WRITE ina
    USEION k WRITE ik
    USEION ca READ cai, cao WRITE ica
    USEION glut READ gluti VALENCE 0
    
    RANGE taur, taud
    RANGE inmda
    
    RANGE P, mg, Pmax
    RANGE  mgb, ica, Area, mgb_k, mg_ref
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
    taur=5 (ms) <1e-9,1e9>
    taud = 50 (ms) <1e-9,1e9>
    cai = 100e-6(mM)	: 100nM
    cao = 2		(mM)
    nai = 18	(mM)	: Set for a reversal pot of +55mV
    nao = 140	(mM)
    ki = 140	(mM)	: Set for a reversal pot of -90mV
    ko = 5		(mM)
    celsius		(degC)
    mg = 1		(mM) : 2 mM in the Johnston et al. 2010, extracellula [MgCl2] = 1 mM in Edelman et al. 2015
    Pmax=1e-6   (cm/s)	: According to Canavier, PNMDAs default value is
    : 1e-6 for 10uM, 1.4e-6 cm/s for 30uM of NMDA
    Area = 1 (cm2)
    mgb_k = 0.062 (/mV)
    mg_ref = 3.57 (mM)
}

ASSIGNED {
    ina     (nA)
    ik      (nA)
    ica     (nA)
    v (mV)
    P (cm/s)
    factor
    mgb	(1)
    inmda	(nA)
    
    
    gluti (mM)
}

STATE {
    A (cm/s)
    B (cm/s)
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
    : Area=1
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    P=B-A
    mgb = mgblock(v)
    
    : Area is just for unit conversion of ghk output
    
    ina = P*mgb*ghk(v, nai, nao,1)*Area	
    ica = P*10.6*mgb*ghk(v, cai, cao,2)*Area
    ik = P*mgb*ghk(v, ki, ko,1)*Area
    inmda = ica + ik + ina
    : printf("nmda%g\t",gluti)
    : printf("nmdaP %g\t",P)
}

DERIVATIVE state {
    A' = -A/taur
    B' = -B/taud
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

FUNCTION mgblock(v(mV)) (1){
    TABLE 
    DEPEND mg
    FROM -140 TO 80 WITH 1000 
    
    : from Jahr & Stevens, JNS, 1990
    
    mgblock = 1 / (1 + exp(mgb_k * -v) * (mg / mg_ref)) 
    : remove the background activation at -70 mV
    : if (mgblock < 0.036 ) {
    : 	mgblock = 0
    : }
}

NET_RECEIVE(weight (uS)) { 	: No use to weight, can be used instead of Pmax,
    : if you want NetCon access to the synaptic
    : conductance.
    : printf("nmda_sp%g\t",gluti)
    if (flag == 0 ) {
	net_send(.01,2)
    }
    if (flag == 2 ) {
	A = A + Pmax*factor * gluti
	B = B + Pmax*factor * gluti
    }
}
