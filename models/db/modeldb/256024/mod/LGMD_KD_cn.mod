TITLE cAMP inactivated K-D channel from RBD


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _KD

        SUFFIX KD_cn
        USEION k READ ek WRITE ik
        USEION cn READ cni VALENCE 1
        RANGE gmax, g, taun, taul
        GLOBAL tnmax, tlmax
}

PARAMETER {
	gmax=0.01 (mho/cm2)
	vhalfn=-52	(mV)
	vn2=-65		(mV)
	zn=6.0		(mV)
	tnmax=20	(ms)
	tnmin=1.5	(ms)
	tns=-7.5	(mV)
	np=2
	
	vhalfl=-67	(mV)
	zl=-3.6		(mV)
	tlmax=2000	(ms)
	tlmin=50	(ms)
	vl2=-75		(mV)
	tls=14		(mV)
	
	cntm=50		(1)
	lcp=2		(1)
	kD=3e-4		(mM)
}

STATE {
	n
	l
}

ASSIGNED {
    v (mV)
    ek (mV)
    cni (mM)

	ik (mA/cm2)
	ninf
	linf
	tf	(1)
	taul (ms)
	taun (ms)
	g (S/cm2)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gmax*n*l
	ik = g*(v-ek)
}

INITIAL {
	rates(v)
	n = ninf
	l = linf
	tf = 1
}

FUNCTION alpn(v(mV)) {
  alpn = exp((vhalfn-v)/zn)
}

FUNCTION betn(v(mV)) {
  betn = exp((vn2-v)/tns) 
}

FUNCTION alpl(v(mV)) {
  alpl = exp((vhalfl-v)/zl)
}

DERIVATIVE states {  
	rates(v)
	lci(cni)
	n' = (ninf - n)/taun
	l' = (linf - l)/(taul/tf)
}

PROCEDURE rates(v (mV)) { :callable from hoc
	LOCAL a
	TABLE ninf, taun, linf, taul DEPEND vhalfn, tlmax, tnmax, tnmin
          FROM -100 TO 50 WITH 600
    
	a = alpn(v)
	ninf = 1/(1 + a)
	taun = 4*(tnmax-tnmin)/(1+betn(v))*ninf+tnmin
	a = alpl(v)
	linf = (1/(1+ a))
	taul = 2*tlmax/( exp((v-vl2)/tls) + exp((vl2-v)/tls) ) + tlmin
}


PROCEDURE lci(cni (mM)) { :callable from hoc
	TABLE tf DEPEND lcp, kD
          FROM 0 TO 0.001 WITH 500
    
	tf = 1+cntm-cntm/(1+(cni/kD)^lcp)
    
}
