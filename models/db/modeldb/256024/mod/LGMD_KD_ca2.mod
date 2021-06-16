TITLE K-D channel from RBD

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(S) = (siemens)
}

NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _KD

        SUFFIX KD_ca2
        USEION ca READ cai
        USEION k READ ek WRITE ik
        RANGE gmax, g, taun, taul, kD_ca, l
        GLOBAL vhalfn, tnmax,tlmax, tnmin, tauc, cap
}

PARAMETER {
	gmax=0.01 (mho/cm2)
	vhalfn=-45	(mV)
	vn2=-60		(mV)
	zn=7.0		(mV)
	tnmax=100	(ms)
	tnmin=3.0	(ms)
	tns=-7.5	(mV)
	
	vhalfl=-66	(mV)
	zl=-2.0		(mV)
	tlmax=820	(ms)
	tlmin=20	(ms)
	tauc=30		(ms)
	vl2=-63		(mV)
	tls=17		(mV)
	kD_ca = 0.0001 (mM)
	cap = 0.6
	lcp=6
}

STATE {
        n
        lv
        lc
}

ASSIGNED {
    v (mV)
    ek (mV)
    cai (mM)
    
	ik (mA/cm2)
	l
	ninf
	linf
	lcinf
	taul (ms)
	taun (ms)
	g (S/cm2)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	l= (1-(1-lc)*(1-lv))*cap+(1-cap)*lv
	g = gmax*n*l
	ik = g*(v-ek)
}

INITIAL {
	rates(v)
	n=ninf
	lv = linf
	lci(cai)
	lc = lcinf
	l= (1-(1-lc)*(1-lv))*cap+(1-cap)*lv
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
	n' = (ninf - n)/taun
	lv' = (linf - lv)/taul
	lci(cai)
	lc' = (lcinf - lc)/tauc
}

PROCEDURE rates(v (mV)) { :callable from hoc
	LOCAL a
	TABLE ninf, taun, linf, taul DEPEND vhalfn, tlmax, tnmax, tnmin
          FROM -100 TO 50 WITH 600
    
	a = alpn(v)
	ninf = 1/(1 + a)
	:taun = tnmax/(1+betn(v))*ninf+tnmin
	taun = 4*(tnmax-tnmin)/(1+betn(v))*ninf+tnmin
	a = alpl(v)
	linf = (1/(1+ a))
	:taul = 4*(tlmax-tlmin)/(1+exp((vhalfl-v)/tls))*ninf + tlmin
	taul = 2*tlmax/( exp((v-vl2)/tls) + exp((vl2-v)/tls) ) + tlmin
}

PROCEDURE lci(cai (mM)) { :callable from hoc
	TABLE lcinf DEPEND lcp, kD_ca
          FROM 0 TO 0.01 WITH 1000
    
    lcinf = (1-(cai/(cai+kD_ca))^(lcp/2))^lcp
    
}
