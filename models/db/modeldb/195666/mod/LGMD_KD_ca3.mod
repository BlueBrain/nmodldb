TITLE Ca2+ deinactivated K-D channel from RBD


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

        SUFFIX KD_ca3
        USEION k READ ek WRITE ik
        USEION ca READ cai
        RANGE gmax, g, taun, taul
        GLOBAL tnmax, tlmax
}

PARAMETER {
	: all values can be adjusted in hoc files
	gmax=0.01 	(mho/cm2)
	vhalfn=-49	(mV)
	vn2=-60		(mV)
	zn=7.0		(mV)
	tnmax=10	(ms)
	tnmin=1.5	(ms)
	tns=-5.5	(mV)
	np=2
	
	vhalfl=-70	(mV)
	zl=-3.5		(mV)
	tlmax=1050	(ms)
	tlmin=70	(ms)
	vl2=-75		(mV)
	tls=15		(mV)
	
	cavm=18		(mV)
	lcp=1		(1)
	kD=4e-4		(mM)
	tauca=20	(ms)
}

STATE {
	n
	l
	ov	(mV)
	ovs	(mV)
}

ASSIGNED {
    v (mV)
    ek (mV)
    cai (mM)

	ik (mA/cm2)
	ninf
	linf
	vs	(mV)
	taul (ms)
	taun (ms)
	g (S/cm2)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	l = 1/(1+exp((vhalfl-ov+ovs)/zl))
	g = gmax*n^np*l
	ik = g*(v-ek)
}

INITIAL {
	lf(v)
	rates(v)
	n = ninf
	l = 1/(1+exp((vhalfl-v)/zl)*exp(vs/zl))
	ov = v
	ovs = vs
}

FUNCTION alpn(v(mV)) {
  alpn = exp((vhalfn-v)/zn)
}

FUNCTION betn(v(mV)) {
  betn = exp((vn2-v)/tns) 
}

DERIVATIVE states {
	lci(cai)
	rates(v)
	n' = (ninf - n)/taun
	:l' = (linf - l)/taul
	:l' = (1/zl)*exp((vhalfl-)/zl)*exp(vs/zl)/(1+exp((vhalfl-v)/zl)*exp(vs/zl))^2
	:l' = linf - l
	ov' = (v-ov)/taul
	ovs' = (vs-ovs)/tauca
}

PROCEDURE lf(v(mV)) {
	LOCAL dvdt, dsdt
	lci(cai)
	:linf = 1/(1+exp((vhalfl-v)/zl)*exp(vs/zl))
	linf = 1/(1+exp((vhalfl-ov+ovs)/zl))
}

PROCEDURE rates(v (mV)) { :callable from hoc
	LOCAL a
	TABLE ninf, taun, taul DEPEND vhalfn, tlmax, tnmax, tnmin
          FROM -100 TO 50 WITH 600
    
	a = alpn(v)
	ninf = 1/(1 + a)
	taun = 4*(tnmax-tnmin)/(1+betn(v))*ninf+tnmin
	taul = 2*tlmax/( exp((v-vl2)/tls) + exp((vl2-v)/tls) ) + tlmin
}


PROCEDURE lci(cai (mM)) { :callable from hoc
	TABLE vs DEPEND lcp, kD, cavm
          FROM 0 TO 0.02 WITH 500
    
	vs = cavm-cavm/(1+(cai/kD)^lcp)
    
}
