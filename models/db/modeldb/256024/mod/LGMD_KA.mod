TITLE K-A channel from RBD

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S) = (siemens)
}

NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _KA

	SUFFIX KA
	USEION k READ ek WRITE ik
	RANGE gmax,g
	GLOBAL tnmax,tlmax
}

PARAMETER {
	gmax=0.01 (mho/cm2)
	vhalfn=-40	(mV)
	vn2=-50		(mV)
	tnmax=3.5	(ms)
	tnmin=0.5	(ms)
	tns=-8		(mV)
	zn=7.5		(mV)
	np=1		(1)

	vhalfl=-45	(mV)
	vl2=-50		(mV)
	tlmax=35	(ms)
	tlmin=3.0	(ms)
	tls=12		(mV)
	zl=-5		(mV)
}

STATE {
        n
        l
}

ASSIGNED {
    v (mV)
    ek (mV)
	
	ik (mA/cm2)
	ninf (1)
	linf (1)  
	taul (ms)
	taun (ms)
    i (mA/cm2)
	g (S/cm2)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gmax*n^np*l
	ik = g*(v-ek)

}

INITIAL {
	rates(v)
	n=ninf
	l=linf
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
	l' = (linf - l)/taul
}

PROCEDURE rates(v (mV)) { :callable from hoc
	LOCAL a
	TABLE ninf, taun, linf, taul DEPEND vhalfn, vhalfl, tlmax, tnmax, tnmin, tlmin, zn, zl
		FROM -100 TO 50 WITH 600
        
		a = alpn(v)
		ninf = 1/(1 + a)
		taun = 4*(tnmax-tnmin)/(1+betn(v))*ninf+tnmin
		a = alpl(v)
		linf = 1/(1+ a)
		taul = 4*(tlmax-tlmin)/(1+exp((vl2-v)/tls))*linf + tlmin
}

