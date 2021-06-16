TITLE K-D channel from RBD

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S) = (siemens)
        (mS) = (millisiemens)
}

NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will be having the suffix _KD

        SUFFIX KD
        USEION k READ ek WRITE ik
        RANGE gmax, g, taun, taul
        GLOBAL tnmax, tlmax
}

PARAMETER {
	gmax=0.01 (mho/cm2)
	vhalfn=-45	(mV)
	vn2=-60		(mV)
	zn=7.8		(mV)
	tnmax=25	(ms)
	tnmin=1.0	(ms)
	tns=-6.0	(mV)
	np=1
	
	vhalfl=-68	(mV)
	zl=-5.1		(mV)
	tlmax=2500	(ms)
	tlmin=800	(ms)
	vl2=-67		(mV)
	tls=14		(mV)
}

STATE {
        n
        l
}

ASSIGNED {
    v (mV)
    ek (mV)

	ik		(mA/cm2)
	ninf
	linf
	taul	(ms)
	taun	(ms)
	g		(S/cm2)
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
	TABLE ninf, taun, linf, taul DEPEND vhalfn, vhalfl, tlmax, tnmax, zn, zl
          FROM -125 TO 25 WITH 750
    
	a = alpn(v)
	ninf = 1/(1 + a)
	:taun = tnmax/(1+betn(v))*ninf+tnmin
	taun = 4*(tnmax-tnmin)/(1+betn(v))*ninf+tnmin
	a = alpl(v)
	linf = 1/(1+ a)
	:taul = 4*(tlmax-tlmin)/(1+exp((vhalfl-v)/tls))*ninf + tlmin
	taul = 2*tlmax/( exp((v-vl2)/tls) + exp((vl2-v)/tls) ) + tlmin
}

