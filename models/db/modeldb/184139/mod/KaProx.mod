TITLE K-A channel from Klee Ficker and Heinemann
: modified to account for Dax A Current ----------
: M.Migliore Jun 1997
: The original model results in too big window current, which is robably contamination of D-type current.
: To reduce the window current component, LSH changed vhalfn, zeatan and zetal as follows
: 	1) vhalfn, 0 -> -15 mV
:	2) vhalfl, -60 -> -70 mV
: 	3) zetan, -1.5 -> -3
: 	4) zetal, 3 -> 5 
: To fit these parameters to Kim Jonas (2012), they should be -5, -65, -1.8, 3.7, respectively.


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	celsius = 24	(degC)
	gbar = .008 (mho/cm2)
    vhalfn = -15    (mV)
    vhalfl = -70  (mV)
    a0l = 0.05    (/ms)
    a0n = 0.05    (/ms)
    zetan = -3    (1)
    zetal = 5  (1)
    gmn = 0.55   (1)
    gml = 1   (1)
	lmin = 6.7 (mS)
	nmin = 0.1 (mS)
	pw = -1    (1)
	tq = -40
	qq = 5
	q10 = 5
	qtl = 0.5
    FRT = 39 (coulombs/joule) 	
}


NEURON {
	SUFFIX KaProx
	USEION k WRITE ik
    RANGE gbar,gka,ik
    GLOBAL ninf,linf,taul,taun,lmin
}

STATE {
	n
    l
}

ASSIGNED {
	ik (mA/cm2)
    ninf
    linf      
    taul
    taun
    gka
}

INITIAL {
	rates(v)
	n=ninf
	l=linf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gka = gbar*n*l
	ik = gka*(v+90.0)
}


FUNCTION alpn(v(mV)) {
LOCAL zeta
  zeta = zetan + pw/(1+exp((v-tq)/qq))
  alpn = exp(1.e-3*zeta*(v-vhalfn)*FRT) 
}

FUNCTION betn(v(mV)) {
LOCAL zeta
  zeta = zetan + pw/(1+exp((v-tq)/qq))
  betn = exp(1.e-3*zeta*gmn*(v-vhalfn)*FRT) 
}

FUNCTION alpl(v(mV)) {
  alpl = exp(1.e-3*zetal*(v-vhalfl)*FRT) 
}

FUNCTION betl(v(mV)) {
  betl = exp(1.e-3*zetal*gml*(v-vhalfl)*FRT) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v)
        n' = (ninf - n)/taun
        l' = (linf - l)/taul
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-24)/10)
        a = alpn(v)
        ninf = 1/(a+1)
        taun = betn(v)/(qt*a0n*(1+a))
		if (taun<nmin) {taun=nmin}
        a = alpl(v)
        linf = 1/(a+1)
		taul = 0.26*(v+50)/qtl
		: taul = 0.26*v + 13.2, Hoffman Nat 98
		if (taul<lmin/qtl) {taul=lmin/qtl}
}
