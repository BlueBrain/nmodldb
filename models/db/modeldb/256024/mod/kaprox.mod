TITLE K-A channel from Klee Ficker and Heinemann
: modified to account for Dax A Current --- M.Migliore Jun 1997
: modified to be used with cvode  M.Migliore 2001

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	gkabar=.008 (mho/cm2)
	vhalfn=11	(mV)
	vhalfl=-56	(mV)
	a0l=0.05	(/ms)
	a0n=0.05	(/ms)
	zetan=-1.5	(degC/mV)
	zetal=3		(degC/mV)
	gmn=0.55	(1)
	gml=1		(1)
	lmin=2		(mV)
	nmin=0.1	(ms)
	pw=-1		(degC/mV)
	tq=-40		(mV)
	qq=5		(mV)
	q10=5		(1)
	qtl=1		(mV/ms)
}


NEURON {
    THREADSAFE

	SUFFIX kap
	USEION k READ ek WRITE ik
	RANGE gkabar,gka
	GLOBAL lmin
}

STATE {
	n
	l
}

ASSIGNED {
	v (mV)
	celsius	(degC)
	ek		(mV)
	ik		(mA/cm2)
	ninf	(1)
	linf	(1)
	taul	(ms)
	taun	(ms)
	gka		(mho/cm2)
}

INITIAL {
	rates(v)
	n=ninf
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	gka = gkabar*n*l
	ik = gka*(v-ek)

}


FUNCTION alpn(v(mV)) {
LOCAL zeta
  zeta=zetan+pw/(1+exp((v-tq)/qq))
  alpn = exp(1.e-3*zeta*(v-vhalfn)*9.648e4/(8.315*(273.16(degC)+celsius))) 
}

FUNCTION betn(v(mV)) {
LOCAL zeta
  zeta=zetan+pw/(1+exp((v-tq)/qq))
  betn = exp(1.e-3*zeta*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16(degC)+celsius))) 
}

FUNCTION alpl(v(mV)) {
  alpl = exp(1.e-3*zetal*(v-vhalfl)*9.648e4/(8.315*(273.16(degC)+celsius))) 
}

FUNCTION betl(v(mV)) {
  betl = exp(1.e-3*zetal*gml*(v-vhalfl)*9.648e4/(8.315*(273.16(degC)+celsius))) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v)
        n' = (ninf - n)/taun
        l' =  (linf - l)/taul
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,qt
	qt=q10^((celsius-24(degC))/10(degC))
	a = alpn(v)
	ninf = 1/(1 + a)
	taun = betn(v)/(qt*a0n*(1+a))
	if (taun<nmin) {taun=nmin}
	a = alpl(v)
	linf = 1/(1+ a)
	taul = 0.26*(v+50(mV))/qtl
	if (taul<lmin/qtl) {taul=lmin/qtl}
}




