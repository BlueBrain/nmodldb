TITLE K-A channel from Beck Ficker and Heinemann (1992)
: M.Migliore 2001

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
	celsius		(degC)
	gkabar=.005 (mho/cm2)
        vhalfn=-73.1   (mV)
        vl=-73.1   (mV)
	vn=11	(mV)
	kn=3
	th=-55	(mV)
        vhalfl=-73.1   (mV)
        a0l=0.02      (/ms)
        a0n=0.3    (/ms)
        zetan=-1.5    (1)
        zetal=2    (1)
        gmn=0.7   (1)
        gml=0.65   (1)
	lmin=7.5  (mS)
	nmin=0.5  (mS)
	q10=3
	ek
}


NEURON {
	SUFFIX ka
	USEION k READ ek WRITE ik
        RANGE gkabar,gka
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
	gka = gkabar*n*l
	ik = gka*(v-ek)

}


FUNCTION alpn(v(mV)) {
  alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV)) {
  betn = exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alpl(v(mV)) {
  alpl = exp(1.e-3*zetal*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betl(v(mV)) {
  betl = exp(1.e-3*zetal*gml*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v)
        n' = (ninf - n)/taun
        l' =  (linf - l)/taul
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-22)/10)
	if (v<=th) {ninf=0} else {ninf = (2*(v-th)^kn)/((vn-th)^kn+ (v-th)^kn)}
        taun = betn(v)/(qt*a0n*(1+alpn(v)))
	if (taun<nmin/qt) {taun=nmin/qt}
        linf = 1/(1+ exp((vl-v)/-6.3))
        taul = betl(v)/(qt*a0l*(1+alpl(v)))
	if (taul<lmin/qt) {taul=lmin/qt}
}
