TITLE HH style channels for spiking retinal ganglion cells
:
: Modified from Fohlmeister et al, 1990, Brain Res 510, 343-345
: by TJ Velte March 17, 1995
: must be used with calcium pump mechanism, i.e. capump.mod
:
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX spike
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	USEION ca READ cai, eca, cao WRITE ica
	RANGE gnabar, gkbar, gabar, gcabar, gkcbar,gkc
	RANGE m_inf, h_inf, n_inf, p_inf, q_inf, c_inf
	RANGE tau_m, tau_h, tau_n, tau_p, tau_q, tau_c,tau_slowh
	RANGE m_exp, h_exp, n_exp, p_exp, q_exp, c_exp,slowh_exp
    RANGE idrk, iak, icak
    GLOBAL vshift_ca,ca50,slowh
	RANGE vshift_na_k
	GLOBAL tausv,tausd,tausb,Nm,Nh,Nsh
	RANGE nkm,  gkmbar
	RANGE km_inf, tau_km

}


UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	gnabar	= 0.04	(mho/cm2)
	gkbar	= 0.012 (mho/cm2)
	gabar	= 0.036	(mho/cm2)
	gcabar	= 0.002	(mho/cm2)
	gkcbar	= 0.00005 (mho/cm2)
	ena	= 65	(mV)
	ek	= -100	(mV)
	eca		(mV)
	cao	= 1.8	(mM)
	cai     = 0.0001 (mM)
	dt              (ms)
	v               (mV)
	vshift_na_k=5	(mV)
	vshift_ca=5		(mV)
	ca50=0.2		(mM)
	slowh=1
	gkmbar = .002   	(mho/cm2)	: 0.03 mho/cm2
	tausb=0.5
	taus=50
	tausv=30
	tausd=1
	Nm=3
	Nh=1
	Nsh=1
}

STATE {
	m h n p q c nkm sh
}

INITIAL {
: The initial values were determined at a resting value of -66.3232 mV in a single-compartment
:	m = 0.0155
:	h = 0.9399
:	n = 0.0768
:	p = 0.0398
:	q = 0.4526
:	c = 0.0016
: at -60 mV
        m = 0.0345
        h = 0.8594
		sh=	1
        n = 0.1213
        p = 0.0862
        q = 0.2534
        c = 0.0038
		nkm=0.0344

}

ASSIGNED {
	ina	(mA/cm2)
	ik	(mA/cm2)
         idrk    (mA/cm2)
         iak     (mA/cm2)
         icak    (mA/cm2)
		 ikm	 (mA/cm2)
	ica	(mA/cm2)
	m_inf h_inf n_inf p_inf q_inf c_inf km_inf
	tau_m tau_h tau_n tau_p tau_q tau_c tau_km :tau_slowh
	m_exp h_exp n_exp p_exp q_exp c_exp km_exp :slowh_exp
	gkc

}

BREAKPOINT {
	SOLVE states
	ina = gnabar * m^Nm*h^Nh*sh^Nsh * (v - ena)
	idrk = gkbar * n*n*n*n * (v - ek)
	iak =  gabar * p*p*p*q * (v - ek)
	gkc=((cai / ca50)/ (1 + (cai / ca50)))
	icak = gkcbar * gkc * (v - ek)
	ikm = 2.3*gkmbar*nkm*(v-ek)
	ik = idrk + iak + icak+ikm
	ica = gcabar * c*c*c * (v - eca)

}

PROCEDURE states() {	: exact when v held constant
	LOCAL sigmas
	evaluate_fct_nak(v+vshift_na_k)
	evaluate_fct_ca(v+vshift_ca)
	m = m + m_exp * (m_inf - m)
	h = h + h_exp * (h_inf - h)
	:sh = sh + slowh_exp * (h_inf - sh)
	sigmas=1/(1+exp((v+tausv+vshift_na_k)/tausd))
	sh = sh + (1 - exp(-dt/(taus*sigmas+tausb)))*(1 / (1 + exp((v + 44+vshift_na_k)/3))  - sh)

	n = n + n_exp * (n_inf - n)
	p = p + p_exp * (p_inf - p)
	q = q + q_exp * (q_inf - q)
	c = c + c_exp * (c_inf - c)
	nkm = nkm + km_exp*(km_inf-nkm)
	VERBATIM
	return 0;
	ENDVERBATIM

}

UNITSOFF
PROCEDURE evaluate_fct_ca(v(mV)) { LOCAL a,b
:CA channel
	a = (-0.3 * (v+13)) / ((exp(-0.1*(v+13))) - 1)
	b = 10 * (exp((-1*(v + 38))/18))
	tau_c = 1 / (a + b)
	c_inf = a * tau_c
	c_exp = 1 - exp(-dt/tau_c)	
}
PROCEDURE evaluate_fct_nak(v(mV)) { LOCAL a,b
	
:NA m
	a = (-0.6 * (v+30)) / ((exp(-0.1*(v+30))) - 1)
	b = 20 * (exp((-1*(v+55))/18))
	tau_m = 1 / (a + b)
	m_inf = a * tau_m

:NA h
	a = 0.4 * (exp((-1*(v+50))/20))
	b = 6 / ( 1 + exp(-0.1 *(v+20)))
	tau_h = 1 / (a + b)
	:tau_slowh=tau_h*slowh
	h_inf = a * tau_h

:K n (non-inactivating, delayed rectifier)
	a = (-0.02 * (v+40)) / ((exp(-0.1*(v+40))) - 1)
	b = 0.4 * (exp((-1*(v + 50))/80))
	tau_n = 1 / (a + b)
	n_inf = a * tau_n

:K (inactivating)
	a = (-0.006 * (v+90)) / ((exp(-0.1*(v+90))) - 1)
	b = 0.1 * (exp((-1*(v + 30))/10))
	tau_p = 1 / (a + b)
	p_inf = a * tau_p

	a = 0.04 * (exp((-1*(v+70))/20))
	b = 0.6 / ( 1 + exp(-0.1 *(v+40)))	
	tau_q = 1 / (a + b)
	q_inf = a * tau_q
:Km
    a = 0.001 * (v +30) / (1 - exp(-(v +30)/9))
    b = -0.001 * (v +30) / (1 - exp((v +30)/9))
    tau_km = 1/(a+b)
	km_inf = a*tau_km


: State vars to inifinity
	m_exp = 1 - exp(-dt/tau_m)
	h_exp = 1 - exp(-dt/tau_h)
	n_exp = 1 - exp(-dt/tau_n)
	p_exp = 1 - exp(-dt/tau_p)
	q_exp = 1 - exp(-dt/tau_q)
	km_exp= 1 - exp(-dt*2.3/tau_km)

}

UNITSON