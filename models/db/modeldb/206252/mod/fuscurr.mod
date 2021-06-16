TITLE fuscurr.mod   DCN pyramidal cell model  
 
COMMENT

Fast transient sodium (Na), noninactivating potassium (K) from Kanold and Manis (2001)
Ih from Destexhe et al 1993, Leao et al (2012) and Nagtegaal and Borst (2010)
Persisten sodium (NaP) and sensitive inwardly rectifying potassium (Kir) currents
based on mathematical model from Leao et al (2012)

ENDCOMMENT
 
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

 
? interface
NEURON {
        SUFFIX fuscurr
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
	USEION h READ eh WRITE ih VALENCE 1
        USEION nap READ enap WRITE inap VALENCE 1
        USEION kir READ ekir WRITE ikir VALENCE 1
	NONSPECIFIC_CURRENT il

	RANGE gna, gk, minf, hinf, ninf		: Transient sodium channels and potassium delayed rectifier
	GLOBAL mtau, htau, ntau, gnabar, gkbar	: time constants for sodium channels and delayed rectifier
		
	RANGE gh, kh_m_inf, kh_n_inf, aih	: Ih channels
	GLOBAL kh_m_tau, kh_n_tau, ghbar	: time constants for Ih channels

	RANGE gnap, gkir, minf_p, ninf_ir		: Persistent Na channels and KIR channels 
	GLOBAL mtau_p, ntau_ir, gnapbar, gkirbar	: time constants for Na channels and K channels 

	GLOBAL gl, el					: leak (linear, ohmic with reversal at el)
}

 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius (degC)
        dt (ms)
	ek (mV)
        ena (mV)
        gnabar = 0.08 (mho/cm2)	<0,1e9>
        gkbar = 0.02 (mho/cm2)	<0,1e9>
        ghbar = 0.00054 (mho/cm2) <0,1e9>
	eh = -43 (mV)
	gl = 0.000150 (mho/cm2)	<0,1e9>
        el = -51.32 (mV) 
	mtau = 0.05 (ms) <0.01,100>
	htau = 0.5 (ms) <0.1,100>
	ntau = 0.5 (ms) <0.1,100>
	ekir (mV)
        enap (mV)
        gnapbar = 0.0001 (mho/cm2) <0,1e9>
        gkirbar = 0.0005 (mho/cm2)	<0,1e9>
	mtau_p = 5.3 (ms) <0,100>
	ntau_ir = 0.5 (ms) <0,100>
}


STATE {
        m h n khm khn mp nir
}
 
ASSIGNED {
    gna (mho/cm2)
    gk (mho/cm2)
    gh (mho/cm2)
    gnap (mho/cm2)
    gkir (mho/cm2)

    ina (mA/cm2)
    ik (mA/cm2)
    ih (mA/cm2)
    il (mA/cm2)
    inap (mA/cm2)
    ikir (mA/cm2)

    minf hinf
    ninf 
    kh_m_inf kh_n_inf
    kh_m_tau kh_n_tau	
    aih
    minf_p 
    ninf_ir 
}


LOCAL mexp, hexp, nexp, kh_m_exp, kh_n_exp, mexp_p, nexp_ir
 
? currents
BREAKPOINT {
    SOLVE states METHOD cnexp

    	gna = gnabar*m*m*h
	ina = gna*(v - ena)

   	gk = gkbar*n*n
	ik = gk*(v - ek)

	aih = 0.5*khm+0.5*khn
	gh = ghbar*aih
 	ih = gh*(v - eh)
	
	il = gl*(v - el)

    	gnap = gnapbar*mp*mp*mp
	inap = gnap*(v - enap)

   	gkir = gkirbar*nir
	ikir = gkir*(v - ekir)
}
? currents

UNITSOFF 
 

INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
	khm = kh_m_inf
	khn = kh_n_inf
	mp = minf_p
	nir = ninf_ir
}


? states
DERIVATIVE states {  
    rates(v)
    m' = (minf - m) / mtau
    h' = (hinf - h) / htau
    n' = (ninf - n) / ntau
    khm' = (kh_m_inf - khm) / kh_m_tau
    khn' = (kh_n_inf - khn) / kh_n_tau
    mp' = (minf_p - mp) / mtau_p
    nir' = (ninf_ir - nir) / ntau_ir
}
 
LOCAL q10

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
	                      :Call once from HOC to initialize inf at resting v.
	LOCAL  alpha, beta, sum
	TABLE minf, mtau, hinf, htau, ninf, ntau, kh_m_inf, kh_n_inf, kh_m_tau, kh_n_tau, minf_p, ninf_ir DEPEND celsius FROM -200 TO 100 WITH 400

UNITSOFF
	q10 = 3^((celsius - 32)/10)

	:"m" sodium activation system
		minf = na_m(v)


	:"h" sodium inactivation system
		hinf = na_h(v)

	:"n" potassium activation system
        ninf = kd_m(v)

	:"kh" adaptation of Destexhe hyp-activated cation current by Patrick Kanold
		kh_m_inf = kh_m(v) 
		kh_n_inf = kh_n(v)
		kh_m_tau = kh_mt(v)
		kh_n_tau = kh_nt(v) 

	:"m" sodium activation system
	minf_p = nap_m(v)

	:"n" potassium activation system
        ninf_ir = kird_m(v)
}

 
FUNCTION na_m(x) { : sodium activation
	na_m = 1/(1+exp(-(x+38)/3.0))	
}

FUNCTION na_h(x) { : sodium inactivation
	na_h = 1/(1+exp((x+43)/3.0))	
}

FUNCTION kd_m(x) { : potassium inactivation
	kd_m = 1/(1+exp(-(x+40)/3))		
}

FUNCTION kh_m(x) { : h activation
	kh_m = 1/(1+exp((x+87)/8.9))
}

FUNCTION kh_n(x) { : h activation
	kh_n = 1/(1+exp((x+87)/8.9)) : same as kh_m
}

FUNCTION kh_mt(v) { : h fast time constant
	kh_mt = 100 + exp((v+183.6)/30.48)
}

FUNCTION kh_nt(v) { : h slow time constant
	kh_nt = 700 + exp((v+188.6)/11.2)/(1+exp((v+105)/5.5))
}

FUNCTION nap_m(x) { : persitent sodium activation
	nap_m = 1/(1+exp(-(x+55)/2.8))	
}

FUNCTION kird_m(x) { : KIR activation
	kird_m = 1/(1+exp((x+85.48)/12)) 
}

UNITSON
