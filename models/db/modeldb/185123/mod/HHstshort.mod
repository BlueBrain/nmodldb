
TITLE Stochastic Hodgkin and Huxley model & M-type potassium & T-and L-type Calcium channels incorporating channel noise .

COMMENT

Based on - Accurate and fast simulation of channel noise in conductance-based model neurons. Linaro, D., Storace, M., and Giugliano, M. 
Added: 
	Km T L channels
	fixed minor bugs and grouped variables into arrays 

ENDCOMMENT

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
    (pS) = (picosiemens)
    (um) = (micrometer)
} : end UNITS


NEURON {
    SUFFIX HH
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT ileak
	RANGE eleak, gleak,NFleak
    RANGE gnabar, gkbar
    RANGE gna, gk
	RANGE nm,  gkmbar
	RANGE m_exp, h_exp, n_exp, km_exp
	
    GLOBAL gamma_na, gamma_k,gamma_km
	RANGE km_inf, tau_km
    RANGE m_inf, h_inf, n_inf
    RANGE tau_m, tau_h, tau_n
    RANGE tau_zn,tau_zm,tau_zkm
    RANGE var_zn,var_zm,var_zkm
    RANGE noise_zn,noise_zm,noise_zkm
    RANGE Nna, Nk,Nkm
    GLOBAL seed    ,hslow,hfast
    RANGE mu_zn,mu_zm,mu_zkm
	GLOBAL vshift

	GLOBAL taukm,NF
    THREADSAFE
} : end NEURON


PARAMETER {
    gnabar  = 0.12   (S/cm2)
    gkbar   = 0.036  (S/cm2)
	gkmbar = .002   	(S/cm2)
    gleak=0.00001		(S/cm2)
	eleak=-60 			(mV)
	NFleak=1
    gamma_na = 10  (pS)		: single channel sodium conductance
    gamma_k  = 10  (pS)		: single channel potassium conductance
    gamma_km  = 10  (pS)		: single channel potassium conductance
    seed = 1              : always use the same seed
	vshift=0				:Voltage shift of the recorded memebrane potential (to offset for liquid junction potential
	taukm=1					:speedup of Km channels
	NF=1					:Noise Factor (set to 0 to zero the noise part)
	hslow=100
	hfast=0.3
} : end PARAMETER


STATE {
    m h n nm 
    zn[3]
    zm[6]
	zkm
	zleak
} : end STATE


ASSIGNED {
    ina        (mA/cm2)
	ikdr	(mA/cm2)
	ikm	 (mA/cm2)
	ik		(mA/cm2)
	ileak	(mA/cm2)
    gna        (S/cm2)
    gk         (S/cm2)
	gkm         (S/cm2)
    ena        (mV)
    ek         (mV)
   
    dt         (ms)
    area       (um2)
    celsius    (degC)
    v          (mV)
        
    Nna  (1) : number of Na channels
    Nk   (1) : number of K channels
    Nkm   (1) : number of Km channels

    m_exp h_exp n_exp km_exp 
    m_inf h_inf n_inf km_inf 
    noise_zn[3]
    noise_zm[6]
	noise_zkm
    var_zn[3]
    var_zm[6]
	var_zkm
    tau_m (ms) tau_h (ms) tau_n (ms) tau_km (ms) 
    tau_zn[3]
    tau_zm[6]
	tau_zkm

    mu_zn[3]
    mu_zm[6]
	mu_zkm

} : end ASSIGNED

INITIAL {
    Nna = ceil(((1e-8)*area)*(gnabar)/((1e-12)*gamma_na))   : area in um2 -> 1e-8*area in cm2; gnabar in S/cm2; gamma_na in pS -> 1e-12*gamma_na in S
    Nk = ceil(((1e-8)*area)*(gkbar)/((1e-12)*gamma_k))   : area in um2 -> 1e-8*area in cm2; gkbar in S/cm2; gamma_k in pS -> 1e-12*gamma_k in S
    Nkm = ceil(((1e-8)*area)*(gkmbar)/((1e-12)*gamma_km))   : area in um2 -> 1e-8*area in cm2; gkbar in S/cm2; gamma_k in pS -> 1e-12*gamma_k in S
    
    rates(v)
    m = m_inf
    h = h_inf
    n = n_inf
	nm=km_inf
	
	zn[0] = 0
	zn[1] = 0
	zn[2] = 0
	zn[3] = 0

    zm[0] = 0.
    zm[1] = 0.
    zm[2] = 0.
    zm[3] = 0.
    zm[4] = 0.
    zm[5] = 0.
    zm[6] = 0.
	
	zkm=0
	
	
	zleak=0
    set_seed(seed)
} : end INITIAL


BREAKPOINT {
    SOLVE states
    gna = gnabar * (m*m*m*h + (zm[0]+zm[1]+zm[2]+zm[3]+zm[4]+zm[5]+zm[6]))
:    if (gna < 0) {
		:gna = 0
:    }
    ina = gna * (v - ena)
    gk = gkbar * (n*n*n*n + (zn[0]+zn[1]+zn[2]+zn[3]))
:   if (gk < 0) {
:		gk = 0
:    }
    ikdr  = gk * (v - ek)
	gkm=gkmbar*(nm+zkm)
:	if (gkm < 0) {
:		gkm= 0
:    }
	ikm = gkm*(v-ek)
	ik=ikm+ikdr
	

	ileak  = gleak*(1+zleak) * (v - eleak)
:	ileak  = (gleak) * (v - eleak)
} : end BREAKPOINT


PROCEDURE states() {
    rates(v+vshift)
	m = m + m_exp * (m_inf - m)
	h = h + h_exp * (h_inf - h)
    n = n + n_exp * (n_inf - n)
	nm = nm + km_exp*(km_inf-nm)


    VERBATIM
    return 0;
    ENDVERBATIM
} : end PROCEDURE states()


PROCEDURE rates(vm (mV)) { 
    LOCAL a,b,m3_inf,n4_inf,sum,one_minus_m,one_minus_h,one_minus_n,i
    
    UNITSOFF
    
    
 :NA m
	a =-.6 * vtrap((vm+30),-10)	
	b = 20 * (exp((-1*(vm+55))/18))
	tau_m = 1 / (a + b)
	m_inf = a * tau_m

    one_minus_m = 1. - m_inf
    m3_inf = m_inf*m_inf*m_inf

:NA h
	a = 0.4 * (exp((-1*(vm+50))/20))
	b = 6 / ( 1 + exp(-0.1 *(vm+20)))
	:b = 6 / ( 1 + exp(-(vm-50)/5))
	:tau_h = 1 / (a + b)
	:h_inf = a * tau_h
	:tau_h=1/(a+6 / ( 1 + exp(-(vm-50)/5)))
	:tau_h=hslow/(1+exp((vm+30)/4))+hfast
	tau_h=hslow/((1+exp((vm+30)/4))+(exp(-(vm+50)/2)))+hfast
	h_inf= 1/(1 + exp((vm + 44)/4))
	:h_inf = a / (a+b)

	one_minus_h = 1. - h_inf
    
    tau_zm[0] = tau_h
    tau_zm[1] = tau_m
    tau_zm[2] = tau_m/2
    tau_zm[3] = tau_m/3
    tau_zm[4] = tau_m*tau_h/(tau_m+tau_h)
    tau_zm[5] = tau_m*tau_h/(tau_m+2*tau_h)
    tau_zm[6] = tau_m*tau_h/(tau_m+3*tau_h)
    var_zm[0] = 1.0 / numchan(Nna) * m3_inf*m3_inf*h_inf * one_minus_h
    var_zm[1] = 3.0 / numchan(Nna) * m3_inf*m_inf*m_inf*h_inf*h_inf * one_minus_m
    var_zm[2] = 3.0 / numchan(Nna) * m3_inf*m_inf*h_inf*h_inf * one_minus_m*one_minus_m
    var_zm[3] = 1.0 / numchan(Nna) * m3_inf*h_inf*h_inf * one_minus_m*one_minus_m*one_minus_m
    var_zm[4] = 3.0 / numchan(Nna) * m3_inf*m_inf*m_inf*h_inf * one_minus_m*one_minus_h
    var_zm[5] = 3.0 / numchan(Nna) * m3_inf*m_inf*h_inf * one_minus_m*one_minus_m*one_minus_h
    var_zm[6] = 1.0 / numchan(Nna) * m3_inf*h_inf * one_minus_m*one_minus_m*one_minus_m*one_minus_h
	
	FROM i=0 TO 6 {
		mu_zm[i] = exp(-dt/tau_zm[i])
        noise_zm[i] = sqrt(var_zm[i] * (1-mu_zm[i]*mu_zm[i])) * normrand(0,NF)
		zm[i] = zm[i]*mu_zm[i] + noise_zm[i]
    }
  
   
:K n (non-inactivating, delayed rectifier)
	a = -0.02 * vtrap((vm+40),-10)
	b = 0.4 * (exp((-1*(vm + 50))/80))
	tau_n = 1 / (a + b)
	n_inf = a * tau_n
    one_minus_n = 1. - n_inf
    n4_inf = n_inf * n_inf * n_inf * n_inf
    tau_zn[0] = tau_n
    tau_zn[1] = tau_n/2
    tau_zn[2] = tau_n/3
    tau_zn[3] = tau_n/4
    var_zn[0] = 4.0/numchan(Nk) * n4_inf*n_inf*n_inf*n_inf * one_minus_n
    var_zn[1] = 6.0/numchan(Nk) * n4_inf*n_inf*n_inf * one_minus_n*one_minus_n
    var_zn[2] = 4.0/numchan(Nk) * n4_inf*n_inf * one_minus_n*one_minus_n*one_minus_n
    var_zn[3] = 1.0/numchan(Nk) * n4_inf * one_minus_n*one_minus_n*one_minus_n*one_minus_n

	FROM i=0 TO 3 {
		mu_zn[i] = exp(-dt/tau_zn[i])
        noise_zn[i] = sqrt(var_zn[i] * (1-mu_zn[i]*mu_zn[i])) * normrand(0,NF)
		zn[i] = zn[i]*mu_zn[i] + noise_zn[i]
    }
:Km
    a = -.001/taukm * vtrap((vm+30),-9)
    b =.001/taukm * vtrap((vm+30),9) 
    tau_km = 1/(a+b)
	km_inf = a*tau_km
	tau_zkm = tau_km
	var_zkm=km_inf*(1-km_inf)/numchan(Nkm)
	mu_zkm = exp(-dt/tau_zkm)
    noise_zkm = sqrt(var_zkm * (1-mu_zkm*mu_zkm)) * normrand(0,NF)
	zkm = zkm*mu_zkm + noise_zkm
	
	zleak=normrand(0,NFleak)
	m_exp = 1 - exp(-dt/tau_m)
	h_exp = 1 - exp(-dt/tau_h)
	n_exp = 1 - exp(-dt/tau_n)
	km_exp= 1 - exp(-dt/tau_km)	

	
   UNITSON
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
    if (fabs(exp(x/y) - 1) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{
        vtrap = x/(exp(x/y) - 1)
    }
}
FUNCTION mulnoise(mean,sd,power){
	LOCAL i,avg
	avg=1
	FROM i=1 TO power {
		avg=avg*normrand(mean,sd)
	}
	mulnoise=avg
}

FUNCTION numchan(Nchannels){
	if (Nchannels>0){
		numchan=(Nchannels)
	}else{
		numchan=1
	}
}