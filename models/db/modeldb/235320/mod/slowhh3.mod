
TITLE  squid sodium, potassium delayed rectifier, and potassium A channels
: Anna's version adding some shift values
	
 
UNITS {
        (molar) = (1/liter)
        (S) = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
         F = (faraday) (coulomb)
         R = (mole k) (mV-coulomb/degC)
        (mM) =  (millimolar)
}
 
NEURON {
        SUFFIX hh3
        USEION na READ nai WRITE ina
        USEION k WRITE ik
        RANGE  sshift, shift,gnabar,gkhhbar,gkabar,ina,ikhh,ika,ik,ena,gk,gna,gka
        GLOBAL minf,hinf,ninf
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v   (mV)
        dt  (ms)
	nai (mM)
	celsius = 35.0 (degC)
        gnabar  = 550.0e-6 (S/cm2)
        gkhhbar = 665.0e-6 (S/cm2)
        gkabar  = 266.0e-6  (S/cm2)
        ek  = -90.0  (mV)
        nao =  145  (mM)
        qv = 60.0 (mV)
        qs = 5.0  (1)
        shift=0 (mV)
        sshift=0 (mV)
 	
}
 
STATE {
        m <1e-4> h <1e-4> n <1e-4> p <1e-4> q <1e-4> 
}
 
ASSIGNED {
        ina (mA/cm2)
        ik (mA/cm2)
        ika (mA/cm2)
        ikhh (mA/cm2)
        ena (mV)
	minf hinf ninf qinf pinf 
	gna 
	gk 
	gka
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
	nai = 0.81
        ena = R*(celsius+273.15)/F*log(nao/nai)
	: ena=77.0
        ina = gnabar*m*m*m*h*(v - ena)
        ikhh = gkhhbar*n*n*n*(v - ek)      
        ika = gkabar*p*p*p*q*(v - ek)      
        ik = ika + ikhh
	gk = n*n*n
	gna = m*m*m*h
	gka = p*p*p*q
}
 
UNITSOFF
 
INITIAL {
        m = boltz(v,-34.6-shift,6.0-sshift)
        h = boltz(v,-56.8-shift,-7.8)
        n = boltz(v,-25.0,12.0)
        p = boltz(v,-43.0,24.0)
        q = boltz(v,-qv,-qs)
}

DERIVATIVE states {  :Computes state variables m, h, and n 
LOCAL minf,hinf,ninf,pinf,qinf,mtau,htau,ntau,ptau,qtau
        minf = boltz(v,-34.6-shift,6.0-sshift)
        hinf = boltz(v,-56.8-shift,-7.8)
        ninf = boltz(v,-25,12.0)
        pinf = boltz(v,-43.0,24.0)
        qinf = boltz(v,-qv,-qs)
        mtau = 1*boltz(v,-45.0,-1.5) - 1*boltz(v,-65.0,-0.5) +0.04
        htau = 3*56.0*boltz(v,-29.0,-4.5) - 3*56.0*boltz(v,-49.0,-2.0) +2.0
        ntau = 19.0/(1.0 + exp((v+ 39.0)/8.0)) - 0.0/(1.0 + exp((v+ 59.0)/20.0))+1.0
        ptau = 2.*exp(-(v+50)*(v+50)/550)+1.1
        qtau = 20.0
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
        p' = (pinf-p)/ptau
        q' = (qinf-q)/qtau
}
 
 
 
FUNCTION boltz(x,y,z) {
                boltz = 1/(1 + exp(-(x - y)/z))
}
 
UNITSON

