
TITLE  squid sodium, potassium delayed rectifier, and potassium A channels

	
 
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
	gk = n*n*n
	gna = m*m*m*h
	gka = p*q
        ina = gnabar*gna*(v - ena)
        ikhh = gkhhbar*gk*(v - ek)      
        ika = gkabar*gka*(v - ek) 
        ik = ika + ikhh
}
 
UNITSOFF
 
INITIAL {
        m = boltz(v,-30.0-shift,9.0-sshift)
        h = boltz(v,-35.0-shift,-10.0)
        n = boltz(v,-25.0,2.0)
        p = boltz(v,-27.0,6.5)
        q = 1:boltz(v,-64.0,-4.9)
}

DERIVATIVE states {  :Computes state variables m, h, and n 
LOCAL minf,hinf,ninf,pinf,qinf,mtau,htau,ntau,ptau,qtau
        minf = boltz(v,-30.0-shift,9.0-sshift)
        hinf = boltz(v,-35.0-shift,-10.0)
        ninf = boltz(v,-25,2.0)
        pinf = boltz(v,-27.0,6.5)
        qinf = boltz(v,-64.0,-4.9)
        mtau = 1*boltz(v,-45.0,-1.5) - 0.5*boltz(v,-65.0,-0.5) +0.01
        htau = 3*56.0*boltz(v,-29.0,-4.5) - 3*56.0*boltz(v,-49.0,-2.0) +0.4
       
	if (v > -40){
        ntau = 1.0 + 10.0*exp(-(log(1.0+0.05*(v+40))/0.05*log(1.0+0.05*(v+40))/0.05)/100) 
        } 
	else {
        ntau = 7 + 4*exp(-((v+40)/10.0)*((v+40)/10.0))
	}
              
	ptau = (1.029 + 4.83*boltz(v,-56.7,-6.22))
        qtau = (58.6+117.57*boltz(v,-68.5,-5.95))
	
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

