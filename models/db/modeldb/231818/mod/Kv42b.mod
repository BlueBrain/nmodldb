TITLE Kv4.2 with auxilliary subunits

COMMENT
This is the Kv4.2 model of Amarillo et al 2008: Ternary Kv4.2 channels recapitulate voltage-dependent inactivation kinetics of A-type K+ channels in cerebellar granule neurons. The Journal of Physiology
: implemented by M.Beining; Beining et al (2016), "A novel comprehensive and consistent electrophysiologcal model of dentate granule cells"


ENDCOMMENT

NEURON {
	SUFFIX Kv42b
	USEION k READ ek WRITE ik
    RANGE  ik, gk, gkbar
	GLOBAL 	f,	g,a0,za,b0,zb, c0, zc, d0, zd, e0, ze, f0, zf ,kci,kic, kappa1, lambda1, kappa2, lambda2, vshift
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
    (S)  = (siemens)
	
	(molar) = (1/liter)
	(mM) = (millimolar)
	(uM) = (micromolar)
	FARADAY = (faraday) (kilocoulombs)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v 		(mV)
	gkbar  = 0.00015                (S/cm2) : to be fitted     	
	f = 0.045  : The allosteric factor f defines the voltage-dependent coupling between activation and inactivation pathways for closed channels. deactivation
	g = 1.03
	a0 = 1.589				(/ms)
	za = 0.64 				
	b0 = 0.0184		(/ms)
	zb = -1.31		
	c0 = 6.668			(/ms)
	zc = 0.15
	d0 = 2.381			(/ms)
	zd = -1.21
	e0 = 0.503			(/ms)
	ze = 0.07		
	f0 = 0.174			(/ms)
	zf = -0.25	
	kci = 0.047			(/ms)
	kic = 0.00003			(/ms)
	kappa1 = 0.229				(/ms)
	lambda1 = 0.151			(/ms)
	kappa2 = 0.0487			(/ms)
	lambda2 = 0.0065		(/ms)
	vshift = 0 					(mV)
}

STATE {
        O C0 C1 C2 C3 C4 I0 I1 I2 I3 I4  IO1 IO2 I5 C5
}

ASSIGNED {
        ik                             (mA/cm2)
        gk                            (S/cm2)
        ek                            (mV)
		alpha   					(/ms)
		beta						(/ms)
		gamma						(/ms)
		delta   					(/ms)
		epsilon 					(/ms)
		phi							(/ms)
		celsius (degC)
}

INITIAL {
	rate(v)
	SOLVE kin STEADYSTATE sparse
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	gk = gkbar * O
    ik = gk * ( v - ek )
}


KINETIC kin {
:LOCAL alpha2, alpha3, beta2, beta3
rate(v)

~ I0 <-> I1 (4*alpha/f,beta*f)
~ I1 <-> I2 (3*alpha/f,2*beta*f)
~ I2 <-> I3 (2*alpha/f,3*beta*f)
~ I3 <-> I4 (alpha/f,4*beta*f)
~ I4 <-> I5 (gamma*g,delta/g)

~ C0 <-> I0 (kci*f^4,kic/f^4)
~ C1 <-> I1 (kci*f^3,kic/f^3)
~ C2 <-> I2 (kci*f^2,kic/f^2)
~ C3 <-> I3 (kci*f,kic/f)
~ C4 <-> I4 (kci,kic)
~ C5 <-> I5 (kci*g,kic/g)

~ C0 <-> C1 (4*alpha,beta)
~ C1 <-> C2 (3*alpha,2*beta)
~ C2 <-> C3 (2*alpha,3*beta)
~ C3 <-> C4 (alpha,4*beta)
~ C4 <-> C5 (gamma,delta)
~ C5 <-> O (epsilon,phi)
~ O <-> IO1 (kappa1,lambda1)
~ IO1 <-> IO2 (kappa2,lambda2)


CONSERVE O + C0 + C1 + C2 + C3 + C4  + I0 + I1 + I2 + I3 + I4  + IO1 + IO2  + I5 + C5 = 1 
}


PROCEDURE rate(v (mV)) { :callable from hoc
	alpha = exponential(a0,za,v,vshift)
	beta = exponential(b0,zb,v,vshift)
	gamma = exponential(c0,zc,v,vshift)
	delta = exponential(d0,zd,v,vshift)
	epsilon  = exponential(e0,ze,v,vshift)
	phi  = exponential(f0,zf,v,vshift)
}


FUNCTION exponential(A(/ms), z , v (mV), D (mV)) (/ms) {
	exponential = A* exp(z*(v-D)*FARADAY/(R*(celsius+273.15)))
}