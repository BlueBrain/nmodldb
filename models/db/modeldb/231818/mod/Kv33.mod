
:Comment :
:from Beining et al (2016), "A novel comprehensive and consistent electrophysiologcal model of dentate granule cells"


NEURON	{
	SUFFIX Kv33
	USEION k READ ek WRITE ik
	RANGE gkbar, gk, ik
	GLOBAL Rinact
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gkbar = 0.0008 (S/cm2) 
	Rinact = 0.1
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gk	(S/cm2)
	mInf
	mTau (ms)
	hInf
	hTau (ms)
	am (/ms)
	bm (/ms)
	ah (/ms)
	bh (/ms)
}

STATE	{ 
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gk = gkbar*m*m*m*m*h
	ik = gk*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
	UNITSOFF 
	am = 7.344 /(1+exp(-0.0807*(v-61.6)))
	bm = 0.611 /(1+exp(0.08625*(v+33)))
	mTau = 1/(am + bm )
	mInf = am * mTau

	
	hInf = Rinact + (1-Rinact)/(1+exp(((v -(-29.7))/(25))))   : this is taken for the model
	
	ah = 0.0066/(1+exp(((v -(10))/(-6))))   : A is 1/tau at very positive values
	bh = 0.01/(1+exp(((v -( -20))/(8 ))))    : A is 1/tau at very negative values
	hTau = 1/(ah+bh)   :    % this is taken for the model

	UNITSON
}