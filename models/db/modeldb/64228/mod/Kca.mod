:Calcium Activated Potassium Channels

NEURON 
{
	SUFFIX Kca
	USEION Ca READ Cai VALENCE 2
	USEION Kca WRITE iKca VALENCE 1
	RANGE infmKcaV,taumKcaV,eKca,gKcabar,Cahalf
}

UNITS
{
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millimho)
	(mol)= (1)
	(M)  = (mol/liter)
	(uM) = (micro M)
}

PARAMETER
{
       
 
       :Ca-dependent K current
       eKca=-80 (mV)
       gKcabar = 5 (mS/cm2)
       Cai
	Cahalf=1 (uM)	 
 

}

STATE
{

	:mKcaV
	mKcaCa
	
}

ASSIGNED
{
	iKca (mA/cm^2)
	v (mV)
           :Ca-dependent potassium channel, Kca
	:infmKcaV
	:taumKcaV  (ms)
	gKca (mho/cm2)

}

INITIAL
{      LOCAL Cas
	:rate(v)
	Cas=Cai*1000 :uM

	:mKcaV= infmKcaV
        :Cas=Cai
        mKcaCa=1/(1+(Cahalf/Cas)^4)
}




BREAKPOINT
{
	LOCAL Cas
	:SOLVE states METHOD cnexp
	Cas=Cai*1000 :uM
	mKcaCa=1/(1+(Cahalf/Cas)^4 )
	:gKca=(0.001)*gKcabar*(mKcaV^2)*mKcaCa
        gKca=(0.001)*gKcabar*mKcaCa^4
	iKca=gKca*(v-eKca) 
	: the current is in the unit of mA/cm2
	
}

