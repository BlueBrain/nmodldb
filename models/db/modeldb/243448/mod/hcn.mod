TITLE HCN current for bladder small DRG neuron soma model  
: Adapted from Kouranova et al., 2008

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
		(S) = (siemens)
}
 

NEURON {
        SUFFIX hcn
		USEION h READ eh WRITE ih VALENCE 1
        RANGE gbarfast, gbarslow, g, ih
        RANGE mtauf, mtausl, minf
		THREADSAFE
}
 
PARAMETER {
	   gbarfast = 1.352e-5 (S/cm2)
	   gbarslow = 6.7615e-5(S/cm2)
       eh = -30 (mV)
}
 
STATE {
        mf msl
}
 
ASSIGNED {
	v (mV)
	
	g (S/cm2)
	ih (mA/cm2)
     
    minf
	mtauf (ms)
	mtausl (ms)
}
 

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = gbarfast*mf+gbarslow*msl
		ih = g*(v - eh)
}
 

INITIAL {
	rates(v)
	mf = minf
	msl = minf
	
}


DERIVATIVE states {  
        rates(v)
        mf'  =  (minf-mf)/mtauf
		msl' =  (minf-msl)/mtausl
}


PROCEDURE rates(v(mV)) {  
		LOCAL q10
UNITSOFF
		minf = 1/(1+exp((v+87.2)/9.7)) : Kouronova 2008

		if (v < -70){
			mtauf = 250 + 12*exp((v+240)/50)
		}
		else{
			mtauf = 140 + 50*exp((v+25)/-20)
		}
		
		if (v < -70){
			mtausl = 2500 + 100*exp((v+240)/50)
		}
		else{
			mtausl = 300 + 542*exp((v+25)/-20)
		}
}
 
 
UNITSON
