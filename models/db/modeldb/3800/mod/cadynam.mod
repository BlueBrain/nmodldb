COMMENT
	 Cardiac intracellular calcium accumulation, modified from Courtemanche et al Am J Physiol 1998 275:H301 
ENDCOMMENT

NEURON {
	SUFFIX Cadynam
	USEION ca READ ica, cai WRITE cai
	USEION cu READ cui VALENCE 2
	USEION cr READ icr, cri  VALENCE 2
	RANGE Vi, Vrel,Vup, Trpn,Kmtrpn,Cmdn, Kmcmdn, cupmax, Kup
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F = (faraday) (coulombs)
}
ASSIGNED{
	icr 		(mA/cm2)
	ica    		(mA/cm2)
	cui		(mm)
	cri              (mM)
}
PARAMETER {
	Vi =  13668e-12   (cm3)
	Vup = 1109.52e-12 (cm3)
	Vrel = 96.48e-12   (cm3)
	Trpn = 0.07    (mM) 
	Kmtrpn = 0.0005  (mM)
	Cmdn = 0.05       (mM)
	Kmcmdn = 0.00238 (mM)
	Iupmax = 0.005   (mM/ms)
	cupmax = 15      (mM)
	   Kup = 0.00092  (mM)
	
}
 

STATE {
	cai START 0.0001 (mM)  <1e-6>
}
LOCAL ViF, VuP, VeR
INITIAL {
	VERBATIM
	cai = _ion_cai;
	
	ENDVERBATIM
 	ViF = Vi*F*2e4
	VuP = Vup/Vi
	VeR = Vrel/Vi
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 
	
	cai' = (-0.5*ica/(ViF) + VuP*(cui*Iupmax/cupmax - Iupmax/(1 + (Kup/cai))) + 0.5*icr*VeR)/(1 + Trpn*Kmtrpn/((cai + Kmtrpn)*(cai + Kmtrpn)) + Cmdn*Kmcmdn/((cai + Kmcmdn)*(cai + Kmcmdn)))
}
