TITLE calcium dynamics modeled by simple first-order decay 

COMMENT Equations from
   McCormick DA, Huguenard JR (1992) A model of the electrophysiological
   properties of thalamocortical relay neurons. J Neurophys 68(4):
   1384-1400.

>< Written by Arthur Houweling for MyFirstNEURON.
ENDCOMMENT

NEURON {
        SUFFIX cadyn
        USEION ca READ cai, ica WRITE cai 
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F    = (faraday) (coul)
}

PARAMETER {
	depth= .1 (micron)	: depth of shell for diffusion
        tau= 1 (ms)		
	cainf= 5e-5 (mM)	
}

ASSIGNED {
	ica	(mA/cm2)
	diam	(micron)
	A	(/coul/cm)
}

STATE { cai (mM) <1e-5> }

BREAKPOINT { 
        SOLVE states METHOD cnexp
}

DERIVATIVE states {
	cai'= -A* ica- (cai- cainf)/ tau
}

INITIAL {
	A= (1e4)* diam/ (2* F* depth* (diam- depth))
}

