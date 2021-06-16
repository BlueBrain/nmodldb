: Template for transient Na+ -type of Hodgkin-Huxley channel
: Formalism from Buhry et al. 2013: "Global parameter estimation of an Hodgkin-Huxley formalism using membrane voltage
:                                    recordings: Application to neuro-mimetic analog integrated circuits", Neurocomputing 81 (2012) 75-85
: NEURON implementation by Tuomo Maki-Marttunen, 2013-2017 (CC-BY 4.0)

NEURON {
       SUFFIX I1
       NONSPECIFIC_CURRENT i
       RANGE i,E,g,m,h,Voffa,Vsloa,taua,Voffi,Vsloi,taui
}

UNITS { 
      (mV) = (millivolt)
      (S) = (siemens) 
}

STATE { m h }


PARAMETER {
	  g = 30 (S/cm2)
	  E = 55 (mV)
          Voffa = -46 (mV)
          Vsloa = 6.5 (mV)
          taua = 0.08 (ms)
          Voffi = -50 (mV)
          Vsloi = 10 (mV)
          taui = 0.3 (ms)
}

ASSIGNED {
	 i (mA/cm2)
	 v (mV)
}

BREAKPOINT {
	   SOLVE states METHOD cnexp
	   i = g*(m^3*h*(v-E))
}

INITIAL {
	m = 1/(1+exp(-(v-Voffa)/Vsloa))
        h = 1/(1+exp((v-Voffi)/Vsloi))
}

DERIVATIVE states {
	   m' = 1/taua*(1/(1+exp(-(v-Voffa)/Vsloa))-m)
	   h' = 1/taui*(1/(1+exp((v-Voffi)/Vsloi))-h)
}

