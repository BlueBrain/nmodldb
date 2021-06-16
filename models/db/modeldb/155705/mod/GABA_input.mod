TITLE first_GABA_recept_model (first order kinetics )
COMMENT
This program just explore basics charactieristics of GABAergic synapses
in the frame of  a simple integrate and fire model. Modified from original of R. van Elburg and based on
Bush PC, Prince DA, Miller KD (1999),J Neurophysiol 82:1748-58

    This is mod-file describes a point process with first order kinetics, that can act as several 
    synapses described by first order kinetics.
     
	Author: 
	Taken from  model for pyramidal cells in J. Tegner, A. Compte, 
	X.J. Wang, J. Neurosci. 22(20): 9053-9062, 2002	
	X.J. Wang, J. Neurosci. 21(19):9587-9603

	Modifications:
	
	ID	Date		Authors				Email			Description
	M_001
	
ENDCOMMENT

NEURON {
	POINT_PROCESS Gaba_syn
	RANGE e, i, g
	RANGE tau_d, frac_rec 
	RANGE area_cell
	RANGE tau_1 : artificion pero realmente pertenece a gabafacil es el decay time dort
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA)   = (nanoamp)
	(mV)   = (millivolt)
	(uS)   = (microsiemens)
	(um)   = (micron)
}

PARAMETER {
	: e = -90 mV for inhibitory synapses,
	:     0 mV for excitatory
	:e = 	:-90 
	e=-80	(mV): value of boergers and koppel's paper



	  tau_1 = 3 (ms):ojo esto es un artificio, no se usa, es solo para incluir el archivo facilmente en hoc!!!!!!!!!!!
	: tau_d decay time, frac_rec usage fraction of receptors
		tau_d = 10 (ms) 	< 1e-9, 1e9 > : 
	    frac_rec = 0.9 (1) 	<0,1>
	
	: maximum coductance
	   g=1e-5 (uS/um2) : this must be multiplied by the area to give actual walue
                   : since the conductance density is  gie=5-7mS/cm2 gii=15-20mS/cm2
                   : Jensen et al. NeuroImage 26 (2005) 347-355
                   
    :cell surface area
        area_cell= 1 (um2)
}

ASSIGNED {
	v (mV)
	i (nA)
	g_eff (uS)
	
}

STATE {
	s
}

INITIAL {
	s=0 :
	g_eff=g*area_cell
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i =g_eff*s*(v - e) :
}

DERIVATIVE state {
	s' = 	-s/tau_d
}

NET_RECEIVE(weight,s_tp,tp(ms)) {:tp time of previous spike
	: Calculate current value single synapse state variable at t-epsilon
	UNITSOFF
	:printf("%g  %g\t",t,s_tp)
      s_tp =s_tp*exp(-(t-tp)/tau_d)
     :s_tp =s_tp*weight*exp(-(t-tp)/tau_d) 
    :printf("%g\t",s_tp)
	UNITSON
	
	: To make sure that we add the same amount to the exact single synapse state variable and the summed
	: state variable we should first update the summed variable and then the single synapse state variable
	 s=s + frac_rec*weight*(1-s_tp)	
    s_tp = s_tp + frac_rec*(1-s_tp)
	:printf("%g %g %g %g\n",weight,s_tp,s,frac_rec)
    tp=t   
}
