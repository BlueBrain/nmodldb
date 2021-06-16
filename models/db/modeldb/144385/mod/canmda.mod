TITLE Ca current through NMDA receptors 

: We use this workaround mechanism to calculate the Ca current through the NMDA receptors 
: separatly from the non specific ion current through the NMDA receptors in the nmda.mod file
: It contains:
: 
: 1.A mechanism to caculate the Ca current through the NMDA receptor
:   the Ca current through the NMDA receptor is added to the total Ca current "ica(mA/cm2)" 
:
: 2.A balance current "i_canmda(mA/cm2)" (the NONSPECIFIC_CURRENT i in the 
:   code above) to the Ca current through NMDA receptors (an inward current) 
:   The balance current is needed because it has already been caculated once as a part of the 
:   total current through NMDA receptors "i" in the "nmda.mod"
: 
: 3.Area (spine head surface area)is declared as a Global variable, and will be used in ampa.mod, nmda.mod, car.mod.
:
: 4.ampa, nmda and R_type current are all sent to this file as current density with the same direction of i_canmda. 
:   The itotal is just the sum of Inmda Iampa and I R_type
:
: Written by Lei Tian on 04/12/06 

NEURON {
	SUFFIX canmda 		
		:will be given to the variables in this file as their family name 
	
	USEION ca WRITE ica
	NONSPECIFIC_CURRENT i 
	RANGE g, i, mg, inmda, gnmda, iampa, gampa, itotal, irtype, Pca, P, f
	
	GLOBAL Area			
		:global varible, will be read by other files as a external one
	
	EXTERNAL i2_nmda, g2_nmda, i2_ampa, g2_ampa, irtype_car
		:declare the external variables which has been declared as Global ones in nmda.mod, ampa.mod and car.mod
	}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {                     : parameters that can be entered when function is called in cell-setup
        dt			(ms)
       	
		mg   = 1	(mM)        :Mg++ concentration
			
		Area = 1.11e-8  (cm2)	:spine head area 1.11e-8  (cm2)
		k = 1e-06   (mA/nA)		:transform the current from in 'nA' to in 'mA'

		P           (cm/s/uS) 	:a factor to convert NMDA conductance to permeability by considering the fraction of ca current at -65mV of NMDAr is about 10% normailize it at -65mV
}  


ASSIGNED {	: parameters needed to solve DE
	ica (mA/cm2)	:calcium current, which will be add to the total Ca current together with ica in 'car.mod'
	v (mV)          :spine head membrane potential
	i (mA/cm2)		:balance current to the ica through NMDA
	g (uS)          :conductance of nmda(not include the effect of Mg block)
	Pca (cm/s)		:Ca permeability of NMDA, it's obtained from gnmda by multiplied with P=0.1*gnmda*(v-e_nmda)/GHK at -65mV
	
	inmda (mA/cm2)	:equal to -i2_nmda which is the total nmda current's density, the direction is changed to be easier compared with i_canmda in this file.  
	gnmda	(uS)	:cunduction of nmda(include the Mg effect), to be easier plot out by just click the 'plot what'button 
	iampa	(mA/cm2):total current of ampa, the direction is changed to be easier compared with i_canmda in this file.  
	gampa	(uS)	:cunductance of ampa
	itotal (mA/cm2)	:total current flow into spinehead (only the aciviated channel current is considered),the direction is chosen the same as i_canmda in this file.  
	irtype (mA/cm2)	:r_type current
	f               :Ca current fraction in nmda current
}

INITIAL {

	P  = (1-exp(-65*-0.0755))/(10*Area*14564*(50e-09-(2e-03*exp(-65*-0.0755))))*k	:converting conductance to permaebility 
}


BREAKPOINT {
	g = g2_nmda	:[uS]
	Pca = P*g	:[cm/s]
	ica = Pca*14564*v*(50e-09-(2e-03*exp(v*-0.0755)))/(1-exp(v*-0.0755))*1/(1+(exp(0.08(/mV) * -v)*(mg / 0.69)))	:ca current density through NMDAr in [mA/cm2]
	i = -Pca*14564*v*(50e-09-(2e-03*exp(v*-0.0755)))/(1-exp(v*-0.0755))*1/(1+(exp(0.08(/mV) * -v)*(mg / 0.69)))	:balance current density of ca current through nmda
	
	:14564=(z^2*F^2)/(R*T); -0.0755 = -z*F/RT in [1/mV] where z=2,F=96500 in[C/mol], R=8.31 in[J/K*mol], T=308 in[K] 
	:and everything should be normalizied to [mV], 0.088 and 0.7474 is from our blocking experiment data fitting.
	

	gnmda=g2_nmda*1/(1+(exp(0.08(/mV) * -v)*(mg / 0.69)))	:[uS]cunduction of nmda(include the Mg effect), to be easier plot out by just click the 'plot what'button 
	gampa=g2_ampa	:[uS]total current of ampa, the direction is changed to be easier compared with i_canmda in this file
	inmda=-i2_nmda	:equal to -i2_nmda which is the total nmda current's density, the direction is changed to be easier compared with i_canmda in this file
	iampa=-i2_ampa	:total current of ampa, the direction is changed to be easier compared with i_canmda in this file.  
	
	irtype=irtype_car	:R-type current,the direction is chosen to be easier compared with i_canmda in this file.  
	itotal=i2_nmda+i2_ampa+irtype_car		:total current flow into spinehead (only the aciviated channel current is considered),the direction is chosen the same as i_canmda in this file. 
	f=i/inmda		:Ca current fraction in nmda current
	
}

	
