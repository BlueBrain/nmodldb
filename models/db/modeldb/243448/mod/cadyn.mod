TITLE Calcium Difffusion, Buffering, ER Mechs- SERCA, IP3R & CICR(RYR), and Mitochondrial Influx(MCU) and Eflux(MNCX) for bladder small DRG neuron soma model

: Author: Darshan Mandge (darshanmandge@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

NEURON{
	SUFFIX cadyn
	USEION ca READ ica, cai, cao WRITE cai, ica VALENCE 2 		: writing cai and ica
	USEION na READ nai VALENCE 1								: for mncx
	USEION caer READ caeri WRITE caeri VALENCE 2 				: caeri = internal ER calcium. WRITING caeri
	USEION camt READ camti WRITE camti VALENCE 2 				: camti = internal mitochondrial calcium
	USEION caip3r READ caip3ri WRITE caip3ri VALENCE 2			: caip3ri is the ca change by ca release by outermost shell's IP3Rs. It is read by CACC
	USEION ip3 READ ip3i VALENCE 1								: ip3i = ip3 conc. 
   
	GLOBAL DCa, cai0, caeri0, camti0							: Diffusion constant of Ca, calcium conc. in cytoplasm, ER and mitochondria 
	RANGE k1, k2, k3, k4, ica_pmp, pump0 			 		 	: Pump Parameters: RANGE as can be different for different segments and sections
	GLOBAL bbr													: Buffer Parameters: GLOBAL as they are props. of buffer and are constant         
	:RANGE kmdye, Bmdye, Dbufm, bbrdye							: Dye correction: fura-2 in O'mullane-2013
	
	RANGE vmaxsr, kpsr		                        			: SERCA parameters
	RANGE kactip3, konip3, kinhip3,  kip3, jmaxsr			  	: IP3 Parameters  
	RANGE ktcicr, kcicr, vcicr                           		: CICR Parameters
    RANGE vmcu, kmcu, nmcu, vncx, kna, kncx			    	   	: Mitochondrial Parameters: Uniporter and MNCX
	RANGE jer, jserca, jip3, jcicr								: Flux Parameters
	RANGE jmcu, jmncx, jmito

	RANGE Kmmt, Bmmt, fmmt										: Mitochondrial Buffer Paramters
	RANGE Kmer, Bmer, fmer										: ER Buffer Paramters
	
    GLOBAL vol
	THREADSAFE
}
DEFINE NANN  12    :IF YOU CHANGE THIS NUMBER DONT FORGET to change the NANN value in ip3dif.mod file.

UNITS{
	(mV)	= (millivolt)
	(um)    = (micron)
	(mM)    = (milli/liter)
    (nM)    = (nano/liter)
	(mA)    = (milliamp)
	F       = (faraday) (coulombs)
	PI      = (pi) (1)
	R 		= (k-mole) (joule/degC)
    (mol)   = (1) :The term mole cannot be used here because it is already defined in NEURON's units database as 6.022169*10^23
}

PARAMETER{
	diam	(um)
	L		(um)
:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	:Ca Parameters
	cai0 	= 136e-6	(mM)
	cao0 	= 2 		(mM)
	caeri0 	= 0.4		(mM)
   	camti0 	= 2e-4 		(mM)
	DCa 	= 0.6		(um2/ms)     : 0.6 in McHugh et al., 2004
:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	: Buffers
	bbr = 370 		:Zeilhofer 1996 (J. neurophysiology) endogenous buffer binding ratio = 370 (Rat DRG)
	
	: Dye paramters
	: Fura-2 O'Mullane et al., 2013 Data for bladder cai measurement
	: kmdye = 224e-6 (mM)		: Grynkiewicz et al., 1985 224 nM in Lu and Gold 2006, 2008
	: Bmdye = 5e-3   (mM)		: O'Mullane 2013 5 uM Fura-2
	: Dbufm = 0.1 (um2/ms)      : Blatter and Wier, 1990 
	
	: Indo-1 Data for cai 
	: kmdye = 250e-6 (mM)	    : 250 nM Grynkiewicz et al., 1985		
	: Bmdye = 100e-3 (mM)     	: 100 uM   Benham et al.,1992
	: Dbufm =  0.1 (um2/ms)     : Blatter and Wier, 1990
	
	: Fura-FF shutov
	: kmdye = 5500e-6 (mM)	    : 5.5 uM Shutov et al., 2013
	: Bmdye = 200e-3 (mM)     	: 200 uM Shutov et al., 2013
	: Dbufm =  0.075 (um2/ms)   : Assumed similar to Fura-2 and indo-1
	
	: mt-pericam
	: kmmitodye = 1.7e-3 (mM)	: 1.7 Nagai et al., 2001
	: Bmmitodye = 15 (mM)
	
:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	:SERCA
	vmaxsr 	= 0.00027	(mM/ms)
	kpsr 	= 3.75e-6 	(mM) : Fink et al., 2000
	
:------------------------------------------------------------------------------------------------------------------------------------------------------------------
	: IP3R
    jmaxsr 	= 3.5e-6	(mM/ms)
    kip3 	= 0.0008    (mM) 	: Fink et al., 2000
	kactip3 = 0.0003	(mM) 	: Fink et al., 2000  0.3e-3
	konip3 	= 2.7		(/mM-ms): Fink et al., 2000  2.7/mM-ms
	kinhip3 = 0.0002	(mM)    : Fink et al., 2000  0.2 e-3 mM

:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	:CICR(RYR) Parameters
	:Schutter and Smolen, 1998
	kcicr 	= 0.00198	(mM) :Lokuta et al., 2002
	ktcicr 	= 0.0006  	(mM)
	vcicr 	= 5e-7     	(/ms)

	: SER Buffer
	Kmer = 0.5 (mM)
	Bmer = 10 (mM)
	
:----------------
:------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
	:: Mitochondrial Modelling
	
	vmcu = 1.4468e-6	(mM/ms)  
	kmcu = 606e-6 	(mM) 	     : 606 nM
    nmcu = 2.3 (1) 				 : Shutov et al., 2013
	
	:: Mitochondrial NCX
	vncx = 6e-5 	(mM/ms)
	kna = 	8		(mM)	:Boyman et al., 2013  8  	(mM)
	kncx = 	35e-3	(mM)	:Boyman et al., 2013  13e-3	(mM)

:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	: Mito Buffer
	Kmmt = 0.01e-3 (mM) : Faville et al., 2008
	Bmmt = 0.065 (mM)
	
:------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
	:PMCA Parameters
	k1 = 3.74e7        (/mM-s) 
	k2 = .25e6      (/s)
	k3 = .5e3       (/s)
	k4 = 5e0        (/mM-s)
	pump0 = 1.3725e-13 (mol/cm2)  : set to 0 in hoc if this pump not wanted
}

ASSIGNED{
	celsius		(degC)
	ica			(mA/cm2)
	
	cai			(mM)
	cao			(mM)
	caeri		(mM)
	camti 		(mM)
	
	vol[NANN]	(1)
:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	:PMCA
	ica_pmp (mA/cm2)
	last_ica_pmp (mA/cm2)
	parea    (um)
	c1      (1+8 um4/ms)
	c2      (1-10 um/ms)
	c3      (1-10 um/ms)
	c4      (1+8 um4/ms)

:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	:: IP3 Ashhad et al.	
	ip3i		(mM)

:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	: Mitochondrial NCX (MNCX)
	nai			(mM)
:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	: Net mitochondrial flux into the cell	
	jmcu[NANN] 			(mM/ms)
	jmncx[NANN]			(mM /ms)
	jmito[NANN]  		(mM /ms)    
:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	::Fluxes from ER
	jer[NANN]		(mM /ms)
	jip3[NANN]		(mM /ms)
	jserca[NANN]	(mM /ms)
	jcicr[NANN]     (mM /ms) 
	
	: Buffer Binding Ratio of Dye
	: bbrdye[NANN]
	
	: SER and Mito Buffers
	fmer[NANN]
	fmmt[NANN]
}

CONSTANT{
volo = 1e10 (um2) }


STATE{
	:PMCA (Calcium ATPase) Pump on the membrane
	pump            (mol/cm2) <1e-16> 
	pumpca          (mol/cm2) <1e-16>
:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	ca[NANN]			(mM)     <1e-8>
	caer[NANN]			(mM)
	camt[NANN]			(mM)
	caip3ri				(mM)
:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	hc[NANN]        (1)		:IP3 channels in closed state
	ho[NANN]		(1)     :IP3 channels in open state
:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	Ln[NANN]		(mM/ms)	
}

LOCAL factors_done

INITIAL{ LOCAL total
	
	if (factors_done==0) {
		factors_done= 1
		factors()
	}
:-------------------------------------------------------------------------------------------------------------------------------------------------------------------
	: initializing intracellular Ca concentration
	cai = cai0	
:-------------------------------------------------------------------------------------------------------------------------------------------------------------------	   
	:CaATPase Pump
	parms() 
	parea = PI*diam
	pump = pump0
	pumpca = cai*pump*k1/k2
	total = pumpca + pump
	if (total > 1e-9) {
		pump = pump*(pump/total)
		pumpca = pumpca*(pump/total)
	}
	ica_pmp = 0
	last_ica_pmp = 0

:------------------------------------------------------------------------------------------------------------------------------------------------------------------

	FROM i=0 TO NANN-1{
		ca[i]=cai0          :Intracellular      Ca Shell Concentration initialization
		caer[i] = caeri0	:Intracellular ER   Ca Shell Concentration initialization
		camt[i] = camti0	:Intracellular Mito Ca Shell Concentration initialization
		caip3ri = cai0
		
		: ER
		jserca[i] = 0
		jip3[i] = 0
		jcicr[i] = 0
		
		:Mito
		jmcu[i] = 0
		jmncx[i] = 0
		
		:Mitochondirial Buffer parameter
		fmmt[i] = 1/(1+(Kmmt*Bmmt)/(Kmmt+camt[i])^2)
		
		:SER Buffer parameter
		fmer[i] = 1/(1+(Kmer*Bmer)/(Kmer+caer[i])^2)
	}
	caeri= caer[0]
	camti = camt[0]
:-------------------------------------------------------------------------------------------------------------------------------------------------------------------	   
	:Balancing SR Fluxes modified from Fink et al., 2000	
	FROM i=0 TO NANN-1 {
    		 ho[i] = kinhip3/(ca[i]+kinhip3) 	        : Intial open IP3 channels 
    		 hc[i] = 1-ho[i] 	   			            : Intial closed IP3 channels

			   jserca[i] = (-vmaxsr*ca[i]^2 / (ca[i]^2 + kpsr^2))
			   
			   jip3[i] = (jmaxsr*(1-(ca[i]/caer[i])) * ( (ip3i/(ip3i+kip3)) * (ca[i]/(ca[i]+kactip3)) * ho[i] )^3 )
			   
			   if(ca[i] > ktcicr){			
					jcicr[i] = (vcicr* (ca[i]/(kcicr+ca[i])) * (caer[i]-ca[i]) ) 
				} else {
					jcicr[i] = 0
				} 			
		
			   
			   jer[i] = jserca[i]+jip3[i]+jcicr[i]
			   
			   UNITSOFF
			   jmcu[i] = (-vmcu*ca[i]^nmcu / (ca[i]^nmcu + kmcu^nmcu))
			   UNITSON
			   jmncx[i] = vncx*(nai^3/(kna^3 + nai^3))*(camt[i]/(kncx+camt[i]))
			   
			   jmito[i] = jmcu[i]+jmncx[i]
			   
			   Ln[i] = -(jserca[i]+jip3[i]+jcicr[i])/(1 - (ca[i]/caeri0))
    		 }
}

BREAKPOINT{
	SOLVE state METHOD sparse
	 last_ica_pmp = ica_pmp
     ica = ica_pmp
}

LOCAL frat[NANN]

PROCEDURE factors(){
	LOCAL r, dr2
	r = 1/2		        :starts at edge (half diam)
	dr2 = r/(NANN-1)/2	:half thickness of annulus
	vol[0] = 0
	frat[0] = 2*r
	FROM i=0 TO NANN-2{
		vol[i] = vol[i] + PI*(r-dr2/2)*2*dr2 :interior half
		r = r - dr2
		frat[i+1] = 2*PI*r/(2*dr2)		:exterior edge of annulus
                                        :divided by distance between centers
		r = r - dr2
		vol[i+1] = PI*(r+dr2/2)*2*dr2   :outer half of annulus
		}
}

LOCAL dsq, dsqvol,dsqvolmt,dsqvoler

KINETIC state {
	COMPARTMENT ii, (1+bbr)*diam*diam*vol[ii]*0.81 {ca} : cytoplasmic volume is 0.81 of total volume in each shell. IF CHANGED, ALSO CHANGE in ip3_diff.mod
	:COMPARTMENT ii, (1+bbr+bbrdye[ii])*diam*diam*vol[ii]*0.81 {ca} : cytoplasmic volume is 0.81 of total cytoplasmic volume in each shell. IF CHANGED, ALSO CHANGE in ip3_diff.mod file
	COMPARTMENT     (1+bbr)*diam*diam*vol[0]*0.81 {caip3ri}: caip3ri is the ca change by ca release from ip3rs in outermost shell
	: COMPARTMENT     (1+bbr+bbrdye[0])*diam*diam*vol[0] *0.81 {caip3ri}: caip3ri is the ca release from ip3rs in outermost shell
	
	COMPARTMENT jj,	(1/fmer[jj])*diam*diam*vol[jj]*0.12 {caer} :SER volume in the cell : > 10 % of total cytosolic volume Verkhratsky 2002. IF CHANGED, ALSO CHANGE in ip3_diff.mod file
	COMPARTMENT kk, (1/fmmt[kk])*diam*diam*vol[kk]*0.07 {camt} :Mitochondrial volume in the cell : 6.96% of total cytosolic volume-Yilmaz 2017. IF CHANGED, ALSO CHANGE in ip3_diff.mod file
	COMPARTMENT (1e10)*parea {pump pumpca}  :Calcium ATPase Pump on the membrane
	COMPARTMENT volo {cao}

:-------------------------------------------------------------------------------------------------------------------------------------------------------------------     :all currents except pump 
	~ ca[0] << (-(ica - last_ica_pmp)*PI*diam*(1e4)*frat[0]/(2*F))
	
	:PMCA
	:Calcium ATPase Pump on the membrane
	 ~ ca[0] + pump <-> pumpca  (c1,c2)  
	 ~ pumpca <-> pump + cao    (c3,c4)
	  
	ica_pmp = (1e-4) * 2*F*(f_flux - b_flux)/parea 
	  
   
	: Diffusion
	 FROM i=0 TO NANN-2{
		::Diffusion of Cytoplasmic Ca
		 ~ ca[i] <-> ca[i+1]	(DCa*frat[i+1], DCa*frat[i+1])
}
         
	 dsq = diam*diam
     	 
	 FROM i=0 TO NANN-1{
		 dsqvol = dsq*vol[i]*0.81
		 dsqvoler = dsq*vol[i]*0.12
:---------------------------------------------------------------------------------------------------------------------------------------------
		:::SERCA pump, IP3R, CICR(RYR)
	
		:: SERCA pump
		jserca[i] = ((-vmaxsr*ca[i]^2 / (ca[i]^2 + kpsr^2)))
		
		~ ca[i] << (dsqvol*jserca[i])
		~ caer[i] << (-dsqvoler*jserca[i])
		
		:: IP3 channel
		~ hc[i] <-> ho[i]  (konip3*kinhip3, konip3*ca[i])
		jip3[i] = (jmaxsr*(1-(ca[i]/caer[i])) * ( (ip3i/(ip3i+kip3)) * (ca[i]/(ca[i]+kactip3)) * ho[i] )^3 )
        
		~ ca[i] << (dsqvol*jip3[i])
		~ caer[i] << (-dsqvoler*jip3[i])
		
		:Calcium release by ip3r at the periphery. For coupling with calcium activated chloride channel (CACC). See Jin et al., 2013,2015
		if (i==0) {
			~ caip3ri << (dsqvol*jip3[0])
		}
		
		:CICR
		if(ca[i] > ktcicr){			
			jcicr[i] = (vcicr* (ca[i]/(kcicr+ca[i])) * (caer[i]-ca[i]) )
			~ ca[i] << (dsqvol * jcicr[i])
			~ caer[i] << (-dsqvoler * jcicr[i])
		} else {
			jcicr[i] = 0
			~ ca[i] << (dsqvol*jcicr[i])
			~ caer[i] << (-dsqvoler*jcicr[i])
		} 			
		
		:Leak Channels ER
		~ ca[i] << (Ln[i]*(1-ca[i]/caeri0)*dsqvol)
		~ caer[i] << (-Ln[i]*(1-ca[i]/caeri0)*dsqvoler)
		
		jer[i] = jserca[i]+jip3[i]+jcicr[i]
		
		:SER Buffer
		fmer[i] = 1/(1+(Kmer*Bmer)/(Kmer+caer[i])^2) 
:---------------------------------------------------------------------------------------------------------------------------------------------------
		: Mitochondria
		dsqvolmt = dsq*vol[i]*0.07
		
		:::Influx - Jinmt  via mcu
		UNITSOFF
		jmcu[i] = ((-vmcu*ca[i]^nmcu / (ca[i]^nmcu + kmcu^nmcu))*1/(camt[i]*1e3))
		UNITSON
		~ ca[i] << (dsqvol*jmcu[i])
		
		::Outflux - Joutmt  via mncx
		jmncx[i] = vncx*(nai^3/(kna^3 + nai^3))*(camt[i]/(kncx+camt[i]))
		~ ca[i] << (jmncx[i]*dsqvol)
		
		:Total flux
		jmito[i] = jmcu[i]+jmncx[i]
		
		:Mitochondrial Ca change
		~ camt[i] <<  (-(jmncx[i]+jmcu[i])*dsqvolmt) 
		
		:Mito Buffer
		fmmt[i] = 1/(1+(Kmmt*Bmmt)/(Kmmt+camt[i])^2)
		
		:: Fluorescent Dye binding ratio
		:bbrdye[i] = kmdye*Bmdye/(kmdye+ca[i])^2
	}
	
	cai = ca[0]
	caeri= caer[0]
	camti = camt[0]
}

PROCEDURE parms() {
	parea = 2*PI*(diam/2)
        c1 = (1e7)*parea * k1
        c2 = (1e7)*parea * k2
        c3 = (1e7)*parea * k3
        c4 = (1e7)*parea * k4
}


COMMENT
The combination of voltage independent current and calcium
accumulation is more difficult because care must be taken not to count
the pump current twice in the computation of the change ca[0].  Hence
the usage of last_ica_pmp to subtract the pump portion of the total
calcium current in ica so that its effect can be calculated implicitly
via the reaction "pumpca <-> pump + cao".  This artifice makes the
pumping much more stable than the assumption of constant pump current
during the step.  Otherwise, ca[0] is prone to become negative and that
crashes the simulation (especially the automatic computation of eca). 
Calcium currents that are inward are generally safe to compute in
separate models. 

ENDCOMMENT