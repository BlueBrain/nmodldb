:Interneuron Cells to Pyramidal Cells GABA with local Ca2+ pool and read public soma Ca2+ pool

NEURON {
	POINT_PROCESS interV2pyrDNE_STFD
	USEION ca READ eca,ica
	NONSPECIFIC_CURRENT igaba
	RANGE initW
	RANGE Cdur_gaba, AlphaTmax_gaba, Beta_gaba, Erev_gaba, gbar_gaba, W, on_gaba, g_gaba
	RANGE eca, tauCa, Icatotal
	RANGE ICag, P0g, fCag
	RANGE Cainf, pooldiam, z
	RANGE lambda1, lambda2, threshold1, threshold2
	RANGE fmax, fmin, Wmax, Wmin, maxChange, normW, scaleW, srcid, destid
	RANGE pregid,postgid, thr_rp
	RANGE F, f, tauF, D1, d1, tauD1, D2, d2, tauD2
	RANGE facfactor
}
 
UNITS {
	(mV) = (millivolt)
        (nA) = (nanoamp)
	(uS) = (microsiemens)
	FARADAY = 96485 (coul)
	pi = 3.141592 (1)
}

PARAMETER {

	srcid = -1 (1)
	destid = -1 (1)
	
	Cdur_gaba = 0.7254 (ms)
	AlphaTmax_gaba = 7.2609 (/ms)
	Beta_gaba = 0.2667 (/ms)
	Erev_gaba = -75 (mV)
	gbar_gaba = 0.6e-3 (uS)

	Cainf = 50e-6 (mM)
	pooldiam =  1.8172 (micrometer)
	z = 2

	k = 0.01	
	
	tauCa = 50 (ms)
	
	P0g = .01
	fCag = .024
	
	lambda1 = 4 : 3 : 2.0 : 2.0
	lambda2 = .01
	threshold1 = 0.47 :  0.48 : 0.45 : 0.4 : 0.95 : 1.35 :0.75 :0.55 (uM)
	threshold2 = 0.52 :  0.53 : 0.5 : 0.45 : 1.0 : 1.4 : 0.8 : 0.65 :0.70 (uM)

	:GABA Weight
	initW = 4.5 :  :  3 :  2.5 : 3 : 5 : 6.25 : 5
	fmax = 4.2
	fmin = .8
	
	GAPstart1 = 96000
	GAPstop1 = 196000
	
	thr_rp = 1 : .7
	
	facfactor = 1
	: the (1) is needed for the range limits to be effective
        f = 0 (1) < 0, 1e9 > : 1.3 (1) < 0, 1e9 >    : facilitation
        tauF = 20 (ms) < 1e-9, 1e9 >
        d1 = 0.95 (1) < 0, 1 >     : fast depression
        tauD1 = 40 (ms) < 1e-9, 1e9 >
        d2 = 0.9 (1) < 0, 1 >     : slow depression
        tauD2 = 70 (ms) < 1e-9, 1e9 >	
	
	NEstart1 = 39500
	NEstop1 = 40000	
	NEstart2 = 35900
	NEstop2 = 36000	

	NE_t1 = 1 : 0.95
	NE_t2 = 1 :0.7 : 0.8
	NE_S = 1 : 0.4	
	Beta1 = 0.001  (/ms) : 1/decay time for neuromodulators
	Beta2 = 0.0001  (/ms)			
}

ASSIGNED {
	v (mV)
	eca (mV)
	ica (nA)
	
	igaba (nA)
	g_gaba (uS)
	on_gaba
	W

	t0 (ms)

	ICan (mA)
	ICag (mA)
	Afactor	(mM/ms/nA)
	Icatotal (mA)

	dW_gaba
	Wmax
	Wmin
	maxChange
	normW
	scaleW
	
	pregid
	postgid

	rp
	tsyn
	
	fa
	F
	D1
	D2	
}

STATE { r_nmda r_gaba capoolcon }

INITIAL {

	on_gaba = 0
	r_gaba = 0
	W = initW

	t0 = -1

	Wmax = fmax*initW
	Wmin = fmin*initW
	maxChange = (Wmax-Wmin)/10
	dW_gaba = 0

	capoolcon = Cainf
	Afactor	= 1/(z*FARADAY*4/3*pi*(pooldiam/2)^3)*(1e6)

	fa =0
	F = 1
	D1 = 1
	D2 = 1	
}

BREAKPOINT {
	SOLVE release METHOD cnexp
}

DERIVATIVE release {
	if (t0>0) {
		if (rp < thr_rp) {
			if (t-t0 < Cdur_gaba) {
				on_gaba = 1
			} else {
				on_gaba = 0
			}
		} else {
			on_gaba = 0
		}
	}
	if (t0>0) {
		if (rp < thr_rp) {
			if (t-t0 < Cdur_gaba) {
				on_gaba = 1
			} else {
				on_gaba = 0
			}
		} else {
			on_gaba = 0
		}
	}

	r_gaba' = AlphaTmax_gaba*on_gaba*(1-r_gaba)-Beta_gaba*r_gaba

	dW_gaba = eta(capoolcon)*(lambda1*omega(capoolcon, threshold1, threshold2)-lambda2*GAP1(GAPstart1, GAPstop1)*W)*dt

	: Limit for extreme large weight changes
	if (fabs(dW_gaba) > maxChange) {
		if (dW_gaba < 0) {
			dW_gaba = -1*maxChange
		} else {
			dW_gaba = maxChange
		}
	}

	:Normalize the weight change
	normW = (W-Wmin)/(Wmax-Wmin)
	if (dW_gaba < 0) {
		scaleW = sqrt(fabs(normW))
	} else {
		scaleW = sqrt(fabs(1.0-normW))
	}

	W = W + dW_gaba*scaleW
	
	:Weight value limits
	if (W > Wmax) { 
		W = Wmax
	} else if (W < Wmin) {
 		W = Wmin
	}


	g_gaba = gbar_gaba*r_gaba*facfactor*NEn(NEstart1,NEstop1)*NE2(NEstart2,NEstop2)   : Norepinephrine effect on GABA		    
	igaba = W*g_gaba*(v - Erev_gaba)

	ICag = P0g*g_gaba*(v - eca)	
	Icatotal = ICag + k*ica*4*pi*((15/2)^2)*(0.01)    :  icag+k*ica*Area of soma*unit change
	capoolcon'= -fCag*Afactor*Icatotal + (Cainf-capoolcon)/tauCa
}

NET_RECEIVE(dummy_weight) {
	t0 = t
	rp = unirand()	
	
	:F  = 1 + (F-1)* exp(-(t - tsyn)/tauF)
	D1 = 1 - (1-D1)*exp(-(t - tsyn)/tauD1)
	D2 = 1 - (1-D2)*exp(-(t - tsyn)/tauD2)
 :printf("%g\t%g\t%g\t%g\t%g\t%g\n", t, t-tsyn, F, D1, D2, facfactor)
	:printf("%g\t%g\t%g\t%g\n", F, D1, D2, facfactor)
	tsyn = t
	
	facfactor = F * D1 * D2

	::F = F+f  :F * f
	
	if (F > 3) { 
	F=3	}	
	if (facfactor < 0.7) { 
	facfactor=0.7
	}
	D1 = D1 * d1
	D2 = D2 * d2
:printf("\t%g\t%g\t%g\n", F, D1, D2)
}

:::::::::::: FUNCTIONs and PROCEDUREs ::::::::::::

FUNCTION eta(Cani (mM)) {
	LOCAL taulearn, P1, P2, P4, Cacon
	P1 = 0.1
	P2 = P1*1e-4
	P4 = 1
	Cacon = Cani*1e3
	taulearn = P1/(P2+Cacon*Cacon*Cacon)+P4
	eta = 1/taulearn*0.001
}

FUNCTION omega(Cani (mM), threshold1 (uM), threshold2 (uM)) {
	LOCAL r, mid, Cacon
	Cacon = Cani*1e3
	r = (threshold2-threshold1)/2
	mid = (threshold1+threshold2)/2
	if (Cacon <= threshold1) { omega = 0}
	else if (Cacon >= threshold2) {	omega = 1/(1+50*exp(-50*(Cacon-threshold2)))}
	else {omega = -sqrt(r*r-(Cacon-mid)*(Cacon-mid))}
}
FUNCTION NEn(NEstart1 (ms), NEstop1 (ms)) {
	LOCAL NEtemp1, NEtemp2, NEtemp3, NEtemp4, NEtemp5, NEtemp6, NEtemp7, NEtemp8, NEtemp9, NEtemp10, NEtemp11, NEtemp12, NEtemp13, NEtemp14, NEtemp15, NEtemp16, NEtemp17, NEtemp18,s
	NEtemp1 = NEstart1+4000
	NEtemp2 = NEtemp1+4000
	NEtemp3 = NEtemp2+4000
	NEtemp4 = NEtemp3+4000
	NEtemp5 = NEtemp4+4000
	NEtemp6 = NEtemp5+4000
	NEtemp7 = NEtemp6+4000
	NEtemp8 = NEtemp7+4000
	NEtemp9 = NEtemp8+4000
	NEtemp10 = NEtemp9+4000
	NEtemp11 = NEtemp10+4000
	NEtemp12 = NEtemp11+4000
	NEtemp13 = NEtemp12+4000
	NEtemp14 = NEtemp13+4000
	NEtemp15 = NEtemp14 + 4000 + 100000     : 100sec Gap
	NEtemp16 = NEtemp15 + 4000 
	NEtemp17 = NEtemp16 + 4000
	NEtemp18 = NEtemp17 + 4000

	if (t <= NEstart1) { NEn = 1.0}
	else if (t >= NEstart1 && t <= NEstop1) {NEn = NE_t1}					: 2nd tone in early conditioning (EC)
		else if (t > NEstop1 && t < NEtemp1) {NEn = 1.0 + (NE_t1-1)*exp(-Beta1*(t-NEstop1))}  		: Basal level
	else if (t >= NEtemp1 && t <= NEtemp1+500) {NEn = NE_t1}					: 3rd tone EC
		else if (t > NEtemp1+500 && t < NEtemp2) {NEn = 1.0 + (NE_t1-1)*exp(-Beta1*(t-NEstop1))}  	: Basal level
	else if (t >= NEtemp2 && t <= NEtemp2+500) {NEn = NE_t1}					: 4th tone EC
		else if (t > NEtemp2+500 && t < NEtemp3) {NEn = 1.0 + (NE_t1-1)*exp(-Beta1*(t-(NEtemp2+500)))}  	: Basal level	
	else if (t >= NEtemp3 && t <= NEtemp3+500) {NEn = NE_t1}					: 5th tone EC
		else if (t > NEtemp3+500 && t < NEtemp4) {NEn = 1.0 + (NE_t1-1)*exp(-Beta1*(t-(NEtemp3+500)))}  	: Basal level
	else if (t >= NEtemp4 && t <= NEtemp4+500) {NEn = NE_t1}					: 6th tone EC
		else if (t > NEtemp4+500 && t < NEtemp5) {NEn = 1.0 + (NE_t1-1)*exp(-Beta1*(t-(NEtemp4+500)))}  		: Basal level
	else if (t >= NEtemp5 && t <= NEtemp5+500) {NEn = NE_t1}					: 7th tone EC
		else if (t > NEtemp5+500 && t < NEtemp6) {NEn = 1.0 + (NE_t1-1)*exp(-Beta1*(t-(NEtemp5+500)))}  	: Basal level
	else if (t >= NEtemp6 && t <= NEtemp6+500) {NEn = NE_t1}					: 8th tone EC
		else if (t > NEtemp6+500 && t < NEtemp7) {NEn = 1.0 + (NE_t1-1)*exp(-Beta1*(t-(NEtemp6+500)))}  	: Basal level
	
	else if (t >= NEtemp7 && t <= NEtemp7+500) {NEn = NE_t2}					: 9th tone	- Second Step late cond (LC)
		else if (t > NEtemp7+500 && t < NEtemp8) {NEn = 1.0 + (NE_t2-1)*exp(-Beta2*(t-(NEtemp7+500)))}  		: Basal level
	else if (t >= NEtemp8 && t <= NEtemp8+500) {NEn = NE_t2}					: 10th tone  LC
		else if (t > NEtemp8+500 && t < NEtemp9) {NEn = 1.0 + (NE_t2-1)*exp(-Beta2*(t-(NEtemp8+500)))}		: Basal level	
	else if (t >= NEtemp9 && t <= NEtemp9+500) {NEn = NE_t2}					: 11th tone  LC 
		else if (t > NEtemp9+500 && t < NEtemp10) {NEn = 1.0 + (NE_t2-1)*exp(-Beta2*(t-(NEtemp9+500)))}  	: Basal level	
	else if (t >= NEtemp10 && t <= NEtemp10+500) {NEn = NE_t2}					: 12th tone  LC
		else if (t > NEtemp10+500 && t < NEtemp11) {NEn = 1.0 + (NE_t2-1)*exp(-Beta2*(t-(NEtemp10+500)))}  	: Basal level
	else if (t >= NEtemp11 && t <= NEtemp11+500) {NEn = NE_t2}					: 13th tone  LC
		else if (t > NEtemp11+500 && t < NEtemp12) {NEn = 1.0 + (NE_t2-1)*exp(-Beta2*(t-(NEtemp11+500)))}  	: Basal level
	else if (t >= NEtemp12 && t <= NEtemp12+500) {NEn = NE_t2}					: 14th tone  LC
		else if (t > NEtemp12+500 && t < NEtemp13) {NEn = 1.0 + (NE_t2-1)*exp(-Beta2*(t-(NEtemp12+500)))}  	: Basal level
	else if (t >= NEtemp13 && t <= NEtemp13+500) {NEn = NE_t2}					: 15th tone  LC
		else if (t > NEtemp13+500 && t < NEtemp14) {NEn = 1.0 + (NE_t2-1)*exp(-Beta2*(t-(NEtemp13+500)))}  	: Basal level
	else if (t >= NEtemp14 && t <= NEtemp14+500) {NEn = NE_t2}					: 16th tone  LC
		else if (t > NEtemp14+500 && t < NEtemp15) {NEn = 1.0 + (NE_t2-1)*exp(-Beta2*(t-(NEtemp14+500)))}  	: Basal level
	
	else if (t >= NEtemp15 && t <= NEtemp15+500) {NEn = NE_t2}					: 1st tone in Early Extinction EE
		else if (t > NEtemp15+500 && t < NEtemp16) {NEn = 1.0 + (NE_t2-1)*exp(-Beta2*(t-(NEtemp15+500)))} 		: Basal level
	else if (t >= NEtemp16 && t <= NEtemp16+500) {NEn = NE_t2}					: 2nd tone EE
		else if (t > NEtemp16+500 && t < NEtemp17) {NEn = 1.0 + (NE_t2-1)*exp(-Beta2*(t-(NEtemp16+500)))}  	: Basal level
	else if (t >= NEtemp17 && t <= NEtemp17+500) {NEn = NE_t2}					: 3rd tone EE
		else if (t > NEtemp17+500 && t < NEtemp18) {NEn = 1.0 + (NE_t2-1)*exp(-Beta2*(t-(NEtemp17+500)))}  	: Basal level	
	else if (t >= NEtemp18 && t <= NEtemp18+500) {NEn = NE_t2}					: 4th tone EE	
		else  {	NEn = 1.0}
}
FUNCTION NE2(NEstart2 (ms), NEstop2 (ms)) {
	LOCAL NE2temp1, NE2temp2, NE2temp3, NE2temp4, NE2temp5, NE2temp6, NE2temp7, NE2temp8, NE2temp9, NE2temp10, NE2temp11, NE2temp12, NE2temp13, NE2temp14, NE2temp15, NE2temp16,s
	NE2temp1 = NEstart2 + 4000
	NE2temp2 = NE2temp1 + 4000
	NE2temp3 = NE2temp2 + 4000
	NE2temp4 = NE2temp3 + 4000
	NE2temp5 = NE2temp4 + 4000
	NE2temp6 = NE2temp5 + 4000
	NE2temp7 = NE2temp6 + 4000
	NE2temp8 = NE2temp7 + 4000
	NE2temp9 = NE2temp8 + 4000
	NE2temp10 = NE2temp9 + 4000
	NE2temp11 = NE2temp10 + 4000
	NE2temp12 = NE2temp11 + 4000 
	NE2temp13 = NE2temp12 + 4000
	NE2temp14 = NE2temp13 + 4000
	NE2temp15 = NE2temp14 + 4000
	
	if (t <= NEstart2) { NE2 = 1.0}
	else if (t >= NEstart2 && t <= NEstop2) {NE2 = NE_S }					: 1st shock
		else if (t > NEstop2 && t < NE2temp1) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NEstop2+500)))} 
	else if (t >= NE2temp1 && t <= NE2temp1+100) {NE2=NE_S}					: 2nd shock
		else if (t > NE2temp1+100 && t < NE2temp2) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp1+100)))}   				 
	else if (t >= NE2temp2 && t <= NE2temp2+100) {NE2=NE_S}					: 3rd shock
		else if (t > NE2temp2+100 && t < NE2temp3) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp2+100)))}  				 
	else if (t >= NE2temp3 && t <= NE2temp3+100) {NE2=NE_S}					: 4th shock
		else if (t > NE2temp3+100 && t < NE2temp4) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp3+100)))}  				 
	else if (t >= NE2temp4 && t <= NE2temp4+100) {NE2=NE_S}					: 5th shock
		else if (t > NE2temp4+100 && t < NE2temp5) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp4+100)))}  				 
	else if (t >= NE2temp5 && t <= NE2temp5+100) {NE2=NE_S}					: 6th shock
		else if (t > NE2temp5+100 && t < NE2temp6) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp5+100)))} 				 
	else if (t >= NE2temp6 && t <= NE2temp6+100) {NE2=NE_S}					: 7th shock
		else if (t > NE2temp6+100 && t < NE2temp7) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp6+100)))}  				 
	else if (t >= NE2temp7 && t <= NE2temp7+100) {NE2=NE_S}					: 8th shock
		else if (t > NE2temp7+100 && t < NE2temp8) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp7+100)))}  				    
	else if (t >= NE2temp8 && t <= NE2temp8+100) {NE2=NE_S }					: 9th shock
		else if (t > NE2temp8+100 && t < NE2temp9) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp8+100)))}  				    
	else if (t >= NE2temp9 && t <= NE2temp9+100) {NE2=NE_S }					: 10th shock
		else if (t > NE2temp9+100 && t < NE2temp10) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp9+100)))}  				    
	else if (t >= NE2temp10 && t <= NE2temp10+100) {NE2=NE_S}					: 11th shock
		else if (t > NE2temp10+100 && t < NE2temp11) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp10+100)))}  				 
	else if (t >= NE2temp11 && t <= NE2temp11+100) {NE2=NE_S }					: 12th shock
		else if (t > NE2temp11+100 && t < NE2temp12) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp11+100)))}  				 
	else if (t >= NE2temp12 && t <= NE2temp12+100) {NE2=NE_S}					: 13th shock
		else if (t > NE2temp12+100 && t < NE2temp13) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp12+100)))} 				 
	else if (t >= NE2temp13 && t <= NE2temp13+100) {NE2=NE_S }					: 14th shock
		else if (t > NE2temp13+100 && t < NE2temp14) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp13+100)))}   				 
	else if (t >= NE2temp14 && t <= NE2temp14+100) {NE2=NE_S}					: 15th shock
		else if (t > NE2temp14+100 && t < NE2temp15) {NE2 = 1.0 + (NE_S-1)*exp(-Beta2*(t-(NE2temp14+100)))}  				 
	else if (t >= NE2temp15 && t <= NE2temp15+100) {NE2=NE_S}					: 16th shock
		else  {	NE2 = 1.0}
}

FUNCTION GAP1(GAPstart1 (ms), GAPstop1 (ms)) {
	LOCAL s
	if (t <= GAPstart1) { GAP1 = 1}
	else if (t >= GAPstart1 && t <= GAPstop1) {GAP1 = 1}					: During the Gap, apply lamda2*2
	else  {	GAP1 = 1}
}
FUNCTION unirand() {    : uniform random numbers between 0 and 1
        unirand = scop_random()
}
