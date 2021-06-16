TITLE Cav2.3 voltage-gated calcium channel with kinetic scheme and two Zn2+ binding sites

COMMENT

Markov model for Cav2.3 channel function and Zn2+-induced modulation based on ionic and gating currents recorded from HEK-293 cells stably transfected with human Cav2.3+ß3-subunits.
Transition rates for movement of sensor 1 in channels with Zn2+ bound to site 1 are shifted by the voltage-offset voff1 and slowed by the slowing factor aoff1.
Channels with Zn2+ bound to site 2 are non-conductive (blocked) and their transition rates for opening and closing are slowed by the factor f. 

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)} 


NEURON {
	SUFFIX CaR
	USEION ca READ cai, cao WRITE ica
	RANGE p, n, igate, iion, zno
	GLOBAL vmin, vmax
}

UNITS {
	F    = (faraday) (coulomb)
	R    = (k-mole) (joule/degC)
	e    = 0.00000000000000000016 (coulomb)
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
}

PARAMETER {
	v 			(mV)
	
	n	= 1.6e010	(/cm2)	: number of channels per cm2
	p  	= 7.5e-014	(cm3/s)	: single channel permeability
	vmin	= -200		(mV)
	vmax	=  200		(mV)
	
	zno	= 0		(mM)	  : free Zn2+ concentration in mM

	

			
}

ASSIGNED {
	ica  	(mA/cm2)	: sum of ionic (iion) and gating (igate) current
	igate	(mA/cm2)	: gating current
	iion	(mA/cm2)	: ionic current
	celsius	(degC)		
	cao	(mM)		
	cai  	(mM)		
	fw1 	(/ms)		
	bw1 	(/ms)
	fw2 	(/ms)
	bw2 	(/ms)
	fw3	(/ms)
	bw3	(/ms)
	fw4	(/ms)
	bw4	(/ms)
	ko	(/ms)
	kc	(/ms)
	fw1b 	(/ms)
	bw1b 	(/ms)
	kob	(/ms)
	kcb	(/ms)	
}

STATE {
	C0
	C1
	C3
	C4	
	C12
	C13
	C14
	C34
	C123
	C124
	C134
	C1234
	O12 
	O123 
	O124 
	O1234 
	IC34 
	IC123 
	IC124 
	IC134 
	IC1234 
	IO123 
	IO124 
	IO1234 
	ICS123 
	ICS124 
	ICS134 
	ICS1234 
	IOS123 
	IOS124 
	IOS1234 
	CB0 
	CB1 
	CB3 
	CB4 
	CB12 
	CB13 
	CB14 
	CB34 
	CB123 
	CB124 
	CB134 
	CB1234 
	OB12 
	OB123 
	OB124 
	OB1234 
	ICB34 
	ICB123 
	ICB124 
	ICB134 
	ICB1234 
	IOB123 
	IOB124 
	IOB1234 
	ICSB123 
	ICSB124 
	ICSB134 
	ICSB1234 
	IOSB123 
	IOSB124 
	IOSB1234 
	OBB12 
	OBB123 
	OBB124 
	OBB1234 
	OBBB12 
	OBBB123 
	OBBB124 
	OBBB1234 
	CBB12 
	CBB123 
	CBB124 
	CBB1234 
	CBBB12 
	CBBB123 
	CBBB124 
	CBBB1234

}


INITIAL {
	SOLVE kin STEADYSTATE sparse

}

BREAKPOINT {
	SOLVE kin METHOD sparse
	iion = (O12 + O123 + O124 + O1234 + OB12 + OB123 + OB124 + OB1234) * n * p * ghk(v, cai, cao)
	igate = (1e+006)*n*e*(1.5716*(fw1*(C0+C3+C4+C34+IC34)+fw1b*(CB0+CB3+CB4+CB34+ICB34)-bw1*(C1+C13+C14+C134+IC134)-(bw1b)*(CB1+CB13+CB14+CB134+ICB134))
+0.19249*(fw2*(C1+C13+C14+C134+IC134+ICS134+CB1+CB13+CB14+CB134+ICB134+ICSB134)-bw2*(C12+C123+C124+C1234+IC1234+ICS1234+CB12+CB123+CB124+CB1234+ICB1234+ICSB1234))
+0.90262*(fw3*(C0+C1+C4+C12+C14+C124+O12+O124+IC124+IO124+ICS124+IOS124+CB0+CB1+CB4+CB12+CB14+CB124+OB12+OB124+ICB124+IOB124+ICSB124+IOSB124+CBB12+CBB124+CBBB12+CBBB124+OBB12+OBB124+OBBB12+OBBB124)-bw3*(C3+C13+C34+C123+C134+C1234+O123+O1234+IC1234+IO1234+ICS1234+IOS1234+CB3+CB13+CB34+CB123+CB134+CB1234+OB123+OB1234+ICB1234+IOB1234+ICSB1234+IOSB1234+CBB123+CBB1234+CBBB123+CBBB1234+OBB123+OBB1234+OBBB123+OBBB1234))
+1.5327*(fw4*(C0+C1+C3+C12+C13+C123+O12+O123+IC123+IO123+ICS123+IOS123+CB0+CB1+CB3+CB12+CB13+CB123+OB12+OB123+ICB123+IOB123+ICSB123+IOSB123+CBB12+CBB123+CBBB12+CBBB123+OBB12+OBB123+OBBB12+OBBB123)-bw4*(C4+C14+C34+C124+C134+C1234+O124+O1234+IC1234+IO1234+ICS1234+IOS1234+CB4+CB14+CB34+CB124+CB134+CB1234+OB124+OB1234+ICB1234+IOB1234+ICSB1234+IOSB1234+CBB124+CBB1234+CBBB124+CBBB1234+OBB124+OBB1234+OBBB124+OBBB1234))
+1.5541*(ko*(C12+C123+C124+C1234+IC123+IC124+IC1234+ICS123+ICS124+ICS1234+CB12+CB123+CB124+CB1234+ICB123+ICB124+ICB1234+ICSB123+ICSB124+ICSB1234)+(kob)*(CBB12+CBB123+CBB124+CBB1234+CBBB12+CBBB123+CBBB124+CBBB1234)-kc*(O12+O123+O124+O1234+IO123+IO124+IO1234+IOS123+IOS124+IOS1234+OB12+OB123+OB124+OB1234+IOB123+IOB124+IOB1234+IOSB123+IOSB124+IOSB1234)-(kcb)*(OBB12+OBB123+OBB124+OBB1234+OBBB12+OBBB123+OBBB124+OBBB1234)))
	ica = iion+igate
}

KINETIC kin {
	rates(v)

	~ C0		<-> 	C1		(fw1, bw1)
	~ C3		<-> 	C13		(fw1, bw1)
	~ C4		<-> 	C14		(fw1, bw1)
	~ C34		<->	C134		(fw1, bw1)

	~ C1		<-> 	C12		(fw2, bw2)
	~ C13		<->	C123		(fw2, bw2)
	~ C14		<->	C124		(fw2, bw2)
	~ C134		<->	C1234		(fw2, bw2)

	~ C0		<-> 	C3		(fw3, bw3)
	~ C1		<-> 	C13		(fw3, bw3)
	~ C4		<-> 	C34		(fw3, bw3)
	~ C12		<->	C123		(fw3, bw3)
	~ C14		<->	C134		(fw3, bw3)
	~ C124		<->	C1234		(fw3, bw3)
	~ O12		<->	O123		(fw3, bw3)
	~ O124		<->	O1234		(fw3, bw3)

	~ C0		<->	C4		(fw4, bw4)
	~ C1		<-> 	C14		(fw4, bw4)
	~ C3		<-> 	C34		(fw4, bw4)
	~ C12		<->	C124		(fw4, bw4)
	~ C13		<->	C134		(fw4, bw4)
	~ C123		<->	C1234		(fw4, bw4)
	~ O12		<->	O124		(fw4, bw4)
	~ O123		<->	O1234		(fw4, bw4)


	~ C12		<->	O12		(ko, kc)
	~ C123		<->	O123		(ko, kc)
	~ C124		<->	O124		(ko, kc)
	~ C1234		<->	O1234		(ko, kc)


	~ C34		<->	IC34		(0.016, 0.0027)
	~ C123		<->	IC123		(0.016, 0.0027)
	~ C124		<->	IC124		(0.016, 0.0027)
	~ C134		<->	IC134		(0.016, 0.0027)
	~ C1234		<->	IC1234		(0.016, 0.0027)
	~ O123		<->	IO123		(0.016, 0.0027)
	~ O124		<->	IO124		(0.016, 0.0027)
	~ O1234		<->	IO1234		(0.016, 0.0027)

	~ IC123		<->	IO123		(ko, kc)
	~ IC124		<->	IO124		(ko, kc)
	~ IC1234	<->	IO1234		(ko, kc)

	~ IC34		<->	IC134		(fw1, bw1)
	~ IC134		<->	IC1234		(fw2, bw2)
	~ IC124		<->	IC1234		(fw3, bw3)
	~ IO124		<->	IO1234		(fw3, bw3)
	~ IC123		<->	IC1234		(fw4, bw4)
	~ IO123		<->	IO1234		(fw4, bw4)


	~ IC123		<->	ICS123		(0.008, 0.00064)
	~ IC124		<->	ICS124		(0.008, 0.00064)
	~ IC134		<->	ICS134		(0.008, 0.00064)
	~ IC1234	<->	ICS1234		(0.008, 0.00064)
	~ IO123		<->	IOS123		(0.008, 0.00064)
	~ IO124		<->	IOS124		(0.008, 0.00064)
	~ IO1234	<->	IOS1234		(0.008, 0.00064)


	~ ICS123	<->	IOS123		(ko, kc)
	~ ICS124	<->	IOS124		(ko, kc)
	~ ICS1234	<->	IOS1234		(ko, kc)

	~ ICS134	<->	ICS1234		(fw2, bw2)
	~ ICS124	<->	ICS1234		(fw3, bw3)
	~ IOS124	<->	IOS1234		(fw3, bw3)
	~ ICS123	<->	ICS1234		(fw4, bw4)
	~ IOS123	<->	IOS1234		(fw4, bw4)


	~ C0		<->	CB0		(zno*100, 0.003*100)
	~ C3		<->	CB3		(zno*100, 0.003*100)		
	~ C4		<->	CB4		(zno*100, 0.003*100)	
	~ C34		<->	CB34		(zno*100, 0.003*100)	
	~ IC34		<->	ICB34		(zno*100, 0.003*100)	
	
	~ C1		<->	CB1		(zno*100, 0.003*100)
	~ C12		<->	CB12		(zno*100, 0.003*100)
	~ C13		<->	CB13		(zno*100, 0.003*100)
	~ C14		<->	CB14		(zno*100, 0.003*100)
	~ C123		<->	CB123		(zno*100, 0.003*100)
	~ C124		<->	CB124		(zno*100, 0.003*100)
	~ C134		<->	CB134		(zno*100, 0.003*100)	
	~ C1234		<->	CB1234		(zno*100, 0.003*100)
	~ O12		<->	OB12		(zno*100, 0.003*100)
	~ O123		<->	OB123		(zno*100, 0.003*100)
	~ O124		<->	OB124		(zno*100, 0.003*100)
	~ O1234		<->	OB1234		(zno*100, 0.003*100)
	~ IC123		<->	ICB123		(zno*100, 0.003*100)	
	~ IC124		<->	ICB124		(zno*100, 0.003*100)	
	~ IC134		<->	ICB134		(zno*100, 0.003*100)	
	~ IC1234	<->	ICB1234		(zno*100, 0.003*100)	
	~ IO123		<->	IOB123		(zno*100, 0.003*100)
	~ IO124		<->	IOB124		(zno*100, 0.003*100)
	~ IO1234	<->	IOB1234		(zno*100, 0.003*100)
	~ ICS123	<->	ICSB123		(zno*100, 0.003*100)
	~ ICS124	<->	ICSB124		(zno*100, 0.003*100)
	~ ICS134	<->	ICSB134		(zno*100, 0.003*100)
	~ ICS1234	<->	ICSB1234	(zno*100, 0.003*100)	
	~ IOS123	<->	IOSB123		(zno*100, 0.003*100)
	~ IOS124	<->	IOSB124		(zno*100, 0.003*100)
	~ IOS1234	<->	IOSB1234	(zno*100, 0.003*100)

	~ CBB12		<->	CBBB12		(zno*100, 0.003*100)
	~ CBB123	<->	CBBB123		(zno*100, 0.003*100)
	~ CBB124	<->	CBBB124		(zno*100, 0.003*100)
	~ CBB1234	<->	CBBB1234	(zno*100, 0.003*100)

	~ OBB12		<->	OBBB12		(zno*100, 0.003*100)
	~ OBB123	<->	OBBB123		(zno*100, 0.003*100)
	~ OBB124	<->	OBBB124		(zno*100, 0.003*100)
	~ OBB1234	<->	OBBB1234	(zno*100, 0.003*100)	


	~ CB0		<-> 	CB1		(fw1b, bw1b)
	~ CB3		<-> 	CB13		(fw1b, bw1b)
	~ CB4		<-> 	CB14		(fw1b, bw1b)
	~ CB34		<->	CB134		(fw1b, bw1b)

	~ CB1		<-> 	CB12		(fw2, bw2)
	~ CB13		<->	CB123		(fw2, bw2)
	~ CB14		<->	CB124		(fw2, bw2)
	~ CB134		<->	CB1234		(fw2, bw2)

	~ CB0		<-> 	CB3		(fw3, bw3)
	~ CB1		<-> 	CB13		(fw3, bw3)
	~ CB4		<-> 	CB34		(fw3, bw3)
	~ CB12		<->	CB123		(fw3, bw3)
	~ CB14		<->	CB134		(fw3, bw3)
	~ CB124		<->	CB1234		(fw3, bw3)
	~ OB12		<->	OB123		(fw3, bw3)
	~ OB124		<->	OB1234		(fw3, bw3)

	~ CB0		<->	CB4		(fw4, bw4)
	~ CB1		<-> 	CB14		(fw4, bw4)
	~ CB3		<-> 	CB34		(fw4, bw4)
	~ CB12		<->	CB124		(fw4, bw4)
	~ CB13		<->	CB134		(fw4, bw4)
	~ CB123		<->	CB1234		(fw4, bw4)
	~ OB12		<->	OB124		(fw4, bw4)
	~ OB123		<->	OB1234		(fw4, bw4)

	~ CB12		<->	OB12		(ko, kc)
	~ CB123		<->	OB123		(ko, kc)
	~ CB124		<->	OB124		(ko, kc)
	~ CB1234	<->	OB1234		(ko, kc)


	~ CB34		<->	ICB34		(0.016, 0.0027)
	~ CB123		<->	ICB123		(0.016, 0.0027)
	~ CB124		<->	ICB124		(0.016, 0.0027)
	~ CB134		<->	ICB134		(0.016, 0.0027)
	~ CB1234	<->	ICB1234		(0.016, 0.0027)
	~ OB123		<->	IOB123		(0.016, 0.0027)
	~ OB124		<->	IOB124		(0.016, 0.0027)
	~ OB1234	<->	IOB1234		(0.016, 0.0027)

	~ ICB123	<->	IOB123		(ko, kc)
	~ ICB124	<->	IOB124		(ko, kc)
	~ ICB1234	<->	IOB1234		(ko, kc)

	~ ICB34		<->	ICB134		(fw1b, bw1b)
	~ ICB134	<->	ICB1234		(fw2, bw2)
	~ ICB124	<->	ICB1234		(fw3, bw3)
	~ IOB124	<->	IOB1234		(fw3, bw3)
	~ ICB123	<->	ICB1234		(fw4, bw4)
	~ IOB123	<->	IOB1234		(fw4, bw4)


	~ ICB123	<->	ICSB123		(0.008, 0.00064)
	~ ICB124	<->	ICSB124		(0.008, 0.00064)
	~ ICB134	<->	ICSB134		(0.008, 0.00064)
	~ ICB1234	<->	ICSB1234	(0.008, 0.00064)
	~ IOB123	<->	IOSB123		(0.008, 0.00064)
	~ IOB124	<->	IOSB124		(0.008, 0.00064)
	~ IOB1234	<->	IOSB1234	(0.008, 0.00064)


	~ ICSB123	<->	IOSB123		(ko, kc)
	~ ICSB124	<->	IOSB124		(ko, kc)
	~ ICSB1234	<->	IOSB1234	(ko, kc)

	~ ICSB134	<->	ICSB1234	(fw2, bw2)
	~ ICSB124	<->	ICSB1234	(fw3, bw3)
	~ IOSB124	<->	IOSB1234	(fw3, bw3)
	~ ICSB123	<->	ICSB1234	(fw4, bw4)
	~ IOSB123	<->	IOSB1234	(fw4, bw4)






	~ C12		<->	CBB12		(zno*100, 0.1*100)
	~ C123		<->	CBB123		(zno*100, 0.1*100)
	~ C124		<->	CBB124		(zno*100, 0.1*100)
	~ C1234		<->	CBB1234		(zno*100, 0.1*100)

	~ O12		<->	OBB12		(zno*100, 0.1*100)
	~ O123		<->	OBB123		(zno*100, 0.1*100)
	~ O124		<->	OBB124		(zno*100, 0.1*100)
	~ O1234		<->	OBB1234		(zno*100, 0.1*100)


	~ CBB12		<->	CBB123		(fw3, bw3)
	~ CBB124	<->	CBB1234		(fw3, bw3)
	~ OBB12		<->	OBB123		(fw3, bw3)
	~ OBB124	<->	OBB1234		(fw3, bw3)

	~ CBB12		<->	CBB124		(fw4, bw4)
	~ CBB123	<->	CBB1234		(fw4, bw4)
	~ OBB12		<->	OBB124		(fw4, bw4)
	~ OBB123	<->	OBB1234		(fw4, bw4)

	~ CBB12		<->	OBB12		(kob, kcb)
	~ CBB123	<->	OBB123		(kob, kcb)
	~ CBB124	<->	OBB124		(kob, kcb)
	~ CBB1234	<->	OBB1234		(kob, kcb)



	~ CB12		<->	CBBB12		(zno*100, 0.1*100)
	~ CB123		<->	CBBB123		(zno*100, 0.1*100)
	~ CB124		<->	CBBB124		(zno*100, 0.1*100)
	~ CB1234	<->	CBBB1234	(zno*100, 0.1*100)

	~ OB12		<->	OBBB12		(zno*100, 0.1*100)
	~ OB123		<->	OBBB123		(zno*100, 0.1*100)
	~ OB124		<->	OBBB124		(zno*100, 0.1*100)
	~ OB1234	<->	OBBB1234	(zno*100, 0.1*100)


	~ CBBB12	<->	CBBB123		(fw3, bw3)
	~ CBBB124	<->	CBBB1234	(fw3, bw3)
	~ OBBB12	<->	OBBB123		(fw3, bw3)
	~ OBBB124	<->	OBBB1234	(fw3, bw3)


	~ CBBB12	<->	CBBB124		(fw4, bw4)
	~ CBBB123	<->	CBBB1234	(fw4, bw4)
	~ OBBB12	<->	OBBB124		(fw4, bw4)
	~ OBBB123	<->	OBBB1234	(fw4, bw4)

	~ CBBB12	<->	OBBB12		(kob, kcb)
	~ CBBB123	<->	OBBB123		(kob, kcb)
	~ CBBB124	<->	OBBB124		(kob, kcb)
	~ CBBB1234	<->	OBBB1234	(kob, kcb)

	CONSERVE C0+C1+C3+C4+C12+C13+C14+C34+C123+C124+C134+C1234+O12+O123+O124+O1234+IC34+IC123+IC124+IC134+IC1234+IO123+IO124+IO1234+ICS123+ICS124+ICS134+ICS1234+IOS123+IOS124+IOS1234+CB0+CB1+CB3+CB4+CB12+CB13+CB14+CB34+CB123+CB124+CB134+CB1234+OB12+OB123+OB124+OB1234+ICB34+ICB123+ICB124+ICB134+ICB1234+IOB123+IOB124+IOB1234+ICSB123+ICSB124+ICSB134+ICSB1234+IOSB123+IOSB124+IOSB1234+OBB12+OBB123+OBB124+OBB1234+OBBB12+OBBB123+OBBB124+OBBB1234+CBB12+CBB123+CBB124+CBB1234+CBBB12+CBBB123+CBBB124+CBBB1234=1

}


FUNCTION rates1(v, keq, z, x, veq) {
	rates1 = keq * exp(z*x*(v-veq)/25)
}

FUNCTION rates2(v, keq, z, x, veq) {
	rates2 = keq * exp(-z*(1-x)*(v-veq)/25)
}

FUNCTION rates3(v, keq, aoff, z, x, veq, voff) {
	rates3 = keq * aoff * exp(z*x*(v-veq-voff)/25)
}

FUNCTION rates4(v, keq, aoff, z, x, veq, voff) {
	rates4 = keq * aoff * exp(-z*(1-x)*(v-veq-voff)/25)
}


PROCEDURE rates(v(mV)) {
UNITSOFF
	fw1 	=(rates1(v, 1.4, 1.5716, 0.001002, 24.354))
	bw1 	=(rates2(v, 1.4, 1.5716, 0.001002, 24.354))
	fw2 	=(rates1(v, 1.7, 0.19249, 0.53623, -96.152))
	bw2 	=(rates2(v, 1.7, 0.19249, 0.53623, -96.152))
	fw3 	=(rates1(v, 1.5135, 0.90262, 0.57543, -28.922))
	bw3 	=(rates2(v, 1.5135, 0.90262, 0.57543, -28.922))
	fw4 	=(rates1(v, 0.02, 1.5327, 0.33817, -30.987))
	bw4 	=(rates2(v, 0.02, 1.5327, 0.33817, -30.987))
	ko 	=(rates1(v, 2.1461, 1.5541, 0.61063, -5.8045))
	kc 	=(rates2(v, 2.1461, 1.5541, 0.61063, -5.8045))
	fw1b 	=(rates3(v, 1.4, 0.05, 1.5716, 0.001002, 24.354, 50))
	bw1b 	=(rates4(v, 1.4, 0.05, 1.5716, 0.001002, 24.354, 50))
	kob 	=(rates3(v, 2.1461, 0.3, 1.5541, 0.61063, -5.8045, 0))
	kcb 	=(rates4(v, 2.1461, 0.3, 1.5541, 0.61063, -5.8045, 0))
UNITSON
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (0.001 coul/cm3) {
	LOCAL z

	z = (0.001)*2*F*v/(R*(celsius+273.15))
	ghk = (0.001)*2*F*(ci*efun(-z) - co*efun(z))
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}