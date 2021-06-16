: Calcium ion accumulation with radial and longitudinal diffusion + Pump + Buffers
NEURON {
	SUFFIX cadifusnpumpOGBenddif
	USEION ca READ cai,cao,ica WRITE cai,ica
	RANGE ica_pmp, TotalPump 
	RANGE CaBuf01,CaBufav1,caiav 
	RANGE CaBuf02,CaBufav2	
	RANGE Bufferav1,Bufferav2 
	RANGE k1,k2,k3,k4 
	RANGE ica 
	GLOBAL vrat
	RANGE TotalBuffer1,TotalBuffer2,k1buf1,k2buf1,k1buf2,k2buf2   
	RANGE cai0,cai1,cai2,cai3
	GLOBAL DCa,Dbuf1,Dcabuf1 
	GLOBAL Dbuf2,Dcabuf2 
	
	}
DEFINE Nannuli 4 : must be >=2 (i.e. at least shell and core)
UNITS {
	(molar)=(1/liter)
	(mM)=(millimolar)
	(um)=(micron)
	(mA)=(milliamp)
	FARADAY=(faraday) (10000 coulomb)
	PI=(pi) (1)
	(mol) = (1)
	}

PARAMETER {
	DCa = 0.6 (um2/ms)

	Dbuf1=0.015 (um2/ms)
	Dcabuf1=0.015 (um2/ms) 

	Dbuf2=0.015 (um2/ms)
	Dcabuf2=0.015 (um2/ms) 
	
	k1buf1 = 100 (/mM-ms) 
	k2buf1 = 0.1 (/ms)
	TotalBuffer1= 0.003 (mM)

	k1buf2 = 100 (/mM-ms) 
	k2buf2 = 0.1 (/ms)
	TotalBuffer2= 0.003 (mM)

	k1 = 1	   (/mM-ms)
	k2 = 0.005 (/ms)
	k3 = 1	   (/ms)
	k4 = 0.005 (/mM-ms)
	
	TotalPump = 1e-11 (mol/cm2)  
	}

ASSIGNED {
	diam (um)
	ica (mA/cm2)
	cai (mM)

	cai0 (mM)
	cai1 (mM)
	cai2 (mM)
	cai3 (mM)

	CaBuf01 (mM) 
	CaBuf02 (mM)
	caiav (mM) 
	CaBufav1 (mM)
	CaBufav2 (mM) 
	Bufferav1 (mM)
	Bufferav2 (mM)

	vrat[Nannuli] (1) 

	Kd1 (/mM)
	Kd2 (/mM)
	B01 (mM)
	B02 (mM) 

	cao	(mM)
	ica_pmp	(mA/cm2)
	parea	(um) 

	}

CONSTANT { volo = 1e10  (um2) }

STATE {

	ca[Nannuli] (mM) <1e-10> 
	CaBuffer1[Nannuli] (mM)
	Buffer1[Nannuli] (mM)
	CaBuffer2[Nannuli] (mM)
	Buffer2[Nannuli] (mM)

	pump	(mol/cm2)
	pumpca	(mol/cm2)
	}

BREAKPOINT { SOLVE state METHOD sparse
	     ica = ica_pmp }

LOCAL factors_done
INITIAL {
	if (factors_done == 0) { : flag becomes 1 in the first segment
	factors_done = 1 : all subsequent segments will have
	factors() : vrat = 0 unless vrat is GLOBAL
				}
	Kd1 = k1buf1/k2buf1
	Kd2 = k1buf2/k2buf2 
	B01 = TotalBuffer1/(1 + Kd1*cai)
	B02 = TotalBuffer2/(1 + Kd2*cai) 
	FROM i=0 TO Nannuli-1 {
			ca[i] = cai
			Buffer1[i] = B01
			CaBuffer1[i] = TotalBuffer1 - B01
			Buffer2[i] = B02
			CaBuffer2[i] = TotalBuffer2 - B02
			}
	
	parea = PI*diam
	pump = TotalPump/(1 + (cai*k1/k2))
	pumpca = TotalPump - pump :
	}

LOCAL frat[Nannuli] : scales the rate constants for model geometry

PROCEDURE factors() {
		LOCAL r, dr2
		r = 1/2 : starts at edge (half diam)
		dr2 = r/(Nannuli-1)/2 : full thickness of outermost annulus
					: half thickness of all other annuli
		vrat[0] = 0
		frat[0] = 2*r
		FROM i=0 TO Nannuli-2 {
				vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2 : interior half
				r = r - dr2
				frat[i+1] = 2*PI*r/(2*dr2) : outer radius of annulus
				: div by distance between centers
				r = r - dr2
				vrat[i+1] = PI*(r+dr2/2)*2*dr2 : outer half of annulus
					}
			}

LOCAL dsq, dsqvol 

KINETIC state {
		COMPARTMENT i, diam*diam*vrat[i] {ca CaBuffer1 Buffer1}
		COMPARTMENT i, diam*diam*vrat[i] {ca CaBuffer2 Buffer2}

		COMPARTMENT (1e10)*parea {pump pumpca}
		COMPARTMENT volo {cao}
		
		LONGITUDINAL_DIFFUSION i, DCa*diam*diam*vrat[i] {ca}

		LONGITUDINAL_DIFFUSION i, Dcabuf1*diam*diam*vrat[i] {CaBuffer1}
		LONGITUDINAL_DIFFUSION i, Dbuf1*diam*diam*vrat[i] {Buffer1}
		
		LONGITUDINAL_DIFFUSION i, Dcabuf2*diam*diam*vrat[i] {CaBuffer2}
		LONGITUDINAL_DIFFUSION i, Dbuf2*diam*diam*vrat[i] {Buffer2}
	
		~ ca[0] + pump <-> pumpca (k1*parea*(1e10), k2*parea*(1e10))
		~ pumpca <-> pump + cao
		(k3*parea*(1e10), k4*parea*(1e10))
		CONSERVE pump + pumpca = TotalPump * parea * (1e10)
		ica_pmp = 2*FARADAY*(f_flux - b_flux)/parea
	
		~ ca[0] << (-(ica - ica_pmp)*PI*diam/(2*FARADAY)) 
		FROM i=0 TO Nannuli-2 {
				~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
			
				~ Buffer1[i] <-> Buffer1[i+1] (Dbuf1*frat[i+1], Dbuf1*frat[i+1])
				~ CaBuffer1[i] <-> CaBuffer1[i+1] (Dcabuf1*frat[i+1], Dcabuf1*frat[i+1])
			
				~ Buffer2[i] <-> Buffer2[i+1] (Dbuf2*frat[i+1], Dbuf2*frat[i+1])
				~ CaBuffer2[i] <-> CaBuffer2[i+1] (Dcabuf2*frat[i+1], Dcabuf2*frat[i+1])
				}
		dsq = diam*diam
		FROM i=0 TO Nannuli-1 {
				dsqvol = dsq*vrat[i]
				~ ca[i] + Buffer1[i] <-> CaBuffer1[i] (k1buf1*dsqvol, k2buf1*dsqvol)
				~ ca[i] + Buffer2[i] <-> CaBuffer2[i] (k1buf2*dsqvol, k2buf2*dsqvol) :this migth cause error
				
				}
			
		cai = ca[0]
		cai0 = ca[0]
		cai1 = ca[1]
		cai2 = ca[2]
		cai3 = ca[3]

		caiav = 0.25*(ca[0]+ca[1]+ca[2]+ca[3])

		CaBuf01 = CaBuffer1[0] 
		CaBufav1= 0.25*(CaBuffer1[0]+CaBuffer1[1]+CaBuffer1[2]+CaBuffer1[3]) 
		Bufferav1=0.25*(Buffer1[0]+Buffer1[1]+Buffer1[2]+Buffer1[3]) 

		CaBuf02 = CaBuffer2[0] 
		CaBufav2= 0.25*(CaBuffer2[0]+CaBuffer2[1]+CaBuffer2[2]+CaBuffer2[3]) 
		Bufferav2=0.25*(Buffer2[0]+Buffer2[1]+Buffer2[2]+Buffer2[3]) 
		}
