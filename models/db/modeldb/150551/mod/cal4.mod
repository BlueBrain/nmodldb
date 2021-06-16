:Modified from NEURON implementation of Fink et al., 2000


NEURON {
	SUFFIX cal4
	USEION ca READ cao, ica WRITE cai, ica
	USEION ip3 READ ip3i VALENCE 1
	RANGE ica_pmp,ca1, ca2,alpha,beta,ca3,gamma,ip3ca, DCa,jip3
 	GLOBAL vrat, TBufs, KDs, TBufs, TBufm, KDm
}

DEFINE Nannuli 4

UNITS {
	(mol)	= (1)
 	(molar) = (1/liter)
  	(uM)    = (micromolar)
  	(mM)    = (millimolar)
  	(um)    = (micron)
  	(mA)    = (milliamp)
  	FARADAY = (faraday)  (10000 coulomb)
  	PI      = (pi)       (1)
}

PARAMETER {
	cai0 = 50e-6(mM)
	cath = 0.2e-3 (mM) : threshold for ca pump activity
	gamma = 8 (um/s) : ca pump flux density
	jmax = 3.5e-3 (mM/ms)
	caer = 0.400 (mM)
	Kip3 = 0.8e-3 (mM)
	Kact = 0.3e-3 (mM)
	kon = 2.7 (/mM-ms)
	Kinh = 0.2e-3 (mM)
	alpha = 1 (1) : relative abundance of ER mechanisms : alpha only specific for ip3 receptor,
	beta  = 1(1)           :introducing beta to take care of other ER mechanisms(SERCA and leak channel density)

	vmax =1e-4   (mM/ms) :1e-4 revised
	Kp = 0.27e-3 (mM)
	DCa = 0.22 (um2/ms) :Fink et al 2000
	TBufs = 0.45 (mM)
        kfs = 1000 (/mM-ms) : try these for now
        KDs = 10 (uM)
	TBufm = 0.075 (mM)
	kfm = 1000 (/mM-ms) : try these for now
        KDm = 0.24 (uM)
        DBufm = 0.050 (um2/ms)


}

ASSIGNED {
	diam      (um)
  	ica       (mA/cm2)
        cai       (mM)
	jip3	  (mM/ms)
        ca1	  (mM)
        ca2       (mM)
        ca3       (mM)
        ica_pmp   (mA/cm2)
	ica_pmp_last   (mA/cm2)
        parea     (um)    :pump area peer unit length
        sump      (mM)
        cao       (mM)
        ip3i      (mM)
        jchnl    (mM/ms)
        vrat[Nannuli]  (1)
	L[Nannuli] (mM/ms)  adjusted so that
                         : jchnl + jpump + jleak = 0  when  ca = 0.05 uM and h = Kinh/(ca + Kinh)
        bufs_0 (mM)
	bufm_0 (mM)
	
	ip3ca	(mM)
} 


CONSTANT { volo = 1e10 (um2) }

STATE {
     	ca[Nannuli]     (mM) <1e-7>
     	hc[Nannuli]    
     	ho[Nannuli]
     	bufs[Nannuli]    (mM) <1e-3>
     	cabufs[Nannuli]  (mM) <1e-7>
	bufm[Nannuli]    (mM) <1e-4>
        cabufm[Nannuli]  (mM) <1e-8>
	ip3cas [Nannuli] (mM)

}



BREAKPOINT {
     	SOLVE state METHOD sparse
     	ica_pmp_last = ica_pmp
     	ica = ica_pmp

}
LOCAL factors_done, jx
INITIAL {
	
    	if (factors_done==0) {
		factors_done= 1
		factors()
    	}
 
        cai = cai0
	jip3=0
	bufs_0 = KDs*TBufs/(KDs + (1000)*cai0)
	bufm_0 = KDm*TBufm/(KDm + (1000)*cai0)

	FROM i=0 TO Nannuli-1 {    
     		ca[i] = cai
		bufs[i] = bufs_0
       		cabufs[i] = TBufs - bufs_0
		bufm[i] = bufm_0
    		cabufm[i] = TBufm - bufm_0

   	}
	
   	ica=0
   	ica_pmp = 0 
   	ica_pmp_last = 0


	FROM i=0 TO Nannuli-1 {
    		ho[i] = Kinh/(ca[i]+Kinh)
    		hc[i] = 1-ho[i]
    		jx = (-vmax*ca[i]^2 / (ca[i]^2 + Kp^2))
    		jx = jx + jmax*(1-(ca[i]/caer)) * ( (ip3i/(ip3i+Kip3)) * (ca[i]/(ca[i]+Kact)) * ho[i] )^3
     	   	L[i] = -jx/(1 - (ca[i]/caer))
    	}

    	sump = cath
    	parea = PI*diam   
}

LOCAL frat[Nannuli]

PROCEDURE factors() {
	LOCAL r, dr2
  	r = 1/2                : starts at edge (half diam)
  	dr2 = r/(Nannuli-1)/2  : full thickness of outermost annulus,
                               : half thickness of all other annuli
  	vrat[0] = 0
  	frat[0] = 2*r

  	FROM i=0 TO Nannuli-2 {
    		vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2  : interior half
   		 r = r - dr2
   		 frat[i+1] = 2*PI*r/(2*dr2)  : outer radius of annulus
                                             : div by distance between centers
   		 r = r - dr2
    		vrat[i+1] = PI*(r+dr2/2)*2*dr2  : outer half of annulus
  	}
}


LOCAL dsq, dsqvol

KINETIC state {
  	COMPARTMENT i, diam*diam*vrat[i] {ca  bufs cabufs bufm cabufm sump}
  	COMPARTMENT volo {cao}
  	LONGITUDINAL_DIFFUSION i, DCa*diam*diam*vrat[i] {ca}
  	LONGITUDINAL_DIFFUSION i, DBufm*diam*diam*vrat[i] {bufm cabufm}




        :cell membrane ca pump
  	~ ca[0] <-> sump  ((0.001)*parea*gamma*u(ca[0]/(1 (mM)), cath/(1 (mM))), (0.001)*parea*gamma*u(ca[0]/(1 (mM)), cath/(1 (mM))))
  	ica_pmp = 2*FARADAY*(f_flux - b_flux)/parea

  	: all currents except cell membrane ca pump
  	~ ca[0] << (-(ica - ica_pmp_last)*PI*diam/(2*FARADAY))  : ica is Ca efflux

 	 : radial diffusion
   	FROM i=0 TO Nannuli-2 {
   		~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
 	}

	: buffering
   	dsq = diam*diam
   
   	FROM i=0 TO Nannuli-1 {
	 	dsqvol = dsq*vrat[i]
     	 	~ ca[i] + bufs[i] <-> cabufs[i]  (kfs*dsqvol, (0.001)*KDs*kfs*dsqvol)
		:~ ca[i] + bufm[i] <-> cabufm[i]  (kfm*dsqvol, (0.001)*KDm*kfm*dsqvol) :to simulate high affinity dyes, used only for the simplified 3 cylinder model in the paper

    	}


       	:SERCA pump, channel
  	FROM i=0 TO Nannuli-1 {
    		dsqvol = dsq*vrat[i]

   	 	: pump
   	 	~ ca[i] << (-dsqvol*beta*vmax*ca[i]^2 / (ca[i]^2 + Kp^2))

    		: channel
   	 	~ hc[i] <-> ho[i]  (kon*Kinh, kon*ca[i])
   	 	~ ca[i] << ( dsqvol*alpha*jmax*(1-(ca[i]/caer)) * ( (ip3i/(ip3i+Kip3)) * (ca[i]/(ca[i]+Kact)) * ho[i] )^3 )
 	 	: leak
   	 	~ ca[i] << (dsqvol*beta*L[i]*(1 - (ca[i]/caer)))
		~ ip3cas[i] << (dsqvol*alpha*jmax*(1-(ca[i]/caer)) * ( (ip3i/(ip3i+Kip3)) * (ca[i]/(ca[i]+Kact)) * ho[i] )^3 )
  	}
	
	jip3 = ( dsqvol*alpha*jmax*(1-(ca[0]/caer)) * ( (ip3i/(ip3i+Kip3)) * (ca[0]/(ca[0]+Kact)) * ho[0] )^3 )

:	ip3ca=0
:	FROM i=0 TO Nannuli-1 {
:		ip3ca=ip3ca+ip3cas[i]
:	}

	ip3ca=ip3cas[0]

  	cai = ca[0]
  	ca1 = ca[1]
  	ca2 = ca[2]
  	ca3 = ca[3]
}


FUNCTION u (x, th) {
  	if (x>th) {
    		u = 1
  	} else {
    		u = 0
  	}
}
