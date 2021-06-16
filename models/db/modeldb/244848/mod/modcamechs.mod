NEURON {
	SUFFIX modcamechs

	USEION ca READ cao, ica,cai WRITE cai, ica

	RANGE ica_pmp,ca1,ca2,ca3,beta,gamma, DCa
    RANGE DBufm,TBufm,KDm,kfm
    RANGE TBufs,KDs,kfs
    RANGE DCaM,TCaM,KDC,kfC,DCaMCa4,cammax, cammin
    RANGE TCAMKII,kfcamkii,KDcamkii,DCAMKII,DCAMKII_CamCa4,camkiimax, camkiimin

    RANGE Kf_Cam_Camkii,Kd_Cam_Camkii
    RANGE Kf_CamCa1_Camkii,Kd_CamCa1_Camkii
    RANGE Kf_CamCa2_Camkii,Kd_CamCa2_Camkii
    RANGE Kf_CamCa3_Camkii,Kd_CamCa3_Camkii
    RANGE Kf_CamCa4_Camkii,Kd_CamCa4_Camkii
    
    RANGE Kf_CamkiiCam_Ca,Kd_CamkiiCam_Ca
    RANGE Kf_CamkiiCamCa1_Ca,Kd_CamkiiCamCa1_Ca
    RANGE Kf_CamkiiCamCa2_Ca,Kd_CamkiiCamCa2_Ca
    RANGE Kf_CamkiiCamCa3_Ca,Kd_CamkiiCamCa3_Ca

    RANGE DpCAMKII_CamCa4,kfpcamkii,pcamkiimax,pcamkiimin
    RANGE TPP1,kfPP1
    RANGE tmax, tmin, vmax
 	GLOBAL vrat, itocfactor
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
  	PI      = (pi)      (1)
}

PARAMETER {
	cai0 = 50e-6(mM)
	cath = 0.2e-3 (mM) : threshold for ca pump activity
	gamma =  8 (um/s) :  ca pump flux density
	jmax = 3.5e-3 (mM/ms)
    caer = 0.400 (mM)
	Kact = 0.3e-3 (mM)
	kon = 2.7 (/mM-ms)
	Kinh = 0.2e-3 (mM)
    beta  = 1(1)     :introducing beta to take care of other ER mechanisms(SERCA and leak channel density)
    vmax = 1e-4 (mM/ms)
	Kp = 0.27e-3 (mM)
    DCa = 0.22 (um2/ms): 0.22 (um2/ms) Fink et al 2000

	TBufs = 0.45 (mM)
        kfs = 1000 (/mM-ms)
        KDs = 10 (uM)
    TBufm = 0.1 (mM)
	    kfm = 1000 (/mM-ms)
        KDm = 0.2 (uM)
        DBufm = 0.050 (um2/ms)
    TCaM = 1e-3(mM)
        kfC = 8.4848(/mM-ms)
        KDC = 1.0001 (uM)
        DCaM = 4e-3 (um2/ms)
        DCaMCa4 = 4e-3 (um2/ms)
    TCAMKII = 70e-3 (mM)
        Kf_Cam_Camkii = 0.2(/uM-s)
        Kd_Cam_Camkii = 13500 (uM)
        Kf_CamCa1_Camkii =  0.8 (/uM-s)
        Kd_CamCa1_Camkii =  3375.0084 (uM)
        Kf_CamCa2_Camkii = 2 (/uM-s)
        Kd_CamCa2_Camkii = 1125 (uM)
        Kf_CamCa3_Camkii = 100.004 (/uM-s)
        Kd_CamCa3_Camkii = 0.225(uM)
        Kf_CamCa4_Camkii = 100.004 (/uM-s)
        Kd_CamCa4_Camkii = 0.045(uM)

        Kf_CamkiiCam_Ca = 4(/uM-s)
        Kd_CamkiiCam_Ca = 5(uM)
        Kf_CamkiiCamCa1_Ca = 100.004 (/uM-s)
        Kd_CamkiiCamCa1_Ca = 0.2 (uM)
        Kf_CamkiiCamCa2_Ca = 100.004 (/uM-s)
        Kd_CamkiiCamCa2_Ca = 0.02 (uM)
        Kf_CamkiiCamCa3_Ca = 100.004 (/uM-s)
        Kd_CamkiiCamCa3_Ca = 1 (uM)

        DCAMKII = 1.6e-3 (um2/ms)
        DCAMKII_CamCa4 = 1.6e-3 (um2/ms)
        kfpcamkii = 10e-3(/ms)
        DpCAMKII_CamCa4 = 1.6e-3 (um2/ms)
        pcamkii_camca4_0 = 0 (mM)
    TPP1 = 10e-3 (mM)
        kfPP1 = 1.72e-3 (/ms)
        itocfactor=120
}

ASSIGNED {
    	diam      (um)
        ica       (mA/cm2)
        cai       (mM)
        ca1	      (mM)
        ca2       (mM)
        ca3       (mM)
        ica_pmp   (mA/cm2)
        ica_pmp_last   (mA/cm2)
        parea     (um)    :pump area peer unit length
        sump      (mM)
        cao       (mM)
        vrat[Nannuli]  (1)
        L[Nannuli] (mM/ms) : 0.1e-6 mM/ms nominally, but adjusted so that
                           : jchnl + jpump + jleak = 0  when  ca = 0.05 uM and h = Kinh/(ca + Kinh)
        bufs_0 (mM)
        bufm_0 (mM)
        cam_0  (mM)
        camkii_0 (mM)
        cammax (mM)
        cammin (mM)
        camkiimax (mM)
        camkiimin (mM)
        pcamkiimax (mM)
        pcamkiimin (mM)
        tmax  (ms)
        tmin  (ms)
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

        cam[Nannuli]    (mM) <1e-8>
        camca1[Nannuli] (mM) <1e-8>
        camca2[Nannuli] (mM) <1e-8>
        camca3[Nannuli] (mM) <1e-8>
        camca4[Nannuli] (mM) <1e-8>

        camkii[Nannuli] (mM) <1e-8>
        camkii_cam[Nannuli] (mM) <1e-8>
        camkii_camca1[Nannuli] (mM) <1e-8>
        camkii_camca2[Nannuli] (mM) <1e-8>
        camkii_camca3[Nannuli] (mM) <1e-8>
        camkii_camca4[Nannuli] (mM) <1e-8>

        pcamkii_camca4[Nannuli] (mM) <1e-8>
        PP1 (mM) <1e-8>
        dPP1 (mM) <1e-8>
}

LOCAL factors_done, jx
INITIAL {
	
    	if (factors_done==0) {
		factors_done= 1
		factors()
    	}
 
        cai = cai0
	
	bufs_0 = KDs*TBufs/(KDs + (1000)*cai0)

	FROM i=0 TO Nannuli-1 {    
          ca[i] = cai
          bufs[i] = bufs_0
          cabufs[i] = TBufs - bufs_0
          bufm[i] = bufm_0
          cabufm[i] = TBufm - bufm_0
          cam[i] = TCaM
          camca1[i]= 0
          camca2[i]= 0
          camca3[i]= 0
          camca4[i]= 0
          camkii[i] = TCAMKII
          camkii_cam[i]=0
          camkii_camca1[i]=0
          camkii_camca2[i]=0
          camkii_camca3[i]=0
          camkii_camca4[i]=0
          pcamkii_camca4[i]=0
          PP1 = TPP1
          dPP1 = 0
   	}

   	ica=0
   	ica_pmp = 0 
   	ica_pmp_last = 0
    FROM i=0 TO Nannuli-1 {
    		jx = (-vmax*ca[i]^2 / (ca[i]^2 + Kp^2))
            L[i] = -jx/(1 - (ca[i]/caer))
    	}

    sump = cath
    parea = PI*diam
    cammax = 0
    cammin = 1e8
    camkiimax = 0
    camkiimin = 1e8
    pcamkiimax = 0
    pcamkiimin = 1e8
    tmax = t
    tmin = t
    
}

BREAKPOINT {
    SOLVE state METHOD sparse
    

    ica_pmp_last = ica_pmp
    ica = ica_pmp

  FROM i=0 TO Nannuli-1 {
   cam[i] = TCaM - (camca1[i]+camca2[i]+camca3[i]+camca4[i]+camkii_cam[i]+camkii_camca1[i]+camkii_camca2[i]+camkii_camca3[i]+camkii_camca4[i]+pcamkii_camca4[i])
   camkii[i]  = TCAMKII - (camkii_cam[i]+camkii_camca1[i]+camkii_camca2[i]+camkii_camca3[i]+camkii_camca4[i]+pcamkii_camca4[i])
   

  }

    VERBATIM
        if (camca4[0] > cammax) {
            cammax = camca4[0];
            tmax = t;
        }

        if (camca4[0] < cammin) {
            cammin = camca4[0];
            tmin = t;
        }
        
        if (camkii_camca4[0] > camkiimax) {
            camkiimax = camkii_camca4[0];
            tmax = t;
        }

        if (camkii_camca4[0] < camkiimin) {
            camkiimin = camkii_camca4[0];
            tmin = t;
        }
        if (pcamkii_camca4[0] > pcamkiimax) {
            pcamkiimax = pcamkii_camca4[0];
            tmax = t;
        }

        if (pcamkii_camca4[0] < pcamkiimin) {
            pcamkiimin = pcamkii_camca4[0];
            tmin = t;
        }
    ENDVERBATIM

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
  	COMPARTMENT i, diam*diam*vrat[i] {ca sump bufs cabufs cam camca1 camca2 camca3 camca4 bufm cabufm}
    COMPARTMENT i, diam*diam*vrat[i] {camkii camkii_cam camkii_camca1 camkii_camca2 camkii_camca3 camkii_camca4 pcamkii_camca4}
  	COMPARTMENT volo {cao}
  	LONGITUDINAL_DIFFUSION i, DCa*diam*diam*vrat[i] {ca}
    LONGITUDINAL_DIFFUSION i, DCaM*diam*diam*vrat[i] {cam camca1 camca2 camca3 camca4}
    LONGITUDINAL_DIFFUSION i, DCAMKII*diam*diam*vrat[i] {camkii camkii_cam camkii_camca1 camkii_camca2 camkii_camca3 camkii_camca4}
    LONGITUDINAL_DIFFUSION i, DpCAMKII_CamCa4*diam*diam*vrat[i] {pcamkii_camca4}



        :cell membrane ca pump
    ~ ca[0] <-> sump  ((0.001)*parea*gamma*u(ca[0]/(1 (mM)), cath/(1 (mM))), (0.001)*parea*gamma*u(ca[0]/(1 (mM)), cath/(1 (mM))))
  	ica_pmp = 2*FARADAY*(f_flux - b_flux)/parea

  	: all currents except cell membrane ca pump
  	~ ca[0] << (-(ica - ica_pmp_last)*PI*diam/(2*itocfactor*FARADAY))  : ica is Ca efflux
                                                                : converts ica from vgccs to cai


 	 : radial diffusion
   	FROM i=0 TO Nannuli-2 {
   		~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
 	}

   : buffering
   	dsq = diam*diam
   
   	FROM i=0 TO Nannuli-1 {
	 	dsqvol = dsq*vrat[i]
        ~ ca[i] + bufs[i] <-> cabufs[i]  (kfs*dsqvol, (0.001)*KDs*kfs*dsqvol)
	}

    FROM i=0 TO Nannuli-1 {
        dsqvol = dsq*vrat[i]
     :pump
   		~ ca[i] << (-dsqvol*beta*vmax*ca[i]^2 / (ca[i]^2 + Kp^2))
     :leak
   	 	~ ca[i] << (dsqvol*beta*L[i]*(1 - (ca[i]/caer)))
    }

    :Calmodulin

   	FROM i=0 TO Nannuli-1 {
          dsqvol = dsq*vrat[i]
          ~ ca[i] + cam[i] <-> camca1[i]  (kfC*dsqvol, (0.001)*KDC*kfC*dsqvol)
          ~ ca[i] + camca1[i] <-> camca2[i]  (kfC*dsqvol, (0.001)*KDC*kfC*dsqvol)
          ~ ca[i] + camca2[i] <-> camca3[i]  (kfC*dsqvol, (0.001)*KDC*kfC*dsqvol)
          ~ ca[i] + camca3[i] <-> camca4[i]  (kfC*dsqvol, (0.001)*KDC*kfC*dsqvol)

	}

   	FROM i=0 TO Nannuli-2 {
   		~ cam[i] <-> cam[i+1] (DCaM*frat[i+1], DCaM*frat[i+1])
   		~ camca1[i] <-> camca1[i+1] (DCaM*frat[i+1], DCaM*frat[i+1])
   		~ camca2[i] <-> camca2[i+1] (DCaM*frat[i+1], DCaM*frat[i+1])
   		~ camca3[i] <-> camca3[i+1] (DCaM*frat[i+1], DCaM*frat[i+1])
   		~ camca4[i] <-> camca4[i+1] (DCaM*frat[i+1], DCaM*frat[i+1])
 	}

    :CAMKII
   	FROM i=0 TO Nannuli-1 {
          dsqvol = dsq*vrat[i]
          ~ cam[i] + camkii[i] <-> camkii_cam[i] (Kf_Cam_Camkii*dsqvol,(0.001)*Kd_Cam_Camkii*Kf_Cam_Camkii*dsqvol)
          ~ camca1[i] + camkii[i] <-> camkii_camca1[i] (Kf_CamCa1_Camkii*dsqvol, (0.001)*Kd_CamCa1_Camkii*Kf_CamCa1_Camkii*dsqvol)
          ~ camca2[i] + camkii[i] <-> camkii_camca2[i] (Kf_CamCa2_Camkii*dsqvol, (0.001)*Kd_CamCa2_Camkii*Kf_CamCa2_Camkii*dsqvol)
          ~ camca3[i] + camkii[i] <-> camkii_camca3[i] (Kf_CamCa3_Camkii*dsqvol, (0.001)*Kd_CamCa3_Camkii*Kf_CamCa3_Camkii*dsqvol)
          ~ camca4[i] + camkii[i] <-> camkii_camca4[i] (Kf_CamCa4_Camkii*dsqvol, (0.001)*Kd_CamCa4_Camkii*Kf_CamCa4_Camkii*dsqvol)
          
          ~ camkii_cam[i] + ca[i] <-> camkii_camca1[i] (Kf_CamkiiCam_Ca*dsqvol,(0.001)*Kd_CamkiiCam_Ca*Kf_CamkiiCam_Ca*dsqvol)
          ~ camkii_camca1[i] + ca[i] <-> camkii_camca2[i] (Kf_CamkiiCamCa1_Ca*dsqvol,(0.001)*Kd_CamkiiCamCa1_Ca*Kf_CamkiiCamCa1_Ca*dsqvol)
          ~ camkii_camca2[i] + ca[i] <-> camkii_camca3[i] (Kf_CamkiiCamCa2_Ca*dsqvol,(0.001)*Kd_CamkiiCamCa2_Ca*Kf_CamkiiCamCa2_Ca*dsqvol)
          ~ camkii_camca3[i] + ca[i] <-> camkii_camca4[i] (Kf_CamkiiCamCa3_Ca*dsqvol,(0.001)*Kd_CamkiiCamCa3_Ca*Kf_CamkiiCamCa3_Ca*dsqvol)
	}

   	FROM i=0 TO Nannuli-2 {
   	  ~ camkii[i] <-> camkii[i+1] (DCAMKII*frat[i+1], DCAMKII*frat[i+1])
      ~ camkii_cam[i] <-> camkii_cam[i+1] (DCAMKII*frat[i+1], DCAMKII*frat[i+1])
   	  ~ camkii_camca1[i] <-> camkii_camca1[i+1] (DCAMKII*frat[i+1], DCAMKII*frat[i+1])
   	  ~ camkii_camca2[i] <-> camkii_camca2[i+1] (DCAMKII*frat[i+1], DCAMKII*frat[i+1])
   	  ~ camkii_camca3[i] <-> camkii_camca3[i+1] (DCAMKII*frat[i+1], DCAMKII*frat[i+1])
   	  ~ camkii_camca4[i] <-> camkii_camca4[i+1] (DCAMKII*frat[i+1], DCAMKII*frat[i+1])
 	}

    : pCAMKII
   	FROM i=0 TO Nannuli-1 {
         dsqvol = dsq*vrat[i]
         ~ camkii_camca4[i] <-> pcamkii_camca4[i]  (kfpcamkii*dsqvol,0)
	}

   	FROM i=0 TO Nannuli-2 {
   	  ~ pcamkii_camca4[i] <-> pcamkii_camca4[i+1] (DpCAMKII_CamCa4*frat[i+1], DpCAMKII_CamCa4*frat[i+1])
 	}

    : PP1(dephosphorylation)
   	FROM i=0 TO Nannuli-1 {
          dsqvol = dsq*vrat[i]
          ~ pcamkii_camca4[i] + PP1 <-> camkii_camca4[i] + dPP1  (kfPP1*dsqvol,0)
    }

    ~ dPP1 <-> PP1 (1e8,0)

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
