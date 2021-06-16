TITLE SBML model: glia_2013 

? Source: Jolivet ... Magistretti, 2015 (original NGV model)
? Source: Xu...Gomez, 2011 (glycogen metabolism expantion module); 
? Boras ... McCulloch, 2014 (cAMP-dependent kinase rate constants)



UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  
  (mM) 	= (millimole)
}


NEURON {

  SUFFIX glia_2013
  RANGE unit_compartment
  RANGE V
  RANGE Fin
  RANGE Jstimg
  RANGE Jstimn
  RANGE ADPn
  RANGE ADPg
  RANGE O2meanc
  RANGE Vv0
  RANGE JleakNan
  RANGE JleakNag
  RANGE Jpumpn
  RANGE Jpumpg
  RANGE JGLCce
  RANGE JGLCcg
  RANGE JGLCeg
  RANGE JGLCen
  RANGE JLACec
  RANGE JLACcg
  RANGE JLACeg
  RANGE JLACen
  RANGE JGLCec
  RANGE JGLCgc
  RANGE JGLCge
  RANGE JGLCne
  RANGE JLACce
  RANGE JLACgc
  RANGE JLACge
  RANGE JLACne
  RANGE JHKPFKn
  RANGE JHKPFKg
  RANGE JPGKn
  RANGE JPGKg
  RANGE JPKn
  RANGE JPKg
  RANGE JLDHn
  RANGE JLDHg
  RANGE Jmitoing
  RANGE Jmitoinn
  RANGE Jmitooutg
  RANGE Jmitooutn
  RANGE Jshuttleg
  RANGE Jshuttlen
  RANGE JCKg
  RANGE JCKn
  RANGE JO2mcg
  RANGE JO2mcn
  RANGE JO2c
  RANGE JGLCc
  RANGE JLACc
  RANGE IL
  RANGE INa
  RANGE IK
  RANGE ICa
  RANGE ImAHP
  RANGE Ipump
  RANGE Fout
  RANGE dAMPdATPn
  RANGE dAMPdATPg
  RANGE Ve
  RANGE Vcap
  RANGE Vg
  RANGE Vn
  RANGE zeta
  RANGE qAK
  RANGE A
  RANGE un
  RANGE ug
  RANGE ren
  RANGE reg
  RANGE rce
  RANGE rcn
  RANGE rcg
  RANGE Isyn
  RANGE kCKplusn
  RANGE kCKplusg
  RANGE Nan0
  RANGE F0
  RANGE SmVn
  RANGE SmVg
  RANGE R
  RANGE F
  RANGE RToverF
  RANGE psig
  RANGE Nae
  RANGE KtGLCen
  RANGE KtGLCeg
  RANGE KtGLCcg
  RANGE KtGLCce
  RANGE KtLACen
  RANGE KtLACeg
  RANGE KtLACge
  RANGE KtLACgc
  RANGE KtLACcg
  RANGE KtLACec
  RANGE KIATP
  RANGE nH1
  RANGE Kg
  RANGE KO2
  RANGE HbOP
  RANGE nh
  RANGE KO2mito
  RANGE Cm
  RANGE gL
  RANGE gNa
  RANGE gK
  RANGE gCa
  RANGE gmAHP
  RANGE KD
  RANGE tauCa
  RANGE Ca0
  RANGE EK
  RANGE ECa
  RANGE phih
  RANGE phin
  RANGE tauv
  RANGE alphav
  RANGE O2a
  RANGE GLCa
  RANGE gNan
  RANGE gNag
  RANGE gKpas
  RANGE kpumpn
  RANGE kpumpg
  RANGE Kmpump
  RANGE C
  RANGE N
  RANGE Kmmito
  RANGE kLDHnplus
  RANGE kLDHnminus
  RANGE kLDHgplus
  RANGE kLDHgminus
  RANGE Mcyton
  RANGE Mcytog
  RANGE Mmiton
  RANGE Mmitog
  RANGE KmADPn
  RANGE KmADPg
  RANGE KmNADn
  RANGE KmNADg
  RANGE KmNADHn
  RANGE KmNADHg
  RANGE kCKn
  RANGE kCKg
  RANGE mNADg
  RANGE TmaxGLCen
  RANGE TmaxGLCce
  RANGE TmaxGLCeg
  RANGE TmaxGLCcg
  RANGE TmaxLACgc
  RANGE TmaxLACcg
  RANGE TmaxLACen
  RANGE TmaxLACne
  RANGE TmaxLACge
  RANGE TmaxLACeg
  RANGE TmaxLACec
  RANGE kHKPFKn
  RANGE kHKPFKg
  RANGE PScapoverVn
  RANGE PScapoverVg
  RANGE Vmaxoutn
  RANGE Vmaxoutg
  RANGE Vmaxinn
  RANGE Vmaxing
  RANGE kPGKn
  RANGE kPGKg
  RANGE kPKn
  RANGE kPKg
  RANGE JATPasesn
  RANGE JATPasesg
  RANGE kCKminusn
  RANGE kCKminusg
  RANGE TNADHn
  RANGE TNADHg
  RANGE LACa
  RANGE Rplusg
  RANGE Rminusg
  RANGE Rplusn
  RANGE Rminusn
  RANGE hh_alphan
  RANGE hh_alpham
  RANGE hh_betam
  RANGE hh_alphah
  RANGE hh_betah
  RANGE hh_betan
  RANGE hh_minfinity
  RANGE hh_ninfinity
  RANGE hh_hinfinity
  RANGE hh_taun
  RANGE hh_mCa
  RANGE hh_tauh
  RANGE hh_EL
  RANGE Finprime
  RANGE JGLUeg
  RANGE vPumpg0 
  RANGE CC	
  RANGE kd
  RANGE ka
  RANGE vL1 	
  RANGE vL2
  RANGE v_L1 		
  RANGE v_L2
  RANGE vtL1 			
  RANGE vtF1			
  RANGE ktF1					
  RANGE vtS1			
  RANGE ktS1					
  RANGE vd_Bgluc		
  RANGE kd_Bgluc		
  RANGE vglucgn		
  RANGE vGgluc		
  RANGE vd_Bglucgn	
  RANGE vins			
  RANGE vd_Bins		
  RANGE B_vip
  RANGE B_ne			
  RANGE kc1		
  RANGE kc2		
  RANGE kcm1		
  RANGE kcm2		
  RANGE kgc1		
  RANGE kgc2		
  RANGE k_gc1	
  RANGE k_gc2	
  RANGE k11
  RANGE k22
  RANGE ni		
  RANGE ng	
  RANGE nn	
  RANGE st		
  RANGE k_a 		
  RANGE kDins		
  RANGE kDins2
  RANGE kDglucgn	
  RANGE kDne	
  RANGE kfeed
  RANGE en1
  RANGE en2
  RANGE en5
  RANGE en6
  RANGE ep1
  RANGE ep2
  RANGE ep3
  RANGE ep5
  RANGE ep8
  RANGE ep9
  RANGE ep12
  RANGE ep13
  RANGE g6p		
  RANGE kt	
  RANGE ki13
  RANGE kL1		
  RANGE kL2 
  RANGE kL3
  RANGE kL4
  RANGE k_L1		
  RANGE k_L2		
  RANGE kmL1		
  RANGE km_L1	
  RANGE kmL2		
  RANGE km_L2
  RANGE kmL3
  RANGE kmL4
  RANGE ktL1		
  RANGE kgi		
  RANGE kg		
  RANGE kg2		
  RANGE kg3		
  RANGE kg4		
  RANGE kg5		
  RANGE kg6		
  RANGE kg7		
  RANGE kg8		
  RANGE kmg3		
  RANGE kmg4		
  RANGE kmg7		
  RANGE kmg8		
  RANGE kins		
  RANGE k1ins		
  RANGE kmIns	
  RANGE kglucgn	
  RANGE kmGlgn	
  RANGE k1glucgn
  RANGE kd_Bglucgn
  RANGE kd_Bins	
  RANGE ka
  RANGE PP1_GPa
  RANGE v_L2
  RANGE glyc 
  RANGE kmg6 
  RANGE kmg5
  RANGE kDvip
  RANGE kDvip2
  RANGE kd_mg
  RANGE kmind
  RANGE kmaxd
  RANGE gluco
  RANGE pt 
  RANGE vigluc 
  RANGE vins
  RANGE vd_Bglucgn 
  RANGE vGgluc
  RANGE vglucgn
  RANGE vd_Bgluc
  RANGE vtS1
  RANGE vtF1
  RANGE v_L2
  RANGE ne
  RANGE vip
  RANGE vip2
  RANGE viptotal
  RANGE ne0
  RANGE vip0
  RANGE cAMP_ne
  RANGE cAMP_vip
  RANGE cAMP_vip2
  RANGE cAMPtotal
  RANGE cAMPtotsyn
  RANGE tauvip
  RANGE taune  
  RANGE sigma
  RANGE tauvip2
  RANGE tauvip3
  RANGE taune2 
  RANGE vL3
  RANGE v_L3
  RANGE vL4
  RANGE vL13
  RANGE pep	
  RANGE oa_c
  RANGE alan
  RANGE R2C2
  RANGE R2CcAMP2
  RANGE R2cAMP4
  RANGE blah
  RANGE kcamp 
  RANGE r2
  RANGE taucamp
  RANGE vL2_sub
  RANGE vL1_sub
  RANGE vL3_sub
  RANGE peptotal 
  RANGE PYRgtotal 
  RANGE vnstim 
  RANGE vgstim
  RANGE ha		 			
  RANGE hb					
  RANGE hc		
  RANGE hd
  RANGE hk
  RANGE cAMP_K

}



PARAMETER {

 unit_compartment = 1.0
  V = 1e-07
  Fin = 0.012
  Jstimg = 0
  Jstimn = 0
  ADPn = 1e-10
  ADPg = 1e-10
  O2meanc = 7.0
  Vv0 = 0.021 ?
  JleakNan = 0.0
  JleakNag = 0.0
  Jpumpn = 0.0
  Jpumpg = 0.0
  JGLCce = 1e-10
  JGLCcg = 1e-10
  JGLCeg = 1e-10
  JGLCen = 1e-10
  JLACec = 1e-10
  JLACcg = 1e-10
  JLACeg = 1e-10
  JLACen = 1e-10
  JGLCec = 1e-10
  JGLCgc = 1e-10
  JGLCge = 1e-10
  JGLCne = 1e-10
  JLACce = 1e-10
  JLACgc = 1e-10
  JLACge = 1e-10
  JLACne = 1e-10
  JHKPFKn = 0.0
  JHKPFKg = 0.0
  JPGKn = 0.0
  JPGKg = 0.0
  JPKn = 0.0
  JPKg = 0.0
  JLDHn = 0.0
  JLDHg = 0.0
  Jmitoing =  0
  Jmitoinn =  0
  Jmitooutg =  0
  Jmitooutn =  0
  Jshuttleg = 0.0
  Jshuttlen = 0.0
  JCKg = 0.0
  JCKn = 0.0
  JO2mcg = 1e-10
  JO2mcn = 1e-10
  JO2c = 0.0
  JGLCc = 0.0
  JLACc = 0.0
  IL = 0.0
  INa = 1.0
  IK = 0.0
  ICa = 0.0
  ImAHP = 1e-10
  Ipump = 1e-10
  Fout = 0.012
  dAMPdATPn = 0.0
  dAMPdATPg = 0.0
  Ve = 0.2
  Vcap = 0.0055
  Vg = 0.25
  Vn = 0.45
  zeta = 0.07  
  qAK = 0.92
  A = 2.212
  un = 0.0
  ug = 0.0
  ren = 0.4444
  reg = 0.8
  rce = 0.0275
  rcn = 0.0122
  rcg = 0.0220
  Isyn = 0.0
  kCKplusn =   0.0433
  kCKplusg = 0.00135
  Nan0 = 8.0
  F0 = 0.012
  SmVn = 25000.0
  SmVg = 25000.0
  R = 8.314151
  F = 96485.3
  RToverF = 26.73
  psig = -70.0
  Nae = 150.0
  KtGLCen = 8.0
  KtGLCeg = 8.0
  KtGLCcg = 8.0
  KtGLCce = 8.0
  KtLACen = 0.74
  KtLACeg = 3.5
  KtLACge = 3.5
  KtLACgc = 1.0
  KtLACcg = 1.0
  KtLACec = 1.0
  KIATP = 1.0
  nH1 = 4.0 ? eredeti
  Kg = 0.05 ?eredeti
  KO2 = 0.0361
  HbOP = 8.6
  nh = 2.73
  KO2mito = 0.001
  Cm = 0.001
  gL = 0.02
  gNa = 40.0
  gK = 18.0
  gCa = 0.02
  gmAHP = 6.5
  KD = 0.03
  tauCa = 0.15
  Ca0 = 5e-05
  EK = -80.0
  ECa = 120.0
  phih = 4.0
  phin = 4.0
  tauv = 35.0
  alphav = 0.5
  O2a = 8.35
  GLCa = 4.75
  gNan = 0.0136
  gNag = 0.0061
  gKpas = 0.2035
  kpumpn = 2.2000e-06
  kpumpg = 4.5000e-07
  Kmpump = 0.5
  C = 10.0
  N = 0.212
  Kmmito = 0.04
  kLDHnplus =   72.3000
  kLDHnminus =    0.7200
  kLDHgplus =    1.5900 
  kLDHgminus =    0.0710 
  Mcyton =   4.9000e-08
  Mcytog =   2.5000e-04
  Mmiton = 393000
  Mmitog =       10600
  KmADPn =    0.00341
  KmADPg =   4.8300e-04
  KmNADHn =    0.0444
  KmNADHg =    0.0269
  kCKn =    0.0433
  kCKg =    0.00135
  KmNADn =    0.4090
  KmNADg =   40.3000
  TmaxGLCen =    0.0410
  TmaxGLCne =    0.0410
  TmaxGLCce =    0.2390
  TmaxGLCec =    0.2390
  TmaxGLCeg =    0.1470
  TmaxGLCge =    0.1470
  TmaxGLCcg =    0.0016
  TmaxLACcg =   0.00243
  TmaxLACne =   24.3000
  TmaxLACen =   24.3000
  TmaxLACge =  106.1000
  TmaxLACeg =  106.1000
  TmaxLACec =    0.2500
  TmaxLACce =    0.2500
  kHKPFKn = 0.050435
  kHKPFKg =  0.1850 
  PScapoverVn =    1.6600
  PScapoverVg =    0.8700
  Vmaxoutn = 0.1640
  Vmaxoutg = 0.0640
  Vmaxinn =0.1303
  Vmaxing =5.7000
  kPGKn =    3.9700
  kPGKg =  401.7000 
  kPKn =   36.7000
  kPKg =  135.2000  
  JATPasesn =    0.1695
  JATPasesg =    0.1404
  kCKminusn =   2.8000e-04
  kCKminusg =   1.0000e-05
  TNADHn =       10330
  TNADHg =   150
  LACa =    0.5060
  Rplusg = 0.0
  Rminusg = 0.0
  Rplusn = 0.0
  Rminusn = 0.0
  hh_alphan = 0.0
  hh_alpham = 0.0
  hh_betam = 0.0
  hh_alphah = 0.0
  hh_betah = 0.0
  hh_betan = 0.0
  hh_minfinity = 0.0
  hh_ninfinity = 0.0
  hh_hinfinity = 0.0
  hh_taun = 1.0
  hh_mCa = 0.0
  hh_tauh = 1.0
  hh_EL = 0.0
  Finprime=0
  JGLUeg=0
  vPumpg0 =    0.0687
  e		= 2.71828
  pi		= 3.1415
  kc1		= 1 			?
  kc2		= 1e-5	 		? mM/sec
  k11		= 0.000043		? mM
  k22		= 0.0007		? mM
  kcm1		= 4e-8 			? mM 			
  kcm2		= 1e-6 			? mM 			
  kgc1		= 1e-6 			? adjusted as per Boras ... McCulloch 2014.
  kgc2		= 1e-6 			? adjusted as per Boras ... McCulloch 2014.
  k_gc1		= 1e-2 			? adjusted as per Boras ... McCulloch 2014.
  k_gc2		= 1e-2 			? adjusted as per Boras ... McCulloch 2014.
  kd 		= 2e-06					 
  ni		= 10						
  ng		= 10			 				
  kDins		= 1e-6			? mM		
  kDins2		= 0.75e-6		?  mM	
  kDglucgn	= 4e-8			? mM 
  kDvip		= 4e-7
  kDvip2 		= 4e-3			 
  kDne		= 3.6e-5		? mM
  kfeed		= 0.01			? mM/sec 
  en1		= 10			
  en2		= 10			
  en5		= 10			
  en6		= 10			
  ep1		= 10			 	
  ep2		= 10			   	
  ep3		= 10			
  ep7		= 10			
  ep8		= 10			
  ep9		= 10 				
  ep12		= 10			 	
  ep13		= 10			 	
  ki13		= 2			 	?mM	
  kL1 		= .05			?	mM/sec 	  
  kL2 		= 3.4			?	/sec	   
  kL3		= .002			? 	mM/sec	  
  kL4		= .034			?	mM/sec 	 
  kL13		= 0.0083		?	/sec 
  k_L1		=  .07			?	mM/sec    
  k_L2		= .34			?	/min  	  
  k_L3 		= .002 			?	mM/sec	  
  kmL1		= 7.7			?	mM	
  km_L1		= 1.3			?	mM	
  kmL2		= 0.57			?	mM			
  km_L2		= 1.4			?	mM	
  kmL3		= 0.01 			?	mM	
  km_L3		= 0.0034		?	mM	
  kmL4		= 0.18			? 	mM 	
  ktL1		= 1.6 			?	/sec
  ktS1		= 0.001 		?  	/sec  
  ktF1		= 0.001 		? 	/sec
  kgi		= 10 			?  	mM	
  kg2		= 0.5			?	sec
  kg3		= 20 			?	/sec	
  kg4		= 5 			?	/sec   
  kg5		= 20 			?	/sec   
  kg6		= 5 			?	/sec	 	
  kg7		= 20 			?	/sec      
  kg8		= 5 		 	?	/min	 					 
  kmg3		= .004 			?	mM	 
  kmg4		= .0011 		?	mM	 
  kmg6		= .005 			?	mM  		
  kmg5		= 0.01			?	mM  
  kmg7		= 0.015			?	mM	 
  kmg8		= 0.00012 		?	mM
  kins		= 7e-4			?	mM/sec
  k1ins		= 6e-4			?	mm/sec
  kmIns		= 8			?	mM  
  kglucgn	= 2e-9			?	mM/sec
  kmGlgn	= 8			?	mM
  k1glucgn	= 5e-9			?	mM/sec
  kd_Bglucgn	= 0.015			?	/sec
  kd_Bins	= 0.015			?	/sec
  kd_Bgluc	= 0.00025		?	/sec   
  nn		= 1			? 	
  k_a 		= 1 			?	/sec*mM	
  kd_mg		= 1?5 			?	mM   
  kmind		= 2e-06			?	mM	            
  kmaxd		= 3.2e-03		?	mM	          
  s1		= 100 			? 	   
  s2		= .001 			?	
  st		= .003 			?	mM  
  pt		= .07 			?	mM  	
  kt		= .0025 		?	mM  
  taucamp	= 50000 		? 	(~ 50000 ms, M. Conti et al., 2014)
  ha		= 0  			? mM
  hb		= 9 			? mM
  hc		= 3  			? mM
  hd		= 1  			? unitless

}


STATE {

  Nan
  Nag
  GLCn
  GLCg
  GAPn
  GAPg
  PEPn
  PEPg
  PYRn
  PYRg
  LACn
  LACg
  NADHcyton
  NADHcytog
  NADHmiton
  NADHmitog
  ATPn
  ATPg
  PCrn
  PCrg
  O2n
  O2g
  O2c
  GLCc
  LACc
  Vv
  Hb
  GLCe
  LACe
  psin
  h
  n
  Ca
  Vvp
  glucgn
  glyc		
  GPa			
  GSa 			
  cAMP
  PP1		
  PP1_GPa 		
  PKa			
  gluc	
  gluco
  R2CcAMP2
  R2cAMP4
  CC	
  R2C2		
  ?B_glyc
  B_ins	
  B_gluc
  B_glucgn
  B_ne
  B_vip
  ka
  V_L1
  k_cg2
  kcg2
  vigluc 
  vfeed
  vd_Bins
  vins
  vd_Bglucgn
  vGgluc
  vglucgn
  vd_Bgluc
  vtS1
  vtF1
  vtL1 
  v_L1
  v_L2
  v_L3
  vL1
  vL2
  vL3
  vL4
  vL13
  cAMPtotal
  hk
  ne
  vip
  sigma
  g6p
  pep
  oa_c
  alan
  vnstim
  vgstim
  ge
  cAMP_K
  cAMP_vip2
  cAMP_ne
  cAMP_vip
  PYRgtotal
  peptotal
  vip2
  viptotal
  blah
  cAMPtotsyn

}



INITIAL {

V= -73.5001
Nan= 7.97118
Nag= 15.106
GLCn= 1.19525
GLCg= 1 ?1.1212
GAPn= 0.00460529
GAPg= 0.0574386
PEPn= 0.0164124
PEPg= 0.0279754
PYRn= 0.173885
PYRg= 0.202024
LACn= 0.599024
LACg= 0.602393
NADHcyton= 0.00625245
NADHcytog= 0.160762
ATPn= 2.19735
ATPg= 2.18717
PCrn= 4.94615
PCrg= 4.92419
O2n= 0.0268577
O2g= 0.027065
O2c= 6.9694
GLCc= 4.50504
LACc= 0.549196
Vv= 0.0210003
Vvp= -5.22878e-008
Hb= 0.0575034
GLCe= 2.49831
LACe= 0.600618
psin= -73.5001
h= 0.993603
n= 0.0187065
Ca= 5.10258e-005
NADHmiton= 0.123493
NADHmitog= 0.126012
Fout= 0.0119993
ADPn= 0.0145487
ADPg= 0.0245268
Rminusn = 0.0303889
Rplusn = 0.71669
Rminusg = 3.13758
Rplusg = 0.682384
JleakNag = 0.207621
JleakNan = 0.53544
Jpumpg = 0.0691609
Jpumpn = 0.178574
JHKPFKn = 0.0043752
JHKPFKg = 0.155087
ADPn = 0.0145487
ADPg = 0.0245268
JPGKn = 0.00875299
JPGKg = 0.180365
JPKn = 0.00876317
JPKg = 0.149144
JLDHn = -0.0101331
JLDHg = 0.0494482
JLACec = 0.00518442
JLACcg = -5.2074e-005
JLACeg = -0.0391787
JLACen = 0.0159702
Jmitoing = 0.0101304
Jmitoinn = 0.0188453
Jmitooutn = 0.0942165
Jmitooutg = 0.0498797
Jshuttlen = 0.0188381
Jshuttleg = 0.00965499
JCKn = 6.45101e-006
JCKg = 5.20289e-005
JO2mcn = 0.0574387
JO2mcg = 0.0299231
JO2c = 6.02442
O2meanc = 5.58881
hh_EL = -70.0742
IL = -0.0685182
hh_alpham = 0.0718109
hh_betam = 14.5556
hh_alphah = 0.734001
hh_betah = 0.00472564
hh_alphan = 0.00775484
hh_betan = 0.406799
hh_mCa = 0.00261347
hh_tauh = 0.00135368
hh_minfinity = 0.00490935
hh_ninfinity = 0.0187065
hh_hinfinity = 0.993603
hh_taun = 0.00241223
INa = -0.00071456
IK = 1.43266e-005
ICa = -2.64329e-005
ImAHP = 0.0717377
Ipump = -0.00249191
ren = 0.444444
reg = 0.8
rce = 0.0275
rcn = 0.0122222
rcg = 0.022
un= 0.870941
ug= 0.88817
dAMPdATPn= -0.0142363
dAMPdATPg= -0.0239314
INa = -0.00071456
IK = 1.43266e-005
ICa = -2.64329e-005
kd = 0.00161907
CC = 1.4339
PKa = 0.00899953
R2C2 = 0.0823673
R2CcAMP2 = 0.3584
R2cAMP4 = 1.55948
pep = 0.0170014
g6p = 0.541282
glyc = 4.88824
gluc = 1.1212?3.8787
cAMPtotal = 0.0449285
GSa = 0.0111569
GPa = 0.0699071
peptotal= 0
PYRgtotal=0
glyc 		= 20 
GPa		= .07 
GSa 		=  0.003 
PP1		= 0.00025			
PP1_GPa 	= 0.000025			
PKa		= 0.0025			
gluc 		= .1212 		
R2CcAMP2	= 1
R2cAMP4		= 1
CC		= 0.0018			
R2C2		= 0.00025			
B_ins		= 0.0				
B_gluc		= 0			
B_glucgn	= 0.0								
g6p		= .7 	?mM   	 	
alan		= 0.1				 
oa_c		= 0.1				
ne              = 0
vip             = 0
hk		= 0 			

}



BREAKPOINT {
   SOLVE states METHOD euler
  ? Need to check order in which assignments/event assignments should be updated!!!

  ? Assignment rule here: Rminusn = NADHcyton / (N - NADHcyton)
  Rminusn = NADHcyton / (N - NADHcyton)

  ? Assignment rule here: Rplusn = (N - NADHmiton) / NADHmiton
  Rplusn = (N - NADHmiton) / NADHmiton

  ? Assignment rule here: Rminusg = NADHcytog / (N - NADHcytog)
  Rminusg = NADHcytog / (N - NADHcytog)

  ? Assignment rule here: Rplusg = (N - NADHmitog) / NADHmitog
  Rplusg = (N - NADHmitog) / NADHmitog

  ? Assignment rule here: JleakNag = SmVg / F * gNag * (RToverF * log(Nae / Nag) - psig)
  JleakNag = SmVg / F * gNag * (RToverF * log(Nae / Nag) - psig)

  ? Assignment rule here: JleakNan = SmVn / F * gNan * (RToverF * log(Nae / Nan) - psin)
  JleakNan = SmVn / F * gNan * (RToverF * log(Nae / Nan) - psin)

  ? Assignment rule here: Jpumpg = SmVg * kpumpg * ATPg * Nag * pow(1 + ATPg / Kmpump, -1)
  Jpumpg = SmVg * kpumpg * ATPg * Nag * pow(1 + ATPg / Kmpump, -1)

  ? Assignment rule here: Jpumpn = SmVn * kpumpn * ATPn * Nan * pow(1 + ATPn / Kmpump, -1)
  Jpumpn = SmVn * kpumpn * ATPn * Nan * pow(1 + ATPn / Kmpump, -1)

  ? Assignment rule here: JGLCce = TmaxGLCce * (GLCc / (GLCc + KtGLCce) - GLCe / (GLCe + KtGLCce))
  JGLCce = TmaxGLCce * (GLCc / (GLCc + KtGLCce) - GLCe / (GLCe + KtGLCce))

  ? Assignment rule here: JGLCcg = TmaxGLCcg * (GLCc / (GLCc + KtGLCcg) - GLCg / (GLCg + KtGLCcg))
  JGLCcg = TmaxGLCcg * (GLCc / (GLCc + KtGLCcg) - GLCg / (GLCg + KtGLCcg))   
  
  ? Assignment rule here: JGLCeg = TmaxGLCeg * (GLCe / (GLCe + KtGLCeg) - GLCg / (GLCg + KtGLCeg))
  JGLCeg = TmaxGLCeg * (GLCe / (GLCe + KtGLCeg) - GLCg / (GLCg + KtGLCeg))  ? ere
  
  ? Assignment rule here: JGLCen = TmaxGLCen * (GLCe / (GLCe + KtGLCen) - GLCn / (GLCn + KtGLCen))
  JGLCen = TmaxGLCen * (GLCe / (GLCe + KtGLCen) - GLCn / (GLCn + KtGLCen))

  ? Assignment rule here: JHKPFKn = kHKPFKn * ATPn * (GLCn / (GLCn + Kg)) * pow(1 + pow(ATPn / KIATP, nH1), -1)
  JHKPFKn = kHKPFKn * ATPn * (GLCn / (GLCn + Kg)) * pow(1 + pow(ATPn / KIATP, nH1), -1)

  ? Assignment rule here: JHKPFKg = kHKPFKg * ATPg * (GLCg / (GLCg + Kg)) * pow(1 + pow(ATPg / KIATP, nH1), -1)
  JHKPFKg = kHKPFKg * ATPg * ((GLCg+g6p) / ((GLCg+g6p) + Kg)) * pow(1 + pow(ATPg / KIATP, nH1), -1)
  
  ? Assignment rule here: ADPn = ATPn / 2 * (pow(qAK * qAK + 4 * qAK * (A / ATPn - 1), 0.5) - qAK)
  ADPn = ATPn / 2 * (pow(qAK * qAK + 4 * qAK * (A / ATPn - 1), 0.5) - qAK)

  ? Assignment rule here: ADPg = ATPg / 2 * (pow(qAK * qAK + 4 * qAK * (A / ATPg - 1), 0.5) - qAK)
  ADPg = ATPg / 2 * (pow(qAK * qAK + 4 * qAK * (A / ATPg - 1), 0.5) - qAK)

  ? Assignment rule here: JPGKn = kPGKn * GAPn * ADPn * ((N - NADHcyton) / NADHcyton)
  JPGKn = kPGKn * GAPn * ADPn * ((N - NADHcyton) / NADHcyton)

  ? Assignment rule here: JPGKg = kPGKg * GAPg * ADPg * ((N - NADHcytog) / NADHcytog)
  JPGKg = kPGKg * GAPg * ADPg * ((N - NADHcytog) / NADHcytog) ?default 
  
  ? Assignment rule here: JPKn = kPKn * PEPn * ADPn
  JPKn = kPKn * PEPn * ADPn

  ? Assignment rule here: JPKg = kPKg * PEPg * ADPg
  JPKg = kPKg * PEPg * ADPg  

  ? Assignment rule here: JLDHn = kLDHnplus * PYRn * NADHcyton - kLDHnminus * LACn * (N - NADHcyton)
  JLDHn = kLDHnplus * PYRn * NADHcyton - kLDHnminus * LACn * (N - NADHcyton)

  ? Assignment rule here: JLDHg = kLDHgplus *   * NADHcytog - kLDHgminus * LACg * (N - NADHcytog)
  JLDHg = 2*(kLDHgplus * PYRg * NADHcytog - kLDHgminus * LACg * (N - NADHcytog)) ? 

  ? Assignment rule here: JLACec = TmaxLACec * (LACe / (LACe + KtLACec) - LACc / (LACc + KtLACec))
  JLACec = TmaxLACec * (LACe / (LACe + KtLACec) - LACc / (LACc + KtLACec))

  ? Assignment rule here: JLACcg = TmaxLACcg * (LACc / (LACc + KtLACcg) - LACg / (LACg + KtLACcg))
  JLACcg = TmaxLACcg * (LACc / (LACc + KtLACcg) - LACg / (LACg + KtLACcg))

  ? Assignment rule here: JLACeg = TmaxLACeg * (LACe / (LACe + KtLACeg) - LACg / (LACg + KtLACeg))
  JLACeg = TmaxLACeg * (LACe / (LACe + KtLACeg) - LACg / (LACg + KtLACeg))

  ? Assignment rule here: JLACen = TmaxLACen * (LACe / (LACe + KtLACen) - LACn / (LACn + KtLACen))
  JLACen = TmaxLACen * (LACe / (LACe + KtLACen) - LACn / (LACn + KtLACen))

  ? Assignment rule here: Jmitoing = Vmaxing * PYRg / (PYRg + Kmmito) * ((N - NADHmitog) / (N - NADHmitog + KmNADHg))
  Jmitoing = Vmaxing * PYRg / (PYRg + Kmmito) * ((N - NADHmitog) / (N - NADHmitog + KmNADg))

  ? Assignment rule here: Jmitoinn = Vmaxinn * PYRn / (PYRn + Kmmito) * ((N - NADHmiton) / (N - NADHmiton + KmNADHn))
  Jmitoinn = Vmaxinn * PYRn / (PYRn + Kmmito) * ((N - NADHmiton) / (N - NADHmiton + KmNADn))

  ? Assignment rule here: Jmitooutn = Vmaxoutn * O2n / (O2n + KO2mito) * ADPn / (ADPn + KmADPn) * (NADHmiton / (NADHmiton + KmNADHn))
  Jmitooutn = Vmaxoutn * O2n / (O2n + KO2mito) * ADPn / (ADPn + KmADPn) * (NADHmiton / (NADHmiton + KmNADHn))

  ? Assignment rule here: Jmitooutg = Vmaxoutg * O2g / (O2g + KO2mito) * ADPg / (ADPg + KmADPg) * (NADHmitog / (NADHmitog + KmNADHg))
  Jmitooutg = Vmaxoutg * O2g / (O2g + KO2mito) * ADPg / (ADPg + KmADPg) * (NADHmitog / (NADHmitog + KmNADHg))

  ? Assignment rule here: Jshuttlen = TNADHn * Rminusn / (Mcyton + Rminusn) * (Rplusn / (Mmiton + Rplusn))
  Jshuttlen = TNADHn * Rminusn / (Mcyton + Rminusn) * (Rplusn / (Mmiton + Rplusn))

  ? Assignment rule here: Jshuttleg = TNADHg * Rminusg / (Mcytog + Rminusg) * (Rplusg / (Mmitog + Rplusg))
  Jshuttleg = TNADHg * Rminusg / (Mcytog + Rminusg) * (Rplusg / (Mmitog + Rplusg))

  ? Assignment rule here: JCKn = kCKplusn * PCrn * ADPn - kCKminusn * (C - PCrn) * ATPn
  JCKn = kCKplusn * PCrn * ADPn - kCKminusn * (C - PCrn) * ATPn

  ? Assignment rule here: JCKg = kCKplusg * PCrg * ADPg - kCKminusg * (C - PCrg) * ATPg
  JCKg = kCKplusg * PCrg * ADPg - kCKminusg * (C - PCrg) * ATPg

  ? Assignment rule here: JO2mcn = PScapoverVn * (KO2 * pow(HbOP / O2c - 1, -1 / nh) - O2n)
  JO2mcn = PScapoverVn * (KO2 * pow(HbOP / O2c - 1, -1 / nh) - O2n)

  ? Assignment rule here: JO2mcg = PScapoverVg * (KO2 * pow(HbOP / O2c - 1, -1 / nh) - O2g)
  JO2mcg = PScapoverVg * (KO2 * pow(HbOP / O2c - 1, -1 / nh) - O2g)

  ? Assignment rule here: JO2c = 2 * Fin / Vcap * (O2a - O2c)
  JO2c = 2 * Fin / Vcap * (O2a - O2c)

  ? Assignment rule here: JGLCc = 2 * Fin / Vcap * (GLCa - GLCc)
  JGLCc = 2 * Fin / Vcap * (GLCa - GLCc)

  ? Assignment rule here: JLACc = 2 * Fin / Vcap * (LACa - LACc)
  JLACc = 2 * Fin / Vcap * (LACa - LACc)

  ? Assignment rule here: O2meanc = 2 * O2c - O2a
  O2meanc = 2 * O2c - O2a

  ? Assignment rule here: hh_EL = (gKpas * EK + gNan * RToverF * log(Nae / Nan)) / (gKpas + gNan)
  hh_EL = (gKpas * EK + gNan * RToverF * log(Nae / Nan)) / (gKpas + gNan)

  ? Assignment rule here: IL = gL * (psin - hh_EL)
  IL = gL * (psin - hh_EL)

  ? Assignment rule here: hh_alpham = -0.1 * ((psin + 33) / (exp(-0.1 * (psin + 33)) - 1))
  hh_alpham = -0.1 * ((psin + 33) / (exp(-0.1 * (psin + 33)) - 1))

  ? Assignment rule here: hh_betam = 4 * exp((psin + 58) / -12)
  hh_betam = 4 * exp((psin + 58) / -12)

  ? Assignment rule here: hh_alphah = 0.07 * exp((0 - (psin + 50)) / 10)
  hh_alphah = 0.07 * exp((0 - (psin + 50)) / 10)

  ? Assignment rule here: hh_betah = 1 / (exp(-0.1 * (psin + 20)) + 1)
  hh_betah = 1 / (exp(-0.1 * (psin + 20)) + 1)

  ? Assignment rule here: hh_alphan = -0.01 * ((psin + 34) / (exp(-0.1 * (psin + 34)) - 1))
  hh_alphan = -0.01 * ((psin + 34) / (exp(-0.1 * (psin + 34)) - 1))

  ? Assignment rule here: hh_betan = 0.125 * exp((0 - (psin + 44)) / 25)
  hh_betan = 0.125 * exp((0 - (psin + 44)) / 25)

  ? Assignment rule here: hh_mCa = 1 / (exp((psin + 20) / -9) + 1)
  hh_mCa = 1 / (exp((psin + 20) / -9) + 1)

  ? Assignment rule here: hh_tauh = 0.001 / (hh_alphah + hh_betah)
  hh_tauh = 0.001 / (hh_alphah + hh_betah)

  ? Assignment rule here: hh_minfinity = hh_alpham / (hh_alpham + hh_betam)
  hh_minfinity = hh_alpham / (hh_alpham + hh_betam)

  ? Assignment rule here: hh_ninfinity = hh_alphan / (hh_alphan + hh_betan)
  hh_ninfinity = hh_alphan / (hh_alphan + hh_betan)

  ? Assignment rule here: hh_hinfinity = hh_alphah / (hh_alphah + hh_betah)
  hh_hinfinity = hh_alphah / (hh_alphah + hh_betah)

  ? Assignment rule here: hh_taun = 0.001 / (hh_alphan + hh_betan)
  hh_taun = 0.001 / (hh_alphan + hh_betan)

  ? Assignment rule here: INa = gNa * pow(hh_minfinity, 3) * h * (psin - RToverF * log(Nae / Nan))
  INa = gNa * pow(hh_minfinity, 3) * h * (psin - RToverF * log(Nae / Nan))

  ? Assignment rule here: IK = gK * pow(n, 4) * (psin - EK)
  IK = gK * pow(n, 4) * (psin - EK)

  ? Assignment rule here: ICa = gCa * pow(hh_mCa, 2) * (psin - ECa)
  ICa = gCa * pow(hh_mCa, 2) * (psin - ECa)

  ? Assignment rule here: ImAHP = gmAHP * Ca / (Ca + KD) * (psin - EK)
  ImAHP = gmAHP * Ca / (Ca + KD) * (psin - EK)

  ? Assignment rule here: Ipump = F * kpumpn * ATPn * (Nan - Nan0) * pow(1 + ATPn / Kmpump, -1)
  Ipump = F * kpumpn * ATPn * (Nan - Nan0) * pow(1 + ATPn / Kmpump, -1)

  ? Assignment rule here: Fout = F0 * pow(Vv / Vv0, 1 / alphav) + F0 * (tauv / Vv0) * pow(Vv / Vv0, -(1 / 2))
  Fout =  F0 * pow(Vv / Vv0, 1 / alphav) + F0 * (tauv / Vv0) * Vvp * pow( Vv / Vv0 , -(1 / 2))

  ? Assignment rule here: ren = Ve / Vn
  ren = Ve / Vn

  ? Assignment rule here: reg = Ve / Vg
  reg = Ve / Vg

  ? Assignment rule here: rce = Vcap / Ve
  rce = Vcap / Ve

  ? Assignment rule here: rcn = Vcap / Vn
  rcn = Vcap / Vn

  ? Assignment rule here: rcg = Vcap / Vg
  rcg = Vcap / Vg

  ? Assignment rule here: JGLCec = 0 - JGLCce
  JGLCec = 0 - JGLCce

  ? Assignment rule here: JGLCgc = 0 - JGLCcg
  JGLCgc = 0 - JGLCcg

  ? Assignment rule here: JGLCge = 0 - JGLCeg
  JGLCge = 0 - JGLCeg

  ? Assignment rule here: JGLCne = 0 - JGLCen
  JGLCne = 0 - JGLCen

  ? Assignment rule here: JLACce = 0 - JLACec
  JLACce = 0 - JLACec

  ? Assignment rule here: JLACgc = 0 - JLACcg
  JLACgc = 0 - JLACcg

  ? Assignment rule here: JLACge = 0 - JLACeg
  JLACge = 0 - JLACeg

  ? Assignment rule here: JLACne = 0 - JLACen
  JLACne = 0 - JLACen

  ? Assignment rule here: un = qAK * qAK + 4 * qAK * (A / ATPn - 1)
  un = qAK * qAK + 4 * qAK * (A / ATPn - 1)

  ? Assignment rule here: ug = qAK * qAK + 4 * qAK * (A / ATPg - 1)
   ug = qAK * qAK + 4 * qAK * (A / ATPg - 1)

  ? Assignment rule here: dAMPdATPn = -1 + qAK / 2 + -0.5 * pow(un, 0.5) + qAK * A / (ATPn * pow(un, 0.5))
  dAMPdATPn = -1 + qAK / 2 + -0.5 * pow(un, 0.5) + qAK * A / (ATPn * pow(un, 0.5))

  ? Assignment rule here: dAMPdATPg = -1 + qAK / 2 + -0.5 * pow(ug, 0.5) + qAK * A / (ATPg * pow(ug, 0.5))
  dAMPdATPg = -1 + qAK / 2 + -0.5 * pow(ug, 0.5) + qAK * A / (ATPg * pow(ug, 0.5))
  
kd	 		= (kmaxd - kmind)*(1/( 1+ pow(glyc/kd_mg,nn))) + kmind     
ka 			= k_a / kd 
vL1 			= (kL1*gluco / (kmL1+gluco) )   
vL2 			= kL2*GSa*g6p / (kmL2+g6p)				
v_L1 			= (k_L1*g6p / (km_L1+g6p) ) 
v_L2			= k_L2*GPa *glyc / (km_L2+glyc)			
vL3			= (kL3*g6p/(kmL3+g6p))
v_L3                    = (k_L3*PEPg/(km_L3+PEPg))  
vL4			= (kL4*pep/(kmL4+pep))
vL13			= kL13*oa_c 
vtL1			= 0 	
vtF1			= ktF1*B_gluc *(1+ ( pow(B_ins, ep12) / (pow(kDins, ep12)+ pow(B_ins,ep12)))) 
vtS1			= ktS1*B_gluc *(1+ ( pow(B_ins, ep13) /  (pow(kDins, ep13)+ pow(B_ins,ep13)))) 
vd_Bgluc		= kd_Bgluc*B_gluc
vglucgn 		= kglucgn
vGgluc			=(k1glucgn*pow(B_gluc, ng)) / ( pow(kmGlgn, ng) + pow(B_gluc, ng)) ?
vd_Bglucgn		= kd_Bglucgn*B_glucgn
vins			= kins
vigluc			= (k1ins*pow(B_gluc, ni)) / ( pow(kmIns, ni) + pow(B_gluc, ni))   
vd_Bins		 	= kd_Bins*B_ins 
gluco			= gluc + GLCg
cAMP_K			= ha + ((hb-ha)/(1+ pow(hc/hk,hd))) ?is this right?
?cAMP_ne		=  (ne/(kDne + ne))  ?.001 is vvlow, .1 is vlow, 1 is norm, 10 is hi, 100 is vhi?
?cAMP_vip		=  (vip/(kDvip + vip))  
?cAMP_vip2		=  (vip2/(kDvip2 + vip2)) 
blah			= sin(t/100000)
?cAMPtotsyn		= cAMP + 2.5*(cAMP_ne + cAMP_vip + cAMP_vip2) ? per Magistretti (1984)

}


 
DERIVATIVE states {

  Nag' = 1/1000*( JleakNag - 3 * Jpumpg + Jstimg + 3 * JGLUeg)
  Nan' = 1/1000*( JleakNan - 3 * Jpumpn + Jstimn -0.310880 *INa)
  GLCn' =1/1000*( JGLCen - JHKPFKn)
  GLCg' =1/1000*( JGLCcg + JGLCeg - JHKPFKg) 
  GAPg' =1/1000*( 2 * JHKPFKg - JPGKg) ? 
  GAPn' =1/1000*( 2 * JHKPFKn - JPGKn) 
  PEPg' =1/1000*( JPGKg - JPKg)
  PEPn' =1/1000*( JPGKn - JPKn)
  PYRg' = 1/1000*(JPKg - JLDHg - Jmitoing)
  PYRn' =1/1000*( JPKn - JLDHn - Jmitoinn)
  LACn' =1/1000*( JLDHn - JLACne) 
  LACg' =1/1000*( JLDHg - JLACge - JLACgc)
  NADHcyton' =1/1000*( pow(1 - zeta, -1) * (JPGKn - JLDHn - Jshuttlen))
  NADHcytog' =1/1000*( pow(1 - zeta, -1) * (JPGKg - JLDHg - Jshuttleg))
  NADHmiton' =1/1000*( pow(zeta, -1) * (4 * Jmitoinn - Jmitooutn + Jshuttlen))
  NADHmitog' = 1/1000*(pow(zeta, -1) * (4 * Jmitoing - Jmitooutg + Jshuttleg))
  ATPn' = 1/1000*((-2 * JHKPFKn + JPGKn + JPKn - JATPasesn - 1 * Jpumpn + 3.6 * Jmitooutn + JCKn) * pow(1 - dAMPdATPn, -1))
  ATPg' = 1/1000*((-2 * JHKPFKg + JPGKg + JPKg - JATPasesg - 1.75 * Jpumpg + 0.75 * vPumpg0+3.6 * Jmitooutg + JCKg) * pow(1 - dAMPdATPg, -1))
  PCrg' = 1/1000*(-JCKg)
  PCrn' = 1/1000*(-JCKn)
  O2n' = 1/1000*(JO2mcn - 0.625 * Jmitooutn)
  O2g' = 1/1000*(JO2mcg - 0.625 * Jmitooutg)
  O2c' = 1/1000*(JO2c - 1 / rcn * JO2mcn - 1 / rcg * JO2mcg)
  GLCc' =1/1000*( JGLCc - 1 / rce * JGLCce - 1 / rcg * JGLCcg)
  LACc' =1/1000*(JLACc + 1 / rce * JLACec + 1 / rcg * JLACgc)
  Vv' = 1/1000*((Fin -  Fout))
  Hb' =  1/1000*(Fin * (O2a - O2meanc) * Vv' /(Fout+Vv'))
  GLCe' =1/1000*(  JGLCce - 1 / reg * JGLCeg -  1 / ren * JGLCen)
  LACe' = 1/1000*( 1 / ren * JLACne + 1 / reg * JLACge -  JLACec)
  psin' =1/1000*( 1 / Cm * (-IL - INa - IK - ICa - ImAHP - Ipump + Isyn))
  h' = 1/1000*(phih / hh_tauh * (hh_hinfinity - h))
  n' = 1/1000*( phin / hh_taun * (hh_ninfinity - n))
  Ca' =1/1000*( -SmVn * ICa / F - 1 / tauCa * (Ca - Ca0))
  ?Vvp should be the same as Vv', for constant Fin. Vv'=Fin-Fout, so Vv'=Fout'. We worked thought this and solved for Vv'', which is set to Vvp in the following equation
  ?Vvp'=(F0 *Vv' *(-2 *Vv0 *pow(Vv/Vv0,(1/2 + 1/alphav)) + alphav* tauv* Vv'))/(2 *alphav *Vv *(F0 *tauv + Vv0* pow(Vv/Vv0,1/2))) )
  Vvp'=1/1000*((2 *alphav *Finprime *Vv0 *Vv *pow((Vv/Vv0),1/2) -  2 *F0 *Vv0 *pow((Vv/Vv0),(1/2 + 1/alphav)) *Vv' + alphav *F0 *tauv *Vv' *Vv')/(2 *alphav *Vv *(F0 * tauv + Vv0 * pow(Vv/Vv0,1/2) ) 
))

alan'			= 0   ?  vtL6 + vL21 - v_L21 - vL20   
oa_c'			= 0   ?  vL12 - vL13 + vL14		
cAMPtotal'		= ( .0001*(ne/(kDne + ne)) + .5*(1/taucamp * ((hb-ha)/(1+ pow(hc/hk,hd))) + ha/taucamp ) - (2*kgc1*R2C2*pow(cAMPtotal,2)) + (2*k_gc1*R2CcAMP2*CC) - (2*kgc2*R2CcAMP2 * pow(cAMPtotal,2)) + (2*k_gc2 *R2cAMP4*CC) ) + cAMPtotal*(-1/taucamp) ? combo
glyc' 			= vL2 - v_L2  
GPa'			= ((kg5*PKa*(pt-GPa))/(kmg5*(1+s1*g6p/kg2) + (pt-GPa))) - ((kg6*(PP1+PP1_GPa)*GPa)/(kmg6/(1+s2*gluc/kgi)+GPa)) - (ka*PP1*GPa) + (k_a*PP1_GPa) 
GSa'			= ((kg8*PP1*(st-GSa))/((kmg8/(1+s1*g6p/kg2)) + (st-GSa))) - ((kg7*(PKa+CC)*GSa) /(kmg7*(1+s1*g6p/kg2)+GSa))
PP1' 			= 1*((-ka*PP1*GPa) + (k_a*PP1_GPa))
PP1_GPa' 		= (ka*PP1*GPa) - (k_a*PP1_GPa) 	
PKa'			= ((kg3*CC*(kt-PKa))/(kmg3+kt-PKa)) - (((kg4*(PP1+PP1_GPa)*PKa))/(kmg4+PKa))
R2CcAMP2'		= kgc1*R2C2*pow(cAMPtotal, 2) - k_gc1*R2CcAMP2*CC - kgc2*R2CcAMP2*pow(cAMPtotal, 2) +  k_gc2*R2cAMP4*CC 
R2cAMP4'		= kgc2*R2CcAMP2*pow(cAMPtotal, 2) - k_gc2*R2cAMP4*CC ? default
R2C2'			= -kgc1*R2C2*pow(cAMPtotal, 2) + k_gc1*R2CcAMP2*CC 
CC'			= kgc1*R2C2*pow(cAMPtotal, 2) - k_gc1*R2CcAMP2*CC + kgc2*R2CcAMP2*pow(cAMPtotal, 2) -  k_gc2*R2cAMP4*CC 
gluc' 			= vtL1 - vL1 + v_L1 
g6p'			= vL1 - v_L1 - vL2 + v_L2 
?B_ins'			= vins + vigluc - vd_Bins  
B_ins'			= 0					
?B_gluc'		= vfeed - vtL1 - vtF1 - vtS1 - vd_Bgluc		
B_gluc'			= 0					
?B_glucgn'		= vglucgn - vGgluc - vd_Bglucgn 
B_glucgn'		= 0 				
  
}

