TITLE SBML model: glia_2013 generated from file: Glia.xml

? this mod file represents a convertion to NEURON from the Matlab model originally published by Jolivet et al., PLoS Comp Biol (2015)

UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
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
RANGE vnstim
RANGE vgstim

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
  Vv0 = 0.021 
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
  nH1 = 4.0
  Kg = 0.05
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
  gNan =    0.0136
  gNag =    0.0061
  gKpas =    0.2035
  kpumpn =   2.2000e-06
  kpumpg =   4.5000e-07
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
  kCKg =    0.0014
  KmNADn =    0.4090
  KmNADg =   40.3000
  TmaxGLCen =    0.0410
  TmaxGLCne =    0.0410
  TmaxGLCce =    0.2390
  TmaxGLCec =    0.2390
  TmaxGLCeg =    0.1470
  TmaxGLCge =    0.1470
  TmaxGLCcg =    0.0016
  TmaxLACcg = 0.0026
  TmaxLACne =   24.3000
  TmaxLACen =   24.3000
  TmaxLACge =  106.1000
  TmaxLACeg =  106.1000
  TmaxLACec =    0.2500
  TmaxLACce =    0.2500
  kHKPFKn = 0.050435
  kHKPFKg =    0.1850
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

}

INITIAL {

 Fout= 0.012
 Vvp= -5.22695e-008
 Nan=   7.9717
 Nag=   15.1060
 GLCn=   1.1952
 GLCg=    1.1917
 GAPn=   0.0046
 GAPg=    0.0046
 PEPn=    0.0164
 PEPg=    0.0142
 PYRn=   0.1739
 PYRg=   0.1630
 LACn=   0.5981
 LACg=   0.6001
 NADHcyton=   0.0062
 NADHcytog=   0.1039
 ATPn=    2.1974
 ATPg=  2.1951
 PCrn=    4.9460
 PCrg=    4.9242
 O2n=   0.0276
 O2g=   0.0276
 O2c=   6.9809
 GLCc=   4.5004
 LACc=    0.5489
 Vv=    0.0210
 Hb=   0.0575
 GLCe=   2.4797
 LACe=    0.5991
 psin= -73.5928
 h=    0.9937
 n=    0.0185
 Ca=    0.000051
 NADHmiton=   0.1235
 NADHmitog=   0.1248
 Nan=     7.971123690827819
 Nag=   15.105962700710011
 GLCn=   1.195242817701208
 GLCg= 1.191683486104078
 GAPn=  0.004604883633382
 GAPg=  0.004592187611272
 PEPn=  0.016413198236988
 PEPg=   0.014203938186866
 PYRn=  0.173879012167164
 PYRg= 0.162979341543002
 LACn=   0.598139671517830
 LACg=   0.600154193323652
 NADHcyton= 0.006244551091528
 NADHcytog= 0.103869031492734
 ATPn=   2.197379714751472
 ATPg=  2.195100303412994
 PCrn=   4.946145277066859
 PCrg=   4.924205857343780
 O2n=   0.027586286704408
 O2g=  0.027594056563463
 O2c=  6.980872113185740
 GLCc=   4.500386628495174
 LACc=  0.548850780742499
 Vv=   0.021000000000000
 Hb=   0.057503371246198
 GLCe=  2.479691528660544
 LACe=   0.599146935390395
 psin= -73.578843708782614
 h=   0.993702020434757
 n=   0.018538671314650
 Ca=    0.000051009989197
 NADHmiton=  0.123456504032627
 NADHmitog=    0.124833635085897

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
  JGLCeg = TmaxGLCeg * (GLCe / (GLCe + KtGLCeg) - GLCg / (GLCg + KtGLCeg))

  ? Assignment rule here: JGLCen = TmaxGLCen * (GLCe / (GLCe + KtGLCen) - GLCn / (GLCn + KtGLCen))
  JGLCen = TmaxGLCen * (GLCe / (GLCe + KtGLCen) - GLCn / (GLCn + KtGLCen))

  ? Assignment rule here: JHKPFKn = kHKPFKn * ATPn * (GLCn / (GLCn + Kg)) * pow(1 + pow(ATPn / KIATP, nH1), -1)
  JHKPFKn = kHKPFKn * ATPn * (GLCn / (GLCn + Kg)) * pow(1 + pow(ATPn / KIATP, nH1), -1)

  ? Assignment rule here: JHKPFKg = kHKPFKg * ATPg * (GLCg / (GLCg + Kg)) * pow(1 + pow(ATPg / KIATP, nH1), -1)
  JHKPFKg = kHKPFKg * ATPg * (GLCg / (GLCg + Kg)) * pow(1 + pow(ATPg / KIATP, nH1), -1)

  ? Assignment rule here: ADPn = ATPn / 2 * (pow(qAK * qAK + 4 * qAK * (A / ATPn - 1), 0.5) - qAK)
  ADPn = ATPn / 2 * (pow(qAK * qAK + 4 * qAK * (A / ATPn - 1), 0.5) - qAK)

  ? Assignment rule here: ADPg = ATPg / 2 * (pow(qAK * qAK + 4 * qAK * (A / ATPg - 1), 0.5) - qAK)
  ADPg = ATPg / 2 * (pow(qAK * qAK + 4 * qAK * (A / ATPg - 1), 0.5) - qAK)

  ? Assignment rule here: JPGKn = kPGKn * GAPn * ADPn * ((N - NADHcyton) / NADHcyton)
  JPGKn = kPGKn * GAPn * ADPn * ((N - NADHcyton) / NADHcyton)

  ? Assignment rule here: JPGKg = kPGKg * GAPg * ADPg * ((N - NADHcytog) / NADHcytog)
  JPGKg = kPGKg * GAPg * ADPg * ((N - NADHcytog) / NADHcytog)

  ? Assignment rule here: JPKn = kPKn * PEPn * ADPn
  JPKn = kPKn * PEPn * ADPn

  ? Assignment rule here: JPKg = kPKg * PEPg * ADPg
  JPKg = kPKg * PEPg * ADPg

  ? Assignment rule here: JLDHn = kLDHnplus * PYRn * NADHcyton - kLDHnminus * LACn * (N - NADHcyton)
  JLDHn = kLDHnplus * PYRn * NADHcyton - kLDHnminus * LACn * (N - NADHcyton)

  ? Assignment rule here: JLDHg = kLDHgplus * PYRg * NADHcytog - kLDHgminus * LACg * (N - NADHcytog)
  JLDHg = kLDHgplus * PYRg * NADHcytog - kLDHgminus * LACg * (N - NADHcytog)

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
}

DERIVATIVE states {
  Nag' =1/1000*(  JleakNag - 3 * Jpumpg + Jstimg + 3 * JGLUeg)
  Nan' = 1/1000*( JleakNan - 3 * Jpumpn + Jstimn -0.310880 *INa)
  GLCn' =1/1000*( JGLCen - JHKPFKn)
  GLCg' =1/1000*( JGLCcg + JGLCeg - JHKPFKg)
  GAPg' =1/1000*( 2 * JHKPFKg - JPGKg)
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
  Vvp'=1/1000*((2 *alphav *Finprime *Vv0 *Vv *pow((Vv/Vv0),1/2) -  2 *F0 *Vv0 *pow((Vv/Vv0),(1/2 + 1/alphav)) *Vv' + alphav *F0 *tauv *Vv' *Vv')/(2 *alphav *Vv *(F0 * tauv + Vv0 * pow(Vv/Vv0,1/2) ) ))
}


