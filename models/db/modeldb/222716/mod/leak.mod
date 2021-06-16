TITLE passive NA, K, and Ca channel leak conductance

COMMENT
        Assembled for MyFirstNEURON by Arthur Houweling
        Modified by William L. Kath
ENDCOMMENT

UNITS {
        (mV) = (millivolt)
        (mA) = (milliamp)
        (mM) = (milli/liter)
        (S) = (siemens)
        FARADAY = 96480 (coul)
        R       = 8.314 (volt-coul/degC)
}

NEURON {
        SUFFIX lc
        USEION na READ nai,nao WRITE ina
        USEION k READ ki,ko WRITE ik
        USEION ca READ cai,cao WRITE ica
	 : USEION cl READ cli,clo WRITE icl VALENCE -1
        RANGE ina,ik,ica,icl
        RANGE gna,gk,gca,gcl
        RANGE cnar,ckr,ccar,cclr
}

PARAMETER {
        v               (mV)    
        celsius         (degC)
        nai             (mM)
        nao             (mM)
        pna= 2.07e-7    (cm/s)
        gna=0.000005     (S/cm2)
        cnar=38.38      (mM)        :=150*15*log(150/15)/(150-15) 
        ki              (mM)
        ko              (mM)
        pk= 3.45e-6     (cm/s)
        gk=0.00005      (S/cm2)
        ckr=15.77       (mM)        :=5*100*log(5/100)/(5-100) 
        cai             (mM)
        cao             (mM)
        pca= 0          (cm/s)
        gca=0           (S/cm2)
        ccar=0.00184    (mM)        :=2*0.0002*log(2/0.0002)/(2-0.0002) 
	:cli		(mM)
	:clo		(mM)
	:pcl=1.5525e-6  (cm/s)
	:gcl=0.00001	(S/cm2)
	:cclr=21.123	(mM)	    :=120*7*log(120/7)/(120-7)
}

ASSIGNED {
        ina     (mA/cm2)
        ik      (mA/cm2)
        ica     (mA/cm2)
	:icl	(mA/cm2)
        v0      (mV)
        vnar    (mV)
        vkr     (mV)
        vcar    (mV)
	vclr	(mV)
}

BREAKPOINT {
        v0 = (1000)*R*(celsius+273.16)/(FARADAY)
        vnar = v0*log(nao/nai)
        ina = gna * nao*nai*log(nao/nai)/(cnar*(nao-nai))* ghkfact(v,nai,nao,1)*(v-vnar)
        vkr = v0*log(ko/ki)
        ik = gk * ko*ki*log(ko/ki)/(ckr*(ko-ki))* ghkfact(v,ki,ko,1)*(v-vkr)
        vcar = v0*log(cao/cai)/2
        ica = gca * cao*cai*log(cao/cai)/(ccar*(cao-cai))*ghkfact(v,cai,cao,2)*(v-vcar)
	:vclr = v0*log(cli/clo)
	:icl = -gcl * clo*cli*log(clo/cli)/(cclr*(clo-cli))*ghkfact(v,cli,clo,-1)*(v-vclr)
}

FUNCTION expfun( v(mV), v0(mV)) (1) { 
        LOCAL e, w
        w=v/v0
        if (fabs(w)>1e-4) 
          { e = w / (exp(w)-1) }
        else 
        : denominator is small -> Taylor series
          { e = 1-w/2+w^2/12 }
        expfun=e
}

FUNCTION ghkfact(v(mV), ci(mM), co(mM), z)  (1) { 
        LOCAL v0, vr
        v0 = (1000)*R*(celsius+273.16)/(z*FARADAY)
        vr = v0*log(co/ci)
        ghkfact = expfun(v,v0)/(expfun(v-vr,v0)*expfun(vr,v0))
}



