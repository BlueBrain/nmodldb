TITLE passive sodium, potassium, calcium, chloride and magnesium channels

COMMENT
        Assembled for MyFirstNEURON by Arthur Houweling
ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(mM) = (milli/liter)
        FARADAY = 96480 (coul)
        R       = 8.314 (volt-coul/degC)
}

NEURON {
	SUFFIX leak
	USEION na READ nai,nao WRITE ina
	USEION k READ ki,ko WRITE ik
	USEION ca READ cai,cao WRITE ica
        USEION cl READ cli,clo WRITE icl VALENCE 1 
	: negative valence not accepted
	: if 'ecl' is needed in a MOD-file use '-ecl' instead
        USEION mg READ mgi,mgo WRITE img VALENCE 2
	RANGE pna,pk,pca,pcl,pmg,ina,ik,ica,icl,img
}

PARAMETER {
	v		(mV)	
        celsius		(degC)
	nai		(mM)
	nao		(mM)
	pna= 2.07e-7	(cm/s)
        ki      	(mM)
        ko      	(mM)
        pk= 3.45e-6	(cm/s)
        cai     	(mM)
        cao     	(mM)
        pca= 0		(cm/s)
        cli     	(mM)
        clo     	(mM)
        pcl= 0		(cm/s)
        mgi     	(mM)
        mgo     	(mM)
        pmg= 0		(cm/s)
}

ASSIGNED {
	ina	(mA/cm2)
	ik      (mA/cm2)
	ica	(mA/cm2)
        icl     (mA/cm2)
	img	(mA/cm2)
}

BREAKPOINT {
	ina = pna * ghk(v,nai,nao,1)
        ik = pk * ghk(v,ki,ko,1)
        ica = pca * ghk(v,cai,cao,2)
        icl = pcl * ghk(v,cli,clo,-1)
        img = pmg * ghk(v,mgi,mgo,2)
}

FUNCTION ghk( v(mV), ci(mM), co(mM), z)  (millicoul/cm3) { LOCAL e, w
        w = v * (.001) * z*FARADAY / (R*(celsius+273.16))
        if (fabs(w)>1e-4) 
          { e = w / (exp(w)-1) }
        else 
	: denominator is small -> Taylor series
          { e = 1-w/2 }
        ghk = - (.001) * z*FARADAY * (co-ci*exp(w)) * e
}

