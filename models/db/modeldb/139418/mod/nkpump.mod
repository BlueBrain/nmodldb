TITLE Sodium Potassium Pump

UNITS {
      (molar) = (1/liter)
      (mV) = (millivolts)
      (mA) = (milliamp)
      (mM) = (millimolar)
}

NEURON {
       SUFFIX pump
       USEION k READ ko WRITE ik
       RANGE ipumpmax, kbath
}

PARAMETER {
	  kbath	= 9		 (mM)
	  ko		         (mM)
	  ipumpmax = 0.0005     (mA/cm2)
	  ipump			 (mA/cm2)
}

ASSIGNED {
	 ik	 (mA/cm2)
}

BREAKPOINT {
	   ipump =  ipumpmax/((1+(kbath/ko))^2)
	   if (ipump <=0) {
	     ipump=0
	     }
	   ik = -2*ipump

}

