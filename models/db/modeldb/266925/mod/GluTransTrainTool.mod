TITLE GlutTransTrainTool

: Generates a decaying Glut time course at a given rate
: using various functions.
: 


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX glurel
	RANGE T, A1, tau1, A2, tau2, A3, tau3, D, Twait, aA, atau, expfunctrain, PulseInterval, PulseNum, alphafunc, difffunc, distance, diffCoeff, numMolecules, avogadro, cleftWidth, diff3dfunc, distance3d, diffCoeff3d, numMolecules3d, lambda3d, alpha3d, diff3dfunctrain
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
    PI = (pi) (1)
}

PARAMETER {
	 Twait = 1
     PulseInterval = 20
     PulseNum = 10
     
     expfunctrain = 0
     A1 = 5
	 tau1 = 1
	 A2 = .10
	 tau2 = 10
     A3 = .01
	 tau3 = 100
	 
     D = 0.001          
     
     alphafunc = 1 
     aA = .0002
     atau = 148
     
     difffunc = 1
     distance = 1100
     diffCoeff = 0.78
     numMolecules = 4700
     avogadro =  6.023e23
     cleftWidth = 20
     
     diff3dfunc = 1
     diff3dfunctrain = 1
     distance3d = 1100
     diffCoeff3d = 0.78       
     numMolecules3d = 4700
     lambda3d = 1.55
     alpha3d  = 0.21
     
     
}

ASSIGNED {
	T (mM)
    tplexp 
    alpha
    diff
    diff3d
}


INITIAL {
	T = 0 
    tplexp = 0
    alpha = 0
    diff = 0
    diff3d = 0
   
}


BREAKPOINT {
    if (t <= Twait) {         :if the transient hasn't started, apply the ambient concentration
        T = D
    }
    if (t > Twait ) {      
        if (expfunctrain == 1){
            if (PulseNum > 0){
                if (t > (Twait+PulseInterval*0) ) {
                    tplexp = A1*exp(-(1/tau1)*(t - Twait-PulseInterval*0)) + A2*exp(-(1/tau2)*(t - Twait-PulseInterval*0)) + A3*exp(-(1/tau3)*(t - Twait-PulseInterval*0)) 
                }         
            }
            if (PulseNum > 1){
                if (t > (Twait+PulseInterval*1) ) {      
                    tplexp = tplexp + A1*exp(-(1/tau1)*(t - Twait-PulseInterval*1)) + A2*exp(-(1/tau2)*(t - Twait-PulseInterval*1)) + A3*exp(-(1/tau3)*(t - Twait-PulseInterval*1))  
                 }
            }
             if (PulseNum > 2){
                if (t > (Twait+PulseInterval*2) ) {
                    tplexp = tplexp + A1*exp(-(1/tau1)*(t - Twait-PulseInterval*2)) + A2*exp(-(1/tau2)*(t - Twait-PulseInterval*2)) + A3*exp(-(1/tau3)*(t - Twait-PulseInterval*2)) 
                }
            }
            if (PulseNum > 3){
                if (t > (Twait+PulseInterval*3) ) {
                    tplexp = tplexp + A1*exp(-(1/tau1)*(t - Twait-PulseInterval*3)) + A2*exp(-(1/tau2)*(t - Twait-PulseInterval*3)) + A3*exp(-(1/tau3)*(t - Twait-PulseInterval*3)) 
                }
            }
            if (PulseNum > 4){
                if (t > (Twait+PulseInterval*4) ) {
                    tplexp = tplexp + A1*exp(-(1/tau1)*(t - Twait-PulseInterval*4)) + A2*exp(-(1/tau2)*(t - Twait-PulseInterval*4)) + A3*exp(-(1/tau3)*(t - Twait-PulseInterval*4)) 
                }
            }
            if (PulseNum > 5){
                if (t > (Twait+PulseInterval*5) ) {
                    tplexp = tplexp + A1*exp(-(1/tau1)*(t - Twait-PulseInterval*5)) + A2*exp(-(1/tau2)*(t - Twait-PulseInterval*5)) + A3*exp(-(1/tau3)*(t - Twait-PulseInterval*5)) 
                }
            }
            if (PulseNum > 6){
                if (t > (Twait+PulseInterval*6) ) {
                    tplexp = tplexp + A1*exp(-(1/tau1)*(t - Twait-PulseInterval*6)) + A2*exp(-(1/tau2)*(t - Twait-PulseInterval*6)) + A3*exp(-(1/tau3)*(t - Twait-PulseInterval*6)) 
                }
            }
            if (PulseNum > 7){
                if (t > (Twait+PulseInterval*7) ) {
                    tplexp = tplexp + A1*exp(-(1/tau1)*(t - Twait-PulseInterval*7)) + A2*exp(-(1/tau2)*(t - Twait-PulseInterval*7)) + A3*exp(-(1/tau3)*(t - Twait-PulseInterval*7)) 
                }
            }
            if (PulseNum > 8){
                if (t > (Twait+PulseInterval*8) ) {
                    tplexp = tplexp + A1*exp(-(1/tau1)*(t - Twait-PulseInterval*8)) + A2*exp(-(1/tau2)*(t - Twait-PulseInterval*8)) + A3*exp(-(1/tau3)*(t - Twait-PulseInterval*8)) 
                }
            }
            if (PulseNum > 9){
                if (t > (Twait+PulseInterval*9) ) {
                    tplexp = tplexp + A1*exp(-(1/tau1)*(t - Twait-PulseInterval*9)) + A2*exp(-(1/tau2)*(t - Twait-PulseInterval*9)) + A3*exp(-(1/tau3)*(t - Twait-PulseInterval*9)) 
                } 
            }
        }
        if (alphafunc == 1 ) {
            if (PulseNum > 0){       
	           if (t > (Twait+PulseInterval*0) ) {
                    alpha = aA*(t -Twait-PulseInterval*0)/atau*exp(-(t - Twait-PulseInterval*0-atau)/atau) 
               }
            }
            if (PulseNum > 1){       
	           if (t > (Twait+PulseInterval*1) ) {
                    alpha = alpha + aA*(t -Twait-PulseInterval*1)/atau*exp(-(t - Twait-PulseInterval*1-atau)/atau) 
               }
            }
            if (PulseNum > 2){       
	           if (t > (Twait+PulseInterval*2) ) {
                    alpha = alpha + aA*(t -Twait-PulseInterval*2)/atau*exp(-(t - Twait-PulseInterval*2-atau)/atau) 
               }
            }
            if (PulseNum > 3){       
	           if (t > (Twait+PulseInterval*3) ) {
                    alpha = alpha + aA*(t -Twait-PulseInterval*3)/atau*exp(-(t - Twait-PulseInterval*3-atau)/atau) 
               }
            }
            if (PulseNum > 4){       
	           if (t > (Twait+PulseInterval*4) ) {
                    alpha = alpha + aA*(t -Twait-PulseInterval*4)/atau*exp(-(t - Twait-PulseInterval*4-atau)/atau) 
               }
            }
            if (PulseNum > 5){       
	           if (t > (Twait+PulseInterval*5) ) {
                    alpha = alpha + aA*(t -Twait-PulseInterval*5)/atau*exp(-(t - Twait-PulseInterval*5-atau)/atau) 
               }
            }
            if (PulseNum > 6){       
	           if (t > (Twait+PulseInterval*6) ) {
                    alpha = alpha + aA*(t -Twait-PulseInterval*6)/atau*exp(-(t - Twait-PulseInterval*6-atau)/atau) 
               }
            }
            if (PulseNum > 7){       
	           if (t > (Twait+PulseInterval*7) ) {
                    alpha = alpha + aA*(t -Twait-PulseInterval*7)/atau*exp(-(t - Twait-PulseInterval*7-atau)/atau) 
               }
            }
            if (PulseNum > 8){       
	           if (t > (Twait+PulseInterval*8) ) {
                    alpha = alpha + aA*(t -Twait-PulseInterval*8)/atau*exp(-(t - Twait-PulseInterval*8-atau)/atau) 
               }
            }
            if (PulseNum > 9){       
	           if (t > (Twait+PulseInterval*9) ) {
                    alpha = alpha + aA*(t -Twait-PulseInterval*9)/atau*exp(-(t - Twait-PulseInterval*9-atau)/atau) 
               }
            }
        }
        if (difffunc == 1) {
            if (PulseNum > 0){
                if (t > (Twait+PulseInterval*0) ) {
                    diff = ( exp(-((1e-8*distance)^2)/(4*diffCoeff*1e-10*(t - Twait-PulseInterval*0)))*((numMolecules/avogadro)/(4*PI*cleftWidth*1e-8*diffCoeff*1e-10*(t -  Twait-PulseInterval*0))) ) *1e3
                }
            }
            if (PulseNum > 1){
                if (t > (Twait+PulseInterval*1) ) {
                    diff = diff + ( exp(-((1e-8*distance)^2)/(4*diffCoeff*1e-10*(t - Twait-PulseInterval*1)))*((numMolecules/avogadro)/(4*PI*cleftWidth*1e-8*diffCoeff*1e-10*(t -  Twait-PulseInterval*1))) ) *1e3
                }
            }        
            if (PulseNum > 2){
                if (t > (Twait+PulseInterval*2) ) {
                    diff = diff + ( exp(-((1e-8*distance)^2)/(4*diffCoeff*1e-10*(t - Twait-PulseInterval*2)))*((numMolecules/avogadro)/(4*PI*cleftWidth*1e-8*diffCoeff*1e-10*(t -  Twait-PulseInterval*2))) ) *1e3
                }
            }     
            if (PulseNum > 3){
                if (t > (Twait+PulseInterval*3) ) {
                    diff = diff + ( exp(-((1e-8*distance)^2)/(4*diffCoeff*1e-10*(t - Twait-PulseInterval*3)))*((numMolecules/avogadro)/(4*PI*cleftWidth*1e-8*diffCoeff*1e-10*(t -  Twait-PulseInterval*3))) ) *1e3
                }
            }     
            if (PulseNum > 4){
                if (t > (Twait+PulseInterval*4) ) {
                    diff = diff + ( exp(-((1e-8*distance)^2)/(4*diffCoeff*1e-10*(t - Twait-PulseInterval*4)))*((numMolecules/avogadro)/(4*PI*cleftWidth*1e-8*diffCoeff*1e-10*(t -  Twait-PulseInterval*4))) ) *1e3
                }
            } 
            if (PulseNum > 5){
                if (t > (Twait+PulseInterval*5) ) {
                    diff = diff + ( exp(-((1e-8*distance)^2)/(4*diffCoeff*1e-10*(t - Twait-PulseInterval*5)))*((numMolecules/avogadro)/(4*PI*cleftWidth*1e-8*diffCoeff*1e-10*(t -  Twait-PulseInterval*5))) ) *1e3
                }
            } 
            if (PulseNum > 6){
                if (t > (Twait+PulseInterval*6) ) {
                    diff = diff + ( exp(-((1e-8*distance)^2)/(4*diffCoeff*1e-10*(t - Twait-PulseInterval*6)))*((numMolecules/avogadro)/(4*PI*cleftWidth*1e-8*diffCoeff*1e-10*(t -  Twait-PulseInterval*6))) ) *1e3
                }
            } 
            if (PulseNum > 7){
                if (t > (Twait+PulseInterval*7) ) {
                    diff = diff + ( exp(-((1e-8*distance)^2)/(4*diffCoeff*1e-10*(t - Twait-PulseInterval*7)))*((numMolecules/avogadro)/(4*PI*cleftWidth*1e-8*diffCoeff*1e-10*(t -  Twait-PulseInterval*7))) ) *1e3
                }
            } 
            if (PulseNum > 8){
                if (t > (Twait+PulseInterval*8) ) {
                    diff = diff + ( exp(-((1e-8*distance)^2)/(4*diffCoeff*1e-10*(t - Twait-PulseInterval*8)))*((numMolecules/avogadro)/(4*PI*cleftWidth*1e-8*diffCoeff*1e-10*(t -  Twait-PulseInterval*8))) ) *1e3
                }
            } 
            if (PulseNum > 9){
                if (t > (Twait+PulseInterval*9) ) {
                    diff = diff + ( exp(-((1e-8*distance)^2)/(4*diffCoeff*1e-10*(t - Twait-PulseInterval*9)))*((numMolecules/avogadro)/(4*PI*cleftWidth*1e-8*diffCoeff*1e-10*(t -  Twait-PulseInterval*9))) ) *1e3
                }
            }
        }
                                                                    
        if (diff3dfunctrain == 1){
            if (PulseNum > 0){
                if (t > (Twait+PulseInterval*0) ) {
                    diff3d = ((exp(-((1e-8*distance3d)^2)/(4*diffCoeff3d*1e-10*(t - Twait -PulseInterval*0)/(lambda3d^2))))*(numMolecules3d/avogadro)/(8*alpha3d*((PI*diffCoeff3d*1e-10*(t - Twait -PulseInterval*0)/(lambda3d^2))^1.5) ) )*1e3 
                }
            }
            if (PulseNum > 1){
                if (t > (Twait+PulseInterval*1) ) {      
                    diff3d = diff3d + ((exp(-((1e-8*distance3d)^2)/(4*diffCoeff3d*1e-10*(t - Twait -PulseInterval*1)/(lambda3d^2))))*(numMolecules3d/avogadro)/(8*alpha3d*((PI*diffCoeff3d*1e-10*(t - Twait -PulseInterval*1)/(lambda3d^2))^1.5) ) )*1e3  
                }
            }
            if (PulseNum > 2){
                if (t > (Twait+PulseInterval*2) ) {
                    diff3d = diff3d + ((exp(-((1e-8*distance3d)^2)/(4*diffCoeff3d*1e-10*(t - Twait -PulseInterval*2)/(lambda3d^2))))*(numMolecules3d/avogadro)/(8*alpha3d*((PI*diffCoeff3d*1e-10*(t - Twait -PulseInterval*2)/(lambda3d^2))^1.5) ) )*1e3 
                }
            }
            if (PulseNum > 3){
                if (t > (Twait+PulseInterval*3) ) {
                    diff3d = diff3d + ((exp(-((1e-8*distance3d)^2)/(4*diffCoeff3d*1e-10*(t - Twait -PulseInterval*3)/(lambda3d^2))))*(numMolecules3d/avogadro)/(8*alpha3d*((PI*diffCoeff3d*1e-10*(t - Twait -PulseInterval*3)/(lambda3d^2))^1.5) ) )*1e3 
                }
            }
            if (PulseNum > 4){
                if (t > (Twait+PulseInterval*4) ) {
                    diff3d = diff3d + ((exp(-((1e-8*distance3d)^2)/(4*diffCoeff3d*1e-10*(t - Twait -PulseInterval*4)/(lambda3d^2))))*(numMolecules3d/avogadro)/(8*alpha3d*((PI*diffCoeff3d*1e-10*(t - Twait -PulseInterval*4)/(lambda3d^2))^1.5) ) )*1e3 
                }
            }
            if (PulseNum > 5){
                if (t > (Twait+PulseInterval*5) ) {
                    diff3d = diff3d + ((exp(-((1e-8*distance3d)^2)/(4*diffCoeff3d*1e-10*(t - Twait -PulseInterval*5)/(lambda3d^2))))*(numMolecules3d/avogadro)/(8*alpha3d*((PI*diffCoeff3d*1e-10*(t - Twait -PulseInterval*5)/(lambda3d^2))^1.5) ) )*1e3 
                }
            }
            if (PulseNum > 6){
                if (t > (Twait+PulseInterval*6) ) {
                    diff3d = diff3d + ((exp(-((1e-8*distance3d)^2)/(4*diffCoeff3d*1e-10*(t - Twait -PulseInterval*6)/(lambda3d^2))))*(numMolecules3d/avogadro)/(8*alpha3d*((PI*diffCoeff3d*1e-10*(t - Twait -PulseInterval*6)/(lambda3d^2))^1.5) ) )*1e3 
                }
            }
            if (PulseNum > 7){
                if (t > (Twait+PulseInterval*7) ) {
                    diff3d = diff3d + ((exp(-((1e-8*distance3d)^2)/(4*diffCoeff3d*1e-10*(t - Twait -PulseInterval*7)/(lambda3d^2))))*(numMolecules3d/avogadro)/(8*alpha3d*((PI*diffCoeff3d*1e-10*(t - Twait -PulseInterval*7)/(lambda3d^2))^1.5) ) )*1e3 
                }
            }
            if (PulseNum > 8){
                if (t > (Twait+PulseInterval*8) ) {
                    diff3d = diff3d + ((exp(-((1e-8*distance3d)^2)/(4*diffCoeff3d*1e-10*(t - Twait -PulseInterval*8)/(lambda3d^2))))*(numMolecules3d/avogadro)/(8*alpha3d*((PI*diffCoeff3d*1e-10*(t - Twait -PulseInterval*8)/(lambda3d^2))^1.5) ) )*1e3 
                }
            }
             if (PulseNum > 9){
                if (t > (Twait+PulseInterval*9) ) {
                    diff3d = diff3d + ((exp(-((1e-8*distance3d)^2)/(4*diffCoeff3d*1e-10*(t - Twait -PulseInterval*9)/(lambda3d^2))))*(numMolecules3d/avogadro)/(8*alpha3d*((PI*diffCoeff3d*1e-10*(t - Twait -PulseInterval*9)/(lambda3d^2))^1.5) ) )*1e3 
                }
            }
        } 
       
    }
    T=tplexp+alpha+diff+diff3d+D   
}
