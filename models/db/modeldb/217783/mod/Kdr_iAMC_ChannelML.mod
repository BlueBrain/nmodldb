COMMENT

   **************************************************
   File generated by: neuroConstruct v1.7.1 
   **************************************************

   This file holds the implementation in NEURON of the Cell Mechanism:
   Kdr_iAMC_ChannelML (Type: Channel mechanism, Model: ChannelML based process)

   with parameters: 
   /channelml/@units = Physiological Units 
   /channelml/notes = ChannelML file containing a single Channel description 
   /channelml/channel_type/@name = Kdr_iAMC_ChannelML 
   /channelml/channel_type/@density = yes 
   /channelml/channel_type/status/@value = in_progress 
   /channelml/channel_type/status/comment = Delayed rectifier potassium channel from an inrinsically oscillating AOB mitral cell from parameters recorded in the lab of M. Spehr RWTH Aachen 
   /channelml/channel_type/status/contributor/name = Simon O'Connor 
   /channelml/channel_type/notes = Mitral cell K DR channel 
   /channelml/channel_type/neuronDBref/modelName = K channels 
   /channelml/channel_type/neuronDBref/uri = http://senselab.med.yale.edu/neuronDB/channelGene2.aspx#table3 
   /channelml/channel_type/current_voltage_relation/@cond_law = ohmic 
   /channelml/channel_type/current_voltage_relation/@ion = k 
   /channelml/channel_type/current_voltage_relation/@default_gmax = 1.0 
   /channelml/channel_type/current_voltage_relation/@default_erev = -86.5 
   /channelml/channel_type/current_voltage_relation/@charge = 1 
   /channelml/channel_type/current_voltage_relation/gate[1]/@name = m 
   /channelml/channel_type/current_voltage_relation/gate[1]/@instances = 1 
   /channelml/channel_type/current_voltage_relation/gate[1]/closed_state/@id = m0 
   /channelml/channel_type/current_voltage_relation/gate[1]/open_state/@id = m 
   /channelml/channel_type/current_voltage_relation/gate[1]/open_state/@fraction = 1 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@name = tau 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@from = m0 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@to = m 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@expr_form = generic 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@expr = v &lt; -66.67 ? 230.256604897 : 29.156*exp(-(v+55)/5.8842) + 18.394 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@name = inf 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@from = m0 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@to = m 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@expr_form = generic 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@expr = v &lt; -66.67 ? 0 : 0.005*v + 0.35 
   /channelml/channel_type/current_voltage_relation/gate[2]/@name = n 
   /channelml/channel_type/current_voltage_relation/gate[2]/@instances = 1 
   /channelml/channel_type/current_voltage_relation/gate[2]/closed_state/@id = n0 
   /channelml/channel_type/current_voltage_relation/gate[2]/open_state/@id = n 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@name = tau 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@from = n0 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@to = n 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@expr_form = generic 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@expr = 1.5 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@name = inf 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@from = n0 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@to = n 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@expr_form = generic 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@expr = v &lt; 5 ? 0 : ((-0.006*v) + 0.40) 
   /channelml/channel_type/impl_prefs/table_settings/@max_v = 100 
   /channelml/channel_type/impl_prefs/table_settings/@min_v = -100 
   /channelml/channel_type/impl_prefs/table_settings/@table_divisions = 400 

// File from which this was generated: /home/Simon/NML2_Test/AOB_MC_neuroConstruct/cellMechanisms/Kdr_iAMC_ChannelML/KChannel.xml

// XSL file with mapping to simulator: /home/Simon/NML2_Test/AOB_MC_neuroConstruct/cellMechanisms/Kdr_iAMC_ChannelML/ChannelML_v1.8.1_NEURONmod.xsl

ENDCOMMENT


?  This is a NEURON mod file generated from a ChannelML file

?  Unit system of original ChannelML file: Physiological Units

COMMENT
    ChannelML file containing a single Channel description
ENDCOMMENT

TITLE Channel: Kdr_iAMC_ChannelML

COMMENT
    Mitral cell K DR channel
ENDCOMMENT


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (um) = (micrometer)
    (molar) = (1/liter)
    (mM) = (millimolar)
    (l) = (liter)
}


    
NEURON {
      

    SUFFIX Kdr_iAMC_ChannelML
    USEION k READ ek WRITE ik VALENCE 1 ? reversal potential of ion is read, outgoing current is written
           
        
    RANGE gmax, gion
    
    RANGE minf, mtau
    
    RANGE ninf, ntau
    
}

PARAMETER { 
      

    gmax = 0.001 (S/cm2)  ? default value, should be overwritten when conductance placed on cell
    
}



ASSIGNED {
      

    v (mV)
    
    celsius (degC)
          

    ? Reversal potential of k
    ek (mV)
    ? The outward flow of ion: k calculated by rate equations...
    ik (mA/cm2)
    
    
    gion (S/cm2)
    minf
    mtau (ms)
    ninf
    ntau (ms)
    
}

BREAKPOINT { 
                        
    SOLVE states METHOD cnexp
        
    gion = gmax * ((1*m)
^1) * (n
^1)      

    ik = gion*(v - ek)
            

}



INITIAL {
    
    ek = -86.5
        
    rates(v)
    m = minf
        n = ninf
        
    
}
    
STATE {
    m
    n
    
}



DERIVATIVE states {
    rates(v)
    m' = (minf - m)/mtau
            n' = (ninf - n)/ntau
            

}

PROCEDURE rates(v(mV)) {  
    
    ? Note: not all of these may be used, depending on the form of rate equations
    LOCAL  alpha, beta, tau, inf, gamma, zeta
, temp_adj_m
, temp_adj_n
    
    TABLE minf, mtau,ninf, ntau
 DEPEND celsius FROM -100 TO 100 WITH 400
    
    UNITSOFF
    temp_adj_m = 1
    temp_adj_n = 1
    
            
                
           

        
    ?      ***  Adding rate equations for gate: m  ***
         
    ? Found a generic form of the rate equation for tau, using expression: v < -66.67 ? 230.256604897 : 29.156*exp(-(v+55)/5.8842) + 18.394
    
    
    if (v < -66.67 ) {
        tau =  230.256604897 
    } else {
        tau =  29.156*exp(-(v+55)/5.8842) + 18.394
    }
    mtau = tau/temp_adj_m
     
    ? Found a generic form of the rate equation for inf, using expression: v < -66.67 ? 0 : 0.005*v + 0.35
    
    
    if (v < -66.67 ) {
        inf =  0 
    } else {
        inf =  0.005*v + 0.35
    }
    minf = inf
    


    ?     *** Finished rate equations for gate: m ***
    

    
            
                
           

        
    ?      ***  Adding rate equations for gate: n  ***
         
    ? Found a generic form of the rate equation for tau, using expression: 1.5
    tau = 1.5
        
    ntau = tau/temp_adj_n
     
    ? Found a generic form of the rate equation for inf, using expression: v < 5 ? 0 : ((-0.006*v) + 0.40)
    
    
    if (v < 5 ) {
        inf =  0 
    } else {
        inf =  ((-0.006*v) + 0.40)
    }
    ninf = inf
    


    ?     *** Finished rate equations for gate: n ***
    

         

}


UNITSON


