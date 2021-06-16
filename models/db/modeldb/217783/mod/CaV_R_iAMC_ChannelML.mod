COMMENT

   **************************************************
   File generated by: neuroConstruct v1.7.1 
   **************************************************

   This file holds the implementation in NEURON of the Cell Mechanism:
   CaV_R_iAMC_ChannelML (Type: Channel mechanism, Model: ChannelML based process)

   with parameters: 
   /channelml/@units = Physiological Units 
   /channelml/notes = AOB Mitral Cell R-Type Calcium Channel 
   /channelml/channel_type/@name = CaV_R_iAMC_ChannelML 
   /channelml/channel_type/status/@value = stable 
   /channelml/channel_type/status/comment = R-Type Calcium Channel from an intrinsically oscillating AOB mitral cell from parameters recorded in the lab of M. Spehr RWTH Aachen 
   /channelml/channel_type/status/contributor/name = Simon O'Connor 
   /channelml/channel_type/notes = A High Voltage Activated Ca2+ channel 
   /channelml/channel_type/authorList/modelTranslator/name = Simon O'Connor 
   /channelml/channel_type/authorList/modelTranslator/institution = UH 
   /channelml/channel_type/authorList/modelTranslator/email = simon.oconnor - at - btinternet.com 
   /channelml/channel_type/neuronDBref/modelName = Calcium channels 
   /channelml/channel_type/neuronDBref/uri = http://senselab.med.yale.edu/neurondb/channelGene2.aspx#table1 
   /channelml/channel_type/current_voltage_relation/@cond_law = ohmic 
   /channelml/channel_type/current_voltage_relation/@ion = ca 
   /channelml/channel_type/current_voltage_relation/@default_gmax = 0.15 
   /channelml/channel_type/current_voltage_relation/@default_erev = 80 
   /channelml/channel_type/current_voltage_relation/@charge = 2 
   /channelml/channel_type/current_voltage_relation/@fixed_erev = yes 
   /channelml/channel_type/current_voltage_relation/gate[1]/@name = m 
   /channelml/channel_type/current_voltage_relation/gate[1]/@instances = 2 
   /channelml/channel_type/current_voltage_relation/gate[1]/closed_state/@id = m0 
   /channelml/channel_type/current_voltage_relation/gate[1]/open_state/@id = m 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@name = tau 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@from = m0 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@to = m 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@expr_form = generic 
   /channelml/channel_type/current_voltage_relation/gate[1]/time_course/@expr = v &lt; -30 ? 28.4118 : 3.1738 + (25.238 * (exp(-1 * ((v + 30)/17.498)))) 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@name = inf 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@from = m0 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@to = m 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@expr_form = sigmoid 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@rate = 1 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@scale = -2.0914 
   /channelml/channel_type/current_voltage_relation/gate[1]/steady_state/@midpoint = -38.037 
   /channelml/channel_type/current_voltage_relation/gate[2]/@name = h 
   /channelml/channel_type/current_voltage_relation/gate[2]/@instances = 1 
   /channelml/channel_type/current_voltage_relation/gate[2]/closed_state/@id = h0 
   /channelml/channel_type/current_voltage_relation/gate[2]/open_state/@id = h 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@name = tau 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@from = h0 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@to = h 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@expr_form = generic 
   /channelml/channel_type/current_voltage_relation/gate[2]/time_course/@expr = v &lt; -30 ? 21.0638148543 : 10.8 + (3.0 * (exp(-1 * ((v+20)/8.13)))) 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@name = inf 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@from = h0 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@to = h 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@expr_form = generic 
   /channelml/channel_type/current_voltage_relation/gate[2]/steady_state/@expr = ((1/(1+(exp(-1 * (v-(-38.037))/-2.0914)))) + (0.6928/(1+(exp(-1 * (v-(-38.037))/2.0914))))) 

// File from which this was generated: /home/Simon/NML2_Test/AOB_MC_neuroConstruct/cellMechanisms/CaV_R_iAMC_ChannelML/CaHVA_Chan.xml

// XSL file with mapping to simulator: /home/Simon/NML2_Test/AOB_MC_neuroConstruct/cellMechanisms/CaV_R_iAMC_ChannelML/ChannelML_v1.8.1_NEURONmod.xsl

ENDCOMMENT


?  This is a NEURON mod file generated from a ChannelML file

?  Unit system of original ChannelML file: Physiological Units

COMMENT
    AOB Mitral Cell R-Type Calcium Channel
ENDCOMMENT

TITLE Channel: CaV_R_iAMC_ChannelML

COMMENT
    A High Voltage Activated Ca2+ channel
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
      

    SUFFIX CaV_R_iAMC_ChannelML
    USEION ca WRITE ica VALENCE 2 ?  outgoing current is written
           
        
    RANGE gmax, gion
    
    RANGE minf, mtau
    
    RANGE hinf, htau
    
}

PARAMETER { 
      

    gmax = 0.00015 (S/cm2)  ? default value, should be overwritten when conductance placed on cell
    
}



ASSIGNED {
      

    v (mV)
    
    celsius (degC)
          

    ? Reversal potential of ca
    eca (mV)
    ? The outward flow of ion: ca calculated by rate equations...
    ica (mA/cm2)
    
    
    gion (S/cm2)
    minf
    mtau (ms)
    hinf
    htau (ms)
    
}

BREAKPOINT { 
                        
    SOLVE states METHOD cnexp
        
    gion = gmax * (m
^2) * (h
^1)      

    ica = gion*(v - eca)
            

}



INITIAL {
    
    eca = 80
        
    rates(v)
    m = minf
        h = hinf
        
    
}
    
STATE {
    m
    h
    
}



DERIVATIVE states {
    rates(v)
    m' = (minf - m)/mtau
            h' = (hinf - h)/htau
            

}

PROCEDURE rates(v(mV)) {  
    
    ? Note: not all of these may be used, depending on the form of rate equations
    LOCAL  alpha, beta, tau, inf, gamma, zeta
, temp_adj_m,
         A_inf_m, B_inf_m, Vhalf_inf_m
, temp_adj_h
    
    TABLE minf, mtau,hinf, htau
 DEPEND celsius FROM -100 TO 100 WITH 400
    
    UNITSOFF
    temp_adj_m = 1
    temp_adj_h = 1
    
            
                
           

        
    ?      ***  Adding rate equations for gate: m  ***
         
    ? Found a generic form of the rate equation for tau, using expression: v < -30 ? 28.4118 : 3.1738 + (25.238 * (exp(-1 * ((v + 30)/17.498))))
    
    
    if (v < -30 ) {
        tau =  28.4118 
    } else {
        tau =  3.1738 + (25.238 * (exp(-1 * ((v + 30)/17.498))))
    }
    mtau = tau/temp_adj_m
    
    ? Found a parameterised form of rate equation for inf, using expression: A / (1 + exp((v-Vhalf)/B))
    A_inf_m = 1
    B_inf_m = -2.0914
    Vhalf_inf_m = -38.037 
    inf = A_inf_m / (exp((v - Vhalf_inf_m) / B_inf_m) + 1)
    
    minf = inf
    


    ?     *** Finished rate equations for gate: m ***
    

    
            
                
           

        
    ?      ***  Adding rate equations for gate: h  ***
         
    ? Found a generic form of the rate equation for tau, using expression: v < -30 ? 21.0638148543 : 10.8 + (3.0 * (exp(-1 * ((v+20)/8.13))))
    
    
    if (v < -30 ) {
        tau =  21.0638148543 
    } else {
        tau =  10.8 + (3.0 * (exp(-1 * ((v+20)/8.13))))
    }
    htau = tau/temp_adj_h
     
    ? Found a generic form of the rate equation for inf, using expression: ((1/(1+(exp(-1 * (v-(-38.037))/-2.0914)))) + (0.6928/(1+(exp(-1 * (v-(-38.037))/2.0914)))))
    inf = ((1/(1+(exp(-1 * (v-(-38.037))/-2.0914)))) + (0.6928/(1+(exp(-1 * (v-(-38.037))/2.0914)))))
        
    hinf = inf
    


    ?     *** Finished rate equations for gate: h ***
    

         

}


UNITSON


