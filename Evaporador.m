function [varout]=Evaporador(Lp,W,H,T3,T3ls,T2r,T5,T5vs,T5ls,T6,Phigh,fm3,fm5,Heat,fluid,Pco2)
%     [varout,LMTD_eco,LMTD_evap,LMTD_sob]=Evaporador(Lp,W,H,T3,T3ls,T2r,T5,T5vs,T5ls,T6,Phigh,fm3,fm5,Heat,fluid)
    Qeco=Heat(1);
    Qevap=Heat(2);
    Qsob=Heat(3);
%Parametros intercambiador
    A_placa=W*H*1e-6;    
    p_thick=1;%0.6; %mm
    plate_tick=2.3; %Longitud espacio entre placas
    head_tick=11; %Longitud de las placas de entrada y salida
    beta=30;
    R=beta/30;
    F_rf= 0.183*R^2-0.275*R+1.10;
    Lp=Lp/1000; %m Path length=Plate length*Number of passes
    Dport=25.4/1000; %m
%U asumidos
    Ueco_a=400; Uevap_a=400; Usob_a=400;  %W/m2-K
%Condiciones para iteraci?n
    swich=1; iter=1; tol=1e-2;
%Temperaturas       
    Pevap=Phigh; %Pevap*1000
    t3=T3-273.15; %t3
    tevap=T3ls-273.15; %tevap
    t2=T2r-273.15;    %t2
    T2=T5-273.15; %T2
    T2p=T5vs-273.15; %T2p
    T1p=T5ls-273.15; %T1p
    T1=T6-273.15;    %T1
%Flujos masicos
    fmcs=fm3; %0.284; %kg/s
    fmHS=fm5; %0.5  %kg/s
    
%Calculo de propiedades

      v_eco_HS=(refpropm('$','T',(T1p+T1)/2+273.15,'P',Pco2,'co2'))*100;   %[m^2/s] kinematic viscosivity
      v_evap_HS=(refpropm('$','T',(T2p+T1p)/2+273.15,'P',Pco2,'co2'))*100; %[m^2/s] kinematic viscosivity
      v_sob_HS=(refpropm('$','T',(T2+T2p)/2+273.15,'P',Pco2,'co2'))*100;   %[m^2/s] kinematic viscosivity
      
%     v_eco_HS=exp(503.471/((T1p+T1)/2+100)-2.25076)/1000000;   %[m^2/s] 
%     v_evap_HS=exp(503.471/((T2p+T1p)/2+100)-2.25076)/1000000; %[m^2/s]
%     v_sob_HS=exp(503.471/((T2+T2p)/2+100)-2.25076)/1000000;   %[m^2/s]

     rho_eco_HS=(refpropm('D','T',(T1p+T1)/2+273.15,'P',Pco2,'co2'))*100;   %[kg/m^3] Density
     rho_evap_HS=(refpropm('D','T',(T2p+T1p)/2+273.15,'P',Pco2,'co2'))*100; %[kg/m^3] Density
     rho_sob_HS=(refpropm('D','T',(T2+T2p)/2+273.15,'P',Pco2,'co2'))*100;   %[kg/m^3] Density
      


%     rho_eco_HS= 989.086-0.692391*(T1p+T1)/2-0.000300218*((T1p+T1)/2)^2;    %[kg/m^3]
%     rho_evap_HS= 989.086-0.692391*(T2p+T1p)/2-0.000300218*((T2p+T1p)/2)^2; %[kg/m^3]
%     rho_sob_HS= 989.086-0.692391*(T2+T2p)/2-0.000300218*((T2+T2p)/2)^2;    %[kg/m^3]
    
     k_eco_HS=(refpropm('L','T',(T1p+T1)/2+273.15,'P',Pco2,'co2'))*100;    %[W/(m-K)] Thermal conductivity
     k_evap_HS=(refpropm('L','T',(T2p+T1p)/2+273.15,'P',Pco2,'co2'))*100;  %[W/(m-K)] Thermal conductivity
     k_sob_HS=(refpropm('L','T',(T2+T2p)/2+273.15,'P',Pco2,'co2'))*100;    %[W/(m-K)] Thermal conductivity
    
%     k_eco_HS=(0.122684 -(0.0000647634)*(T1p+T1)/2-(0.000000136278)*((T1+T1p)/2)^2);    %[W/(m-K)] Thermal conductivity
%     k_evap_HS=(0.122684 -(0.0000647634)*(T2p+T1p)/2-(0.000000136278)*((T2p+T1p)/2)^2); %[W/(m-K)] Thermal conductivity
%     k_sob_HS=(0.122684 -(0.0000647634)*(T2+T2p)/2-(0.000000136278)*((T2+T2p)/2)^2);    %[W/(m-K)] Thermal conductivity
    
     Cp_eco_HS=(refpropm('C','T',(T1p+T1)/2+273.15,'P',Pco2,'co2'))*100;    %[J/(kg K)] cp
     Cp_evap_HS=(refpropm('C','T',(T2p+T1p)/2+273.15,'P',Pco2,'co2'))*100;  %[[J/(kg K)] cp
     Cp_sob_HS=(refpropm('C','T',(T2+T2p)/2+273.15,'P',Pco2,'co2'))*100;    %[J/(kg K)] cp  

% 
%     Cp_eco_HS=(1.61665 + 0.00318359*(T1+T1p)/2+(0.000000546009)*((T1+T1p)/2)^2)*1000; %[J/(kg-K)]
%     Cp_evap_HS=(1.61665 + 0.00318359*(T2p+T1p)/2+(0.000000546009)*((T2p+T1p)/2)^2)*1000; %[J/(kg-K)]
%     Cp_sob_HS=(1.61665 + 0.00318359*(T2+T2p)/2+(0.000000546009)*((T2+T2p)/2)^2)*1000; %[J/(kg-K)]
%     
    
    % Propiadades working fluid wf
    mu_eco_wf=refpropm('V','T',(t2+tevap)/2+273.15,'P',Pevap,fluid); %[N/(m^2-s)]
    mu_sob_wf=refpropm('V','T',(t3+tevap)/2+273.15,'P',Pevap,fluid); %[N/(m^2-s)]
    rho_eco_wf=refpropm('D','T',(t2+tevap)/2+273.15,'P',Pevap,fluid); %[kg/m^3]
    rho_sob_wf=refpropm('D','T',(t3+tevap)/2+273.15,'P',Pevap,fluid); %[kg/m^3]
    k_eco_wf=refpropm('L','T',(t2+tevap)/2+273.15,'P',Pevap,fluid); %[W/(m-K)]
    k_evap_wf=refpropm('L','T',tevap+273.15,'Q',0,fluid); %[W/(m-K)]
    k_sob_wf=refpropm('L','T',(t3+tevap)/2+273.15,'P',Pevap,fluid); %[W/(m-K)]
    rho_f_wf=refpropm('D','T',tevap+273.15,'Q',0,fluid); %[kg/m^3]
    rho_g_wf=refpropm('D','T',tevap+273.15,'Q',1,fluid); %[kg/m^3]
    mu_f_wf=refpropm('V','T',tevap+273.15,'Q',0,fluid); %[N/(m^2-s)]
    mu_g_wf=refpropm('V','T',tevap+273.15,'Q',1,fluid); %[N/(m^2-s)]
    rho_fg_wf=rho_f_wf-rho_g_wf; %[kg/m^3]
    rho_evap_m=1/(0.5/rho_g_wf+0.5/rho_f_wf); %[kg/m^3]
    mu_evap_m=(rho_evap_m*(0.5*mu_g_wf/rho_g_wf+0.5*mu_f_wf/rho_f_wf)); %[N/(m^2-s)]
    h_fg_wf=refpropm('H','T',tevap+273.15,'Q',1,fluid)-refpropm('H','T',tevap+273.15,'Q',0,fluid); %[J/kg]
    alpha_evap_wf=refpropm('%','T',tevap+273.15,'Q',0,fluid)/10000; %[m^2/s]
    k_l_wf=refpropm('L','T',tevap+273.15,'Q',0,fluid); %[W/(m-K)]
%Calculo diferencia de temperatura media logaritmica
    LMTD_eco=((T1p-tevap)-(T1-t2))/log((T1p-tevap)/(T1-t2));
    LMTD_evap=((T2p-tevap)-(T1p-tevap))/log((T2p-tevap)/(T1p-tevap));
    LMTD_sob=((T2-t3)-(T2p-tevap))/log((T2-t3)/(T2p-tevap));
%Calculo del U
while (swich==1)    
    %Calculo de ?reas requeridas
        Aeco=Qeco/Ueco_a/LMTD_eco;     %m2
        Aevap=Qevap/Uevap_a/LMTD_evap; %m2
        if Qsob==0
            Asob=0;
        else
            Asob=Qsob/Usob_a/LMTD_sob;     %m2
        end
    %Calculo temperatura de pared
%         Ts=(T2p+tevap)/2;
        qevap_asum=Uevap_a*(T1p-tevap);
    %Numero de platos
        Np_eco=floor(Aeco/A_placa)+1;
        Np_evap=floor(Aevap/A_placa)+1;
        if Qsob==0
            Np_sob=0;
        else
            Np_sob=floor(Asob/A_placa)+1;
        end
        Np=Np_eco+Np_evap+Np_sob;
        Nch_pass=floor((Np-1)/2)+1;
        L=(head_tick+plate_tick*Np)/1000;
        pitch=L/Np;
        ch_gap=pitch-p_thick/1000; %m
        Dh=2*ch_gap; %m        
    %Velocidades
        Vch_eco_HS=fmHS/(rho_eco_HS*ch_gap*W/1000*Nch_pass); %m/s
        Vch_evap_HS=fmHS/(rho_evap_HS*ch_gap*W/1000*Nch_pass); %m/s            
        Vch_eco_wf=fmcs/(rho_eco_wf*ch_gap*W/1000*Nch_pass); %m/s
        Vch_evap_wf=fmcs/(rho_evap_m*ch_gap*W/1000*Nch_pass); %m/s
        if Qsob==0
            Vch_sob_HS=0;
            Vch_sob_wf=0;   
        else
            Vch_sob_HS=fmHS/(rho_sob_HS*ch_gap*W/1000*Nch_pass); %m/s
            Vch_sob_wf=fmcs/(rho_sob_wf*ch_gap*W/1000*Nch_pass); %m/s
        end
    %Reynolds
        Re_eco_HS=Vch_eco_HS*Dh/v_eco_HS;
        Re_evap_HS=Vch_evap_HS*Dh/v_evap_HS;            
        Re_eco_wf=Vch_eco_wf*Dh*rho_eco_wf/mu_eco_wf;            
        Re_tp=Vch_evap_wf*Dh*rho_evap_m/mu_evap_m;
        if Qsob==0
            Re_sob_HS=0;
            Re_sob_wf=0;
        else
            Re_sob_HS=Vch_sob_HS*Dh/v_sob_HS;
            Re_sob_wf=Vch_sob_wf*Dh*rho_sob_wf/mu_sob_wf;
        end
    %Prandtl
        Pr_eco_HS=v_eco_HS*Cp_eco_HS*rho_eco_HS/k_eco_HS;
        Pr_evap_HS=v_evap_HS*Cp_evap_HS*rho_evap_HS/k_evap_HS;            
        Pr_eco_wf=refpropm('^','T',(t2+tevap)/2+273.15,'P',Pevap,fluid);
        Pr_evap_wf=refpropm('^','T',tevap+273.15,'Q',0,fluid);            
        if Qsob==0
            Pr_sob_wf=0;
            Pr_sob_HS=0;
        else
            Pr_sob_wf=refpropm('^','T',(t3+tevap)/2+273.15,'P',Pevap,fluid);
            Pr_sob_HS=v_sob_HS*Cp_sob_HS*rho_sob_HS/k_sob_HS;
        end
    %Nusselt
        Nu_eco_HS=0.78*Re_eco_HS^0.5*Pr_eco_HS^(1/3);
        Nu_evap_HS=0.78*Re_evap_HS^0.5*Pr_evap_HS^(1/3);
        if strcmp(fluid,'r245fa')==1
            if Re_eco_wf<400
                Nu_eco_wf=1.29*Re_eco_wf^0.5*Pr_eco_wf^(1/3);
            else
                Nu_eco_wf=1.29*Re_eco_wf^0.7*Pr_eco_wf^(1/3);
            end
        else
            Nu_eco_wf=0.78*Re_eco_wf^0.5*Pr_eco_wf^(1/3);
        end
        if Qsob==0
            Nu_sob_wf=0;
            Nu_sob_HS=0;
        else            
            Nu_sob_HS=0.78*Re_sob_HS^0.5*Pr_sob_HS^(1/3);
            if strcmp(fluid,'r245fa')==1
                if Re_sob_wf<400
                    Nu_sob_wf=1.29*Re_sob_wf^0.5*Pr_sob_wf^(1/3);
                else
                    Nu_sob_wf=1.29*Re_sob_wf^0.7*Pr_sob_wf^(1/3);
                end
            else
                Nu_sob_wf=0.78*Re_sob_wf^0.5*Pr_sob_wf^(1/3);
            end
        end
        if strcmp(fluid,'r245fa')==1
            Gp_evap_wf=fmcs/(ch_gap*W/1000*Nch_pass); %kg/(m2*s)
            Re_f_wf=Gp_evap_wf*Dh/mu_f_wf;
            if Re_f_wf<400
                Nu_l=1.29*Re_f_wf^0.5*Pr_evap_wf^(1/3); %Reynolds liquid phase
            else
                Nu_l=1.29*Re_f_wf^0.7*Pr_evap_wf^(1/3); %Reynolds liquid phase
            end
            hl=Nu_l*k_l_wf/Dh;
            Bo=qevap_asum/(h_fg_wf*Gp_evap_wf);
            h_evap_wf=1.98*hl*Bo^0.5;
        else
            do=0.0146*35*(2*1.52/1000/(9.81*rho_fg_wf))^0.5;
            Nu_evap_wf=(0.00187)*((qevap_asum*do/(k_evap_wf*(tevap+273.15)))^0.56)*(h_fg_wf*do/(alpha_evap_wf^2))^0.31*(Pr_evap_wf^(1/3));
            h_evap_wf=Nu_evap_wf*k_evap_wf/Dh;
        end
    %h de transferencia de calor
        h_eco_HS=Nu_eco_HS*k_eco_HS/Dh; %[W/(m^2-K)]
        h_evap_HS=Nu_evap_HS*k_evap_HS/Dh; %[W/(m^2-K)]
        h_sob_HS=Nu_sob_HS*k_sob_HS/Dh; %[W/(m^2-K)]
        h_eco_wf=Nu_eco_wf*k_eco_wf/Dh;            
        h_sob_wf=Nu_sob_wf*k_sob_wf/Dh;            
    %U
        Ueco_c=1/(1/h_eco_HS+1/h_eco_wf);
        Uevap_c=1/(1/h_evap_HS+1/h_evap_wf);
        if Qsob==0
            Usob_c=0;
        else
            Usob_c=1/(1/h_sob_HS+1/h_sob_wf);                
        end            
    %Verificacion
        if (abs(Ueco_c-Ueco_a)<tol)&&(abs(Uevap_c-Uevap_a)<tol)&&(abs(Usob_c-Usob_a)<tol)
            swich=0;
        else
            Ueco_a=Ueco_c;
            Uevap_a=Uevap_c;
            Usob_a=Usob_c;
        end
        iter=iter+1;
        if iter==1000
            swich=0;
        end
end
    Aeco=Qeco/Ueco_a/LMTD_eco;     %m2
    Aevap=Qevap/Uevap_a/LMTD_evap; %m2
    if Qsob==0
        Asob=0;
    else
        Asob=Qsob/Usob_a/LMTD_sob;     %m2
    end
% Fanning friction factor
    jf_eco_HS=0.6/Re_eco_HS^0.3;
    jf_evap_HS=0.6/Re_evap_HS^0.3;
    jf_eco_wf=0.6/Re_eco_wf^0.3;
    jf_evap_wf=(3.81e4)*F_rf/(Re_tp^0.9*(rho_f_wf/rho_g_wf)^0.16);
    if Qsob==0
        jf_sob_wf=0;
        jf_sob_HS=0;
    else
        jf_sob_wf=0.6/Re_sob_wf^0.3;
        jf_sob_HS=0.6/Re_sob_HS^0.3;
    end            
    %Flow pressure drop
    Dpp_eco_HS=8*jf_eco_HS*(Lp/Dh)*rho_eco_HS/2*Vch_eco_HS^2/100; %mbar
    Dpp_evap_HS=8*jf_evap_HS*(Lp/Dh)*rho_evap_HS/2*Vch_evap_HS^2/100; %mbar
    Dpp_eco_wf=8*jf_eco_wf*(Lp/Dh)*rho_eco_wf/2*Vch_eco_wf^2/100; %mbar
    Dpp_evap_wf=jf_evap_wf*(Lp/Dh)*rho_evap_m/2*Vch_evap_wf^2/100; %mbar
    if Qsob==0
        Dpp_sob_wf=0;
        Dpp_sob_HS=0;
    else
        Dpp_sob_wf=8*jf_sob_wf*(Lp/Dh)*rho_sob_wf/2*Vch_sob_wf^2/100; %mbar
        Dpp_sob_HS=8*jf_sob_HS*(Lp/Dh)*rho_sob_HS/2*Vch_sob_HS^2/100; %mbar
    end
%Port velocity
    Vpt_eco_HS=(fmHS/rho_eco_HS)/(pi/4*Dport^2);
    Vpt_evap_HS=(fmHS/rho_evap_HS)/(pi/4*Dport^2);
    Gpt_evap_wf=fmcs/(pi/4*Dport^2);
    Vpt_eco_wf=(fmcs/rho_eco_wf)/(pi/4*Dport^2);
    if Qsob==0
        Vpt_sob_wf=0;
        Vpt_sob_HS=0;
    else
        Vpt_sob_wf=(fmcs/rho_sob_wf)/(pi/4*Dport^2);
        Vpt_sob_HS=(fmHS/rho_sob_HS)/(pi/4*Dport^2);
    end            
%Port pressure drop
    Dpt_eco_HS=1.3*rho_eco_HS/2*Vpt_eco_HS^2/100; %mbar
    Dpt_evap_HS=1.3*rho_evap_HS/2*Vpt_evap_HS^2/100; %mbar
    Dpt_eco_wf=1.3*rho_eco_wf/2*Vpt_eco_wf^2/100; %mbar
    Dpt_evap_wf=0.75*(Gpt_evap_wf^2/(2*rho_f_wf)+Gpt_evap_wf^2/(2*rho_evap_m))/100; %mbar
    if Qsob==0
        Dpt_sob_wf=0;
        Dpt_sob_HS=0;
    else
        Dpt_sob_wf=1.3*rho_sob_wf/2*Vpt_sob_wf^2/100; %mbar
        Dpt_sob_HS=1.3*rho_sob_HS/2*Vpt_sob_HS^2/100; %mbar
    end            
%Total pressure drop
    P_drop_HS=Dpp_eco_HS+Dpp_evap_HS+Dpp_sob_HS+Dpt_eco_HS+Dpt_evap_HS+Dpt_sob_HS;  %mbar            
    P_drop_wf=Dpp_eco_wf+Dpp_sob_wf+Dpt_eco_wf+Dpt_sob_wf+Dpp_evap_wf+Dpt_evap_wf;  %mbar
%Exergy destroyed by pressure drop
    ED_precal_DP=fm3*(Dpp_eco_wf+Dpt_eco_wf)/rho_eco_wf+fm5*(Dpp_eco_HS+Dpt_eco_HS)/rho_eco_HS;
    ED_evap_DP=fm3*(Dpp_evap_wf+Dpt_evap_wf)/rho_evap_m+fm5*(Dpp_evap_HS+Dpt_evap_HS)/rho_evap_HS;
    if Qsob==0
        ED_sob_DP=0;
    else                
        ED_sob_DP=fm3*(Dpp_sob_wf+Dpt_sob_wf)/rho_sob_wf+fm5*(Dpp_sob_HS+Dpt_sob_HS)/rho_sob_HS;
    end            
%Output
    varout(1)=ED_precal_DP/10; %kW
    varout(2)=ED_evap_DP/10; %kW
    varout(3)=ED_sob_DP/10; %kW
    varout(4)=P_drop_wf; %mbar
    varout(5)=P_drop_HS; %mbar
    varout(6)=Aeco+Aevap+Asob; %m2
    varout(7)=Np;
    varout(8)=Ueco_c; %W/(m^2-K)
    varout(9)=Uevap_c; %W/(m^2-K)
    varout(10)=Usob_c; %W/(m^2-K)
    varout(11)=Re_eco_wf;
    varout(12)=Re_sob_wf;
    varout(13)=Re_tp;
    varout(14)=Re_eco_HS;
    varout(15)=Re_evap_HS;
    varout(16)=Re_sob_HS;
end