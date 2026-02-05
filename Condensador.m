function [varout]=Condensador(Lpc,Wc,Hc,Qcond,Qenf,T4r,T1p,T1,T7,T7p,T8,fm3,fm7,Plow,x4r,wfluid)
T4r=T4r-273.15;
T1p=T1p-273.15;
T1=T1-273.15;
T7=T7-273.15;
T7p=T7p-273.15;
T8=T8-273.15;
fmrf=fm7;
fmwf=fm3;
% Parametros intercambiador
    W=Wc;
    H=Hc;
    A_placa=W*H*1e-6;    
    p_thick=0.6; %mm
    plate_tick=2.3; %Longitud espacio entre placas
    head_tick=11; %Longitud de las placas de entrada y salida    
    Lp=Lpc/1000; %m
    Dport=25.4/1000; %m
% Solo condensacion
if Qenf==0
%Calculo diferencia de temperatura media logaritmica
    LMTD_cond=((T1-T7)-(T1p-T7p))/log((T1-T7)/(T1p-T7p));    
%Calculo de propiedades fluido de trabajo
    k_l_cond=refpropm('L','P',Plow,'Q',0,wfluid); % [W/(m K)]
    rho_l_cond=refpropm('D','P',Plow,'Q',0,wfluid); % [kg/m^3]
    rho_g_cond=refpropm('D','P',Plow,'Q',1,wfluid); % [kg/m^3]    
    xm=x4r/2;
    rho_m_cond=xm*(rho_g_cond-rho_l_cond)+rho_l_cond;
    Pr_l_cond=refpropm('^','P',Plow,'Q',0,wfluid); % [dimmensionless]
    mu_l_cond=refpropm('V','P',Plow,'Q',0,wfluid); %[N*s/m^2]
    i_g_cond=refpropm('H','P',Plow,'Q',1,wfluid); % [J/kg]
    i_f_cond=refpropm('H','P',Plow,'Q',0,wfluid); % [J/kg]
    i_fg_cond=i_g_cond-i_f_cond;
%Calculo de propiedades fluido de enfriamiento
    nu_rf_cond=refpropm('$','T',(T7+T7p)/2+273.15,'P',101.3,'water')/10000; %[m^2/s]
    k_rf_cond=refpropm('L','T',(T7+T7p)/2+273.15,'P',101.3,'water'); % [W/(m K)]
    rho_rf_cond=refpropm('D','T',(T7+T7p)/2+273.15,'P',101.3,'water'); % [kg/m^3]
    Pr_rf_cond=refpropm('^','T',(T7+T7p)/2+273.15,'P',101.3,'water'); % [dimmensionless]
%U asumidos
    Ucond_a=400; %W/m2-K
%Condiciones para iteraci?n
    swich=1; iter=1; tol=1e-2;
%Calculo del U
while (swich==1)    
    %Calculo de ?reas requeridas
        Acond=Qcond/Ucond_a/LMTD_cond; %m2
    %Calculo flujo de calor kW/m2
        qcond=Ucond_a*LMTD_cond;
    %Numero de platos
        Np_cond=floor(Acond/A_placa)+1;    
        Np=Np_cond;
        Nch_pass=floor((Np-1)/2)+1;
        L=(head_tick+plate_tick*Np)/1000;
        pitch=L/Np;
        ch_gap=pitch-p_thick/1000; %m
        Dh=2*ch_gap; %m
    %Refrigeration fluid (water)
    %Velocidades
        Vch_cond_rf=fmrf/(rho_rf_cond*ch_gap*W/1000*Nch_pass); %m/s        
    %Reynolds
        Re_cond_rf=Vch_cond_rf*Dh/nu_rf_cond;        
    %Nusselt
        Nu_cond_rf=0.78*Re_cond_rf^0.5*Pr_rf_cond^(1/3);        
    %h de transferencia de calor
        h_cond_rf=Nu_cond_rf*k_rf_cond/Dh; %W/(m^2-K)        
    %Working fluid
    %Velocidades  
        xm=x4r/2;
        G=fmwf/(ch_gap*W/1000*Nch_pass); %kg/(m^2*s)
        Geq=G*(1-xm+xm*(rho_l_cond/rho_g_cond)^0.5); %kg/(m^2*s)
    %Reynolds
        Re_eq_wf=Geq*Dh/mu_l_cond;
    % Nusselt
    if strcmp(wfluid,'r245fa')
        Re_wf=G*Dh/mu_l_cond;
        if Re_wf<400
            Nu_l=0.1248*Re_wf^0.5*Pr_l_cond^(1/3);
        elseif Re_wf>=400
            Nu_l=0.1248*Re_wf^0.7*Pr_l_cond^(1/3);
        end
        hl=Nu_l*k_l_cond/Dh;
        Co=(rho_g_cond/rho_l_cond)*((1-xm)/xm)^0.8;
        Fr_l=(G/rho_l_cond)^2/(9.81*Dh);
        Bo_wf=qcond/(G*i_fg_cond);
        Nu_cond_wf=45.08*hl*(0.25*Co^-0.45*Fr_l^0.25+75*Bo_wf^0.75);
    else
        Nu_cond_wf=4.118*Re_eq_wf^0.4*Pr_l_cond^0.3;
    end
    %h de transferencia de calor
        h_cond_wf=Nu_cond_wf*k_l_cond/Dh;
    %U
        Ucond_c=1/(1/h_cond_rf+1/h_cond_wf);
    %Verificacion
        if (abs(Ucond_c-Ucond_a)<tol)
            swich=0;
        else
            Ucond_a=Ucond_c;
        end
        iter=iter+1;
        if iter==1000
            swich=0;
        end
end
        Acond=Qcond/Ucond_a/LMTD_cond; %m2
        varout(1)=Acond;
        varout(2)=0;
        varout(3)=Np;
        varout(4)=Ucond_a;
        varout(5)=0;
    % Pressure drop
    % Working fluid side  
        Bo_wf=qcond/(G*i_fg_cond);
        Re_wf=G*Dh/mu_l_cond;
        Pcrit=refpropm('P','C',0,' ',0,wfluid);
        ftp=94.75*Re_eq_wf^(-0.0467)/((Re_wf^0.4)*(Bo_wf^(-0.5))*(Plow/Pcrit)^(-0.8));
        Dpp_cond_wf=2*(G^2)*ftp*(Lp/Dh)*(1/rho_m_cond); %mbar
        Vpt_cond_wf=fmwf/rho_m_cond/(pi/4*Dport^2);
        Dpt_cond_wf=1.3*rho_m_cond/2*Vpt_cond_wf^2;
    % Refrigeration fluid side
    % jf
        jf_cond_rf=0.6/Re_cond_rf^0.3;    
    % Caida de presion flujo
        Dpp_cond_rf=8*jf_cond_rf*(Lp/Dh)*rho_rf_cond/2*Vch_cond_rf^2; %mbar
        Vpt_cond_rf=fmrf/rho_rf_cond/(pi/4*Dport^2);
        Dpt_cond_rf=1.3*rho_rf_cond/2*Vpt_cond_rf^2;
        varout(6)=Dpp_cond_wf;
        varout(7)=0;
        varout(8)=Dpp_cond_rf;
        varout(9)=0;
        ED_cond_DP=fmwf*(Dpp_cond_wf+Dpt_cond_wf)/rho_m_cond+fmrf*(Dpp_cond_rf+Dpt_cond_rf)/rho_rf_cond;
        ED_enf_DP=0;
        
        varout(10)=ED_cond_DP+ED_enf_DP; %W
else
%Calculo diferencia de temperatura media logaritmica
    LMTD_cond=((T1-T7)-(T1p-T7p))/log((T1-T7)/(T1p-T7p));
    LMTD_enf=((T1p-T7p)-(T4r-T8))/log((T1p-T7p)/(T4r-T8));
%Calculo de propiedades fluido de trabajo
    k_wf_enf=refpropm('L','T',(T1p+T4r)/2+273.15,'P',Plow,wfluid); % [W/(m K)]
    rho_wf_enf=refpropm('D','T',(T1p+T4r)/2+273.15,'P',Plow,wfluid); % [kg/m^3]
    Pr_wf_enf=refpropm('^','T',(T1p+T4r)/2+273.15,'P',Plow,wfluid); % [dimmensionless]
    mu_wf_enf=refpropm('V','T',(T1p+T4r)/2+273.15,'P',Plow,wfluid); %[N*s/m^2]        
    k_l_cond=refpropm('L','P',Plow,'Q',0,wfluid); % [W/(m K)]
    rho_l_cond=refpropm('D','P',Plow,'Q',0,wfluid); % [kg/m^3]
    rho_g_cond=refpropm('D','P',Plow,'Q',1,wfluid); % [kg/m^3]
    xm=0.5;
    rho_m_cond=xm*(rho_g_cond-rho_l_cond)+rho_l_cond;
    Pr_l_cond=refpropm('^','P',Plow,'Q',0,wfluid); % [dimmensionless]
    mu_l_cond=refpropm('V','P',Plow,'Q',0,wfluid); %[N*s/m^2]
    i_g_cond=refpropm('H','P',Plow,'Q',1,wfluid); % [J/kg]
    i_f_cond=refpropm('H','P',Plow,'Q',0,wfluid); % [J/kg]
    i_fg_cond=i_g_cond-i_f_cond;
%Calculo de propiedades fluido de enfriamiento
    nu_rf_enf=refpropm('$','T',(T8+T7p)/2+273.15,'P',101.3,'water')/10000; %[m^2/s]
    k_rf_enf=refpropm('L','T',(T8+T7p)/2+273.15,'P',101.3,'water'); % [W/(m K)]
    rho_rf_enf=refpropm('D','T',(T8+T7p)/2+273.15,'P',101.3,'water'); % [kg/m^3]
    Pr_rf_enf=refpropm('^','T',(T8+T7p)/2+273.15,'P',101.3,'water'); % [dimmensionless]    
    nu_rf_cond=refpropm('$','T',(T7+T7p)/2+273.15,'P',101.3,'water')/10000; %[m^2/s]
    k_rf_cond=refpropm('L','T',(T7+T7p)/2+273.15,'P',101.3,'water'); % [W/(m K)]
    rho_rf_cond=refpropm('D','T',(T7+T7p)/2+273.15,'P',101.3,'water'); % [kg/m^3]
    Pr_rf_cond=refpropm('^','T',(T7+T7p)/2+273.15,'P',101.3,'water'); % [dimmensionless]
%U asumidos
    Uenf_a=400; Ucond_a=400; %W/m2-K
%Condiciones para iteraci?n
    swich=1; iter=1; tol=1e-2;
%Calculo del U
while (swich==1)    
    %Calculo de ?reas requeridas
        Acond=Qcond/Ucond_a/LMTD_cond; %m2
        Aenf=Qenf/Uenf_a/LMTD_enf; %m2
    %Calculo flujo de calor kW/m2
        qcond=Ucond_a*LMTD_cond;
    %Numero de platos
        Np_cond=floor(Acond/A_placa)+1;
        Np_enf=floor(Aenf/A_placa)+1;        
        Np=Np_cond+Np_enf;
        Nch_pass=floor((Np-1)/2)+1;
        L=(head_tick+plate_tick*Np)/1000;
        pitch=L/Np;
        ch_gap=pitch-p_thick/1000; %m
        Dh=2*ch_gap; %m
    %Refrigeration fluid (water)
    %Velocidades
        Vch_enf_rf=fmrf/(rho_rf_enf*ch_gap*W/1000*Nch_pass); %m/s
        Vch_cond_rf=fmrf/(rho_rf_cond*ch_gap*W/1000*Nch_pass); %m/s        
    %Reynolds
        Re_enf_rf=Vch_enf_rf*Dh/nu_rf_enf;
        Re_cond_rf=Vch_cond_rf*Dh/nu_rf_cond;        
    %Nusselt
        Nu_enf_rf=0.78*Re_enf_rf^0.5*Pr_rf_enf^(1/3);
        Nu_cond_rf=0.78*Re_cond_rf^0.5*Pr_rf_cond^(1/3);        
    %h de transferencia de calor
        h_enf_rf=Nu_enf_rf*k_rf_enf/Dh; %W/(m^2-K)
        h_cond_rf=Nu_cond_rf*k_rf_cond/Dh; %W/(m^2-K)        
    %Working fluid
    %Velocidades
        xm=0.5;
        Vch_enf_wf=fmwf/(rho_wf_enf*ch_gap*W/1000*Nch_pass); %m/s
        G=fmwf/(ch_gap*W/1000*Nch_pass); %kg/(m^2*s)
        Geq=G*(1-xm+xm*(rho_l_cond/rho_g_cond)^0.5); %kg/(m^2*s)
    %Reynolds
        Re_enf_wf=Vch_enf_wf*Dh*rho_wf_enf/mu_wf_enf;
        Re_eq_wf=Geq*Dh/mu_l_cond;
    % Nusselt    
        Nu_enf_wf=0.78*Re_enf_wf^0.5*Pr_wf_enf^(1/3);        
    if strcmp(wfluid,'r245fa')
        Re_wf=G*Dh/mu_l_cond;
        if Re_wf<400
            Nu_l=0.1248*Re_wf^0.5*Pr_l_cond^(1/3);
        elseif Re_wf>=400
            Nu_l=0.1248*Re_wf^0.7*Pr_l_cond^(1/3);
        end
        hl=Nu_l*k_l_cond/Dh;
        Co=(rho_g_cond/rho_l_cond)*((1-xm)/xm)^0.8;
        Fr_l=(G/rho_l_cond)^2/(9.81*Dh);
        Bo_wf=qcond/(G*i_fg_cond);
        Nu_cond_wf=45.08*hl*(0.25*Co^-0.45*Fr_l^0.25+75*Bo_wf^0.75);
    else
        Nu_cond_wf=4.118*Re_eq_wf^0.4*Pr_l_cond^0.3;
    end
    %h de transferencia de calor
        h_enf_wf=Nu_enf_wf*k_wf_enf/Dh;
        h_cond_wf=Nu_cond_wf*k_l_cond/Dh;
    %U
        Uenf_c=1/(1/h_enf_rf+1/h_enf_wf);
        Ucond_c=1/(1/h_cond_rf+1/h_cond_wf);
    %Verificacion
        if (abs(Uenf_c-Uenf_a)<tol)&&(abs(Ucond_c-Ucond_a)<tol)
            swich=0;
        else
            Uenf_a=Uenf_c;
            Ucond_a=Ucond_c;
        end
        iter=iter+1;
        if iter==1000
            swich=0;
        end
end
        Acond=Qcond/Ucond_a/LMTD_cond; %m2
        Aenf=Qenf/Uenf_a/LMTD_enf; %m2
        varout(1)=Acond;
        varout(2)=Aenf;
        varout(3)=Np;
        varout(4)=Ucond_a;
        varout(5)=Uenf_a;
    % Pressure drop
        Bo_wf=qcond/(G*i_fg_cond);
        Re_wf=G*Dh/mu_l_cond;
        Pcrit=refpropm('P','C',0,' ',0,wfluid);
        ftp=94.75*Re_eq_wf^(-0.0467)/((Re_wf^0.4)*(Bo_wf^(-0.5))*(Plow/Pcrit)^(-0.8));
        Dpp_cond_wf=2*(G^2)*ftp*(Lp/Dh)*(1/rho_m_cond)/100;
        Vpt_cond_wf=fmwf/rho_m_cond/(pi/4*Dport^2);
        Dpt_cond_wf=1.3*rho_m_cond/2*Vpt_cond_wf^2;
    % jf
        jf_cond_rf=0.6/Re_cond_rf^0.3;
        jf_enf_rf=0.6/Re_enf_rf^0.3;
        jf_enf_wf=0.6/Re_enf_wf^0.3;
    % Caida de presion flujo
        Dpp_cond_rf=8*jf_cond_rf*(Lp/Dh)*rho_rf_cond/2*Vch_cond_rf^2; %N/m2
        Vpt_cond_rf=fmrf/rho_rf_cond/(pi/4*Dport^2);
        Dpt_cond_rf=1.3*rho_rf_cond/2*Vpt_cond_rf^2; %N/m2
        
        Dpp_enf_rf=8*jf_enf_rf*(Lp/Dh)*rho_rf_enf/2*Vch_enf_rf^2; %N/m2
        Vpt_enf_rf=fmrf/rho_rf_enf/(pi/4*Dport^2);
        Dpt_enf_rf=1.3*rho_rf_enf/2*Vpt_enf_rf^2;  %N/m2
        
        Dpp_enf_wf=8*jf_enf_wf*(Lp/Dh)*rho_wf_enf/2*Vch_enf_wf^2; %N/m2
        Vpt_enf_wf=fmwf/rho_wf_enf/(pi/4*Dport^2);
        Dpt_enf_wf=1.3*rho_wf_enf/2*Vpt_enf_wf^2;   %N/m2
        
        varout(6)=Dpp_cond_wf;
        varout(7)=Dpp_enf_wf;
        varout(8)=Dpp_cond_rf;
        varout(9)=Dpp_enf_rf;
        
        ED_cond_DP=fmwf*(Dpp_cond_wf+Dpt_cond_wf)/rho_m_cond+fmrf*(Dpp_cond_rf+Dpt_cond_rf)/rho_rf_cond;
        ED_enf_DP=fmwf*(Dpp_enf_wf+Dpt_enf_wf)/rho_wf_enf+fmrf*(Dpp_enf_rf+Dpt_enf_rf)/rho_rf_enf;
        
        varout(10)=ED_cond_DP+ED_enf_DP; %W
end        
end