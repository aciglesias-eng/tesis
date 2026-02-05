function [varout]=Cooler(Lpr,Wr,Hr,T2r,T2,T4,T4r,Phigh,Plow,Qreg,fluid,fm_h, fm_c)        
    %Heat exchanger dimmensions
        W=Wr; H=Hr;
        A_placa=W*H*1e-6;    
        p_thick=0.6;%0.6; %mm
        plate_tick=2.3; %Plate lenght
        head_tick=11; %Head Plate lenght
        Lp=Lpr/1000; %m Path length=Plate length*Number of passes
        Dport=25.4/1000; %m
    %U assumed
        Ureg_a=400;
    %Parameters
        swich=1; iter=1; tol=1e-2;
    %Properties
        Tfilm_c = 0.5*(T4r+T4);

        if ~isfinite(Tfilm_c) || ~isfinite(Plow)
            error('Cooler: Tfilm_c o Plow no finitos. T4=%g T4r=%g Plow=%g', T4, T4r, Plow);
        end

        if Tfilm_c < 150 || Tfilm_c > 2000
            error('Cooler: Tfilm_c fuera de rango. T4=%g T4r=%g Tfilm_c=%g K', T4, T4r, Tfilm_c);
        end

        if Plow <= 0 || Plow > 1e6
            error('Cooler: Plow sospechosa (kPa). Plow=%g', Plow);
        end

        mu_hf=refpropm('V','T',0.5*(T2r+T2),'P',Phigh,'air.ppf');   %[N/(m^2-s)]
        mu_cf=refpropm('V','T',0.5*(T4r+T4),'P',Plow,fluid);        %[N/(m^2-s)]
        rho_hf=refpropm('D','T',0.5*(T2r+T2),'P',Phigh,'air.ppf');  %[kg/m^3]
        rho_cf=refpropm('D','T',0.5*(T4r+T4),'P',Plow,fluid);       %[kg/m^3]
        k_hf=refpropm('L','T',0.5*(T2r+T2),'P',Phigh,'air.ppf');    %[W/(m-K)]
        k_cf=refpropm('L','T',0.5*(T4r+T4),'P',Plow,fluid);         %[W/(m-K)]
    %Log mean temperature difference
        LMTD=((T4-T2r)-(T4r-T2))/log((T4-T2r)/(T4r-T2));        
    %Calculo del U
    while (swich==1)           
        %Required heat transfer area
            Areg=Qreg/Ureg_a/LMTD;     %m2        
        %Number of plates
            Np=floor(Areg/A_placa)+1;            
            Nch_pass=floor((Np-1)/2)+1;
            L=(head_tick+plate_tick*Np)/1000;
            pitch=L/Np;
            ch_gap=pitch-p_thick/1000; %m
            Dh=2*ch_gap; %m        
        %Channel velocity
            Vch_hf=fm_h/(rho_hf*ch_gap*W/1000*Nch_pass); %m/s
            Vch_cf=fm_c/(rho_cf*ch_gap*W/1000*Nch_pass); %m/s
        %Reynolds
            Re_hf=Vch_hf*Dh*rho_hf/mu_hf;
            Re_cf=Vch_cf*Dh*rho_cf/mu_cf;
        %Prantdl
            Tfilm_h = 0.5*(T2r+T2);
            Tfilm_c = 0.5*(T4r+T4);
            Cp_hf = refpropm('C','T',Tfilm_h,'P',Phigh,'air.ppf');  % J/kg-K
            Cp_cf = refpropm('C','T',Tfilm_c,'P',Plow,fluid);      % J/kg-K
            Pr_hf = Cp_hf*mu_hf/k_hf;
            Pr_cf = Cp_cf*mu_cf/k_cf;

        %Nusselt    
            Nu_hf=0.78*Re_hf^0.5*Pr_hf^(1/3);
            Nu_cf=0.78*Re_cf^0.5*Pr_cf^(1/3);
        %h
            h_hf=Nu_hf*k_hf/Dh;
            h_cf=Nu_cf*k_cf/Dh;
        %U
            Ureg_c=1/(1/h_hf+1/h_cf);
        %Verifying U
            if (abs(Ureg_c-Ureg_a)<tol)
                swich=0;
            else
                Ureg_a=Ureg_c;                
            end
            iter=iter+1;
            if iter==1000
                swich=0;
            end
    end
            Areg=Qreg/Ureg_a/LMTD;     %m2
        % Friction factor
            jf_hf=0.6/Re_hf^0.3;
            jf_cf=0.6/Re_cf^0.3;
        %Channel pressure drop
            Dpp_hf=8*jf_hf*(Lp/Dh)*rho_hf/2*Vch_hf^2; %N/m2
            Dpp_cf=8*jf_cf*(Lp/Dh)*rho_cf/2*Vch_cf^2; %N/m2
        %Port pressure drop
            Vpt_hf=(fm_h/rho_hf)/(pi/4*Dport^2); %m/2
            Vpt_cf=(fm_c/rho_cf)/(pi/4*Dport^2); %m/2        
            Dpt_hf=1.3*rho_hf/2*Vpt_hf^2; %N/m2
            Dpt_cf=1.3*rho_cf/2*Vpt_cf^2; %N/m2
        %Exergy destroyed because of pressure drop
            ExDP_hf=fm_h*(Dpp_hf+Dpt_hf)/rho_hf; %W
            ExDP_cf=fm_c*(Dpp_cf+Dpt_cf)/rho_cf; %W
        %Output
            varout(1)=ExDP_hf+ExDP_cf; %W
            varout(2)=Areg;
            varout(3)=Np;
            varout(4)=Ureg_c; %W/(m^2-K)
            
end