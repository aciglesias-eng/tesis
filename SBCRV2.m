function[LCOH, c_el_bus,LCOEn, WNETO, W_Rankine, n_ciclo, n_ex, n_ex_2, n_th, n_en_PEM,n_ex_PEM, T20, m_H2, m_H2O_in, XD, In, Ytotal, Ytotal_PEM]=SBCRV2(nt,nc,Thigh,rc,pmaire)
format long

%ENTRADAS A LA FUNCIÓN-------------------------
% nt            = Eficiencia turbinas [0.7-0.9];
% nc            = Eficiencia comprensores [0.7-0.9];
% Thigh         = Temperatura de alta del ciclo Brayton dado por el recibidor de la torre solar [650 -800°C];
% rc            = Relacion de presion del Brayton
% pmaire        = porcentaje de flujo másico de aire que entra al cooler

% SALIDAS DE LA FUNCIÓN-------------------------
% W_neto        = Trabajo neto del sistema global (Brayton + ORC)
% W_Rankine     = Trabajo neto del sistema global (ORC)
% n_ciclo       = Eficiencia térmica total de los ciclos (Brayton + ORC)
% n_ex          = Eficiencia exergética global del sistema
%                   (1 – Exergía destruida total / exergía de entrada solar)
% n_th          = Eficiencia térmica total de los ciclos (Brayton + ORC + PEM)
% n_en_PEM      = Eficiencia energética total de la PEM
% n_ex_PEM      = Eficiencia exergética total de la PEM
% T20           = Temperatura de entrada de la PEM
% m_H2          = Flujo másico de hidrógeno producido
% m_H2O_in      = Flujo másico de agua que requiere la PEM
% XD            = Vector de exergías destruidas y aportadas por componente
% In            = Vector de indicadores de sostenibilidad exergética:
%                   [EWR (Exergy Waste Ratio), EEF (Environmental Effect Factor), 
%                   ESI (Exergetic Sustainability Index)]
% Ytotal        = Vector de impacto ambiental total asociado al sistema


%Entradas SISTEMA BRAYTON----------------------------------------------------------------------------------------------------------------------
    fluido = 'co2';          % Fluido de trabajo Brayton 
    T0 = 25 + 273.15;        % Temperatura ambiente (K)
    P0 = 101.325;            % Presión atmosférica (kPa) 
    Tlow = 55 + 273.15;      % Temp low (K)--> Guajira
    m_CO2 = 1;              % Flujo másico (kg/s)
    eff = 0.95;              % Eficacia de los recuperadores (aprox para el punto 6)
    Nr = 40;                 % Nr = Número de segmentos para discretizar los intercambiadores
    r1 = rc;                 % Relación de presión1 --> r1 = Phigh/Plow (turbinas completas)
    r2 = 1.43;               % Relación de presión2 --> r2 = Phigh/Pintermedia (Turbina 2)
    Phigh = 25000;           % Presion alta (Salida compresores)
    Plow = Phigh/r1;         % Presión después de salir de la turbina 2
    Prr = Phigh/r2;          % Presión después de salir de la turbina 1
    Tr = Thigh + 150;        % Temperatura del recibidor
    I = 914;                 % Radiacion solar = 950 W/m2
    tao = 8000;              % 
    interes = 0.12;
    years = 20;
    gamma = 0.06;
    HX = [460,160,336,125,5; % Librería de datos de intercambiadores de calor
        799,160,675,150,12;
        837,310,590,200,28;
        1066,310,819,200,40;
        735,310,494,200,20;
        940,310,694,200,30;
        1135,310,894,200,45;
        1080,440,650,500,100;
        1160,480,719,500,105;
        1332,480,894,500,150;
        1579,480,1141,500,200;
        1826,480,1388,500,250;
        1470,610,941.4,700,280;
        1835,610,1306.2,700,420;
        2200,610,1671,700,560;
        1380,760,770,700,300;
        1740,760,1130,700,455;
        2100,760,1490,700,700;
        2460,760,1850,700,910;
        1930,980,1100,700,585;
        2320,980,1490,700,875;
        2710,980,1879,700,1120;
        3100,980,2267,700,1330;
        2855,1370,1822,700,1400;
        3211,1370,2178,700,1750;
        3567,1370,2534,700,2100];
    
% Estado de referencia para exergía (estado muerto)
    [h0, s0] = refpropm('HS','T',T0,'P',P0,fluido);

% Puntos --------------------------------------------------------------------------------------------------------------------------------------------------
    % Punto 1
        T1 = Thigh; P1 = Phigh;
        [h1, s1] = refpropm('HS','T',T1,'P',P1,fluido);
        x1 = (h1-h0)-T0*(s1-s0);

    % Punto 2
        P2 = Prr;
        [h2, T2, s2] = turbina(T1, P1, nt, P2, fluido);
        x2 = (h2-h0)-T0*(s2-s0);

    % Punto 3
        T3 = T1; P3 = P2;   %Se asume que el recalentador sube la temperatura de punto 3 igual a punto 1 1
        [h3,s3] = refpropm('HS','T',T3,'P',P3,fluido);
        x3 = (h3-h0)-T0*(s3-s0);
        
    % Punto 4
        P4 = Plow;
        [h4,T4,s4] = turbina(T3,P3,nt,P4,fluido);
        x4 = (h4-h0)-T0*(s4-s0);
      
    % Punto 7
        T7 = Tlow; P7 = P4;
        [h7,s7] = refpropm('HS','T',T7,'P',P7,fluido);
        x7 = (h7-h0)-T0*(s7-s0);

    % Punto 8
        P8 = Phigh;
        [h8,T8,s8] = compresor(T7,P7,nc,P8,fluido);
        x8 = (h8-h0)-T0*(s8-s0);
        
    % Punto 6
        P6 = P4;
        h6_p = refpropm('H','T',T8,'P',P6,fluido);%Por caso ideal T8=T6
        h6 = h4-(eff*(h4-h6_p)); 
        [T6,s6] = refpropm('TS','P',P6,'h',h6,fluido);
        x6 = (h6-h0)-T0*(s6-s0);
        
    % Punto 10
        P10 = Phigh;
        [h10,T10,s10] = compresor(T6,P6,nc,P10,fluido); 
        x10 = (h10-h0)-T0*(s10-s0);
        
    % Punto 11
        T11 = T10; P11 = P10;                                                                   
        h11 = h10; s11 = s10;
        x11 = x10;

    % Punto 5
        [~,T5_htr,T12_htr,~,Q_htr] = intercambiador_V1(T4,P4,T11,P11,fluido,eff,Nr,m_CO2,m_CO2);
                                    %Tmin: Temperatura minima HTR 
                                    %T5_htr:Tlinea caliente 
                                    %T12_htr: Tlinea fria
                                    %UA_htr: Conductancia global HTR
                                    %Q_htr: Calor recuperado HTR
        T5 = T5_htr; P5 = P4;
        [h5,s5] = refpropm('HS','T',T5,'P',P5,fluido);
        x5 = (h5-h0)-T0*(s5-s0);
        
    % Punto 6AI
        T6AI = T6; P6AI = P6;
        h6AI = h6; s6AI = s6;
        x6AI = x6;
                     
    % Punto 9
        T9 = T10; P9 = P10;                                                                   
        h9 = h10; s9 = s10;
        x9 = x10;       
        
% --------------------------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------
        % Flujo másico (línea 6AI)
            splitA = (h5-h6)/(h9-h8); %split = m6AI/mtotal= fraccion masica de linea c1
            
        % Flujo másico (línea 6BI)
            splitB = 1-splitA ;
            
% --------------------------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------         
    
% Punto 6AO1/ACOPLE CON ORC
        P6AO1 = P6AI;
        [T6AO3,A_ORC,W_t3,W_b2,Xd_orc,TPhsmx_orc,Q_cond_vec,Q_evap_vec]=ORCV2(T6AI,P6AO1,m_CO2*splitA,'RORC', 0.85);
        T6AO1 = TPhsmx_orc(12,1);
        h6AO1 = TPhsmx_orc(12,3); s6AO1 = TPhsmx_orc(12,4);
        x6AO1 = TPhsmx_orc(12,6);
        
     % Punto 6AO2/ACOPLE CON ORC
        T6AO2 = TPhsmx_orc(13,1);P6AO2 = P6AI;
        h6AO2 = TPhsmx_orc(13,3); s6AO2 = TPhsmx_orc(13,4);
        x6AO2 = TPhsmx_orc(13,6);
        
     % Punto 6AO3/ACOPLE CON ORC
        P6AO3 = P6AI;
        h6AO3 = TPhsmx_orc(14,3); s6AO3 = TPhsmx_orc(14,4);
        x6AO3 = TPhsmx_orc(14,6);      
        
    % Punto 6BI
        T6BI = T6; P6BI = P6;
        h6BI = h6; s6BI = s6;
        x6BI = x6;
        
    % Punto 12
        T12 = T12_htr;
        P12 = Phigh;
        [h12,s12] = refpropm('HS','T',T12,'P',P12,fluido);
        x12 = (h12-h0)-T0*(s12-s0);
        
% --------------------------------------------------------------------------------------------------------------------------------------------------
% CALCULO ENERGÉTICO Y EXERGÉTICO ---------------------------------------------------------------------    
% HEATER AND REHEATER (Puntos: 12, 1, 2, 3)
    % HEATER (Puntos: 12, 1)
        Q_heater = m_CO2*(h1-h12);

    % REHEATER (Puntos: 2, 3)
        Q_reheater = m_CO2*(h3-h2); 

    % Calor total
        Q_tt = Q_heater+Q_reheater; % Q_tt = calor total aportado al Brayton: heater + reheater

%---------------------------------------------------------------------                                    
% ANÁLISIS DE LA TORRE SOLAR (Puntos: 12, 1, 2, 3)
    [nsr, Fi, nfield] = solarReceiver(Thigh, I);  % nsr: eficiencia térmica del recibidor solar
                                                  % Fi: fracción exergética de la radiación solar (Petela)
                                                  % nfield: eficiencia del campo heliostato
    % Vector con parámetros solares para exportar/reportar
        vector_field = [nsr, Fi, nfield]; 

    % Exergía solar requerida para entregar q_tt al ciclo
        X_solar = Q_tt*Fi/(nsr*nfield);

    % Exergía perdida en el recibidor por pérdidas térmicas (Carnot)
        Xloss_sr = Q_tt*(0.95/nsr-1)*(1-T0/Tr); 

    % Exergía perdida en el campo solar (pérdidas ópticas del campo)
        Xloss_field = Q_tt*Fi*(1-nfield)/(nsr*nfield);
        
    % Exergía destruida en recibidor + heater/reheater (balance exergético)
        Xd_sr = (X_solar*nfield-Xloss_sr+m_CO2*(x12-x1+x2-x3));
     
%---------------------------------------------------------------------
    % ANÁLISIS DE LA TURBINA 1 (Puntos: 1, 2)
        % Potencia de la turbina 1
            W_t1 = m_CO2*(h1-h2);

        % Exergía destruida en la turbina 1
            Xd_t1 = -W_t1+m_CO2*(x1-x2);
 
%---------------------------------------------------------------------
    % ANÁLISIS DE LA TURBINA 2 (Puntos: 3, 4)
        % Potencia de la turbina 2
            W_t2 = m_CO2*(h3-h4);

        % Exergía destruida en la turbina 2
            Xd_t2 = -W_t2+m_CO2*(x3-x4);
            
%---------------------------------------------------------------------
    % ANÁLISIS DEL COMPRESOR 1 (Puntos: 7, 8)
        % Potencia del compresor 1
            W_c1 = m_CO2*splitA*(h7-h8);

        % Exergía destruida en el compresor 1
            Xd_c1 = -W_c1+m_CO2*splitA*(x7-x8);
            
%---------------------------------------------------------------------
    % ANÁLISIS DEL COMPRESOR 2 (Puntos: 6BI, 10)
        % Potencia del compresor 2
            W_c2 = m_CO2*splitB*(h6BI-h10);

        % Exergía destruida en el compresor 2
            Xd_c2 = -W_c2+m_CO2*splitB*(x6BI-x10);
            
%---------------------------------------------------------------------
    % ANÁLISIS DEL INTERCAMBIADOR HTR. ---> High recuperator (Puntos: 4, 5, 11, 12)                                            
        % Área HTR
            for ii=1:26
                Hr=HX(ii,1); Wr=HX(ii,2); Lpr=HX(ii,3);
                [varout]=HTR(Lpr,Wr,Hr,T12,T11,T4,T5,Phigh,Phigh,m_CO2,Q_htr,fluido);
                A_HTR=varout(2);
                if ii==1
                    Amin_HTR=A_HTR;
                else
                    if (A_HTR<Amin_HTR)&&(varout(3)<HX(ii,4))&&(A_HTR<HX(ii,5))
                        Amin_HTR=A_HTR;
                    end
                end
            end
    
        % Exergía destruida en HTR
            Xd_htr = m_CO2*((x4-x5)+(x11-x12));
            
%---------------------------------------------------------------------
    % ANÁLISIS DEL INTERCAMBIADOR LTR. ---> Low recuperator (Puntos: 5, 6, 8, 9) 
        % Calculo de efectividad
            h6_max = refpropm('H','T',T8,'P',P6,fluido);
            h9_max = refpropm('H','T',T5,'P',P9,fluido);
            eff_ltr = (h5-h6)/ (min((h5-h6_max),splitA*(h9_max-h8)));
        % Propiedades LTR
            [~,~,~,~,Q_ltr] = ...
            intercambiador_V1(T5,P5,T8,P8,fluido,eff_ltr,Nr,m_CO2,m_CO2*splitA);

        % Área LTR
            for ii=1:26
                Hr=HX(ii,1); Wr=HX(ii,2); Lpr=HX(ii,3);
                [varout]=LTR(Lpr,Wr,Hr,T9,T8,T5,T6,Phigh,Phigh,Q_ltr,fluido,m_CO2,m_CO2*splitA);
                A_LTR=varout(2);
                if ii==1
                    Amin_LTR=A_LTR;
                else
                    if (A_LTR<Amin_LTR)&&(varout(3)<HX(ii,4))&&(A_LTR<HX(ii,5))
                        Amin_LTR=A_LTR;
                    end
                end
            end
        
        % Exergía destruida en LTR
            Xd_ltr = m_CO2*(x5-x6)+m_CO2*splitA*(x8-x9);
        
%---------------------------------------------------------------------
    % ANÁLISIS DEL COOLER (Puntos: 6AO3, 6AO, 13, 14) ACOPLE 
        [T6AO, A_cooler, Q_cooler, Xd_cooler, x_loss_cooler, TPhsmx_cooler] = ...
        Cooler_it(T6AO3, Plow, m_CO2*splitA, pmaire);

%--------------------------------------------------------------------------
    % ANÁLISIS DE LA PEM (Proton Exchange Membrane) (Puntos: 19, 20, 6AO, 7) ACOPLE
        % Potencias -------------------------
            W_Brayton = W_t1+W_t2+W_c1+W_c2;    % W
            W_Rankine = W_t3+W_b2;              % W
            WNETO = W_Brayton+W_Rankine;        % W
            
        [TPhsmx_PEM, Xd_PEM_TOT, A_HE, Q_HE, ~, ~, m_H2O_in, m_H2] = PEM_OFICIAL(T6AO, Plow, m_CO2*splitA, T7, W_Rankine);
        T19 = TPhsmx_PEM(3,1);
        T20 = TPhsmx_PEM(4,1);
        LHV_H2 = 120100*1000;                   % J/kg
        Cp_H2O = 4.186;                         % kJ/kgK
        Q_H2O = m_H2O_in*Cp_H2O*(T20-T19);
%--------------------------------------------------------------------------

    % Eficiencias ciclo --------------------------------
        % Eficiencia térmica Rankine
            n_th_ORC = W_Rankine/Q_evap_vec(4)*100;
            
        % Eficiencia térmica Brayton
            n_th_Brayton = W_Brayton/Q_tt*100;
            
        % Eficiencia energética PEM
            n_en_PEM = (m_H2*LHV_H2/((W_Rankine+Q_H2O*1000)))*100;
            
        % Eficiencia térmica
            n_th = ((W_Brayton+m_H2*LHV_H2)/Q_tt)*100;
            
        % Eficiencia energética SIN PEM
            n_ciclo = (WNETO/Q_tt)*100;
            
        % Eficiencia exergetica Rankine
            n_ex_ORC = (W_Rankine/(m_CO2*splitA*(x6AI-x6AO3)))*100;
            
        % Eficiencia exergetica Brayton
            n_ex_Brayton = W_Brayton/X_solar*100;
            
        % Eficiencia exergetica PEM --------------------
            n_ex_PEM = (Xd_PEM_TOT(2)/(W_Rankine/1000+Xd_PEM_TOT(3)))*100;
            
        % Eficiencia exergética ------------------------
            n_ex = ((W_Brayton/1000+Xd_PEM_TOT(2))/(X_solar/1000))*100;
         
% --------------------------------------------------------------------------------------------------------------------------------------------------
% RESULTADOS
    % Exergia destruida
        XD(1,1) = Xd_c1/1000;           % Comprensor 1
        XD(1,2) = Xd_c2/1000;           % Comprensor 2
        XD(1,3) = Xd_t1/1000;           % Turbina 1
        XD(1,4) = Xd_t2/1000;           % Turbina 2
        XD(1,5) = Xd_orc(2)/1000;       % Turbina 3
        XD(1,6) = Xd_htr/1000;          % HTR
        XD(1,7) = Xd_ltr/1000;          % LTR
        XD(1,8) = Xloss_sr/1000;        % Loss receiver
        XD(1,9) = Xloss_field/1000;     % Loss field
        XD(1,10) = Xd_sr/1000;          % Destruida Re_Heater 
        XD(1,11) = Xd_cooler/1000;      % Destruida cooler
        XD(1,12) = x_loss_cooler/1000;  % Perdida en el aire
        XD(1,13) = Xd_orc(4)/1000;      % Destruida evaporador
        XD(1,14) = Xd_orc(3)/1000;      % Destruida condensador
        XD(1,15) = Xd_orc(1)/1000;      % Destruida pump 2
        XD(1,16) = Xd_orc(5)/1000;      % Loss en el agua
        XD(1,17) = Xd_PEM_TOT(1)/1000;  % Destruida HE
        XD(1,18) = X_solar/1000;        % Exergia aportada por el sistema solar
        XD(1,19) = Xd_PEM_TOT(4);       % Exergía destruida por la PEM
        Xd_total = sum(XD(1,1:19));  
        XD(1,20) = Xd_total-X_solar/1000;% Exergia destruida total
        
        % Eficiencia exergética ------------------------
            n_ex_2 = (1-(XD(1,20)/(X_solar/1000)))*100;   
            
    % Áreas de los intercambiadores
        Areas(1,1) = A_HTR;          % Área HTR     
        Areas(1,2) = A_LTR;          % Área LTR
        Areas(1,3) = A_ORC(1);       % Área Evaporador 
        Areas(1,4) = A_ORC(2);       % Área Condensador
        Areas(1,5) = A_cooler;       % Área Cooler
        Areas(1,6) = A_HE;           % Área Heat Exchanger
        
    % Calor
        Calor(1,1) = Q_heater/1000000;      % +Calor Heater      
        Calor(1,2) = Q_reheater/1000000;    % +Calor reheater
        Calor(1,3) = Q_tt/1000000;          % +Calor total entrada Brayton
        Calor(1,4) = Q_htr/1000000;         % -Calor linea caliente a fria en htr
        Calor(1,5) = Q_ltr/1000000;         % -Calor linea caliente a fria en ltr
        Calor(1,6) = Q_evap_vec(1)/1000000; % +Calor recibido etapa 1 evaporador
        Calor(1,7) = Q_evap_vec(2)/1000000; % +Calor recibido etapa 2 evaporador
        Calor(1,8) = Q_evap_vec(3)/1000000; % +Calor recibido etapa 3 evaporador
        Calor(1,9) = Q_evap_vec(4)/1000000; % +Calor recibido total evaporador
        Calor(1,10) = Q_cond_vec(1)/1000000;% +Calor recibido etapa 1 condensador
        Calor(1,11) = Q_cond_vec(2)/1000000;% +Calor recibido etapa 2 condensador
        Calor(1,12) = Q_cond_vec(3)/1000000;% +Calor recibido total condensador
        Calor(1,13) = Q_cooler/1000000;     % +Calor recibido cooler
        Calor(1,14) = Q_HE/1000000;         % -Calor cedido HE 
    
    % Propiedades BRAYTON
        T = [T1,T2,T3,T4,T5,T6,T6AI,T6AO1,T6AO2,T6AO3,T6AO,T6BI,T7,T8,T9,T10,T11,T12];      % Temperaturas
        P = [P1,P2,P3,P4,P5,P6,P6AI,P6AO1,P6AO2,P6AO3,TPhsmx_cooler(4,2)*1000,P6BI,P7,P8,P9,P10,P11,P12];      % Presiones
        H = [h1,h2,h3,h4,h5,h6,h6AI,h6AO1,h6AO2,h6AO3,TPhsmx_cooler(4,3)*1000,h6BI,h7,h8,h9,h10,h11,h12];      % Entalpias
        S = [s1,s2,s3,s4,s5,s6,s6AI,s6AO1,s6AO2,s6AO3,TPhsmx_cooler(4,4)*1000,s6BI,s7,s8,s9,s10,s11,s12];      % Entropias
        m = [m_CO2,m_CO2,m_CO2,m_CO2,m_CO2,m_CO2,m_CO2*splitA,m_CO2*splitA,m_CO2*splitA,m_CO2*splitA,m_CO2*splitA,m_CO2*splitB,m_CO2*splitA,m_CO2*splitA,m_CO2*splitA,m_CO2*splitB,m_CO2,m_CO2];       % J/kg
        x = [x1,x2,x3,x4,x5,x6,x6AI,x6AO1,x6AO2,x6AO3,TPhsmx_cooler(4,6)*1000,x6BI,x7,x8,x9,x10,x11,x12];      % kg/s
          
        TPhsmx_brayton = [T',(P./1000)',(H./1000)',(S./1000)',m',(x./1000)'];
        Propierties=[TPhsmx_brayton; TPhsmx_orc; TPhsmx_cooler; TPhsmx_PEM];
        
    % Potencias
        Potencias(1,1) = W_t1/1000;
        Potencias(1,2) = W_t2/1000;
        Potencias(1,3) = W_t3/1000;
        Potencias(1,4) = W_c1/1000;
        Potencias(1,5) = W_c2/1000;
        Potencias(1,6) = W_b2/1000;
        Potencias(1,7) = W_Brayton/1000;
        Potencias(1,8) = W_Rankine/1000;
        Potencias(1,9) = WNETO/1000;    
        
    % Eficiencias
        Eficiencias(1,1) = n_th_ORC;
        Eficiencias(1,2) = n_th_Brayton;
        Eficiencias(1,3) = n_en_PEM;
        Eficiencias(1,4) = n_th;
        Eficiencias(1,5) = n_ciclo;
        Eficiencias(1,6) = n_ex_ORC;
        Eficiencias(1,7) = n_ex_Brayton;
        Eficiencias(1,8) = n_ex_PEM;
        Eficiencias(1,9) = n_ex;
        Eficiencias(1,10)= n_ex_2;
        
% -------------------------------------------------------------------------
% CALCULO AMBIENTAL -------------------------------------------------------
    EWR = XD(1,20)*1000/(X_solar);   % Exergy waste ratio EWR
    EEF = EWR/(n_ex/100);       % Environmental effect factor (EEF)
    ESI = 1/EEF;                % Exergetic sustainability index (ESI)

    In(1,1)=EWR;
    In(1,2)=EEF;
    In(1,3)=ESI;
    
    W1    = abs(W_t1)/1000;
    W2    = abs(W_t2)/1000;
    Wc1   = abs(W_c1)/1000;
    Wc2   = abs(W_c2)/1000;
    Wt3   = abs(W_t3)/1000;
    Wb2   = abs(W_b2)/1000;
    W_PEM = W_Rankine/1000; 
    
    [YCOPROD, Ytotal, Yco_total, Yde_total, Yom_total, Y_componentes, Ytotal_PEM,Yco_total_PEM, Yde_total_PEM, Yom_total_PEM] = EnvironmentalV22(W1,W2,Wt3,Wc1,Wc2,Wb2,Areas(1,1), Areas(1,2), Areas(1,5), Areas(1,3) ,Areas(1,4), Areas(1,6));
    
    

% -------------------------------------------------------------------------
% CALCULO EXERGOECONOMICO
         % TODO el trabajo neto disponible se envia a la PEM

    % Z--------------------------------------------------------------------
    % Actualización de costos por año
    CEPCI_2022 = 816;
    F_TC   = CEPCI_2022/381.7; 
    F_HXC  = CEPCI_2022/318.4;  
    F_HXS  = CEPCI_2022/575.5; 
    F_T3   = CEPCI_2022/584.6;  
    F_PUMP = CEPCI_2022/468.2;  
    F_PEM  = CEPCI_2022/368.1;  

    Z_t1       = F_TC  * 479.34*m_CO2*(1/(0.93-nt))*log(r1)*(1+exp(0.036*T1-54.4));
    Z_t2       = F_TC  * 479.34*m_CO2*(1/(0.93-nt))*log(r2)*(1+exp(0.036*T3-54.4));
    Z_c1       = F_TC  * 71.1*m_CO2*splitA*(1/(0.92-nc))*r1*log(r1);
    Z_c2       = F_TC  * 71.1*m_CO2*splitB*(1/(0.92-nc))*r1*log(r1);
    Z_HTR      = F_HXC * 2681*(A_HTR^0.59);
    Z_LTR      = F_HXC * 2681*(A_LTR^0.59);
    Z_ORC_EVAP = F_HXC * 2681*(A_ORC(1)^0.59);
    Z_t3       = F_T3  * 4405*(Wt3^0.7);
    Z_b2       = F_PUMP* 1120*(Wb2^0.8);
    Z_ORC_COND = F_HXS * 2143*(A_ORC(2)^0.514);
    Z_COOLER   = F_HXS * 2143*(A_cooler^0.514);
    Z_HE       = F_HXS * 2143*(A_HE^0.514);
    Z_PEM      = F_PEM * 1000*W_PEM;

    % Tasa de costo de inversión [USD/h] ----------------------------------
    CRF = (interes*(1+interes)^years)/((1+interes)^years - 1);
    f_Z = (CRF + gamma)/tao;

    Zdot_t1       = Z_t1*f_Z;
    Zdot_t2       = Z_t2*f_Z;
    Zdot_c1       = Z_c1*f_Z;
    Zdot_c2       = Z_c2*f_Z;
    Zdot_HTR      = Z_HTR*f_Z;
    Zdot_LTR      = Z_LTR*f_Z;
    Zdot_ORC_EVAP = Z_ORC_EVAP*f_Z;
    Zdot_t3       = Z_t3*f_Z;
    Zdot_b2       = Z_b2*f_Z;
    Zdot_ORC_COND = Z_ORC_COND*f_Z;
    Zdot_COOLER   = Z_COOLER*f_Z;
    Zdot_HE       = Z_HE*f_Z;
    Zdot_PEM      = Z_PEM*f_Z;

    % Costos unitarios ----------------------------------------------------
    c_solar = 0.02;   % USD/kWh_ex 
    c_water = 0.002;  % USD/kWh_ex

    % Exergias de corrientes [kW] -----------------------------------------
    E1    = (m_CO2*x1)/1000;
    E2    = (m_CO2*x2)/1000;
    E3    = (m_CO2*x3)/1000;
    E4    = (m_CO2*x4)/1000;
    E5    = (m_CO2*x5)/1000;
    E6    = (m_CO2*x6)/1000;
    E7    = (m_CO2*splitA*x7)/1000;
    E8    = (m_CO2*splitA*x8)/1000;
    E9    = (m_CO2*splitA*x9)/1000;
    E10   = (m_CO2*splitB*x10)/1000;
    E11   = (m_CO2*x11)/1000;
    E12   = (m_CO2*x12)/1000;
    E6BI  = (m_CO2*splitB*x6BI)/1000;
    E6AI  = (m_CO2*splitA*x6AI)/1000;
    E6AO3 = (m_CO2*splitA*x6AO3)/1000;
    E6AO  =  m_CO2*splitA*TPhsmx_cooler(4,6);   % kW
    E19   =  m_H2O_in*TPhsmx_PEM(3,6);          % kW
    E20   =  m_H2O_in*TPhsmx_PEM(4,6);          % kW
    E15   = TPhsmx_orc(5,5)*TPhsmx_orc(5,6);
    E16   = TPhsmx_orc(6,5)*TPhsmx_orc(6,6);
    E17   = TPhsmx_orc(1,5)*TPhsmx_orc(1,6);
    E18   = TPhsmx_orc(2,5)*TPhsmx_orc(2,6);

    EH2   = Xd_PEM_TOT(2);   % kW
    EO2   = Xd_PEM_TOT(5);   % kW

    names = {'c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12', ...
             'c6BI','c6AI','c6AO3','c6AO','c19','c20','c15','c16','c17','c18', ...
             'cH2','cO2','cel'};
    n = numel(names);
    idx = containers.Map(names,1:n);
    A = zeros(n,n);
    b = zeros(n,1);

    % (1) Heater solar
    A(1,idx('c1'))  =  E1;
    A(1,idx('c12')) = -E12;
    b(1) = c_solar*(E1 - E12);

    % (2) Reheater solar
    A(2,idx('c3')) =  E3;
    A(2,idx('c2')) = -E2;
    b(2) = c_solar*(E3 - E2);

    % (3) Electricidad común: T1 + T2 + T3
    A(3,idx('c2'))  =  E2;
    A(3,idx('c1'))  = -E1;
    A(3,idx('c4'))  =  E4;
    A(3,idx('c3'))  = -E3;
    A(3,idx('c16')) =  E16;
    A(3,idx('c15')) = -E15;
    A(3,idx('cel')) =  (W1 + W2 + Wt3);
    b(3) = Zdot_t1 + Zdot_t2 + Zdot_t3;

    % Reglas turbinas
    A(4,idx('c1'))  = 1;  A(4,idx('c2'))  = -1;  b(4) = 0;
    A(5,idx('c3'))  = 1;  A(5,idx('c4'))  = -1;  b(5) = 0;
    A(6,idx('c15')) = 1;  A(6,idx('c16')) = -1;  b(6) = 0;

    % (7) Compresor 1
    A(7,idx('c8'))  =  E8;
    A(7,idx('c7'))  = -E7;
    A(7,idx('cel')) = -Wc1;
    b(7) = Zdot_c1;

    % (8) Compresor 2
    A(8,idx('c10')) =  E10;
    A(8,idx('c6BI'))= -E6BI;
    A(8,idx('cel')) = -Wc2;
    b(8) = Zdot_c2;

    % (9) Bomba ORC
    A(9,idx('c18')) =  E18;
    A(9,idx('c17')) = -E17;
    A(9,idx('cel')) = -Wb2;
    b(9) = Zdot_b2;

    % (10) HTR
    A(10,idx('c12')) =  E12;
    A(10,idx('c5'))  =  E5;
    A(10,idx('c11')) = -E11;
    A(10,idx('c4'))  = -E4;
    b(10) = Zdot_HTR;

    % (11) LTR
    A(11,idx('c6')) =  E6;
    A(11,idx('c9')) =  E9;
    A(11,idx('c5')) = -E5;
    A(11,idx('c8')) = -E8;
    b(11) = Zdot_LTR;

    % (12) Evaporador ORC (acople CO2-ORC)
    A(12,idx('c6AO3')) =  E6AO3;
    A(12,idx('c15'))   =  E15;
    A(12,idx('c6AI'))  = -E6AI;
    A(12,idx('c18'))   = -E18;
    b(12) = Zdot_ORC_EVAP;

    % (13) PEM
    % C_H2O,out = C_20, el costo del agua se cancela en la PEM.
    % La PEM carga costo al H2 a partir del trabajo electrico consumido.
    A(13,idx('cH2')) =  EH2;
    A(13,idx('cO2')) =  EO2;
    A(13,idx('cel')) = -W_PEM;
    b(13) = Zdot_PEM;

    % (14) Cooler CO2-ambiente
    A(14,idx('c6AO'))  =  E6AO;
    A(14,idx('c6AO3')) = -E6AO3;
    b(14) = Zdot_COOLER;

    % (15) HE de calentamiento de agua para PEM
    A(15,idx('c7'))   =  E7;
    A(15,idx('c20'))  =  E20;
    A(15,idx('c6AO')) = -E6AO;
    A(15,idx('c19'))  = -E19;
    b(15) = Zdot_HE;

    % (16) Condensador ORC 
    A(16,idx('c17')) =  E17;
    A(16,idx('c16')) = -E16;
    b(16) = Zdot_ORC_COND;

    % (17) Mezclador
    A(17,idx('c11')) =  E11;
    A(17,idx('c9'))  = -E9;
    A(17,idx('c10')) = -E10;
    b(17) = 0;

    % (18)-(23) Ecuaciones auxiliares 
    A(18,idx('c4'))   = 1; A(18,idx('c5'))   = -1; b(18) = 0;   % HTR, lado combustible
    A(19,idx('c5'))   = 1; A(19,idx('c6'))   = -1; b(19) = 0;   % LTR, lado combustible
    A(20,idx('c6'))   = 1; A(20,idx('c6AI')) = -1; b(20) = 0;   % split hacia ORC
    A(21,idx('c6'))   = 1; A(21,idx('c6BI')) = -1; b(21) = 0;   % split hacia compresor 2
    A(22,idx('c6AI')) = 1; A(22,idx('c6AO3'))= -1; b(22) = 0;   % evaporador ORC, lado combustible
    A(23,idx('c6AO')) = 1; A(23,idx('c7'))   = -1; b(23) = 0;   % HE agua-PEM, lado combustible

    % (24) Agua de reposicion
    A(24,idx('c19')) = 1;
    b(24) = c_water;

    % (25) Oxigeno sin valorizacion economica
    A(25,idx('cO2')) = 1;
    b(25) = 0;

    x_costos = A\b;
    
    sol = struct();
    for k = 1:n
        sol.(names{k}) = x_costos(k);
    end
    % Flujos de costo [USD/h] ---------------------------------------------
    Vector_E = [E1; E2; E3; E4; E5; E6; E7; E8; E9; E10; E11; E12; ...
                E6BI; E6AI; E6AO3; E6AO; E19; E20; E15; E16; E17; E18; ...
                EH2; EO2; W_PEM];
    Cdot = x_costos .* Vector_E;

    % Parametros exergoeconomicos por componente --------------------------
    Comp_Names = {'Turbina 1'; 'Turbina 2'; 'Compresor 1'; 'Compresor 2'; ...
                  'HTR'; 'LTR'; 'Evaporador ORC'; 'Turbina ORC'; 'Bomba ORC'; ...
                  'Condensador ORC'; 'Cooler'; 'HE-PEM'; 'PEM'};

    Zdot = [Zdot_t1; Zdot_t2; Zdot_c1; Zdot_c2; Zdot_HTR; Zdot_LTR; ...
               Zdot_ORC_EVAP; Zdot_t3; Zdot_b2; Zdot_ORC_COND; Zdot_COOLER; ...
               Zdot_HE; Zdot_PEM];

    ED_k = [XD(1,3); XD(1,4); XD(1,1); XD(1,2); XD(1,6); XD(1,7); ...
            XD(1,13); XD(1,5); XD(1,15); XD(1,14); XD(1,11); XD(1,17); XD(1,19)];

    cF_k = [sol.c1; sol.c3; sol.cel; sol.cel; sol.c4; sol.c5; ...
            sol.c6AI; sol.c15; sol.cel; sol.c16; sol.c6AO3; sol.c6AO; sol.cel];

    CD_k = cF_k .* ED_k;
    f_k  = 100*(Zdot ./(Zdot + CD_k));

    Cdot_fuel = c_solar * (X_solar / 1000);   % [USD/h]
    E_H2_LHV   = m_H2*LHV_H2/1000;                              % kW
    W_grid = max(W_Brayton,0) / 1000;                           % [kW]

    Zdot_total = Z_t1*f_Z + Z_t2*f_Z + Z_c1*f_Z + Z_c2*f_Z + ...
                 Z_HTR*f_Z + Z_LTR*f_Z + Z_ORC_EVAP*f_Z + ...
                 Z_t3*f_Z + Z_b2*f_Z + Z_ORC_COND*f_Z + ...
                 Z_COOLER*f_Z + Z_HE*f_Z + Z_PEM*f_Z;      % [USD/h]
             

    LCOEn = (Cdot_fuel + Zdot_total) / max(W_grid + E_H2_LHV, eps);   % [USD/kWh]
    LCOH = sol.cH2;           % [USD/kWh_ex]
    c_el_bus = sol.cel;       % [USD/kWh_ex]
    
% Salida resultados ---------------------------------------------------

    % Escribir resultados previos
    xlswrite('Salidas.xlsx',Propierties,'Energy_Exergy','C3')
    xlswrite('Salidas.xlsx',Calor','Energy_Exergy','L4')
    xlswrite('Salidas.xlsx',Areas','Energy_Exergy','L20')
    xlswrite('Salidas.xlsx',vector_field','Energy_Exergy','O27')        
    xlswrite('Salidas.xlsx',Potencias','Energy_Exergy','O4')   
    xlswrite('Salidas.xlsx',Eficiencias','Energy_Exergy','O15')  
    xlswrite('Salidas.xlsx',Yco_total,'Energy_Exergy','R23')  
    xlswrite('Salidas.xlsx',Yde_total,'Energy_Exergy','R24')  
    xlswrite('Salidas.xlsx',Yom_total,'Energy_Exergy','R22')  
    xlswrite('Salidas.xlsx',Y_componentes','Energy_Exergy','R4')
    xlswrite('Salidas.xlsx',Yco_total_PEM,'Energy_Exergy','T22') 
    xlswrite('Salidas.xlsx',Yde_total_PEM,'Energy_Exergy','T23')         
    xlswrite('Salidas.xlsx',Yom_total_PEM,'Energy_Exergy','T24')         
    xlswrite('Salidas.xlsx',YCOPROD,'Energy_Exergy','T25')         

    Titulos_Flujos = {'Fluido', 'Costo Unit. (USD/kWh_ex)', 'Cdot (USD/h)'};
    Datos_Flujos = [names', num2cell(x_costos), num2cell(Cdot)];
    xlswrite('Salidas.xlsx', Titulos_Flujos, 'Exergoeconomico', 'A1');
    xlswrite('Salidas.xlsx', Datos_Flujos, 'Exergoeconomico', 'A2');

 
    Datos_Comp = [Comp_Names, num2cell(Zdot), num2cell(ED_k), num2cell(CD_k), num2cell(f_k)];
    xlswrite('Salidas.xlsx', Datos_Comp, 'Exergoeconomico', 'E2');

    Titulos_Globales = {'LCOEn sistema [USD/kWh]', ...
                        'c_H2 [USD/kWh_ex]', ...
                        'c_el_bus interno [USD/kWh_ex]'};

    Datos_Globales = [LCOEn, LCOH, c_el_bus];
    xlswrite('Salidas.xlsx', Titulos_Globales, 'Exergoeconomico', 'K1');
    xlswrite('Salidas.xlsx', Datos_Globales,   'Exergoeconomico', 'K2');
end