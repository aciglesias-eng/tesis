function [T6AO3, A_ORC, W_t3, W_b2, Xd_orc, TPhsmx_orc, Q_cond_vec, Q_evap_vec] = ORCV2(Ts, P_CO2, m_CO2)
%function[T6,A,W_t3,W_b2,Xd,TPhsmx_orc]=ORC_fastV2(Ts,P_CO2,m_CO2)

%ENTRADAS A LA FUNCIÓN-------------------------
% Ts = Temperatura de la fuente caliente (6AI)
% P_CO2 = Presión de la línea caliente (Brayton)
% m_CO2 = Flujo másico de la línea caliente

% SALIDAS DE LA FUNCIÓN-------------------------
% T6AO3 = Temperatura de salida de la línea caliente
% A = Vector de áreas requeridas de intercambiadores
% W_t3 = Potencia generada por la turbina del ORC
% W_b2 = Potencia consumida por la bomba
% XD_orc = Vector con la exergía destruida por componente
% TPhsmx_orc = Matriz resumen de estados del ciclo: T, P, h, s, Ex y flujo másico en cada punto
    
%Entradas SISTEMA ORC----------------------------------------------------------------------------------------------------------------------
    fluid = 'toluene';      % Fluido de trabajo ORC
    np = 0.85;              % Eficiencia isentrópica de la bomba
    nt = 0.85;              % Eficiencia isentrópica de la turbina
    Ap = 39;                % Temperatura minima de aproximacion en el evaporador
    pinch = 10.0;           % Diferencia mínima de temperatura (pinch) en intercambiadores (K)
    Tc = 40+273.15;         % Temperatura mínima del ciclo
    T0 = 25+273.15;         % Temperatura del refrigerante/agua que entra al condensador (°C)
    r_ORC = 4;              % Relación de presión en ORC
    Plow = refpropm('P','T',Tc,'Q',0,fluid); % Presion minima del ciclo ORC
    Phigh = Plow*r_ORC;                      % Presion máxima del ciclo ORC
    Pcrit = refpropm('P','C',0,' ',0,fluid); % Presion crítica del ciclo ORC
    P0 = 101.325;           % Presión atmosférica (kPa) 
    
    % Estado de referencia para exergía (estado muerto) Tolueno
        [h0t,s0t] = refpropm('HS','T',25+273.15,'P',101.3,fluid);
    % Estado de referencia para exergía (estado muerto) CO2
        [h0c,s0c] = refpropm('HS','T',25+273.15,'P',101.3,'co2');
    % Estado de referencia para exergía (estado muerto) agua
        [h0w,s0w] = refpropm('HS','T',25+273.15,'P',101.3,'water');
    display(Pcrit);
    if Phigh >= Pcrit
        T6AO3=0; A_ORC=[0;0]; W_t3=0; W_b2=0; Xd_orc=zeros(1,5);
        Q_cond_vec=zeros(3,1); Q_evap_vec=zeros(4,1);
        TPhsmx_orc = zeros(14,6);
        return
    end

   
% Puntos --------------------------------------------------------------------------------------------------------------------------------------------------
    % ORC
        % Punto 17
            T17 = Tc; P17 = Plow;
            [h17, s17] = refpropm('HS','T',T17,'P',P17,fluid);
            x17 = (h17-h0t)-T0*(s17-s0t);

        % Punto 18
            P18 = Phigh;
            s18 = s17;
            h18s = refpropm('H','P',P18,'S',s18,fluid);
            h18 = h17+(h18s-h17)/np;
            T18 = refpropm('T','P',P18,'H',h18,fluid);
            x18 = (h18-h0t)-T0*(s18-s0t);

        % Punto 15ls
            P15ls = P18;
            [T15ls,h15ls] = refpropm('TH','P',P15ls,'Q',0,fluid);
            s15ls = refpropm('S','P',Phigh,'Q',0,fluid);
            x15ls = (h15ls-h0t)-T0*(s15ls-s0t);

        % Punto 15vs
            T15vs = T15ls; P15vs = P15ls;
            [h15vs,s15vs] = refpropm('HS','P',P15vs,'Q',1,fluid);
            x15vs = (h15vs-h0t)-T0*(s15vs-s0t);

        % Punto 15
            P15 = P15vs;
            T15 = Ts-Ap;
            [h15,s15] = refpropm('HS','T',T15,'P',P15,fluid);
            x15 = (h15-h0t)-T0*(s15-s0t);
            if h15 < h15vs
                T6AO3=0; A_ORC=[0;0;0]; W_t3=0; W_b2=0; Xd_orc=zeros(1,5); TPhsmx_orc=zeros(14,6);
                Q_cond_vec=zeros(3,1); Q_evap_vec=zeros(4,1);
                return
            end
        % Punto 16
            P16 = Plow;
            h16s = refpropm('H','P',P16,'S',s15,fluid);
            h16 = h15-(h15-h16s)*nt;
            [T16,s16] = refpropm('TS','P',P16,'H',h16,fluid);
            x16 = (h16-h0t)-T0*(s16-s0t);
            
        % Punto 16vs
           T16vs = T17; P16vs = P17;
           [h16vs,s16vs] = refpropm('HS','T',T16vs,'Q',1,'toluene');
           x16vs = (h16vs-h0t)-T0*(s16vs-s0t);
           
% --------------------------------------------------------------------------------------------------------------------------------------------------        
    % Brayton (CO2) y Condensador
        % Punto 6AI B
            P6AI = P_CO2;
            T6AI = Ts;
            [h6AI,s6AI] = refpropm('HS','T',Ts,'P',P6AI,'co2');
            x6AI = (h6AI-h0c) - T0*(s6AI-s0c);

        % Punto 6AO2  B
            T6AO2 = T15ls+pinch;                     
            P6AO2 = P_CO2;
            [h6AO2,s6AO2] = refpropm('HS','T',T6AO2,'P',P6AO2,'co2');  
            x6AO2 = (h6AO2-h0c)-T0*(s6AO2-s0c);
        
        % Punto C1 C
            TC1 = T0; PC1 = P0;
            [hC1,sC1] = refpropm('HS','T',TC1,'P',PC1,'water');
            xC1 = (hC1-h0w)-T0*(sC1-s0w);
            
        % Punto C2 C
            PC2 = P0;
            TC2 = T17-pinch;
            [hC2,sC2] = refpropm('HS','T',TC2,'P',PC2,'water');
            xC2 = (hC2-h0w)-T0*(sC2-s0w);
            

    %---------------------------------------------------------------------
        % Flujo másico (línea ORC)
            T6AOCp = 0.5*(T6AO2+Ts); % Cálculo de una temperatura promedio para evaluar Cp
            CpHS = refpropm('C','T',T6AOCp,'P',P_CO2,'co2');
            m_ORC = m_CO2*CpHS*(Ts-T6AO2)/(h15-h15ls);
            
        % Flujo másico (línea Condensador)
            h17p = refpropm('H','T',T17,'Q',1,fluid);     % Entalpía de vapor saturado a Plow 
            if h16<h17p                                   % Entra al condensador como mezcla húmeda
                m_water = m_ORC*(h16-h17)/(hC2-hC1);      % Agua calcula calor total hasta líquido saturado
            else                                          % Entra como vapor o sobrecalentado
                m_water = m_ORC*(h17p-h17)/(hC2-hC1);     % Agua dimensionada con calor latente
            end

    %---------------------------------------------------------------------   
        % Punto C3 C
            PC3 = P0;
            hC3 = hC1+m_ORC*(h16-h17)/m_water;
            [TC3,sC3] = refpropm('TS','P',PC3,'H',hC3,'water');  
            xC3 = (hC3-h0w)-T0*(sC3-s0w);
            xloss_cond = (xC1-xC3)*m_water;
    
        % Punto 6AO1 B
            T6AO1 = Ts-m_ORC*(h15-h15vs)/(m_CO2*CpHS);
            P6AO1 = P_CO2;
            [h6AO1,s6AO1] = refpropm('HS','T',T6AO1,'P',P6AO1,'co2');  
            x6AO1 = (h6AO1-h0c)-T0*(s6AO1-s0c);
            
        % Punto 6AO3 B
            T6AO3 = Ts-m_ORC*(h15-h18)/(m_CO2*CpHS);  
            P6AO3 = P_CO2;
            [h6AO3,s6AO3] = refpropm('HS','T',T6AO3,'P',P6AO3,'co2');        
            x6AO3 = (h6AO3-h0c)-T0*(s6AO3-s0c);

% --------------------------------------------------------------------------------------------------------------------------------------------------        
%---------------------------------------------------------------------
    % ANÁLISIS DE LA BOMBA 2 (Puntos: 17, 18)
        % Potencia de la bomba 2
            W_b2 = m_ORC*(h17-h18);

        % Exergía destruida en la bomba 2
            Xd_b2 = -W_b2+m_ORC*(x17-x18);
            
%---------------------------------------------------------------------
    % ANÁLISIS DE LA TURBINA 3 (Puntos: 15, 16)
        % Potencia de la turbina 3
            W_t3 = m_ORC*(h15-h16);

        % Exergía destruida en la turbina 3
            Xd_t3 = -W_t3+m_ORC*(x15-x16);

%---------------------------------------------------------------------
    % ANÁLISIS DEL CONDENSADOR (Puntos: C1, C2, C3, 16, 16vs, 17)
        % Calor del condensador
            Q_condensador_e1 = m_water*(hC2-hC1);
            if h16<h17p
                Q_condensador_e2 = 0;
            else
                Q_condensador_e2 = m_water*(hC3-hC2);
            end
        % Exergía destruida en el condensador
            Xd_condensador = m_ORC*(x16-x17)+m_water*(xC1-xC3);
            
%---------------------------------------------------------------------
    % ANÁLISIS DEL EVAPORADOR (Puntos: 6AI, 6AO1, 6AO2, 6AO3, 18, 15ls, 15vs, 15)
        % Calor del evaporador por etapas
            Q_evaporador_e1 = m_ORC*(h15ls-h18);
            Q_evaporador_e2 = m_ORC*(h15vs-h15ls);
            Q_evaporador_e3 = m_ORC*(h15-h15vs);
            Q_evaporador    = m_ORC*(h15-h18);
            
        % Exergía destruida en el evaporador
            Xd_evaporador = m_CO2*(x6AI-x6AO3)+m_ORC*(x18-x15);
%---------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------------------------
    Heat(1)=Q_evaporador_e1;
    Heat(2)=Q_evaporador_e2;
    Heat(3)=Q_evaporador_e3;
    HX=[460,160,336,125,5;
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
        1470,610,941.400000000000,700,280;
        1835,610,1306.20000000000,700,420;
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
    
          

    for ii=1:26
         H=HX(ii,1); W=HX(ii,2); Lp=HX(ii,3);
         [varout] = Evaporador(Lp,W,H,T15,T15ls,T18,Ts,T6AO1,T6AO2,T6AO3, Phigh, m_ORC, m_CO2, Heat, fluid, P_CO2);
         Aevap=varout(6);
         if ii==1
            Amin=varout(6);
         else                    
             if (varout(6)<Amin)&&(varout(7)<HX(ii,4))&&(varout(6)<HX(ii,5))
                 Amin=varout(6);
             end
        end
    end
    

    Qcond = Q_condensador_e1;
    Qenf  = Q_condensador_e2;
    x16q  = refpropm('Q','P',P16,'H',h16,fluid);
    

    for ii=1:26
         Hc=HX(ii,1); Wc=HX(ii,2); Lpc=HX(ii,3);
         [varout] = Condensador(Lpc,Wc,Hc,Qcond,Qenf,T16,T17,T17,TC1,TC2,TC3,m_ORC,m_water,Plow, x16q, fluid);
         Acond=varout(1)+varout(2);
         if ii==1
             Amin=Acond;
         else
             if (Acond<Amin)&&(varout(3)<HX(ii,4))&&(Acond<HX(ii,5))
                 Amin=Acond;
             end
         end
    end
    
    A_ORC = [Aevap;Acond];
%-------------------------------------------------------------------------------------------------------------------------------------------------
% RESULTADOS
    Xd_orc = [Xd_b2, Xd_t3,Xd_condensador, Xd_evaporador, -xloss_cond];

    %A = [Aevap;Aregen;Acond];
    T  = [T17, T18, T15ls, T15vs, T15, T16, T16vs, TC1, TC2, TC3, T6AI, T6AO1, T6AO2, T6AO3];  % K
    P  = [P17, P18, P15ls, P15vs, P15, P16, P16vs,  PC1, PC2, PC3, P6AI, P6AO1, P6AO2, P6AO3];  % kPa
    H  = [h17, h18, h15ls, h15vs, h15, h16, h16vs, hC1, hC2, hC3, h6AI, h6AO1, h6AO2, h6AO3];  % J/kg
    S  = [s17, s18, s15ls, s15vs, s15, s16, s16vs, sC1, sC2, sC3, s6AI, s6AO1, s6AO2, s6AO3];  % J/kg-k
    m  = [m_ORC, m_ORC, m_ORC, m_ORC, m_ORC, m_ORC, m_ORC, m_water,m_water, m_water, m_CO2, m_CO2, m_CO2, m_CO2]; % kg/s
    x  = [x17, x18, x15ls, x15vs, x15, x16, x16vs, xC1, xC2, xC3, x6AI, x6AO1, x6AO2, x6AO3];  % J/kg

    TPhsmx_orc = [T', P', (H./1000)', (S./1000)', m',(x./1000)'];
    Q_evap_vec = [Q_evaporador_e1; Q_evaporador_e2; Q_evaporador_e3; Q_evaporador];
    Q_cond_vec = [Q_condensador_e1; Q_condensador_e2; (Q_condensador_e1 + Q_condensador_e2)];
end





            
            
