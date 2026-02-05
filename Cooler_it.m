function [T6AO, A_Cooler, Q_cooler, Xd_cooler,x_loss_cooler, TPhsmx_cooler] = Cooler_it( ...
    T6AO3, P_CO2, m_CO2, apertura_pct)

% ENTRADAS
%   T6AO3        	Temp CO2 entrada cooler
%   P_CO2           Presión CO2
%   m_CO2           Flujo másico CO2
%   apertura_pct    (escala el área "efectiva" del cooler)

% SALIDAS
%   A_cooler        área efectiva (A_nominal * apertura)
%   T6AO            temp CO2 salida cooler (punto 6AO)
%   Q_cooler        calor removido del CO2 (positivo si enfría CO2)
%   Xd_cooler       exergía destruida (incluye ExDP de Cooler())
%   TPhsmx_cooler   Matriz resumen de estados del ciclo: T, P, h, s, Ex y flujo másico en cada punto

%Entradas COOLER----------------------------------------------------------------------------------------------------------------------
    
    m_air_nominal = 1.16;   
    phi = max(0, min(100, apertura_pct))/100;
    m_air_equiv = phi * m_air_nominal;
    
    T0 = 25 + 273.15;        % Temperatura ambiente (K)
    P0 = 101.325;            % Presión atmosférica (kPa)     

    % Estado de referencia para exergía (estado muerto) aire
        [h0a,s0a] = refpropm('HS','T',25+273.15,'P',101.3,'air.ppf');   
    % Estado de referencia para exergía (estado muerto) CO2
        [h0c,s0c] = refpropm('HS','T',25+273.15,'P',101.3,'co2');
% Puntos --------------------------------------------------------------------------------------------------------------------------------------------------
     % Punto 13
        T13 = T0;
        P13 = P0;
        [h13,s13] = refpropm('HS','T',T0,'P',P0,'air.ppf');
        x13 = (h13-h0a)-T0*(s13-s0a);
        
    % Punto 14
        T14 = T0-20;
        P14 = P0;
        [h14,s14] = refpropm('HS','T',T6AO3-20,'P',P0,'air.ppf');
        x14 = (h14-h0a)-T0*(s14-s0a);
        
    % Punto 6AO3
        [h6AO3, s6AO3] = refpropm('HS','T',T6AO3,'P',P_CO2,'co2');
        x6AO3 = (h6AO3-h0c) - T0*(s6AO3-s0c);
       

% Catálogo HX --------------------------------------------------------------------------------------------
    HX = [460,160,336,125,5;
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
% --------------------------------------------------------------------------------------------------------------------------------------------------
    if phi == 0
        A_Cooler = 0;
        T6AO = T6AO3;
        P6AO = P_CO2;
        Q_cooler = 0;
        Xd_cooler = 0;
        h6AO = h6AO3; s6AO = s6AO3; x6AO = x6AO3;
        T  = [T13, T14, T6AO3, T6AO];  % K
        P  = [P13, P14, P_CO2, P6AO];  % kPa
        H  = [h13, h14, h6AO3, h6AO];  % J/kg
        S  = [s13, s14, s6AO3, s6AO];  % J/kg-k
        m  = [m_CO2, m_CO2, 0, 0]; % kg/s
        x  = [x13, x14, x6AO3, x6AO];  % J/kg
        x_loss_cooler = 0;
        TPhsmx_cooler = [T', P', (H./1000)', (S./1000)', m',(x./1000)'];
        return
    end
    
%---------------------------------------------------------------------
% ANÁLISIS DEL COOLER (Puntos: 6AO3, 6AO, 13, 14)
    % Calor del cooler
        Q_cooler = m_air_equiv*(h14 - h13); 

    % Punto 6AO
        h6AO = h6AO3 - Q_cooler/m_CO2;
        P6AO = P_CO2;
        [T6AO,s6AO] = refpropm('TS','H',h6AO,'P',P6AO,'co2');
        x6AO = (h6AO-h0c)-T0*(s6AO-s0c);
        
    % Exergía destruida en el cooler    
        Xd_cooler = m_CO2*(x6AO3-x6AO)+m_air_equiv*(x13-x14);
        x_loss_cooler = -m_air_equiv*(x13-x14);
               


    for ii = 1:size(HX,1)
        Hr  = HX(ii,1);
        Wr  = HX(ii,2);
        Lpr = HX(ii,3);

        varout = Cooler(Lpr,Wr,Hr, T14,T13, T6AO3,T6AO, P0,P_CO2, Q_cooler, 'co2', m_air_equiv, m_CO2);
        A_Cooler=varout(2);
        if ii==1
            Amin_Cooler=A_Cooler;
        else
            if (A_Cooler<Amin_Cooler)&&(varout(3)<HX(ii,4))&&(A_Cooler<HX(ii,5))
                Amin_Cooler=A_Cooler;
            end
        end
    end

%---------------------------------------------------------------------
% RESULTADOS
        T  = [T13, T14, T6AO3, T6AO];  % K
        P  = [P13, P14, P_CO2, P6AO];  % Pa
        H  = [h13, h14, h6AO3, h6AO];  % J/kg
        S  = [s13, s14, s6AO3, s6AO];  % J/kg-k
        m  = [m_air_equiv, m_air_equiv, m_CO2, m_CO2]; % kg/s
        x  = [x13, x14, x6AO3, x6AO];  % J/kg

        TPhsmx_cooler = [T', (P./1000)', (H./1000)', (S./1000)', m',(x./1000)'];
end

        
        
        
        