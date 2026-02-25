function [TPhsmx_PEM, Xd_PEM_TOT, A_HE, Q_HE, J, V, m_H2O_in, m_H2] = PEM_OFICIAL(T6AO, P_CO2, m_CO2, T7, W_PEM)
format long

% ENTRADAS
    %   T_PEM= Temperatura de la celda
    %   W_PEM = Potencia eléctrica entregada a la celda
    %   T7 = Temperatura de 7
    %   P_CO2 = Presión de la línea caliente (Brayton)
    %   m_CO2 = Flujo másico de la línea caliente
    
% SALIDAS
    %   J= Densidad de corriente [A/m^2]
    %   V= Voltaje de celda [V]
    %   mdotH2= Flujo másico de H2 producido [kg/s]
    %   mdotH2O_in= Flujo másico mínimo de H2O consumido [kg/s]

        % CONSTANTES FÍSICAS
            R = 8.314;             % J/mol-K
            F = 96486;             % C/mol
            M_H2  = 2.01588e-3;    % kg/mol
            M_O = 31.998e-3;       % kg/mol
            M_H2O = 18.01528e-3;   % kg/mol
            x_ch_H2 = 236.09;      % kJ/mol
            x_ch_O = 3.97;         % kJ/mol
            x_ch_H2O = 0.9;       % kJ/mol
        % Entradas del ciclo combinado a la PEM
            T0 = 25 + 273.15;      % Temp low (K)--> Guajira
            fluido = 'water';
            P0 = 101.325;           % Presión atmosférica (kPa)
            
        % Estado de referencia para exergía (estado muerto) CO2
            [h0c,s0c] = refpropm('HS','T',25+273.15,'P',101.3,'co2');
        % Estado de referencia para exergía (estado muerto) agua
            [h0w,s0w] = refpropm('HS','T',25+273.15,'P',101.3,'water');
            
        %PARÁMETROS DEL PEM 
             A_PEM = 1;     % m^2 1*10-3
             L_mem = 50e-6;     % m 127

            % Hidratación 
                 lambda_an = 14;
                 lambda_ca = 10;

            Eact_an = 76000;       % K/mol
            Eact_ca = 18000;       % K/mol
            Jref_an = 170000;      % A/m^2
            Jref_ca = 4600;        % A/m^2
        

% Puntos ------------------------------------------------------------------------------------------------
        % Punto 7 Brayton
            P7 = P_CO2;
            [h7,s7] = refpropm('HS','T',T7,'P',P7,'co2');
            x7 = (h7-h0c)-T0*(s7-s0c);
            
        % Punto 6AO (entrada CO2 a HE)
            P6AO = P_CO2;
            [h6AO, s6AO] = refpropm('HS','T',T6AO,'P',P6AO,'co2');
            x6AO = (h6AO - h0c) - T0*(s6AO - s0c);

        % Punto 19 Agua
            T19 = T0; P19 = P0;
            [h19,s19] = refpropm('HS','T',T19,'P',P19,fluido);
            x19 = (h19-h0w)-T0*(s19-s0w);
        
% --------------------------------------------------------------------------------------------------------------------------------------------------
% Catálogo HX --------------------------------------------------------------------------------------------
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


 %-------------------------------------------------------------------------
    % ANÁLISIS DEL HEAT EXCHANGER (Puntos: 6AO, 7, 19, 20)
        % Calor del HEAT EXCHANGER 
            Q_HE = m_CO2*(h7-h6AO);
            
            T20 = max(T19+1.0, min(T6AO-1.0,T7+20.0));
            tolT = 1e-6;
            maxit = 200;
            
        % Cálculo densidad de corriente y voltaje/PRODUCCIÓN de H2 y H2O --
            for ii = 1:maxit
                
                [g, J, V, m_H2, m_H2O_in, m_O] = balance_T20(T20);

                if abs(g) < tolT
                    break
                end

                % derivada numérica central
                dT = max(1e-3, 1e-4*abs(T20));
                gp = funtemp(T20 + dT);
                gm = funtemp(T20 - dT);
                dg = (gp - gm)/(2*dT);

                % protección dg ~ 0
                if ~isfinite(dg) || abs(dg) < 1e-12
                    % no hay pendiente confiable -> sal o fuerza un paso pequeño
                    break
                end

                Tnew = T20 - g/dg;

                Tnew = max(T19 + 1.0, min(T6AO - 1.0, Tnew));

                T20 = Tnew;

            end

        % Punto 20 final (agua a PEM)
            P20 = P0;
            [h20, s20] = refpropm('HS','T',T20,'P',P20,fluido);
            x20 = (h20 - h0w) - T0*(s20 - s0w);
            
        % Exergía destruida en el HEAT EXCHANGER
            Xd_HE = m_CO2*(x6AO-x7)+m_H2O_in*(x19-x20);  
    
    
% --------------------------------------------------------------------------------------------------------------------------------------------------

Heat = [-Q_HE; 0; 0];

for ii=1:26
    Hhx = HX(ii,1); Whx = HX(ii,2); Lphx = HX(ii,3);

    varout = Evaporador(Lphx,Whx,Hhx, ...
        T20, (0.5*(T19 + T20)), T19, ...
        T6AO, (0.5*(T6AO + T7)), T7, T7, ...
        P0, m_H2O_in, m_CO2, Heat, fluido, P_CO2);

    Atest = varout(6);

    if ii==1
        Amin = Atest;
    else
        if (Atest < Amin) && (varout(7) < HX(ii,4)) && (Atest < HX(ii,5))
            Amin = Atest;
        end
    end
end

A_HE = Amin;
%--------------------------------------------------------------------------
xH2ch = x_ch_H2/M_H2;              % kJ/kg
XH2 = m_H2*xH2ch;                   % kJ/s

XH2O_in_ph = m_H2O_in*x20/1000;
XH2O_in_ch = m_H2O_in*x_ch_H2O/M_H2O;

xOch = x_ch_O/M_O;              % kJ/kg
XO = m_O*xOch;               % kJ/s

Xd_PEM = W_PEM/1000+XH2O_in_ph+XH2O_in_ch-(XH2+XO);
Xd_PEM_TOT = [Xd_HE,XH2,XH2O_in_ph,Xd_PEM];

% RESULTADOS
T = [T6AO, T7,  T19, T20];      % K
P = [P6AO, P7,  P19, P20];      % kPa
H = [h6AO, h7,  h19, h20];      % J/kg
S = [s6AO, s7,  s19, s20];      % J/kg-K
m = [m_CO2, m_CO2, m_H2O_in, m_H2O_in]; % kg/s
x = [x6AO, x7,  x19, x20];      % J/kg

TPhsmx_PEM = [T', P', (H./1000)', (S./1000)', m', (x./1000)'];

%--------------------------------------------------------------------------------------------------------------------------------------------------
% Funciones ---------------------------------------------------------------
    function [g, Jloc, Vloc, mH2loc, mH2Oloc, mOloc] = balance_T20(Tk)

        % PEM -> J,V -> mH2, mH2O
        [Jloc, Vloc] = solv_potencia_variable(Tk, W_PEM, ...
            A_PEM, L_mem, lambda_an, lambda_ca, ...
            Jref_an, Jref_ca, Eact_an, Eact_ca, R, F);

        Iloc = Jloc*A_PEM;
        N_H2 = Iloc/(2*F);
        N_H2O = Iloc/(2*F);
        N_O = Iloc/(4*F);
        mH2loc  = N_H2*M_H2;
        mH2Oloc = N_H2O*M_H2O;
        mOloc = N_O*M_O;
        
        % Balance HE
        h20loc = refpropm('H','T',Tk,'P',P0,'water');
        g = mH2Oloc*(h20loc - h19) + Q_HE;
    end

    function g = funtemp(Tk)
        [g,~,~,~,~] = balance_T20(Tk);
    end

    function Vcell = Voltaje_PEM(Tk, Jk, Lm, lam_an, lam_ca, ...
                                 Jr_an, Jr_ca, Ea_an, Ea_ca, R, F)
        % Voltaje Reversible 
        V0 = 1.229 - 8.5e-4*(Tk - 298.15);

        % Corrientes de intercambio
        J0_an = Jr_an * exp(-Ea_an/(R*Tk));
        J0_ca = Jr_ca * exp(-Ea_ca/(R*Tk));

        % Resistencia de membrana
        R_PEM = ResArea_PEM(Tk, Lm, lam_an, lam_ca);

        % Sobrepotenciales
        K = R*Tk/F;
        V_act_anca = K*(asinh(Jk/(2*J0_ca)) + asinh(Jk/(2*J0_an)) );
        V_ohm = Jk * R_PEM;

        Vcell = V0 + V_act_anca + V_ohm;
    end

    function R_PEM = ResArea_PEM(Tk, Lm, lam_an, lam_ca)
        n = 1001; % impar para Simpson 1/3
        x = linspace(0, Lm, n);

        lam_x = ((lam_an - lam_ca)/Lm).*x + lam_ca;

        sigma = (0.5139.*lam_x - 0.326) .* exp(1268.*(1/303 - 1./Tk));
        sigma = max(sigma, 1e-12);

        R_PEM = simpson13(x, 1./sigma);
    end

    function Int = simpson13(x, y)
        n = numel(x);
        h = (x(end)-x(1))/(n-1);
        Int = y(1) + y(end) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2));
        Int = Int*h/3;
    end

    function [Jsol, Vsol] = solv_potencia_variable(Tk, Wk, A, Lm, lam_an, lam_ca, ...
                                                   Jr_an, Jr_ca, Ea_an, Ea_ca, Rg, Fc)
        J_in = 1000;     
        tol  = 1e-6;

        f  = @(Jk) A*Jk*Voltaje_PEM(Tk, Jk, Lm, lam_an, lam_ca, ...
                                    Jr_an, Jr_ca, Ea_an, Ea_ca, Rg, Fc) - Wk;

        df = @(Jk) (f(Jk + max(1e-6,1e-4*Jk)) - f(Jk - max(1e-6,1e-4*Jk))) ...
                   / (2*max(1e-6,1e-4*Jk));
        % Newton-Rapshon
        for i = 1:200
            J_new = J_in - f(J_in)/df(J_in);

            if abs(f(J_new)) < tol
                Jsol = J_new;
                Vsol = Wk/(A*Jsol);
                return;
            end

            J_in = J_new;
        end
    end
end


