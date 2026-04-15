function [T6AO3, A_ORC, W_t3, W_b2, Xd_orc, TPhsmx_orc, Q_cond_vec, Q_evap_vec, Q_regen] = ORCV2(Ts, P_CO2, m_CO2, cycle_type, ehr)
% ORCV2 adaptado para ORC simple y RORC
%
% ENTRADAS:
%   Ts         = Temperatura de la fuente caliente [K]
%   P_CO2      = Presion de la linea caliente (Brayton/CO2) [kPa]
%   m_CO2      = Flujo masico de la linea caliente [kg/s]
%   cycle_type = 'ORC' o 'RORC'  (opcional, por defecto 'ORC')
%   ehr        = Efectividad del regenerador (0-1) (opcional, por defecto 0.85)
%
% SALIDAS:
%   T6AO3      = Temperatura de salida de la linea caliente [K]
%   A_ORC      = [Aevap; Aregen; Acond]
%   W_t3       = Potencia generada por la turbina [W]
%   W_b2       = Potencia consumida por la bomba [W] (negativa)
%   Xd_orc     = [Xd_bomba, Xd_turbina, Xd_regenerador, Xd_condensador, Xd_evaporador, Ex_rechazada_agua]
%   TPhsmx_orc = Matriz de estados: T, P, h, s, m, ex
%   Q_cond_vec = [Q_cond_1; Q_cond_2; Q_cond_total]
%   Q_evap_vec = [Q_evap_1; Q_evap_2; Q_evap_3; Q_evap_total]
%   Q_regen    = Calor recuperado en el regenerador [W]

    if nargin < 4 || isempty(cycle_type)
        cycle_type = 'ORC';
    end
    if nargin < 5 || isempty(ehr)
        ehr = 0.85;
    end
    cycle_type = upper(string(cycle_type));
    use_regen = (cycle_type == "RORC");

    %----------------------------------------------------------------------
    % ENTRADAS / PARAMETROS DEL SISTEMA ORC
    %----------------------------------------------------------------------
    fluid = 'cyclohex';
    np = 0.85;
    nt = 0.85;
    Ap = 39;
    pinch = 10.0;
    Tc = 40 + 273.15;
    T0 = 25 + 273.15;
    r_ORC = 4;

    Plow  = refpropm('P','T',Tc,'Q',0,fluid);
    Phigh = Plow*r_ORC;
    Pcrit = refpropm('P','C',0,' ',0,fluid);
    P0 = 101.325;

    % Estados muertos para exergia
    [h0t,s0t] = refpropm('HS','T',25+273.15,'P',101.3,fluid);
    [h0c,s0c] = refpropm('HS','T',25+273.15,'P',101.3,'co2');
    [h0w,s0w] = refpropm('HS','T',25+273.15,'P',101.3,'water');

    % Inicializacion por seguridad
    T6AO3 = 0;
    A_ORC = [0;0;0];
    W_t3 = 0;
    W_b2 = 0;
    Xd_orc = zeros(1,6);
    Q_cond_vec = zeros(3,1);
    Q_evap_vec = zeros(4,1);
    Q_regen = 0;
    TPhsmx_orc = zeros(16,6);

    if Phigh >= Pcrit
        return
    end

    %----------------------------------------------------------------------
    % PUNTOS ORC
    %----------------------------------------------------------------------
    % Punto 17: salida de condensador / entrada a bomba
    T17 = Tc;
    P17 = Plow;
    [h17, s17] = refpropm('HS','T',T17,'P',P17,fluid);
    x17 = (h17-h0t)-T0*(s17-s0t);

    % Punto 18: salida de bomba
    P18 = Phigh;
    s18s = s17;
    h18s = refpropm('H','P',P18,'S',s18s,fluid);
    h18  = h17 + (h18s-h17)/np;
    [T18,s18] = refpropm('TS','P',P18,'H',h18,fluid);
    x18 = (h18-h0t)-T0*(s18-s0t);

    % Punto 15ls: liquido saturado a alta presion
    P15ls = P18;
    [T15ls,h15ls] = refpropm('TH','P',P15ls,'Q',0,fluid);
    s15ls = refpropm('S','P',Phigh,'Q',0,fluid);
    x15ls = (h15ls-h0t)-T0*(s15ls-s0t);

    % Punto 15vs: vapor saturado a alta presion
    T15vs = T15ls;
    P15vs = P15ls;
    [h15vs,s15vs] = refpropm('HS','P',P15vs,'Q',1,fluid);
    x15vs = (h15vs-h0t)-T0*(s15vs-s0t);

    % Punto 15: salida de evaporador / entrada a turbina
    P15 = P15vs;
    T15 = Ts - Ap;
    [h15,s15] = refpropm('HS','T',T15,'P',P15,fluid);
    x15 = (h15-h0t)-T0*(s15-s0t);

    if h15 < h15vs
        return
    end

    % Punto 16: salida de turbina
    P16 = Plow;
    h16s = refpropm('H','P',P16,'S',s15,fluid);
    h16 = h15-(h15-h16s)*nt;
    [T16,s16] = refpropm('TS','P',P16,'H',h16,fluid);
    x16 = (h16-h0t)-T0*(s16-s0t);

    % Punto 16vs: vapor saturado a baja presion
    T16vs = T17;
    P16vs = P17;
    [h16vs,s16vs] = refpropm('HS','T',T16vs,'Q',1,fluid);
    x16vs = (h16vs-h0t)-T0*(s16vs-s0t);

    %----------------------------------------------------------------------
    % REGENERADOR (si aplica)
    % Estados:
    %   18r = salida lado frio del regenerador (antes del evaporador)
    %   16r = salida lado caliente del regenerador (antes del condensador)
    %----------------------------------------------------------------------
    h18r = h18; T18r = T18; s18r = s18; x18r = x18;
    h16r = h16; T16r = T16; s16r = s16; x16r = x16;
    Aregen = 0;
    Xd_regenerador = 0;
    Q_regen = 0;

    if use_regen && (ehr > 0) && (T16 > T18)
        h18r_max = refpropm('H','T',T16,'P',Phigh,fluid);
        h16r_min = refpropm('H','T',T18,'P',Plow,fluid);

        Qmax1 = h18r_max - h18;
        Qmax2 = h16 - h16r_min;
        Qmax = min(Qmax1,Qmax2);

        if Qmax > 0
            Q_regen = ehr * Qmax;
            h18r = h18 + Q_regen;
            h16r = h16 - Q_regen;

            [T18r,s18r] = refpropm('TS','P',Phigh,'H',h18r,fluid);
            [T16r,s16r] = refpropm('TS','P',Plow,'H',h16r,fluid);

            x18r = (h18r-h0t)-T0*(s18r-s0t);
            x16r = (h16r-h0t)-T0*(s16r-s0t);
        end
    end

    %----------------------------------------------------------------------
    % BRAYTON (CO2) Y CONDENSADOR
    %----------------------------------------------------------------------
    % Punto 6AI
    P6AI = P_CO2;
    T6AI = Ts;
    [h6AI,s6AI] = refpropm('HS','T',Ts,'P',P6AI,'co2');
    x6AI = (h6AI-h0c) - T0*(s6AI-s0c);

    % Punto 6AO2
    T6AO2 = T15ls + pinch;
    P6AO2 = P_CO2;
    [h6AO2,s6AO2] = refpropm('HS','T',T6AO2,'P',P6AO2,'co2');
    x6AO2 = (h6AO2-h0c)-T0*(s6AO2-s0c);

    % Agua de enfriamiento
    TC1 = T0;  PC1 = P0;
    [hC1,sC1] = refpropm('HS','T',TC1,'P',PC1,'water');
    xC1 = (hC1-h0w)-T0*(sC1-s0w);

    PC2 = P0;
    TC2 = T17 - pinch;
    [hC2,sC2] = refpropm('HS','T',TC2,'P',PC2,'water');
    xC2 = (hC2-h0w)-T0*(sC2-s0w);

    %----------------------------------------------------------------------
    % FLUJOS MASICOS
    % En RORC el evaporador arranca en 18r; en ORC arranca en 18
    %----------------------------------------------------------------------
    T6AOCp = 0.5*(T6AO2 + Ts);
    CpHS = refpropm('C','T',T6AOCp,'P',P_CO2,'co2');

    m_ORC = m_CO2*CpHS*(Ts-T6AO2)/(h15-h15ls);

    h17p = refpropm('H','T',T17,'Q',1,fluid);
    if h16r < h17p
        m_water = m_ORC*(h16r-h17)/(hC2-hC1);
    else
        m_water = m_ORC*(h17p-h17)/(hC2-hC1);
    end

    %----------------------------------------------------------------------
    % PUNTOS RESTANTES
    %----------------------------------------------------------------------
    PC3 = P0;
    hC3 = hC1 + m_ORC*(h16r-h17)/m_water;
    [TC3,sC3] = refpropm('TS','P',PC3,'H',hC3,'water');
    xC3 = (hC3-h0w)-T0*(sC3-s0w);
    xloss_cond = (xC1-xC3)*m_water;

    % Punto 6AO1
    T6AO1 = Ts - m_ORC*(h15-h15vs)/(m_CO2*CpHS);
    P6AO1 = P_CO2;
    [h6AO1,s6AO1] = refpropm('HS','T',T6AO1,'P',P6AO1,'co2');
    x6AO1 = (h6AO1-h0c)-T0*(s6AO1-s0c);

    % Punto 6AO3: salida total de la fuente caliente
    T6AO3 = Ts - m_ORC*(h15-h18r)/(m_CO2*CpHS);
    P6AO3 = P_CO2;
    [h6AO3,s6AO3] = refpropm('HS','T',T6AO3,'P',P6AO3,'co2');
    x6AO3 = (h6AO3-h0c)-T0*(s6AO3-s0c);

    %----------------------------------------------------------------------
    % ANALISIS DE COMPONENTES
    %----------------------------------------------------------------------
    % Bomba
    W_b2 = m_ORC*(h17-h18);   % negativa: consumo
    Xd_b2 = -W_b2 + m_ORC*(x17-x18);

    % Turbina
    W_t3 = m_ORC*(h15-h16);
    Xd_t3 = -W_t3 + m_ORC*(x15-x16);

    % Regenerador
    if use_regen && (Q_regen > 0)
        Xd_regenerador = m_ORC*((x16-x16r) - (x18r-x18));
    else
        Xd_regenerador = 0;
    end

    % Condensador
    Q_condensador_e1 = m_water*(hC2-hC1);
    if h16r < h17p
        Q_condensador_e2 = 0;
    else
        Q_condensador_e2 = m_water*(hC3-hC2);
    end
    Xd_condensador = m_ORC*(x16r-x17) + m_water*(xC1-xC3);

    % Evaporador
    Q_evaporador_e1 = m_ORC*(h15ls-h18r);
    Q_evaporador_e2 = m_ORC*(h15vs-h15ls);
    Q_evaporador_e3 = m_ORC*(h15-h15vs);
    Q_evaporador    = m_ORC*(h15-h18r);

    Xd_evaporador = m_CO2*(x6AI-x6AO3) + m_ORC*(x18r-x15);

    %----------------------------------------------------------------------
    % DIMENSIONAMIENTO DE INTERCAMBIADORES
    %----------------------------------------------------------------------
    Heat = zeros(1,3);
    Heat(1) = Q_evaporador_e1;
    Heat(2) = Q_evaporador_e2;
    Heat(3) = Q_evaporador_e3;

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

    %-------------------------
    % Evaporador
    %-------------------------
    Aevap = 0;
    Evap = [];
    iioptevap = 0;
    Amin = inf;

    for ii = 1:26
        H = HX(ii,1); W = HX(ii,2); Lp = HX(ii,3);
        varout = Evaporador(Lp, W, H, T15, T15ls, T18r, Ts, T6AO1, T6AO2, T6AO3, Phigh, m_ORC, m_CO2, Heat, fluid, P_CO2);
        Atest = varout(6);

        if ii == 1 || ((Atest < Amin) && (varout(7) < HX(ii,4)) && (Atest < HX(ii,5)))
            Amin = Atest;
            Aevap = Atest;
            Evap = varout;
            iioptevap = ii;
        end
    end

    %-------------------------
    % Condensador
    %-------------------------
    Qcond = Q_condensador_e1;
    Qenf  = Q_condensador_e2;
    x16q  = refpropm('Q','P',P16,'H',h16r,fluid);

    Acond = 0;
    Cond = [];
    iioptcond = 0;
    Amin = inf;

    for ii = 1:26
        Hc = HX(ii,1); Wc = HX(ii,2); Lpc = HX(ii,3);
        varout = Condensador(Lpc, Wc, Hc, Qcond, Qenf, T16r, T17, T17, TC1, TC2, TC3, m_ORC, m_water, Plow, x16q, fluid);
        Atest = varout(1) + varout(2);

        if ii == 1 || ((Atest < Amin) && (varout(3) < HX(ii,4)) && (Atest < HX(ii,5)))
            Amin = Atest;
            Acond = Atest;
            Cond = varout;
            iioptcond = ii;
        end
    end

    %-------------------------
    % Regenerador
    %-------------------------
    Regen = zeros(1,4);
    iioptreg = 0;

    if use_regen && (Q_regen > 0)
        Amin = inf;
        for ii = 1:26
            Hr = HX(ii,1); Wr = HX(ii,2); Lpr = HX(ii,3);
            varout = Regenerador(Lpr, Wr, Hr, T18r, T18, T16, T16r, Phigh, Plow, m_ORC, Q_regen, fluid);
            Atest = varout(2);

            if ii == 1 || ((Atest < Amin) && (varout(3) < HX(ii,4)) && (Atest < HX(ii,5)))
                Amin = Atest;
                Aregen = Atest;
                Regen = varout;
                iioptreg = ii;
            end
        end
    else
        Aregen = 0;
        Xd_regenerador = 0;
        Q_regen = 0;
    end

    A_ORC = [Aevap; Aregen; Acond];

    %----------------------------------------------------------------------
    % RESULTADOS DE EXERGIA
    %----------------------------------------------------------------------
    Xd_orc = [Xd_b2, Xd_t3, Xd_regenerador, Xd_condensador, Xd_evaporador, -xloss_cond];

    %----------------------------------------------------------------------
    % MATRIZ DE ESTADOS
    % Orden:
    % 17, 18, 18r, 15ls, 15vs, 15, 16, 16r, 16vs, C1, C2, C3, 6AI, 6AO1, 6AO2, 6AO3
    %----------------------------------------------------------------------
    T = [T17, T18, T18r, T15ls, T15vs, T15, T16, T16r, T16vs, TC1, TC2, TC3, T6AI, T6AO1, T6AO2, T6AO3];
    P = [P17, P18, P18,  P15ls, P15vs, P15, P16, P16,  P16vs, PC1, PC2, PC3, P6AI, P6AO1, P6AO2, P6AO3];
    H = [h17, h18, h18r, h15ls, h15vs, h15, h16, h16r, h16vs, hC1, hC2, hC3, h6AI, h6AO1, h6AO2, h6AO3];
    S = [s17, s18, s18r, s15ls, s15vs, s15, s16, s16r, s16vs, sC1, sC2, sC3, s6AI, s6AO1, s6AO2, s6AO3];
    m = [m_ORC, m_ORC, m_ORC, m_ORC, m_ORC, m_ORC, m_ORC, m_ORC, m_ORC, m_water, m_water, m_water, m_CO2, m_CO2, m_CO2, m_CO2];
    x = [x17, x18, x18r, x15ls, x15vs, x15, x16, x16r, x16vs, xC1, xC2, xC3, x6AI, x6AO1, x6AO2, x6AO3];

    TPhsmx_orc = [T', P', (H./1000)', (S./1000)', m', (x./1000)'];

    %----------------------------------------------------------------------
    % VECTORES DE CALOR
    %----------------------------------------------------------------------
    Q_evap_vec = [Q_evaporador_e1; Q_evaporador_e2; Q_evaporador_e3; Q_evaporador];
    Q_cond_vec = [Q_condensador_e1; Q_condensador_e2; (Q_condensador_e1 + Q_condensador_e2)];

    %----------------------------------------------------------------------
    % VALIDACIONES BASICAS
    %----------------------------------------------------------------------
    Wnet = W_t3 + W_b2;
    Qin = Q_evaporador;

    if (m_ORC < 0) || (Wnet < 0) || (Qin < 0) || (T6AO3 < T18r)
        T6AO3 = 0;
        A_ORC = [0;0;0];
        W_t3 = 0;
        W_b2 = 0;
        Xd_orc = zeros(1,6);
        Q_cond_vec = zeros(3,1);
        Q_evap_vec = zeros(4,1);
        Q_regen = 0;
        TPhsmx_orc = zeros(16,6);
        return
    end
end