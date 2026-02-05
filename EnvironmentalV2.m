function [Ytotal,Yco_total,Yde_total,Yom_total,Y_componentes,Ytotal_PEM] = EnvironmentalV2(W_t1,W_t2,W_t3,W_c1,W_c2,W_b2,A_HTR,A_LTR,A_cooler,A_evap,A_cond,A_HE)
format long
% %ENTRADAS A LA FUNCIÓN-------------------------
%   W_t1,W_t2,W_t3                          = Potencias de turbinas      (kW) 
%   W_c1,W_c2                               = Potencias de compresores   (kW)
%   W_b2                                    = Potencia de bomba          (kW)
%   A_HTR,A_LTR,Acooler,Aevap,Acond,Aregen  = Areas de intercambiadores (m^2)

% SALIDAS DE LA FUNCIÓN-------------------------
%   Ytotal                                  = Impacto total normalizado [kgCO2-eq/kWh]



%--------------------------------------------------------------------------
% PARAMETROS GENERALES ----------------------------------------------------
    years           = 20;        % años
    disponibilidad  = 0.85;      % factor de disponibilidad
    horas           = 24*365;    % h/año
    peed            = 0.05;      % perdida en decomisionamiento (5%)

% Material de intercambiadores 
    rho   = 8900;                % kg/m3 (Cobre:rho=8900, Aluminio:rho=2698.4, Acero:rho=7930)
    delta = 0.002;               % m     Espesor intercambiador
    mu_eq = 2.20;                % kgCO2-eq/kWh (Cobre:mu_eq=2.20, Aluminio:mu_eq=2.9, Acero:mu_eq=1.77)
%--------------------------------------------------------------------------
% Potencia neta y energía de vida útil (kWh)
    WNETO = W_t1+W_t2+W_c1+W_c2;          % kW (IMPORTANTE: en kW)
    E_util = WNETO*horas*years*disponibilidad;      % kWh Energía de vida útil
    
% Acople PEM
[Ytotal_PEM,Yco_total_PEM,Yde_total_PEM,Yom_total_PEM,~] = EnvironmentalPEM_V2(W_t3,W_b2);
%--------------------------------------------------------------------------
% Factores de impacto equipos (mu_eq COBRE)
    mu_Evap   = mu_eq;       % kgCO2-eq/kWh
    mu_Cond   = mu_eq;       % kgCO2-eq/kWh
    mu_HTR    = mu_eq;       % kgCO2-eq/kWh
    mu_LTR    = mu_eq;       % kgCO2-eq/kWh
    mu_Cooler = mu_eq;       % kgCO2-eq/kWh
    mu_HE     = mu_eq;       % kgCO2-eq/kWh

    mu_T1     = mu_eq;       % kgCO2-eq/kWh
    mu_T2     = mu_eq;       % kgCO2-eq/kWh
    mu_T3     = mu_eq;       % kgCO2-eq/kWh
    mu_C1     = mu_eq;       % kgCO2-eq/kWh
    mu_C2     = mu_eq;       % kgCO2-eq/kWh
    mu_b2     = mu_eq;       % kgCO2-eq/kWh

    mu_CSP    = 9.456e-3;    % kgCO2-eq/kWh

% Fluidos
    mu_tolueno  = 32.3;      % kgCO2-eq/kWh
    mu_CO2      = 1.0;       % kgCO2-eq/kWh
    
%--------------------------------------------------------------------------
% Calidad para equipos 
    a_turb = 35.03;  % kg/kW (turbinas)
    a_pump = 15.71;  % kg/kW (bombas/compresores)
    
%--------------------------------------------------------------------------
% Masa equivalente equipos (M_eq)
    M_Evap   = A_evap*delta*rho;        % kg
    M_Cond   = A_cond*delta*rho;        % kg
    M_HTR    = A_HTR*delta*rho;         % kg
    M_LTR    = A_LTR*delta*rho;         % kg
    M_Cooler = A_cooler*delta*rho;      % kg
    M_HE     = A_HE*delta*rho;          % kg

    M_T1   = a_turb*W_t1;               % kg
    M_T2   = a_turb*W_t2;               % kg
    M_T3   = a_turb*W_t3;               % kg

    M_C1   = a_pump*W_c1;               % kg
    M_C2   = a_pump*W_c2;               % kg
    M_B2   = a_pump*W_b2;               % kg

% Masa equivalente fluidos (M_eq)
    %Tolueno
        wf_c    = 5.57*W_t3;                % kg (carga inicial)
        wf_op   = wf_c*0.10;                % kg (reposicion O&M)
        wf_de   = (wf_c-wf_op)*0.03;        % kg (remanente a decommission)
    %CO2
        co2_c   = 5.57*(W_t1+W_t2);                         % kg (carga inicial)
        co2_op  = co2_c*0.10;                % kg (reposicion O&M)
        co2_de  = (co2_c-co2_op)*0.03;      % kg (remanente a decommission)

%--------------------------------------------------------------------------
% Decomisionamiento
    Yde.Evap   = peed*mu_Evap*M_Evap;
    Yde.Cond   = peed*mu_Cond*M_Cond;
    Yde.HTR    = peed*mu_HTR*M_HTR;
    Yde.LTR    = peed*mu_LTR*M_LTR;
    Yde.Cooler = peed*mu_Cooler*M_Cooler;
    Yde.HE     = peed*mu_HE*M_HE;

    Yde.T1     = peed*mu_T1*M_T1;
    Yde.T2     = peed*mu_T2*M_T2;
    Yde.T3     = peed*mu_T3*M_T3;

    Yde.C1     = peed*mu_C1*M_C1;
    Yde.C2     = peed*mu_C2*M_C2;
    Yde.B2     = peed*mu_b2*M_B2;

    Yde.wf     = mu_tolueno*wf_de;
    Yde.co2    = mu_CO2*co2_de;

    Yde.CSP    = peed*(E_util*mu_CSP);
   
    Yde_total  = Yde.Evap+Yde.HE+Yde.Cond+Yde.HTR+Yde.LTR+Yde.Cooler+Yde.T1...
                +Yde.T2+Yde.T3+Yde.C1+Yde.C2+Yde.B2+Yde.wf+Yde.co2+Yde.CSP+Yde_total_PEM; %
    
%--------------------------------------------------------------------------
% Comisionamiento
    Yco.Evap   = mu_Evap*M_Evap+Yde.Evap;
    Yco.Cond   = mu_Cond*M_Cond+Yde.Cond;
    Yco.HTR    = mu_HTR*M_HTR+Yde.HTR;
    Yco.LTR    = mu_LTR*M_LTR+Yde.LTR;
    Yco.Cooler = mu_Cooler*M_Cooler+Yde.Cooler;
    Yco.HE     = mu_HE*M_HE+Yde.HE;

    Yco.T1     = mu_T1*M_T1+Yde.T1;
    Yco.T2     = mu_T2*M_T2+Yde.T2;
    Yco.T3     = mu_T3*M_T3+Yde.T3;

    Yco.C1     = mu_C1*M_C1+Yde.C1;
    Yco.C2     = mu_C2*M_C2+Yde.C2;
    Yco.B2     = mu_b2*M_B2+Yde.B2;

    Yco.wf     = mu_tolueno*wf_c+Yde.wf;
    Yco.co2    = mu_CO2*co2_c+Yde.co2;

    Yco.CSP    = mu_CSP*E_util+Yde.CSP;
    
    Yco_total  = Yco.Evap+Yco.HE+Yco.Cond+Yco.HTR+Yco.LTR+Yco.Cooler+Yco.T1...
                +Yco.T2+Yco.T3+Yco.C1+Yco.C2+Yco.B2+Yco.wf+Yco.co2+Yco.CSP+Yco_total_PEM; %
%--------------------------------------------------------------------------
% Operación y mantenimiento
    Yom.Evap   = 0;
    Yom.Cond   = 0;
    Yom.HTR    = 0;
    Yom.LTR    = 0;
    Yom.Cooler = 0;
    Yom.HE  = 0;

    Yom.T1     = 0;
    Yom.T2     = 0;
    Yom.T3     = 0;

    Yom.C1     = 0;
    Yom.C2     = 0;
    Yom.B2     = 0;

    Yom.wf     = 0.10 * Yco.wf;
    Yom.co2    = 0.10 * Yco.co2;
    
    Yom.CSP    = 0;
    
    Yom_total = Yom.wf+Yom.co2+Yom_total_PEM;
%--------------------------------------------------------------------------
% Impacto total por componente
    Y.Evap   = Yco.Evap   + Yde.Evap   + Yom.Evap;
    Y.Cond   = Yco.Cond   + Yde.Cond   + Yom.Cond;
    Y.HTR    = Yco.HTR    + Yde.HTR    + Yom.HTR;
    Y.LTR    = Yco.LTR    + Yde.LTR    + Yom.LTR;
    Y.Cooler = Yco.Cooler + Yde.Cooler + Yom.Cooler;
    Y.HE     = Yco.HE  + Yde.HE  + Yom.HE;

    Y.T1     = Yco.T1 + Yde.T1 + Yom.T1;
    Y.T2     = Yco.T2 + Yde.T2 + Yom.T2;
    Y.T3     = Yco.T3 + Yde.T3 + Yom.T3;

    Y.C1     = Yco.C1 + Yde.C1 + Yom.C1;
    Y.C2     = Yco.C2 + Yde.C2 + Yom.C2;
    Y.B2     = Yco.B2 + Yde.B2 + Yom.B2;

    Y.wf     = Yco.wf + Yde.wf + Yom.wf;
    Y.co2    = Yco.co2 + Yde.co2 + Yom.co2;

    Y.CSP    = Yco.CSP + Yde.CSP + Yom.CSP;
    Y.PEM    = Yco_total_PEM + Yde_total_PEM + Yom_total_PEM; % kgCO2-eq (ABSOLUTO)
    
%--------------------------------------------------------------------------
% Resultados
    Y_sum = Y.Evap + Y.HE + Y.Cond + Y.HTR + Y.LTR + Y.Cooler + Y.CSP + ...
            Y.T1 + Y.T2 + Y.T3 + Y.C1 + Y.C2 + Y.B2 + Y.wf + Y.co2 + Y.PEM;%
    
    Ytotal = Y_sum / E_util; %(kgCO2-eq/kWh)
    
    Y_componentes = [Y.Evap,Y.Cond,Y.HTR,Y.LTR,Y.Cooler,Y.HE,Y.T1,Y.T2,Y.T3, ...
                 Y.C1,Y.C2,Y.B2,Y.wf,Y.co2,Y.CSP,Y.PEM];
end


function [Ytotal_PEM,Yco_total_PEM,Yde_total_PEM,Yom_total_PEM,Y_componentes_PEM] = EnvironmentalPEM_V2(W_t3,W_b2)
format long
% =========================================================================
% EnvironmentalPEM_V2
% Impacto ambiental (kgCO2-eq/kWh) de un PEM electrolyzer alimentado por
% el trabajo neto del ciclo Rankine (T3 + B2).
%
% IMPORTANTE:
% - El LCA del Rankine (electricidad) YA está calculado afuera.
% - Por lo tanto, aquí NO se contabiliza impacto por electricidad (mu_el = 0)
%   y se eliminan todos los términos Y*.El para evitar doble conteo.
%
% ESTRUCTURA (igual a EnvironmentalV2):
%   Parametros -> E_util -> factores -> masas -> Yde -> Yco -> Yom -> Y -> Ytotal
%
% ENTRADAS:
%   W_t3 = Potencia disponible desde Turbina 3 hacia PEM [kW]
%   W_b2 = Potencia asociada a Bomba 2 hacia PEM [kW]
%
% SALIDAS:
%   Ytotal = Impacto total normalizado [kgCO2-eq/kWh] (normalizado por E_util)
% =========================================================================
% REFERENCIAS (en comentarios tipo paper)
% - DOE (U.S. Dept. of Energy), PEM Electrolysis Technical Targets:
%   consumo típico de sistema ~55 kWh/kgH2 (estado ~2022) y metas.
% - Bareiß et al., Applied Energy (2019): LCA de PEMWE; dominancia de electricidad
%   cuando se incluye mix eléctrico (aquí se excluye por doble conteo).
% - Iyer et al., Int. J. Hydrogen Energy (2024): LCA de electrólisis; manufactura,
%   O&M y sensibilidad (electricidad usualmente domina si se incluye).
% - Estequiometría de electrólisis: ~9 kg H2O/kg H2 (mínimo).
% =========================================================================

%--------------------------------------------------------------------------
% PARAMETROS GENERALES ----------------------------------------------------
    years           = 20;        % años vida útil PEM
    disponibilidad  = 0.85;      % factor disponibilidad [-]
    horas           = 24*365;    % h/año
    peed            = 0.05;      % penalidad decomisionamiento (5%)

%--------------------------------------------------------------------------
% Energía de vida útil que alimenta a la PEM (Rankine -> PEM) --------------
% (se respeta el estilo de tu código: suma de potencias en kW)
    WNETO  = W_t3 + W_b2;                                % kW hacia PEM
    E_util = WNETO * horas * years * disponibilidad;     % kWh vida útil (base de normalización)

%--------------------------------------------------------------------------
% PARAMETROS DE PROCESO PEM (incluidos en el código) ----------------------
% Consumo eléctrico específico del sistema PEM [kWh/kgH2] (DOE targets)
    E_H2 = 55 ;     % kWh/kgH2  (ejemplo típico de sistema; ajustable)

% Agua mínima estequiométrica [kgH2O/kgH2]
    W_H2 = 9;      % kgH2O/kgH2 

% Impacto del agua (tratamiento) [kgCO2-eq/m3]
    mu_w  = 0.30;  % kgCO2-eq/m3 (ejemplo; ajusta a tu inventario)

% IMPACTO ELECTRICIDAD DESACTIVADO (Rankine ya calculado)
    mu_el = 0.0;   %#ok<NASGU>  % kgCO2-eq/kWh (no se usa)

%--------------------------------------------------------------------------
% FRACCIONES DE BOP (para repartir WNETO en equipos) -----------------------
% Estas fracciones son supuestos de ingeniería para asignar potencia a BOP.
% Ajusta con datos del proveedor o tu modelo de planta.
    alpha_rect   = 0.03;  % [-] rectificador / power electronics
    alpha_pump   = 0.02;  % [-] bombas
    alpha_comp   = 0.06;  % [-] compresor H2
    alpha_cooler = 0.02;  % [-] enfriamiento
    alpha_wt     = 0.01;  % [-] tratamiento agua

% Potencias por equipo (kW)
    W_rect   = alpha_rect   * WNETO;   % kW
    W_pump   = alpha_pump   * WNETO;   % kW
    W_comp   = alpha_comp   * WNETO;   % kW
    W_cooler = alpha_cooler * WNETO;   % kW
    W_wt     = alpha_wt     * WNETO;   % kW

% Potencia equivalente del stack (kW)
    W_stack = WNETO;                   % kW (stack equivalente)

%--------------------------------------------------------------------------
% Factores de impacto de materiales/equipos (kgCO2-eq/kg) ------------------
% Sustituir por base LCA (ecoinvent/GREET) o inventario propio.
    mu_stack_mat  = 20.0;   % kgCO2-eq/kg
    mu_rect_mat   = 8.0;    % kgCO2-eq/kg
    mu_pump_mat   = 6.0;    % kgCO2-eq/kg
    mu_comp_mat   = 7.0;    % kgCO2-eq/kg
    mu_cooler_mat = 5.0;    % kgCO2-eq/kg
    mu_wt_mat     = 6.0;    % kgCO2-eq/kg

%--------------------------------------------------------------------------
% Escalamiento de masa "a_*" (kg/kW) --------------------------------------
% Similar a a_turb y a_pump en tu código original.
    a_stack  = 12.0;   % kg/kW
    a_rect   = 3.0;    % kg/kW
    a_pump   = 1.0;    % kg/kW
    a_comp   = 2.5;    % kg/kW
    a_cooler = 1.5;    % kg/kW
    a_wt     = 0.8;    % kg/kW

%--------------------------------------------------------------------------
% Producción total de H2 y consumo de agua en vida útil --------------------
    H2_tot   = E_util / E_H2;          % kgH2 en la vida útil
    W_tot_kg = H2_tot * W_H2;          % kgH2O consumida en vida útil
    W_tot_m3 = W_tot_kg / 1000.0;      % m3

%--------------------------------------------------------------------------
% Masas equivalentes (M_*) ------------------------------------------------
    M_stack  = a_stack  * W_stack;     % kg
    M_rect   = a_rect   * W_rect;      % kg
    M_pump   = a_pump   * W_pump;      % kg
    M_comp   = a_comp   * W_comp;      % kg
    M_cooler = a_cooler * W_cooler;    % kg
    M_wt     = a_wt     * W_wt;        % kg

%--------------------------------------------------------------------------
% Decomisionamiento (Yde) -------------------------------------------------
% Operación: peed * (impacto equipo)
    Yde.Stack  = peed * mu_stack_mat  * M_stack;
    Yde.Rect   = peed * mu_rect_mat   * M_rect;
    Yde.Pump   = peed * mu_pump_mat   * M_pump;
    Yde.Comp   = peed * mu_comp_mat   * M_comp;
    Yde.Cooler = peed * mu_cooler_mat * M_cooler;
    Yde.WT     = peed * mu_wt_mat     * M_wt;

% Agua al final de vida (simple): fracción peed del impacto de agua total
    Yde.Water  = peed * (W_tot_m3 * mu_w);

% Total decomisionamiento
    Yde_total_PEM  = Yde.Stack + Yde.Rect + Yde.Pump + Yde.Comp + Yde.Cooler + ...
                 Yde.WT + Yde.Water;

%--------------------------------------------------------------------------
% Comisionamiento (Yco) ---------------------------------------------------
% Operación: fabricación del equipo + su término Yde (igual que tu patrón)
    Yco.Stack  = mu_stack_mat  * M_stack  + Yde.Stack;
    Yco.Rect   = mu_rect_mat   * M_rect   + Yde.Rect;
    Yco.Pump   = mu_pump_mat   * M_pump   + Yde.Pump;
    Yco.Comp   = mu_comp_mat   * M_comp   + Yde.Comp;
    Yco.Cooler = mu_cooler_mat * M_cooler + Yde.Cooler;
    Yco.WT     = mu_wt_mat     * M_wt     + Yde.WT;

% Puesta en marcha: no se añade carga inicial de agua (puedes activarla si quieres)
    Yco.Water  = 0 + Yde.Water;

% Total comisionamiento
    Yco_total_PEM  = Yco.Stack + Yco.Rect + Yco.Pump + Yco.Comp + Yco.Cooler + ...
                 Yco.WT + Yco.Water;

%--------------------------------------------------------------------------
% Operación y mantenimiento (Yom) ------------------------------------------
% NO se contabiliza electricidad (Rankine ya calculado).
% O&M incluye agua y reemplazos de stack (si aplica).

% Agua en operación
    Yom.Water  = W_tot_m3 * mu_w;

% Reemplazos de stack según vida (h)
% Referencia: durabilidad objetivo/valores reportados en DOE targets y literatura LCA
    L_stack = 60000;  % h vida stack (ejemplo típico; ajustar)
    h_tot   = horas * years * disponibilidad;      % h totales de operación
    N_stack = max(0, ceil(h_tot / L_stack) - 1);   % número de reemplazos

% Impacto por reemplazos (fabricación adicional del stack)
    Yom.Stack  = N_stack * (mu_stack_mat * M_stack);

% Otros equipos en O&M (simplificado = 0)
    Yom.Rect   = 0;
    Yom.Pump   = 0;
    Yom.Comp   = 0;
    Yom.Cooler = 0;
    Yom.WT     = 0;

% Total O&M
    Yom_total_PEM = Yom.Stack + Yom.Rect + Yom.Pump + Yom.Comp + Yom.Cooler + ...
                Yom.WT + Yom.Water;

%--------------------------------------------------------------------------
% Impacto total por componente (Y) ----------------------------------------
    Y.Stack  = Yco.Stack  + Yde.Stack  + Yom.Stack;
    Y.Rect   = Yco.Rect   + Yde.Rect   + Yom.Rect;
    Y.Pump   = Yco.Pump   + Yde.Pump   + Yom.Pump;
    Y.Comp   = Yco.Comp   + Yde.Comp   + Yom.Comp;
    Y.Cooler = Yco.Cooler + Yde.Cooler + Yom.Cooler;
    Y.WT     = Yco.WT     + Yde.WT     + Yom.WT;
    Y.Water  = Yco.Water  + Yde.Water  + Yom.Water;

%--------------------------------------------------------------------------
% Resultados ---------------------------------------------------------------
    Y_sum  = Y.Stack + Y.Rect + Y.Pump + Y.Comp + Y.Cooler + Y.WT + Y.Water;
    Ytotal_PEM = Y_sum / E_util;  % kgCO2-eq/kWh

% Vector de componentes (orden estilo tu Y_componentes)
    Y_componentes_PEM = [Y.Stack,Y.Rect,Y.Pump,Y.Comp,Y.Cooler,Y.WT,Y.Water];

end

