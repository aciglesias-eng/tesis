function [YCOPROD, GWP_elec, Yco_total, Yde_total, Yom_total, Y_componentes, Ytotal_PEM,Yco_total_PEM, Yde_total_PEM, Yom_total_PEM] = EnvironmentalV22( ...
    W_t1, W_t2, W_t3, W_c1, W_c2, W_b2, A_HTR, A_LTR, A_cooler, A_evap, A_cond, A_HE)

format long

% ------------------ PARAMETROS GENERALES ------------------
years          = 20;      % [años]
disponibilidad = 0.85;     % [-]
horas          = 24*365;  % [h/año]
peed           = 0.05;    % 

% ------------------ GEOMETRIA / MATERIAL HEX ------------------
rho   = 8900;    % [kg/m3]
delta = 0.002;   % [m]

% Factores
mu_eq_mass = 2.20;      % [kgCO2-eq/kg]
mu_CSP     = 9.456e-3;  % [kgCO2-eq/kWh_e]
mu_tolueno_mass = 32.3; % [kgCO2-eq/kg]
mu_CO2_mass     = 1.0;  % [kgCO2-eq/kg]

% ------------------ POTENCIA NETA Y ENERGIA------------------
% Convención física: turbinas producen (+), compresores/bombas consumen (-)
WNETO  = (W_t1 + W_t2 + W_t3) - (W_c1 + W_c2 + W_b2); % [kW]
E_util = WNETO * horas * years * disponibilidad;      % [kWh_e]

% ------------------ MASA EQUIVALENTE EQUIPOS ------------------
a_turb = 35.03;  % [kg/kW]
a_pump = 15.71;  % [kg/kW]

M_Evap   = abs(A_evap   * delta * rho);
M_Cond   = abs(A_cond   * delta * rho);
M_HTR    = abs(A_HTR    * delta * rho);
M_LTR    = abs(A_LTR    * delta * rho);
M_Cooler = abs(A_cooler * delta * rho);
M_HE     = abs(A_HE     * delta * rho);

M_T1 = a_turb * W_t1;
M_T2 = a_turb * W_t2;
M_T3 = a_turb * W_t3;

M_C1 = a_pump * abs(W_c1);
M_C2 = a_pump * abs(W_c2);
M_B2 = a_pump * abs(W_b2);

% ------------------ FLUIDOS ------------------
wf_c  = 5.57 * W_t3;
wf_op = wf_c * 0.10;
wf_de = (wf_c - wf_op) * 0.03;

co2_c  = 5.57 * (W_t1 + W_t2);
co2_op = co2_c * 0.10;
co2_de = (co2_c - co2_op) * 0.03;
%--------------------------------------------------------------------------
% Decomisionamiento
    Yde.Evap   = peed*mu_eq_mass*M_Evap;
    Yde.Cond   = peed*mu_eq_mass*M_Cond;
    Yde.HTR    = peed*mu_eq_mass*M_HTR;
    Yde.LTR    = peed*mu_eq_mass*M_LTR;
    Yde.Cooler = peed*mu_eq_mass*M_Cooler;
    Yde.HE     = peed*mu_eq_mass*M_HE;

    Yde.T1     = peed*mu_eq_mass*M_T1;
    Yde.T2     = peed*mu_eq_mass*M_T2;
    Yde.T3     = peed*mu_eq_mass*M_T3;

    Yde.C1     = peed*mu_eq_mass*M_C1;
    Yde.C2     = peed*mu_eq_mass*M_C2;
    Yde.B2     = peed*mu_eq_mass*M_B2;

    Yde.wf     = mu_tolueno_mass*wf_de;
    Yde.co2    = mu_CO2_mass*co2_de;

    Yde.CSP    = peed*(E_util*mu_CSP);
   
    Yde_block  = Yde.Evap+Yde.HE+Yde.Cond+Yde.HTR+Yde.LTR+Yde.Cooler+Yde.T1...
                +Yde.T2+Yde.T3+Yde.C1+Yde.C2+Yde.B2+Yde.wf+Yde.co2+Yde.CSP; %
%--------------------------------------------------------------------------
% Comisionamiento
    Yco.Evap   = mu_eq_mass*M_Evap;
    Yco.Cond   = mu_eq_mass*M_Cond;
    Yco.HTR    = mu_eq_mass*M_HTR;
    Yco.LTR    = mu_eq_mass*M_LTR;
    Yco.Cooler = mu_eq_mass*M_Cooler;
    Yco.HE     = mu_eq_mass*M_HE;

    Yco.T1     = mu_eq_mass*M_T1;
    Yco.T2     = mu_eq_mass*M_T2;
    Yco.T3     = mu_eq_mass*M_T3;

    Yco.C1     = mu_eq_mass*M_C1;
    Yco.C2     = mu_eq_mass*M_C2;
    Yco.B2     = mu_eq_mass*M_B2;

    Yco.wf     = mu_tolueno_mass*wf_c;
    Yco.co2    = mu_CO2_mass*co2_c+2;

    Yco.CSP    = mu_CSP*E_util;

    Yco_block  = Yco.Evap+Yco.HE+Yco.Cond+Yco.HTR+Yco.LTR+Yco.Cooler+Yco.T1...
                +Yco.T2+Yco.T3+Yco.C1+Yco.C2+Yco.B2+Yco.wf+Yco.co2+Yco.CSP; %
%--------------------------------------------------------------------------
% Operación y mantenimiento
    Yom.Evap   = 0;
    Yom.Cond   = 0;
    Yom.HTR    = 0;
    Yom.LTR    = 0;
    Yom.Cooler = 0;
    Yom.HE     = 0;

    Yom.T1     = 0;
    Yom.T2     = 0;
    Yom.T3     = 0;

    Yom.C1     = 0;
    Yom.C2     = 0;
    Yom.B2     = 0;

    Yom.wf     = 0.10 * Yco.wf;
    Yom.co2    = 0.10 * Yco.co2;
    
    Yom.CSP    = 0;
    
    Yom_block = Yom.wf+Yom.co2;
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
    
% ------------------ FACTOR ELECTRICO TOTAL (bloque de potencia) ------------------
GWP_elec = (Yco_block + Yom_block + Yde_block)/E_util; % [kgCO2-eq/kWh_e]

% ------------------ ACOPLE PEM ------------------
[Ytotal_PEM,Yco_total_PEM,Yde_total_PEM,Yom_total_PEM]= PEM_LCA(GWP_elec, WNETO, years, disponibilidad, horas);
Y_PEME = Yom_total_PEM+Yde_total_PEM+Yco_total_PEM;
% ------------------ TOTALES SISTEMA COMPLETO ------------------
Yco_total = Yco_block;
Yde_total = Yde_block;
Yom_total = Yom_block;



Y_componentes = [Y.Evap,Y.Cond,Y.HTR,Y.LTR,Y.Cooler,Y.HE,Y.T1,Y.T2,Y.T3, ...
                 Y.C1,Y.C2,Y.B2,Y.wf,Y.co2,Y.CSP,Y_PEME];
YCOPROD = Y.Evap+Y.Cond+Y.HTR+Y.LTR+Y.Cooler+Y.HE+Y.T1+Y.T2+Y.T3+ ...
                 Y.C1+Y.C2+Y.B2+Y.wf+Y.co2+Y.CSP+Y_PEME;
             
end

function [Ytotal_PEM,Yco_total_PEM,Yde_total_PEM,Yom_total_PEM]=PEM_LCA(GWP_elec, WNETO, years, disponibilidad, horas)

    % O&M PEM
    e_PEM = 54.81; % [kWh/kgH2]

    % Manufactura PEM kg material/kgH2
    m_low_alloy   = 5.30E-3;
    m_high_alloy  = 2.10E-3;
    m_stainless   = 1.11E-4;
    m_aluminum    = 1.40E-4;
    m_copper      = 3.36E-4;
    m_titanium    = 5.83E-4;
    m_iridium     = 8.30E-7;
    m_platinum    = 8.00E-8;
    m_act_carbon  = 9.94E-6;
    m_nafion      = 1.77E-5;
    m_plastics    = 1.11E-4;
    m_cement      = 6.19E-3;


    EF_steel_low  = 1.7; %kgCO2/kg
    EF_steel_high = 1.7;
    EF_stainless  = 6.1;
    EF_aluminum   = 8.32;
    EF_copper     = 2.7;
    EF_titanium   = 20.6;
    EF_iridium    = 12000;
    EF_platinum   = 35000;
    EF_act_carbon = 2.0;
    EF_nafion     = 6.0;
    EF_plastics   = 3.31;
    EF_cement     = 0.84;

    GWP_mfg_perkg = ...
        m_low_alloy  * EF_steel_low  + ...
        m_high_alloy * EF_steel_high + ...
        m_stainless  * EF_stainless  + ...
        m_aluminum   * EF_aluminum   + ...
        m_copper     * EF_copper     + ...
        m_titanium   * EF_titanium   + ...
        m_iridium    * EF_iridium    + ...
        m_platinum   * EF_platinum   + ...
        m_act_carbon * EF_act_carbon + ...
        m_nafion     * EF_nafion     + ...
        m_plastics   * EF_plastics   + ...
        m_cement     * EF_cement;

    GWP_op_perkg    = e_PEM * GWP_elec;

    E_util = WNETO * horas * years * disponibilidad;
    H2_life       = E_util/e_PEM;

    Yco_total_PEM    = 0;
    Yom_total_PEM    = GWP_op_perkg*H2_life + GWP_mfg_perkg*H2_life;
    Yde_total_PEM    = 0;
    Ytotal_PEM       = GWP_op_perkg + GWP_mfg_perkg;
end