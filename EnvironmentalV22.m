function [Ytotal, Yco_total, Yde_total, Yom_total, Y_componentes, Ytotal_PEM] = EnvironmentalV22( ...
    W_t1, W_t2, W_t3, W_c1, W_c2, W_b2, A_HTR, A_LTR, A_cooler, A_evap, A_cond, A_HE)

format long

% ------------------ PARAMETROS GENERALES ------------------
years          = 20;      % [años]
disponibilidad = 0.7;     % [-]
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
WNETO  = (W_t1 + W_t2 + W_t3) + (W_c1 + W_c2 + W_b2); % [kW]
E_util = WNETO * horas * years * disponibilidad;      % [kWh_e]

% ------------------ MASA EQUIVALENTE EQUIPOS ------------------
a_turb = 35.03;  % [kg/kW]
a_pump = 15.71;  % [kg/kW]

M_Evap   = A_evap   * delta * rho;
M_Cond   = A_cond   * delta * rho;
M_HTR    = A_HTR    * delta * rho;
M_LTR    = A_LTR    * delta * rho;
M_Cooler = A_cooler * delta * rho;
M_HE     = A_HE     * delta * rho;

M_T1 = a_turb * W_t1;
M_T2 = a_turb * W_t2;
M_T3 = a_turb * W_t3;

M_C1 = a_pump * W_c1;
M_C2 = a_pump * W_c2;
M_B2 = a_pump * W_b2;

% ------------------ FLUIDOS ------------------
wf_c  = 5.57 * W_t3;
wf_op = wf_c * 0.10;
wf_de = (wf_c - wf_op) * 0.03;

co2_c  = 5.57 * (W_t1 + W_t2);
co2_op = co2_c * 0.10;
co2_de = (co2_c - co2_op) * 0.03;

% ------------------ DECOMISIONAMIENTO (bloque potencia) ------------------
Yde_block = peed * mu_eq_mass * (M_Evap + M_Cond + M_HTR + M_LTR + M_Cooler + M_HE) +...
            peed * mu_eq_mass * (M_T1 + M_T2 + M_T3 + M_C1 + M_C2 + M_B2) +...
            (mu_tolueno_mass * wf_de) + (mu_CO2_mass * co2_de);

% ------------------ COMISIONAMIENTO / CAPEX (bloque potencia) ------------------
Yco_block = mu_eq_mass * (M_Evap + M_Cond + M_HTR + M_LTR + M_Cooler + M_HE) +...
            mu_eq_mass * (M_T1 + M_T2 + M_T3 + M_C1 + M_C2 + M_B2) +...
            (mu_tolueno_mass * wf_c) + (mu_CO2_mass * co2_c);

% ------------------ O&M (bloque potencia) ------------------
Yom_block = (mu_tolueno_mass * wf_op) + (mu_CO2_mass * co2_op);

% ------------------ CSP (impacto por energía) ------------------
Y_CSP_abs  = mu_CSP * E_util;        % [kgCO2-eq]
Yco_CSP    = Y_CSP_abs;
Yde_CSP    = peed * Y_CSP_abs;  

% ------------------ FACTOR ELECTRICO TOTAL (bloque de potencia) ------------------
GWP_elec = (Yco_block + Yom_block + Yde_block) / E_util; % [kgCO2-eq/kWh_e]

% ------------------ ACOPLE PEM ------------------
Ytotal_PEM = PEM_LCA_from_paper(GWP_elec, WNETO, years, disponibilidad, horas);

% ------------------ TOTALES SISTEMA COMPLETO ------------------
Yco_total = Yco_block + Yco_CSP + Ytotal_PEM.Yco_abs;
Yde_total = Yde_block + Yde_CSP + Ytotal_PEM.Yde_abs;
Yom_total = Yom_block + Yom_CSP + Ytotal_PEM.Yom_abs;

Y_sum_abs = Yco_total + Yde_total + Yom_total;
Ytotal    = Y_sum_abs / E_util; % [kgCO2-eq/kWh_e]

% Componentes: [bloque potencia, CSP, PEM]
Y_componentes = [ (Yco_block + Yom_block + Yde_block), Y_CSP_abs, Ytotal_PEM.Ytotal_abs ];

end

% ------------------ PEM: materiales + electricidad (como el paper) ------------------
function PEM = PEM_LCA_from_paper(GWP_elec, WNETO, years, disponibilidad, horas)

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
    EF_titanium   = 13.5;
    EF_iridium    = 12000;
    EF_platinum   = 35000;
    EF_act_carbon = 2.0;
    EF_nafion     = 6.0;
    EF_plastics   = 3.0;
    EF_cement     = 0.12;

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
    GWP_total_perkg = GWP_mfg_perkg + GWP_op_perkg;

    E_to_PEM_life = WNETO * horas * years * disponibilidad;
    H2_life       = E_to_PEM_life / e_PEM;


    PEM.Yco_abs    = GWP_mfg_perkg * H2_life;
    PEM.Yom_abs    = GWP_op_perkg  * H2_life;
    PEM.Yde_abs    = 0;
    PEM.Ytotal_abs = PEM.Yco_abs + PEM.Yom_abs + PEM.Yde_abs;

    PEM.GWP_mfg_perkg   = GWP_mfg_perkg;
    PEM.GWP_op_perkg    = GWP_op_perkg;
    PEM.GWP_total_perkg = GWP_total_perkg;
end