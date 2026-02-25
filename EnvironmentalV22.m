function [Ytotal, Yco_total, Yde_total, Yom_total, Y_componentes, Ytotal_PEM] = EnvironmentalV2( ...
    W_t1, W_t2, W_t3, W_c1, W_c2, W_b2, A_HTR, A_LTR, A_cooler, A_evap, A_cond, A_HE)

format long

% ------------------ PARAMETROS GENERALES ------------------
years          = 20;      % [años]
disponibilidad = 0.7;     % [-]
horas          = 24*365;  % [h/año]
peed           = 0.05;    % [-] proxy EoL si decides usarlo

% ------------------ GEOMETRIA / MATERIAL HEX ------------------
rho   = 8900;    % [kg/m3] (cobre)
delta = 0.002;   % [m]

% Factor 
mu_eq_mass = 2.20;      % [kgCO2-eq/kg]
mu_CSP    = 9.456e-3;   % [kgCO2-eq/kWh]
% Fluidos (por masa)
mu_tolueno_mass = 32.3; % [kgCO2-eq/kg]
mu_CO2_mass     = 1.0;  % [kgCO2-eq/kg]

% ------------------ POTENCIA NETA Y ENERGIA VIDA ------------------
% Convención física: turbinas generan (+), compresores/bombas consumen (-)
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

% ------------------ FLUIDOS (tu formulación) ------------------
wf_c  = 5.57 * W_t3;
wf_op = wf_c * 0.10;
wf_de = (wf_c - wf_op) * 0.03;

co2_c  = 5.57 * (W_t1 + W_t2);
co2_op = co2_c * 0.10;
co2_de = (co2_c - co2_op) * 0.03;

% ------------------ DECOMISIONAMIENTO (separado) ------------------
Yde.Evap   = peed * mu_eq_mass * M_Evap;
Yde.Cond   = peed * mu_eq_mass * M_Cond;
Yde.HTR    = peed * mu_eq_mass * M_HTR;
Yde.LTR    = peed * mu_eq_mass * M_LTR;
Yde.Cooler = peed * mu_eq_mass * M_Cooler;
Yde.HE     = peed * mu_eq_mass * M_HE;

Yde.T1     = peed * mu_eq_mass * M_T1;
Yde.T2     = peed * mu_eq_mass * M_T2;
Yde.T3     = peed * mu_eq_mass * M_T3;

Yde.C1     = peed * mu_eq_mass * M_C1;
Yde.C2     = peed * mu_eq_mass * M_C2;
Yde.B2     = peed * mu_eq_mass * M_B2;

Yde.wf     = mu_tolueno_mass * wf_de;
Yde.co2    = mu_CO2_mass     * co2_de;

% ------------------ COMISIONAMIENTO / CAPEX ------------------
Yco.Evap   = mu_eq_mass * M_Evap;
Yco.Cond   = mu_eq_mass * M_Cond;
Yco.HTR    = mu_eq_mass * M_HTR;
Yco.LTR    = mu_eq_mass * M_LTR;
Yco.Cooler = mu_eq_mass * M_Cooler;
Yco.HE     = mu_eq_mass * M_HE;

Yco.T1     = mu_eq_mass * M_T1;
Yco.T2     = mu_eq_mass * M_T2;
Yco.T3     = mu_eq_mass * M_T3;

Yco.C1     = mu_eq_mass * M_C1;
Yco.C2     = mu_eq_mass * M_C2;
Yco.B2     = mu_eq_mass * M_B2;

Yco.wf     = mu_tolueno_mass * wf_c;
Yco.co2    = mu_CO2_mass     * co2_c;

% ------------------ O&M ------------------
Yom.wf  = mu_tolueno_mass * wf_op;
Yom.co2 = mu_CO2_mass     * co2_op;

% ------------------ TOTALES DEL BLOQUE DE POTENCIA (sin PEM) ------------------
Yco_block = Yco.Evap + Yco.Cond + Yco.HTR + Yco.LTR + Yco.Cooler + Yco.HE + ...
            Yco.T1 + Yco.T2 + Yco.T3 + Yco.C1 + Yco.C2 + Yco.B2 + Yco.wf + Yco.co2;

Yde_block = Yde.Evap + Yde.Cond + Yde.HTR + Yde.LTR + Yde.Cooler + Yde.HE + ...
            Yde.T1 + Yde.T2 + Yde.T3 + Yde.C1 + Yde.C2 + Yde.B2 + Yde.wf + Yde.co2;

Yom_block = Yom.wf + Yom.co2;

% ---- ESTE ES EL FACTOR CORRECTO PARA ELECTRICIDAD QUE ALIMENTA LA PEM ----
GWP_elec = (Yco_block + Yom_block + Yde_block) / E_util; % [kgCO2-eq/kWh_e]

% ------------------ ACOPLE PEM (artículo: electricidad + materiales) ------------------
Ytotal_PEM = PEM_LCA_from_paper(GWP_elec, WNETO, years, disponibilidad, horas);

% CSP
    Yde_CSP = peed*(E_util*mu_CSP);
    Yco_CSP = mu_CSP*E_util+Yde.CSP;

% ------------------ TOTALES CON PEM Y CSP------------------
Yco_total = Yco_block + Ytotal_PEM.Yco_abs + Yco_CSP;
Yde_total = Yde_block + Ytotal_PEM.Yde_abs + Yde_CSP;
Yom_total = Yom_block + Ytotal_PEM.Yom_abs;

Y_sum_abs = Yco_total + Yde_total + Yom_total;
Ytotal    = Y_sum_abs / E_util; % [kgCO2-eq/kWh_e] del sistema completo (bloque + PEM)

% Vector componentes (ahora sí incluye PEM)
Y_componentes = [Yco_block+Yde_block+Yom_block, Ytotal_PEM.Ytotal_abs];

end

% ------------------ PEM: LCI del artículo (materiales + electricidad) ------------------
function PEM = PEM_LCA_from_paper(GWP_elec, WNETO, years, disponibilidad, horas)

    % O&M PEM del artículo (Tabla 10)
    e_PEM = 54.81; % [kWh/kgH2]  (paper)

    % Manufactura PEM del artículo (Tabla 7): kg material/kgH2
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

    % EF por masa (pon aquí tus EF reales; estos son placeholders)
    EF_steel_low  = 1.7;
    EF_steel_high = 1.7;
    EF_stainless  = 6.1;
    EF_aluminum   = 9.2;
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

    % Energía vida disponible al PEM (si el PEM consume toda la salida neta)
    E_to_PEM_life = WNETO * horas * years * disponibilidad;
    H2_life       = E_to_PEM_life / e_PEM;

    % En el paper: construcción PEM excluida y decom no cuantificado
    PEM.Yco_abs    = GWP_mfg_perkg * H2_life; % lo tratamos como CAPEX
    PEM.Yom_abs    = GWP_op_perkg  * H2_life; % operación
    PEM.Yde_abs    = 0;
    PEM.Ytotal_abs = PEM.Yco_abs + PEM.Yom_abs + PEM.Yde_abs;

    PEM.GWP_mfg_perkg   = GWP_mfg_perkg;
    PEM.GWP_op_perkg    = GWP_op_perkg;
    PEM.GWP_total_perkg = GWP_total_perkg;

end