clear; clc; close all

% Fluido
fluid = 'butane';
% Ejemplos:
% fluid = 'WATER';
% fluid = 'CO2';
% fluid = 'AMMONIA';
% fluid = 'PROPANE';

% Nombre del archivo Excel de salida
excelFile = ['campana_TS_REFPROP_', fluid, '.xlsx'];

% Número de puntos
N = 600;


try
    T_test = refpropm('T','P',101.325,'Q',0,'WATER');
    fprintf('REFPROP funciona. T_sat agua a 101.325 kPa = %.4f K\n', T_test);
catch ME
    error(['No se pudo llamar a REFPROP. Revisa la ruta, refpropm.m, ', ...
           'REFPRP64.dll y REFPROP.h. Error: ', ME.message])
end

%% ==========================
%  TEMPERATURA CRÍTICA Y TRIPLE
%  ==========================

% Temperatura crítica
Tcrit = refpropm('T','C',0,' ',0,fluid);

% Temperatura triple o límite inferior
try
    Tmin = refpropm('T','R',0,' ',0,fluid);
catch
    % Si falla el punto triple para algún fluido, coloca un valor manual.
    % Para R134a, por ejemplo, puedes usar aproximadamente 169.85 K.
    warning('No se pudo obtener Ttriple. Usando Tmin manual.')
    Tmin = 170;
end

% Evitar exactamente los extremos para no tener problemas numéricos
T = linspace(Tmin + 0.01, Tcrit - 0.01, N)';

%% ==========================
%  INICIALIZAR VECTORES
%  ==========================

P_sat_kPa = nan(N,1);

s_liq_JkgK = nan(N,1);
s_vap_JkgK = nan(N,1);

h_liq_Jkg = nan(N,1);
h_vap_Jkg = nan(N,1);

rho_liq_kgm3 = nan(N,1);
rho_vap_kgm3 = nan(N,1);

%% ==========================
%  CALCULAR CURVA DE SATURACIÓN
%  ==========================

for i = 1:N

    Ti = T(i);

    try
        % Presión de saturación
        P_sat_kPa(i) = refpropm('P','T',Ti,'Q',0,fluid);

        % Líquido saturado, Q = 0
        s_liq_JkgK(i)   = refpropm('S','T',Ti,'Q',0,fluid);
        h_liq_Jkg(i)    = refpropm('H','T',Ti,'Q',0,fluid);
        rho_liq_kgm3(i) = refpropm('D','T',Ti,'Q',0,fluid);

        % Vapor saturado, Q = 1
        s_vap_JkgK(i)   = refpropm('S','T',Ti,'Q',1,fluid);
        h_vap_Jkg(i)    = refpropm('H','T',Ti,'Q',1,fluid);
        rho_vap_kgm3(i) = refpropm('D','T',Ti,'Q',1,fluid);

    catch
        % Si REFPROP no converge en algún punto, se deja como NaN
        P_sat_kPa(i) = NaN;
        s_liq_JkgK(i) = NaN;
        s_vap_JkgK(i) = NaN;
        h_liq_Jkg(i) = NaN;
        h_vap_Jkg(i) = NaN;
        rho_liq_kgm3(i) = NaN;
        rho_vap_kgm3(i) = NaN;
    end
end

%% ==========================
%  LIMPIAR DATOS INVÁLIDOS
%  ==========================

valid = ~isnan(P_sat_kPa) & ...
        ~isnan(s_liq_JkgK) & ...
        ~isnan(s_vap_JkgK);

T = T(valid);
P_sat_kPa = P_sat_kPa(valid);

s_liq_JkgK = s_liq_JkgK(valid);
s_vap_JkgK = s_vap_JkgK(valid);

h_liq_Jkg = h_liq_Jkg(valid);
h_vap_Jkg = h_vap_Jkg(valid);

rho_liq_kgm3 = rho_liq_kgm3(valid);
rho_vap_kgm3 = rho_vap_kgm3(valid);

%% ==========================
%  CONVERSIONES
%  ==========================

T_C = T - 273.15;

s_liq_kJkgK = s_liq_JkgK / 1000;
s_vap_kJkgK = s_vap_JkgK / 1000;

h_liq_kJkg = h_liq_Jkg / 1000;
h_vap_kJkg = h_vap_Jkg / 1000;

%% ==========================
%  CREAR TABLA PARA EXCEL
%  ==========================

tabla = table( ...
    T, ...
    T_C, ...
    P_sat_kPa, ...
    s_liq_kJkgK, ...
    s_vap_kJkgK, ...
    h_liq_kJkg, ...
    h_vap_kJkg, ...
    rho_liq_kgm3, ...
    rho_vap_kgm3, ...
    'VariableNames', { ...
        'T_K', ...
        'T_C', ...
        'P_sat_kPa', ...
        's_liq_kJ_kgK', ...
        's_vap_kJ_kgK', ...
        'h_liq_kJ_kg', ...
        'h_vap_kJ_kg', ...
        'rho_liq_kg_m3', ...
        'rho_vap_kg_m3' ...
    });


writetable(tabla, excelFile, 'Sheet', 'Campana_TS');

fprintf('Archivo Excel creado: %s\n', excelFile);

