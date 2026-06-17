clear%;clc
Ts=[72 74 76 78 80 82 84 112 114 116 118 120 122 124 126 128 130 132 134]+273.15;
P_CO2 = 9124;
m_CO2 = 180;
cycle_type = 'RORC';
ehr = 0.80;

% Pre-asignacion: los casos invalidos quedan como NaN en vez de 0/0
varout1 = nan(1,length(Ts));
varout2 = nan(1,length(Ts));

for i=1:length(Ts)
    display(i);
    try
        [T6AO3, A_ORC, W_t3, W_b2, Xd_orc, TPhsmx_orc, Q_cond_vec, Q_evap_vec, Q_regen] = ORCV2(Ts(i), P_CO2, m_CO2, cycle_type, ehr);
        if (T6AO3 > 0) && (Q_evap_vec(4) > 0)
            varout1(i)   = (W_t3+W_b2)/Q_evap_vec(4)*100; % MW
            varout2(i)   = (W_t3+W_b2)/1000; %MW
        end
    catch ME
        warning('Caso %d (Ts=%.2f C) fallo: %s', i, Ts(i)-273.15, ME.message);
    end
end
resultados = [ ...
    (Ts-273.15)', varout1', varout2'];
xlswrite('Salidas2.xlsx',resultados ,'T','B3')