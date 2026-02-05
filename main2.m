clear%;clc
T=[500 550 600 650 700 750 800]+273.15;
nc=0.89;
nt=0.93;
rc=2.74;
pmaire=100; % EN PORCENTAJE

varout1  = nan(1,length(T)); varout2  = nan(1,length(T)); varout3  = nan(1,length(T)); varout4  = nan(1,length(T));
varout5  = nan(1,length(T)); varout6  = nan(1,length(T)); varout7  = nan(1,length(T)); varout8  = nan(1,length(T));
varout9  = nan(1,length(T)); varout10 = nan(1,length(T)); varout11 = nan(1,length(T)); varout12 = nan(1,length(T));
varout13 = nan(1,length(T)); varout14 = nan(1,length(T)); varout15 = nan(1,length(T)); varout16 = nan(1,length(T));
varout17 = nan(1,length(T)); varout18 = nan(1,length(T)); varout19 = nan(1,length(T)); varout20 = nan(1,length(T));
varout21 = nan(1,length(T)); varout22 = nan(1,length(T)); varout23 = nan(1,length(T)); varout24 = nan(1,length(T));
varout25 = nan(1,length(T)); varout26 = nan(1,length(T)); varout27 = nan(1,length(T)); varout28 = nan(1,length(T));
varout29 = nan(1,length(T)); varout30 = nan(1,length(T)); varout31 = nan(1,length(T)); varout32 = nan(1,length(T));
varout33 = nan(1,length(T)); varout34 = nan(1,length(T)); varout35 = nan(1,length(T)); varout36 = nan(1,length(T));


for i=1:length(T)   
[WNETO, W_Rankine, n_ciclo, n_ex, n_ex_2, n_th, n_en_PEM,n_ex_PEM, T20, m_H2, m_H2O_in, XD, In, Ytotal, Ytotal_PEM] = SBCRV2(nt, nc, T(i), rc, pmaire);
varout1(i)   = WNETO/1000000; % MW
varout2(i)   = W_Rankine/1000000; %MW
varout3(i)   = n_ciclo;    % Eficiencia energética ciclo
varout4(i)   = n_ex;       % Eficiencia exergética 
varout5(i)   = n_th;       % Eficiencia termica
varout6(i)   = n_en_PEM;   % Eficiencia termica PEM
varout7(i)   = n_ex_PEM;   % Eficiencia exergetica PEM 
varout8(i)   = T20-273.15; % Temperatura entrada PEM
varout9(i)   = m_H2*3600;       % Flujo masico hidrogeno
varout10(i)  = m_H2O_in*3600;   % Flujo masico de entrada agua
varout11(i)  = XD(1,1);    % Comprensor 1
varout12(i)  = XD(1,2);    % Comprensor 2
varout13(i)  = XD(1,3);    % Turbina 1
varout14(i)  = XD(1,4);    % Turbina 2
varout15(i)  = XD(1,5);    % Turbina 3
varout16(i)  = XD(1,6);    % HTR
varout17(i)  = XD(1,7);    % LTR
varout18(i)  = XD(1,8);    % Loss solar receiver
varout19(i)  = XD(1,9);    % Loss field
varout20(i)  = XD(1,10);   % Destruida CSP
varout21(i)  = XD(1,11);   % Destruida cooler
varout22(i)  = XD(1,12);   % Perdida en el aire
varout23(i)  = XD(1,13);   % Destruida evaporador
varout24(i)  = XD(1,14);   % Destruida condensador
varout25(i)  = XD(1,15);   % Destruida bomba 2
varout26(i)  = XD(1,16);   % Perdida en el agua
varout27(i)  = XD(1,17);   % Destruida HE
varout28(i)  = XD(1,18);   % Exergia aportada por el sistema solar
varout29(i)  = XD(1,19);   % Exergia destruida PEM
varout30(i)  = XD(1,20);   % Exergia destruida total

varout31(i) = In(1,1);     % Exergy waste ratio EWR
varout32(i) = In(1,2);     % Environmental effect factor (EEF)
varout33(i) = In(1,3);     % Exergetic sustainability index (ESI)

varout34(i) = Ytotal;     % Impacto ambiental total;
varout35(i) = Ytotal_PEM; % Impacto ambiental PEM respecto a Rankine
varout36(i) = n_ex_2;
end


resultados = [ ...
    (T-273.15)', varout1', varout2', varout3', varout4', varout5', varout6', ...
    varout7', varout8', varout9', varout10'
];
resultados2 = [(T-273.15)', varout11', varout12', varout13', ...
    varout14', varout15', varout16', varout17', varout18', varout19', varout20', ...
    varout21', varout22', varout23', varout24', varout25', varout26', varout27', ...
    varout28', varout29', varout30'
];
resultados3 = [(T-273.15)', varout31', varout32', varout33', varout34', varout35', varout36'
    
];



xlswrite('Salidas.xlsx',resultados,'T','B3')
xlswrite('Salidas.xlsx',resultados2,'T','B16')
xlswrite('Salidas.xlsx',resultados3,'T','B29')








