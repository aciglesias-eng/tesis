%{
%% 
clear%;clc
rc = [2.50 2.70 2.90 3.10 3.30 3.50];
nc = 0.89;
nt = 0.929;
T = 900 + 273.15;
pmaire=100; % EN PORCENTAJE


for i=1:length(rc)   
[LCOH, c_el_bus,LCOEn, WNETO, W_Rankine, n_ciclo, n_ex, n_ex_2, n_th, n_en_PEM,n_ex_PEM, T20, m_H2, m_H2O_in, XD, In, Ytotal, Ytotal_PEM] = SBCRV2(nt, nc, T, rc(i), pmaire);
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
varout35(i) = Ytotal_PEM; % GWP PEM respecto a total kgCO2eq/kgH2
varout36(i) = n_ex_2;

varout37(i) = LCOEn;
varout38(i) = LCOH;
varout39(i) = c_el_bus;
end


resultados = [ ...
    (rc)', varout1', varout2', varout3', varout4', varout5', varout6', ...
    varout7', varout8', varout9', varout10'
];
resultados2 = [(rc)', varout11', varout12', varout13', ...
    varout14', varout15', varout16', varout17', varout18', varout19', varout20', ...
    varout21', varout22', varout23', varout24', varout25', varout26', varout27', ...
    varout28', varout29', varout30'
];
resultados3 = [(rc)', varout31', varout32', varout33', varout34', varout35',...
    varout36', varout37', varout38', varout39'
    ];



xlswrite('Salidas.xlsx',resultados ,'rc','B3')
xlswrite('Salidas.xlsx',resultados2 ,'rc','B16')
xlswrite('Salidas.xlsx',resultados3 ,'rc','B29')

%% 
clear%;clc
nt =[0.60 0.65 0.70 0.75 0.80 0.85 0.90];
nc = 0.89;
T = 900 + 273.15;
rc = 2.74;
pmaire=100; % EN PORCENTAJE


for i=1:length(nt)   
[LCOH, c_el_bus,LCOEn, WNETO, W_Rankine, n_ciclo, n_ex, n_ex_2, n_th, n_en_PEM,n_ex_PEM, T20, m_H2, m_H2O_in, XD, In, Ytotal, Ytotal_PEM] = SBCRV2(nt(i), nc, T, rc, pmaire);
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
varout35(i) = Ytotal_PEM; % GWP PEM respecto a total kgCO2eq/kgH2
varout36(i) = n_ex_2;

varout37(i) = LCOEn;
varout38(i) = LCOH;
varout39(i) = c_el_bus;
end


resultados = [ ...
    (nt)', varout1', varout2', varout3', varout4', varout5', varout6', ...
    varout7', varout8', varout9', varout10'
];
resultados2 = [(nt)', varout11', varout12', varout13', ...
    varout14', varout15', varout16', varout17', varout18', varout19', varout20', ...
    varout21', varout22', varout23', varout24', varout25', varout26', varout27', ...
    varout28', varout29', varout30'
];
resultados3 = [(nt)', varout31', varout32', varout33', varout34', varout35',...
    varout36', varout37', varout38', varout39'
];



xlswrite('Salidas.xlsx',resultados ,'nt','B3')
xlswrite('Salidas.xlsx',resultados2 ,'nt','B16')
xlswrite('Salidas.xlsx',resultados3 ,'nt','B29')


%% 
clear%;clc
nc = [0.60 0.65 0.70 0.75 0.80 0.85 0.90];
nt = 0.929;
T = 900 + 273.15;
rc = 2.74;
pmaire=1; % EN PORCENTAJE




for i=1:length(nc)   
[LCOH, c_el_bus,LCOEn, WNETO, W_Rankine, n_ciclo, n_ex, n_ex_2, n_th, n_en_PEM,n_ex_PEM, T20, m_H2, m_H2O_in, XD, In, Ytotal, Ytotal_PEM] = SBCRV2(nt, nc(i), T, rc, pmaire);
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
varout35(i) = Ytotal_PEM; % GWP PEM respecto a total kgCO2eq/kgH2
varout36(i) = n_ex_2;

varout37(i) = LCOEn;
varout38(i) = LCOH;
varout39(i) = c_el_bus;
end


resultados = [ ...
    (nc)', varout1', varout2', varout3', varout4', varout5', varout6', ...
    varout7', varout8', varout9', varout10'
];
resultados2 = [(nc)', varout11', varout12', varout13', ...
    varout14', varout15', varout16', varout17', varout18', varout19', varout20', ...
    varout21', varout22', varout23', varout24', varout25', varout26', varout27', ...
    varout28', varout29', varout30'
];
resultados3 = [(nc)', varout31', varout32', varout33', varout34', varout35',...
    varout36', varout37', varout38', varout39'
];



xlswrite('Salidas.xlsx',resultados ,'nc','B3')
xlswrite('Salidas.xlsx',resultados2 ,'nc','B16')    
xlswrite('Salidas.xlsx',resultados3 ,'nc','B29')
%}

%%
clear%;clc
T=[500 550 600 650 700 750 800 850 900 800]+273.15;
nc=0.89;
nt=0.929;
rc=2.74;
pmaire=100; % EN PORCENTAJE


for i=1:length(T)   
    display(i);
[LCOH, c_el_bus,LCOEn, WNETO, W_Rankine, n_ciclo, n_ex, n_ex_2, n_th, n_en_PEM,n_ex_PEM, T20, m_H2, m_H2O_in, XD, In, Ytotal, Ytotal_PEM] = SBCRV2(nt, nc, T(i), rc, pmaire);
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
varout30(i)  = XD(1,20);   % Exergia destruida regenerador

varout31(i) = In(1,1);     % Exergy waste ratio EWR
varout32(i) = In(1,2);     % Environmental effect factor (EEF)
varout33(i) = In(1,3);     % Exergetic sustainability index (ESI)

varout34(i) = Ytotal;     % Impacto ambiental total;
varout35(i) = Ytotal_PEM; % GWP PEM respecto a total kgCO2eq/kgH2
varout36(i) = n_ex_2;
varout37(i) = LCOEn;
varout38(i) = LCOH;
varout39(i) = c_el_bus;
varout40(i) = XD(1,21);   % Exergia destruida total
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
resultados3 = [(T-273.15)', varout31', varout32', varout33', varout34', varout35',...
    varout36', varout37', varout38', varout39'
    
];



xlswrite('Salidas.xlsx',resultados ,'T','B3')
xlswrite('Salidas.xlsx',resultados2 ,'T','B16')
xlswrite('Salidas.xlsx',resultados3 ,'T','B29')

