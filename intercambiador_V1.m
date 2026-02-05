 function [Tmin, T_hN, T_c0, sumUA,Qp_real] = intercambiador_V1(Th,Ph,Tc,Pc,fluido,eff,N,mh,mc)


    %   Función de análisis termodinámico del intercambiador de calor
    %   Nos permite obtener las condiciones de salida del intercambiador en
    %   función de las condiciones de entrada V1.1
    
    %---------------------------------------------------------
    % Condiciones de Operación
    % Las propiedades en las entradas de las corrientes caliente y fría se
    % determinan en función de la temperatura y la presión en estos puntos,
    % mientras que para las propiedades en las salidas de estos flujos, se
    % supone que lo máximo que se puede calentar el flujo frío es hasta la
    % tempreratura del flujo caliente, mientras que lo máximo que se puede
    % enfríar el flujo caliente es hasta la temperatura del flujo frío. 
    
    
    % Para el flujo caliente
    h_hi = refpropm('H','T',Th,'P',Ph,fluido);      %Entalpia Entrada corriente caliente
    h_ho = refpropm('H','T',Tc,'P',Ph,fluido);      %Entalpia Salida corriente caliente
        
    % Para el flujo frío
    h_ci = refpropm('H','T',Tc,'P',Pc,fluido);      %Entalpia Entrada corriente fría
    h_co = refpropm('H','T',Th,'P',Pc,fluido);      %Entalpia Salida corriente fría
    
    % Calor transferido
    Qp_max = min((mh*(h_hi - h_ho)),(mc*(h_co - h_ci)));
    Qp_real = Qp_max*eff; 
    
    %Entalpias en función del calor real
    h_ho = h_hi - (Qp_real/mh); 
    h_co = h_ci + (Qp_real/mc);
    
    q_nodes = Qp_real/N;
    T_co = refpropm('T','P',Pc,'H',h_co,fluido);
    nodes = N+1;
    
    % Vectores de propiedades
    h_h = zeros(1,nodes);        %Entalpía de los nodos calientes
    h_c = zeros(1,nodes);        %Entalpía de los nodos fríos
    T_h = zeros(1,nodes);        %Temperatura de los nodos calientes
    T_c = zeros(1,nodes);        %Temperatura de los nodos fríos
    % cmin = zeros(1,N);
    UA = zeros(1,N);
    
    
    h_h(1) = h_hi; 
    h_h(N+1) = h_ho;
    T_h(1) = Th; 
    T_c(1) = T_co;
    h_c(N+1) = h_ci;
    h_c(1) = h_co;
    T_c(N+1) = Tc;
    
    
    for i=1:N+1 
        
        h_h(i+1) = h_h(i)-(q_nodes/mh);
        h_c(i+1) = h_c(i)-(q_nodes/mc);
        T_h(i+1) = refpropm('T','P',Ph,'H',h_h(i+1),fluido);    
        T_c(i+1) = refpropm('T','P',Pc,'H',h_c(i+1),fluido);
        
        dTh = T_h(i) - T_h(i+1);
        dTc = T_c(i) - T_c(i+1);
        dh_h = h_h(i) - h_h(i+1);
        dh_c = h_c(i) - h_c(i+1);
        
        ch = abs(mh*dh_h/dTh);
        cc = abs(mc*dh_c/dTc);
        cmin = min(ch,cc);
        cmax = max(ch,cc);
        cr = cmin/cmax;
        qmax_i = cmin*(T_h(i) - T_c(i+1));
        effi = q_nodes/qmax_i;
        
        
        if abs(cr-1) < 1e-5
            NTU_i = effi/(1-effi);
            % NTU_i=0;

        else 
            NTU_i = (1/(cr-1))*log((effi-1) / ((effi*cr) - 1));
            % NTU_i =0;
            
        end
        
        UA(i) = NTU_i*cmin;    
        
    end
    
    Tmin = min(T_h - T_c);
    % dTlm = zeros(N);
    % dTL = Th-Tc;
    
    T_hN = T_h(N);
    T_c0 = T_c(1);
    sumUA = sum(UA);
    
 end

 