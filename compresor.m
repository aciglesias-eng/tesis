function [h_out,T_out,s_out] = compresor(T,P,nc,Phigh,fluido)
    %   Función de análisis termodinámico del compresor 
    %   Nos permite obtener las condiciones de salida del compresor en
    %   función de las condiiones de entrada
    
    %---------------------------------------------------------
    %Condiciones de Entrada
    
    [h,s] = refpropm('HS','T',T,'P',P,fluido);
    %---------------------------------------------------------
    
    
    %---------------------------------------------------------
    %Condiciones de Salida
    
    s_out = s; %Suposición Isentropica
    
    HS = refpropm('H','P',Phigh,'S',s_out,fluido);  %Entalpía Isentrópica
    
    h_out = h+((HS-h)/(nc));    %Entalpía de salida, usando eficiencia. 
    
    
    [T_out,s_out] = refpropm('TS','P',Phigh,'H',h_out,fluido);
    %---------------------------------------------------------
end

