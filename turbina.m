function [h_out, T_out, s_out] = turbina(T,P,nt,P_out,fluido)
    %   Función de análisis termodinámico de la turbina
    %   Nos permite obtener las condiciones de salida de la turbina en
    %   función de las condiciones de entrada
    
    %---------------------------------------------------------
    %Condiciones de Entrada
    
    [h,s] = refpropm('HS','T',T,'P',P,fluido);
    %---------------------------------------------------------
   
    %---------------------------------------------------------
    %Condiciones de Salida PAPAYA
    s_outs = s; %Suposición Isentropica
    
    HS = refpropm('H','P',P_out,'S',s_outs,fluido);  %Entalpía Isentropica
    
    h_out = h - (nt*(h-HS));    %Entalpía de salida, usando eficiencia. 
    
    [T_out,s_out] = refpropm('TS','P',P_out,'H',h_out,fluido);
    %---------------------------------------------------------
end

