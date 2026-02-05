Caso = 2;
% ENTRADAS A LA FUNCIÓN
    % Constantes físicas
        R = 8.314;     % J/mol-K
        F = 96485;     % C/mol
        M_H2 = 2.01588e-3;  % kg/mol

    % Membrana
        A_PEM = 1.0;   % m^2 (área activa)
        % Hidratación
            lambda_an = 14;
            lambda_ca = 10;
        % Energías de activación (del set que vienes usando)
            Eact_an = 76000; % J/mol
            Eact_ca = 18000; % J/mol
            Jref_an = 1744271.56;              % A/m^2  (pre-exponencial)
            Jref_ca = 1744271559894.8206;      % A/m^2  (pre-exponencial)



% ---------- Caso 1 ---------- T variable, J constante
T_vec_caso1 = [300 310 320 330 340 350]; % K (como tu tabla Ni)
J_fijo      = 5000;                     % A/m^2 (Ni validación T vs V)
L_mem_caso1  = 178e-6;                   % m (Ni: L=178e-6 para Fig 2c)

% ---------- Caso 2 ---------- T constante, Wdot variable
T_fijo      = 353;                      % K (típico Ni para J-V)
L_mem_caso2  = 50e-6;                    % m (Ni: L=50e-6 para Fig 2a)

% Potencia eléctrica de entrada (elige el barrido que quieras)
Wdot_vec_caso2 = [28.5782 82.733 169.04 884.3 1811.81 3717.94 5670.15 7647.88 9652.6 11660.04];  % W  (ejemplo)


if Caso == 1
    V_out   = zeros(size(T_vec_caso1));
    Wdot_out= zeros(size(T_vec_caso1));
    mdotH2  = zeros(size(T_vec_caso1));

    for k = 1:numel(T_vec_caso1)
        T = T_vec_caso1(k);

        V = Voltaje_PEM(T, J_fijo, L_mem_caso1, lambda_an, lambda_ca, ...
                            Jref_an, Jref_ca, Eact_an, Eact_ca, R, F);

        Wdot = A_PEM * J_fijo * V;
        
        NdotH2 = (J_fijo*A_PEM)/(2*F);
        mdot   = NdotH2 * M_H2;

        V_out(k)    = V;
        Wdot_out(k) = Wdot;
        mdotH2(k)   = mdot;
    end

elseif Caso == 2
    J_out   = zeros(size(Wdot_vec_caso2));
    V_out   = zeros(size(Wdot_vec_caso2));
    mdotH2  = zeros(size(Wdot_vec_caso2));

    for k = 1:numel(Wdot_vec_caso2)
        Wdot = Wdot_vec_caso2(k);

        [J, V] = solv_potencia_variable(T_fijo, Wdot, A_PEM, L_mem_caso2, lambda_an, lambda_ca, ...
                                    Jref_an, Jref_ca, Eact_an, Eact_ca, R, F);

        NdotH2 = (J*A_PEM)/(2*F);
        mdot   = NdotH2 * M_H2;

        J_out(k)  = J;
        V_out(k)  = V;
        mdotH2(k) = mdot;
    end
else
    error('Caso debe ser 1 o 2');
end


% --- Ni: Variación de temperatura (tabla izquierda de tu imagen)
T_Ni = [300 310 320 330 340 350];
V_Ni_T = [2.1728 2.1378 2.10406 2.07483 2.04819 2.02255];

% --- Ni: Variación de densidad de corriente (tabla derecha de tu imagen)
J_Ni = [0 20 50 100 500 1000 2000 3000 4000 5000 6000];
V_Ni_J = [1.22379 1.42891 1.65466 1.6904 1.7686 1.81181 1.85897 1.89005 1.91197 1.93052 1.94334];


% FUNCIONES
    function V = Voltaje_PEM(T, J, L_mem, lambda_an, lambda_ca, ...
                                 Jref_an, Jref_ca, Eact_an, Eact_ca, R, F)
        % Potencial reversible (forma simplificada usada en muchos modelos basados en Ni)
        V0 = 1.229 - 8.5e-4*(T - 298.15);

        % Corrientes de intercambio (Ni: J0 = Jref*exp(-Eact/RT))
        J0_an = Jref_an * exp(-Eact_an/(R*T));
        J0_ca = Jref_ca * exp(-Eact_ca/(R*T));

        % Resistencia de membrana: R_PEM = integral dx/sigma(lambda(x))
        R_PEM = ResArea_PEM(T, L_mem, lambda_an, lambda_ca);

        % Sobrepotenciales
        K = R*T/F;
        eta_act = K*( asinh(J/(2*J0_ca)) + asinh(J/(2*J0_an)) );
        eta_ohm = J * R_PEM;

        V = V0 + eta_act + eta_ohm;
    end

    function R_PEM = ResArea_PEM(T, L_mem, lambda_an, lambda_ca)
        n = 1001; % impar para Simpson 1/3
        x = linspace(0, L_mem, n);

        lambda_x = ((lambda_an - lambda_ca)/L_mem).*x + lambda_ca;

        sigma = (0.5139*lambda_x - 0.326).*exp(1268*(1/303 - 1/T)); % S/m

        inv_sigma = 1./sigma;

        R_PEM = simpson13(x, inv_sigma); % ohm*m^2
    end
    
    
    function I = simpson13(x, y)
        n = numel(x);
        if mod(n,2)==0
            error('Simpson 1/3 requiere numero IMPAR de puntos');
        end
        h = (x(end)-x(1))/(n-1);
        I = y(1) + y(end) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2));
        I = I*h/3;
    end

    function [J, V] = solv_potencia_variable(T, Wdot, A, L_mem, lambda_an, lambda_ca, ...
                                          Jref_an, Jref_ca, Eact_an, Eact_ca, R, F)
        J_in = 1000;
        tol = 1e-6;
        f  = @(J) A*J*Voltaje_PEM(T, J, L_mem, lambda_an, lambda_ca, ...
                                     Jref_an, Jref_ca, Eact_an, Eact_ca, R, F) - Wdot;

        % Derivada numérica estable
        df = @(J) (f(J + max(1e-6,1e-4*J)) - f(J - max(1e-6,1e-4*J))) / (2*max(1e-6,1e-4*J));

        while true
            raiz = J_in - f(J_in)/df(J_in);
            if abs(f(raiz)) < tol
                J = raiz;
                break;
            else
                J_in = raiz;
            end
        end
         V = Wdot/(A*J);
    end
