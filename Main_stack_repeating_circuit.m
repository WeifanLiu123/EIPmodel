clear all; clc;
%% Constants definition
F         = 96485;                  % Faraday constant
RT        = 8.3144 * 298;           % gas constant*temperature
V_T       = RT / F;                 % thermal voltage
E         = 300;                    % the pore ion correlation energy
C_St_vol  = 1.45e8;                 % F m-3, Zero-charge capacitance
Alfa      = 30;                     % Factor for charge dependence of Stern capacitance
L_elec    = 450e-6;                 % Electrode thickness
L_aem     = 500e-6;                 % Electrode thickness
p_pore    = 0.3;                    % pore volume fraction
p_IEP     = 0.6;                    % polymer volume fraction
X_poly    = 2500;                   % charge density
Eta_cse   = 0.1*(0.6)^(1.5);        % diffusion reduction factor
Eta_aem   = 0.1;                    % diffusion reduction factor
Di_Na_b   = 1.33 * 10^-9;           % m2 s-1, Na+ effective diffusion coefficient
Di_Cl_b   = 2.03 * 10^-9;           % m2 s-1, Cl- effective diffusion coefficient
Di_Na_cse = Eta_cse * Di_Na_b;      % m2 s-1, Na+ effective diffusion coefficient
Di_Cl_cse = Eta_cse * Di_Cl_b;      % m2 s-1, Cl- effective diffusion coefficient
Di_Na_aem = Eta_aem * Di_Na_b;      % m2 s-1, Na+ effective diffusion coefficient
Di_Cl_aem = Eta_aem * Di_Cl_b;      % m2 s-1, Cl- effective diffusion coefficient
z_Na      = 1;                      % Na+ valance
z_Cl      = -1;                     % Cl- valance
% operational conditions
c_di_0    = 100;                    % diluate concentration
c_con_0   = 100;                    % concentrate concentration
c_feed_di = c_di_0;                 % initial concentration in diluate tank
c_feed_con= c_con_0;                % initial concentration in concentrate tank
c_sp_di   = c_feed_di;              % initial concentration in diluate spacer channel
c_sp_con  = c_feed_con;             % initial concentration in concentrate spacer channel
n         = 2;                      % cycle number 
TT        = 20;                     % half cycle time
np        = 100;                    % position grid number
mp        = 100;                    % time grid number
dx        = L_elec / np;            % position grid
dTT       = TT / mp;                % time grid
L_spacer  = 5e-3;                   % m, thickness of spacer channel
dsp       = L_spacer/np;            % position grid
daem      = L_aem/np;               % position grid
Area      = 6e-4;                   % m2, electrode area
flowrate  = 5e-6/60;                % m3 s-1, flowrate
V_sp      = L_spacer*Area;          % m3, spacer channel volume
V_feed    = 30e-6-V_sp;             % m3, feed tank volume
tau_sp    = V_sp/flowrate;          % s, hydraulic retention time in the spacer
tau_feed  = V_feed /flowrate;       % s, hydraulic retention time in the tank
EER       = 0.0035;                 % external resistance
i_list    = 20*0.85;
cell_unit_list = zeros(1, length(i_list));
for mm = 1:length(i_list)           
    I_max_0   = i_list(mm);                     % A m-2, current density
        %% create list to store results
        % CSE 1
        c_IEP_Na_list1            = zeros(2*n*mp,np+1);
        c_IEP_Cl_list1            = zeros(2*n*mp,np+1);
        c_pore_Na_list1           = zeros(2*n*mp,np+1);
        c_pore_Cl_list1           = zeros(2*n*mp,np+1);
        c_v_list1                 = zeros(2*n*mp,np+1);
        phi_IEP_list1             = zeros(2*n*mp,np+1);
        phi_interfacial_list1     = zeros(2*n*mp,np+1);
        phi_Donnan_list1          = zeros(2*n*mp,np+1);
        phi_Stern_list1           = zeros(2*n*mp,np+1);
    
        % concentrate 1
        c_sp_con_list1            = zeros(2*n*mp,np+1);
        phi_sp_con_list1          = zeros(2*n*mp,np+1);
    
        % AEM 1
        c_AEM_Na_list1            = zeros(2*n*mp,np+1);
        c_AEM_Cl_list1            = zeros(2*n*mp,np+1);
        phi_AEM_di_list1          = zeros(2*n*mp,np+1);
    
        % diluate 1
        c_sp_di_list1             = zeros(2*n*mp,np+1);
        phi_sp_di_list1           = zeros(2*n*mp,np+1);
    
        % CSE 2
        c_IEP_Na_list2            = zeros(2*n*mp,np+1);
        c_IEP_Cl_list2            = zeros(2*n*mp,np+1);
        c_pore_Na_list2           = zeros(2*n*mp,np+1);
        c_pore_Cl_list2           = zeros(2*n*mp,np+1);
        c_v_list2                 = zeros(2*n*mp,np+1);
        phi_IEP_list2             = zeros(2*n*mp,np+1);
        phi_interfacial_list2     = zeros(2*n*mp,np+1);
        phi_Donnan_list2          = zeros(2*n*mp,np+1);
        phi_Stern_list2           = zeros(2*n*mp,np+1);
    
        % concentrate 2
        c_sp_con_list2            = zeros(2*n*mp,np+1);
        phi_sp_con_list2          = zeros(2*n*mp,np+1);
    
        % AEM 2
        c_AEM_Na_list2            = zeros(2*n*mp,np+1);
        c_AEM_Cl_list2            = zeros(2*n*mp,np+1);
        phi_AEM_di_list2          = zeros(2*n*mp,np+1);
    
        % diluate 2
        c_sp_di_list2             = zeros(2*n*mp,np+1);
        phi_sp_di_list2           = zeros(2*n*mp,np+1);
    
        phi_carbon_list1          = zeros(2*n*mp,1);
        j_ion_flux_list1          = zeros(2*n*mp,1);
        Aver_sigma_list1          = zeros(2*n*mp,1);
        j_Na_left_flux_list1      = zeros(2*n*mp,1);
        j_Na_right_flux_list1     = zeros(2*n*mp,1);
        j_Cl_left_flux_list1      = zeros(2*n*mp,1);
        j_Cl_right_flux_list1     = zeros(2*n*mp,1);
    
        phi_carbon_list2          = zeros(2*n*mp,1);
        j_ion_flux_list2          = zeros(2*n*mp,1);
        Aver_sigma_list2          = zeros(2*n*mp,1);
    
        j_Na_con_flux_list        = zeros(2*n*mp,1);
        j_Cl_con_flux_list        = zeros(2*n*mp,1);
        j_Na_di_flux_list        = zeros(2*n*mp,1);
        j_Cl_di_flux_list        = zeros(2*n*mp,1);
    
        j_Na_left_flux_list2      = zeros(2*n*mp,1);
        j_Na_right_flux_list2     = zeros(2*n*mp,1);
        j_Cl_left_flux_list2      = zeros(2*n*mp,1);
        j_Cl_right_flux_list2     = zeros(2*n*mp,1);
    
        cell_voltage_list        = zeros(2*n*mp,1);
    
        c_feed_di_list           = zeros(2*n*mp,1);
        c_feed_con_list          = zeros(2*n*mp,1);
        %%
        % Initialize
        c_pore_Na_ini = (0:np)*(c_con_0-c_di_0)/np+c_di_0;
        c_pore_Cl_ini = (0:np)*(c_con_0-c_di_0)/np+c_di_0;
        c_IEP_Na_ini = c_pore_Na_ini.*exp(-z_Na.*asinh(-X_poly./(2.*c_pore_Na_ini)));
        c_IEP_Cl_ini = c_pore_Cl_ini.*exp(-z_Cl*asinh(-X_poly./(2.*c_pore_Cl_ini)));
        c_v_ini =  (0:np)*(c_con_0-c_di_0)/np+c_di_0;
        phi_IEP_ini = asinh(-X_poly./(2.*c_pore_Na_ini));
        phi_interfacial_ini = asinh(X_poly./(2.*c_v_ini));
        phi_Donnan_ini = asinh(-(c_pore_Na_ini-c_pore_Cl_ini)./(2.*c_pore_Na_ini));
        phi_Stern_ini = -(c_pore_Na_ini-c_IEP_Cl_ini)./(C_St_vol.*V_T./F+Alfa.*(c_pore_Na_ini-c_IEP_Cl_ini).^2.*V_T./F);
    
        x0 = [c_IEP_Na_ini , c_IEP_Cl_ini, phi_IEP_ini,   0,0,0];
        options = optimoptions('fsolve', 'MaxFunEvals', 10000000, 'Maxiter', 10000000, 'Display', 'off','Algorithm', 'trust-region','FunctionTolerance',1e-10,'OptimalityTolerance',1e-10);
    
        [x, fval, exitflag] = fsolve(@(x) Donnandialysis(x, x0, np, dx, X_poly, Di_Na_cse, Di_Cl_cse, c_di_0,c_con_0,  z_Na, z_Cl), x0, options);
    
        c_IEP_Na_ini = x(1:np+1);
        c_IEP_Cl_ini = x(np+2:2*np+2);
        phi_IEP_ini =x(2*np+3:3*np+3);
        Na_ion_flux_donnan_dialysis = x(3*(np+1)+2);
        c_v_ini = sqrt( (c_IEP_Na_ini+c_IEP_Cl_ini).^2 -X_poly^2)./2;
        phi_interfacial_ini = asinh(X_poly./2./c_v_ini);
    
        x0 = [c_pore_Na_ini, c_pore_Cl_ini,  phi_Donnan_ini, phi_Stern_ini,0];% Initial conditions vector
    
        I_max =0;
        charging_discharge = 1;
        [x, fval, exitflag] = fsolve(@(x) deal_equations_ini(x, x0, np, V_T, F, C_St_vol, Alfa, z_Na, z_Cl,E,phi_IEP_ini,c_v_ini, phi_interfacial_ini), x0, options);
        disp(exitflag)

        c_pore_Na_ini = x(1:np+1);
        c_pore_Cl_ini = x(np+2:2*np+2);
        phi_Donnan_ini =  x(2*np+3:3*np+3);
        phi_Stern_ini  = x(3*np+4:4*np+4);
        phi_carbon_ini = x(4*(np+1)+1);
    
        j_ion_flux_ini = 0;
        Aver_sigma_ini = 0;
    
        c_con_ini   = c_con_0*ones(1,np+1);
        c_di_ini    = c_di_0*ones(1,np+1);
        phi_con_ini = zeros(1,np+1);
        phi_di_ini  = zeros(1,np+1);
        c_aem_Na_ini = c_IEP_Cl_ini;
        c_aem_Cl_ini = c_IEP_Na_ini;
        phi_aem_ini  = -phi_IEP_ini;

        c_di_interface  = c_di_0;
        c_con_interface  = c_di_0;
        c_di_interface1  = c_di_0;
        c_con_interface1 = c_con_0;
        c_di_interface2  = c_di_0;
        c_con_interface2 = c_con_0;

        %% main
        x0 = [c_IEP_Na_ini , c_IEP_Cl_ini, c_pore_Na_ini, c_pore_Cl_ini, c_v_ini, phi_IEP_ini, phi_interfacial_ini, phi_Donnan_ini, phi_Stern_ini, c_IEP_Na_ini , c_IEP_Cl_ini, c_pore_Na_ini, c_pore_Cl_ini, c_v_ini, phi_IEP_ini, phi_interfacial_ini, phi_Donnan_ini, phi_Stern_ini, phi_carbon_ini, j_ion_flux_ini,Aver_sigma_ini,0,0,0,0 ,phi_carbon_ini, j_ion_flux_ini,Aver_sigma_ini,0,0,0,0 ];% Initial conditions vector
        x1 = [c_con_ini,c_con_ini,phi_con_ini,c_aem_Na_ini,c_aem_Cl_ini,phi_aem_ini,c_di_ini,c_di_ini ,phi_di_ini,0,0,0,0,0,0 ];
        x2 = [c_con_ini,c_con_ini,phi_con_ini,c_aem_Na_ini,c_aem_Cl_ini,phi_aem_ini,c_di_ini,c_di_ini ,phi_di_ini,0,0,0,0,0,0 ];
        kk = 18;
       %% Initialize
        cell_voltage_inital = zeros(mp/2,1);
        phi_carbon_list1_inital = zeros(mp/2,1);
        phi_carbon_list2_inital = zeros(mp/2,1);
        c_IEP_Na_list1_inital            = zeros(mp/2,np+1);
        c_IEP_Cl_list1_inital            = zeros(mp/2,np+1);
        c_IEP_Na_list2_inital            = zeros(mp/2,np+1);
        c_IEP_Cl_list2_inital            = zeros(mp/2,np+1);
        c_pore_Na_list1_inital           = zeros(mp/2,np+1);
        c_pore_Cl_list1_inital           = zeros(mp/2,np+1);
        c_v_list1_inital                 = zeros(mp/2,np+1);
        c_pore_Na_list2_inital           = zeros(mp/2,np+1);
        c_pore_Cl_list2_inital           = zeros(mp/2,np+1);
        c_v_list2_inital                 = zeros(mp/2,np+1);

        charging_discharge = 1;
        for j = 1:mp/2
            I_max = I_max_0 *(1-exp(-j*30/mp*2)-exp(-(mp/2-j)*30/mp*2));
            %% CSE unit
                [x, fval, exitflag] = fsolve(@(x) stack_unit(x, x0, charging_discharge, np, dx,dsp,  daem, dTT,X_poly, V_T, F, C_St_vol, p_pore, p_IEP, Di_Na_cse, Di_Cl_cse,Di_Na_b, Di_Cl_b, I_max, L_elec, c_di_0,c_con_0,  Alfa, z_Na, z_Cl,E,tau_sp, kk,c_di_interface1,c_con_interface1,c_di_interface2,c_con_interface2 ), x0, options);
                disp(exitflag)
                x0 = x;
                c_IEP_Na_list1_inital(j,1:np+1)            = x(0*(np+1)+1:1*(np+1));
                c_IEP_Cl_list1_inital(j,1:np+1)            = x(1*(np+1)+1:2*(np+1));
                c_pore_Na_list1_inital(j,1:np+1)           = x(2*(np+1)+1:3*(np+1));
                c_pore_Cl_list1_inital(j,1:np+1)           = x(3*(np+1)+1:4*(np+1));
                c_v_list1_inital(j,1:np+1)                 = x(4*(np+1)+1:5*(np+1));

                c_IEP_Na_list2_inital(j,1:np+1)            = x(9*(np+1)+1:10*(np+1));
                c_IEP_Cl_list2_inital(j,1:np+1)            = x(10*(np+1)+1:11*(np+1));
                c_pore_Na_list2_inital(j,1:np+1)           = x(11*(np+1)+1:12*(np+1));
                c_pore_Cl_list2_inital(j,1:np+1)           = x(12*(np+1)+1:13*(np+1));
                c_v_list2_inital(j,1:np+1)                 = x(13*(np+1)+1:14*(np+1));


                phi_carbon_list1_inital(j)                = x(kk*(np+1)+1)+x2(9*(np+1));            
                phi_carbon_list2_inital(j)                = x(kk*(np+1)+8);
                j_Na_left_flux                      =  x(kk*(np+1)+4);
                j_Na_right_flux                     =  x(kk*(np+1)+5);
                j_Cl_left_flux                      =  x(kk*(np+1)+6);
                j_Cl_right_flux                     =  x(kk*(np+1)+7);

                cell_voltage_inital(j)                    = (phi_carbon_list2_inital(j) - phi_carbon_list1_inital(j)+2*EER*I_max_0*F/RT)*RT/F;

            %% Spacer channel unit
                [x, fval, exitflag] = fsolve(@(x) spacer_unit1(x, x1, charging_discharge, np, dx,dsp,  daem, dTT,X_poly, V_T, F, C_St_vol, p_pore, p_IEP, Di_Na_aem, Di_Cl_aem,Di_Na_b, Di_Cl_b, I_max, L_elec, c_feed_di,c_feed_con,  Alfa, z_Na, z_Cl,E,tau_sp, kk, j_Na_left_flux, j_Na_right_flux,j_Cl_left_flux, j_Cl_right_flux), x1, options);
                disp(exitflag)
                x1 = x;
                c_di_interface1                                 = x(7*(np+1));
                c_con_interface1                                = x(1);

                [x, fval, exitflag] = fsolve(@(x) spacer_unit2(x, x2, charging_discharge, np, dx,dsp,  daem, dTT,X_poly, V_T, F, C_St_vol, p_pore, p_IEP, Di_Na_aem, Di_Cl_aem,Di_Na_b, Di_Cl_b, I_max, L_elec, c_feed_di,c_feed_con,  Alfa, z_Na, z_Cl,E,tau_sp, kk, j_Na_left_flux, j_Na_right_flux,j_Cl_left_flux, j_Cl_right_flux), x2, options);
                disp(exitflag)
                x2 = x;
                c_di_interface2                                 = x(7*(np+1));
                c_con_interface2                                = x(1);

        end

        for  ii = 1:n
            charging_discharge = 2;
    
            for j = 1:mp
                I_max = I_max_0 *(1-exp(-j*30/mp)-exp(-(100-j)*30/mp));
                %% CSE unit
                [x, fval, exitflag] = fsolve(@(x) stack_unit(x, x0, charging_discharge, np, dx,dsp,  daem, dTT,X_poly, V_T, F, C_St_vol, p_pore, p_IEP, Di_Na_cse, Di_Cl_cse,Di_Na_b, Di_Cl_b, I_max, L_elec, c_di_0,c_con_0,  Alfa, z_Na, z_Cl,E,tau_sp, kk,c_di_interface1,c_con_interface1,c_di_interface2,c_con_interface2 ), x0, options);
                % disp(exitflag)
                x0 = x;
                % CSE 1
                c_IEP_Na_list1((ii-1)*2*mp+j,1:np+1)            = x(0*(np+1)+1:1*(np+1));
                c_IEP_Cl_list1((ii-1)*2*mp+j,1:np+1)            = x(1*(np+1)+1:2*(np+1));
                c_pore_Na_list1((ii-1)*2*mp+j,1:np+1)           = x(2*(np+1)+1:3*(np+1));
                c_pore_Cl_list1((ii-1)*2*mp+j,1:np+1)           = x(3*(np+1)+1:4*(np+1));
                c_v_list1((ii-1)*2*mp+j,1:np+1)                 = x(4*(np+1)+1:5*(np+1));
                phi_IEP_list1((ii-1)*2*mp+j,1:np+1)             = x(5*(np+1)+1:6*(np+1));
                phi_interfacial_list1((ii-1)*2*mp+j,1:np+1)     = x(6*(np+1)+1:7*(np+1));
                phi_Donnan_list1((ii-1)*2*mp+j,1:np+1)          = x(7*(np+1)+1:8*(np+1));
                phi_Stern_list1((ii-1)*2*mp+j,1:np+1)           = x(8*(np+1)+1:9*(np+1));
                % CSE 2
                c_IEP_Na_list2((ii-1)*2*mp+j,1:np+1)            = x(9*(np+1)+1:10*(np+1));
                c_IEP_Cl_list2((ii-1)*2*mp+j,1:np+1)            = x(10*(np+1)+1:11*(np+1));
                c_pore_Na_list2((ii-1)*2*mp+j,1:np+1)           = x(11*(np+1)+1:12*(np+1));
                c_pore_Cl_list2((ii-1)*2*mp+j,1:np+1)           = x(12*(np+1)+1:13*(np+1));
                c_v_list2((ii-1)*2*mp+j,1:np+1)                 = x(13*(np+1)+1:14*(np+1));
                phi_IEP_list2((ii-1)*2*mp+j,1:np+1)             = x(14*(np+1)+1:15*(np+1))+x1(9*(np+1));
                phi_interfacial_list2((ii-1)*2*mp+j,1:np+1)     = x(15*(np+1)+1:16*(np+1));
                phi_Donnan_list2((ii-1)*2*mp+j,1:np+1)          = x(16*(np+1)+1:17*(np+1));
                phi_Stern_list2((ii-1)*2*mp+j,1:np+1)           = x(17*(np+1)+1:18*(np+1));
    
    
                phi_carbon_list1((ii-1)*2*mp+j)                 = x(kk*(np+1)+1);
                j_ion_flux_list1((ii-1)*2*mp+j)                 = x(kk*(np+1)+2);
                Aver_sigma_list1((ii-1)*2*mp+j)                 = x(kk*(np+1)+3);
                j_Na_left_flux_list1((ii-1)*2*mp+j)             = x(kk*(np+1)+4);
                j_Na_right_flux_list1((ii-1)*2*mp+j)            = x(kk*(np+1)+5);
                j_Cl_left_flux_list1((ii-1)*2*mp+j)             = x(kk*(np+1)+6);
                j_Cl_right_flux_list1((ii-1)*2*mp+j)            = x(kk*(np+1)+7);
    
                phi_carbon_list2((ii-1)*2*mp+j)                 = x(kk*(np+1)+8)+x1(9*(np+1));
                j_ion_flux_list2((ii-1)*2*mp+j)                 = x(kk*(np+1)+9);
                Aver_sigma_list2((ii-1)*2*mp+j)                 = x(kk*(np+1)+10);
                j_Na_left_flux_list2((ii-1)*2*mp+j)             = x(kk*(np+1)+11);
                j_Na_right_flux_list2((ii-1)*2*mp+j)            = x(kk*(np+1)+12);
                j_Cl_left_flux_list2((ii-1)*2*mp+j)             = x(kk*(np+1)+13);
                j_Cl_right_flux_list2((ii-1)*2*mp+j)            = x(kk*(np+1)+14);
        
                j_Na_left_flux                      =  x(kk*(np+1)+4);
                j_Na_right_flux                     =  x(kk*(np+1)+5);
                j_Cl_left_flux                      =  x(kk*(np+1)+6);
                j_Cl_right_flux                     =  x(kk*(np+1)+7);

                %% Spacer channel unit
                [x, fval, exitflag] = fsolve(@(x) spacer_unit1(x, x1, charging_discharge, np, dx,dsp,  daem, dTT,X_poly, V_T, F, C_St_vol, p_pore, p_IEP, Di_Na_aem, Di_Cl_aem,Di_Na_b, Di_Cl_b, I_max, L_elec, c_feed_di,c_feed_con,  Alfa, z_Na, z_Cl,E,tau_sp, kk, j_Na_left_flux, j_Na_right_flux,j_Cl_left_flux, j_Cl_right_flux), x1, options);
                % disp(exitflag)
                x1 = x;
    
                % concentrate 1
                c_sp_con_list1((ii-1)*2*mp+j,1:np+1)            = x(0*(np+1)+1:1*(np+1));
                phi_sp_con_list1((ii-1)*2*mp+j,1:np+1)          = x(2*(np+1)+1:3*(np+1));
    
                % AEM 1
                c_AEM_Na_list1((ii-1)*2*mp+j,1:np+1)            = x(3*(np+1)+1:4*(np+1));
                c_AEM_Cl_list1((ii-1)*2*mp+j,1:np+1)            = x(4*(np+1)+1:5*(np+1));
                phi_AEM_di_list1((ii-1)*2*mp+j,1:np+1)          = x(5*(np+1)+1:6*(np+1));
                %
                % diluate 1
                c_sp_di_list1((ii-1)*2*mp+j,1:np+1)             = x(6*(np+1)+1:7*(np+1));
                phi_sp_di_list1((ii-1)*2*mp+j,1:np+1)           = x(8*(np+1)+1:9*(np+1));
                %
    
                c_di_interface1                                 = x(7*(np+1));
                c_con_interface1                                = x(1);
    
                [x, fval, exitflag] = fsolve(@(x) spacer_unit2(x, x2, charging_discharge, np, dx,dsp,  daem, dTT,X_poly, V_T, F, C_St_vol, p_pore, p_IEP, Di_Na_aem, Di_Cl_aem,Di_Na_b, Di_Cl_b, I_max, L_elec, c_feed_di,c_feed_con,  Alfa, z_Na, z_Cl,E,tau_sp, kk, j_Na_left_flux, j_Na_right_flux,j_Cl_left_flux, j_Cl_right_flux), x2, options);
                % disp(exitflag)
                x2 = x;
    
                % concentrate 2
                c_sp_con_list2((ii-1)*2*mp+j,1:np+1)            = x(0*(np+1)+1:1*(np+1));
                phi_sp_con_list2((ii-1)*2*mp+j,1:np+1)          = x(2*(np+1)+1:3*(np+1));
    
                % AEM 2
                c_AEM_Na_list2((ii-1)*2*mp+j,1:np+1)            = x(3*(np+1)+1:4*(np+1));
                c_AEM_Cl_list2((ii-1)*2*mp+j,1:np+1)            = x(4*(np+1)+1:5*(np+1));
                phi_AEM_di_list2((ii-1)*2*mp+j,1:np+1)          = x(5*(np+1)+1:6*(np+1));
    
                % diluate 2
                c_sp_di_list2((ii-1)*2*mp+j,1:np+1)             = x(6*(np+1)+1:7*(np+1));
                phi_sp_di_list2((ii-1)*2*mp+j,1:np+1)           = x(8*(np+1)+1:9*(np+1));
    
                c_di_interface2                                 = x(7*(np+1));
                c_con_interface2                                = x(1);
    
                
    
                cell_voltage_list((ii-1)*2*mp+j)                = phi_carbon_list1((ii-1)*2*mp+j) - phi_carbon_list2((ii-1)*2*mp+j)+2*EER*I_max_0*F/RT;
                % concentrate 2
                phi_sp_con_list2((ii-1)*2*mp+j,1:np+1)          = phi_sp_con_list2((ii-1)*2*mp+j,1:np+1) + x0(15*(np+1))+x1(9*(np+1))+asinh(X_poly/(2*c_sp_con_list2((ii-1)*2*mp+j,1)));
                % AEM 2
                phi_AEM_di_list2((ii-1)*2*mp+j,1:np+1)          = phi_AEM_di_list2((ii-1)*2*mp+j,1:np+1) + x0(15*(np+1))+x1(9*(np+1))+asinh(X_poly/(2*c_sp_con_list2((ii-1)*2*mp+j,1)));
                % diluate 2
                phi_sp_di_list2((ii-1)*2*mp+j,1:np+1)           = phi_sp_di_list2((ii-1)*2*mp+j,1:np+1) + x0(15*(np+1))+x1(9*(np+1))+asinh(X_poly/(2*c_sp_con_list2((ii-1)*2*mp+j,1)));
            end
    
            charging_discharge = 1;
    
            for j = mp+1:2*mp
                I_max = - I_max_0 *(exp(-(j-100)*30/mp)+exp(-(200-j)*30/mp)-1);
                %% CSE unit
                [x, fval, exitflag] = fsolve(@(x) stack_unit(x, x0, charging_discharge, np, dx,dsp,  daem, dTT,X_poly, V_T, F, C_St_vol, p_pore, p_IEP, Di_Na_cse, Di_Cl_cse,Di_Na_b, Di_Cl_b, I_max, L_elec, c_di_0,c_con_0,  Alfa, z_Na, z_Cl,E,tau_sp, kk,c_di_interface1,c_con_interface1,c_di_interface2,c_con_interface2 ), x0, options);
                % disp(exitflag)
                x0 = x;
                % CSE 1
                c_IEP_Na_list1((ii-1)*2*mp+j,1:np+1)            = x(0*(np+1)+1:1*(np+1));
                c_IEP_Cl_list1((ii-1)*2*mp+j,1:np+1)            = x(1*(np+1)+1:2*(np+1));
                c_pore_Na_list1((ii-1)*2*mp+j,1:np+1)           = x(2*(np+1)+1:3*(np+1));
                c_pore_Cl_list1((ii-1)*2*mp+j,1:np+1)           = x(3*(np+1)+1:4*(np+1));
                c_v_list1((ii-1)*2*mp+j,1:np+1)                 = x(4*(np+1)+1:5*(np+1));
                phi_IEP_list1((ii-1)*2*mp+j,1:np+1)             = x(5*(np+1)+1:6*(np+1))+x2(9*(np+1));
                phi_interfacial_list1((ii-1)*2*mp+j,1:np+1)     = x(6*(np+1)+1:7*(np+1));
                phi_Donnan_list1((ii-1)*2*mp+j,1:np+1)          = x(7*(np+1)+1:8*(np+1));
                phi_Stern_list1((ii-1)*2*mp+j,1:np+1)           = x(8*(np+1)+1:9*(np+1));
    
                % CSE 2
                c_IEP_Na_list2((ii-1)*2*mp+j,1:np+1)            = x(9*(np+1)+1:10*(np+1));
                c_IEP_Cl_list2((ii-1)*2*mp+j,1:np+1)            = x(10*(np+1)+1:11*(np+1));
                c_pore_Na_list2((ii-1)*2*mp+j,1:np+1)           = x(11*(np+1)+1:12*(np+1));
                c_pore_Cl_list2((ii-1)*2*mp+j,1:np+1)           = x(12*(np+1)+1:13*(np+1));
                c_v_list2((ii-1)*2*mp+j,1:np+1)                 = x(13*(np+1)+1:14*(np+1));
                phi_IEP_list2((ii-1)*2*mp+j,1:np+1)             = x(14*(np+1)+1:15*(np+1));
                phi_interfacial_list2((ii-1)*2*mp+j,1:np+1)     = x(15*(np+1)+1:16*(np+1));
                phi_Donnan_list2((ii-1)*2*mp+j,1:np+1)          = x(16*(np+1)+1:17*(np+1));
                phi_Stern_list2((ii-1)*2*mp+j,1:np+1)           = x(17*(np+1)+1:18*(np+1));
    
                phi_carbon_list1((ii-1)*2*mp+j)                 = x(kk*(np+1)+1)+x2(9*(np+1));
                j_ion_flux_list1((ii-1)*2*mp+j)                 = x(kk*(np+1)+2);
                Aver_sigma_list1((ii-1)*2*mp+j)                 = x(kk*(np+1)+3);
                j_Na_left_flux_list1((ii-1)*2*mp+j)             = x(kk*(np+1)+4);
                j_Na_right_flux_list1((ii-1)*2*mp+j)            = x(kk*(np+1)+5);
                j_Cl_left_flux_list1((ii-1)*2*mp+j)             = x(kk*(np+1)+6);
                j_Cl_right_flux_list1((ii-1)*2*mp+j)            = x(kk*(np+1)+7);
    
                phi_carbon_list2((ii-1)*2*mp+j)                 = x(kk*(np+1)+8);
                j_ion_flux_list2((ii-1)*2*mp+j)                 = x(kk*(np+1)+9);
                Aver_sigma_list2((ii-1)*2*mp+j)                 = x(kk*(np+1)+10);
                j_Na_left_flux_list2((ii-1)*2*mp+j)             = x(kk*(np+1)+11);
                j_Na_right_flux_list2((ii-1)*2*mp+j)            = x(kk*(np+1)+12);
                j_Cl_left_flux_list2((ii-1)*2*mp+j)             = x(kk*(np+1)+13);
                j_Cl_right_flux_list2((ii-1)*2*mp+j)            = x(kk*(np+1)+14);

                j_Na_left_flux                      =  x(kk*(np+1)+4);
                j_Na_right_flux                     =  x(kk*(np+1)+5);
                j_Cl_left_flux                      =  x(kk*(np+1)+6);
                j_Cl_right_flux                     =  x(kk*(np+1)+7);
                %% Spacer channel unit
                [x, fval, exitflag] = fsolve(@(x) spacer_unit1(x, x1, charging_discharge, np, dx,dsp,  daem, dTT,X_poly, V_T, F, C_St_vol, p_pore, p_IEP, Di_Na_aem, Di_Cl_aem,Di_Na_b, Di_Cl_b, I_max, L_elec, c_feed_di,c_feed_con,  Alfa, z_Na, z_Cl,E,tau_sp, kk, j_Na_left_flux, j_Na_right_flux,j_Cl_left_flux, j_Cl_right_flux), x1, options);
                % disp(exitflag)
                x1 = x;
    
                % concentrate 1
                c_sp_con_list1((ii-1)*2*mp+j,1:np+1)            = x(0*(np+1)+1:1*(np+1));
                phi_sp_con_list1((ii-1)*2*mp+j,1:np+1)          = x(2*(np+1)+1:3*(np+1));
    
                % AEM 1
                c_AEM_Na_list1((ii-1)*2*mp+j,1:np+1)            = x(3*(np+1)+1:4*(np+1));
                c_AEM_Cl_list1((ii-1)*2*mp+j,1:np+1)            = x(4*(np+1)+1:5*(np+1));
                phi_AEM_di_list1((ii-1)*2*mp+j,1:np+1)          = x(5*(np+1)+1:6*(np+1));
    
                % diluate 1
                c_sp_di_list1((ii-1)*2*mp+j,1:np+1)             = x(6*(np+1)+1:7*(np+1));
                phi_sp_di_list1((ii-1)*2*mp+j,1:np+1)           = x(8*(np+1)+1:9*(np+1));
                %
    
                c_di_interface1                                 = x(7*(np+1));
                c_con_interface1                                = x(1);
    
                [x, fval, exitflag] = fsolve(@(x) spacer_unit2(x, x2, charging_discharge, np, dx,dsp,  daem, dTT,X_poly, V_T, F, C_St_vol, p_pore, p_IEP, Di_Na_aem, Di_Cl_aem,Di_Na_b, Di_Cl_b, I_max, L_elec, c_feed_di,c_feed_con,  Alfa, z_Na, z_Cl,E,tau_sp, kk, j_Na_left_flux, j_Na_right_flux,j_Cl_left_flux, j_Cl_right_flux), x2, options);
                % disp(exitflag)
                x2 = x;
 
                % concentrate 2
                c_sp_con_list2((ii-1)*2*mp+j,1:np+1)            = x(0*(np+1)+1:1*(np+1));
                phi_sp_con_list2((ii-1)*2*mp+j,1:np+1)          = x(2*(np+1)+1:3*(np+1));
    
                % AEM 2
                c_AEM_Na_list2((ii-1)*2*mp+j,1:np+1)            = x(3*(np+1)+1:4*(np+1));
                c_AEM_Cl_list2((ii-1)*2*mp+j,1:np+1)            = x(4*(np+1)+1:5*(np+1));
                phi_AEM_di_list2((ii-1)*2*mp+j,1:np+1)          = x(5*(np+1)+1:6*(np+1));
    
                % diluate 2
                c_sp_di_list2((ii-1)*2*mp+j,1:np+1)             = x(6*(np+1)+1:7*(np+1));
                phi_sp_di_list2((ii-1)*2*mp+j,1:np+1)           = x(8*(np+1)+1:9*(np+1));
    
                c_di_interface2                                 = x(7*(np+1));
                c_con_interface2                                = x(1);
               
                cell_voltage_list((ii-1)*2*mp+j)                = phi_carbon_list2((ii-1)*2*mp+j) - phi_carbon_list1((ii-1)*2*mp+j)+2*EER*I_max_0*F/RT;
    
                % concentrate 1
                phi_sp_con_list1((ii-1)*2*mp+j,1:np+1)          = phi_sp_con_list1((ii-1)*2*mp+j,1:np+1) + phi_sp_di_list1((ii-1)*2*mp+j,np+1) - phi_sp_con_list1((ii-1)*2*mp+j,1)  + phi_IEP_list2((ii-1)*2*mp+j,1) +asinh(X_poly/(2*c_sp_di_list1((ii-1)*2*mp+j,np+1)));
                % AEM 1
                phi_AEM_di_list1((ii-1)*2*mp+j,1:np+1)          = phi_AEM_di_list1((ii-1)*2*mp+j,1:np+1) + phi_sp_di_list1((ii-1)*2*mp+j,np+1) - phi_sp_con_list1((ii-1)*2*mp+j,1) +phi_IEP_list2((ii-1)*2*mp+j,1) +asinh(X_poly/(2*c_sp_di_list1((ii-1)*2*mp+j,np+1)));
                % diluate 1
                phi_sp_di_list1((ii-1)*2*mp+j,1:np+1)           = phi_sp_di_list1((ii-1)*2*mp+j,1:np+1) + phi_sp_di_list1((ii-1)*2*mp+j,np+1) - phi_sp_con_list1((ii-1)*2*mp+j,1) + phi_IEP_list2((ii-1)*2*mp+j,1) +asinh(X_poly/(2*c_sp_di_list1((ii-1)*2*mp+j,np+1)));
    
            end
    
        end
    cell_unit = cell_voltage_list*RT/F;
    cell_unit_list(mm) = mean(cell_voltage_list)*RT/F;
end




