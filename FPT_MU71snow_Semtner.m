% This code solves the FPT equation using the Semtner model for ice growth.
% A fully implicit method is used for solving the equation.
% There is a mix of non-dimensional and dimensional quantities -- Semtner
% Three-Layer model is solved using dimenional variables, and the scaled result is
% used in the FPT equation. (SI units wherever dimensional.)

clear all; clc;

T = 20;             % Total duration of the run
xi = 0.17;          % Fraction of net shortwave radiation penetrating ice.
snowlayer = input('Enter 1 for snow layer and 0 for none.');
    
delt = 0.025; Tyear = 30; T = T * Tyear;  Nt = T/delt; Nyear = Tyear/delt;                        % Non-dimensional time variables
Hmin = 0; Hmax = 10; Nh = 1200; delh = (Hmax - Hmin)/Nh; h = linspace(Hmin,Hmax,Nh);              % Non-dimensional thickness variables
Heq = 1.5; hd = h * Heq; delhd = Heq * delh; Nh1 = round(0.5/delhd);                              % Thickness scale and dimensional thickness
t_plot = delt; Nplot = round(T/t_plot); Navg1 = Tyear/t_plot; Navg2 = Nplot - Navg1;

NAug = round(19.17/delt); NOct = round(25/delt); 
NNov = NOct + 1; NDec = round(30/delt); NJan = 2; 
NApr = round(10/delt); NMay1 = NApr + 1; NMay2 = round(12.5/delt);

delt_dim = delt * (12 * 24 * 3600);             % Dimensional time step for integrating equation for open water energy balance.

Hcd = 0.1; Hc = Hcd/Heq;                    % Cut-off thickness to determine open water. Currently 10 cm.
delhc = Hc;
Nc = round(Hc/delh);

% Properties of sea ice (all dimensional)
rhoi = 900; ki = 2.03; rhocpi = 1882800;                 % lambda is (extinction coefficient)^-1 for Beer's law 
qI = 301248000; kappa = ki/(rhocpi);
qB = 267776000;

f0 = kappa/Heq;                                                                 % Scale for growth rate

% Properties of snow
rhos = 330; ks = 0.31; rhocps = 690360; qS = 109620800;

% Time scales (all dimensional)
L = 10^5; U = 0.1; tm = L/U;                                                    % "Inertial" time scale
td = Heq^2/kappa;                                                               % Diffusion time scale

tau = tm/td;                                                                    % Ratio of the two time scales. Should be < 1.

% Albedo using EW09
ai = 0.68; aw = 0.2; lambda = 0.67;
a = (ai + aw)/2 + (aw - ai)/2 * tanh(-hd/lambda);

% Various fluxes from MU-71
Fs = 15.9 * [0 0 1.9 9.9 17.7 19.2 13.6 9.0 3.7 0.4 0 0];                           % Shortwave 
Fl = 15.9 * [10.4 10.3 10.3 11.6 15.1 18.0 19.1 18.7 16.5 13.9 11.2 10.9];          % Incoming longwave
Fsh = 15.9 * [1.18 0.76 0.72 0.29 -0.45 -0.39 -0.30 -0.40 -0.17 0.10 0.56 0.79];    % Turbulent specific heat
Flh = 15.9 * [0 -0.02 -0.03 -0.09 -0.46 -0.70 -0.64 -0.66 -0.39 -0.19 -0.01 -0.01]; % Turbulent latent heat
as = [0 0 0.83 0.81 0.82 0.78 0.64 0.69 0.84 0.85 0 0];
Fb = input('Bottom heat flux in W/m^2');
DelF = input('Greenhouse-gas forcing in W/m^2');

% Interpolating the fluxes
t0 = (1.25:2.5:28.75);                                                             % The energy values are taken to be at 15th of each month...
ti = linspace(0,Tyear,Nyear);
Fsi = interp1([0 t0 Tyear], Fs([end 1:end 1]), ti, 'pchip');
Fli = interp1([0 t0 Tyear], Fl([end 1:end 1]), ti, 'pchip');
Fshi = interp1([0 t0 Tyear], Fsh([end 1:end 1]), ti, 'pchip');
Flhi = interp1([0 t0 Tyear], Flh([end 1:end 1]), ti, 'pchip');
asi = interp1([0 t0 Tyear], as([end 1:end 1]), ti, 'pchip');

Fli = Fli + DelF;

% Stefan-Boltzmann law (all dimensional, from EW09)
sigma = 5.67 * 10^(-8);                                                         % These are constants for linearized Stefan-Boltzmann law

% Initial condition
epsilon = 0.046; k1 = 0.048; k2 = 0.025;           % The values of k1 and k2 are dimensionless.
q = epsilon/k2; H = k2/k1;
g0 = h.^q.*exp(-h/H); g0 = g0/trapz(h,g0);
g = g0';

% Initial fraction of open water
Af = 1 - trapz(h,g);

% Arrays
A = zeros(Nh,Nh); hmean = zeros(Nplot+1,1);
Tsrf = zeros(Nh,1); f = zeros(Nh,1);
fmean = zeros(Nplot+1,1); Fsw = zeros(Nh,1); Flw = zeros(Nh,1);
a_mean = zeros(Nplot+1,1); gplot = zeros(Nplot+1,Nh);
f_plot = zeros(Nplot+1,Nh); Tplot = zeros(Nplot+1,1); Twater = zeros(Nplot+1,1);

cgrid = delt/delh; cgridc = delt/delhc;

% Temperature variables

Tsnow = 273 * ones(Nh,1); T0 = Tsnow; TI = Tsnow;
T1 = Tsnow; T2 = Tsnow;

TB = 271; % Equilibrium temperature at the ice-ocean interface

count = 1;

% Snow accumulation from MU71

hsp = zeros(Nyear,1); hsnow = zeros(Nh,1);

%20th Aug -- 30th Oct

if(snowlayer == 1)
    
    snowfactor = 1; %input('Factor by which the snowfall is reduced in comparison with the MU71 standard case');
    
    slope1 = 0.3/(NOct - NAug + 1)/snowfactor;
    slope2 = 0.05/(NNov - NApr + 1)/snowfactor;
    slope3 = 0.05/(NMay2 - NMay1 + 1)/snowfactor;

    for i = NAug:NOct
    
        hsp(i) = slope1 * (i - NAug + 1);
    
    end

    % 1st Nov -- 30th Dec

    for i = NNov:NDec
    
        hsp(i) = hsp(NOct) + slope2 * (i - NNov + 1);

    end

    % 1 Jan -- 30th Apr

    hsp(1) = hsp(Nyear) + slope2;

    for i = NJan:NApr
    
        hsp(i) = hsp(1) + slope2 * (i - NJan + 1);
    
    end

    %1st May -- 30th May

    for i = NMay1:NMay2
    
        hsp(i) = hsp(NApr) + slope3 * (i - NMay1 + 1);
    
    end
    
end

snow_count = zeros(Nyear,1);

% Time loop begins

for t = 0:Nt
    
    t1 = delt * t;
    tp = mod(t,Nyear);
    tp = tp + 2;
    
    if(tp == Nyear+1)
        
        tp = 1;
        
    elseif(tp == Nyear+2)
        
        tp = 2;
        
    end
    
    Ti_avg = dot(TI,g) * delh;
    Ts_avg = mean(Tsnow);
    
    snow_count(tp) = 0;
    
    if(hsp(tp) > 0)
        
        if((Ti_avg < 273) || (Ts_avg < 273))
        
            hsnow = hsp(tp) * ones(Nh,1);               % Snow is accumulated only if surface temperature is below freezing.
            
            snow_count(tp) = 1;
            
        end
        
    end
    
    asnow = asi(tp);
    
    % This loop is for ice of thickness upto 50 cm.
    
    for i = 1:Nh1
        
        hstemp = hsnow(i);
        
        if(hstemp == 0)
            
            Fa = -(1 - a(i)) * (1 - xi) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * TI(i)^4;
            Tdiff = (-hd(i) * Fa + ki * (TB - TI(i)))/(4 * sigma * TI(i)^3 * hd(i) + ki);
            
            Ttemp = TI(i) + Tdiff;
            
            if(Ttemp < 273)
                
                TI(i) = Ttemp;
                Fs = -(1 - a(i)) * (1 - xi) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * TI(i)^4;
                f(i) = (Fs - Fb)/qB; 
                
            elseif(Ttemp >= 273)
                
                TI(i) = 273;
                Fa = -(1 - a(i)) * (1 - xi) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * TI(i)^4;
                
                if(hd(i) == 0)
                    
                    Fs = 0;
                    f(i) = (Fa - Fs)/qB;
                    
                else
                
                    Fs = ki * (TB - TI(i))/hd(i);
                    f(i) = (Fa - Fs)/qI + (Fs - Fb)/qB;
                    
                end
            end
            
        elseif(hstemp > 0)
            
            beta1 = ks * hd(i)/(ki * hstemp + ks * hd(i)); beta2 = ki * hstemp/(ki * hstemp + ks * hd(i));
            Fa = -(1 - asnow) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * Tsnow(i)^4;
            Tdiff = (-Fa * hstemp + ks * ((beta1 - 1) * Tsnow(i) + beta2 * TB))/(4 * sigma * hstemp * Tsnow(i)^3 + ks * (1 - beta1));
            Ttemp = Tsnow(i) + Tdiff;
            
            if(Ttemp < 273)
                
                Tsnow(i) = Ttemp;
                Fs = -(1 - asnow) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * Tsnow(i)^4;
                f(i) = (Fs - Fb)/qB;
                
            elseif(Ttemp >= 273)
                
                Tsnow(i) = 273;
                Fs = ks * ((beta1 - 1) * Tsnow(i) + beta2 * TB)/hstemp;
                f(i) = (Fs - Fb)/qB;
                
                Fa = -(1 - asnow) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * Tsnow(i)^4;
                hsnow(i) = hsnow(i) + delt_dim * (Fa - Fs)/qS;
                
            end
            
            if(hsnow(i) < 0)
                    
                hsnow(i) = 0;
                
            elseif(hsnow(i) > 0)
                
                TI(i) = (ks * hd(i) * Tsnow(i) + ki * hstemp * TB)/(ki * hsnow(i) + ks * hd(i));
                    
            end
            
        end
    end
    
    %This loop is for ice of thickness greater than 50 cm.
    
    for i = Nh1+1:Nh
        
        hstemp = hsnow(i);
        
        if(hstemp == 0)
            
            % Conductive fluxes in the ice interior
            F1 = 2 * ki * (T2(i) - T1(i))/hd(i);
            F2 = 4 * ki * (TB - T2(i))/hd(i);
            
            Fa = -(1 - a(i)) * (1 - xi) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * TI(i)^4;
            Tdiff = (-Fa * hd(i) + 4 * ki * (T1(i) - TI(i)))/(4 * sigma * TI(i)^3 * hd(i) + 4 * ki);
            Ttemp = TI(i) + Tdiff;
            
            if(Ttemp < 273)
                
                TI(i) = Ttemp;
                F0 = -(1 - a(i)) * (1 - xi) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * TI(i)^4;
                f(i) = (F2 - Fb)/qB;
                
            elseif(Ttemp >= 273)
                
                TI(i) = 273;
                Fa = -(1 - a(i)) * (1 - xi) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * TI(i)^4;
                F0 = 4 * ki * (T1(i) - TI(i))/hd(i);
                f(i) = (Fa - F0)/qI + (F2 - Fb)/qB;
                
            end
            
            % Temperature updates
            
            T1(i) = T1(i) + 2 * delt_dim * (F1 - F0)/(rhocpi * hd(i));
            T2(i) = T2(i) + 2 * delt_dim * (F2 - F1)/(rhocpi * hd(i));
            
        end
         
        % This loop is for snow thickness between 0 and 15 cm.
        
        if((hstemp > 0) && (hstemp <= 0.15))
            
            % Conductive fluxes in the ice and snow interior
            F1 = 2 * ki * (T2(i) - T1(i))/hd(i);
            F2 = 4 * ki * (TB - T2(i))/hd(i);
            
            Fa = -(1 - asnow) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * Tsnow(i)^4;
            beta1 = ks * hd(i)/(4 * ki * hstemp + ks * hd(i)); 
            beta2 = 4 * ki * hstemp/(4 * ki * hstemp + ks * hd(i));
            Tdiff = (-Fa * hstemp + ks * (beta2 * T1(i) + (beta1 - 1) * Tsnow(i)))/(4 * sigma * Tsnow(i)^3 * hstemp + ks * (1 - beta1));
            Ttemp = Tsnow(i) + Tdiff;
            
            if(Ttemp < 273)
                
                Tsnow(i) = Ttemp;
                Fs = -(1 - asnow) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * Tsnow(i)^4;
                
            elseif(Ttemp >= 273)
                
                Tsnow(i) = 273;
                Fa = -(1 - asnow) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * Tsnow(i)^4;
                Fs = ks * ((beta1 - 1) * Tsnow(i) + beta2 * T1(i))/hstemp;
                hsnow(i) = hsnow(i) + delt_dim * (Fa - Fs)/qS;
                
            end
            
            if(hsnow(i) <= 0)
                    
                hsnow(i) = 0;
                    
            end
            
            f(i) = (F2 - Fb)/qB;
            
            % Temperature updates
            
            T1(i) = T1(i) + 2 * delt_dim * (F1 - Fs)/(rhocpi * hd(i));
            T2(i) = T2(i) + 2 * delt_dim * (F2 - F1)/(rhocpi * hd(i));
            
        end
        
        % This loop is for snow thickness above 15 cm.
        
        if(hstemp > 0.15)
            
            % Conductive fluxes in the ice interior
            F0 = 4 * ki * (T1(i) - TI(i))/hd(i);
            F1 = 2 * ki * (T2(i) - T1(i))/hd(i);
            F2 = 4 * ki * (TB - T2(i))/hd(i);
            
            Fa = -(1 - asnow) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * Tsnow(i)^4;
            Tdiff = (-Fa * hstemp + 2 * ks * (T0(i) - Tsnow(i)))/(4 * sigma * Tsnow(i)^3 * hstemp + 2 * ks);
            Ttemp = Tsnow(i) + Tdiff;
            
            if(Ttemp < 273)
                
                Tsnow(i) = Ttemp;
                Fs = -(1 - asnow) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * Tsnow(i)^4;
                
            elseif(Ttemp >= 273)
                
                Tsnow(i) = 273;
                Fa = -(1 - asnow) * Fsi(tp) - Fli(tp) - Fshi(tp) - Flhi(tp) + sigma * Tsnow(i)^4;
                Fs = 2 * ks * (T0(i) - Tsnow(i))/hstemp;
                hsnow(i) = hsnow(i) + delt_dim * (Fa - Fs)/qS;
                
            end
            
            % Temperature updates
            
            T1(i) = T1(i) + 2 * delt_dim * (F1 - F0)/(rhocpi * hd(i));
            T2(i) = T2(i) + 2 * delt_dim * (F2 - F1)/(rhocpi * hd(i));
            
            if(hsnow(i) > 0.15)
                
                T0(i) = T0(i) + delt_dim * (F0 - Fs)/(hsnow(i) * rhocps);
                TI(i) = (2 * ks * hd(i) * T0(i) + 4 * ki * hstemp * T1(i))/(4 * ki * hstemp + 2 * ks * hd(i));
                
            end
            
            f(i) = (F2 - Fb)/qB;
            
        end
    end

    f = f/f0;
    
    % Here the Fokker-Planck equation is solved.
    
    phi01 = k1 - tau * f(1); phi02 = k1 - tau * f(2);
    
    J = f(1);
    
    A(1,1) = 1; A(Nh,Nh) = 1;
    
    if(J > 0)
        
        A(1,2) = 0;
        
    else
        
        A(1,2) = - cgridc * (0.5 * phi02 + k2/delh)/(1 - cgridc * (0.5 * phi01 - k2/delh));
        
    end
    
    for i = 2:Nh-1
        
        phi1 = k1 - tau * f(i+1); phi2 = k1 - tau * f(i-1);
        
        A(i,i-1) = 0.5 * phi2 - k2/delh;
        
        A(i,i) = 2 * k2/delh + 1/cgrid;
        
        A(i,i+1) = - (0.5 * phi1 + k2/delh);
        
    end
    
    rhs = g/cgrid;
    
    if(J > 0)
        
        rhs(1) = 0;
        
    else
        
        rhs(1) = cgrid * rhs(1)/(1 - cgridc * (0.5 * phi01 - k2/delh));
        
    end
    
    g = A\rhs;
    
    if(J > 0)
        
        g = g/trapz(h,g);
        
        Af = 0;
        
    else
        
        Af = 1 - trapz(h,g); 
        
    end

    if(rem(t1,t_plot) == 0)
        
        f_plot(count,:) = f;
        hs_plot(count,:) = hsnow;
        gplot(count,:) = g;
        hmean(count) = dot(h,g)*delh;
        fmean(count) = J;
        a_mean(count) = dot(a,g) * delh;
        Af_plot(count) = Af;
        
    end
    
    count = count + 1;
    
end

tplot = linspace(0,360,Navg1+1);
t1 = linspace(0,T,Nplot+1)/Tyear;

fyear = f_plot(Navg2+1:end,:);
ayear = a_mean(Navg2+1:end);
hsm = hs_plot(Navg2+1:end,:);
gm = gplot(Navg2+1:end,:);
Af_year = Af_plot(Navg2+1:end);

hm = hmean(Navg2+1:end); hm1 = mean(hm);
disp(['Mean thickness =' num2str(hm1)])

fm = fmean(Navg2+1:end);

% Plotting

figure('color','white');
plot(t1, hmean,'linewidth',4);
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',32,'FontName','Times')
xlabel('$t (year)$','FontUnits','points','interpreter','latex','FontSize',32,'FontName','Times')
ylabel('$h (m)$','FontUnits','points','interpreter','latex','FontSize',32,'FontName','Times')

figure('color','white');
plot(tplot,hm,'linewidth',4);
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',32,'FontName','Times')
set(gca,'xtick',15:30:360), datetick('x',4,'keepticks'), set(gca,'xlim',[0,360])
ylabel('$\left<h\right> (m)$','FontUnits','points','interpreter','latex','FontSize',32,'FontName','Times')

figure('color','white');
plot(tplot,fm,'linewidth',4);
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',32,'FontName','Times')
set(gca,'xtick',15:30:360), datetick('x',4,'keepticks'), set(gca,'xlim',[0,360])
ylabel('$\left<f\right>$','FontUnits','points','interpreter','latex','FontSize',32,'FontName','Times')

figure('color','white');
plot(tplot,Af_year,'linewidth',4); 
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',32,'FontName','Times')
set(gca,'xtick',15:30:360), datetick('x',4,'keepticks'), set(gca,'xlim',[0,360])
ylabel('$A(t)$','FontUnits','points','interpreter','latex','FontSize',32,'FontName','Times')











