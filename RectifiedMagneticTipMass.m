%% 
clear;
%% ============================ USER INPUT =========================
Mag_force = input('Magnetic force included? (1 = Yes): ');
calc_dis = 1;
min_freq = input('Minimum frequency to be assessed: ');
max_freq = input('Maximum frequency to be assessed: ');
int_freq = input('Frequency interval (1 if single value): ');

Frequency_range = min_freq:int_freq:max_freq;
w_range = Frequency_range * 2 * pi;

Resistance = input('Resistances in ohms to be assessed (vector): ');
C = input('Capacitance in nF to be assessed: ')*10^-9;

clear min_freq max_freq;

%% ========================== BEAM & MATERIAL PROPERTIES ===========
% Beam Dimensions (units: meters) %Checked
L = 60e-3;
hp = 0.267e-3;
hs = 0.3e-3;
b = 25e-3;

% Material Properties %Checked (NOTE: Tip mass still present as original)
Ys = 72e9; % Young's modulus (substrate)
Yp = 62e9; % Young's modulus (piezoelectric)
rhop = 7800; % Density (piezoelectric)
rhos = 2700; % Density (substrate)
d31 = -320e-12; % Piezoelectric constant
ep33t = 1800 * 8.854187817e-12; % Permittivity at constant stress %%
%Relative dielectric constant should be 1800 not 3800 - mistake from
%previous student
ep33s = ep33t - (d31^2) * Yp; % Permittivity at constant strain

%%%%CORRECTION (18/2/26)^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%ep33s = 1.3555e-08;%1.83369e-8; %Updated 14/02/26
%%NB: there was a misprint in Sensors paper; This value agrees with Dalzell Bonello 2012 paper in 
%%Smart Mater. Struct. 21 (2012) 10502 (page 5)
%%You can also uncheck lines 31,32 and find that the expression in line 32 yields the same value (1.3555e-08)
%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

xi = 1.09e-2; % Damping ratio for identical beam with no tip mass
Omega = -(d31 * Yp * b * (hp + hs)) / 2;
m = b * (2 * hp * rhop + hs * rhos) % Distributed mass
mt = 4.3e-3; % Tip magnet mass

n = 1.984; % Emission factor
k = 1.38e-23; % Boltzmann const
T = 300; % Reference temp
Q = 1.6e-19; % Charge of an electron

VT = (k*T)/Q; % Thermal voltage
IS = 1.411e-8; % Saturation current

% Derived Beam Constants %Checked
YI = b * (((Ys * hs^3) / 12) + Yp * (((hp * hs^2) / 2) + (hp^2 * hs) + (2/3) * hp^3)); % flexuaral stiffness of entire beam
Cp = (ep33s * b * L) / hp; %Capacitance
zc = (hp + hs) / 2; % midpoint depth 


%% ========================= MODAL ANALYSIS ======================== %Checked 14/02/26
lambdar0= 1.87510407; %1.87510407 is lamdar for a clamped beam with no tip mass
lambdar = fzero(@findlamdar, lambdar0, [], mt, m*L); % With tip mass

% 23/02/26 sigmar now for tip mass included
sigmar = (sinh(lambdar) - sin(lambdar) + lambdar * (mt / (m *L)) *(cosh(lambdar) - cos(lambdar))) /(cosh(lambdar) + cos(lambdar) + lambdar * (mt / (m * L)) *(sinh(lambdar) - sin(lambdar)));
%sigmar = (cosh(lambdar) + cos(lambdar))/(sinh(lambdar) + sin(lambdar)); %Updated mode shape coefficient w/o tip mass check if other form makes a difference 15/02/26


% Symbolic Calculation for Mode Shape and Integrals
syms x;
psir_x = cosh((lambdar/L)*x) - cos((lambdar/L)*x)- sigmar * (sinh((lambdar/L)*x) - sin((lambdar/L)*x));
psir_x_sqr = psir_x^2;
int_psi_r_sqr = double(int(psir_x_sqr, x, 0, L));

psir_x_sqr_L = (cosh((lambdar/L)*L) - cos((lambdar/L)*L)- sigmar * (sinh((lambdar/L)*L) - sin((lambdar/L)*L)))^2;
Mr=m*int_psi_r_sqr + mt*psir_x_sqr_L; % 23/02/26 updated with tip mass 
phir_x = sqrt(1 / Mr) * psir_x; % Mass-normalized mode shape
phir_L = sqrt(1 / Mr) * (cosh((lambdar/L)*L) - cos((lambdar/L)*L)- sigmar * (sinh((lambdar/L)*L) - sin((lambdar/L)*L)));
int_phi_r = double(int(phir_x, x, 0, L));
diff_phi_r_L = double(subs(diff(phir_x, x), x, L));


%% ============================ MODEL CONSTANTS =================== %Checked 14/0/26
wr = lambdar^2 * sqrt(YI / (m * L^4)); % Natural frequency (undamped)
wr0= lambdar0^2 * sqrt(YI / (m * L^4)); % Natural frequency (undamped) with no mass
xi_mass= wr*(xi/wr0); % damping ratio with mass, wr = wr0 for this anyway 23/02/26 updated for with mass
%xi_mass = xi; %New damping ratio w/o mass
Xr = Omega * diff_phi_r_L; % Piezoelectric coupling term
alph_r = -Yp * d31 * b * zc * diff_phi_r_L;
% Constants for equation of motion
A1 = -2 * xi_mass * wr;
A2 = -wr^2;
A3 = -m * int_phi_r -mt * phir_L; % 23/02/26 updated for with mass
A4 = phir_L;

B1 = 2 / Cp;
B2 = IS;
B3 = n*VT;

%% ========== QZS Optimization via Magnetic Force Tuning ===========
if Mag_force == 1
 % Magnetic Properties
 Sa = 200*1e-6; % Surface area (m^2)
 Sb = Sa;
 Ma = 0.96*1e6; % Magnetization (A/m)
 Mb = Ma;
 u0 = 4 * pi * 1e-7;
 a = (2.75 / 2) * 1e-3;

 wa = linspace(-0.02, 0.02, 200); % Deflection range
 best_d = NaN; % Distance at which QZS is achieved is initialized
 min_derivative = Inf;

 for dk = 5:0.05:15
     d = dk * 1e-3;
     Fmagy = zeros(size(wa)); % Initializing the magnetic tip force
     Fbeam = zeros(size(wa)); % Initializing the restoring elastic tip force
     Fmag_1 = zeros(size(wa));
     Fmag_2 = zeros(size(wa));

     for c = 1:length(wa)
         w = wa(c);
         du = diff_phi_r_L * w / phir_L;

         Fmag_1 = (Ma * Sa * (-Mb * Sb * w / ((2*a + d)^2 + w^2)^(3/2) + Mb * Sb * w / ((4*a + d)^2 + w^2)^(3/2)));

         Fmag_2 = (Ma * Sa * (-Mb * Sb * (w + 2*a*du) / (d^2 + (w + 2*a*du)^2)^(3/2) + Mb * Sb * (w + 2*a*du) / ((2*a + d)^2 + (w +2*a*du)^2)^(3/2)));

         Fmagy(c) = (u0 / (4*pi)) * (Fmag_1 - Fmag_2);

         Fbeam(c) = 3 * (YI / L^3) * w;

      end

      Ftotal = -Fmagy + Fbeam;
      % Potential energy: Numerical method to find the integral of the magnetic force over the deflection
      U_mag = -cumtrapz(wa, Fmagy);
      dw = wa(2) - wa(1);
      derivative_values = gradient(Ftotal, dw); % Gradient for a step dw

      [~, idx_zero] = min(abs(wa-0)); % finds index value at wa=0
      derivative_at_zero = derivative_values(idx_zero); % finds derivative at wa = 0


 % if function finds the distance with smallest derivative throught the loop

        if abs(derivative_at_zero) < abs(min_derivative)
            min_derivative = derivative_at_zero;
            best_d = d;
            best_Fmagy = Fmagy;
            best_U_mag = U_mag;
        end
 end

 % Plotting Restoring force and static potential energy
 U_beam = ((wa*(1/phir_L)).^2)*0.5 * YI* double(int((diff(diff(phir_x,x),x))^2,x,0,L)); %static mechanical potential energy
 U_total = U_beam + best_U_mag;
 Ftotal= -best_Fmagy+Fbeam;

 fprintf('Optimal magnet gap (d): %.4f m\n', best_d);
 fprintf('Derivative at wa=0: %.4e\n', min_derivative);

 figure;
 plot(wa, Fbeam, 'k','LineWidth', 0.5)
 hold on
 plot(wa, best_Fmagy, 'k--', wa, Ftotal, 'k', 'LineWidth', 1.5);
 hold off
 xlabel('Tip Deflection (m)', 'FontWeight', 'bold','FontSize', 15, 'FontName','Cambria');
 ylabel('F_{restoring} (N)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName','Cambria');
 legend('Beam', 'Magnet', 'Total','FontWeight', 'bold', 'FontSize', 13,'FontName', 'Cambria')
 set(gca, 'FontSize', 12, 'FontName', 'Cambria'); % axes ticks font

 figure;
 plot(wa, U_beam, 'k','LineWidth', 0.5)
 hold on
 plot(wa, best_U_mag, 'k--', wa, U_total, 'k', 'LineWidth', 1.5);
 hold off
 legend('Beam', 'Magnet', 'Total', 'FontWeight', 'bold','FontSize', 13,'FontName', 'Cambria');
 xlabel('Tip Deflection (m)', 'FontWeight', 'bold','FontSize', 15, 'FontName','Cambria');
 ylabel('Potential Energy (N)', 'FontWeight', 'bold','FontSize', 15,'FontName', 'Cambria');
 set(gca, 'FontSize', 12, 'FontName', 'Cambria'); % axes ticks font
else
        best_d=0;
end

%% ====================== RUNNING MODEL ===========================
volt = zeros(length(Frequency_range), length(Resistance));
voltcap = zeros(length(Frequency_range), length(Resistance));
voltDC = zeros(length(Frequency_range), length(Resistance));
voltcapAC = zeros(length(Frequency_range), length(Resistance));
if calc_dis == 1
 tip_dis = zeros(length(Frequency_range), length(Resistance));
end

% finds the respective positions of a chosen frequency for deriving time histories
[~, idx_zero] = min(abs(w_range - (4*2*pi)));
[~, idx_zero2] = min(abs(w_range - (70*2*pi)));

w_range_back = flip(w_range);      % Reverse frequency vector
nFreq = length(w_range_back);

for a = 1:length(Resistance)
    y0=[0,0,0,0];
    R = Resistance(a);
    for c = 1:nFreq
        w = w_range_back(c);
        %tf=t0+ncycles*2*pi/w;
        [t, x] = ode15s(@(t, q) piezoelecbeamresistor(t, q, A1, A2, A3, B1, B2, B3, R, C, w, alph_r, Xr, A4, diff_phi_r_L, Mag_force, best_d), [0 4], y0);
        y0=x(end,:); 
        
        %d = sum(t <= tf-5*2*pi/w);1
        d = sum(t <= 3.8);

        volt(c, a) = (max(x(d:end, 3)) - min(x(d:end, 3))) / 2;   % Dominant AC component of Vpiezo
     
        voltcap(c, a) = mean(x(d:end,4)); % Dominant DC component of Vcap This was taking the DC component of PIEZO
        %voltcap(c,a) = (1/4) * sum(x(1:end-1, 3).*(t(2:end)-t(1:end-1)));
     
        voltDC(c, a) = mean(x(d:end, 3)); % DC component of Vpiezo
        %voltDC(c, a) = (1/4) * sum(x(1:end-1, 3).*(t(2:end)-t(1:end-1)));

        voltcapAC(c, a) = (max(x(d:end, 4)) - min(x(d:end, 4))) / 2;  % AC component of Vcap Returns faulty for some reason

        if calc_dis == 1
             tip_dis(c, a) = (max(x(d:end, 1)) - min(x(d:end, 1))) / 2; 
             f = w / (2 * pi);
             disp(f)
        end

        if ismember(c, [idx_zero, idx_zero2])
            xdot = x(:, 2) * phir_L;
            x_dis= x(:, 1) * phir_L;
            figure;
            plot(t, x_dis*1e3);
            xlabel('Time (s)','FontWeight', 'bold','FontSize', 15, 'FontName','Cambria'); 
            ylabel('Tip Displacement (mm)','FontWeight', 'bold','FontSize', 15,'FontName', 'Cambria');
            set(gca, 'FontWeight', 'bold','FontSize', 12, 'FontName', 'Cambria')
            figure;
            plot(x_dis*1e3, xdot);
            xlabel('Tip Displacement (mm)','FontWeight', 'bold','FontSize', 15,'FontName', 'Cambria');
            ylabel('Velocity (m/s)','FontWeight', 'bold','FontSize', 15,'FontName', 'Cambria');
            set(gca, 'FontWeight', 'bold','FontSize', 12, 'FontName', 'Cambria')
        

            eta=x(:,1);    
            etadot=x(:,2);
            vol=x(:,3);
            vcap = x(:,4);

            % figure('Color',w)
            % plot(t, vol, 'b', 'LineWidth', 1.5)
            % hold on
            % plot(t, vcap, 'r', 'LineWidth', 1.5)
            % hold off
            % xlabel('Time (s)')
            % ylabel('Voltage (V)')
            % legend('V_p','V_{cap}')
            % grid on
            figure('Color','w');   % White background
hold on

% Plot signals
h1 = plot(t, vol, ...
    'Color', [0 0.447 0.7410], ...
    'LineWidth', 1.5);

h2 = plot(t, vcap, ...
    'Color', [0.8 0 0], ...
    'LineWidth', 1);

hold off
grid on
box on

% --- Title and Labels (bigger, clean) ---
title('Unrectified Piezoelectric Voltage and Rectified External Load Voltage vs Time', ...
    'FontSize', 28, ...
    'FontWeight', 'bold', ...
    'FontName', 'Arial');

xlabel('Time (s)', ...
    'FontSize', 18, ...
    'FontWeight', 'bold', ...
    'FontName', 'Arial');

ylabel('Voltage (V)', ...
    'FontSize', 18, ...
    'FontWeight', 'bold', ...
    'FontName', 'Arial');

% Legend
legend([h1 h2], {'V_{piezo}','V_{cap}'}, ...
    'Location','best', ...
    'FontSize', 20, ...
    'FontWeight','bold',...
    'Box','on');

% --- Axes styling ---
ax = gca;
ax.FontSize = 18;     
ax.FontWeight = 'bold';
ax.LineWidth = 1.2;
ax.TickDir = 'in';      
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

% --- X axis ---
xlim([0 4])
xticks(0:1:4)

% --- PROFESSIONAL AUTO Y-SCALING ---
yall = [vol(:); vcap(:)];
yall = yall(isfinite(yall));

if ~isempty(yall)
    ymin = min(yall);
    ymax = max(yall);

    if ymin == ymax
        ymin = ymin - 1;
        ymax = ymax + 1;
    end

    pad = 0.05 * (ymax - ymin);
    ylim([ymin - pad, ymax + pad]);
end

yticks('auto')

            wa= eta*phir_L;
            d=best_d;
            du=diff_phi_r_L*eta;
            Fmag_y = zeros(size(wa)); % Initializing the magnetic tip force
            if Mag_force == 1%%%%%added 27 Nov 2025
              for k = 1:length(wa)
                   w = wa(k);

                   Fmag1 = (Ma * Sa * (-Mb * Sb * w / ((2*a + d)^2 + w^2)^(3/2) + Mb * Sb * w / ((4*a + d)^2 + w^2)^(3/2)));

                   Fmag2 = (Ma * Sa * (-Mb * Sb * (w + 2*a*du(k)) / (d^2 + (w + 2*a*du(k))^2)^(3/2)+ Mb * Sb * (w + 2*a*du(k)) / ((2*a + d)^2 + (w +2*a*du(k))^2)^(3/2)));

                   Fmag_y(k) = (u0 / (4*pi)) * (Fmag1 - Fmag2);

              end
            end%if Mag_force == 1%%%%%added 27 Nov 2025
            Umag = -cumtrapz(wa, Fmag_y);
            vdot=gradient(x(:,3), t);        %vdot for piezo voltage not for vcap
            syms x
            int_2diff_phir_L=double(int((diff(diff(phir_x,x),x))^2,x,0,L));
            int_diff_phir_L=double(int(diff(phir_x,x),x,0,L));
            Ubeam = (eta.^2)*0.5 * YI*int_2diff_phir_L+(eta.*(alph_r*etadot-((Cp/2)*gradient(vol, t))))*2*Omega*R*int_diff_phir_L - ((Cp/2)*(R*(alph_r*etadot-((Cp/2)*vdot))).^2);
            Utotal= Umag+Ubeam;
            figure
            plot(wa*1e3, Utotal, 'b','LineWidth', 1.5)
            hold on
            plot(wa*1e3, Umag, 'r--','LineWidth', 1.5);
            plot(wa*1e3, Ubeam, 'g--', 'LineWidth', 1.5);
            hold off
            legend('Total', 'Magnet', 'Beam', 'FontWeight', 'bold','FontSize', 13,'FontName', 'Cambria');
            xlabel('Tip Deflection (mm)', 'FontWeight', 'bold','FontSize', 15,'FontName', 'Cambria');
            ylabel('Potential Energy-dynamic (N)', 'FontWeight','bold','FontSize', 15, 'FontName', 'Cambria');
            set(gca, 'FontSize', 12, 'FontName', 'Cambria'); % axes ticks font
        end
    end
    disp(['R = ', num2str(R), ' Ohms']);
end

if calc_dis == 1
    tip_dis = tip_dis * phir_L;
end

%% ========================= PLOTTING RESULTS ======================
% % Voltage FRF
% if input('Plot Vpiezo AC component FRF (1 = Yes): ') == 1
%     figure;
%     for a = 1:length(Resistance)
%         semilogy(Frequency_range, volt(:, a),'-');
%         hold on;
%     end
%     xlabel('Frequency (Hz)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName','Cambria');
%     ylabel('Voltage (V)', 'FontWeight', 'bold','FontSize', 15, 'FontName','Cambria');
%     legend(arrayfun(@(R) sprintf('R = %.2f \\Omega', R), Resistance,'UniformOutput', false), 'FontWeight', 'bold','FontSize', 13, 'FontName','Cambria');
%     set(gca, 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Cambria'); % axes ticks font
% end

if input('Plot Vpiezo AC component FRF (1 = Yes): ') == 1
    figure('Color','w');
    hold on

    for a = 1:length(Resistance)
        semilogy(Frequency_range, volt(:, a), '-', 'LineWidth', 3);
    end

    grid on
    box on

    title('Piezoelectric Voltage Frequency Response Plot for a Rectified Harvester', ...
        'FontSize', 28, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial');

    xlabel('Frequency (Hz)', ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'FontName','Arial');

    ylabel('Voltage (V)', ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'FontName','Arial');

    legend(arrayfun(@(R) sprintf('R = %.2f \\Omega', R), Resistance, ...
        'UniformOutput', false), ...
        'FontSize', 18, ...
        'Box','on');

    % --- Axes styling (MATCHED) ---
    ax = gca;
    ax.FontSize = 18;
    ax.FontWeight = 'bold';
    ax.LineWidth = 1.2;
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';

    % --- Auto Y scaling (log-safe) ---
    yall = volt(:);
    yall = yall(isfinite(yall) & yall > 0);

    if ~isempty(yall)
        ymin = min(yall);
        ymax = max(yall);

        if ymin == ymax
            ymin = ymin * 0.9;
            ymax = ymax * 1.1;
        end

        ylim([ymin ymax]);
    end

    yticks('auto')
end

% % Voltage FRF
% if input('Plot Vcap DC component FRF (1 = Yes): ') == 1
%     figure;
%     for a = 1:length(Resistance)
%         semilogy(Frequency_range, voltcap(:, a),'-');
%         hold on;
%     end
%     xlabel('Frequency (Hz)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName','Cambria');
%     ylabel('Voltage (V)', 'FontWeight', 'bold','FontSize', 15, 'FontName','Cambria');
%     legend(arrayfun(@(R) sprintf('R = %.2f \\Omega', R), Resistance,'UniformOutput', false), 'FontWeight', 'bold','FontSize', 13, 'FontName','Cambria');
%     set(gca, 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Cambria'); % axes ticks font
% end

if input('Plot Vcap DC component FRF (1 = Yes): ') == 1
    figure('Color','w');
    hold on

    for a = 1:length(Resistance)
        semilogy(Frequency_range, voltcap(:, a), '-', 'LineWidth', 2);
    end

    grid on
    box on

    title('External Load Voltage Frequency Response Plot for a Rectified Harvester ', ...
        'FontSize', 28, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial');

    xlabel('Frequency (Hz)', ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'FontName','Arial');

    ylabel('Voltage (V)', ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'FontName','Arial');

    legend(arrayfun(@(R) sprintf('R = %.2f \\Omega', R), Resistance, ...
        'UniformOutput', false), ...
        'FontSize', 18, ...
        'Box','on');

    % --- Axes styling (MATCHED) ---
    ax = gca;
    ax.FontSize = 18;
    ax.FontWeight = 'bold';
    ax.LineWidth = 1.2;
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';

    % --- Auto Y scaling (log-safe) ---
    yall = voltcap(:);
    yall = yall(isfinite(yall) & yall > 0);

    if ~isempty(yall)
        ymin = min(yall);
        ymax = max(yall);

        if ymin == ymax
            ymin = ymin * 0.9;
            ymax = ymax * 1.1;
        end

        ylim([ymin ymax]);
    end

    yticks('auto')
end

if input('Plot FRF comparison (1 = Yes): ') == 1
    figure('Color','w');
    hold on

    % --- Plot Piezo ---
    for a = 1:length(Resistance)
        h1 = semilogy(Frequency_range, volt(:, a), '-', ...
            'Color', [0 0.447 0.7410], ...
            'LineWidth', 2);
    end

    % --- Plot Capacitor ---
    for a = 1:length(Resistance)
        h2 = semilogy(Frequency_range, voltcap(:, a), '-', ...
            'Color', [0.8 0 0], ...
            'LineWidth', 2);
    end

    grid on
    box on

    title('Voltage Frequency Response Plot of a Rectified Harvester', ...
        'FontSize', 28, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial');

    xlabel('Frequency (Hz)', ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'FontName','Arial');

    ylabel('Voltage (V)', ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'FontName','Arial');

    % --- Clean Legend (NO AUTO UPDATE) ---
    lgd = legend([h1 h2], {'V_{piezo}', 'V_{cap}'}, ...
        'FontSize', 18, ...
        'Box', 'on');
    lgd.AutoUpdate = 'off';

    % --- Axes styling ---
    ax = gca;
    ax.FontSize = 18;
    ax.FontWeight = 'bold';
    ax.LineWidth = 1.21;
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';

    % --- Auto Y scaling (log safe) ---
    yall = [volt(:); voltcap(:)];
    yall = yall(isfinite(yall) & yall > 0);

    if ~isempty(yall)
        ylim([min(yall) max(yall)]);
    end

    drawnow   % lock axis before using limits

    %% ===== Resonance Markers =====

    xl = xlim;
    yl = ylim;

    x_offset = xl(1) + 0.12*(xl(2)-xl(1));  % shift labels right

    % -------- PIEZO --------
    [Vp_max, idx_p] = max(volt(:));
    [row_p, ~] = ind2sub(size(volt), idx_p);
    f_res_p = Frequency_range(row_p);

    line([f_res_p f_res_p], [yl(1) Vp_max], ...
        'Color', [0 0 0], ...
        'LineStyle', '--', ...
        'LineWidth', 1.8, ...
        'HandleVisibility','off');

    line([xl(1) f_res_p], [Vp_max Vp_max], ...
        'Color', [0 0 0], ...
        'LineStyle', '--', ...
        'LineWidth', 1.8, ...
        'HandleVisibility','off');

    text(x_offset, Vp_max, sprintf('%.3f V', Vp_max), ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'Color', 'k', ...
        'VerticalAlignment', 'bottom');


    % -------- CAPACITOR --------
    [Vc_max, idx_c] = max(voltcap(:));
    [row_c, ~] = ind2sub(size(voltcap), idx_c);
    f_res_c = Frequency_range(row_c);

    line([f_res_c f_res_c], [yl(1) Vc_max], ...
        'Color', [0 0 0], ...
        'LineStyle', '--', ...
        'LineWidth', 1.8, ...
        'HandleVisibility','off');

    line([xl(1) f_res_c], [Vc_max Vc_max], ...
        'Color', [0 0 0], ...
        'LineStyle', '--', ...
        'LineWidth', 1.8, ...
        'HandleVisibility','off');

    text(x_offset, Vc_max, sprintf('%.3f V', Vc_max), ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'Color', 'k', ...
        'VerticalAlignment', 'bottom');


    % -------- Single Resonant Frequency Label --------
    if Vp_max >= Vc_max
        f_label = f_res_p;
    else
        f_label = f_res_c;
    end

end

% Voltage FRF
if input('Plot Vpiezo DC component FRF (1 = Yes): ') == 1
    figure;
    for a = 1:length(Resistance)
        plot(Frequency_range, voltDC(:, a),'-');
        hold on;
    end
    xlabel('Frequency (Hz)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName','Cambria');
    ylabel('Voltage (V)', 'FontWeight', 'bold','FontSize', 15, 'FontName','Cambria');
    legend(arrayfun(@(R) sprintf('R = %.2f \\Omega', R), Resistance,'UniformOutput', false), 'FontWeight', 'bold','FontSize', 13, 'FontName','Cambria');
    set(gca, 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Cambria'); % axes ticks font
end

if input('Plot Vcap AC component FRF (1 = Yes): ') == 1
    figure;
    for a = 1:length(Resistance)
        semilogy(Frequency_range, voltcapAC(:, a),'-');
        hold on;
    end
    xlabel('Frequency (Hz)', 'FontWeight', 'bold', 'FontSize', 15, 'FontName','Cambria');
    ylabel('Voltage (V)', 'FontWeight', 'bold','FontSize', 15, 'FontName','Cambria');
    legend(arrayfun(@(R) sprintf('R = %.2f \\Omega', R), Resistance,'UniformOutput', false), 'FontWeight', 'bold','FontSize', 13, 'FontName','Cambria');
    set(gca, 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Cambria'); % axes ticks font
end


% Tip Displacement FRF1
if input('Plot Tip Displacement FRF (1 = Yes): ') == 1
    figure;
    for a = 1:length(Resistance)
        plot(Frequency_range, tip_dis(:, a));
        hold on;
    end
    xlabel('Frequency (Hz)', 'FontWeight', 'bold','FontSize', 15, 'FontName','Cambria');
    ylabel('Tip Displacement (m)', 'FontWeight', 'bold','FontSize', 15,'FontName', 'Cambria');
    legend(arrayfun(@(R) sprintf('R = %.2f \\Omega', R), Resistance,'UniformOutput', false),'FontWeight', 'bold','FontSize', 12, 'FontName','Cambria');
    set(gca, 'FontSize', 12, 'FontName', 'Cambria'); % axes ticks font
end

%PIEZOELECTRICRESISTOR

%defining the three 1st odes
function qprime = piezoelecbeamresistor (t,q,A1,A2,A3,B1,B2,B3,R,C,w,alph_r,Xr,phir_L,diff_phi_r_L,Mag_force,best_d)
% Modelling circuit with resistor
%q(1)=eta;
%q(2)=etadot;
%q(3)=vp;
%q(4)=vcap
%qprime(1)=etadot;
%qprime(2)=etadotdot;
%qprime(3)=vdot;
%qprime(4)=vcapdot
if Mag_force==1

    Sa=200*1e-6; %(m^2)SA and SB are the magnet surface area
    Sb= Sa;
    Ma=0.96*1e6; %(A/m) MA and MB are the magnetization of magnets A and B
    Mb=Ma;
    u0=4*pi*1e-7; % (Tm/A) free space permittvity
    a=(2.75/2)*1e-3; %length of magnet(m)
    wa= q(1)*phir_L;
    d=best_d;
    du=diff_phi_r_L*q(1);% wa and du is the tip displacement od tip slope

    Fmag_1 = (Ma * Sa * (-Mb * Sb * wa/((2*a + d)^2 + wa^2)^(3/2) + Mb * Sb * wa / ((4*a + d)^2 + wa^2)^(3/2)));

    Fmag_2 = Ma * Sa * (-Mb * Sb * (wa + (2*a*du))/(d^2 + (wa + (2*a*du))^2)^(3/2) + Mb * Sb * (wa + 2*a*du) / ((2*a + d)^2 + (wa +2*a*du)^2)^(3/2));

    Fmagy = (u0 / (4*pi)) * (Fmag_1 - Fmag_2);
else
    Fmagy=0;
end
qprime(1)=q(2); %q=etadot
ub2dot=2*9.81*sin(w*t); % Currently 1g but can be changed to different accelerations like 0.1
qprime(2)=A1*q(2)+A2*q(1)-Xr*q(3)+(ub2dot*A3)+(phir_L*Fmagy);

%Determining ID1 & ID2:

VD1 = (q(3)-q(4))/2;
VD2 = (-q(3)-q(4))/2;

ID1=B2*(exp(VD1/B3)-1);
ID2=B2*(exp(VD2/B3)-1);

%qprime(3)=B1*(alph_r*q(2)-(q(3)/R));
qprime(3)=-B1*(ID1-ID2)+B1*Xr*q(2);
qprime(4)=(ID1+ID2-(q(4)/R))/C;        %Whereby R and C are resistor & capacitor connected externally to rectifier and inputted at the start by user
qprime=qprime(:);
end

%FINDLAMBDAR
function f = findlamdar(lamdar,mt,mb)
% equationn to find the roots
f=1 + cos(lamdar)*cosh(lamdar) + (mt/mb)*lamdar*(cos(lamdar)*sinh(lamdar)-sin(lamdar)*cosh(lamdar)); % Changed to frequency equation with tip mass 23/02/26
%f = 1 + cos(lamdar)*cosh(lamdar);  %freq. eq. for cantilever w/o tip mass
end