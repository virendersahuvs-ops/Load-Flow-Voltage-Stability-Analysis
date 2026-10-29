%% Load Flow & Voltage Stability Analysis - IEEE 9-Bus (WSCC) System
% Method : Newton-Raphson (polar form, analytic Jacobian)
% Study  : Base-case load flow, stressed loading, PV (nose) curves,
%          critical point of voltage collapse
% Run    : save as load_flow_voltage_stability.m and press Run
%          (local functions at the end require MATLAB R2016b or later)

clear; clc; close all;

%% 1. System data (100 MVA base, per unit)
Sbase = 100;

% bus: [no  type  Pg  Qg  Pd  Qd  |V|0]   type 1=slack, 2=PV, 3=PQ
bus = [ ...
    1  1  0.00  0  0.00  0.00  1.040
    2  2  1.63  0  0.00  0.00  1.025
    3  2  0.85  0  0.00  0.00  1.025
    4  3  0.00  0  0.00  0.00  1.000
    5  3  0.00  0  1.25  0.50  1.000
    6  3  0.00  0  0.90  0.30  1.000
    7  3  0.00  0  0.00  0.00  1.000
    8  3  0.00  0  1.00  0.35  1.000
    9  3  0.00  0  0.00  0.00  1.000];

% line: [from  to  R  X  B_total]
line = [ ...
    1  4  0.0000  0.0576  0.000
    4  5  0.0170  0.0920  0.158
    5  7  0.0320  0.1610  0.306
    4  6  0.0100  0.0850  0.176
    6  9  0.0390  0.1700  0.358
    7  2  0.0000  0.0625  0.000
    7  8  0.0085  0.0720  0.149
    8  9  0.0119  0.1008  0.209
    9  3  0.0000  0.0586  0.000];

nb   = size(bus, 1);
Ybus = build_ybus(line, nb);

fprintf('=== Y-bus (magnitude, pu) ===\n');
disp(round(abs(Ybus), 3));

%% 2. Base-case Newton-Raphson load flow
[V, ok, iters, J, hist] = nr_loadflow(Ybus, bus, 1.0, []);
assert(ok, 'Base case did not converge.');

fprintf('\n=== Base case: converged in %d iterations ===\n', iters);
fprintf('Max mismatch per iteration:\n');
fprintf('  it %d : %.3e\n', [(0:numel(hist)-1); hist]);

S      = V .* conj(Ybus * V);                 % net injection at each bus
Pgen   = real(S) + bus(:,5);                  % generation = injection + load
Qgen   = imag(S) + bus(:,6);

fprintf('\nBus |  |V| (pu) | Angle (deg) | Pg (MW) | Qg (MVAr) | Pd (MW) | Qd (MVAr)\n');
for k = 1:nb
    if bus(k,2) == 3, pg = 0; qg = 0; else, pg = Pgen(k)*Sbase; qg = Qgen(k)*Sbase; end
    fprintf('%3d | %8.4f | %11.3f | %7.2f | %9.2f | %7.2f | %9.2f\n', ...
        k, abs(V(k)), rad2deg(angle(V(k))), pg, qg, bus(k,5)*Sbase, bus(k,6)*Sbase);
end

% Line flows and losses
fprintf('\nFrom To |  P_ij (MW)  Q_ij (MVAr) |  Loss (MW)\n');
Ploss = 0;
for l = 1:size(line,1)
    f = line(l,1); t = line(l,2);
    z = line(l,3) + 1j*line(l,4); y = 1/z; b = 1j*line(l,5)/2;
    Ift = (V(f)-V(t))*y + V(f)*b;
    Itf = (V(t)-V(f))*y + V(t)*b;
    Sft = V(f)*conj(Ift);
    Stf = V(t)*conj(Itf);
    loss = real(Sft + Stf);
    Ploss = Ploss + loss;
    fprintf('%3d %2d | %9.2f %11.2f | %9.3f\n', f, t, real(Sft)*Sbase, imag(Sft)*Sbase, loss*Sbase);
end
fprintf('Total active loss = %.3f MW\n', Ploss*Sbase);

%% 3. Jacobian structure at the base-case solution
type = bus(:,2); pv = find(type==2); pq = find(type==3); pvpq = [pv; pq];
fprintf('\nJacobian size at solution: %d x %d  (J11: dP/dTheta, J12: dP/d|V|, J21: dQ/dTheta, J22: dQ/d|V|)\n', size(J,1), size(J,2));
fprintf('Base-case min singular value of J = %.4f\n', min(svd(J)));

%% 4. Stressed loading -> PV (nose) curves
% All loads scaled by lambda at constant power factor; PV-bus generation
% is held fixed and the slack bus picks up the extra load and losses.
step     = 0.05;
minStep  = 1e-5;
lam      = 0;
Vguess   = [];
lamVec   = [];
VmMat    = [];          % rows = lambda points, cols = buses
sigMin   = [];

while step > minStep
    [Vt, okt, ~, Jt] = nr_loadflow(Ybus, bus, lam + step, Vguess);
    if okt
        lam    = lam + step;
        Vguess = Vt;
        lamVec(end+1,1) = lam;               %#ok<SAGROW>
        VmMat(end+1,:)  = abs(Vt).';         %#ok<SAGROW>
        sigMin(end+1,1) = min(svd(Jt));      %#ok<SAGROW>
    else
        step = step/2;                       % refine towards the nose
    end
end
lamCrit = lam;

% prepend lambda = 0 (no extra load)
[V0, ~, ~, J0] = nr_loadflow(Ybus, bus, 0, []);
lamVec = [0; lamVec];
VmMat  = [abs(V0).'; VmMat];
sigMin = [min(svd(J0)); sigMin];

Pbase  = sum(bus(:,5))*Sbase;
fprintf('\n=== Voltage stability ===\n');
fprintf('Base load                 : %.1f MW\n', Pbase);
fprintf('Critical loading factor   : %.4f x base\n', lamCrit);
fprintf('Maximum loadability       : %.1f MW\n', lamCrit*Pbase);
fprintf('Load margin to collapse   : %.1f MW (%.1f %%)\n', (lamCrit-1)*Pbase, (lamCrit-1)*100);

[~, weak] = min(VmMat(end, pq));              % weakest PQ bus at the nose
weakBus   = pq(weak);
fprintf('Weakest bus at nose point : Bus %d (|V| = %.4f pu)\n', weakBus, VmMat(end, weakBus));

%% 5. Plots
figure('Name','PV (Nose) Curves','Color','w');
plot(lamVec*Pbase, VmMat(:, [5 6 8 4 7 9]), 'LineWidth', 1.5); hold on; grid on;
xline(Pbase, '--k', 'Base load');
xline(lamCrit*Pbase, '--r', 'Collapse point');
xlabel('Total system load (MW)'); ylabel('Bus voltage magnitude (pu)');
title('PV Curves - IEEE 9-Bus System (Newton-Raphson)');
legend('Bus 5','Bus 6','Bus 8','Bus 4','Bus 7','Bus 9','Location','southwest');

figure('Name','Jacobian Singular Value','Color','w');
plot(lamVec*Pbase, sigMin, 'LineWidth', 1.8); grid on; hold on;
xline(lamCrit*Pbase, '--r', 'Collapse point');
xlabel('Total system load (MW)'); ylabel('\sigma_{min}(J)');
title('Minimum singular value of Jacobian -> 0 at voltage collapse');

%% 6. Save figures to ../results
outDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
saveas(findobj('Type','figure','Name','PV (Nose) Curves'),        fullfile(outDir,'pv_curves_matlab.png'));
saveas(findobj('Type','figure','Name','Jacobian Singular Value'), fullfile(outDir,'jacobian_sigma_matlab.png'));

%% ===================== Local functions =====================

function Ybus = build_ybus(line, nb)
    Ybus = zeros(nb);
    for l = 1:size(line,1)
        f = line(l,1); t = line(l,2);
        y  = 1/(line(l,3) + 1j*line(l,4));
        bc = 1j*line(l,5)/2;                 % half line charging at each end
        Ybus(f,f) = Ybus(f,f) + y + bc;
        Ybus(t,t) = Ybus(t,t) + y + bc;
        Ybus(f,t) = Ybus(f,t) - y;
        Ybus(t,f) = Ybus(t,f) - y;
    end
end

function [V, ok, it, J, hist] = nr_loadflow(Y, bus, lam, Vinit)
    tol = 1e-8; maxit = 30;
    nb   = size(bus,1);
    type = bus(:,2);
    pv   = find(type==2); pq = find(type==3); pvpq = [pv; pq];
    npvpq = numel(pvpq);

    Psch = bus(:,3) - lam*bus(:,5);
    Qsch =          - lam*bus(:,6);

    Vm = bus(:,7); Va = zeros(nb,1);
    if ~isempty(Vinit)
        Vm = abs(Vinit); Va = angle(Vinit);
        Vm(type~=3) = bus(type~=3,7);        % slack/PV magnitudes are fixed
    end

    ok = false; hist = []; J = [];
    for it = 0:maxit
        V = Vm .* exp(1j*Va);
        S = V .* conj(Y*V);
        mis = [Psch(pvpq) - real(S(pvpq)); Qsch(pq) - imag(S(pq))];
        hist(end+1) = max(abs(mis));         %#ok<AGROW>

        if hist(end) < tol, ok = true; break; end
        if any(isnan(mis)) || any(Vm(pq) < 0.3) || hist(end) > 1e3, break; end

        % Analytic Jacobian [J11 J12; J21 J22]
        Ibus = Y*V;
        dS_dVa = 1j*diag(V) * conj(diag(Ibus) - Y*diag(V));
        dS_dVm = diag(V) * conj(Y*diag(V./abs(V))) + conj(diag(Ibus))*diag(V./abs(V));
        J11 = real(dS_dVa(pvpq, pvpq));  J12 = real(dS_dVm(pvpq, pq));
        J21 = imag(dS_dVa(pq,   pvpq));  J22 = imag(dS_dVm(pq,   pq));
        J   = [J11 J12; J21 J22];

        dx = J \ mis;
        Va(pvpq) = Va(pvpq) + dx(1:npvpq);
        Vm(pq)   = Vm(pq)   + dx(npvpq+1:end);
    end
    if ok
        % Jacobian at the converged point (used for stability index)
        Ibus = Y*V;
        dS_dVa = 1j*diag(V) * conj(diag(Ibus) - Y*diag(V));
        dS_dVm = diag(V) * conj(Y*diag(V./abs(V))) + conj(diag(Ibus))*diag(V./abs(V));
        J = [real(dS_dVa(pvpq,pvpq)) real(dS_dVm(pvpq,pq)); ...
             imag(dS_dVa(pq,pvpq))   imag(dS_dVm(pq,pq))];
    end
end
