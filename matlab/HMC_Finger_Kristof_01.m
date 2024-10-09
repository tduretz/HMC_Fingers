% 2D hydro-mechanical-chemical model, atg --> fo + en + water (de)hydration reaction at ca. 3 GPa
% isothermal, compressible, porous medium, porosity-dependent viscosity
% initial porosity anomaly defines weak inclusion in center
% dehydration reaction is triggered by rising fluid pressure at the edges of the weak inclusion, due to kinematic BC
% by Kristóf Porkoláb, Evangelos Moulas, Stefan Schmalholz, 07.2024
% Modified to work for Erwan's look-up table for eclogitization

clear all, close all, clc
% ----------------------- Switches ------------------------------ %
SAVE   = 0;                             % 1 for saving model parameters, folder specified below
YES_vp = 0;                             % 1 for using Von Mises plasticity
% ------------------------- Physics ------------------------------------ %
% Physical Parameters
phi_amb                = 0.01;
% dphi0                  = 0.025;
L_phys                 = 0.1;             % Physical Length                     (m)
P_phys                 = 19e8;            % Reaction Pressure                   (Pa)
k_phys                 = 1e-22;           % permeability                        (m2)
eta_f_phys             = 1e-3;            % fluid viscosity                     (Pa.s)
k_eta_phys             = k_phys/eta_f_phys;    % permeability / fluid visc.     (m2/Pa/s)
eta_phys               = 1e22;            % Viscosity Solid                    (Pa.s)
dL_phys                = L_phys/10;       % Lengthscale of Anomaly              (m) (for non-ellipsoidal anomaly)
dP_phys                = 3e8;%0.1*P_phys;      % magnitude of P perturbation         (Pa)
Ebg_phys               = 2.0e-13;         % background strain rate              (1/s)
n                      = 3;               % Karman Kozeny exponent              ()
angle                  = 30;              % param for ellipsoid                 (dgr)
t_yield                = 1e8;             % yield stress for plasticity         (Pa)
% Independent Scales
L_sc                   = L_phys;                     % LenghScale               (m)
P_sc                   = P_phys;                     % Pscale                   (Pa)
t_sc                   = L_phys^2/k_eta_phys/P_phys; % timescale                (sec)
% Depentent Scales (shortcuts)
V_sc                   = L_sc/t_sc;     % velocity scale
eta_sc                 = P_sc*t_sc;     % viscosity scale
edot_sc                = 1/t_sc;        % strain rate scale
% Arbitrary scales (can be eliminated during iterations) !
rho_db                 = 3000;          % To be used only in mass balance
P0                     = 18e8/P_sc;%P_phys/P_sc;   % Ndim background P
% Nondimensional properties-----------------------------------------------
Ly                     = L_phys/L_sc;                   % Height of Model
Lx                     = Ly;                            % Width of Model
a                      = 0.05*Ly;                 % Parameters for ellipsoid anomaly
b                      = 3*a;                             % b = a --> circle
eta0                   = eta_phys/eta_sc;               % NDim initial solid viscosity
k                      = k_phys/L_sc^2;                 % NDim permeability
eta_f                  = eta_f_phys/eta_sc;             % NDim fluid viscosity
k_eta                  = k_eta_phys/(L_sc^2/P_sc/t_sc); % NDim k_eta
dL                     = dL_phys/L_sc;                  % only used for non-ellipsodal anomaly
dP                     = dP_phys/P_sc;                  % NDim P anomaly
tau_yield              = t_yield/P_sc;                  % NDim yield stress
ebg                    = 1*Ebg_phys/edot_sc;            % NDim backgrounhd strain rate
% ------------------------- Mesh setup and Numerics ------------------------- %
nvx                   = 151;       nvy    = 151;        % numerical resolution
ncx                   = nvx-1;     ncy    = nvy-1;
dx                    = Ly/ncy;    dy     = Ly/ncy;  # ERROR : dx
% Making the grids ------------------- %
xv                    = [-Lx/2:dx:Lx/2];            yv     = [-Ly/2:dy:Ly/2]; 
xRPf                  = [-Lx/2+dx:dy:Lx/2-dx];      yRPf   = [-Ly/2+dy:dy:Ly/2-dy];
xce                   = [-Lx/2-dx/2:dy:Lx/2+dx/2];  yce    = [-Ly/2-dy/2:dy:Ly/2+dy/2];
[X, Y]                = ndgrid(xv, yv);                                                         % create 2D grids
[X_vx, Y_vx]          = ndgrid(xce, yv);                                                        % velocity grids larger in one dimension for staggering
[X_vy, Y_vy]          = ndgrid(xv, yce);
[xRPf, yRPf]          = ndgrid(xRPf, yRPf);
X_incl                =  a*cos([0:0.01:2*pi])*cosd(angle)+b*sin([0:0.01:2*pi])*sind(angle);     % define ellispoidal inclusion
Y_incl                = -a*cos([0:0.01:2*pi])*sind(angle)+b*sin([0:0.01:2*pi])*cosd(angle);
incl_XY               = [-X_incl; Y_incl];
damp                  = 2*(Ly)/ncy;                                                             % convergence accelerator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USE ERWANS'S LOOK-UP TABLE FOR ECLOGITIZATION
load DATA_LOOKUP_ERWAN.mat
Pf_1d           = Pf_lt'    / P_sc;
rho_s_1d        = rhos_lt'  / rho_db;
rho_f_1d        = rhof_lt'  / rho_db;
X_s_1d          = Xs_lt';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------- Initial conditions ------------------------- %
% Pf perturbation ---------
Pf_ini                 = P0;
Pf                     = Pf_ini*ones(nvx, nvy);
rotx                                  =  X*cosd(angle)+Y*sind(angle);
roty                                  = -X*sind(angle)+Y*cosd(angle);
rotx2                                 =  xRPf*cosd(angle)+yRPf*sind(angle);
roty2                                 = -xRPf*sind(angle)+yRPf*cosd(angle);
Pf(sqrt(rotx.^2./a^2+roty.^2./b^2)<1) = Pf_ini + dP;
% Smoothing initial perturbation
for s=1:6
    Ii                  = [2:nvx-1];
    Pf(Ii,:)            = Pf(Ii,:) + 0.4*(Pf(Ii+1,:)-2*Pf(Ii,:)+Pf(Ii-1,:));
    Pf(:,Ii)            = Pf(:,Ii) + 0.4*(Pf(:,Ii+1)-2*Pf(:,Ii)+Pf(:,Ii-1));
end
P                       = Pf;

% Initial porosity consistent with initial pressure perturbation
[X_s, rho_s, rho_f]     = interp(Pf_1d,X_s_1d,rho_s_1d,rho_f_1d,Pf);
rhos_amb                = min(rho_s(:)) * ones(size(rho_s));
X_s_amb                 = max(X_s(:))   * ones(size(rho_s));
phi0                    = 1 - ( ( rhos_amb.*X_s_amb*(1-phi_amb) )./( rho_s.*X_s ) );
phi                     = phi0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta_mod                = 1e-11 * P_sc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ks_mod                  = 1/beta_mod;                                % solid bulk modulus
Kd                      = 0.5*Ks_mod;                                            % drained bulk modulus
G                       = Kd*1;                                                    % shear modulus
eta                     = eta0.*(phi0 ./ phi).^2;                                % porosity-dependent shear viscosity
lam                     = 2.*eta;                                                % bulk viscosity
% Characteristic scales
c_tsc         = (a^2*eta_f)/(phi_amb^3*Ks_mod*k);                                   % characteristic time in model
c_tsc_phys    = c_tsc*t_sc;                                                      % char. time in nature
c_Lsc         = sqrt((k*phi_amb^3/eta_f)*(max(max(lam)) + 4/3*max(max(eta))));      % char. length scale
Maxw_tsc      = eta0/Ks_mod;                                                     % Maxwell time
Maxw_tsc_phys = Maxw_tsc*t_sc;
% Physical timestep
dt                     = 1e-7*c_tsc;                                             % physical timestep in model units defined as a fraction of char. diffusion time scale
dt_phys                = dt*t_sc;                                                % phyisical timestep in SI units
t_kin                  = 1e-3*c_tsc;                                             % kinetic time scale defined as a fraction of diff. time scale
% Velocity ini -------------------------------
% Pure shear
Vx              = -ebg * X_vx;
Vy              =  ebg * Y_vy;
% Initialize further
rho_t           = rho_f .* phi + rho_s .* (1-phi);
rho_x           = rho_s .* X_s;
k_darcy         = k_eta*ones(nvx, nvy);
ax              = [-5*(L_phys/Lx) 5*(L_phys/Lx) -5*(L_phys/Ly) 5*(L_phys/Ly)];

it              = 0;        t_mod           = 0;        t_phys          = 0;
Pf_ini          = Pf;       phi_ini         = phi;
max_tii         = [0];      max_Pf          = [0];      strength        = [0];
max_rho_s       = [0];      av_eta_eff      = [0];      av_phi          = [0];
min_Peff        = [0];      itercount       = 0;        ttot            = 101*dt;
% Array allocation ------------------------------------------------------ %
pre_darcy_x = zeros(ncx, nvy); pre_darcy_y = zeros(nvx, ncy);     V_fluid_x = zeros(ncx, nvy); V_fluid_y = zeros(nvx, ncy);
evol        = zeros(nvx, nvy); exx         = zeros(nvx, nvy);     eyy       = zeros(nvx, nvy); exy       = zeros(ncx, ncy);
txx         = zeros(nvx, nvy); tyy         = zeros(nvx, nvy);     txy       = zeros(ncx, ncy); tii       = zeros(nvx, nvy); eii = zeros(nvx, nvy);
R_phi       = zeros(nvx, nvy); R_Pf        = zeros(ncx-1, ncx-1); Rvx       = zeros(ncx, ncy-1); Rvy     = zeros(ncx-1, ncy);
divv_rho_t  = zeros(nvx, nvy); divv_rho_x  = zeros(nvx, nvy);     dVxdtau   = zeros(ncx, ncy-1); dVydtau = zeros(ncx-1, ncy);
dphidtau    = zeros(nvx, nvy); dPfdtau     = zeros(ncx-1, ncy-1); P_eff     = zeros(nvx, nvy);   dPdtau = zeros(nvx, nvy);

divv_rho_x  = rand(size(divv_rho_x)) * 1e-1; % Perturbation to have initial error at first iteration different to zero

% ------------------------- Phyisical time loop ------------------------- %
while t_mod < ttot
    it        = it +1;
    rho_t_old = rho_t;
    rho_x_old = rho_x;
    phi_old   = phi;
    X_s_old   = X_s;
    rho_s_old = rho_s;
    rho_f_old = rho_f;
    txx_old   = txx;
    tyy_old   = tyy;
    txy_old   = txy;
    P_old     = P;
    Pf_old    = Pf;
    dtau      = 0.1*dy^2/max(max(eta))        * 6e4;                                    % pseudo timestep for solid velocity
    dtau_P    = 25.*a.*max(max(eta)).*dx/Lx   / 1e3;
    dtau_Pf   = 0.01*dy^2/max(max((k_darcy.*phi.^3.*Ks_mod)/1e-3)) * 5e4;         % pseudo timestep for fluid pressure
    dtau_phi  = dt /10;                                                       % pseudo timestep for porosity
    tol       = 1e-3;                                                      % tolerance of PT loop
    err_tot   = 1e10;
    iter      = 0;
    % ------------------------- Pseudo-time loop -------------------------  %
    tic
    while iter < 1e5 && err_tot>tol
        iter                 = iter+1;
        itercount            = itercount+1;
        max_Eta   = max(max(eta));
        % Conservation of non-volatile mass - porosity evolution
        rho_x                = rho_s .* X_s;
        R_phi                = ((1-phi).*rho_x - (1-phi_old).*rho_x_old)/dt + divv_rho_x;
        dphidtau             = R_phi + (1.0 - damp) .*dphidtau;
        phi                  = phi + dtau_phi.*dphidtau;
        % Conservation of total mass - fluid pressure evolution
        rho_t                = rho_f .* phi + rho_s .* (1-phi);
        pre_darcy_x          = (rho_f(1:end-1, :) .* k_darcy(1:end-1, :) .* phi(1:end-1, :).^3 + rho_f(2:end, :) .* k_darcy(2:end, :) .* phi(2:end, :).^3)/2;
        pre_darcy_y          = (rho_f(:, 1:end-1) .* k_darcy(:, 1:end-1) .* phi(:, 1:end-1).^3 + rho_f(:, 2:end) .* k_darcy(:, 2:end) .* phi(:, 2:end).^3)/2;
        q_darcy_x            = - pre_darcy_x .*(diff(Pf, 1, 1)/dx);
        q_darcy_y            = - pre_darcy_y .*(diff(Pf, 1, 2)/dy);
        V_fluid_x            = q_darcy_x./((rho_f(1:end-1, :)./phi(1:end-1, :) + rho_f(2:end, :)./phi(2:end, :))/2);
        V_fluid_y            = q_darcy_y./((rho_f(:, 1:end-1)./phi(:, 1:end-1) + rho_f(:, 2:end)./phi(:, 2:end))/2);
        div_q_darcy          = diff(q_darcy_x(:, 2:end-1), 1, 1)/dx + diff(q_darcy_y(2:end-1, :), 1, 2)/dy;
        R_Pf                 = - div_q_darcy - (rho_t(2:end-1, 2:end-1) - rho_t_old(2:end-1, 2:end-1))/dt - divv_rho_t(2:end-1, 2:end-1);
        dPfdtau              = R_Pf + (1.0 - damp) .*dPfdtau;
        Pf(2:end-1, 2:end-1) = Pf(2:end-1, 2:end-1) + dtau_Pf .* dPfdtau;
        % No flux at boundaries in y dimension (to avoid late-stage boundary effects as reaction progresses vertically)
        Pf(:,1)              = Pf(:,2);
        Pf(:, end)            = Pf(:, end-1);
        % Lookup based on Pf and X
        [X_s_eq, rho_s_eq, rho_f]  = interp(Pf_1d,X_s_1d,rho_s_1d,rho_f_1d,Pf);
        % Reaction kinetics
        X_s                  = X_s_old   + dt.*(X_s_eq   - X_s)  ./t_kin;
        rho_s                = rho_s_old + dt.*(rho_s_eq - rho_s)./t_kin;
        % Strain rates, viscosity, deviatoric stresses
        divv_rho_x           = diff(yv2yce((1-phi).*rho_x, 1).*Vx, 1, 1)/dx + diff(yv2yce((1-phi).*rho_x, 2).*Vy, 1, 2)/dy;
        divv_rho_t           = diff(yv2yce(rho_t, 1).*Vx, 1, 1)/dx + diff(yv2yce(rho_t, 2).*Vy, 1, 2)/dy;
        evol                 = diff(Vx, 1, 1)/dx + diff(Vy, 1, 2)/dy;                                       % volumetric strain rate
        exx                  = diff(Vx,1,1)/dx - 1/3*evol;                                                  % viscous deviatoric strain rate in x
        eyy                  = diff(Vy,1,2)/dy - 1/3*evol;
        exy                  = 0.5*(diff(Vx(2:end-1,:),1,2)/dy + diff(Vy(:,2:end-1),1,1)/dx);
        % visco-elastic deviatoric strain rates
        Exx                  = exx + txx_old./2./(dt*G);                                                    % viscoelastic deviatoric strain rate in x
        Eyy                  = eyy + tyy_old./2./(dt*G);
        Exy                  = exy + txy_old./2./(dt*G);
        % porostiy-dependent viscosity
        eta                  = eta0.*(phi0 ./ phi).^2;                                                      % porosity-dependent shear viscosity
        lam                  = 2.*eta;                                                                      % bulk viscosity
        eta_xy               = (eta(2:end,:)    + eta(1:end-1,:))/2;
        eta_xy               = (eta_xy(:,2:end) + eta_xy(:,1:end-1))/2;
        % visco-elastic modulus
        eta_ve               = 1./(1/(dt*G) + 1./eta);                                                      % viscoelastic modulus
        eta_ve_xy            = 1./(1/(dt*G) + 1./eta_xy);
        % visco-elastic deviatoric stresses
        txx                  = 2 .*eta_ve    .* Exx;                                                        % deviatoric stress in x
        tyy                  = 2 .*eta_ve    .* Eyy;
        txy                  = 2 .*eta_ve_xy .* Exy;
        tii                  = sqrt( 0.5.*(txx.^2 + tyy.^2) + yv2yce(yv2yce(txy,1),2).^2);                  % 2nd invariant of dev. stress tensor
        eii                  = sqrt( 0.5.*(exx.^2 + eyy.^2) + yv2yce(yv2yce(exy,1),2).^2);                  % 2nd invariant of dev. strain rate tensor
        % ------------------------ option for using Von Mises Plasticity -----------------------------

        % -----------------------------------------------------------------
        % Total Pressure
        R_P                  = (P - (P_old - dt .* (Kd .* evol - 0.5 .* ((Pf - Pf_old) ./ dt) ...
            - Kd .* Pf ./ ((1 - phi) .* lam))) ./ (1 + dt .* Kd ./ ((1.0 - phi) .* lam)));
        P                    = P - dtau_P .* R_P;
        % P                    = (P_old - dt .* (Kd .* evol - 0.5 .* ((Pf - Pf_old) ./ dt) ...
        %     - Kd .* Pf ./ ((1 - phi) .* lam))) ./ (1 + dt .* Kd ./ ((1.0 - phi) .* lam));
        % Effective Pressure
        P_eff                = P-Pf;
        % Velocity Residuals
        Rvx                  =   -diff(P(:,2:end-1),1,1)/dx + diff(txx(:,2:end-1),1,1)/dx + diff(txy,1,2)/dy;
        Rvy                  =   -diff(P(2:end-1,:),1,2)/dy + diff(tyy(2:end-1,:),1,2)/dy + diff(txy,1,1)/dx;
        % Convergence accelerator
        dPdtau               = R_P - (1.0 - damp) .* dPdtau;
        dVxdtau              = Rvx + (1.0 - damp) .* dVxdtau;
        dVydtau              = Rvy + (1.0 - damp) .* dVydtau;
        % Update Velocities
        Vx(2:end-1, 2:end-1) = Vx(2:end-1, 2:end-1) + dVxdtau * dtau;
        Vy(2:end-1, 2:end-1) = Vy(2:end-1, 2:end-1) + dVydtau * dtau;
        P                    = P + dPdtau * dtau_P;
        err_Pf               = norm(R_Pf(:))  / sqrt(length(R_Pf(:)));     	% Error fluid pressure
        err_P                = norm(R_P(:))   / sqrt(length(R_P(:)));     	% Error fluid pressure
        err_phi              = norm(R_phi(:)) / sqrt(length(R_phi(:)));                   % Error phi
        err_Vx               = norm(Rvx(:))   / sqrt(length(Rvx(:)));        % Error horizontal velocitiy
        err_Vy               = norm(Rvy(:))   / sqrt(length(Rvy(:)));        % Error vertical velocity
        if iter == 1
            err_Pf_0         = err_Pf;     	% Error fluid pressure
            err_P_0          = err_P;     	% Error fluid pressure
            err_phi_0        = err_phi;                   % Error phi
            err_Vx_0         = err_Vx;        % Error horizontal velocitiy
            err_Vy_0         = err_Vy;        % Error vertical velocity
        end

        err_tot_vec          = [err_Pf/err_Pf_0, err_Vx/err_Vx_0, err_Vy/err_Vy_0, err_phi/err_phi_0, err_P/err_P_0];           % Error total
        err_tot              = max(err_tot_vec);           % Error total
        %fprintf('it = %d, iter = %d, err = %1.3e norm[Rvx=%1.3e, Rvy=%1.3e, R_Pf=%1.3e, R_phi=%1.3e] \n', it, iter, err_tot, err_Vx, err_Vy, err_Pf, err_phi);
        % Plot error evolution
        if mod(itercount, 2500)==0
            figure(1)
            plot(itercount,log10(err_tot_vec(2)),'kx'),hold on;plot(itercount,log10(err_tot_vec(3)),'rv');plot(itercount,log10(err_tot_vec(1)),'bd');plot(itercount,log10(err_tot_vec(end-1)),'go');plot(itercount,log10(err_tot_vec(end)),'ms');
            plot([0 itercount+500],[log10(tol) log10(tol)],'--r','linewidth',1)
            legend('Error Vx','Error Vy','Error Pf','Error phi','Error P','Tolerance','Location','Northwest'); grid on;xlabel('total PT iterations');ylabel('log_{10} Error')
            drawnow
            fprintf('it = %d, iter = %d, err = %1.3e norm[Rvx=%1.3e, Rvy=%1.3e, R_Pf=%1.3e, R_phi=%1.3e, R_P=%1.3e] \n', it, iter, err_tot, err_tot_vec(2), err_tot_vec(3), err_tot_vec(1), err_tot_vec(4), err_tot_vec(5));
        end
    end
    toc
    % Track variables
    max_tii(it)          = max(max(tii(5:end-4, 5:end-4)));                                    % exclude a few cells along boundaries to not track possible late-stage boundary effects
    max_Pf(it)           = max(max(Pf(5:end-4, 5:end-4)));
    strength(it)         = mean(mean(((tii(5:end-4, 5:end-4)*P_sc)./(eii(5:end-4, 5:end-4)*edot_sc))));
    max_rho_s(it)        = max(max(rho_s(5:end-4, 5:end-4)));
    av_phi(it)           = mean(mean(phi(5:end-4, 5:end-4)));
    min_Peff(it)         = min(min(P_eff(5:end-4, 5:end-4)));
    % ---------------- Plot rho_s, Pf, phi, evol -------------------------
    if mod(it, 2) == 0
        figure(2)
        subplot(221), pcolor(X, Y, rho_s*rho_db), xlabel('X'), ylabel('Y'),
        title(['rho_s ', 'it = ', num2str(it), ', t phys =', num2str(t_phys/(60*60)), ' hrs']),
        shading interp,
        colorbar, hold on
        plot(incl_XY(1,:),incl_XY(2,:),':r','linewidth',1.3)
        axis equal; axis(ax)
        subplot(222), pcolor(X, Y, Pf*P_sc), xlabel('X'), ylabel('Y'),
        title('Pf'), shading interp, colorbar, hold on
        plot(incl_XY(1,:),incl_XY(2,:),':r','linewidth',1.3)
        axis equal; axis(ax)
        subplot(223), pcolor(X, Y, phi), xlabel('X'), ylabel('Y'),
        title('phi'), shading interp, colorbar, hold on,
        plot(incl_XY(1,:),incl_XY(2,:),':r','linewidth',1.3)
        axis equal; axis(ax)
        subplot(224), pcolor(X, Y, evol*dt/dt_phys), xlabel('X'), ylabel('Y'),
        title('volumetric def'), shading interp, colorbar; hold on,
        plot(incl_XY(1,:),incl_XY(2,:),':r','linewidth',1.3)
        axis equal; axis(ax)
        drawnow
        % ---------------- Plot tii, eii, V_s, V_f ---------------------------
        figure(3)
        subplot(221), pcolor(X, Y, tii.*P_sc), xlabel('X'), ylabel('Y'),
        title('tii'), shading interp
        colorbar, hold on,
        plot(incl_XY(1,:),incl_XY(2,:),':r','linewidth',1.3)
        axis equal; axis(ax)
        subplot(222), pcolor(X, Y, eii*dt/dt_phys), xlabel('X'), ylabel('Y'),
        title('eii'), shading interp, colorbar, hold on,
        plot(incl_XY(1,:),incl_XY(2,:),':r','linewidth',1.3)
        axis equal; axis(ax)
        subplot(223), pcolor(xRPf,yRPf,sqrt(V_fluid_x(2:end,2:end-1).^2+V_fluid_y(2:end-1, 2:end).^2).*V_sc), xlabel('X'), ylabel('Y'),
        title('V fluid'), shading interp
        %colormap(cmap.tokyo)
        colorbar, hold on
        quiver(xRPf, yRPf, V_fluid_x(2:end,2:end-1).*V_sc, V_fluid_y(2:end-1, 2:end).*V_sc, 'b'), hold on
        plot(incl_XY(1,:),incl_XY(2,:),':r','linewidth',1.3)
        axis equal; axis(ax)
        subplot(224), pcolor(X, Y, P_eff*P_sc), xlabel('X'), ylabel('Y'),
        title('P effective'), shading interp, colorbar, hold on
        plot(incl_XY(1,:),incl_XY(2,:),':r','linewidth',1.3)
        axis equal; axis(ax)
        drawnow
        % ------------- Plot solid velocity components, total P, eta ----------
        figure(4)
        subplot(221), pcolor(X_vx, Y_vx, Vx*V_sc), xlabel('X'), ylabel('Y'),
        title('Vx'), shading interp
        colorbar, hold on
        plot(incl_XY(1,:),incl_XY(2,:),':r','linewidth',1.3)
        axis equal; axis(ax)
        subplot(222), pcolor(X_vy, Y_vy, Vy*V_sc), xlabel('X'), ylabel('Y'),
        title('Vy'), shading interp,colorbar, hold on
        plot(incl_XY(1,:),incl_XY(2,:),':r','linewidth',1.3)
        axis equal; axis(ax)
        subplot(223), pcolor(X, Y, P*P_sc), xlabel('X'), ylabel('Y'),
        title('P'), shading interp,colorbar, hold on
        plot(incl_XY(1,:),incl_XY(2,:),':r','linewidth',1.3)
        axis equal; axis(ax)
        subplot(224), pcolor(X, Y, eta*eta_sc), xlabel('X'), ylabel('Y'),
        title('viscosity'), shading interp,colorbar, hold on
        plot(incl_XY(1,:),incl_XY(2,:),':r','linewidth',1.3)
        axis equal; axis(ax)
        drawnow
        % -------------------- Plot tracked variables ------------------------
        figure(5)
        subplot(231), plot([1:1:it], max_tii*P_sc,     'o', Color='k'), title('max tii')
        subplot(232), plot([1:1:it], strength/2,       'o', Color='k'), title('mean tii/eii/2')
        subplot(233), plot([1:1:it], min_Peff*P_sc,    'o', Color='k'), title('min P effective')
        subplot(234), plot([1:1:it], max_Pf*P_sc,      'o', Color='r'), title('max Pf')
        subplot(235), plot([1:1:it], max_rho_s*rho_db, 'o', Color='r'), title('max rho solid')
        subplot(236), plot([1:1:it], av_phi,           'o', Color='r'), title('mean phi')
        drawnow
    end
    % ------------ Save model parameters -----------------------------------
    if SAVE == 1
        if mod(it, 100)==0 || it == 1
            mkdir 'C:\Users\porkolab\Documents\PD 22\Results\paper 1 models\ref_vscelast_dt8em9_tol1em10';
            dir = 'C:\Users\porkolab\Documents\PD 22\Results\paper 1 models\ref_vscelast_dt8em9_tol1em10';
            save(fullfile(dir, ['Pf_', num2str(it), '.mat']), 'Pf', '-mat');
            save(fullfile(dir, ['P_tot_', num2str(it), '.mat']), 'P', '-mat');
            save(fullfile(dir, ['P_eff_', num2str(it), '.mat']), 'P_eff', '-mat');
            save(fullfile(dir, ['rhos_', num2str(it), '.mat']), 'rho_s', '-mat');
            save(fullfile(dir, ['phi_', num2str(it), '.mat']), 'phi', '-mat');
            save(fullfile(dir, ['eii_', num2str(it), '.mat']), 'eii', '-mat');
            save(fullfile(dir, ['eta_', num2str(it), '.mat']), 'eta', '-mat');
            save(fullfile(dir, ['tii_', num2str(it), '.mat']), 'tii', '-mat');
            save(fullfile(dir, ['evol_', num2str(it), '.mat']), 'evol', '-mat');
            save(fullfile(dir, ['Vx_', num2str(it), '.mat']), 'Vx', '-mat');
            save(fullfile(dir, ['Vy_', num2str(it), '.mat']), 'Vy', '-mat');
            save(fullfile(dir, ['V_f_x_', num2str(it), '.mat']), 'V_fluid_x', '-mat');
            save(fullfile(dir, ['V_f_y_', num2str(it), '.mat']), 'V_fluid_y', '-mat');
            save(fullfile(dir, 'av_phi.mat'), 'av_phi', '-mat');
            save(fullfile(dir, 'max_Pf.mat'), 'max_Pf', '-mat');
            save(fullfile(dir, 'max_rho_s.mat'), 'max_rho_s', '-mat');
            save(fullfile(dir, 'max_tii.mat'), 'max_tii', '-mat');
            save(fullfile(dir, 'strength.mat'), 'strength', '-mat');
            save(fullfile(dir, 'min_Peff.mat'), 'min_Peff', '-mat');
        end
    end
    t_mod  = t_mod+dt;
    t_phys = t_phys+dt_phys;
end
% ----------------------------- additional functions -------------------
% Interpolate from yv to yce
function [B] = yv2yce(A, ndim)
if ndim == 1
    B = zeros(size(A,1)+1,size(A,2));
    B(2:end-1,:)   = (A(2:end,:)   + A(1:end-1,:))/2;
    B(1,  :)       = 1.5*A(1,:)   -0.5*A(2,:);
    B(end,:)       = 1.5*A(end,:) -0.5*A(end-1,:);
elseif ndim == 2
    B = zeros(size(A,1),size(A,2)+1);
    B(:,2:end-1)   = (A(:,2:end) + A(:,1:end-1))/2;
    B(:,1  )       = 1.5*A(:,1)   -0.5*A(:,2);
    B(:,end)       = 1.5*A(:,end) -0.5*A(:,end-1);
end
end

% Interpolate from yce to yv
function [B] = yce2yv(A, ndim)
if ndim == 1
    B = zeros(size(A,1)-1,size(A,2));
    B(1:end,:)     = (A(2:end,:)   + A(1:end-1,:))/2;
elseif ndim == 2
    B = zeros(size(A,1),size(A,2)-1);
    B(:,1:end)     = (A(:,2:end) + A(:,1:end-1))/2;
end
end

% Extract 1D coordinates from lookup tables for interpolation
function [Pf_1d,X_s_1d,rho_s_1d,rho_f_1d] = extract1d(X_coord,Pf_lt,X_s_lt,rho_s_lt, rho_f_lt)
Pf_1d       = Pf_lt(X_coord,:);
X_s_1d      = X_s_lt(X_coord,:);
rho_s_1d    = rho_s_lt(X_coord,:);
rho_f_1d    = rho_f_lt(X_coord,:);
end

% 1D interpolation using Pf values from total mass balance equation
function [X_s_eq, rho_s_eq, rho_f] = interp(Pf_1d,X_s_1d,rho_s_1d,rho_f_1d,Pf)
X_s_eq     = interp1(Pf_1d,X_s_1d,Pf,'linear');
rho_s_eq   = interp1(Pf_1d,rho_s_1d,Pf,'linear');
rho_f      = interp1(Pf_1d,rho_f_1d,Pf,'linear');
end