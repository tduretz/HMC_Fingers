using Plots, MAT, Printf
import LinearAlgebra:norm
#using Infiltrator 

# Interpolate from yv to yce
@views function yv2yce!(B, A, ndim)
    if ndim == 1
        @. B[2:end-1,:]   = (A[2:end,:]   + A[1:end-1,:])/2;
        B[1,  :]       .= 1.5*A[1,:]   .-0.5*A[2,:]
        B[end,:]       .= 1.5*A[end,:] .-0.5*A[end-1,:]
    elseif ndim == 2
        @. B[:,2:end-1]   = (A[:,2:end] + A[:,1:end-1])/2;
        B[:,1  ]       .= 1.5*A[:,1]   .-0.5*A[:,2];
        B[:,end]       .= 1.5*A[:,end] .-0.5*A[:,end-1];
    end
end

# Interpolate from yce to yv
@views function yce2yv!(B, A, ndim)
    if ndim == 1
        @. B[1:end,:]     = (A[2:end,:] + A[1:end-1,:])/2;
    elseif ndim == 2
        @. B[:,1:end]     = (A[:,2:end] + A[:,1:end-1])/2;
    end
end

@views function Itp1D( Xc, Pc, Xdata, Pdata, p )
    dP   =  Pdata[2] - Pdata[1]
    Pmin =  Pdata[1]
    # Interpolate in 1D:              X_W -----------o---x_E
    p.iW .= Int64.(round.(trunc.( (Pc .-  Pmin)./dP ) .+ 1))   # index of the west node is database
    p.wW .= 1.0 .- (Pc .- Pdata[p.iW])./dP
    for i in  eachindex(Xc)
        Xc[i] = p.wW[i] * Xdata[p.iW[i]] +  (1.0 - p.wW[i]) * Xdata[p.iW[i] + 1]
    end
    # Xc .= @views wW .* Xdata[iW] .+  (1.0 .- wW) .* Xdata[iW.+1]
end

function main(;nt=1, tol=1e-4, max_its=1e5)
    # Physical Parameters
    ϕ_amb                  = 0.01
    L_phys                 = 0.1             # Physical Length                     (m)
    P_phys                 = 19e8            # Reaction Pressure                   (Pa)
    k_phys                 = 1e-22           # permeability                        (m2)
    η_f_phys               = 1e-3            # fluid viscosity                     (Pa.s)
    k_η_phys               = k_phys/η_f_phys    # permeability / fluid visc.     (m2/Pa/s)
    η_phys                 = 1e22            # Viscosity Solid                     (Pa.s)
    dL_phys                = L_phys/10       # Lengthscale of Anomaly              (m) (for non-ellipsoidal anomaly)
    dP_phys                = 3e8#0.1*P_phys      # magnitude of P perturbation        (Pa)
    Ebg_phys               = 1.0e-14         # background strain rate              (1/s)
    n_darcy                = 3               # Karman Kozeny exponent              ()
    angle                  = -30              # param for ellipsoid                 (dgr)
    t_yield                = 1e8             # yield stress for plasticity         (Pa)
    # Independent Scales
    L_sc                   = L_phys                     # LenghScale               (m)
    P_sc                   = P_phys                     # Pscale                   (Pa)
    t_sc                   = L_phys^2/k_η_phys/P_phys # timescale                (sec)
    # Dependent Scales (shortcuts)
    V_sc                   = L_sc/t_sc     # velocity scale
    η_sc                   = P_sc*t_sc     # viscosity scale
    ε̇_sc                   = 1/t_sc        # strain rate scale
    # Arbitrary scales (can be eliminated during iterations) !
    ρ_db                   = 3000          # To be used only in mass balance
    P0                     = 18e8/P_sc#P_phys/P_sc   # Ndim background P
    # Nondimensional properties-----------------------------------------------
    Ly                     = L_phys/L_sc                   # Height of Model
    Lx                     = Ly                            # Width of Model
    a                      = 0.05*Ly                 # Parameters for ellipsoid anomaly
    b                      = 3*a                             # b = a --> circle
    η0                     = η_phys/η_sc               # NDim initial solid viscosity
    k                      = k_phys/L_sc^2                 # NDim permeability
    η_f                    = η_f_phys/η_sc             # NDim fluid viscosity
    k_η                    = k_η_phys/(L_sc^2/P_sc/t_sc) # NDim k_η
    dP                     = dP_phys/P_sc                  # NDim P anomaly
    ε̇bg                    = 1*Ebg_phys/ε̇_sc            # NDim backgrounhd strain rate
    # ------------------------- Mesh setup and Numerics ------------------------- #
    nv                     = (x=151,       y    = 151)        # numerical resolution
    nc                     = (x=nv.x-1,    y    = nv.y-1)
    Δ                      = (x=Lx/nc.x,   y    = Ly/nc.y)
    # Making the grids ------------------- #
    xv                    = LinRange(-Lx/2, Lx/2, nv.x);   yv      = LinRange(-Ly/2, Ly/2, nv.y);
    xc_ex                    = LinRange(-Lx/2-Δ.x/2, Lx/2+Δ.x/2, nc.x+2);  yc_ex    = LinRange(-Ly/2-Δ.y/2, Ly/2+Δ.y/2, nc.y+2);
    θ                     = LinRange(0., 2π, 100)
    Y_incl                =  a*cos.(θ)*cosd(angle) .+ b*sin.(θ)*sind(angle)     # define ellispoidal inclusion
    X_incl                = -a*cos.(θ)*sind(angle) .+ b*sin.(θ)*cosd(angle)
    damp                  = 2*(Ly)/nc.y
    ### USE ERWANS'S LOOK-UP TABLE FOR ECLOGITIZATION
    file   = matopen("./data/DATA_LOOKUP_ERWAN.mat")
    Pf_1d  = read(file, "Pf_lt")   ./ P_sc # note that this does NOT introduce a variable ``varname`` into scope
    ρ_s_1d = read(file, "rhos_lt") ./ ρ_db
    ρ_f_1d = read(file, "rhof_lt") ./ ρ_db
    X_s_1d = read(file, "Xs_lt")
    close(file)
    #@infiltrate
    # ------------------------- Initial conditions ------------------------- %
    # Pf perturbation ---------
    Pf_ini                 = P0
    Pf                     = Pf_ini*ones(nv...)
    roty                                  =  xv*cosd(angle) .+ yv'*sind(angle)
    rotx                                  = -xv*sind(angle) .+ yv'*cosd(angle)
    @. Pf[sqrt(rotx^2 /a^2 + roty^2/b^2) < 1] = Pf_ini + dP
    # Smoothing initial perturbation
    for s=1:6
        @. Pf[2:end-1,:] += 0.4*(Pf[3:end,:] - 2*Pf[2:end-1,:] + Pf[1:end-2,:]);
        @. Pf[:,2:end-1] += 0.4*(Pf[:,3:end] - 2*Pf[:,2:end-1] + Pf[:,1:end-2]);
    end
    P                       = copy(Pf)
    # Initial porosity consistent with initial pressure perturbation
    itp      = (wW = zeros(nv...), iW = zeros(Int64, nv...) )
    X_s      = zero(Pf)
    ρ_s      = zero(Pf)
    ρ_f      = zero(Pf)
    @views Itp1D( X_s, Pf, X_s_1d, Pf_1d, itp )
    @views Itp1D( ρ_s, Pf, ρ_s_1d, Pf_1d, itp )
    @views Itp1D( ρ_f, Pf, ρ_f_1d, Pf_1d, itp )
    ρ_s_amb = minimum(ρ_s)
    X_s_amb = maximum(X_s)
    ϕ0      = 1.0 .- ( ( ρ_s_amb.*X_s_amb*(1-ϕ_amb) )./( ρ_s.*X_s ) )
    ϕ       = copy(ϕ0)
    #######################################################
    β_mod              = 1e-11 * P_sc
    Ks_mod             = 1/β_mod                                          # solid bulk modulus
    Kd                 = 0.5*Ks_mod                                       # drained bulk modulus
    G                  = Kd                                               # shear modulus
    η                  = η0.*(ϕ0 ./ ϕ).^2                                 # porosity-dependent shear viscosity
    λ                  = 2 .*η                                            # bulk viscosity
    # Characteristic scales
    c_tsc         = (a^2*η_f)/(ϕ_amb^3*Ks_mod*k)                          # characteristic time in model
    c_tsc_phys    = c_tsc*t_sc                                            # char. time in nature
    c_Lsc         = sqrt((k*ϕ_amb^3/η_f)*(maximum(λ) + 4/3*maximum(η)))   # char. length scale
    Maxw_tsc      = η0/Ks_mod                                             # Maxwell time
    Maxw_tsc_phys = Maxw_tsc*t_sc
    # Physical timestep
    dt                     = 2e-7*c_tsc                                   # physical timestep in model units defined as a fraction of char. diffusion time scale
    dt_phys                = dt*t_sc                                      # phyisical timestep in SI units
    t_kin                  = 1e-3*c_tsc                                   # kinetic time scale defined as a fraction of diff. time scale
    # Velocity ini -------------------------------
    # Pure shear
    Vx              = -ε̇bg * xc_ex .+ 0   * yv'
    Vy              =  0   * xv .+ ε̇bg * yc_ex'
    # Initialize further
    ρ_t             = ρ_f .* ϕ   .+ ρ_s .* (1 .- ϕ)
    ρ_x             = ρ_s .* X_s
    k_darcy         = k_η*ones(nv...)
    # Allocate more
    τxx      = zeros(nv...)
    τyy      = zeros(nv...)
    τxy      = zeros(nc...) 
    ρ_t_old  = zero(ρ_t)
    ρ_x_old  = zero(ρ_x)
    ϕ_old    = zero(ϕ)
    X_s_old  = zero(X_s)
    ρ_s_old  = zero(ρ_s)
    ρ_f_old  = zero(ρ_f)
    τxx_old  = zero(τxx)
    τyy_old  = zero(τyy)
    τxy_old  = zero(τxy)
    P_old    = zero(P  )
    Pf_old   = zero(Pf )
    R_ϕ      = zero(ϕ  )
    divv_ρ_x = zero(ϕ  )
    dϕdτ     = zero(ϕ  )
    pre_darcy_x = zeros(nv.x-1, nv.y)
    pre_darcy_y = zeros(nv.x, nv.y-1)
    q_darcy_x   = zeros(nv.x-1, nv.y)
    q_darcy_y   = zeros(nv.x, nv.y-1)
    V_fluid_x   = zeros(nv.x-1, nv.y)
    V_fluid_y   = zeros(nv.x, nv.y-1)
    div_q_darcy = zeros(nv.x-2, nv.y-2)
    R_Pf        = zeros(nv.x-2, nv.y-2)
    divv_ρ_t    = zero(ϕ  )
    dPfdτ       = zeros(nv.x-2, nv.y-2)
    X_s_eq      = zero(X_s)
    ρ_s_eq      = zero(ρ_s)
    Av          = zeros(nv.x, nv.y)
    Av_cx       = zeros(nv.x+1, nv.y)
    Av_cy       = zeros(nv.x, nv.y+1)
    divV        = zeros(nv...)
    exx         = zeros(nv...)
    eyy         = zeros(nv...)
    exy         = zeros(nc...)
    Exx         = zeros(nv...)
    Eyy         = zeros(nv...)
    Exy         = zeros(nc...)
    η_xy        = zeros(nc...)
    η_ve_xy     = zeros(nc...)
    η_ve        = zeros(nv...)
    R_P         = zeros(nv...)
    P_eff       = zeros(nv...)
    Rvx         = zeros(nc.x, nv.y-2)      
    Rvy         = zeros(nv.x-2, nc.y)  
    dPdτ        = zeros(nv...)
    dVxdτ       = zeros(nc.x, nv.y-2)
    dVydτ       = zeros(nv.x-2, nc.y)    
    # Time loop
    for it=1:nt
        @printf("Time step = %06d\n", it)
        ρ_t_old .= ρ_t
        ρ_x_old .= ρ_x
        ϕ_old   .= ϕ
        X_s_old .= X_s
        ρ_s_old .= ρ_s
        ρ_f_old .= ρ_f
        τxx_old .= τxx
        τyy_old .= τyy
        τxy_old .= τxy
        P_old   .= P
        Pf_old  .= Pf
        max_η   = maximum(η)
        @show dτ      = 0.1*Δ.y^2/maximum(η) * 5e6                               # pseudo timestep for solid velocity
        @show dτ_P    = 25. *a.*maximum(η).*Δ.y/Lx   / 2e3
        @show dτ_Pf   = 0.01*Δ.y^2/maximum((k_darcy.*ϕ.^3 .* Ks_mod)/1e-3) * 1e5         # pseudo timestep for fluid pressure
        @show dτ_ϕ    = dt /25 
        errVx0, errVy0 = 1., 1.
        for iter=1:max_its
            # Conservation of non-volatile mass - porosity evolution
            @. ρ_x  = ρ_s * X_s;
            @. R_ϕ  = ((1-ϕ).*ρ_x - (1-ϕ_old)*ρ_x_old)/dt + divv_ρ_x
            @. dϕdτ = R_ϕ + (1.0 - damp) .*dϕdτ
            @. ϕ    += dτ_ϕ.*dϕdτ;
            # Conservation of total mass - fluid pressure evolution
            @. ρ_t                  = ρ_f .* ϕ + ρ_s .* (1-ϕ);
            @. pre_darcy_x          = (ρ_f[1:end-1, :] .* k_darcy[1:end-1, :] .* ϕ[1:end-1, :].^n_darcy + ρ_f[2:end, :] .* k_darcy[2:end, :] .* ϕ[2:end, :].^3)/2
            @. pre_darcy_y          = (ρ_f[:, 1:end-1] .* k_darcy[:, 1:end-1] .* ϕ[:, 1:end-1].^n_darcy + ρ_f[:, 2:end] .* k_darcy[:, 2:end] .* ϕ[:, 2:end].^3)/2;
            @. q_darcy_x            = - pre_darcy_x .*(Pf[2:end,:] - Pf[1:end-1,:])/Δ.x
            @. q_darcy_y            = - pre_darcy_y .*(Pf[:,2:end] - Pf[:,1:end-1])/Δ.y
            @. V_fluid_x            = q_darcy_x./((ρ_f[1:end-1, :]./ϕ[1:end-1, :] + ρ_f[2:end, :]./ϕ[2:end, :])/2)
            @. V_fluid_y            = q_darcy_y./((ρ_f[:, 1:end-1]./ϕ[:, 1:end-1] + ρ_f[:, 2:end]./ϕ[:, 2:end])/2)
            @. div_q_darcy          = (q_darcy_x[2:end, 2:end-1] - q_darcy_x[1:end-1, 2:end-1])/Δ.x + (q_darcy_y[2:end-1,2:end] - q_darcy_y[2:end-1,1:end-1])/Δ.y;
            @. R_Pf                 = - div_q_darcy - (ρ_t[2:end-1, 2:end-1] - ρ_t_old[2:end-1, 2:end-1])/dt - divv_ρ_t[2:end-1, 2:end-1]
            @. dPfdτ                = R_Pf + (1.0 - damp) .*dPfdτ
            @. Pf[2:end-1, 2:end-1] = Pf[2:end-1, 2:end-1] + dτ_Pf .* dPfdτ
            # No flux at boundaries in y dimension (to avoid late-stage boundary effects as reaction progresses vertically)
            @. Pf[:,1]              = Pf[:,2]
            @. Pf[:, end]           = Pf[:, end-1]
            # Thermo
            @views Itp1D( X_s_eq, Pf, X_s_1d, Pf_1d, itp )
            @views Itp1D( ρ_s_eq, Pf, ρ_s_1d, Pf_1d, itp )
            @views Itp1D( ρ_f, Pf, ρ_f_1d, Pf_1d, itp )
            # Reaction kinetics
            @. X_s  = X_s_old + dt.*(X_s_eq - X_s  ) ./t_kin
            @. ρ_s  = ρ_s_old + dt.*(ρ_s_eq - ρ_s  ) ./t_kin
            # Strain rates, viscosity, deviatoric stresses
            @. Av = (1-ϕ)* ρ_x
            yv2yce!(Av_cx, Av, 1)
            yv2yce!(Av_cy, Av, 2)
            @. divv_ρ_x           = (Av_cx[2:end,:] - Av_cx[1:end-1,:])/Δ.x  + (Av_cy[:,2:end] - Av_cy[:,1:end-1])/Δ.y
            yv2yce!(Av_cx, ρ_t, 1)
            yv2yce!(Av_cy, ρ_t, 2)
            @. divv_ρ_t           = (Av_cx[2:end,:] - Av_cx[1:end-1,:])/Δ.x  + (Av_cy[:,2:end] - Av_cy[:,1:end-1])/Δ.y
            @. divV               = (Vx[2:end,:] - Vx[1:end-1,:])/Δ.x + (Vy[:,2:end] - Vy[:,1:end-1])/Δ.y                                       # volumetric strain rate
            @. exx                = (Vx[2:end,:] - Vx[1:end-1,:])/Δ.x - 1/3*divV;                                           # viscous deviatoric strain rate in x
            @. eyy                = (Vy[:,2:end] - Vy[:,1:end-1])/Δ.y - 1/3*divV;
            @. exy                = 0.5*((Vx[2:end-1,2:end] - Vx[2:end-1,1:end-1])/Δ.y + (Vy[2:end,2:end-1] - Vy[1:end-1,2:end-1])/Δ.x)
            # visco-elastic deviatoric strain rates
            @. Exx                = exx + τxx_old./2/(dt*G) # viscoelastic deviatoric strain rate in x
            @. Eyy                = eyy + τyy_old./2/(dt*G)
            @. Exy                = exy + τxy_old./2/(dt*G)
            # porosity-dependent viscosity
            @. η                  = η0.*(ϕ_amb ./ ϕ).^2   # porosity-dependent shear viscosity
            @. λ                  = 2*η;               # bulk viscosity
            @. η_xy               = 0.25*(η[1:end-1,1:end-1] + η[1:end-1,2:end] + η[2:end,1:end-1] + η[2:end,2:end])
            # visco-elastic modulus
            @. η_ve               = 1/(1/(dt*G) + 1/η)    # viscoelastic modulus
            @. η_ve_xy            = 1/(1/(dt*G) + 1/η_xy)
            # visco-elastic deviatoric stresses
            @. τxx                = 2 .*η_ve    .* Exx # deviatoric stress in x
            @. τyy                = 2 .*η_ve    .* Eyy
            @. τxy                = 2 .*η_ve_xy .* Exy
            # @. tii                  = sqrt( 0.5.*(txx.^2 + tyy.^2) + yv2yce(yv2yce(txy,1),2).^2);                  % 2nd invariant of dev. stress tensor
            # @. eii                  = sqrt( 0.5.*(exx.^2 + eyy.^2) + yv2yce(yv2yce(exy,1),2).^2);
            # Total Pressure
            @. R_P                  = (P - (P_old - dt .* (Kd .* divV - 0.5 .* ((Pf - Pf_old) ./ dt)  - Kd .* Pf ./ ((1 - ϕ) .* λ))) ./ (1 + dt .* Kd ./ ((1.0 - ϕ) .* λ)));
            @. P                    = P - dτ_P .* R_P;
            # Effective Pressure
            @. P_eff                = P-Pf;
            # Velocity Residuals
            @. Rvx                  =   -(P[2:end,2:end-1]-P[1:end-1,2:end-1])/Δ.x + (τxx[2:end,2:end-1]-τxx[1:end-1,2:end-1])/Δ.x + (τxy[:,2:end]-τxy[:,1:end-1])/Δ.y
            @. Rvy                  =   -(P[2:end-1,2:end]-P[2:end-1,1:end-1])/Δ.y + (τyy[2:end-1,2:end]-τyy[2:end-1,1:end-1])/Δ.y + (τxy[2:end,:]-τxy[1:end-1,:])/Δ.x
            # Convergence accelerator
            @. dPdτ                = R_P - (1.0 - damp) .* dPdτ
            @. dVxdτ               = Rvx + (1.0 - damp) .* dVxdτ
            @. dVydτ               = Rvy + (1.0 - damp) .* dVydτ
            # Update Velocities
            @. Vx[2:end-1, 2:end-1] = Vx[2:end-1, 2:end-1] + dVxdτ * dτ
            @. Vy[2:end-1, 2:end-1] = Vy[2:end-1, 2:end-1] + dVydτ * dτ
            # Checks
            if iter==1 || mod(iter, 1000)==0
                errVx, errVy = norm(Rvx), norm(Rvy)
                if iter==1
                    errVx0, errVy0 = errVx, errVy 
                end
                @printf("iter = %06d:  Err. Vx = %2.1e --- Err. Vx = %2.1e\n", iter, errVx/errVx0, errVy/errVy0)
                if (max(errVx/errVx0, errVy/errVy0) < tol) break end
            end
        end
        # Viz
        p1 = plot(aspect_ratio=1, xlims=extrema(xv), ylims=extrema(yv))
        p1 = heatmap!(xv, yv, Pf'.*P_sc/1e9, c=:turbo, title="Pf")
        p1 = plot!(X_incl, Y_incl, c=:white, label=:none)
        p2 = plot(aspect_ratio=1, xlims=extrema(xv), ylims=extrema(yv))
        p2 = heatmap!(xv, yv, X_s', c=:turbo, title="X_s")
        p2 = plot!(X_incl, Y_incl, c=:white, label=:none)
        p3 = plot(aspect_ratio=1, xlims=extrema(xv), ylims=extrema(yv))
        p3 = heatmap!(xv, yv, ρ_s'.*ρ_db, c=:turbo, title="ρ_s")
        p3 = plot!(X_incl, Y_incl, c=:white, label=:none)
        p4 = plot(aspect_ratio=1, xlims=extrema(xv), ylims=extrema(yv))
        p4 = heatmap!(xv, yv, ϕ'.*ρ_db, c=:turbo, title="ρ_f")
        p4 = plot!(X_incl, Y_incl, c=:white, label=:none)
        display(plot(p1, p2, p3, p4))
    end
end

main(nt=10, tol=1e-3, max_its=50e3)