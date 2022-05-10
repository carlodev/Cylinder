using Gridap
using GridapGmsh
using Gridap.Fields
using Gridap.CellData
using Gridap.Arrays
using LineSearches: BackTracking, Static, MoreThuente
using FillArrays

"""
LidDrivenCavityFlow
"""

#Parameters

Re = 1000
St = 0.22


ν = 0.001 #0.1 m2/s 
D = 0.1 #2 [m] cylinder diameter
p0 = 0

order = 1 #Order of pressure and velocity
hf = VectorValue(0.0,0.0)

umax = 1.5
u0 = Re*ν/D
f = St*u0/D #[Hz] vortex shedding approximately frequence
dtmin = 1/f
initial_condition = true #print model of initial condition

#ODE settings
dt = 0.01 
t0 = 0
tF = 1000*dt

Ntimestep = (tF-t0)/dt





#MESH DEFINITION
model = GmshDiscreteModel("Cylinder2.msh")
writevtk(model,"model")

labels = get_face_labeling(model)



#u_free(x,t) = VectorValue(4*umax/0.41^2*(0.205^2 - x[2]^2),0)
u_free(x,t) = VectorValue(u0,0)
u_free(t::Real) = x -> u_free(x,t)


u_wall(x,t) = VectorValue(0.0, 0.0)
u_wall(t::Real) = x -> u_wall(x,t)

p_out(x,t) = p0
p_out(t::Real) = x -> p_out(x,t)






reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["Inlet","P_inlet","Limits","P_limits","Cylinder","P_cylinder"])
reffeₚ = ReferenceFE(lagrangian,Float64, order)
#reffeₚ = ReferenceFE(lagrangian,Float64,order-1; space=:P)
#reffeₚ = ReferenceFE(lagrangian, Float64, order - 1)
#Q = TestFESpace(model,reffeₚ, conformity=:L2, constraint=:zeromean)
#Q = TestFESpace(model,reffeₚ, conformity=:L2, dirichlet_tags="interior")
Q = TestFESpace(model,reffeₚ, conformity=:H1, dirichlet_tags=["Outlet", "P_outlet"])

#Since we impose Dirichlet boundary conditions on the entire boundary ∂Ω, the mean value of the pressure is constrained to zero in order have a well posed problem
#Q = TestFESpace(model, reffeₚ)


#Transient is just for the fact that the boundary conditions change with time
U = TransientTrialFESpace(V, [u_free, u_free, u_wall, u_wall, u_wall, u_wall])
#P = TrialFESpace(Q) #?transient
P = TransientTrialFESpace(Q, [p_out,p_out]) #?transient



Y = MultiFieldFESpace([V, Q]) #?transient
X = TransientMultiFieldFESpace([U, P])

degree = 4
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)


h = lazy_map(h->h^(1/2),get_cell_measure(Ω))


# Momentum residual, without the viscous term
Rm(t,(u,p)) = ∂t(u) + u⋅∇(u) + ∇(p) - hf

# Continuity residual
Rc(u) = ∇⋅u


function τ(u,h)
  

    τ₂ = h^2/(4*ν)
    val(x) = x
    val(x::Gridap.Fields.ForwardDiff.Dual) = x.value
    u = val(norm(u))
    
    if iszero(u)
        return τ₂
        
    end
    τ₃ =  dt/2 #h/(2*u) #0  dt/2 #

    τ₁ = h/(2*u) #h/(2*u) #
    return 1/(1/τ₁ + 1/τ₂ + 1/τ₃)
    
end


τb(u,h) = (u⋅u)*τ(u,h)

var_equations(t,(u,p),(v,q)) = ∫(
    ν*∇(v)⊙∇(u) # Viscous term
    + v⊙Rm(t,(u,p)) # Other momentum terms
    + q*Rc(u)
 )dΩ # Continuity


stab_equations(t,(u,p),(v,q)) = ∫(  (τ∘(u,h)*(u⋅∇(v) + ∇(q)))⊙Rm(t,(u,p)) # First term: SUPG, second term: PSPG
    +τb∘(u,h)*(∇⋅v)⊙Rc(u) # Bulk viscosity. Try commenting out both stabilization terms to see what happens in periodic and non-periodic cases
)dΩ


res(t,(u,p),(v,q)) = var_equations(t,(u,p),(v,q)) + stab_equations(t,(u,p),(v,q))


op = TransientFEOperator(res,X,Y)
nls = NLSolver(show_trace=true, method=:newton, linesearch=MoreThuente(), iterations=30)

solver = FESolver(nls)




U0 = U(0.0)
P0 = P(0.0)
X0 = X(0.0)

#uh0 = interpolate_everywhere(u_start00, U0)
uh0 = interpolate_everywhere(VectorValue(u0,0), U0)
ph0 = interpolate_everywhere(p_out(0), P0)

xh0 = interpolate_everywhere([uh0, ph0], X0)

vuh0 = interpolate_everywhere(VectorValue(0,0), U0)
vph0 = interpolate_everywhere(0, P0)


vxh0 = interpolate_everywhere([vuh0, vph0], X0)


ρ∞ = 0.5
ode_solver = GeneralizedAlpha(nls,dt,ρ∞)
sol_t = solve(ode_solver,op,(xh0,vxh0),t0,tF)



_t_nn = t0
iteration = 0
createpvd("TV_2d") do pvd
  for (xh_tn, tn) in sol_t
    global _t_nn
    _t_nn += dt
    global iteration
    iteration += 1
    println("it_num = $iteration\n")
    uh_tn = xh_tn[1]
    ph_tn = xh_tn[2]
    ωh_tn = ∇ × uh_tn
  
    #if mod(iteration, 10)<1
      pvd[tn] = createvtk(Ω, "Results/TV_2d_$_t_nn" * ".vtu", cellfields=["uh" => uh_tn, "ph" => ph_tn, "wh" => ωh_tn])
    #end
  end

end

