using Pkg
Pkg.activate(".")
using Gridap
using GridapGmsh
using Gridap.Fields
using Gridap.CellData
using Gridap.Arrays
using LineSearches: BackTracking, Static, MoreThuente
using FillArrays

using Statistics
using JLD2

using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.CellData
"""
Unsteady Cylinder
add Gridap#generalized_alpha

Parabolic inlet value taken from "A cell-based smoothed finite element method stabilized by implicit SUPG/SPGP/Fractional step
method for incompressible flow", Mingyang, Guanjun, Huifen, Chen

VMS formulation based on bazilevs  
"""

#Parameters

Re = 1e6

St = 0.22 #approximately


ν = 0.000001 #0.001 m2/s 
ρ = 1 #1 kg/m3, not in the equations, just for Computing μ

D = 0.1 # [m] cylinder diameter
p0 = 0

order = 1 #Order of pressure and velocity
hf = VectorValue(0.0,0.0)

#umax = 1.5 #Velocity parameter for parabolic inlet; 
u0 = Re*ν/D

f = St*u0/D #[Hz] vortex shedding approximately frequence
dtmin = 1/f
initial_condition = false #print model of initial condition

#ODE settings #dt is 1/10 of the extimated period of the vortex shedding; total number of timestep 2000
#dt = ceil(dtmin/0.0001) * 0.0001/10 
#dt small to ensure CFL<1, the problem is more unstable at high reynolds. At low Re (1000) it tolerates a CFL>1, not at 1e6
dt = 0.0001

t0 = 0
tF = 2000*dt #Define the number of time step. Usually 5-10 for testing

Ntimestep = (tF-t0)/dt #Just to check how many time step





#MESH DEFINITION
model = GmshDiscreteModel("Cylinder2.msh")
writevtk(model,"model")

labels = get_face_labeling(model)



#u_free(x,t) = VectorValue(4*umax/0.41^2*(0.205^2 - x[2]^2),0) #parabolic inlet
u_free(x,t) = VectorValue(u0,0)
u_free(t::Real) = x -> u_free(x,t)


u_wall(x,t) = VectorValue(0.0, 0.0)
u_wall(t::Real) = x -> u_wall(x,t)

p_out(x,t) = p0
p_out(t::Real) = x -> p_out(x,t)






reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["Inlet","P_inlet","Limits","P_limits","Cylinder","P_cylinder"]) #P_ are the points 
reffeₚ = ReferenceFE(lagrangian,Float64, order)
Q = TestFESpace(model,reffeₚ, conformity=:H1, dirichlet_tags=["Outlet", "P_outlet"])

#TransientTrial because the condition are "time dependent"
U = TransientTrialFESpace(V, [u_free, u_free, u_wall, u_wall, u_wall, u_wall])
P = TransientTrialFESpace(Q, [p_out,p_out]) 


Y = MultiFieldFESpace([V, Q]) 
X = TransientMultiFieldFESpace([U, P])

degree = 4
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

#Integration Ω
Ωint(f) = ∫(f)dΩ


# Momentum residual, without viscous term
Rm(t,(u,p)) = ∂t(u) + u⋅∇(u) + ∇(p) #- ν*Δ(u) -hf

# Continuity residual
Rc(u) = ∇⋅u


Bᴳ(t,(u,p),(v,q)) = Ωint(∂t(u)⋅v) + Ωint((u⋅∇(u))⋅v) - Ωint((∇⋅v)*p) + Ωint(q*(∇⋅u)) + ν*Ωint(∇(v)⊙∇(u)) #Variational equations

B_SUPG(t,(u,p),(v,q)) =  Ωint((u⋅∇(v) + ∇(q))⋅(τm∘(u,G,GG).*Rm(t,(u,p)))) + Ωint((∇⋅v)⋅(τc∘(u,gg,G,GG) *Rc(u))) #SUPG terms 

B_VMS1(t,(u,p),(v,q)) = Ωint((u⋅∇(v)')⊙(τm∘(u,G,GG)*Rm(t,(u,p)))) #first VMS term

TRm(t,(u,p)) = τm∘(u,G,GG)*Rm(t,(u,p)) #Function used in the second VMS term

B_VMS2(t,(u,p),(v,q)) = -1*Ωint(∇(v)⊙(outer(TRm(t,(u,p)),TRm(t,(u,p)))))# second VMS term (To be added, it does not work at the moment)



Bᴹ(t,(u,p),(v,q)) = Bᴳ(t,(u,p),(v,q)) + B_SUPG(t,(u,p),(v,q)) + B_VMS1(t,(u,p),(v,q)) + B_VMS2(t,(u,p),(v,q))


h = lazy_map(h->h^(1/2),get_cell_measure(Ω))

ξₖ = get_cell_map(Ω)

#nn = get_node_coordinates(Ω)

Jt     = lazy_map(Broadcasting(∇),ξₖ)
inv_Jt = lazy_map(Operation(inv),Jt)

#d1 = evaluate(Jt, [Point(0.0,0.0)]) #x2

d = evaluate(inv_Jt, [Point(0.0,0.0)]) #x2

G = copy(d)
for i = 1:1:length(Jt)
q = (d[i] ⋅ d[i]') 
G[i] = q
end


GG = G .⊙ G

gg = zeros(num_cells(Ω))
for i  = 1:1:num_cells(Ω)
gg[i] = (d[i][1] +  d[i][3])^2 + (d[i][2] +  d[i][4])^2
end
gg

h = lazy_map(h->h^(1/2),get_cell_measure(Ω))

#CFL = u0*dt/minimum(h)
#Function for computing an element-wise CFL, to check approximately if <1
function compute_CFL(u,h)
  u = norm(u)
    return u*dt/h
  
end

#NO CFL correction, try with and without

function τm(uu,G,GG)
  Cᵢ = 1

  τ₃ = (ν^2 *GG)
  val(x) = x
  function val(x::Gridap.Fields.ForwardDiff.Dual)
    #print("Here\n")

    x.value
  end
  
  uu1 = val(uu[1])
  uu2 = val(uu[2])
  uu_new = VectorValue(uu1,uu2)
  

  if norm(uu_new) < 0.00001 #you can try if norm(u)<1e-5 or similar
      return (τ₃).^(-1/2)      
  end
  τ₁ = Cᵢ * (20/dt)^2 #add this 20 factor as STABILIZATION corrections, to be verified
  τ₂ = uu_new⋅G⋅uu_new
  return (τ₁ .+  τ₂ .+ τ₃).^(-1/2)     

end



function τc(uu,gg,G,GG)
  return 1/(τm(uu,G,GG)⋅ gg)
end


res(t,(u,p),(v,q)) = Bᴹ(t,(u,p),(v,q))



op = TransientFEOperator(res,X,Y)
nls = NLSolver(show_trace=true, method=:newton, linesearch=MoreThuente(), iterations=30)

solver = FESolver(nls)




U0 = U(0.0)
P0 = P(0.0)
X0 = X(0.0)

#Initial condition
uh0 = interpolate_everywhere(VectorValue(u0,0), U0)
ph0 = interpolate_everywhere(p_out(0), P0)
xh0 = interpolate_everywhere([uh0, ph0], X0)

#initial condition - derivative; It is not ideal, the first iteration are not perfect
vuh0 = interpolate_everywhere(VectorValue(0,0), U0)
vph0 = interpolate_everywhere(0, P0)
vxh0 = interpolate_everywhere([vuh0, vph0], X0)

#Check the initial CFL
CFL = compute_CFL∘(interpolate_everywhere(VectorValue(u0,0), U0),h)
writevtk(Ω, "CFL", cellfields=["CFL" => CFL])

#Compute Drag/Lift
Γ = BoundaryTriangulation(model; tags=["Cylinder","P_cylinder"]) 
dΓ = Measure(Γ,degree)
n_Γ = get_normal_vector(Γ) #Check if the direction is inside or outdside the domain; probably a -1 coefficient needed
F(p) = sum(∫(p⋅n_Γ)dΓ)

μ = ρ*ν
Friction(∇u) = μ*sum(∫(∇u⋅n_Γ)dΓ) #should be equivalent to ∂(u)/∂(n)

Friction(∇(uh0)) 


ρ∞ = 0.8 #ρ∞=1 no dissipation, ρ∞=0 max dissipation, ρ∞=0.5 quite good 
ode_solver = GeneralizedAlpha(nls,dt,ρ∞)
sol_t = solve(ode_solver,op,(xh0,vxh0),t0,tF)



_t_nn = t0
iteration = 1
L = F(ph0)[2] #1 Drag (x direction), 2 Lift (y direction)
F_cylinder = zeros(5,floor(Int, round(Ntimestep+1))) #time-step, Drag, Lift, Friction (x direction), Friction (y direction) 


F_cylinder[:, 1] = [_t_nn, F(ph0)[1], F(ph0)[2], Friction(∇(uh0))[1], Friction(∇(uh0))[2]]


createpvd("Cylinder_2d") do pvd
  for (xh_tn, tn) in sol_t
    global _t_nn
    _t_nn += dt
    global iteration
    iteration += 1
    println("it_num = $iteration\n")
    uh_tn = xh_tn[1]
    ph_tn = xh_tn[2]
    ωh_tn = ∇ × uh_tn
    global F_cylinder
    F_cylinder[:, iteration] = [_t_nn, F(ph_tn)[1], F(ph_tn)[2], Friction(∇(uh_tn))[1], Friction(∇(uh_tn))[2]]


    #if mod(iteration, 10)<1 #Useful when many time-step, reducing the number of file generated
      pvd[tn] = createvtk(Ω, "Results/Cylinder_2d_$_t_nn" * ".vtu", cellfields=["uh" => uh_tn, "ph" => ph_tn, "wh" => ωh_tn])
    #end
  end

end
@save "Results/Cylinder_Re_$Re.jld2" F_cylinder

