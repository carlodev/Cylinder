using Statistics
using JLD2
using Plots
using FFTW
using Gridap
Re = 1000
St = 0.22
@load "Cylinder_Re_1000.jld2"

n=100
idx = 400
idxp = 400+n

D = 0.1
ν=0.001
u0 = Re*ν/D

f = St*u0/D #[Hz] vortex shedding approximately frequence
dtmin = 1/f
dt = 0.005

dtmin/dt

t_cylinder = F_cylinder[1,idx:idxp]
D_cylinder = F_cylinder[2,idx:idxp]
L_cylinder = F_cylinder[3,idx:idxp]
Friction_cylinder = norm.(eachcol(F_cylinder[4:5,idx:idxp]))

plt = plot(title="Unsteady Cylinder Re = $Re", ylabel="Lift[N/m]", xlabel="time [s]")
plot!(t_cylinder, L_cylinder, linestyle=:dot, markershape=:hexagon, label=false)
savefig("L_unsteady.pdf")



plt = plot(title="Unsteady Cylinder Re = $Re", ylabel="Drag[N/m]", xlabel="time [s]")
plot!(t_cylinder, D_cylinder, linestyle=:dot, markershape=:hexagon,label=false)
savefig("D_unsteady.pdf")


#Spectral analysis

nt=length(t_cylinder)
fhat=fft(L_cylinder)

PSD=fhat.*conj(fhat)/(nt)
PSD = real(fftshift(PSD))
freqs = fftshift(fftfreq(nt,1/dt))

high_PSD_idx = findall(real(PSD) .>100) #threshold for selcting meaningful frequencies


plot(title="Unsteady Cylinder Re = $Re", ylabel="PSD - Lift signal", xlabel="f [Hz]",yaxis=:log)
plot!(freqs, PSD, label=false, xlim = [0, 0.5/dt])
for i in high_PSD_idx
    println("$i\n")
    i = Int(i)
    if freqs[i]>0
        fi = freqs[i]
        PSDi = PSD[i]
        println("$fi\n $PSDi\n")
        plot!([fi], [PSDi], seriestype = :scatter, markershape=:hexagon, markersize=5, label="$fi Hz")


    end
end


savefig("PSD_L.pdf")


fhat=fft(D_cylinder)
PSD=fhat.*conj(fhat)/(nt)
PSD = real(fftshift(PSD))
freqs = fftshift(fftfreq(nt,1/dt))

high_PSD_idx = findall(real(PSD) .>10) #threshold for selcting meaningful frequencies

plot(title="Unsteady Cylinder Re = $Re", ylabel="PSD - Drag signal", xlabel="f [Hz]", yaxis=:log)
plot!(freqs, PSD, label=false, xlim = [0, 0.5/dt])
for i in high_PSD_idx
    println("$i\n")
    i = Int(i)
    if freqs[i]>0
        fi = freqs[i]
        PSDi = PSD[i]
        println("$fi\n $PSDi\n")
        plot!([fi], [PSDi], seriestype = :scatter, markershape=:hexagon, markersize=5, label="$fi Hz")


    end
end
savefig("PSD_D.pdf")


