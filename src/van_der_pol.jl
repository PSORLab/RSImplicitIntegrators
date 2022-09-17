include("extended_preamble.jl")

x0(p::Vector{T}) where T = T[1.4; 2.3]
function RHS!(dx::Vector{T}, x::Vector{S}, p::Vector{Q}, t::R) where {T,S,Q,R} 
    dx[1] = x[2]
    dx[2] = p[1]*(1. - x[1]*x[1])*x[2] - x[1]
    nothing
end
function Jx_RHS!(out::VecOrMat{T}, x::Vector{S}, p::Vector{Q}, t::R) where {T,S,Q,R}
    out[1,1] = 0.0
    out[1,2] = 1.0
    out[2,1] = -2.0*p[1]*x[1]*x[2] - 1.0
    out[2,2] = p[1]*(1.0 - x[1]*x[1])
    nothing
end
function Jp_RHS!(out::VecOrMat{T}, x::Vector{S}, p::Vector{Q}, t::R) where {T,S,Q,R}
    out[1,1] = 0.0
    out[2,1] = (1. - x[1]*x[1])*x[2]
    nothing
end

nx = 2
tL = 0.0
tU = 3.0
nt = 200 
nt_plot = 1
support_increment = (tU - tL)/nt

xL = Float64[-3.0; -6.0]
xU = Float64[3.0; 4.0]

discrete_relax_factory(d, t, is_adaptive, is_relax = false) = discrete_relax_factory(d, t, is_adaptive, xL, xU, support_increment, is_relax)

function factory1(d)
    d = DBDisc.Wilhelm2019(d, DBDisc.AM2())
    DBB.set!(d, ConstantStateBounds(xL, xU))
    d.evaluate_interval = true
    return d
end
factory2(d) = discrete_relax_factory(d, DBDisc.AdamsMoulton(2), false, false)
factory3(d) = discrete_relax_factory(d, DBDisc.HermiteObreschkoff(3,3,false), false, false)

get_step_count(integrator::DiscretizeRelax, nt) = integrator.step_count
get_step_count(integrator, nt) = nt + 1

function get_bounds(f, i, nt, is_adaptive = false)
    
    support_increment = (tU - tL)/nt
    tspan = (tL, tU)
    pL = [1.30]
    pU = [1.40]
    
    prob = ODERelaxProb(RHS!, tspan, x0, pL, pU; Jx! = Jx_RHS!, Jp! = Jp_RHS!)

    if !is_adaptive
        support_t = Float64[i for i in 0.0:support_increment:tU]
        DBB.set!(prob, DBB.SupportSet{Float64}(support_t))
    end

    integrator = f(prob)
    @time DBB.relax!(integrator)
    step_count = get_step_count(integrator, nt)

    xL = zeros(nx, step_count)
    xU = zeros(nx, step_count)
    t = !is_adaptive ? support_t : integrator.time
    #@show length(integrator.time), length(support_t), is_adaptive

    DBB.getall!(xL, integrator, DBB.Bound{Lower}())
    DBB.getall!(xU, integrator, DBB.Bound{Upper}())

    integrator2 = f(prob)
    DBB.relax!(integrator2)

    return step_count, t, xL[i,:], xU[i,:]
end

function add_line_plots!(p, x, yL, yU, label, lc, ls, msc, ms)
    plot!(p, x, yL, label=label, lw = LINEWIDTH, lc=lc, ls=ls, markershape = msc, ms = ms, mc = lc)
    plot!(p, x, yU, label="", lw = LINEWIDTH, lc=lc, ls=ls, markershape = msc, ms = ms, mc = lc)
    return 
end

function plot_integrator_result!(p, f, i, l, c, s, msc, ms, nt_plot, nt, is_adaptive = false)
    step_count, t, xL, xU = get_bounds(f, i, nt, is_adaptive)
    tt = is_adaptive ? t[1:nt_plot:end] : t[1:nt_plot:end]
    xLt = is_adaptive ? xL[1:nt_plot:end] : xL[1:nt_plot:end]
    xUt = is_adaptive ? xU[1:nt_plot:end] : xU[1:nt_plot:end]
    add_line_plots!(p, tt, xLt, xUt, l, c, s, msc, ms)
end

function add_local_trajectories_to_plot!(d, i, n)
    p1 = [1.30; 1.32; 1.34; 1.36; 1.38; 1.4]
    tspan = (tL, tU)
    δt = (tU - tL)/n
    for p1i in p1
        x0v = Float64[1.4; 2.3]
        prob = ODEProblem(RHS!, x0v, tspan, [p1i])
        sol = solve(prob, Tsit5(), saveat=tL:δt:tU)
        uv = [sol.u[j][i] for j = 1:length(sol.u)]
        plot!(d, sol.t, uv, lc=:cyan, ls=:dash, lw = 2)
    end
    return
end

function generate_plot4!()

    p1l = plot()
    add_local_trajectories_to_plot!(p1l, 1, 100)
    println("generated plot 1")
    plot_integrator_result!(p1l, factory1, 1, "1", :red, :dot, :diamond, 4, 10, 200)
    println("generated plot 2")
    plot_integrator_result!(p1l, factory2, 1, "2", :blue,  :dot, :rect, 4, 4, 200)
    println("generated plot 3")
    plot_integrator_result!(p1l, factory3, 1, "3", :green, :dot, :circle, 3, 4, 200)
    xlabel!(p1l, L"t")
    ylabel!(p1l, L"x_1")

    p1r = plot()
    add_local_trajectories_to_plot!(p1r, 2, 100)
    println("generated plot 4")
    plot_integrator_result!(p1r, factory1, 2, "1", :red, :dot, :diamond, 4, 10, 200)
    println("generated plot 5")
    plot_integrator_result!(p1r, factory2, 2, "2", :blue,  :dot, :rect, 4, 4, 200)
    println("generated plot 6")
    plot_integrator_result!(p1r, factory3, 2, "3", :green, :dot, :circle, 3, 4, 200)
    xlabel!(p1r, L"t")
    ylabel!(p1r, L"x_2")

    p2 = plot(p1l, p1r, layout=grid(1,2, widths=(4/8,4/8)), size=(800,400), margin=3mm, legend = false, xtickfontsize=10, ytickfontsize=10, xguidefontsize=12, yguidefontsize=12)
    savefig(p2, "van_der_pol.pdf")
    return 
end

generate_plot4!()