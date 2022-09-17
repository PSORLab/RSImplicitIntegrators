include("extended_preamble.jl")

x0(p::Vector{T}) where T = T[1.2; 1.1]
function RHS!(dx::Vector{T}, x, p, t) where T 
    dx[1] = p[1]*x[1]*(1. - x[2])
    dx[2] = p[2]*x[2]*(x[1] - 1.)
    nothing
end
function Jx_RHS!(out::VecOrMat{T}, x, p, t) where T
    out[1,1] = p[1]*(1. - x[2])
    out[1,2] = -p[1]*x[1]
    out[2,1] = p[2]*x[2]
    out[2,2] = p[2]*(x[1] - 1.)
    nothing
end
function Jp_RHS!(out::VecOrMat{T}, x, p, t) where T
    out[1,1] = x[1]*(1. - x[2])
    out[1,2] = zero(T)
    out[2,1] = zero(T)
    out[2,2] = x[2]*(x[1] - 1.)
    nothing
end

nx = 2
tL = 0.0
tU = 4.0
nt = 20 
nt_plot = 1
support_increment = (tU - tL)/nt

xL = Float64[0.1; 0.1]
xU = Float64[2.0; 2.0]

discrete_relax_factory(d, t, is_adaptive, is_relax = false) = discrete_relax_factory(d, t, is_adaptive, xL, xU, support_increment, is_relax)

function factory1(d)
    d = DBDisc.Wilhelm2019(d, DBDisc.AM2())
    DBB.set!(d, ConstantStateBounds(xL, xU))
    d.evaluate_interval = true
    return d
end
factory2(d) = discrete_relax_factory(d, DBDisc.AdamsMoulton(2), false, true)
factory3(d) = discrete_relax_factory(d, DBDisc.HermiteObreschkoff(3,3,false), false, true)

function factory4(d)
    d = DBDisc.Wilhelm2019(d, DBDisc.AM2())
    DBB.set!(d, ConstantStateBounds(xL, xU))
    d.evaluate_interval = true
    return d
end
factory5(d) = discrete_relax_factory(d, DBDisc.AdamsMoulton(2), true, true)
factory6(d) = discrete_relax_factory(d, DBDisc.HermiteObreschkoff(3,3,false), true, true)

get_step_count(integrator::DiscretizeRelax, nt) = integrator.step_count
get_step_count(integrator, nt) = nt + 1

function get_bounds(f, i, nt, is_adaptive = false)
    
    support_increment = (tU - tL)/nt
    tspan = (tL, tU)
    pL = [2.9; 0.9]
    pU = [3.1; 1.1]
    
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
    x0v = Float64[1.2; 1.1]
    tspan = (tL, tU)
    δt = (tU - tL)/n
    p1 = [2.9; 3.0; 3.1]
    p2 = [0.9; 1.0; 1.1]
    for p1i in p1
        for p2i in p2
            prob = ODEProblem(RHS!, x0v, tspan, [p1i; p2i])
            sol = solve(prob, Tsit5(), saveat=tL:δt:tU)
            uv = [sol.u[j][i] for j = 1:length(sol.u)]
            plot!(d, sol.t[26:100], uv[26:100], lc=:cyan, ls=:dash)
        end
    end
    return
end

function generate_plot3!()

    p1l = plot()
    plot_integrator_result!(p1l, factory1, 1, "1", :red, :dot, :diamond, 4, 5, 100)
    plot_integrator_result!(p1l, factory2, 1, "2", :blue,  :dash, :rect, 3, 5, 100)
    plot_integrator_result!(p1l, factory3, 1, "3", :green, :dashdot, :circle, 4, 5, 100)
    add_local_trajectories_to_plot!(p1l, 1, 100)

    xlabel!(p1l, L"t")
    ylabel!(p1l, L"x_1")

    p1r = plot()
    xlabel!(p1r, L"t")
    plot_integrator_result!(p1r, factory4, 1, "4", :red, :dot, :diamond, 4, 5, 100)
    plot_integrator_result!(p1r, factory5, 1, "5", :blue,  :dash, :rect, 3, 2, 56, true)
    plot_integrator_result!(p1r, factory6, 1, "6", :green, :dashdot, :circle, 4, 3, 45, true)
    add_local_trajectories_to_plot!(p1r, 1, 100)

    p2l = plot()
    println("Condition 1")
    plot_integrator_result!(p2l, factory1, 2, "1", :red, :dot, :diamond, 4, 5, 100)
    println("Condition 2")
    plot_integrator_result!(p2l, factory2, 2, "2", :blue,  :dash, :rect, 3, 5, 100)
    println("Condition 3")
    plot_integrator_result!(p2l, factory3, 2, "3", :green, :dashdot, :circle, 4, 5, 100)
    add_local_trajectories_to_plot!(p2l, 2, 100)
    xlabel!(p2l, L"t")
    ylabel!(p2l, L"x_2")

    p2r = plot()
    xlabel!(p2r, L"t")
    plot_integrator_result!(p2r, factory4, 2, "4", :red, :dot, :diamond, 4, 5, 100)
    println("Condition 4")
    plot_integrator_result!(p2r, factory5, 2, "5", :blue,  :dash, :rect, 3, 2, 56, true)
    println("Condition 5")
    plot_integrator_result!(p2r, factory6, 2, "6", :green, :dashdot, :circle, 4, 3, 45, true)
    add_local_trajectories_to_plot!(p2r, 2, 100)

    p2temp = plot()
    plot_integrator_result!(p2temp, factory3, 2, "3", :green, :dashdot, :circle, 4, 5, 100)
    
    p2 = plot(p1l, p1r, p2l, p2r, layout=grid(2,2, widths=(4/8,4/8)), size=(800,800), margin=5mm, legend = false, xtickfontsize=10, ytickfontsize=10, xguidefontsize=12, yguidefontsize=12)
    savefig(p2, "lokta_volterra_plot.pdf")
    return 
end

generate_plot3!()