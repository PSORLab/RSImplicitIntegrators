include("preamble.jl")

nx = 1
tL = 0.0
tU = 1.0
nt = 100
nt_plot = 1
support_increment = (tU - tL)/nt

function discrete_relax_factory(d, t, is_adaptive, h, is_relax = false)
    dr = DBDisc.DiscretizeRelax(d, t, h = h, skip_step2 = false, print_relax_time = false, relax = is_relax)
    xL = [0.01]
    xU = [9.0]
    DBB.set!(dr, ConstantStateBounds(xL, xU))
    dr
end

function wilhelm_factory(d)
    d = DBDisc.Wilhelm2019(d, DBDisc.AM2())
    xL = [0.01]
    xU = [9.0]
    DBB.set!(d, ConstantStateBounds(xL, xU))
    d.evaluate_interval = true
    return d
end


pilms2_factory_100(d) = discrete_relax_factory(d, DBDisc.AdamsMoulton(2), false, 0.01)
pilms3_factory_100(d) = discrete_relax_factory(d, DBDisc.AdamsMoulton(3), false, 0.01)
pilms4_factory_100(d) = discrete_relax_factory(d, DBDisc.AdamsMoulton(4), false, 0.01)

ho11_factory_100(d) = discrete_relax_factory(d, DBDisc.HermiteObreschkoff(1,1,false), false, 0.01)
ho22_factory_100(d) = discrete_relax_factory(d, DBDisc.HermiteObreschkoff(2,2,false), false, 0.01)
ho33_factory_100(d) = discrete_relax_factory(d, DBDisc.HermiteObreschkoff(3,3,false), false, 0.01)

ho11_factory_30(d) = discrete_relax_factory(d, DBDisc.HermiteObreschkoff(1,1,false), false, 0.025)
ho22_factory_30(d) = discrete_relax_factory(d, DBDisc.HermiteObreschkoff(2,2,false), false, 0.025)
ho33_factory_30(d) = discrete_relax_factory(d, DBDisc.HermiteObreschkoff(3,3,false), false, 0.025)

ho11A_factory(d) = discrete_relax_factory(d, DBDisc.HermiteObreschkoff(1,1,false), true)
ho22A_factory(d) = discrete_relax_factory(d, DBDisc.HermiteObreschkoff(2,2,false), true)
ho33A_factory(d) = discrete_relax_factory(d, DBDisc.HermiteObreschkoff(3,3,false), true)

x0(p::Vector{T}) where T = T[9.0]
function RHS!(dx::Vector{T}, x, p, t) where T 
    dx[1] = -x[1]*x[1] + p[1]
    nothing
end
function Jx_RHS!(out::VecOrMat{T}, x, p, t) where T
    out[1,1] = -2.0*x[1]
    nothing
end
function Jp_RHS!(out::VecOrMat{T}, x, p, t) where T
    out[1,1] = one(p[1])
    nothing
end

function get_bounds(f, i, nt)
    
    support_increment = (tU - tL)/nt
    tspan = (tL, tU)
    pL = [-1.0]
    pU = [1.0]
    
    prob = ODERelaxProb(RHS!, tspan, x0, pL, pU; Jx! = Jx_RHS!, Jp! = Jp_RHS!)

    t = Float64[i for i in 0.0:support_increment:tU]
    DBB.set!(prob, DBB.SupportSet{Float64}(t))

    integrator = f(prob)

    println("----------------")
    println("----------------")
    println("BEGIN THE RELAX")
    println("----------------")
    println("----------------")

    DBB.relax!(integrator)

    println("----------------")
    println("----------------")
    println("END THE RELAX")
    println("----------------")
    println("----------------")

    xL = zeros(nx, nt+1)
    xU = zeros(nx, nt+1)


    @show size(t)
    @show size(xL)
    @show size(xU)

    DBB.getall!(xL, integrator, DBB.Bound{Lower}())
    DBB.getall!(xU, integrator, DBB.Bound{Upper}())

    integrator2 = f(prob)
    DBB.relax!(integrator2)

    return t, xL[i,:], xU[i,:]
end

function add_line_plots!(p, x, yL, yU, label, lc, ls, msc, ms)
    plot!(p, x, yL, label=label, lw = LINEWIDTH, lc=lc, ls=ls, markershape = msc, ms = ms, mc = lc)
    plot!(p, x, yU, label="", lw = LINEWIDTH, lc=lc, ls=ls, markershape = msc, ms = ms, mc = lc)
    return 
end

function plot_integrator_result!(p, f, i, l, c, s, msc, ms, nt_plot, nt)
    t, xL, xU = get_bounds(f, i, nt)
    tt = t[1:nt_plot:end]
    xLt = xL[1:nt_plot:end]
    xUt = xU[1:nt_plot:end]
    add_line_plots!(p, tt, xLt, xUt, l, c, s, msc, ms)
end

function generate_plot1!()

    p1l = plot()
    plot_integrator_result!(p1l, pilms2_factory_100, 1, "ORDER 2", :red, :dot, :diamond, 4, 5, 100)
    plot_integrator_result!(p1l, wilhelm_factory, 1, "NUMERICAL EXACT", :green, :dashdot, :circle, 4, 5, 100)

    x0v = Float64[9.0]
    tspan = (tL, tU)
    for p = -1:0.25:1
        prob = ODEProblem(RHS!, x0v, tspan, p)
        sol = solve(prob, Tsit5(), saveat=0.0:0.01:1.0)
        uv = [sol.u[i][1] for i = 1:length(sol.u)]
        plot!(p1l, sol.t, uv, lc=:cyan, ls=:dash) #label="SOLUTION", lw = LINEWIDTH, lc=:cyan, ls=:solid, markershape = :star, ms = :solid, mc = :cyan)
    end

    xlabel!(p1l, L"t")
    ylabel!(p1l, L"x")

    p1r = plot()
    plot_integrator_result!(p1r, pilms2_factory_100, 1, "ORDER 2", :red, :dot, :diamond, 4, 5, 100)
    plot_integrator_result!(p1r, pilms3_factory_100, 1, "ORDER 3", :blue,  :dash, :rect, 3, 5, 100)
    plot_integrator_result!(p1r, pilms4_factory_100, 1, "ORDER 4", :green, :dashdot, :circle, 4, 5, 100)

    x0v = Float64[9.0]
    tspan = (tL, tU)
    for p = -1:0.25:1
        prob = ODEProblem(RHS!, x0v, tspan, p)
        sol = solve(prob, Tsit5(), saveat=0.0:0.01:1.0)
        uv = [sol.u[i][1] for i = 1:length(sol.u)]
        plot!(p1r, sol.t, uv, lc=:cyan, ls=:dash) #label="SOLUTION", lw = LINEWIDTH, lc=:cyan, ls=:solid, markershape = :star, ms = :solid, mc = :cyan)
    end

    xlabel!(p1r, L"t")

    p1 = plot(p1l, p1r, layout=grid(1,2, widths=(4/8,4/8)), size=(800,500), margin=5mm, legend = false, xtickfontsize=10, ytickfontsize=10, xguidefontsize=12, yguidefontsize=12)
    savefig(p1, "exact_bounds_pilms_plot.pdf")
    return 
end

function generate_plot2!()

    p1l = plot()
    plot_integrator_result!(p1l, ho11_factory_100, 1, "HO11", :red, :dot, :diamond, 4, 5, 100)
    plot_integrator_result!(p1l, ho22_factory_100, 1, "HO22", :blue,  :dash, :rect, 3, 5, 100)
    plot_integrator_result!(p1l, ho33_factory_100, 1, "HO33", :green, :dashdot, :circle, 4, 5, 100)
    plot_integrator_result!(p1l, pilms2_factory_100, 1, "ORDER 2 PILMS", :purple, :solid, :diamond, 3, 5, 100)
    plot_integrator_result!(p1l, wilhelm_factory, 1, "NUMERICAL EXACT", :orange, :solid, :circle, 3, 5, 100)
    
    x0v = Float64[9.0]
    tspan = (tL, tU)
    for p = -1:0.25:1
        prob = ODEProblem(RHS!, x0v, tspan, p)
        sol = solve(prob, Tsit5(), saveat=0.0:0.01:1.0)
        uv = [sol.u[i][1] for i = 1:length(sol.u)]
        plot!(p1l, sol.t, uv, lc=:cyan, ls=:dash) #label="SOLUTION", lw = LINEWIDTH, lc=:cyan, ls=:solid, markershape = :star, ms = :solid, mc = :cyan)
    end
    #p2l = plot!(p1l, pt[2], pt[3], pt[4], pt[5], pt[6], pt[7])

    xlabel!(p1l, L"t")
    ylabel!(p1l, L"x")

    p1r = plot()
    xlabel!(p1r, L"t")
    plot_integrator_result!(p1r, ho11_factory_30, 1, "HO11", :red, :dot, :diamond, 4, 1, 40)
    plot_integrator_result!(p1r, ho22_factory_30, 1, "HO22", :blue,  :dash, :rect, 3, 1, 40)
    plot_integrator_result!(p1r, ho33_factory_30, 1, "HO33", :green, :dashdot, :circle, 4, 1, 40)

    for p = -1:0.25:1
        prob = ODEProblem(RHS!, x0v, tspan, p)
        sol = solve(prob, Tsit5(), saveat=0.0:0.04:1.0)
        uv = [sol.u[i][1] for i = 1:length(sol.u)]
        plot!(p1r, sol.t, uv, lc=:cyan, ls=:dash) #label="SOLUTION", lw = LINEWIDTH, lc=:cyan, ls=:solid, markershape = :star, ms = :solid, mc = :cyan)
    end
   
    p2 = plot(p1l, p1r, layout=grid(1,2, widths=(4/8,4/8)), size=(800,500), margin=5mm, legend = false, xtickfontsize=10, ytickfontsize=10, xguidefontsize=12, yguidefontsize=12)
    savefig(p2, "ho_bounds_plot.pdf")
    return 
end

generate_plot1!()
generate_plot2!()