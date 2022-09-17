include("preamble.jl")

function discrete_relax_factory(d, t, is_adaptive, xL, xU, support_increment, is_relax = false)
    h = is_adaptive ? 0.0 : support_increment
    dr = DBDisc.DiscretizeRelax(d, t, h = h, skip_step2 = false, print_relax_time = false, relax = is_relax)
    DBB.set!(dr, ConstantStateBounds(xL, xU))
    dr
end
