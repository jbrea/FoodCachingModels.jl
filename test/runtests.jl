using FoodCachingModels, Unitful
using Test

@testset "FoodCachingModels.jl" begin
    # Write your tests here.
end

@testset "Hunger" begin
    import FoodCachingModels: Hunger, update!, HungerModulatedCaching, duration_above_threshold
    using DifferentialEquations
    function hungerdyn!(du, u, p, t)
        if u[1] > 0.
            du[1] = -1.0/p.digestionduration
            du[2] = -u[2]/p.digestiontimeconstant
        else
            u[1] = 0.
            du[1] = 0.
            du[2] = (1. - u[2])/p.hungertimeconstant
        end
    end
    p = (hungertimeconstant = 6., digestiontimeconstant = 2., digestionduration = 4.)
    u0 = [3.2, .97]
    tspan = (0.0, 30.0)
    prob = ODEProblem(hungerdyn!, u0, tspan, p)
    sol = solve(prob, Tsit5(), abstol = 1e-12, reltol = 1e-6)
    h = Hunger(6.0u"minute", 2.0u"minute", 4.0u"minute", [.97], [3.2],
               HungerModulatedCaching())
    h0 = deepcopy(h)
    update!(h, 3u"minute")
    @test [h.stomach; h.hunger] ≈ sol(3.) atol = 1e-6
    update!(h, 15u"minute")
    @test [h.stomach; h.hunger] ≈ sol(18.) atol = 1e-6
    update!(h, 12u"minute")
    @test [h.stomach; h.hunger] ≈ sol(30) atol = 1e-6
    dur = duration_above_threshold(h0, 30.0u"minute", 0.9, false, 1)
    dur_approx = length(findall(u -> u[2] ≥ 0.9, sol(0:1e-4:30).u))*1e-4
    @test dur/1.0u"minute" ≈ dur_approx atol = 1e-3
end
