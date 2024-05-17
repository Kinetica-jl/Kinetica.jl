using Test
using Kinetica

@testset "Profile Construction" begin
    static = Kinetica.StaticConditionProfile(10.0)
    @test static.value == 10.0

    nulldirect = Kinetica.NullDirectProfile(;
        X_start = 300.0, 
        t_end = 10.0)
    @test nulldirect.X_start == 300.0
    @test nulldirect.t_end == 10.0
    @test nulldirect.f(5.0, nulldirect) ≈ 300.0
    @test length(nulldirect.tstops) == 1
    @test nulldirect.tstops[1] ≈ 10.0

    lineardirect = LinearDirectProfile(;
        rate = 50.0,
        X_start = 300.0,
        X_end = 500.0)
    @test lineardirect.rate == 50.0
    @test lineardirect.X_start == 300.0
    @test lineardirect.X_end == 500.0
    @test lineardirect.t_end ≈ 4.0
    @test lineardirect.f(2.0, lineardirect) ≈ 400.0
    @test length(lineardirect.tstops) == 1
    @test lineardirect.tstops[1] ≈ 4.0

    nullgradient = Kinetica.NullGradientProfile(;
        X = 300.0,
        t_end = 10.0)
    @test nullgradient.X_start == 300.0
    @test nullgradient.t_end == 10.0
    @test nullgradient.grad(5.0, nullgradient) == 0.0
    @test length(nullgradient.tstops) == 1
    @test nullgradient.tstops[1] ≈ 10.0

    lineargradient = LinearGradientProfile(;
        rate = 50.0,
        X_start = 300.0,
        X_end = 500.0)
    @test lineargradient.rate == 50.0
    @test lineargradient.X_start == 300.0
    @test lineargradient.X_end == 500.0
    @test lineargradient.t_end ≈ 4.0
    @test lineargradient.grad(2.0, lineargradient) == 50.0
    @test lineargradient.grad(5.0, lineargradient) == 0.0
    @test length(lineargradient.tstops) == 1
    @test lineargradient.tstops[1] ≈ 4.0

    doublerampgradient = DoubleRampGradientProfile(;
        X_start = 300.0,
        t_start_plateau = 5.0,
        rate1 = 10.0,
        X_mid = 500.0,
        t_mid_plateau = 3.0,
        rate2 = -20.0,
        X_end = 200.0,
        t_end_plateau = 5.0)
    @test doublerampgradient.rate1 == 10.0
    @test doublerampgradient.rate2 == -20.0
    @test doublerampgradient.X_start == 300.0
    @test doublerampgradient.X_mid == 500.0
    @test doublerampgradient.X_end == 200.0
    @test doublerampgradient.t_start_plateau == 5.0
    @test doublerampgradient.t_mid_plateau == 3.0
    @test doublerampgradient.t_end_plateau == 5.0
    @test doublerampgradient.t_blend == 0.0
    @test doublerampgradient.t_end ≈ 48.0
    @test doublerampgradient.tstops ≈ [5.0, 25.0, 28.0, 43.0, 48.0]
    @test doublerampgradient.grad(1.0, doublerampgradient) == 0.0
    @test doublerampgradient.grad(15.0, doublerampgradient) == 10.0
    @test doublerampgradient.grad(27.0, doublerampgradient) == 0.0
    @test doublerampgradient.grad(35.0, doublerampgradient) == -20.0
    @test doublerampgradient.grad(45.0, doublerampgradient) == 0.0
    @test doublerampgradient.grad(100.0, doublerampgradient) == 0.0

    doublerampgradientblended = DoubleRampGradientProfile(;
        X_start = 300.0,
        t_start_plateau = 5.0,
        rate1 = 10.0,
        X_mid = 500.0,
        t_mid_plateau = 3.0,
        rate2 = -20.0,
        X_end = 200.0,
        t_end_plateau = 5.0,
        t_blend = 0.1)
    @test doublerampgradientblended.t_blend == 0.1
    @test doublerampgradientblended.tstops ≈ [4.9, 5.1, 24.9, 25.1, 27.9, 28.1, 42.9, 43.1, 48.0]
end

@testset "ConditionSet Construction" begin
    csc = ConditionSet(Dict(
        :T => LinearDirectProfile(;
            rate = 50.0,
            X_start = 300.0,
            X_end = 500.0),
        :P => DoubleRampGradientProfile(
            X_start = 1e5,
            t_start_plateau = 1.0,
            rate1 = 1e3,
            X_mid = 2e5,
            t_mid_plateau = 10.0,
            rate2 = -1e3,
            X_end = 1e5,
            t_end_plateau = 1.0,
            t_blend = 0.1),
        :V => 1e3
    ))
    @test issetequal(Set(csc.symbols), Set([:T, :P, :V]))
    @test length(csc.profiles) == 3
    @test csc.discrete_updates == false
    @test isnothing(csc.ts_update)

    csd = ConditionSet(Dict(
        :T => LinearDirectProfile(;
            rate = 50.0,
            X_start = 300.0,
            X_end = 500.0),
        :P => DoubleRampGradientProfile(
            X_start = 1e5,
            t_start_plateau = 1.0,
            rate1 = 1e3,
            X_mid = 2e5,
            t_mid_plateau = 10.0,
            rate2 = -1e3,
            X_end = 1e5,
            t_end_plateau = 1.0,
            t_blend = 0.1),
        :V => 1e3),
        ts_update = 1e-3)
    @test csd.discrete_updates == true
    @test csd.ts_update ≈ 1e-3

    @test_throws ArgumentError ConditionSet(Dict(:X => "abc"))
end