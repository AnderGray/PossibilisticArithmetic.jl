@testset "Conservative memberships" begin

    Nsamp = 10^4
    for i in [5, 10, 100, 200, 500, 1000]
        A = FuzzyNumber(0, 1, 2, steps = i)

        U1 = Uniform(0, 1);
        U2 = Uniform(1, 2);

        rand1 = rand(Nsamp)
        U1cdf = cdf.(U1, rand1);

        rand2 = rand(Nsamp) .+ 1;
        U2cdf = 1 .- cdf.(U2, rand2);

        @test all(membership.(A, rand1) .>= U1cdf)
        @test all(membership.(A, rand2) .>= U2cdf)
    end
end


@testset "Mass tests" begin

    types = [Float32, Float64, BigFloat]

    U1 = Uniform(1, 3)
    U2 = Uniform(1, 2)
    U3 = Uniform(2, 3)
    U4 = Uniform(1, 3)

    for t in types

        for j in [5, 10, 100, 200, 500]
            A = FuzzyNumber(1, 2, 3, steps = j);

            @test mass(A, interval(1, 3)).hi == 1
            @test mass(A, interval(2)) == interval(0,1)
            @test mass(A, interval(-10, -5)) == interval(0)

            @test check_inside(A, U1)
            @test check_inside(A, U2)
            @test check_inside(A, U3)
            @test check_inside(A, U4)

            @test check_inside(A, U1, usemass = true)
            @test check_inside(A, U2, usemass = true)
            @test check_inside(A, U3, usemass = true)
            @test check_inside(A, U4, usemass = true)
        end
    end

end

@testset "Levelwise" begin

    A = FuzzyNumber(1, 2, 3);
    B = FuzzyNumber(2, 3, 4);

    @test interval(1, 3) ⊆ A.Membership[1]
    @test interval(2) ⊆ A.Membership[end]

    C = levelwise(A, B)

    @test  interval(3, 7) ⊆ C.Membership[1]
    @test interval(5) ⊆ C.Membership[end]

end

@testset "Floating point error" begin

    types = [Float32, Float64, BigFloat]

    for i = 1:length(types)-1
        t1 = types[i]
        t2 = types[i+1]
        for j in [5, 10, 100, 200, 500]

            A = Fuzzy(t1(1), t1(2), t1(3), steps = j)
            B = Fuzzy(t2(1), t2(2), t2(3), steps = j)

            @test all(B.Membership .⊆ A.Membership)
        end
    end
end
