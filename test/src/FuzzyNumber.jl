@testset "Conservative memberships" begin

    Nsamp = 10^4

    A = FuzzyNumber(0,1,2, steps = 1000)

    U1 = Uniform(0,1);
    U2 = Uniform(1,2);

    rand1 = rand(Nsamp)
    U1cdf = cdf(U1,rand1);

    rand2 = rand(Nsamp) .+1;
    U2cdf = 1 .- cdf(U2,rand2);

    @test all(membership.(A, rand1) .>= U1cdf)
    @test all(membership.(A, rand2) .>= U2cdf)


    A = FuzzyNumber(0,1,2, steps = 200)

    U1 = Uniform(0,1);
    U2 = Uniform(1,2);

    rand1 = rand(Nsamp)
    U1cdf = cdf(U1,rand1);

    rand2 = rand(Nsamp) .+1;
    U2cdf = 1 .- cdf(U2,rand2);

    @test all(membership.(A, rand1) .>= U1cdf)
    @test all(membership.(A, rand2) .>= U2cdf)


    A = FuzzyNumber(0,1,2, steps = 10)

    U1 = Uniform(0,1);
    U2 = Uniform(1,2);

    rand1 = rand(Nsamp)
    U1cdf = cdf(U1,rand1);

    rand2 = rand(Nsamp) .+1;
    U2cdf = 1 .- cdf(U2,rand2);

    @test all(membership.(A, rand1) .>= U1cdf)
    @test all(membership.(A, rand2) .>= U2cdf)

end

@testset "Mass tests" begin
    A = FuzzyNumber(1,2,3);

    @test mass(A,interval(1,3)) == interval(0,1)
    @test mass(A,interval(2)).hi == 1
    @test mass(A,interval(-10, -5)).hi == interval(0)

    A = FuzzyNumber(1,2,3, steps = 10);

    @test mass(A,interval(1,3)) == interval(0,1)
    @test mass(A,interval(2)).hi == 1
    @test mass(A,interval(-10, -5)).hi == interval(0)

end

@testset "Levelwise" begin

    A = FuzzyNumber(1,2,3);
    B = FuzzyNumber(2,3,4);

    @test interval(1,3) ⊆ A.Membership[1]
    @test interval(2) ⊆ A.Membership[end]

    C = levelwise(A,B)

    @test  interval(3,7) ⊆ C.Membership[1]
    @test interval(5) ⊆ C.Membership[end]

end
