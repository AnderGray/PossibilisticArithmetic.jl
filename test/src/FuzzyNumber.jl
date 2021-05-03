@testset "Fuzzy Numbers" begin
    A = FuzzyNumber(1,2,3);
    B = FuzzyNumber(2,3,4);

    @test interval(1,3) ⊆ A.Membership[1]
    @test interval(2) ⊆ A.Membership[end]

    C = levelwise(A,B)

    @test  interval(3,7) ⊆ C.Membership[1]
    @test interval(5) ⊆ C.Membership[end]

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
