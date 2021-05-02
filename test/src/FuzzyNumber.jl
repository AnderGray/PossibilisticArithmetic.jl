@testset "Fuzzy Numbers" begin
    A = FuzzyNumber(1,2,3);
    B = FuzzyNumber(2,3,4);

    @test A.Membership[1] == interval(1,3)
    @test A.Membership[end] == interval(2)

    C = levelwise(A,B)

    @test C.Membership[1] == interval(3,7)
    @test C.Membership[end] == interval(5)

end
