

using MAT

function getMatData()

    f = matopen("D1.mat")
    D = read(f, "D1_data")

    ys = zeros(5001,100)
    for i = 1:100
        ys[:,i] = D[i]["y"]
    end
    return ys
end
