###
#   This file is a part of FuzzyArithmetic.jl package
#
#   Defines plotting functions
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
###


function plot(x :: FuzzyNumber, fill = true; name = missing, col = missing, alpha = 0.2, fontsize = 18)

    edgeCol = "red"; fillcol = "grey"
    if !(ismissing(col)); edgeCol = fillcol = col;end

    if ismissing(name); fig = figure(figsize=(10,10))else; fig = figure(name,figsize=(10,10));end
    ax = fig.add_subplot()

    mems = x.Membership; n = length(mems);
    
    j = range(0, 1, length = n+1);

    lefts = left.(mems)
    PyPlot.step([lefts; lefts[end]], j, color = edgeCol, where = "pre");

    rights = right.(mems)
    PyPlot.step([rights[1]; rights], j, color = edgeCol, where = "post");

    PyPlot.plot([x.Core.lo, x.Core.hi] , [1, 1], color = edgeCol);

    if fill
        ax.fill_between([lefts; reverse(rights)], zeros(2 * n), [j[2:end]; reverse(j[2:end])], alpha=alpha, color =fillcol)
    end

    xticks(fontsize = fontsize); yticks(fontsize = fontsize)
    xlabel("Range",fontsize = fontsize); ylabel("α",fontsize=fontsize);
end



function plot(x :: tnorm; name =missing, pn = 50, fontsize=18, alpha = 0.8)
    A = x.T
    if !ismissing(x.func)
        if x.func == zee; 
            for i in 1:5
                A[end-i,1:end-1] .= A[end,1:end-1]; A[1:end-1,end-i] .= A[1:end-1,end]; 
            end
        end
    end

    m = size(A)[1];
    if m < pn; ppn = m; else ppn = pn; end

    x = y = range(0, stop = 1,length=ppn)
    xgrid = repeat(x',ppn,1)
    ygrid = repeat(y,1,ppn)

    nm = round(m/ppn);

    z = A[1:Int(nm):end,1:Int(nm):end]

    if ismissing(name); fig = figure(figsize=(10,10)) else fig = figure(name,figsize=(10,10));end
    ax = fig.add_subplot(1,1,1,projection="3d")
    #ax = fig.add_subplot(2,1,1)
    plot_surface(xgrid, ygrid, z, rstride=2,edgecolors="k", cstride=2, alpha=alpha, linewidth=0.25, cmap=ColorMap("RdGy"))
    ax.view_init(45-27, 180+ 26)
    xlabel("y",fontsize=fontsize)
    ylabel("x",fontsize=fontsize)
    ax.zaxis.set_rotate_label(false);
    zlabel("T(x,y)", rotation = 0,fontsize=fontsize)
    #xticks(fontsize = fontsize); yticks(fontsize = fontsize);
    #PyPlot.title(title,fontsize=fontsize)
    tight_layout()
end

function plotContourCdf(x; name = "SurfacePlots",fontsize=18)
    A = x.T
    m = size(A)[1];
    if m < 200; ppn = m; else ppn = 200; end

    x = y = range(0,stop = 1,length=ppn)
    xgrid = repeat(x',ppn,1)
    ygrid = repeat(y,1,ppn)

    nm = round(m/ppn);

    z = A[1:Int(nm):end,1:Int(nm):end]
    fig = figure(title,figsize=(10,10))

    cp = contour(xgrid, ygrid, z,cmap=ColorMap("coolwarm"),levels = 20)
    xlabel("y",fontsize=fontsize)
    ylabel("x",fontsize=fontsize)

end



function plotTgen(name =missing)

    x = range(0,1,length = 200)

    Ts = [Tgen(xs,ys) for xs in x, ys in x]

    plot(tnorm(Ts),name = name)

end

function plotTInd(name = missing)

    x = range(0,1,length = 200)

    Ts = [Tind(xs,ys) for xs in x, ys in x]

    plot(tnorm(Ts), name = name)

end

function plotTest(name = missing)

    x = range(0,1,length = 200)

    Ts = TestNorm(x,x)

    plot(tnorm(Ts), name = name)

end