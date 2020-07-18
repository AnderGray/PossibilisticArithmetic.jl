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
    
    j = range(0, 1, length = n);

    lefts = left.(mems)
    PyPlot.step(lefts, j, color = edgeCol, where = "pre");

    rights = right.(mems)
    PyPlot.step(rights, j, color = edgeCol, where = "pre");

    PyPlot.plot([x.Core.lo, x.Core.hi] , [1, 1], color = edgeCol);

    if fill
        ax.fill_between([lefts; reverse(rights)], zeros(2 * n), [j; reverse(j)], alpha=alpha, color =fillcol)
    end

    xticks(fontsize = fontsize); yticks(fontsize = fontsize)
    xlabel("Range",fontsize = fontsize); ylabel("α",fontsize=fontsize);


end