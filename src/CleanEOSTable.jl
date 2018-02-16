module CleanEOSTable

export cleanEOSTable

using HDF5



function clip(range::Range{Int}, mask::Range{Int})::Range{Int}
    imin = max(start(range), start(mask))
    imax = min(last(range), last(mask))
    imin:imax
end

function linterp{T<:AbstractFloat}(x0::T, y0::T, x1::T, y1::T, x::T,
                                   extrap::Bool=false)::T
    if !extrap
        if x0 < x1
            if !(x0 <= x && x <= x1)
                @show x0, x1, x
            end
            @assert x0 <= x && x <= x1
        elseif x0 > x1
            @assert x1 <= x && x <= x0
        else
            # We expect the input to have a non-empty range
            @assert false
            @assert x == x0
        end
    end
    y = (x1 - x) / (x1 - x0) * y0 + (x - x0) / (x1 - x0) * y1
    if !extrap
        if y0 < y1
            @assert y0 <= y && y <= y1
        elseif y0 > y1
            @assert y1 <= y && y <= y0
        else
            @assert y == y0
        end
    end
    y
end

function linterp(x0::Number, y0::Number, x1::Number, y1::Number, x::Number,
                 extrap::Bool=false)::Number
    linterp(promote(x0, y0, x1, y1, x)..., extrap)
end

function linterp0{T<:AbstractFloat}(arr::Array{T,3}, i::Int, j::Int, k::Int)::T
    arr[i,j,k]
end
function linterp1{T<:AbstractFloat}(arr::Array{T,3}, i::Int, j::Int, k::Int,
                                    x::T)::T
    (1 - x) * linterp0(arr, i, j, k) + x * linterp0(arr, i+1, j, k)
end
function linterp2{T<:AbstractFloat}(arr::Array{T,3}, i::Int, j::Int, k::Int,
                                    x::T, y::T)::T
    (1 - y) * linterp1(arr, i, j, k, x) + y * linterp1(arr, i, j+1, k, x)
end
function linterp3{T<:AbstractFloat}(arr::Array{T,3}, i::Int, j::Int, k::Int,
                                    x::T, y::T, z::T)::T
    @assert x>=0 && x<=1
    @assert y>=0 && y<=1
    @assert z>=0 && z<=1
    (1 - z) * linterp2(arr, i, j, k, x, y) + z * linterp2(arr, i, j, k+1, x, y)
end

"""
Flag inconsistent table entries:
The internal energy needs to be a strictly monotonically increasing function
of the temperature
"""
function flag_inconsistencies{T<:AbstractFloat}(table_logenergy::Array{T, 3})::
    BitArray{3}
    println("Flagging inconsistendies in table...")
    (pointsrho, pointstemp, pointsye) = size(table_logenergy)
    outliers = falses(pointsrho, pointstemp, pointsye)
    for k in 1:pointsye, i in 1:pointsrho
        prev_j = 1
        j = 2
        while j <= pointstemp
            @assert !outliers[i,prev_j,k]
            if table_logenergy[i,j,k] <= table_logenergy[i,prev_j,k]
                # We found an inconsistency. Either this point j or the
                # previous point prev_j is an outlier that needs to be
                # removed.
                irange = clip(i-1:i+1, 1:pointsrho)
                jrange = clip(prev_j-1:j+1, 1:pointstemp)
                krange = clip(k-1:k+1, 1:pointsye)
                avg = mean(table_logenergy[irange, jrange, krange])
                # Assume that the point further from the average is the
                # outlier
                if (abs(table_logenergy[i,prev_j,k] - avg)
                    > abs(table_logenergy[i,j,k] - avg))
                    # Disable the previous point
                    outliers[i,prev_j,k] = true
                    # Next, compare the current point to the one before
                    # the previous point
                    prev_j = prev_j - 1
                    while prev_j >= 1 && outliers[i,prev_j,k]
                        prev_j = prev_j - 1
                    end
                    if prev_j >= 1
                        # keep j
                    else
                        # All points before this are disabled: move on to
                        # the next point
                        prev_j = j
                        j = j + 1
                    end
                else
                    # Disable the current point
                    outliers[i,j,k] = true
                    # keep prev_j
                    j = j + 1
                end
            else
                # Move on to the next point pair
                prev_j = j
                j = j + 1
            end
        end
        # Find the first non-disabled point
        prev_j = 1
        while prev_j <= pointstemp && outliers[i,prev_j,k]
            prev_j = prev_j + 1
        end
        # Check consistency for all following points
        for j = prev_j+1:pointstemp
            if !outliers[i,j,k]
                @assert table_logenergy[i,j,k] > table_logenergy[i,prev_j,k]
                prev_j = j
            end
        end
    end
    noutliers = count(outliers)
    npoints = length(table_logenergy)
    noutliers_percent = round(100 * noutliers / npoints, 2)
    println("    Detected $noutliers outliers ($noutliers_percent%)")
    # println("    Outliers:")
    # n = 0
    # for k in 1:pointsye, j in 2:pointstemp, i in 1:pointsrho
    #     if outlier[i,j,k]
    #         n = n + 1
    #         println("        #$n: [$i,$j,$k]")
    #     end
    # end
    return outliers
end

# Correct outliers
function correct_table{T<:AbstractFloat}(old_table_logenergy::Array{T, 3},
                                         outliers::BitArray{3})::Array{T, 3}
    println("Correcting outliers...")
    @assert size(old_table_logenergy) == size(outliers)
    (pointsrho, pointstemp, pointsye) = size(old_table_logenergy)
    table_logenergy = copy(old_table_logenergy)
    for k in 1:pointsye, j in 1:pointstemp, i in 1:pointsrho
        if outliers[i,j,k]
            # Find two valid neighbouring points
            extrap = false
            prev_j = j - 1
            while prev_j >= 1 && outliers[i,prev_j,k]
                prev_j = prev_j - 1
            end
            next_j = j + 1
            while next_j <= pointstemp && outliers[i,next_j,k]
                next_j = next_j + 1
            end
            if prev_j < 1
                extrap = true
                prev_j = next_j
                next_j = next_j + 1
                while next_j <= pointstemp && outliers[i,next_j,k]
                    next_j = next_j + 1
                end
            end
            if next_j > pointstemp
                extrap = true
                next_j = prev_j
                prev_j = prev_j - 1
                while prev_j >= 1 && outliers[i,prev_j,k]
                    prev_j = prev_j - 1
                end
            end
            @assert prev_j >= 1 && prev_j <= pointstemp
            @assert next_j >= 1 && next_j <= pointstemp
            # Linear interpolation
            avg = linterp(prev_j, table_logenergy[i,prev_j,k],
                          next_j, table_logenergy[i,next_j,k], j, extrap)
            table_logenergy[i,j,k] = avg
        end
    end
    return table_logenergy
end

# Check consistency
function check_consistency{T<:AbstractFloat}(table_logenergy::Array{T, 3})::Void
    println("Checking consistency again...")
    (pointsrho, pointstemp, pointsye) = size(table_logenergy)
    for k in 1:pointsye, j in 2:pointstemp, i in 1:pointsrho
        @assert table_logenergy[i,j,k] > table_logenergy[i,j-1,k]
    end
end



function cleanEOSTable(tablename::AbstractString,
                       new_tablename::AbstractString)
    println("CleanEOSTable")

    println("Reading EOS table \"$tablename\"...")
    table = h5open(tablename, "r")

    # The 3D tables have three indices (rho, temp, ye), in this order
    # (in Fortran, Cactus, and Julia)
    pointsrho = read(table, "pointsrho")[1]
    pointstemp = read(table, "pointstemp")[1]
    pointsye = read(table, "pointsye")[1]

    table_logrho = read(table, "logrho")
    @assert size(table_logrho) == (pointsrho,)
    table_logtemp = read(table, "logtemp")
    @assert size(table_logtemp) == (pointstemp,)
    table_ye = read(table, "ye")
    @assert size(table_ye) == (pointsye,)

    table_logenergy = read(table, "logenergy")
    @assert size(table_logenergy) == (pointsrho, pointstemp, pointsye)
    npoints = length(table_logenergy)
    println("    Read table with ($pointsrho,$pointstemp,$pointsye) = $npoints data points")

    outliers = flag_inconsistencies(table_logenergy)
    new_table_logenergy = correct_table(table_logenergy, outliers)
    check_consistency(new_table_logenergy)

    println("Writing corrected table...")
    new_table = h5open(new_tablename, "w")
    write(new_table, "corrected_logenergy", new_table_logenergy)
    close(new_table)

    println("Done.")
end

end # module
