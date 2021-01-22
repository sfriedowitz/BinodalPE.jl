#==============================================================================#
# Binodal plotting tools
#==============================================================================#

@recipe function f(bndl::BinodalData; 
    map = x -> [x[1], x[2]], 
    primary = true, 
    tie = false,
    bulk = false,
    interp = false, 
    composition = :auto,
    omega = nothing
)
    if composition in (:phi, :conc) && !(typeof(omega) <: AbstractVector)
        error("Must specify monomer sizes :omega if composition mapping turned on.")
    end
    if composition == :phi
        bndl = tophi(bndl, omega)
    elseif composition == :conc
        bndl = toconc(bndl, omega)
    end

    # Map the full data state
    bulk_data = hcat([map(state.bulk) for state in bndl]...)'

    sup = hcat([map(state.sup) for state in bndl]...)'
    dense = hcat([map(state.dense) for state in bndl]...)'

    bndl_data = vcat(sup, reverse(dense, dims = 1))
    npoints = size(bndl_data, 1)
    ndim = size(bndl_data, 2)

    legend --> false
    if ndim == 3
        x = bndl_data[:,1]
        y = bndl_data[:,2]
        z = bndl_data[:,3]
            
        seriestype --> :path
        @series begin
            primary := true
            x, y, z
        end

        # Tie line plots if necessary    
        if tie
            mid = Int(npoints/2)

            lx = x[1:mid]
            ly = y[1:mid]
            lz = z[1:mid]
            rx = x[end:-1:mid+1] 
            ry = y[end:-1:mid+1]
            rz = z[end:-1:mid+1]

            primary := false
            linestyle := :dash
            linealpha := 0.5
            linewidth := 1.5
            marker := false

            for i = 1:mid
                tl = transpose([[lx[i], ly[i], lz[i]] [rx[i], ry[i], rz[i]]])
                @series begin
                    seriestype := :line
                    tl[:,1], tl[:,2], tl[:,3]
                end
            end
        end

        # Plot the bulk lines if necessary
        if bulk
            @series begin
                primary := false
                seriestype := :scatter
                markercolor := :white
                markeralpha := 0.75
                markerstrokecolor := :black
                markerstrokewidth := 2
                markerstrokealpha := 1.0
                bulk_data[:,1], bulk_data[:,2]
            end
        end

    elseif ndim == 2
        x = bndl_data[:,1]
        y = bndl_data[:,2]

        if interp
            s = range(0, stop = 1, length = npoints)
            itp = Interpolations.scale(interpolate(hcat(x, y), (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), s, 1:2)

            sf = range(0, stop = 1, length = 250)
            xf, yf = [itp(s,1) for s in sf], [itp(s,2) for s in sf]

            # Plot of interpolated data series w/ fine resolution
            @series begin
                seriestype --> :path
                primary := true
                marker := false
                xf, yf
            end
        else   
            # Standard plot of lines connecting data points
            seriestype --> :path
            primary --> true
            @series begin
                x, y
            end
        end

        # Tie line plots if necessary    
        if tie
            mid = Int(npoints/2)

            lx = x[1:mid]
            ly = y[1:mid]
            rx = x[end:-1:mid+1] 
            ry = y[end:-1:mid+1]

            primary := false
            linestyle := :dash
            linealpha := 0.5
            linewidth := 1.5
            marker := false

            for i = 1:mid
                tl = transpose([[lx[i], ly[i]] [rx[i], ry[i]]])
                @series begin
                    seriestype := :line
                    tl[:,1], tl[:,2]
                end
            end
        end

        # Plot the bulk lines if necessary
        if bulk
            @series begin
                primary := false
                seriestype := :scatter
                markercolor := :white
                markeralpha := 0.75
                markerstrokecolor := :black
                markerstrokewidth := 2
                markerstrokealpha := 1.0
                bulk_data[:,1], bulk_data[:,2]
            end
        end

    end
end