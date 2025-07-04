multi_step_refinement(x_intervals; npoints = nothing, gridsize = nothing) = multi_step_refinement(x_intervals, npoints, gridsize)
multi_step_refinement(::Any, ::Nothing, ::Nothing)      = error("Either `points` or `dx` must be provided.")
multi_step_refinement(x_intervals, npoints, ::Nothing)  = multi_step_refinement_from_npoints(x_intervals, npoints)
multi_step_refinement(x_intervals, ::Nothing, gridsize) = multi_step_refinement_from_gridsize(x_intervals, gridsize)
multi_step_refinement(x_intervals, npoints, gridsize)   = error("Both `points` or `dx` are provided. You can provide only one of them!")

function multi_step_refinement_from_npoints(x_intervals, n_intervals)
    xi = ntuple(Val(length(n_intervals))) do i
        collect(LinRange(x_intervals[i], x_intervals[i+1], n_intervals[i]))
    end
    x = merge_grid_intervals(xi)
    return x, abs.(diff(x))
end


function multi_step_refinement_from_gridsize(x_intervals, dx_intervals)
    xi = ntuple(Val(length(dx_intervals))) do i
        collect(x_intervals[i]:dx_intervals[i]:x_intervals[i+1])
    end
    x = merge_grid_intervals(xi)
    return x, abs.(diff(x))
end
