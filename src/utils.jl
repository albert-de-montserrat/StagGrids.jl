@inline function value_and_derivative(f::F, x::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    f_dual = f(ForwardDiff.Dual{T}(x, one(x)))
    dfdx = ForwardDiff.extract_derivative(T, f_dual)
    return f_dual.value, dfdx
end

function merge_grid_intervals(xi::NTuple{N, AbstractVector{T}}) where {N, T}
    x = xi[1]
    for i in 2:length(xi)
        x2 = xi[i]
        x[end] == x2[1] && popfirst!(x2)  # Remove the first element of the second array if it matches the last of the first
        x  = vcat(x, x2)    # Concatenate excluding the first element of the second array
    end
    return x
end

@generated function genereate_center_grid(xi::NTuple{N, AbstractVector{T}}) where {N, T}
    quote
        Base.@ntuple $N i-> begin
            x = xi[i]
            [(x[i] + x[i+1]) / 2 for i in 1:length(x)-1]
        end
    end
end