macro dxi(ex)
    return quote
        $(esc(dxi(ex)))
    end
end

macro dx(grid, I)
    return quote
        _dx($(esc(grid)), $(esc.(I)))
    end
end

macro dy(grid, I)
    return quote
        _dy($(esc(grid)), $(esc.(I)))
    end
end

macro dz(grid, I)
    return quote
        _dz($(esc(grid)), $(esc.(I)))
    end
end
