function cx = clamp(x, lo, hi)
    if x < lo
        cx = lo;
    elseif x > hi
        cx = hi;
    else
        cx = x;
    end
endfunction