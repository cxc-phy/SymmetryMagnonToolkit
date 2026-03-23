module Tool
export issubarray
export find_range_of_subarray

"""
issubarray(subarray::Vector{Int64}, array::Vector{Int64}) -> Bool
"""
function issubarray(subarray::Vector{Int64}, array::Vector{Int64})::Bool
    if length(array) < length(subarray) @goto no_result end
    idxs = findall(x->x==subarray[1], array)
    for idx in idxs
        if idx+length(subarray)-1>length(array) @goto no_result end
        if array[idx:idx+length(subarray)-1] == subarray return true end
    end
    @label no_result
    return false
end

"""
find_range_of_subarray(subarray::Vector{Int64}, array::Vector{Int64}) -> UnitRange{Int64}
"""
function find_range_of_subarray(subarray::Vector{Int64}, array::Vector{Int64})::UnitRange{Int64}
    if length(array) < length(subarray) @goto no_result end
    idxs = findall(x->x==subarray[1], array)
    for idx in idxs
        if idx+length(subarray)-1>length(array) @goto no_result end
        if array[idx:idx+length(subarray)-1] == subarray return idx:idx+length(subarray)-1 end
    end
    @label no_result
    return 0:-1
end

end #module Tool
