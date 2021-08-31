
#QQQQ - Make a copy for the heap

"""
QQQQ
"""
==(h1::MutableBinaryMaxHeap, h2::MutableBinaryMaxHeap)::Bool = extract_all!(deepcopy(h1)) == extract_all!(deepcopy(h2))

"""
QQQQ
"""
function map!(f, h::MutableBinaryMaxHeap)
    for n in h.nodes
        update!(h, n.handle, f(n.value))
    end
    return h
end

"""
QQQQ
"""
map(f,h::MutableBinaryMaxHeap) = map!(f,deepcopy(h))

#QQQQ - make the iterate

#QQQQ - T