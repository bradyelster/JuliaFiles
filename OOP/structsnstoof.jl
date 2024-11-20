struct Person
    # info about person
    age::Int
    name::String

end

function UpperCamelCase(v::Vector{Float64})
    list1 = []
    for i in eachindex[v]
        list1 += i
    end
    return list1
end
