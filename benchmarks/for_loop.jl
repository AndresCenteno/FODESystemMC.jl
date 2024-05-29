using BenchmarkTools

# try a bunch of matrix multiplications with different type
mydict = Dict{String, Type}(
    "Float16" => Float16,
    "Float32" => Float32,
    "Float64" => Float64,
    "BigFloat" => BigFloat,
    "Int16" => Int16,
    "Int32" => Int32,
    "Int64" => Int64,
    "Int128" => Int128
)

function bench(d::Dict)
    for (key, value) in d
        println("Current type = $key")
        b = @benchmark rand($value, 10, 10) * rand($value, 10, 10)
        display(b)
    end
end

bench(mydict)