using Dates 
mutable struct RADOLAN
    nodata::Float64
    producttype::String
    datetime::DateTime
    radarId::Int64
    datasize::Union{Int32,Missing}
    maxrange::Union{String,Missing}
    radolanversion::Union{String,Missing}
    precision::Union{Float32,Missing}
    intervalseconds::Union{Int32,Missing}
    intervalunit::Union{Int32,Missing}
    nrow::Union{Int32,Missing}
    ncol::Union{Int32,Missing}
    radarlocations::Union{Vector{String},Missing}
    reanalysisversion::Union{String,Missing}
    nlevel::Union{Int32,Missing}
    level::Union{Float64,Missing}
    radardays::Union{Vector{String},Missing}
    indicator::Union{String,Missing}
    imagecount::Union{Int32,Missing}
    predictiontime::Union{Int32,Missing}
    moduleflag::Union{Int32,Missing}
    quantification::Union{Int32,Missing}
    data::Union{Matrix{Union{Float64,Missing}},Missing}
    nodatamask::Union{BitArray,Missing}
    secondary::Union{Matrix{Bool},Missing}
    cluttermask::Union{Matrix{Bool},Missing}
    RADOLAN() = new(-9999, "rw", DateTime(2000), 10000, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing)

end
function print_info(obj::RADOLAN)
    @show obj.producttype
    @show obj.datetime
    @show obj.radarId
    @show obj.datasize
    @show obj.maxrange
    @show obj.radolanversion
    @show obj.precision
    @show obj.intervalseconds
    @show obj.intervalunit
    @show obj.nrow
    @show obj.ncol
    @show obj.radarlocations
    @show obj.reanalysisversion
    @show obj.nlevel
    @show obj.level
    @show obj.radardays
    @show obj.indicator
    @show obj.imagecount
    @show obj.predictiontime
    @show obj.moduleflag
    @show obj.quantification
end
function find_pos(x, string)
    pos=findfirst(x,string)
    if !isnothing(pos)
        return pos[1]
    else 
        return 0 
    end
end
function stop_pos(start_pos,all_pos, header_len)
    pos_gt_start=collect(Iterators.filter(x->x>start_pos,all_pos))
    minimum(push!(pos_gt_start,header_len))
end
function read_header_to_dict(io)
    rad_obj = RADOLAN()
    header_end = findfirst(io .== b"\x03")
    header_string = String(io[1:header_end-1])
    rad_obj.producttype = header_string[1:2]
    rad_obj.datetime = DateTime(header_string[3:8] * header_string[14:15] * "20" * header_string[16:17] * "00", "ddHHMMmmyyyySS") # watchout for 19xx as year
    token_position = Dict{String,Int}("BY" => 0, "VS" => 0, "SW" => 0, "PR" => 0, "INT" => 0, "GP" => 0, "MS" => 0, "LV" => 0, "CS" => 0,
        "MX" => 0, "BG" => 0, "ST" => 0, "VV" => 0, "MF" => 0, "QN" => 0, "VR" => 0, "U" => 0)
    for x in keys(token_position)
        token_position[x] = find_pos(x, header_string)
    end

    all_pos = values(token_position)
    header_len = length(header_string)
    for (token, position) in token_position
        if position == 0
            continue
        end
        stop = stop_pos(position, all_pos, header_len) - 1
        @inbounds rel = @view header_string[position+length(token):stop]

        if token == "BY"
            rad_obj.datasize = parse(Int32, rel) - header_len - 1
        elseif token == "VS"
            vs = rel
            text::String = ""
            if vs == 0
                text = "100 km and 128 km (mixed)"
            elseif vs == 1
                text = "100 km"
            elseif vs == 2
                text = "128 km"
            elseif vs == 3
                text = "150 km"
            end
            rad_obj.maxrange = text
        elseif token == "SW"
            rad_obj.radolanversion = strip(rel)
        elseif token == "PR"
            rad_obj.precision = parse(Float32, "1" * strip(rel))
        elseif token == "INT"
            rad_obj.intervalseconds = parse(Int32, rel) * 60
        elseif token == "GP"
            row_col = split(strip(rel), "x")
            rad_obj.nrow = parse(Int32, row_col[1])
            rad_obj.ncol = parse(Int32, row_col[2])
        elseif token == "MS"
            rad_obj.radarlocations = split(split(split(strip(rel), "<")[2], ">")[1], ",")
        elseif token == "LV"
            rad_obj.nlevel = parse(Int32, rel[1])
            rad_obj.level = parse(Float32, rel[2:end])
        elseif token == "CS"
            text = ""
            if rel == 0
                text = "near ground level"
            elseif rel == 1
                text = "maximum"
            elseif rel == 2
                text = "tops"
            end
            rad_obj.indicator = text
        elseif token == "MX"
            rad_obj.imagecount = parse(Int32, rel)
        elseif token == "BG"
            rad_obj.nrow = rel[1:Int(length(rel) // 2)]
            rad_obj.ncol = rel[Int(length(rel) // 2)+1:end]
        elseif token == "ST"
            rad_obj.radardays = split(split(split(strip(rel), "<")[2], ">")[1], ",")
        elseif token == "VV"
            rad_obj.predictiontime = parse(Int32, rel)
        elseif token == "MF"
            rad_obj.moduleflag = parse(Int32, rel)
        elseif token == "QN"
            rad_obj.quantification = parse(Int32, rel)
        elseif token == "VR"
            rad_obj.reanalysisversion = strip(rel)
        elseif token == "U"
            rad_obj.intervalunit = parse(Int32, rel)
        end
    end
    if (!ismissing(rad_obj.intervalunit)) && (rad_obj.intervalunit == 1)
        if !ismissing(rad_obj.intervalseconds)
            rad_obj.intervalseconds *= 1440
        end
    end
    if (!ismissing(rad_obj.nrow)) && (!ismissing(rad_obj.ncol))
        rad_obj.data = fill(rad_obj.nodata, (rad_obj.nrow, rad_obj.ncol))
    end
    return rad_obj, @view io[header_end+1:end]
end


function read_array(radolan::RADOLAN, io, missing_val)
    dat = @view io[1:radolan.datasize]
    
    if radolan.producttype ∈ ["RX","WX"]
        data_tmp_rx::Vector{UInt8} = reinterpret(UInt8, dat)
        radolan.nodatamask = reshape(data_tmp_rx .==250, (radolan.ncol, radolan.nrow))
        radolan.cluttermask = reshape(data_tmp_rx .==249, (radolan.ncol, radolan.nrow))
        radolan.data = reshape(data_tmp_rx, (radolan.ncol, radolan.nrow)) .* radolan.precision
        radolan.data = radolan.data .* 0.5 .- 32.5

    elseif radolan.producttype ∈ ["RY", "RW","RH","AY","RV","YW"]
        data_tmp::Vector{UInt16} = reinterpret(UInt16, dat)
        radolan.nodatamask = reshape(.!(iszero.(data_tmp .& 0x2000)), (radolan.ncol, radolan.nrow))
        radolan.secondary = reshape(.!(iszero.(data_tmp .& 0x1000)), (radolan.ncol, radolan.nrow))
        radolan.cluttermask = reshape(.!(iszero.(data_tmp .& 0x8000)), (radolan.ncol, radolan.nrow))
        data_tmp .&= 0xFFF
        radolan.data = reshape(data_tmp, (radolan.ncol, radolan.nrow)) .* radolan.precision
        
    end
    @inbounds radolan.data[radolan.nodatamask] .= missing_val
    if !ismissing(missing_val)
        radolan.data = convert(Array{Float64}, radolan.data)
    end
end

function read_radolan_composite(path, filename, missing_val=missing)
    raw = read(joinpath(path, filename))
    rad_obj, raw_without_header = read_header_to_dict(raw)
    read_array(rad_obj, raw_without_header, missing_val)
    return rad_obj
end


function read_radolan(path, filename; missing_val=NaN)
    read_radolan_composite(path, filename,missing_val).data
end

