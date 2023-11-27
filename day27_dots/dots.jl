using Makie.StructArrays,GLMakie,CSV,DataFrames,GeoMakie
include("utils.jl")
function read_data(filename)
    stat=CSV.File(joinpath("day17_flow/data",filename),delim=";")|>DataFrame
    stat.AE_GB_POS_withmissings=replace(stat.AE_GB_POS,-999=>missing)
    stat.AE_GL_POS_withmissings=replace(stat.AE_GL_POS,-999=>missing)
    dropmissing!(stat)
    return stat
end
all_sonds=DataFrame()
for file in readdir("day17_flow/data")
    all_sonds=vcat(all_sonds,read_data(file))
end
points = StructArray{Point2f}((all_sonds.AE_GL_POS[1:end],all_sonds.AE_GB_POS[1:end]))
border=loadGermanBorder()
cmap = to_colormap(:BuPu)
cmap[1] = RGBAf(0, 0, 0, 1)
let
fig=Figure(backgroundcolor=:black)
ax=Axis(fig[1,1],backgroundcolor=:black,height=700,width=700)
plot = datashader!(ax, points;
            colormap=cmap,
            async=false)
p=poly!(ax,GeoMakie.to_multipoly.(GeoMakie.geo2basic.(border."")),color=:transparent,strokecolor=:white,strokewidth=0.5)

hidedecorations!(ax)
hidespines!(ax)
        
Label(fig[2,1],text="Data: © Deutscher Wetterdienst\nGeodata: © GeoBasis-DE / BKG (2023)",color=:white)
Makie.resize_to_layout!(fig)
save("dots_final.png",fig,px_per_unit=2)
fig
end