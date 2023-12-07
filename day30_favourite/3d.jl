using HDF5
using GLMakie,GeometryBasics,LinearAlgebra,Colors,DataFrames,Rasters,ArchGDAL,Proj
cd("day30_favourite")
path="data"

function read_hdf(path,name)
    f=HDF5.h5open(joinpath(path,name),"r")
    a=f["dataset1"]["data1"]["data"]
    undetect=HDF5.h5readattr(joinpath(path,name),"dataset1/data1/what")["undetect"]
    gain=HDF5.h5readattr(joinpath(path,name),"dataset1/data1/what")["gain"]
    offset=HDF5.h5readattr(joinpath(path,name),"dataset1/data1/what")["offset"]
    nodata=HDF5.h5readattr(joinpath(path,name),"dataset1/data1/what")["nodata"]
    data::Matrix{Float64}=read(a).*gain.+offset
    data[isequal.(data,(nodata*gain+offset))].=NaN
    data[isequal.(data,(undetect*gain+offset))].=-32.5
    HDF5.close(f)
    return data
end
dgm=Raster("../dgm200.utm32s.geotiff/dgm200/dgm200_utm32s.tif")
function loadGermanBorder_utm()
    layer=ArchGDAL.getlayer(ArchGDAL.read(joinpath("..","nuts5000_12-31.gk3.shape/nuts5000/NUTS5000_N1.shp")),0)
    bounds=DataFrame(layer);
    ArchGDAL.createcoordtrans(ArchGDAL.importURL("https://sg.geodatenzentrum.de/web_public/gdz/dokumentation/crs/gk3.prj"), ArchGDAL.importURL("http://sg.geodatenzentrum.de/web_public/gdz/dokumentation/crs/utm32s.prj")) do transform
        for x in eachrow(bounds)
        ArchGDAL.transform!(x."", transform)
        end
    end
    return bounds
end
borders=loadGermanBorder_utm()
transform =Proj.Transformation("+proj=longlat +datum=WGS84 +no_defs +type=crs",ArchGDAL.toPROJ4(ArchGDAL.importURL("http://sg.geodatenzentrum.de/web_public/gdz/dokumentation/crs/utm32s.prj"))) 

fig=Figure()
ax=Axis3(fig[1,1],limits=(nothing,nothing,nothing,nothing,0,6000),height=1200,width=1200,zlabel="",xlabel="",ylabel="",viewmode=:fit)
surface!(ax,dgm,colormap=cgrad(:oleron,rev=false,alpha=0.7),colorrange=(-1000,1000),shading=MultiLightShading)
hidespines!(ax)
hidedecorations!(ax)
ss=nothing
#Label(fig[2,1],text="Data: © Deutscher Wetterdienst\nGeodata: © GeoBasis-DE / BKG (2023) ",justification=:left)
for name in readdir("data")
 
    @show name
    if !startswith(name,"ras07-pcp")
        continue
    end
    range_val=HDF5.h5readattr(joinpath(path,name),"dataset1/how/")["range"]
    elevation=HDF5.h5readattr(joinpath(path,name),"dataset1/how/")["startelA"]
    lat,long,height=(HDF5.h5readattr(joinpath(path,name),"where")["lat"],HDF5.h5readattr(joinpath(path,name),"where")["lon"],HDF5.h5readattr(joinpath(path,name),"where")["height"])
    long,lat=    transform([long,lat])
    vals=read_hdf(path,name)
    phi=LinRange(0,2pi,360) #azimuth
    alpha=elevation
    beta=90 .-elevation
    alpha_rad=alpha.*pi./180
    beta_rad=beta.*pi./180
    a=250:250:range_val
    ye=[cos(phi_v)*sin(beta_v)*a_v for a_v in a,(phi_v,beta_v) in zip(phi,beta_rad)] .+lat
    xe=[sin(phi_v)*sin(beta_v)*a_v for a_v in a,(phi_v,beta_v) in zip(phi,beta_rad)] .+long
    ze=[sin(alpha_v)*a_v for a_v in a,(phi_v,alpha_v) in zip(phi,alpha_rad)].+height
    ss=surface!(ax,xe,ye,ze,color=vals,colormap=cgrad(:ice,rev=true,alpha=0.5),colorrange=(0,45),lowclip=RGBA(0.7,0.7,0.7,0.3),shading=NoShading)
end