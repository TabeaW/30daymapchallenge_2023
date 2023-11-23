using HDF5
using GLMakie,GeometryBasics,LinearAlgebra,Colors,DataFrames,Rasters,ArchGDAL,Proj
cd("day23_3d")
path="data"
name="ras07-vol5minng01_sweeph5onem_dbzh_09-2023112021490200-tur-10832-hd5"
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
dgm=Raster("../dgm1000.utm32s.geotiff/dgm1000/dgm1000_utm32s.tif")
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
fig=Figure()
Label(fig[2,1:2],text="Radar scans from Türkheim",fontsize=30)
ax=Axis3(fig[1,1],limits=(nothing,nothing,nothing,nothing,0,10000),height=900,width=900,zlabel="Height in meter",xlabel="",ylabel="",viewmode=:fit)
ax2=Axis3(fig[1,2],limits=(nothing,nothing,nothing,nothing,0,10000),height=900,width=900,zlabel="Height in meter",xlabel="",ylabel="",viewmode=:fit)
surface!(ax,dgm,colormap=cgrad(:oleron,rev=false,alpha=0.7),colorrange=(-1000,1000),shading=true)
surface!(ax2,dgm,colormap=cgrad(:oleron,rev=false,alpha=0.7),colorrange=(-1000,1000),shading=true)
hidespines!(ax)
hidexdecorations!(ax)
hideydecorations!(ax)
hidespines!(ax2)
hidedecorations!(ax2)
ss=nothing
Label(fig[3,1],text="Data: © Deutscher Wetterdienst\nGeodata: © GeoBasis-DE / BKG (2023) ",justification=:left)
Makie.resize_to_layout!(fig)
for name in readdir("data")
    @show name
    range_val=HDF5.h5readattr(joinpath(path,name),"dataset1/how/")["range"]
    elevation=HDF5.h5readattr(joinpath(path,name),"dataset1/how/")["startelA"]
    lat,long,height=(HDF5.h5readattr(joinpath(path,name),"where")["lat"],HDF5.h5readattr(joinpath(path,name),"where")["lon"],HDF5.h5readattr(joinpath(path,name),"where")["height"])
    transform =Proj.Transformation("+proj=longlat +datum=WGS84 +no_defs +type=crs",ArchGDAL.toPROJ4(ArchGDAL.importURL("http://sg.geodatenzentrum.de/web_public/gdz/dokumentation/crs/utm32s.prj"))) 
    long,lat=    transform([long,lat])
    @show lat,long
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
    #points=vec([Point3f(xv,yv,zv) for (xv,yv,zv) in zip(xe,ye,ze)])
    #faces_mesh=decompose(QuadFace{GLIndex},Tesselation(Rect(0,0,1,1),size(ze)))
    #normals=normalize.(points)
    #gb_mesh=GeometryBasics.Mesh(meta(points;normals),faces_mesh)

    ss=surface!(ax,xe,ye,ze,color=vals,colormap=cgrad(:ice,rev=true,alpha=0.5),colorrange=(0,45),lowclip=RGBA(0.7,0.7,0.7,0.06),shading=true)
    surface!(ax2,xe,ye,ze,color=vals,colormap=cgrad(:ice,rev=true,alpha=0.5),colorrange=(0,45),lowclip=RGBA(0.7,0.7,0.7,0.06),shading=true)
    #mesh!(ax,gb_mesh,color=vals,colormap=:fire,colorrange=(-32,40))
    #break
end
Colorbar(fig[3,2],vertical=false,ss,label="Reflectivity in dBZ")
