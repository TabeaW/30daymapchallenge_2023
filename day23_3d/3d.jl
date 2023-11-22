
using POreClim
using HDF5
using CairoMakie,GeometryBasics,LinearAlgebra
path="30daymapchallenge/3d/data"
name="ras07-vol5minng01_sweeph5onem_dbzh_09-2023112021490200-tur-10832-hd5"

fig=Figure()
ax=Axis3(fig[1,1])
for name in readdir("30daymapchallenge/3d/data")
    @show name
    range_val=HDF5.h5readattr(joinpath(path,name),"dataset1/how/")["range"]/1000
    elevation=HDF5.h5readattr(joinpath(path,name),"dataset1/how/")["startelA"]
    lat,long,height=(48.585379,9.782675,768)
    values=read_hdf(path,name)
    phi=LinRange(0,2pi,360) #azimuth
    alpha=elevation
    beta=90 .-elevation
    alpha_rad=alpha.*pi./180
    beta_rad=beta.*pi./180
    a=0:0.25:range_val
    xe=[cos(phi_v)*sin(beta_v)*a_v for (phi_v,beta_v) in zip(phi,beta_rad), a_v in a]
    ye=[sin(phi_v)*sin(beta_v)*a_v for (phi_v,beta_v) in zip(phi,beta_rad), a_v in a]
    ze=[sin(alpha_v)*a_v for (phi_v,alpha_v) in zip(phi,alpha_rad), a_v in a]

    points=vec([Point3f(xv,yv,zv) for (xv,yv,zv) in zip(xe,ye,ze)])
    faces=decompose(QuadFace{GLIndex},Tesselation(Rect(0,0,1,1),size(ze)))
    normals=normalize.(points)
    gb_mesh=GeometryBasics.Mesh(meta(points;normals),faces)


    mesh!(ax,gb_mesh,color=RGBA(0.5,0.5,0.5,0.5))
end