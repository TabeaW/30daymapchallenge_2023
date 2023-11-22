using ArchGDAL
using DataFrames
function loadGermanBorder()
    layer=ArchGDAL.getlayer(ArchGDAL.read(joinpath(@__DIR__,"nuts5000_12-31.gk3.shape/nuts5000/NUTS5000_N1.shp")),0)
    bounds=DataFrame(layer);
    ArchGDAL.createcoordtrans(ArchGDAL.importURL("https://sg.geodatenzentrum.de/web_public/gdz/dokumentation/crs/gk3.prj"), ArchGDAL.importPROJ4("+proj=longlat +datum=WGS84 +no_defs +type=crs")) do transform
        for x in eachrow(bounds)
        ArchGDAL.transform!(x."", transform)
        end
    end
    return bounds
end

using GeoJSON
function plotWorldBorders(ax,kwargs...)
    worldCountries = GeoJSON.read(joinpath(@__DIR__,"ne_110m_admin_0_countries.geojson"))
    poly!(
    ax, worldCountries;
    color= :transparent,
    strokecolor = :black,
    strokewidth = 1,kwargs...)
end

latlon_string="+proj=longlat +datum=WGS84 +no_defs +type=crs"

radolan_string="+proj=stere +lat_0=90 +lat_ts=90 +lon_0=10 +k=0.93301270189 +x_0=0 +y_0=0 +a=6370040 +b=6370040  +no_defs"
radolan_rw_string="+proj=stere +lat_0=90 +lat_ts=90 +lon_0=10 +k=0.93301270189 +x_0=0 +y_0=0 +a=6370040 +b=6370040 +to_meter=1000 +no_defs"