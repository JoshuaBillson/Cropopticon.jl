struct NAIP
    src::String
end

function plot!(ax, naip::NAIP)
    # Read the raster data and select the RGB bands (R, G, B)
    raster = Rasters.Raster(naip.src; lazy=false, raw=true)

    # Create an RGB image from the stretched data
    naip_img = raster_to_image(raster.data, 0x00, 0xff)

    # Plot the image
    plot_image!(ax, naip_img, Rasters.dims(raster))
end

struct CDL
    src::String
end

function plot!(ax, cdl::CDL)
    # Read the CDL raster data
    raster = Rasters.Raster(cdl.src, lazy=false, raw=true)
    raster = permutedims(raster, (2,1))
    legend = CSV.read("cdl_legend.csv", DataFramesMeta.DataFrame)
    labels = filter(!=(0), unique(raster.data)) |> sort!

    # Write colors based on cdl class
    dst = zeros(ImageCore.RGB{ImageCore.N0f8}, size(raster.data))
    for label in labels
        color_code = legend[legend.Code .== label, :Color][1]
        color = parse(Colorant, color_code)
        dst[raster.data .== label] .= color
    end

    # Plot the mask
    plot_image!(ax, dst, Rasters.dims(raster))
end

struct AlphaEarth
    src::String
end

function plot!(ax, alphaearth::AlphaEarth)
    # Load Raster
    raster = Rasters.Raster(alphaearth.src, lazy=false, raw=true)
    
    # Perform PCA on the raster data to reduce it to 3 bands for visualization
    W, H, C = size(raster.data)
    pca_result = @pipe reshape(raster.data, W*H, C) |> pca |> reshape(_, (W, H, 3))

    # Compute the lower and upper bounds for contrast stretching
    lb = mapslices(x -> quantile(vec(x), 0.02), pca_result, dims=(1,2))
    ub = mapslices(x -> quantile(vec(x), 0.98), pca_result, dims=(1,2))
    image = raster_to_image(pca_result, lb, ub)

    # Plot the image
    plot_image!(ax, image, Rasters.dims(raster))
end