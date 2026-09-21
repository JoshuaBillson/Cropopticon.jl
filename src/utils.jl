function linear_stretch(x::AbstractArray{T}, low, high) where {T <: Real}
    return clamp!((x .- low) ./ (high - low), T(0), T(1)) .|> ImageCore.N0f8
end

function raster_to_image(raster::AbstractArray{<:Real,3}, lb, ub)
    # Stretch the data into a (0, 1) range for visualization
    img_stretched = linear_stretch(raster, lb, ub)

    # Create an RGB image from the stretched data
    return ImageCore.colorview(ImageCore.RGB, permutedims(img_stretched, (3,2,1)))
end

function pca(X::AbstractMatrix{<:Real})
    # Randomly sample rows for PCA computation to speed up SVD on large datasets (e.g., 1 million samples)
    sample_size = min(size(X,1), 100_000_000)  # Limit to 100 million samples
    sample_indices = Random.randperm(size(X,1))[1:sample_size]
    sample = X[sample_indices, :]

    # Compute mean of the sample for centering
    μ = mean(sample; dims=1) |> vec

    # Economy SVD: X = U * Diagonal(S) * V'; columns of V are principal components
    V = LinearAlgebra.svd(sample; full=false).V     # (C, min(N,C))
    V = mapslices(x -> x * sign(x[1]), V, dims=1)  # Ensure consistent sign for principal components

    # Apply PCA transformation to the entire dataset
    transformed = (X .- reshape(μ, (1,:))) * V[:, 1:3]
    return transformed
end

function polygonize(mask::Rasters.AbstractRaster{<:Integer,2})
    # Get the coordinate ranges for the raster dimensions
    xdim, ydim = Rasters.dims(mask)
    xcoords = first(xdim):Rasters.span(xdim).step:last(xdim)
    ycoords = first(ydim):Rasters.span(ydim).step:last(ydim)

    # Get unique labels and polygonize each segment
    labels = Int[]
    polygons = GeoInterface.Polygon[]
    for label in unique_labels(mask.data)
        component_mask = ImageMorphology.label_components(mask.data .== label)
        for component_label in unique_labels(component_mask)
            component_polygon = GeometryOps.polygonize(xcoords, ycoords, component_mask .== component_label).geom
            for p in component_polygon
                push!(polygons, GeometryOps.simplify(p, tol=2))
                push!(labels, label)
            end
            #@info length(component_polygon)
            #_, i = findmax(GeometryOps.area.(component_polygon))
            #push!(polygons, GeoInterface.convert(GeometryBasics, GeometryOps.simplify(component_polygon[i], tol=2)))
            #push!(labels, label)
        end
    end

    return labels .=> polygons
end

function extent_to_polygon(dims)
    return GeometryOps.extent_to_polygon(Rasters.extent(dims))
end

function intersects(poly1::GeoInterface.Polygon, poly2::GeoInterface.Polygon)
    return GeometryOps.intersects(poly1, poly2)
end

function intersection(poly1::GeoInterface.Polygon, poly2::GeoInterface.Polygon)
    return GeometryOps.intersection(poly1, poly2, target=GeoInterface.PolygonTrait()) |> first
end

function union(poly1::GeoInterface.Polygon, poly2::GeoInterface.Polygon)
    return GeometryOps.union(poly1, poly2, target=GeoInterface.PolygonTrait()) |> first
end

function difference(poly1::GeoInterface.Polygon, poly2::GeoInterface.Polygon)
    return GeometryOps.difference(poly1, poly2, target=GeoInterface.PolygonTrait())
end

function simplify(poly::GeoInterface.Polygon, tol::Real=2)
    return GeometryOps.simplify(poly, tol=tol)
end

unique_labels(x) = filter(!iszero, unique(x)) |> sort!

function rasterize(polys::AbstractVector, dims, fill=nothing, missingval=0x00)
    fill = isnothing(fill) ? UInt8.(collect(eachindex(polys))) : fill
    return reduce([rasterize(p, dims, UInt8(f), missingval) for (f, p) in zip(fill, polys)]) do acc, x
        return ifelse.(x .== missingval, acc, x)
    end
end
function rasterize(poly, dims, fill=0x01, missingval=0x00)
    return Rasters.rasterize(first, poly, to=dims, fill=fill, missingval=missingval)
end

function remove_boundary_pixels(segments::AbstractMatrix{<:Integer}, buffer=1)
    dst = copy(segments)
    labels = filter(!iszero, unique(segments)) |> sort!
    for label1 in labels
        segment_mask = segments .== label1
        for label2 in neighboring_regions(segments, label1, distance=buffer)
            neighbor_mask = segments .== label2
            boundary_mask = ImageMorphology.dilate(neighbor_mask, r=buffer) .& segment_mask
            dst[boundary_mask] .= 0
        end
    end
    return dst
end

function neighboring_regions(mask, label::Integer; distance=1)
    label_mask = ImageMorphology.dilate(mask .== label, r=distance)
    neighboring_labels = unique(mask[label_mask])
    return setdiff(neighboring_labels, [0, label]) |> sort!
end

function plot_image!(ax, image, dims)
    xdims = dimrange(dims, Rasters.X)
    ydims = dimrange(dims, Rasters.Y)
    GLMakie.image!(ax, xdims, ydims, image')
end

dimrange(d) = first(d.val)..last(d.val)
dimrange(d, dim) = dimrange(Rasters.dims(d, dim))
