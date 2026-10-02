module Cropopticon

import ArchGDAL
import CSV
import DataFramesMeta
import GLMakie
import GeoInterface
import GeometryBasics
import GeometryOps
import ImageCore
import ImageMorphology
import Random
import Rasters
import StatsBase
import LibGEOS
import LinearAlgebra
import YAML

import GLMakie: (..)

using Pipe: @pipe
using Colors
using Statistics
using Dates
using Accessors: @set
using InvertedIndices: Not

const GLM = GLMakie

include("utils.jl")
include("data_sources.jl")

struct State{P,D}
    dims::D
    polygons::Vector{P}
    show_polygons::Bool
end

struct Polygon{G}
    geom::G
    isfield::Bool
    tomerge::Bool
end

intersects(poly1::Polygon, poly2::Polygon) = intersects(poly1.geom, poly2.geom)
intersects(poly1::Polygon, poly2::GeoInterface.Polygon) = intersects(poly1.geom, poly2)
intersects(poly1::GeoInterface.Polygon, poly2::Polygon) = intersects(poly1, poly2.geom)

function intersection(poly1::Polygon, poly2::GeoInterface.Polygon)
    return @set poly1.geom = intersection(poly1.geom, poly2)
end

function union(poly1::Polygon, poly2::GeoInterface.Polygon)
    return @set poly1.geom = union(poly1.geom, poly2)
end

function difference(poly1::Polygon, poly2::GeoInterface.Polygon)
    map(difference(poly1.geom, poly2)) do poly
        return @set poly1.geom = poly
    end
end

function simplify(poly::Polygon, tol::Real=2)
    return @set poly.geom = GeometryOps.simplify(poly.geom, tol=tol)
end

"""
Create a new State object from a raster of segments.

# Arguments
- `segments`: The input raster containing segment labels as integers.

# Returns
- `State`: A new State object containing the polygons, labels, and other metadata for the segments.
"""
function State(segments::Rasters.AbstractRaster{<:Integer,2})
    polygon_results = polygonize(segments)
    polygons = isempty(polygon_results) ? Polygon{Tuple{Int,Int}}[] : [Polygon(p, false, false) for p in last.(polygon_results)]
    return State(Rasters.dims(segments), polygons, true)
end

function polygonindex(state::State, point)
    for (i, polygon) in enumerate(state.polygons)
        if GeometryOps.contains(polygon.geom, point)
            return i
        end
    end
    return 0
end

getpolygon(state::State, segment_index) = state.polygons[segment_index].geom

isfield(state, segment_index) = state.polygons[segment_index].isfield

tomerge(state, segment_index) = state.polygons[segment_index].tomerge

function toggle_field(state::State, segment_index::Int)
    return @set state.polygons[segment_index].isfield = !isfield(state, segment_index)
end

function toggle_merge(state::State, segment_index::Int)
    return @set state.polygons[segment_index].tomerge = !tomerge(state, segment_index)
end

function merge_regions(state::State)
    polys_to_merge_indices = [i for i in eachindex(state.polygons) if tomerge(state, i)]

    if isempty(polys_to_merge_indices)
        return state
    end

    polys_to_merge = [getpolygon(state, i) for i in polys_to_merge_indices]

    rasterized_polys = rasterize(polys_to_merge, state.dims)
    merged_polys = ImageMorphology.closing(rasterized_polys .> 0x00, r=8) # Merge the regions and fill small holes/gaps between them

    # Update the polygons, labels, isfield, and tomerge arrays to reflect the merged regions
    polygon_results = polygonize(merged_polys)
    polygons = [Polygon(p, false, false) for p in last.(polygon_results)]

    return @set state.polygons = vcat(state.polygons[Not(polys_to_merge_indices)], polygons)
end

function split_regions(state, line)
    line_poly = GeometryOps.buffer(line, 1)
    return erase_regions(state, line_poly)
end

function fill_regions(state::State{P}, region) where P
    # Crop region to tile extent
    region = intersection(region, extent_to_polygon(state.dims))

    # Find the indices of polygons that intersect with the region to be filled
    overlapping_polygon_indices = findall(i -> intersects(getpolygon(state, i), region), eachindex(state.polygons))

    # If there are no overlapping polygons, simply append the region as a new polygon
    if isempty(overlapping_polygon_indices)
        return @set state.polygons = vcat(state.polygons, [Polygon(region, false, false)])
    end

    # Determine which overlapping polygon has the maximum overlap with the region to be filled
    updated_polygons = P[]
    max_overlap_index = argmax(i -> GeometryOps.area(intersection(getpolygon(state, i), region)), overlapping_polygon_indices)
    for (i, poly) in zip(overlapping_polygon_indices, state.polygons[overlapping_polygon_indices])
        if i == max_overlap_index
            push!(updated_polygons, union(poly, region))
        else
            for new_poly in difference(poly, region)
                push!(updated_polygons, new_poly)
            end
        end
    end

    # Combine the polygons to keep with the updated polygons
    return @set state.polygons = vcat(state.polygons[Not(overlapping_polygon_indices)], updated_polygons)
end

function erase_regions(state::State{P}, region) where P
    # Find the indices of polygons that intersect with the region to be erased
    overlapping_polygon_indices = findall(i -> intersects(getpolygon(state, i), region), eachindex(state.polygons))

    # Compute the difference between each overlapping polygon and the region to be erased
    updated_polygons = P[]
    for poly1 in state.polygons[overlapping_polygon_indices]
        for poly2 in difference(poly1, region)
            push!(updated_polygons, poly2)
        end
    end
    
    # Combine the polygons to keep with the updated polygons
    return @set state.polygons = vcat(state.polygons[Not(overlapping_polygon_indices)], updated_polygons)
end

function render_state!(ax, state::GLM.Observable{<:State})
    show_polygons = GLM.lift(s -> s.show_polygons, state)
    polygon_alphas = GLM.lift(s -> [(isfield(s, l) || tomerge(s, l)) ? 0.2 : 0.05 for l in eachindex(s.polygons)], state)
    polygon_colors = GLM.lift(s -> [tomerge(s, l) ? :blue : :red for l in eachindex(s.polygons)], state)
    strokes = GLM.lift((colors, show) -> [(color, show * 0.8) for color in colors], polygon_colors, show_polygons)
    colors = GLM.lift((colors, alphas, show) -> [(c, show * a) for (c, a) in zip(colors, alphas)], polygon_colors, polygon_alphas, show_polygons)
    polygons = GLM.lift(s -> [GeoInterface.convert(GeometryBasics, p.geom) for p in s.polygons], state)
    GLM.poly!(ax, polygons, color=colors, strokecolor=strokes, strokewidth=1.0)
end

"""
Pan the given axes by holding Space and moving the mouse (no click needed).

Replaces Makie's default right-drag pan. Because the axes are linked, panning
whichever axis is under the cursor moves all of them.
"""
function enable_space_pan!(fig, axes...)
    # Remove the default right-click-drag pan from every axis
    for ax in axes
        GLM.deregister_interaction!(ax, :dragpan)
    end

    active = Ref{Any}(nothing)      # axis being panned while Space is held
    last_pos = Ref((0.0, 0.0))      # last mouse position in window pixels

    # Latch onto the panel currently under the cursor
    function latch!()
        idx = findfirst(ax -> GLM.is_mouseinside(ax.scene), axes)
        active[] = isnothing(idx) ? nothing : axes[idx]
        last_pos[] = GLM.events(fig).mouseposition[]
    end

    GLM.on(GLM.events(fig).keyboardbutton) do event
        if event.key == GLM.Keyboard.space
            if event.action == GLM.Keyboard.press
                latch!()
            elseif event.action == GLM.Keyboard.release
                active[] = nothing
            end
        end
    end

    GLM.on(GLM.events(fig).mouseposition) do pos
        # Only pan while Space is held
        if !GLM.ispressed(fig, GLM.Keyboard.space)
            active[] = nothing
            return
        end

        # Space is down but the cursor wasn't over a panel yet
        if isnothing(active[])
            latch!()
            return
        end
        ax = active[]

        # Convert the pixel movement into data units for this axis
        dx, dy = pos .- last_pos[]
        last_pos[] = pos
        lims = ax.finallimits[]
        viewport = ax.scene.viewport[]
        xscale = lims.widths[1] / viewport.widths[1]
        yscale = lims.widths[2] / viewport.widths[2]

        # Content follows the cursor, so the view window moves the opposite way
        x0 = lims.origin[1] - dx * xscale
        y0 = lims.origin[2] - dy * yscale
        GLM.limits!(ax, x0, x0 + lims.widths[1], y0, y0 + lims.widths[2])
    end
end

function run_labelling(config::String)
    cfg = YAML.load_file(config)
    sample_id = cfg["sample_id"] == "random" ? random_unlabelled_sample(cfg["segment_dir"], cfg["dst_dir"]) : cfg["sample_id"]
    return run_labelling(
        cfg["naip_dir"], 
        cfg["cdl_dir"], 
        cfg["ae_dir"], 
        cfg["segment_dir"], 
        cfg["dst_dir"], 
        sample_id,
        cfg["figsize"]
    )
end

function run_labelling(naip_dir::String, cdl_dir::String, ae_dir::String, segment_dir::String, dst_dir::String, sample_id::String, figsize)
    # Create figure and layout
    @info "Labelling sample $sample_id"
    fig = GLM.Figure(size=(figsize * 1.40, figsize))

    # Generate mask segments
    segments = Rasters.Raster(joinpath(segment_dir, sample_id * ".tif"), lazy=false, raw=true)
    segments = Rasters.replace_missing(segments, 0)
    segments = Rasters.modify(x -> ImageMorphology.label_components(x), segments)

    # Polygonize segments
    history = GLM.Observable([State(segments)])
    state = GLM.lift(h -> h[end], history)

    # Plot the NAIP image
    naip_src = NAIP(joinpath(naip_dir, sample_id * ".tif"))
    ax1 = GLM.Axis(fig[1:2, 1], aspect=GLM.DataAspect(), yreversed=false)
    GLM.deregister_interaction!(ax1, :rectanglezoom)
    GLM.deregister_interaction!(ax1, :limitreset)
    GLM.hidedecorations!(ax1)
    plot!(ax1, naip_src)

    # Plot CDL mask
    ax2 = GLM.Axis(fig[2, 2], aspect=GLM.DataAspect(), yreversed=false)
    GLM.deregister_interaction!(ax2, :rectanglezoom)
    GLM.hidedecorations!(ax2)
    cdl_src = CDL(joinpath(cdl_dir, sample_id * ".tif"))
    plot!(ax2, cdl_src)

    # Plot AlphaEarth
    ax3 = GLM.Axis(fig[1, 2], aspect=GLM.DataAspect(), yreversed=false)
    ae_src = AlphaEarth(joinpath(ae_dir, sample_id * ".tif"))
    GLM.deregister_interaction!(ax3, :rectanglezoom)
    GLM.hidedecorations!(ax3)
    plot!(ax3, ae_src)

    # Plot the polygons on top of the NAIP image
    render_state!(ax1, state)

    # Interactivity
    points = GLM.Observable(Tuple{Int,Int}[])
    poly = GLM.Observable(Tuple{Float64,Float64}[])
    GLM.on(GLM.events(fig).mousebutton) do event

        # Toggle Segment Label on Left Mouse Button Press
        if event.button == GLM.Mouse.left && event.action == GLM.Mouse.press && GLM.is_mouseinside(ax1.scene) && !GLM.ispressed(fig, GLM.Keyboard.space)

            # Determine selected segment
            world_pos_x, world_pos_y = GLM.mouseposition(ax1.scene)
            current_segment = polygonindex(state[], (world_pos_x, world_pos_y))

            if GLM.ispressed(fig, GLM.Keyboard.left_shift) && current_segment != 0 # Merge selected segment with other segments
                state[] = toggle_merge(state[], current_segment)
            elseif GLM.ispressed(fig, GLM.Keyboard.left_control) # Add point
                points[] = push!(points[], round.(Int, (world_pos_x, world_pos_y)))
            elseif current_segment != 0 # Toggle selected segment
                state[] = toggle_field(state[], current_segment)
            end
        end

        if event.button == GLM.Mouse.right && event.action == GLM.Mouse.press && GLM.is_mouseinside(ax1.scene)
            if GLM.ispressed(fig, GLM.Keyboard.left_control) && length(points[]) >= 3 # Close polygon
                verts = [(Float64(p[1]), Float64(p[2])) for p in points[]]
                poly[] = vcat(verts, verts[1:1])
                points[] = Tuple{Int,Int}[]
            end
        end
    end

    GLM.lines!(ax1, points)
    GLM.poly!(ax1, poly, color=(:blue,0.0), strokecolor=(:blue,1.0), strokewidth=2.0)

    # Live preview, like Serval's "select raster cells by polygon": a rubber-band
    # line from the last vertex to the cursor, a dashed line back to the first
    # vertex, and a shaded polygon showing what would be covered if closed now.
    cursor = GLM.Observable(GeometryBasics.Point2f(0, 0))
    update_cursor!(_) = (cursor[] = GeometryBasics.Point2f(GLM.mouseposition(ax1.scene)))
    GLM.on(update_cursor!, GLM.events(fig).mouseposition)
    GLM.on(update_cursor!, ax1.finallimits)   # stay in sync while panning/zooming

    preview_edge = GLM.lift(points, cursor) do pts, c
        isempty(pts) ? GeometryBasics.Point2f[] : [GeometryBasics.Point2f(pts[end]...), c]
    end
    preview_closing = GLM.lift(points, cursor) do pts, c
        isempty(pts) ? GeometryBasics.Point2f[] : [c, GeometryBasics.Point2f(pts[1]...)]
    end
    preview_fill = GLM.lift(points, cursor) do pts, c
        if length(pts) < 2
            GeometryBasics.Point2f[]
        else
            verts = GeometryBasics.Point2f[GeometryBasics.Point2f(p...) for p in pts]
            push!(verts, c)
            verts
        end
    end

    # xautolimits/yautolimits=false so the preview can never change the view
    GLM.poly!(ax1, preview_fill, color=(:yellow, 0.25), xautolimits=false, yautolimits=false)
    GLM.lines!(ax1, preview_closing, color=(:yellow, 0.6), linewidth=1.5, linestyle=:dash, xautolimits=false, yautolimits=false)
    GLM.lines!(ax1, preview_edge, color=:yellow, linewidth=2, xautolimits=false, yautolimits=false)

    buttongrid = GLM.GridLayout(fig[3, 1], tellwidth=false)
    showbutton = GLM.Button(buttongrid[1, 1], label="Show/Hide")
    GLM.on(showbutton.clicks) do _
        current_state = state[]
        state[] = @set current_state.show_polygons = !current_state.show_polygons
    end

    mergebutton = GLM.Button(buttongrid[1, 2], label="Merge")
    GLM.on(mergebutton.clicks) do _
        newstate = merge_regions(state[])
        history[] = vcat(history[], [newstate])
    end

    splitbutton = GLM.Button(buttongrid[1, 3], label="Split")
    GLM.on(splitbutton.clicks) do _
        line = GeometryBasics.LineString([GeometryBasics.Point2f(p...) for p in points[]])
        newstate = split_regions(state[], line)
        history[] = vcat(history[], [newstate])
        points[] = Tuple{Int,Int}[]
    end

    fillbutton = GLM.Button(buttongrid[1, 4], label="Fill")
    GLM.on(fillbutton.clicks) do _
        if length(poly[]) > 2
            poly_to_fill = GeoInterface.Polygon([map(collect, poly[])])
            newstate = fill_regions(state[], poly_to_fill)
            history[] = vcat(history[], [newstate])
            poly[] = Tuple{Float64,Float64}[]
        end
    end

    erasebutton = GLM.Button(buttongrid[1, 5], label="Erase")
    GLM.on(erasebutton.clicks) do _
        if length(poly[]) > 2
            poly_to_erase = GeoInterface.Polygon([map(collect, poly[])])
            newstate = erase_regions(state[], poly_to_erase)
            history[] = vcat(history[], [newstate])
            poly[] = Tuple{Float64,Float64}[]
        end
    end

    undo = GLM.Button(buttongrid[1, 6], label="Undo")
    GLM.on(undo.clicks) do _
        if length(history[]) > 1
            history[] = history[][1:end-1]
        end
    end

    # Write the current field labels to <dst_dir>/<sample_id>.tif
    function save_labels!()
        current_state = state[]
        polygons = [getpolygon(current_state, i) for i in eachindex(current_state.polygons) if isfield(current_state, i)]
        if isempty(polygons)
            # Nothing marked as field: write an all-zero mask so the sample still counts as done
            mask = Rasters.Raster(zeros(UInt8, map(length, current_state.dims)), current_state.dims; missingval=0x00)
        else
            mask = rasterize(polygons, current_state.dims)
            mask = Rasters.modify(x -> remove_boundary_pixels(x, 2), mask)
        end
        Rasters.write(joinpath(dst_dir, "$sample_id.tif"), mask, force=true)
        @info "Saved $sample_id"
    end

    savebutton = GLM.Button(buttongrid[1, 7], label="Save")
    GLM.on(savebutton.clicks) do _
        save_labels!()
    end

    # Next: save this sample, then load the next unlabelled one in the same window
    loading = Ref(false)   # ignore extra clicks while the next sample is loading
    nextbutton = GLM.Button(buttongrid[1, 8], label="Next")
    GLM.on(nextbutton.clicks) do _
        loading[] && return
        loading[] = true
        @async try
            save_labels!()
            next_id = next_unlabelled_sample(segment_dir, dst_dir)
            if isnothing(next_id)
                @info "No unlabelled samples left"
                loading[] = false
                return
            end
            newfig = run_labelling(naip_dir, cdl_dir, ae_dir, segment_dir, dst_dir, next_id, figsize)
            screen = GLM.Makie.getscreen(fig.scene)
            isnothing(screen) ? display(newfig) : display(screen, newfig)
        catch err
            loading[] = false
            @error "Next failed; staying on this sample" exception=(err, catch_backtrace())
        end
    end

    axes = [ax1]
    if cdl_dir !== nothing
        push!(axes, ax2)
    end
    if ae_dir !== nothing
        push!(axes, ax3)
    end
    GLM.linkaxes!(axes...)
    enable_space_pan!(fig, axes...)

    GLM.colsize!(fig.layout, 1, GLM.Relative(2/3))

    return fig
end

function random_unlabelled_sample(mask_dir::String, dst_dir::String)
    tif_stems(dir) = [splitext(f)[1] for f in readdir(dir) if endswith(f, ".tif")]
    candidates = setdiff(tif_stems(mask_dir), tif_stems(dst_dir))
    @assert !isempty(candidates) "No unlabelled samples found"
    return rand(candidates)
end

# Like random_unlabelled_sample, but returns `nothing` instead of erroring when
# every sample has been labelled (used by the Next button).
function next_unlabelled_sample(mask_dir::String, dst_dir::String)
    tif_stems(dir) = [splitext(f)[1] for f in readdir(dir) if endswith(f, ".tif")]
    candidates = setdiff(tif_stems(mask_dir), tif_stems(dst_dir))
    return isempty(candidates) ? nothing : rand(candidates)
end

end # module Cropopticon