# Cropopticon

Cropopticon is a Julia package for interactively labelling agricultural field segments. It overlays segmentation results on satellite and aerial imagery, and lets you mark, merge, and export field boundaries with a point-and-click UI.

---

## Requirements

- [Julia 1.11 or later](https://julialang.org/downloads/)

---

## Installation

1. **Install Julia.**  
   Download and install the latest stable release from [julialang.org/downloads](https://julialang.org/downloads/).

2. **Clone or open the project.**  
   ```bash
   git clone <repo-url>
   cd Cropopticon
   ```

3. **Install dependencies.**  
   Start Julia with the project environment and instantiate all packages:
   ```bash
   julia --project=.
   ```
   Then in the Julia REPL:
   ```julia
   ] instantiate
   ```
   Press `Backspace` to exit package mode when done.

---

## Usage

Start a Julia REPL with the project environment:

```bash
julia --project=.
```

Import the package and run the labelling tool:

```julia
using Cropopticon

fig = run_labelling("config.yaml")
display(fig)
```

`run_labelling` reads all paths and settings from the YAML config file, loads the imagery and segments for the specified sample, and opens an interactive labelling window. Set `sample_id` to `random` to automatically pick an unlabelled sample from `segment_dir` (one that has no corresponding file yet in `dst_dir`).

---

## Configuration

All settings are stored in a YAML file. See [`config.yaml`](config.yaml) for a full example.

```yaml
figsize: 800
naip_dir: /path/to/naip
cdl_dir: /path/to/cdl
ae_dir: /path/to/alpha_earth
segment_dir: /path/to/sam_segments
dst_dir: /path/to/output
sample_id: california_2018_tile_005_patch_238
```

| Key | Type | Description |
|-----|------|-------------|
| `figsize` | Integer | Width and height (in pixels) of each image panel in the figure. |
| `naip_dir` | String | Directory containing NAIP GeoTIFF files. Each file must be named `<sample_id>.tif`. |
| `cdl_dir` | String | Directory containing USDA Cropland Data Layer (CDL) GeoTIFF files. Each file must be named `<sample_id>.tif`. |
| `ae_dir` | String | Directory containing AlphaEarth GeoTIFF files. Each file must be named `<sample_id>.tif`. |
| `segment_dir` | String | Directory containing SAM segmentation GeoTIFFs. Each file must be named `<sample_id>.tif`. |
| `dst_dir` | String | Directory where the output labelled raster will be written. The output file is named `<sample_id>.tif`. |
| `sample_id` | String | Identifier of the sample to label, or `random` to pick an unlabelled sample automatically. |

---

## Labelling UI

The window displays the NAIP image alongside the Cropland Data Layer (CDL) and AlphaEarth panels, with coloured polygon overlays for each segment drawn on top of the NAIP image.

### Mouse Controls

| Action | Effect |
|--------|--------|
| **Left click** on a segment | Toggle the segment as a field. Field segments are highlighted in red. |
| **Shift + Left click** on a segment | Toggle the segment for merging. Merge-selected segments are highlighted in blue. |
| **Ctrl + Left click** | Add a point to the line/polygon currently being drawn. |
| **Ctrl + Right click** | Close the current set of points into a polygon (used by **Fill** and **Erase**). |
| **Right click + drag** | Pan the image. |
| **Scroll wheel** | Zoom in and out. |

All image panels are linked — panning or zooming one panel moves all others simultaneously.

### Buttons

| Button | Description |
|--------|-------------|
| **Show/Hide** | Toggle the visibility of all polygon overlays. |
| **Merge** | Merge all blue-highlighted segments into a single region. The merged segments are morphologically closed to fill small gaps between them. |
| **Split** | Split segments along the line currently drawn with Ctrl+Left click. |
| **Fill** | Add the drawn polygon to the segment it overlaps most, and cut it out of any other overlapping segments. Requires a closed polygon (Ctrl+Right click). |
| **Erase** | Cut the drawn polygon out of any overlapping segments. Requires a closed polygon (Ctrl+Right click). |
| **Undo** | Revert the last merge, split, fill, or erase operation. |
| **Save** | Write the current field labels to `<dst_dir>/<sample_id>.tif` as a binary raster (1 = field, 0 = non-field). |

### Workflow

1. Scroll and pan to explore the image.
2. Left-click segments that correspond to agricultural fields to mark them red.
3. If two adjacent segments belong to the same field, Shift+click each one to highlight them blue, then press **Merge**.
4. Use Ctrl+Left click to draw a line or polygon boundary, then use **Split**, **Fill**, or **Erase** to adjust segment boundaries.
5. Use **Undo** to revert an edit if needed.
6. Press **Save** when labelling is complete. The output raster is written to `dst_dir`.