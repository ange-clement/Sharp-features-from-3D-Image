# Sharp-features-from-3D-Image

Source code for the paper "Sharp Features Graph Extraction from Segmented 3D Images for MeshGeneration" (Not yet published).

Allows to create sharp feature lines from 3D Segmented images that can be used to generate a mesh.

This repos is under construction.

## Building the executable

Using cmake:

```shell
git clone https://github.com/ange-clement/Sharp-features-from-3D-Image.git
cd Sharp-features-from-3D-Image
mkdir build ; cd build
cmake ..
```

Minimum system requirements: C++11 enabled compiler, [cmake](https://cmake.org), [cgal](https://cgal.org), [dgtal](https://www.dgtal.org/).

## Examples

By default, resulting figures are in `<folder_containing_executable>/<input_file_name>/<measure>/`

---
Extract sharp features from the `cut_cube_sphere.inr` image with Ambrosio-Tortorelli using the default parameters.

```shell
$: extract_sharp_features "input/cut_cube_sphere.inr"
```

---
Extract sharp features from the `cut_cube_sphere.inr` image with four measures (Ambrosio-Tortorelli, Curvature, Voronoi Covariance Measure, and the feedback measure) using the default parameters.

```shell
$: extract_sharp_features -a on -c on -v on -f on "input/cut_cube_sphere.inr"
```

---
Extract sharp features from the `anchor(30).off` mesh, by generating an 512x512x512 image, with Ambrosio-Tortorelli using the default parameters, and measure the line accuracy (extracted with an angle threshold of 30°) and the mesh accuracy.

```shell
$: extract_sharp_features "input/anchor(30).off" --extract_line_angle 30 --mesh_image_resolution 512 --measure_mesh_accuracy true --measure_line_accuracy true
```

Extract sharp features from the `anchor(30).off` mesh, by generating an 256x256x256 image with added Kanungo noise (`k=0.2`), with the Curvature measure only and the selection threshold at 30%, and measure the line accuracy.

```shell
$: extract_sharp_features "input/anchor(30).off" -a off -c on --selection_threshold 0.3 --extract_line_angle 30 --mesh_image_resolution 256 --noise_amount .3 --measure_line_accuracy true
```
