# User Manual for inPsights
After a successful build, different modules of inPsights are available as executables in $BUILD_DIR/bin.
The essential executables are
* `ProcessMaxima`
* `inPsights`.

The input for these executables has to be specified in YAML-format.

## `ProcessMaxima`  Executable
The executable `ProcessMaxima` processes binary files created with `Amolqc` in `rawdata` mode.
Exemplary `.bin` files can be found in the `Amolqc/examples/MAX/RAW` folder for acetone and ethane.

Additionally, an input file in YAML-format is required, where different options can be specified. 
For unspecified options, suitable defaults will be used automatically resulting in notifications in the console output.

Example: `input.yml`
```yaml
MaximaProcessing:
  binaryFileBasename: 'raw'
  samplesToAnalyze: 0         # 0 = all
  minimalClusterWeight: 0.0
  deleteCoreElectrons: false
Clustering:
  PreClusterer:
    radius: 0.01 # [a0]
  DensityBasedClusterer:
    radius: 0.2  # [a0]
VoxelCubeGeneration:
  generateVoxelCubesQ: true
VoxelCubeOverlapCalculation:
  calculateOverlapQ: false
```

The `.bin` files must be located in the same folder as the `input.yml` file.
A successful run of `ProcessMaxima input.yml` produces an `input-out.yml` file containing the processing results in YAML-format.

### Structure of the Input
Four top-level nodes must be present in the YAML input.
* MaximaProcessing (general options)
* Clustering (specification of clustering steps)
* VoxelCubeGeneration (calculation of SEDs)
* VoxelCubeOverlapCalculation (calculation of SED overlaps)


#### MaximaProcessing
General settings have to be given under the top-level YAML node `MaximaProcessing`:
* `binaryFileBasename` (`string`): basename of the `.bin` files (without file-extension and counters) created by an Amolqc run.  E.g. for `raw`, all `raw-00.bin`, `raw-01.bin` [...] files are considered.
* `samplesToAnalyze` (`unsigned integer`): Samples/maxima to analyze. A value of `0` means all samples/maxima from the subsequent `.bin` files are processed.
* `minimalClusterWeight` (`positive float`): minimal weight of clusters to be printed in the output file .
* `deleteCoreElectrons` (`bool`):  If `true`, all core electrons are deleted and thus not considered in the clustering or any statistic.

#### Clustering
The clustering process is specified under the top level YAML node `Clustering` and consists of a list of clusterers.
Subsequently, all specified clusterers are applied in their specified order and can be executed multiple times.

Example:
```yaml
Clustering:
  PreClusterer:      # 1. spherical pre-clustering with the PreClusterer and a small radius of 0.01 a0
    radius: 0.01  # [a0]
  DensityBasedClusterer:  # 2. density-based clustering with the DensityBasedClusterer and a radius of 0.2 a0
    radius: 0.2  # [a0]
```

Each clusterer starts with the clusters of the previous one clustering step and produces another hierachary level.
In the beginning, each maximum constitues an individual cluster. 
Note that only the last two clustering hierarchy levels are stored in the `-out.yml`.

##### PreClusterer
Greedy spherical clusterer employing a spin-agnostic best-match distance metric.
* `radius` (`positive float`,`[a0]`): Radius in which similar maxima (irrespective of spin) are clustered together.
* `valueIncrement`  (`positive float`,`[a0]`): Function value increment used in the greedy clusterer. This value is determined automatically from the standard error of the `-ln(|Ψ|^2)` if not specified.

##### SphericalClusterer
Spherical clusterer employing a spin-agnostic best-match distance metric.
* `radius` (`positive float`,`[a0]`): Radius in which similar maxima (irrespective of spin) are clustered together.
* `local` (`bool`): true unlocks the `ParticleSelection` Options in which the subset of considered electrons can be specified.


##### DensityBasedClusterer
Density-based clusterer employing a spin-agnostic best-match distance metric.
* `radius` (`positive float`,`[a0]`): radius in which similar, density connected maxima (irrespective of spin) are clustered together
* `local` (`bool`): true unlocks the `ParticleSelection` options in which the subset of considered electrons can be specified.


##### Local Clustering with the `ParticleSelection` Option
`ParticleSelection` is A sub-node that can be added to a `SphericalClusterer` or `DensityBasedClusterer` node.
The following options can be specified:
* `maximalCount` (`unsigned int`): Maximal number of electrons that are compared for clustering (the subset).
* `maximalDistance` (`positive float`,`[a0]`): Maximal distance of electrons from the reference positions to be included in the subset for comparison.
* `distanceMode` (`string`): Calculation of the distance for deciding which electrons to compare.
    * `minimum` The distance of an electron is calculated as the minimum of distances to the reference positions.
    * `average` The distance of an electron is calculated as the average of all distances to the reference positions.
* `invertSelection` (`bool`): If `true`, inverts the selection of electrons.
* `valenceOnly` (`bool`): If `true`, ignores non-valence electrons.
* `sortRemainder` (`bool`): If `true`, best-match sorts electrons that are not within the subset.
* `positions`: A list of positions, each of which can be given in one of the following ways:
    * `atAtom`: Give the index of an atom to add its position to the list.
    * `atCoordinates`: Give the cartesian coordinates as a three-membered list of floats.
    * `inBetween`: Give a list of positions (as above) to get the position in between all the given positions. 

Example: Density-based clustering of the four valence electrons closest to the point in between the 4. and 5. atom and sort the remaining electrons.
```yaml
Clustering:
  PreClusterer:
    radius: 0.01
  DensityBasedClusterer:
    radius: 0.2
    local: true
    sortRemainder: true
    ParticleSelection:
      maximalCount: 4
      distanceMode: minimum
      positions: [
        inBetween: [atAtom: 3, atAtom: 4]
      ]
      valenceOnly: true
```

#### VoxelCubeGeneration:
`VoxelCubeGeneration` specifies settings for the calculation of SEDs.
* `generateVoxelCubesQ` (`bool`): If `true`, SED voxel cubes are calculated.
* `centerCubesAtElectronsQ` (`bool`): If `true`, voxel cubes are centered at the electron positions of the maximum. If `false`, the cubes are centered at the origin. 
* `length` (`positive float`,`[a0]`):  side length of the considered volume being split into voxel cubes
* `dimension` (`unsigned int`): `dimension`^3 is the total number voxel cubes 
* `smoothingQ` (`bool`): If `true`, smoothing is enabled and each cube count is averaged with its next neighbors specified by `smoothingNeighbors`.
* ``smoothingNeighbors`` (`unsigned int`): number of next neighbors considered for the smoothing.

#### VoxelCubeOverlapCalculation:
`VoxelCubeOverlapCalculation` specifies settings for the calculation of SED overlaps. 
Since all cube origins have to be in the same center, a new grid has to be calculated compared to `VoxelCubeGeneration`. 
Possible options are:
* `calculateOverlapQ` (`bool`): If `true`, SED overlaps are calculated.
* `length` (`positive float`,`[a0]`):  side length of the considered volume being split into voxel cubes
* `dimension` (`unsigned int`): `dimension`^3 is the total number voxel cubes 

## `inPsights` Executable
The `inPsights` executable is the `Qt5` based GUI to visualize the processed data and requires an `-out.yml` file as an input.
Execution of `inPsights input-out.yml` from the command-line starts the GUI in a separate window.

