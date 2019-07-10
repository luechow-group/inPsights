# User manual for inPsights
* Different modules of inPsights are available as executables in $BUILD_DIR/bin.
* For all of these executables, the input has to be given in the YAML format.
## MaximaProcessing
* Handles the clustering (sorting) of maxima.
* General Settings have to be given under the top level YAML node `MaximaProcessing`:
    * `binaryFileBasename`: basename (without extension and counters) of the .bin files created by an Amolqc run.
    * `samplesToAnalyze`: number of samples to Analyze from the given binary files. `0` means all samples.
    * `minimalClusterWeight`: Minimal weight of clusters to be printed in the output file.
### Cluster
* `Clustering:` A top level node for the MaximaProcessing input. The subsequent called clusterers are applied in the order they are given in the input.
#### ReferencePositionsClusterer
* Clusters maxima only comparing a subset of electrons which are closest to given reference positions.
* `ReferencePositionsClusterer`: A sub-node of `Clustering` with the following settings:
    * `similarityRadius`: Clustering threshold referring to the single particle maximum euclidian distance.
    * `maximalCount`: Maximal number of electrons that are compared for clustering (the subset).
    * `maximalDistance`: Maximal distance of electrons from the reference positions to be included in the subset for comparison.
    * `distanceMode`: Calculation of the distance for deciding which electrons to compare.
        * `minimum` The distance of an electron is calculated as the minimum of distances to the reference positions.
        * `average` The distance of an electron is calculated as the average of all distances to the reference positions.
    * `invertSelection`: `true` or `false`, inverts the selection of electrons.
    * `valenceOnly`: `true` or `false`. If `true`, ignores non-valence electrons.
    * `positions`: A list of positions, each of which can be given in one of the following ways:
        * `atAtom`: Give the index of an atom to add its position to the list.
        * `atCoordinates`: Give the cartesian coordinates as a three-membered list of floats.
        * `inBetween`: Give a list of positions (as above) to get the position in between all the given positions. 
