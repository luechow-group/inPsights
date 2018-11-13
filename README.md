# inPsights [![master status](https://git.rwth-aachen.de/luechow-group/inPsights/badges/master/pipeline.svg)](https://git.rwth-aachen.de/luechow-group/inPsights/commits/master)
inPsights is a toolset for wave function analysis and visualization written in C++. It contains and interfaces the QMC program Amolqc (written in Fortran) as an independent submodule.

## Purpose
It is intended to analyse maxima data generated by Amolqc by providing.

## Installation
### Installing required packages
To install inPsights, Git, CMake, the GNU Compiler Collection (gcc) and the Eigen3 library must be installed. This is most easily achieved by employing a packet manager.

#### MacOS
Make sure that the Xcode Command Line Tools are installed already with: 
```bash
xcode-select --install
```

To install the required packages on MacOS, the [homebrew package manager](https://brew.sh) can be used. 
It can be downloaded and installed from the command line as follows:
```bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```
To download the packages with homebrew` execute the following command in the terminal:
```bash
brew update
brew upgrade
brew install -y git cmake gcc lapack eigen boost qt
```

#### Ubuntu
To install the required packages on Ubuntu, the package manager `aptitude` can be used.
To download the required packages with aptitude` execute the following command in the terminal:
```bash
sudo apt-get -y install \
    build-essential git cmake \
    gcc g++ gfortran \
    libgomp1 libblas-dev liblapack-dev libeigen3-dev libboost-all-dev \ 
    qtbase5-dev qt3d5-dev 

```

#### Cloning inPsights
Make sure you are registered for the GitLab service of https://git.rwth-aachen.de and that you have an SSH-Key for your account on the local machine you are using.
Furthermore, make sure you have permission to all repositiories by asking the git-administrator of the luechow-group.

```bash
git clone git@git.rwth-aachen.de:luechow-group/inpsights.git
```
and check out the branch of interest e.g.`git checkout develop`. Next, execute
```bash
git submodule update --init --recursive
```
to initialize all the submodules.

The next time you checkout a branch e.g.
```bash
git checkout develop
```
the submodules do not need to be initialized again afterwards. Thus
```bash
git submodule update --recursive
```
should be sufficient.

### Setting environment variables

#### Amolqc
The path to `Amolqc` *must* be exported as an environment variable e.g.

```bash
export AMOLQC=/Users/michaelheuer/Projects/inPsights/src/AmolqcInterface/Amolqc
```
to allow for interfacing `Amolqc` in the `AmolqcInterface` module.

#### Compilers
Environment variables for the different compilers and the libraries must be specified e.g. 

for the `GNU Compiler Collection`
```bash
export CC=/usr/local/bin/gcc-8
export CXX=/usr/local/bin/g++-8
export FC=/usr/local/bin/gfortran-8
export CMAKE_PREFIX_PATH=/usr/local/Cellar/qt/5.11.2
```
or `Intel Parallel Studio XE`
```bash
export FC=/opt/intel/bin/ifort
export CC=/opt/intel/bin/icc
export CXX=/opt/intel/bin/icpc
export INTELROOT=/opt/intel/
export MKLROOT=/opt/intel/mkl
export CMAKE_PREFIX_PATH=/usr/local/Cellar/qt/5.11.2
```
Otherwise the default compilers are used.

## Running inPsights from the Command-Line

Create a build directory e.g.
```bash
mkdir cmake-build-release
cd cmake-build-release
```
and conduct an out-of-source-build:

```bash
cmake ..
make
```
To build only a specific target write
```bash
make TargetName
```
instead. 

