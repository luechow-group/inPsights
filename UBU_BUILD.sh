#! /bin/sh

mkdir -p build_DEBUG
mkdir -p build_RELEASE
cd build_DEBUG
cmake -DCMAKE_BUILD_TYPE=Debug ..
make all
./src/tests/LocalSpinMultiplicity_RunAllTests
cd ../build_RELEASE
cmake -DCMAKE_BUILD_TYPE=Release ..
make all
cd ..
