
FROM ubuntu:18.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update 
RUN apt-get -y upgrade

RUN apt-get -y install build-essential git cmake wget apt-utils
RUN apt-get -y install gcc g++ gfortran 
RUN apt-get -y install libgomp1 

# Eigen, Blas, Lapack
RUN apt-get -y install libeigen3-dev libblas-dev liblapack-dev

# Boost
RUN apt-get -y install libboost-all-dev

# Qt
RUN apt-get -y install qtbase5-dev qt3d5-dev
