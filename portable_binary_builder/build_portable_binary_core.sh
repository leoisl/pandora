#!/bin/bash

# Builds portable pandora
# Inspired by Páll Melsted blog (https://pmelsted.wordpress.com/2015/10/14/building-binaries-for-bioinformatics/),
# on how he and the other authors managed to make kallisto (Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter,
# Near-optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, 525–527 (2016), doi:10.1038/nbt.3519)
# portable in different linux distributions.

set -e

# Activate Holy Build Box environment.
source /hbb_exe/activate

set -eux

# install packages needed to compile
yum install wget git -y

# install boost
wget https://sourceforge.net/projects/boost/files/boost/1.62.0/boost_1_62_0.tar.gz -O - | tar xzf -
cd boost_1_62_0
./bootstrap.sh --prefix=/usr/ --with-libraries=system,filesystem,iostreams,log,thread,date_time
./b2 install
cd ..

# compile pandora
cd io
mkdir build_portable_executable
cd build_portable_executable
cmake -DCMAKE_BUILD_TYPE=RELEASE_WITH_ASSERTS ..
make VERBOSE=1 -j 8
ctest -VV

# verify if the binary is portable
set +e
libcheck pandora #TODO: libomp is still a shared library dependency
set -e

# print ldd output for us to check the dependencies also
ldd pandora

# copy binary to host filesystem
cp pandora /io/pandora-linux-precompiled
