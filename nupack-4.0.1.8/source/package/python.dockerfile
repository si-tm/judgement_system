FROM archlinux:base-devel

###############################################################################

RUN pacman -Syy
#RUN pacman-key --refresh
RUN pacman -Syu --noconfirm elinks cmake git wget vim ninja patch clang \
    make openblas libffi curl unzip zip tar pkg-config gcc-fortran lapack

################################################################################

ENV HOME=/root
WORKDIR $HOME
RUN wget -nv https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN /bin/bash Miniconda3-latest-Linux-x86_64.sh -b
RUN /root/miniconda3/bin/conda update -n base -c defaults conda
RUN /root/miniconda3/bin/conda install numpy scipy pandas termcolor matplotlib pip jinja2 jupyterlab nbconvert
ENV PATH /root/miniconda3/bin:$PATH

################################################################################

WORKDIR /root
RUN git clone https://github.com/mfornace/vcpkg.git -b add-gecode-6.2.0
RUN ./vcpkg/bootstrap-vcpkg.sh
RUN ./vcpkg/vcpkg install armadillo tbb nlohmann-json gecode spdlog protobuf \
    boost-core boost-preprocessor boost-align boost-sort boost-variant yaml-cpp \
    boost-iterator boost-fusion boost-container boost-functional libsimdpp taskflow

################################################################################

COPY ./ $HOME/nupack/
RUN mkdir $HOME/nupack-build
WORKDIR $HOME/nupack-build

RUN cmake ../nupack -DLAPACK_LIBRARIES=/usr/lib64/liblapack.so -DNUPACK_PIC=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE=/root/vcpkg/scripts/buildsystems/vcpkg.cmake \
    -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF \
    -DCMAKE_CXX_COMPILER=clang++ \
    -DCMAKE_CXX_COMPILER_AR=/usr/bin/ar \
    -DCMAKE_CXX_COMPILER_RANLIB=/usr/bin/ranlib \
    -G Ninja
    # -DCMAKE_POSITION_INDEPENDENT_CODE=ON \

# The following could crash/segfault if VM doesn't have enough memory to compile.
# Either limit the CPU via "--build-arg CPU=1" or allow more memory via "--memory=16g".
ARG CPU=10
RUN echo "-- Building with $CPU cores" && ninja -v -j $CPU -k 50 nupack-python
RUN pip install -e .

################################################################################
