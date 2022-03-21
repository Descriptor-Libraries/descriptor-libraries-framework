FROM postgres:12.2 AS builder

RUN apt update && apt upgrade -y

RUN apt install -y python3 \
                   python3-pip \
                   build-essential \
                   zip \
                   unzip \
                   cmake \
                   wget \
                   postgresql-client-12 \
                   postgresql-server-dev-12 \
                   postgresql-plpython3-12 \
                   libeigen3-dev \
                   libboost-dev \
                   libboost-iostreams-dev \
                   libboost-python-dev \
                   libboost-regex-dev \
                   libboost-serialization-dev \
                   libboost-system-dev \
                   libboost-thread-dev \
                   python3-numpy \
                   python3-dev

ARG RDKIT_VERSION=2020_03_2

RUN wget -q https://github.com/rdkit/rdkit/archive/Release_${RDKIT_VERSION}.tar.gz && \
    tar zxf Release_${RDKIT_VERSION}.tar.gz && \
    mv rdkit-Release_${RDKIT_VERSION} rdkit && \
    mkdir rdkit/build && cd rdkit/build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr \
             -DRDK_BUILD_PGSQL=on \
             -DRDK_BUILD_SWIG_JAVA_WRAPPER=OFF \
             -DRDK_BUILD_PYTHON_WRAPPERS=ON \
             -DRDK_SWIG_STATIC=OFF \
             -DPostgreSQL_ROOT=/usr \
             -DPostgreSQL_TYPE_INCLUDE_DIR=/usr/include/postgresql/12/server/ \
             -DRDK_BUILD_AVALON_SUPPORT=ON \
             -DRDK_BUILD_CAIRO_SUPPORT=OFF \
             -DRDK_BUILD_CPP_TESTS=OFF \
             -DRDK_BUILD_INCHI_SUPPORT=ON \
             -DRDK_BUILD_FREESASA_SUPPORT=ON \
             -DRDK_INSTALL_INTREE=OFF \
             -DRDK_INSTALL_STATIC_LIBS=OFF \
             -DPYTHON_EXECUTABLE=/usr/bin/python3 && \
    make -j `nproc` && \
    make install

FROM postgres:12.2 as postgres-rdkit

RUN apt update && apt upgrade -y

RUN apt install -y python3 \
                   python3-pip \
                   python3-numpy \
                   zip \
                   unzip \
                   postgresql-client-12 \
                   postgresql-plpython3-12 \
                   libboost-iostreams1.67.0 \
                   libboost-python1.67.0 \
                   libboost-regex1.67.0 \
                   libboost-serialization1.67.0 \
                   libboost-system1.67.0 \
                   libboost-thread1.67.0


# equivalent to pgsql_install.sh
COPY --from=builder /rdkit/build/Code/PgSQL/rdkit/rdkit--3.8.sql /usr/share/postgresql/12/extension
COPY --from=builder /rdkit/Code/PgSQL/rdkit/rdkit.control /usr/share/postgresql/12/extension
COPY --from=builder /rdkit/build/Code/PgSQL/rdkit/librdkit.so /usr/lib/postgresql/12/lib/rdkit.so

COPY --from=builder /usr/lib/libRDKit* /usr/lib/
COPY --from=builder /usr/lib/cmake/rdkit /usr/lib/cmake/rdkit
COPY --from=builder /usr/share/RDKit /usr/share/RDKit
COPY --from=builder /usr/include/rdkit /usr/include/rdkit
COPY --from=builder /usr/lib/python3/dist-packages/rdkit /usr/lib/python3/dist-packages/rdkit
