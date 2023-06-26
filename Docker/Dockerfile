# Based on Dockerfile from Dfam TETools https://github.com/Dfam-consortium/TETools, with added database configuration, conda configuration, and Earl Grey installation
# and configuration

FROM debian:latest AS builder

RUN apt-get -y update && apt-get -y install \
    curl gcc g++ make zlib1g-dev libgomp1 \
    perl \
    python3-h5py \
    libfile-which-perl \
    libtext-soundex-perl \
    libjson-perl liburi-perl libwww-perl

# Install dependencies and some basic utilities
RUN apt-get -y update \
    && apt-get -y install \
        aptitude \
        libgomp1 \
        perl \
        python3-h5py \
        libfile-which-perl \
        libtext-soundex-perl \
        libjson-perl liburi-perl libwww-perl \
        libdevel-size-perl \
    && aptitude install -y ~pstandard ~prequired \
        curl wget \
        vim nano \
        procps strace \
        libpam-systemd-

COPY src/* /opt/src/
WORKDIR /opt/src

# Extract RMBlast
RUN cd /opt \
    && mkdir rmblast \
    && tar --strip-components=1 -x -f src/rmblast-2.13.0+-x64-linux.tar.gz -C rmblast \
    && rm src/rmblast-2.13.0+-x64-linux.tar.gz

# Compile HMMER
RUN tar -x -f hmmer-3.3.2.tar.gz \
    && cd hmmer-3.3.2 \
    && ./configure --prefix=/opt/hmmer && make && make install \
    && make clean

# Compile TRF
RUN tar -x -f trf-4.09.1.tar.gz \
    && cd TRF-4.09.1 \
    && mkdir build && cd build \
    && ../configure && make && cp ./src/trf /opt/trf \
    && cd .. && rm -r build

# Compile RepeatScout
RUN tar -x -f RepeatScout-1.0.6.tar.gz \
    && cd RepeatScout-1.0.6 \
    && sed -i 's#^INSTDIR =.*#INSTDIR = /opt/RepeatScout#' Makefile \
    && make && make install

# Compile and configure RECON
RUN tar -x -f RECON-1.08.tar.gz \
    && mv RECON-1.08 ../RECON \
    && cd ../RECON \
    && make -C src && make -C src install \
    && cp 00README bin/ \
    && sed -i 's#^\$path =.*#$path = "/opt/RECON/bin";#' scripts/recon.pl

# Compile cd-hit
RUN tar -x -f cd-hit-v4.8.1-2019-0228.tar.gz \
    && cd cd-hit-v4.8.1-2019-0228 \
    && make && mkdir /opt/cd-hit && PREFIX=/opt/cd-hit make install

# Compile genometools (for ltrharvest)
RUN tar -x -f gt-1.6.2.tar.gz \
    && cd genometools-1.6.2 \
    && make -i -j4 cairo=no && make -i cairo=no prefix=/opt/genometools install \
    && make cleanup

# Configure LTR_retriever
RUN cd /opt \
    && tar -x -f src/LTR_retriever-2.9.0.tar.gz \
    && mv LTR_retriever-2.9.0 LTR_retriever \
    && cd LTR_retriever \
    && sed -i \
        -e 's#BLAST+=#BLAST+=/opt/rmblast/bin#' \
        -e 's#RepeatMasker=#RepeatMasker=/opt/RepeatMasker#' \
        -e 's#HMMER=#HMMER=/opt/hmmer/bin#' \
        -e 's#CDHIT=#CDHIT=/opt/cd-hit#' \
        paths

# Compile MAFFT
RUN tar -x -f mafft-7.471-without-extensions-src.tgz \
    && cd mafft-7.471-without-extensions/core \
    && sed -i 's#^PREFIX =.*#PREFIX = /opt/mafft#' Makefile \
    && make clean && make && make install \
    && make clean

# Compile NINJA
RUN cd /opt \
    && mkdir NINJA \
    && tar --strip-components=1 -x -f src/NINJA-cluster.tar.gz -C NINJA \
    && cd NINJA/NINJA \
    && make clean && make all

# Move UCSC tools
RUN mkdir /opt/ucsc_tools \
    && mv faToTwoBit twoBitInfo twoBitToFa  /opt/ucsc_tools \
    && chmod +x /opt/ucsc_tools/*
#COPY LICENSE.ucsc /opt/ucsc_tools/LICENSE

# Compile and configure coseg
RUN cd /opt \
    && tar -x -f src/coseg-0.2.2.tar.gz \
    && cd coseg \
    && sed -i 's@#!.*perl@#!/usr/bin/perl@' preprocessAlignments.pl runcoseg.pl refineConsSeqs.pl \
    && sed -i 's#use lib "/usr/local/RepeatMasker";#use lib "/opt/RepeatMasker";#' preprocessAlignments.pl \
    && make -i

# Configure RepeatMasker
RUN cd /opt \
    && tar -x -f src/RepeatMasker-4.1.4.tar.gz \
    && chmod a+w RepeatMasker/Libraries \
    && cd RepeatMasker/Libraries \
    && cp /opt/src/Dfam.h5.gz . \
    && gunzip -f Dfam.h5.gz \
    && cd /opt/RepeatMasker \
    #&& cp /opt/src/RepBaseRepeatMaskerEdition-20181026.tar.gz . \
    #&& tar -zxf RepBaseRepeatMaskerEdition-20181026.tar.gz \
    #&& wget https://github.com/rmhubley/RepeatMasker/blob/development/Libraries/RMRBMeta.embl \
    #&& cd /opt/RepeatMasker \
    && mv famdb.py famdb.py.bak \
    && cp /opt/src/famdb.py . \
    && chmod 755 famdb.py \ 
    && perl configure \
        -hmmer_dir=/opt/hmmer/bin \
        -rmblast_dir=/opt/rmblast/bin \
        -libdir=/opt/RepeatMasker/Libraries \
        -trf_prgm=/opt/trf \
        -default_search_engine=rmblast \
    && echo "RepeatMasker configured..." && cd .. && rm src/RepeatMasker-4.1.4.tar.gz

# Configure RepeatModeler
RUN cd /opt \
    && tar -x -f src/RepeatModeler-2.0.4.tar.gz \
    && mv RepeatModeler-2.0.4 RepeatModeler \
    && cd RepeatModeler \
    && perl configure \
         -cdhit_dir=/opt/cd-hit -genometools_dir=/opt/genometools/bin \
         -ltr_retriever_dir=/opt/LTR_retriever -mafft_dir=/opt/mafft/bin \
         -ninja_dir=/opt/NINJA/NINJA -recon_dir=/opt/RECON/bin \
         -repeatmasker_dir=/opt/RepeatMasker \
         -rmblast_dir=/opt/rmblast/bin -rscout_dir=/opt/RepeatScout \
         -trf_dir=/opt/ \
         -ucsctools_dir=/opt/ucsc_tools

FROM debian:latest

# Install dependencies and some basic utilities
RUN apt update && apt-get -y install git
RUN apt-get -y update \
    && apt-get -y install \
        aptitude \
        libgomp1 \
        perl \
        python3-h5py \
        libfile-which-perl \
        libtext-soundex-perl \
        libjson-perl liburi-perl libwww-perl \
        libdevel-size-perl \
    && aptitude install -y ~pstandard ~prequired \
        curl wget \
        vim nano \
        procps strace \
        libpam-systemd-
RUN cd /opt/ \
        && curl https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh --output anaconda.sh \
        && bash ./anaconda.sh -b -p /anaconda3 \
        && eval "$(/anaconda3/bin/conda shell.bash  hook)"


ENV PATH=$PATH:/opt/RepeatMasker:/opt/RepeatMasker/util:/opt/RepeatModeler:/opt/RepeatModeler/util:/opt/coseg:/opt/ucsc_tools:/opt:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/cd-hit/

COPY --from=builder /opt /opt
RUN echo "PS1='(earlgrey \$(pwd))\\\$ '" >> /etc/bash.bashrc
ENV LANG=C
ENV PYTHONIOENCODING=utf8

RUN apt install build-essential -y --no-install-recommends
RUN apt install libcurl4-openssl-dev -y && apt-get install libbz2-dev liblzma-dev libtiff-dev libfreetype6-dev pkgconf libfontconfig1-dev libharfbuzz-dev libfribidi-dev bc -y

ARG USER_ID
ARG GROUP_ID

RUN addgroup --gid $GROUP_ID user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user
RUN chmod a+rwx /opt/*
USER user

RUN cd ~ \
        && git clone https://github.com/TobyBaril/EarlGrey \
        && cd EarlGrey \
        && chmod +x ./configureForDocker \
        && eval "$(/anaconda3/bin/conda shell.bash  hook)" \
        && ./configureForDocker 

RUN cd ~ \
        && git clone https://github.com/ridgelab/SA-SSR \
        && cd SA-SSR/ \
        && sed -i "s|PREFIX=/usr/local/bin|PREFIX=/home/user/SA-SSR/|g" Makefile && make && make install  

ENV PATH=$PATH:/home/user/EarlGrey/:/home/user/SA-SSR/bin/

#RUN echo "#!/bin/bash" > /work/start.sh \
#        && echo 'eval "$(/anaconda3/bin/conda shell.bash  hook)"' >> /work/start.sh \
#        && echo "conda activate earlGrey" >> /work/start.sh \
#        && chmod +x /work/start.sh
