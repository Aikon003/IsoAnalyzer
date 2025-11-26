FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive

# --- 1. System update, Add R 4.4 Repo, and Install System Packages ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    gpg-agent \
    wget \
    ca-certificates \
    apt-transport-https \
    && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/r-project.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | tee /etc/apt/sources.list.d/r-project.list \
    && apt-get update \
    && apt-get install -y \
        # --- R 4.4  ---
        r-base \
        r-base-dev \
        cutadapt \
        fastqc \
        trim-galore \
        xxd \
        gcc \
        g++ \
        gfortran \
        make \
        cmake \
        unzip \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        gdebi-core \
        gnupg \
        lsb-release \
        libudunits2-dev \
        libgdal-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libcairo2-dev \
        libgmp-dev \
        libmpfr-dev \
        libncurses-dev \
        cython3 \
        git \
        libblas-dev \
        liblapack-dev \
        libgsl-dev \
        parallel \
        python3.11 \
        python3.11-dev \
        python3.11-venv \
        python3-pip \
        python-is-python3 \
    # Configura python3.11 come default
    && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 \
    # Pulisci
    && rm -rf /var/lib/apt/lists/*

# --- 2. Python packages  ---
RUN pip install --no-cache-dir \
        numpy \
        pandas \
        scipy \
        cython \
        pysam \
        future \
        regex \
        matplotlib \
        umi_tools \
        argparse \
        gtftools \
        spladder

# --- 3. STAR aligner  ---
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz && \
    tar -xzf 2.7.11b.tar.gz && cd STAR-2.7.11b/source && make STAR && \
    cp STAR /usr/local/bin && cd / && rm -rf STAR-2.7.11b 2.7.11b.tar.gz

# --- 4. SAMtools  ---
RUN wget -q https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -xjf samtools-1.21.tar.bz2 && cd samtools-1.21 && \
    ./configure --prefix=/usr/local && make && make install && cd / && rm -rf samtools*

# --- 5. Subread (featureCounts)  ---
RUN wget "https://sourceforge.net/projects/subread/files/subread-2.1.1/subread-2.1.1-Linux-x86_64.tar.gz/download" -O subread.tar.gz && \
    tar -xzf subread.tar.gz && cd subread-2.1.1-Linux-x86_64 && \
    cp -r bin/* /usr/local/bin/ && cd / && rm -rf subread*

# --- 6. R packages  ---
RUN Rscript -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')"


RUN Rscript -e "BiocManager::install( \
    c('dplyr','readr','ggplot2','tidyr','stringr','readxl', \
      'IsoformSwitchAnalyzeR', 'org.Hs.eg.db', 'rtracklayer'), \
    ask=FALSE, update=FALSE, dependencies=TRUE, \
    lib = R.home('site-library') \
)"

# --- 7. UMI_parallel 
RUN git clone https://github.com/Milda85/UMI_parallel.git /opt/UMI_parallel \
    && chmod +x /opt/UMI_parallel/para_umi_dedup.sh \
    && chmod +x /opt/UMI_parallel/para_umi_extract.sh
ENV PATH="/opt/UMI_parallel/bin:${PATH}"

# --- 8. rMATS  ---
RUN pip install --upgrade "cython>=3.0.0"
RUN mkdir /rmats_build && cd /rmats_build && \
    git clone https://github.com/Xinglab/rmats-turbo.git && cd rmats-turbo && \
    echo '' > setup_environment.sh && ./build_rmats && \
    mkdir /rmats && cd /rmats && \
    cp /rmats_build/rmats-turbo/rmats.py ./ && \
    cp /rmats_build/rmats-turbo/cp_with_prefix.py ./ && \
    cp /rmats_build/rmats-turbo/*.so ./ && \
    mkdir rMATS_C && cp /rmats_build/rmats-turbo/rMATS_C/rMATSexe ./rMATS_C && \
    mkdir rMATS_P && cp /rmats_build/rmats-turbo/rMATS_P/*.py ./rMATS_P && \
    mkdir rMATS_R && cp /rmats_build/rmats-turbo/rMATS_R/*.R ./rMATS_R && \
    rm -rf /rmats_build
ENV PATH="/rmats:${PATH}"

# --- 9. RSEM  ---
RUN git clone https://github.com/deweylab/RSEM.git /opt/RSEM && cd /opt/RSEM && make && make install prefix=/usr/local && rm -rf /opt/RSEM

# --- 10. fastp  ---
RUN wget http://opengene.org/fastp/fastp -O /usr/local/bin/fastp && chmod a+x /usr/local/bin/fastp 

# --- BSgenome installation ---
RUN Rscript -e "BiocManager::install( \
    'BSgenome.Hsapiens.UCSC.hg38', \
    ask=FALSE, update=FALSE, dependencies=TRUE, \
    lib = R.home('site-library') \
)"

# --- 11. IsoformSwitchAnalyzeR dependencies ---
RUN git clone https://github.com/weishwu/isoformSwitchAnalyzeR.git /opt/isoformSwitchAnalyzeR \
    && chmod +x /opt/isoformSwitchAnalyzeR/split_AA_fasta.py \
    && chmod +x /opt/isoformSwitchAnalyzeR/modify_gtf.py
ENV PATH="/opt/isoformSwitchAnalyzeR:${PATH}"

# --- 11. Cleanup  ---
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*



