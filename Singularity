Bootstrap: docker
From: continuumio/miniconda3@sha256:3a2017213a16daff5bc8dec8571354249c3370d6b0d64ac78e8538257ce42d4c

%post
    apt-get update
    apt-get install -y libxml2 libxml2-dev file coinor-cbc
    # install build tools, zlib, and git
    apt-get install -y build-essential zlib1g-dev git
    apt-get install -y default-jre
    mkdir -p /opt/conda/config
    export CONDARC=/opt/conda/config/.condarc
    touch /opt/conda/config/.condarc
    # Install Mamba
    conda install -c conda-forge mamba

    # Install Snakemake
    mamba install -c bioconda -c conda-forge snakemake coin-or-cbc=2.10.10

    # test if cbc is correctly installed, if not exit with error
    cbc ? || exit 1


    # configure strict channel priority
    conda config --set channel_priority strict
    conda init bash

    # Source the conda.sh script to ensure conda commands are available
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/profile
    echo "conda activate base" >>  /etc/profile

    # Create environments using snakemake
    . /opt/conda/etc/profile.d/conda.sh
    conda activate base
    # Verify that strict channel priority is set
    conda config --show | grep channel_priority

    cd /opt/pipeline
    # Pre-create EVERY conda environment referenced by ANY workflow so the
    # read-only runtime image never needs to build one. Snakefile_build_envs
    # references each envs/*.yaml directly, independent of any workflow DAG or
    # checkpoint (Snakefile_generic's discover_haplotypes checkpoint hides its
    # post-assembly rules from --conda-create-envs-only, so pointing it at that
    # workflow would silently miss those envs). Envs are created at the same
    # hashed paths the workflows resolve at runtime.
    snakemake -s Snakefile_build_envs --use-conda --conda-prefix /opt/conda/envs \
      --conda-create-envs-only --cores 4

    # Clean up
    mamba clean --all

    # make root accessible for everyone
    chmod -R 777 /root
    # remove all temp files

    # install hifiasm
    cd /opt
    git clone https://github.com/kavonrtep/hifiasm
    cd hifiasm && git checkout cf27b25413e01032283f5a030baac5dd3f1a8280 && make
    cp hifiasm /usr/local/bin

    # Stamp image labels from version.py so the release artefact carries
    # the pipeline source version it was built from.
    PIPELINE_VERSION=$(python3 /opt/pipeline/version.py)
    BUILD_DATE=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    mkdir -p /.singularity.d
    cat > /.singularity.d/labels.json <<EOF
{
    "Version": "${PIPELINE_VERSION}",
    "Build-Date": "${BUILD_DATE}"
}
EOF
    echo "Stamped container labels: Version=${PIPELINE_VERSION} Build-Date=${BUILD_DATE}"

%files
    envs /opt/pipeline/envs
    Snakefile /opt/pipeline/Snakefile
    Snakefile_generic /opt/pipeline/Snakefile_generic
    Snakefile_mapping /opt/pipeline/Snakefile_mapping
    Snakefile_build_envs /opt/pipeline/Snakefile_build_envs
    config_template.yaml /opt/pipeline/config.yaml
    config_template.yaml /opt/pipeline/config_template.yaml
    config_generic_template.yaml /opt/pipeline/config_generic_template.yaml
    config_mapping_template.yaml /opt/pipeline/config_mapping_template.yaml
    run_pipeline.py /opt/pipeline/run_pipeline.py
    run_mapping.py /opt/pipeline/run_mapping.py
    version.py /opt/pipeline/version.py
    scripts /opt/pipeline/scripts
    gepard-1.30 /opt/gepard-1.30



%environment
    export PATH=/opt/pipeline/scripts:/opt/conda/bin:$PATH
    export CONDA_ENVS_PATH=/opt/conda/envs
    export CONDA_PREFIX=/opt/conda
    export CONDARC=/opt/conda/config/.condarc
    export HOME=/root
    export GEPARD=/opt/gepard-1.30


%runscript
    # Navigate to the pipeline directory
    # set cache directory

    exec /opt/pipeline/run_pipeline.py "$@"
