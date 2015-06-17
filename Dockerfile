# Dockerfile for the CGAT Code Collection
# http://www.cgat.org/

# Let us use an Ubuntu base image
FROM ubuntu:12.04

# Contact person
MAINTAINER Sebastian Luna Valero, sebastian.luna.valero@gmail.com

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    fixincludes \
    unzip

# Install CGAT code
RUN wget --no-check-certificate https://raw.github.com/CGATOxford/cgat/master/install-CGAT-tools.sh && \
    mkdir /shared && \
    bash install-CGAT-tools.sh --cgat-devel --zip --location /shared

# Set environment variables
ENV PATH=/shared/conda-install/envs/cgat-devel/bin:$PATH

# Add an entry point to the cgat command
ENTRYPOINT ["/shared/conda-install/envs/cgat-devel/bin/cgat"]

# Create a shared folder between docker container and host
VOLUME ["/shared/data"]
