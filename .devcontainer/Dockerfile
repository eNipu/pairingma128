FROM ubuntu:20.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Add retry and timeout options to apt to make it more resilient
RUN echo 'Acquire::Retries "5";' > /etc/apt/apt.conf.d/80retries && \
    echo 'APT::Get::Assume-Yes "true";' > /etc/apt/apt.conf.d/90assumeyes && \
    echo 'Acquire::http::Timeout "180";' > /etc/apt/apt.conf.d/99timeout && \
    echo 'Acquire::https::Timeout "180";' > /etc/apt/apt.conf.d/99timeout

# Install required packages with increased resilience
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    gcc \
    g++ \
    gdb \
    git \
    libgmp-dev \
    make \
    python3 \
    valgrind \
    vim \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create a non-root user
RUN useradd -m developer \
    && echo "developer:developer" | chpasswd \
    && adduser developer sudo \
    && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# Set up the working directory
WORKDIR /workspace

# Switch to the non-root user
USER developer