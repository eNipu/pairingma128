FROM ubuntu:20.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install essential build tools and GMP
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    gcc \
    g++ \
    make \
    libgmp-dev \
    git \
    openssh-client \
    python3 \
    valgrind \
    gdb \
    lcov \
    gcovr \
    clang \
    clang-tools \
    clang-format \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set up a non-root user for development
RUN useradd -m developer && \
    echo "developer:developer" | chpasswd && \
    adduser developer sudo && \
    echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# Switch to the non-root user
USER developer
WORKDIR /home/developer

# Set up environment
ENV PATH="/usr/local/bin:${PATH}"

# Set default working directory
WORKDIR /workspace