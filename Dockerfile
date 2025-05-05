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
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set up a non-root user for development
RUN useradd -m developer && \
    echo "developer:developer" | chpasswd && \
    adduser developer sudo

# Switch to the non-root user
USER developer
WORKDIR /home/developer

# Create directory for mounting the source code
RUN mkdir -p /home/developer/pairing

# Set default working directory to the source code mount
WORKDIR /home/developer/pairing

# Expose port for remote debugging
EXPOSE 2222

# Default command: start an interactive shell
CMD ["/bin/bash"]