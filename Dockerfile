# Image for INTERVENE Synthetic data pipeline
FROM debian:bullseye-slim

RUN set -eux; \
	apt-get update; \
	apt-get install -y --no-install-recommends \
		ca-certificates \
# ERROR: no download agent available; install curl, wget, or fetch
		curl \
	; \
	rm -rf /var/lib/apt/lists/*


# ********************* Build basic julia image ******************************************
ENV JULIA_PATH /usr/local/julia
ENV PATH $JULIA_PATH/bin:$PATH

# https://julialang.org/juliareleases.asc
# Julia (Binary signing key) <buildbot@julialang.org>
ENV JULIA_GPG 3673DF529D9049477F76B37566E3C7DC03D6E495

# https://julialang.org/downloads/
ENV JULIA_VERSION 1.6.4

RUN set -eux; \
	\
	savedAptMark="$(apt-mark showmanual)"; \
	if ! command -v gpg > /dev/null; then \
		apt-get update; \
		apt-get install -y --no-install-recommends \
			gnupg \
			dirmngr \
		; \
		rm -rf /var/lib/apt/lists/*; \
	fi; \
	\
# https://julialang.org/downloads/#julia-command-line-version
# https://julialang-s3.julialang.org/bin/checksums/julia-1.6.4.sha256
	arch="$(dpkg --print-architecture)"; \
	case "$arch" in \
		'amd64') \
			url='https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.4-linux-x86_64.tar.gz'; \
			sha256='52244ae47697e8073dfbc9d1251b0422f0dbd1fbe1a41da4b9f7ddf0506b2132'; \
			;; \
		'armhf') \
			url='https://julialang-s3.julialang.org/bin/linux/armv7l/1.6/julia-1.6.4-linux-armv7l.tar.gz'; \
			sha256='9ad3f6bd71eb6840d4cef1569855da20c0b4931a2bdf77703a64df672b0702a1'; \
			;; \
		'arm64') \
			url='https://julialang-s3.julialang.org/bin/linux/aarch64/1.6/julia-1.6.4-linux-aarch64.tar.gz'; \
			sha256='072daac7229c15fa41d0c1b65b8a3d6ee747323d02f5943da3846b075291b48b'; \
			;; \
		'i386') \
			url='https://julialang-s3.julialang.org/bin/linux/x86/1.6/julia-1.6.4-linux-i686.tar.gz'; \
			sha256='9d43d642174cf22cf0fde18dc2578c84f357d2c619b9d846d3a6da4192ba48cf'; \
			;; \
		*) \
			echo >&2 "error: current architecture ($arch) does not have a corresponding Julia binary release"; \
			exit 1; \
			;; \
	esac; \
	\
	curl -fL -o julia.tar.gz.asc "$url.asc"; \
	curl -fL -o julia.tar.gz "$url"; \
	\
	echo "$sha256 *julia.tar.gz" | sha256sum --strict --check -; \
	\
	export GNUPGHOME="$(mktemp -d)"; \
	gpg --batch --keyserver keyserver.ubuntu.com --recv-keys "$JULIA_GPG"; \
	gpg --batch --verify julia.tar.gz.asc julia.tar.gz; \
	command -v gpgconf > /dev/null && gpgconf --kill all; \
	rm -rf "$GNUPGHOME" julia.tar.gz.asc; \
	\
	mkdir "$JULIA_PATH"; \
	tar -xzf julia.tar.gz -C "$JULIA_PATH" --strip-components 1; \
	rm julia.tar.gz; \
	\
	apt-mark auto '.*' > /dev/null; \
	[ -z "$savedAptMark" ] || apt-mark manual $savedAptMark; \
	apt-get purge -y --auto-remove -o APT::AutoRemove::RecommendsImportant=false; \
	\
# smoke test
	julia --version


# **************************** Install tools ***********************************
ENV TOOLS_DIR /opt/tools
ENV DOWNLOAD_DIR /home/Downloads

ENV PLINK_PATH "$TOOLS_DIR/plink"
ENV PLINK2_PATH "$TOOLS_DIR/plink2"
ENV KING_PATH "$TOOLS_DIR/king"
ENV VCFTOOLS_PATH "$TOOLS_DIR/vcftools"
ENV BCFTOOLS_PATH "$TOOLS_DIR/bcftools"

# Add tools to PATH
ENV PATH $PLINK_PATH:$PLINK2_PATH:$KING_PATH:$VCFTOOLS_PATH/bin:$BCFTOOLS_PATH/bin:$PATH

# Setup compression tools
RUN set -eux; \
    apt-get update; \
	apt-get install -y --no-install-recommends \
        unzip \
        git \
        autoconf automake build-essential pkg-config zlib1g-dev \
        libbz2-dev liblzma-dev

RUN set -eux; \
# Setup tools directory
    mkdir -p $TOOLS_DIR; \
    mkdir -p $DOWNLOAD_DIR; \
# Setup Plink 1.9
    echo "Setting up Plink 1.9"; \
    curl -fL -o $DOWNLOAD_DIR/plink1.zip "http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip"; \
    unzip $DOWNLOAD_DIR/plink1.zip -d "$PLINK_PATH"; \
# Setup Plink2
    echo "Setting up Plink 2.0"; \
    curl -fL -o $DOWNLOAD_DIR/plink2.zip "http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_avx2.zip"; \
    unzip $DOWNLOAD_DIR/plink2.zip -d "$PLINK2_PATH"; \
# Setup KING
    echo "Setting up KING"; \
    curl -fL -o $DOWNLOAD_DIR/king.tar.gz "https://www.kingrelatedness.com/Linux-king.tar.gz"; \
    mkdir "$KING_PATH"; \
    tar -xzf $DOWNLOAD_DIR/king.tar.gz -C "$KING_PATH"; \
# Setup vcftools
    echo "Setting up vcftools"; \
    cd "$DOWNLOAD_DIR"; \
    git clone https://github.com/vcftools/vcftools.git; \
    cd vcftools; \
    ./autogen.sh; \
    ./configure --prefix="$VCFTOOLS_PATH"; \
    make; \
    make install; \
    cd -; \
# Setup bcftools
    echo "Setting up bcftools"; \
    cd "$DOWNLOAD_DIR"; \
    curl -fL -o bcftools.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2; \
    mkdir -p bcftools; \
    tar -xf bcftools.tar.bz2 -C "bcftools" --strip-components 1; \
    cd bcftools; \
    ./configure --prefix="$BCFTOOLS_PATH"; \
    make; \
    make install; \
    cd -; \
# Debug
    # ls -la $TOOLS_DIR; \
    # ls -la $PLINK_PATH; \
    # ls -la $PLINK2_PATH; \
    # ls -la $KING_PATH; \
    # ls -la $BCFTOOLS_PATH; \
    # ls -ls $BCFTOOLS_PATH/bin; \
    # ls -ls $BCFTOOLS_PATH/libexec; \
    \
# Test
    plink --version; \
    plink2 --version; \
    which king; \
    vcftools --version; \
    bcftools --version; \
# Clean up
    rm -rf $DOWNLOAD_DIR


# ************************* Add INTERVENE scripts *********************************
WORKDIR /opt/intervene
ENV SCRIPT_DIR "/opt/intervene/scripts"
ENV DATA_DIR "/data"
RUN mkdir -p $DATA_DIR

# Copy source files
WORKDIR $SCRIPT_DIR
COPY . .

# Install Julia packages
RUN set -eux; \
	echo "$DATA_DIR"; \
	julia --project=$SCRIPT_DIR -e "using Pkg; Pkg.instantiate()"

# Setup path for commands
ENV PATH "$SCRIPT_DIR/commands":$PATH
