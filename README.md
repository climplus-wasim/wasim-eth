This repository contains the sources of the hydrologic model system WaSiM-ETH (Water balance Simulation Model) starting with version 11.01.01. Earlier versions are available via https://www.wasim.ch/en/dialog/contact.htm

The model description (manual) is in [docs/wasim_2025_en.pdf](docs/wasim_2025_en.pdf), also available from https://www.wasim.ch/downloads/doku/wasim/wasim_2025_en.pdf.

In /docs, a readme can be found how to compile and run the model on RHEL 9. An annotated sample control file, with references to the manual, is in [docs/sample_control_file.wasim](docs/sample_control_file.wasim). Known limits on model domain size are described in [docs/domain-size-limits.md](docs/domain-size-limits.md), and known problems are listed in [TODO.md](TODO.md).

## Dependencies

Required:

- a C++11 compiler (g++ 11 or newer on Linux, Apple clang on macOS)
- OpenMP (part of g++; on macOS provided by `libomp`)
- netCDF 4.x C library with headers (e.g. 4.6.1). The netCDF C++ interface is bundled in `src/netcdfcpp_Linux` and does not need to be installed separately.
- make

Optional:

- OpenMPI, for the parallel MPI builds (makefiles with `openmpi` / `mpi` in their name). Without it, use one of the non-MPI makefiles, which build an OpenMP-only binary.

Package names for common systems:

| System | Packages |
| --- | --- |
| RHEL / CentOS | `gcc-c++ make netcdf-devel` (EPEL), optional `openmpi-devel` |
| Ubuntu / Debian | `g++ make libnetcdf-dev`, optional `libopenmpi-dev openmpi-bin` |
| openSUSE | `gcc-c++ make netcdf-devel`, optional `openmpi4-devel` |
| macOS (Homebrew) | `xcode-select --install`, then `brew install netcdf libomp`, optional `brew install open-mpi` |

## Compiling

The makefiles are in `src`, one per platform, e.g.:

    cd src
    make -f makefile_rhel9_openmpi

`make clean` in `src` removes the build products of any of the makefiles.

The Linux makefiles expect the netCDF C++ sources in `src/netcdfcpp`; copy or symlink `netcdfcpp_Linux` to that name first (`ln -s netcdfcpp_Linux netcdfcpp`).

### macOS

`src/makefile_macos` works on Apple Silicon and Intel Macs with Apple clang and Homebrew:

    cd src
    make -f makefile_macos          # OpenMP + MPI (needs open-mpi)
    make -f makefile_macos MPI=0    # OpenMP only

By default the binary is tuned for the build machine's CPU family; pass `ARCHFLAGS=` for a portable build. Homebrew paths are found via `brew --prefix`; they can be overridden with `NETCDF=...` and `LIBOMP=...`. To use Homebrew GCC instead of Apple clang, pass e.g. `CXX=g++-15 OMPFLAGS=-fopenmp OMPLIBS=-fopenmp MPI=0`.
