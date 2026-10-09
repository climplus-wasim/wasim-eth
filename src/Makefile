# general makefile: cleanup for all platform builds
#
# to build, use one of the platform makefiles, e.g.
#   make -f makefile_rhel9_openmpi
#   make -f makefile_macos
#
# make clean removes the objects, archives and binaries produced by any of them

.PHONY: help clean

help:
	@echo "build with one of the platform makefiles, e.g.:"
	@echo "  make -f makefile_rhel9_openmpi"
	@echo "  make -f makefile_macos"
	@echo "available makefiles:"
	@ls makefile_* | sed 's/^/  /'
	@echo "make clean removes all build products"

clean:
	rm -f *.o netcdfcpp*/*.o *.a
	rm -f wasimuzr-* wasimtop-*
