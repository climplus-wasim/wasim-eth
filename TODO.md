# TODO

Known problems, noted but not yet fixed. Keep in mind that there is some idle-ish interest in running this logic on GPU. If the performance and architecture arguments keep piling up in favour of a rewrite (likely in [gridtools](https://github.com/GridTools/gridtools), [gt4py](https://github.com/GridTools/gt4py)), the balance will tip sooner. For now, gather painpoints/todos here; it's hard to distinguish between "want" and "need" for now. We can provide tentative labels on feasibility for tasks, though.

## Documentation

- [ ] How/why to run using mpirun/mpiexec
- [ ] translate german comments (not symbols!) to english, but translate/explain when a symbol uses an abbreviation, be it english or german.

## Build

- [ ] The Linux makefiles list `bodenmod.o` and `infilt.o` in `OBJ`, but `bodenmod.cc` and `infilt.cc` are not in the repository, so `make` should stop with "No rule to make target". The source list in `wasimmpic.vcxproj` and `makefile_macos` does not include them.
- [ ] The non-MPI Linux makefiles (`makefile_ubuntu`, `makefile_rhel9`, …) have a `topmodel` target that needs `wasim.cc` and `modappli.cc`, which are not in the repository either.
- [ ] The different makefiles seem largely redundant, find better ways of distinguishing between distros and architectures
- [ ] The netcdf C++ interface headers are copied from `netcdf-cxx` (homebrew formula name, aka `ncxx4`) and should likely be retrieved on the target system using something like `ncxx4-config --includedir`. Given that we're linking to system netCDF anyways, I'm a little confused why we need more than the headers in there. Would likely resolve the points below.
  - [ ] see if we can find a license
  - [ ] see if this is also a feasible path for windows
- [ ] The Linux makefiles use `-I./netcdfcpp` and `netcdfcpp/*.o`, but the directory is called `src/netcdfcpp_Linux`. Workaround: `ln -s netcdfcpp_Linux netcdfcpp`.
- [ ] The Linux makefiles' `%.o: %.cc` rule compiles `netcdfcpp/ncAtt.cc` to `./ncAtt.o` instead of `netcdfcpp/ncAtt.o` (no `-o $@`), so make never sees those objects as up to date and rebuilds them on every run. The archive step (`ar rcv ... *.o`) depends on this misplacement.
- [ ] There are 21 platform makefiles that differ mostly in include/library paths and CPU flags; they could be merged into one makefile with variables (as in `makefile_macos`).
- [ ] consider simply building a `wasim` executable; since topmodel was dropped, the `wasimuzr` moniker is confusing to new users; there only ever is an "Unsaturated Zone model Richards version".
- [ ] the parser is not only vendored in, but actually change from upstream (see `parser.h`) - should we try to replace it with something more performant, or are there no performance critical routines dependent on this?
  - [ ] see if we can find a license

## Code

- [ ] `delete` used on memory allocated with `new[]` (undefined behaviour), reported by clang `-Wmismatched-new-delete`: `src/unsatzon.cc:14798`, `src/netcdfcpp_Linux/ncGroup.cc:107`.
- [ ] `printf`-style format/argument type mismatches, e.g. `src/unsatzon.cc:10428` (`%d` with a `long`). Build with clang `-Wformat` for the full list.
- [ ] consider building a friendly entrypoint printing a friendly/formatted message when calling with `-h` or `--help` or without arguments.
  - [ ] Do not start MPI or OMP processes before parsing the flags, if possible.
- [ ] convert to a consistent utf-8, LF format
- [ ] check for low-hanging performance fruit
- [ ] autoformat using clang-format?

## Domain size and limits

See [docs/domain-size-limits.md](docs/domain-size-limits.md) for details.

- [ ] netCDF I/O builds flat buffers sized `ncols*nrows` using `int` arithmetic (`src/modefuzr.cc:4063-4104`, `4358-4394`, `6014-6056`; `src/netcdfhlp.cc:880-960`). This overflows above 2,147,483,647 cells per grid (e.g. 46,341 × 46,341), which leads to a crash or heap corruption. Use `size_t`.
- [ ] Zone and region code checks are off by one (`src/defltuzr.cc:1585-1616`, same pattern at `3783ff`): `zzarray` has `maxzonenwert` = 32767 entries, but code 32767 passes the check `zonenid > maxzonenwert`, and `zzone` has `maxzonenzahl` = 16000 entries, but `zonmax == 16000` passes `zonmax > maxzonenzahl`. Both write one element past the end. The comment there also still says 8000 zones.
- [ ] `zonenid = (int)zn[i][j] + 0.5;` truncates instead of rounding (the cast happens before the `+ 0.5`). Same for regions.
- [ ] The netCDF reader stores the `x`/`y` coordinate variables in `int` buffers (`src/netcdfhlp.cc:819`, `864`), so fractional coordinates are truncated. This gives a cellsize of 0 for geographic (degree) grids and wrong results for cells smaller than 1 m. It also reads `ybuf[1]` when `nc_rows == 1`.
- [ ] The netCDF reader stores `add_offset` as `int` (`src/netcdfhlp.cc:785`), so fractional offsets of packed data (e.g. 273.15) are truncated.
- [ ] `memory_handling::matrix()` computes sizes in `long`. On Windows `long` is 32 bit, so grids with more than about 268 million cells (2^31 bytes per grid, e.g. 16,384 × 16,384) overflow there. This is fine on Linux and macOS.
- [ ] MPI grid exchange uses the row number as message tag. The MPI standard only guarantees tags up to 32767; Open MPI allows much more, but other implementations may not. Check `MPI_TAG_UB` if domains have more than 32767 rows.
- [ ] Practical limits on writing out grids: is reading/writing handled in the main thread or delegated, can we tune a batch size for our system/memory availability...
