# TODO

Known problems, noted but not yet fixed.

## Build

- [ ] The Linux makefiles list `bodenmod.o` and `infilt.o` in `OBJ`, but `bodenmod.cc` and `infilt.cc` are not in the repository, so `make` should stop with "No rule to make target". The source list in `wasimmpic.vcxproj` and `makefile_macos` does not include them.
- [ ] The non-MPI Linux makefiles (`makefile_ubuntu`, `makefile_rhel9`, …) have a `topmodel` target that needs `wasim.cc` and `modappli.cc`, which are not in the repository either.
- [ ] The Linux makefiles use `-I./netcdfcpp` and `netcdfcpp/*.o`, but the directory is called `src/netcdfcpp_Linux`. Workaround: `ln -s netcdfcpp_Linux netcdfcpp`.
- [ ] The Linux makefiles' `%.o: %.cc` rule compiles `netcdfcpp/ncAtt.cc` to `./ncAtt.o` instead of `netcdfcpp/ncAtt.o` (no `-o $@`), so make never sees those objects as up to date and rebuilds them on every run. The archive step (`ar rcv ... *.o`) depends on this misplacement.
- [ ] There are 21 platform makefiles that differ mostly in include/library paths and CPU flags; they could be merged into one makefile with variables (as in `makefile_macos`).

## Code

- [ ] `delete` used on memory allocated with `new[]` (undefined behaviour), reported by clang `-Wmismatched-new-delete`: `src/unsatzon.cc:14798`, `src/netcdfcpp_Linux/ncGroup.cc:107`.
- [ ] `printf`-style format/argument type mismatches, e.g. `src/unsatzon.cc:10428` (`%d` with a `long`). Build with clang `-Wformat` for the full list.

## Domain size and limits

See [docs/domain-size-limits.md](docs/domain-size-limits.md) for details.

- [ ] netCDF I/O builds flat buffers sized `ncols*nrows` using `int` arithmetic (`src/modefuzr.cc:4063-4104`, `4358-4394`, `6014-6056`; `src/netcdfhlp.cc:880-960`). This overflows above 2,147,483,647 cells per grid (e.g. 46,341 × 46,341), which leads to a crash or heap corruption. Use `size_t`.
- [ ] Zone and region code checks are off by one (`src/defltuzr.cc:1585-1616`, same pattern at `3783ff`): `zzarray` has `maxzonenwert` = 32767 entries, but code 32767 passes the check `zonenid > maxzonenwert`, and `zzone` has `maxzonenzahl` = 16000 entries, but `zonmax == 16000` passes `zonmax > maxzonenzahl`. Both write one element past the end. The comment there also still says 8000 zones.
- [ ] `zonenid = (int)zn[i][j] + 0.5;` truncates instead of rounding (the cast happens before the `+ 0.5`). Same for regions.
- [ ] The netCDF reader stores the `x`/`y` coordinate variables in `int` buffers (`src/netcdfhlp.cc:819`, `864`), so fractional coordinates are truncated. This gives a cellsize of 0 for geographic (degree) grids and wrong results for cells smaller than 1 m. It also reads `ybuf[1]` when `nc_rows == 1`.
- [ ] The netCDF reader stores `add_offset` as `int` (`src/netcdfhlp.cc:785`), so fractional offsets of packed data (e.g. 273.15) are truncated.
- [ ] `memory_handling::matrix()` computes sizes in `long`. On Windows `long` is 32 bit, so grids with more than about 268 million cells (2^31 bytes per grid, e.g. 16,384 × 16,384) overflow there. This is fine on Linux and macOS.
- [ ] MPI grid exchange uses the row number as message tag. The MPI standard only guarantees tags up to 32767; Open MPI allows much more, but other implementations may not. Check `MPI_TAG_UB` if domains have more than 32767 rows.
