# TODO

Known problems, noted but not yet fixed. Keep in mind that there is some idle-ish interest in running this logic on GPU. If the performance and architecture arguments keep piling up in favour of a rewrite (likely in [gridtools](https://github.com/GridTools/gridtools), [gt4py](https://github.com/GridTools/gt4py)), the balance will tip sooner. For now, gather painpoints/todos here; it's hard to distinguish between "want" and "need" for now. We can provide tentative labels on feasibility for tasks, though.

Experiments that should not go onto branches of the main repo (e.g. due to data size) should make use of <https://github.com/mmyrte/wasim-forclim-evoland-experiments>. It currently has git-lfs enabled on zip files (github standard limits).

## Documentation

- [x] How/why to run using mpirun/mpiexec: README, section Running
- [ ] translate german comments (not symbols!) to english, but translate/explain when a symbol uses an abbreviation, be it english or german.
  - [x] glossary of German terms, abbreviations and file extensions: [docs/glossary.md](docs/glossary.md)
  - [ ] translation of the comments themselves (about 90,000 lines of code; best done module by module, together with other changes to that module)
- [x] do we still need docs/RHEL8und9_make_wasim_readme.txt? what can we integrate into the main readme? Integrated into the README (build, MPI run, debugging) and removed; the VirtualBox mount steps were machine-specific.

## Build

- [x] (superseded by `src/Makefile`; the legacy file itself is unchanged) The Linux makefiles list `bodenmod.o` and `infilt.o` in `OBJ`, but `bodenmod.cc` and `infilt.cc` are not in the repository, so `make` should stop with "No rule to make target". The source list in `wasimmpic.vcxproj` and `makefile_macos` does not include them.
- [x] (superseded by `src/Makefile`; the legacy file itself is unchanged) The non-MPI Linux makefiles (`makefile_ubuntu`, `makefile_rhel9`, …) have a `topmodel` target that needs `wasim.cc` and `modappli.cc`, which are not in the repository either.
- [x] The different makefiles seem largely redundant, find better ways of distinguishing between distros and architectures: `src/Makefile` builds on Linux and macOS (OpenMP or MPI, bundled or system netCDF C++), objects in `src/build/`. The `makefile_*` files are kept but superseded; delete them once nobody uses them.
- [ ] The netcdf C++ interface headers are copied from `netcdf-cxx` (homebrew formula name, aka `ncxx4`) and should likely be retrieved on the target system using something like `ncxx4-config --includedir`. Given that we're linking to system netCDF anyways, I'm a little confused why we need more than the headers in there. Would likely resolve the points below.
  - Findings: `src/netcdfcpp_Linux` is an older, unmodified release of [netcdf-cxx4](https://github.com/Unidata/netcdf-cxx4) (C++ API for netCDF-4, not just headers: the `.cc` files are the library itself). Compared with netcdf-cxx4 4.3.1 (Ubuntu 24.04), its headers differ only by a missing `;` in `ncCheck.h`, two `const` return types, and functions added upstream later. `make NCXX=system` builds against the system library; golden tests pass (see below).
  - [ ] see if we can find a license: netcdf-cxx4 is distributed by Unidata/UCAR under its BSD-style licence (the `COPYRIGHT` file in the netcdf-cxx4 repository); the bundled copy contains no licence file. If the bundled copy stays, add that file.
  - [ ] see if this is also a feasible path for windows: vcpkg and conda-forge both package netcdf-cxx4, so the Visual Studio project could use one of them instead of `netcdfcpp_Windows` (not tested)
- [x] (superseded by `src/Makefile`; the legacy file itself is unchanged) The Linux makefiles use `-I./netcdfcpp` and `netcdfcpp/*.o`, but the directory is called `src/netcdfcpp_Linux`. Workaround: `ln -s netcdfcpp_Linux netcdfcpp`.
- [x] (superseded by `src/Makefile`; the legacy file itself is unchanged) The Linux makefiles' `%.o: %.cc` rule compiles `netcdfcpp/ncAtt.cc` to `./ncAtt.o` instead of `netcdfcpp/ncAtt.o` (no `-o $@`), so make never sees those objects as up to date and rebuilds them on every run. The archive step (`ar rcv ... *.o`) depends on this misplacement.
- [x] (superseded by `src/Makefile`; the legacy file itself is unchanged) There are 21 platform makefiles that differ mostly in include/library paths and CPU flags; they could be merged into one makefile with variables (as in `makefile_macos`).
- [x] consider simply building a `wasim` executable (`src/Makefile` builds `wasim` and `wasim-mpi`); since topmodel was dropped, the `wasimuzr` moniker is confusing to new users; there only ever is an "Unsaturated Zone model Richards version".
- [ ] the parser is not only vendored in, but actually change from upstream (see `parser.h`) - should we try to replace it with something more performant, or are there no performance critical routines dependent on this?
  - Findings: `parser.cc/.h` is the expression parser by Jos de Jong (2007, speqmath.com), merged into one file and adapted by J. Schulla in 2011. It is used only by the silting-up model, method 3, which evaluates user-defined expressions from the control file for every cell and time step (`siltingup.cc:688ff`), so it is performance-critical only there. The Swiss setups do not use the silting-up model, and the parser does not appear in their profiles.
  - [ ] see if we can find a license: the files carry no licence; the original was published on speqmath.com (not checked from here whether it still is). Ask the author or J. Schulla.

## Code

- [x] `delete` used on memory allocated with `new[]` (undefined behaviour), reported by clang `-Wmismatched-new-delete`: `src/unsatzon.cc:14798`, `src/netcdfcpp_Linux/ncGroup.cc:107`.
- [x] `printf`-style format/argument type mismatches, e.g. `src/unsatzon.cc:10428` (`%d` with a `long`). Build with clang `-Wformat` for the full list. Fixed all reported by clang 18, including `sscanf("%ld")` into an `int` (`defltuzr.cc:1361`, `1412`, irrigation table, overwrote 4 bytes next to the variable), non-literal format strings in `parser.cc`, and a missing virtual destructor (`dataassim.h`). Golden tests unchanged.
- [ ] Possible logic errors reported by clang, not changed because a fix changes results:
  - `snowmod.cc:5294`: `refTemperSnow[i][j][...];` has no effect, presumably meant `snowTemp = ...`; with layered snow and heat transfer, `UpdateSWEStack` always gets a snow temperature of 0.
  - `meteoint.cc:2250`: `STAT[...] && 0` is always false, so the check for "no station at all" can never fire.
  - The routine that fills zone codes for missing codes appears three times with an ambiguous `else` (`siltingup.cc:257`, `snowmod.cc:395`, `unsatzon.cc:930`): the `else` binds to the inner `if`, but in `unsatzon.cc` it is indented as if it belonged to the outer one. Check which is intended.
- [x] consider building a friendly entrypoint printing a friendly/formatted message when calling with `-h` or `--help` or without arguments.
  - [x] Do not start MPI or OMP processes before parsing the flags, if possible. (`-h`, `--help`, `--version` and no arguments are handled before `MPI_Init`.)
- [x] convert to a consistent utf-8, LF format: all sources UTF-8 (13 files were Latin-1), LF was already used everywhere; `.editorconfig` added. The commit is listed in `.git-blame-ignore-revs` (`git config blame.ignoreRevsFile .git-blame-ignore-revs`).
- [ ] check for low-hanging performance fruit (profile of the Rhone sample, hourly, 1 thread: Richards zone 60%, of which the van Genuchten table lookups `h_calc`, `get_theta`, `k_calc` 36%; IDW interpolation 11%; snow 6%; I/O about 5%)
  - [x] link-time optimization (`-flto`): 3% faster, identical results
  - [ ] the parallel loops use `schedule(runtime)` and the default static schedule leaves threads idle: `OMP_SCHEDULE=dynamic,4` made a 10-day January run 27% faster on 4 threads, but a full hourly year ran about twice as slow as with the static schedule (still to be understood, possibly through the convergence of the groundwater iteration; see below). Measure per module before changing the default.
  - [ ] IDW recomputes all station distances (`sqrt`, power) for every cell and time step; weights only change when stations have missing values, so they could be cached per cell (memory: cells x stations)
  - [ ] the Richards loop allocates about 30 scratch arrays per grid row and time step (`unsatzon.cc:8130ff`); allocate once per thread
- [ ] autoformat using clang-format? A proposed `.clang-format` close to the current style (tabs, braces on their own line) is in the repository root. Applied to all sources it would change about 40,000 of 91,000 lines, burying `git blame` and conflicting with any upstream change. Suggestion: format only changed lines (`git clang-format` before committing), and reformat whole files only when a module is rewritten.

## Unit/integration tests

- [ ] Unit tests to ensure correctness and performance experiments
- [x] golden tests: [tests/golden](tests/golden/README.md), cases Doubs and Rhone (data in wasim-forclim-evoland-experiments `golden/`)
- [ ] integration or golden tests for interface modernization
  - [ ] control file
  - [ ] wrapper scripts for control file generation, modern/standardized string interpolation? so far Ahmed used regexes to manipulate the files; path construction is fickle
  - [ ] input data formats (geotiff?)
- [ ] suitable abstraction integration tests for any experiments to translate to `gt4py`...
  - [x] method used for the T-index snow prototype: a WaSiM run writes a process's inputs and outputs at every step, the port is driven by the same inputs and compared per step (wasim-forclim-evoland-experiments `gt4py/`)
  - [x] grid dump hook `griddump` in `src/modapuzr.cc`: with `WASIM_DUMP_GRIDS` (grid names, comma-separated) and `WASIM_DUMP_DIR` set, named grids are written as raw float64 at fixed points of the time step (`start`, `after_meteo`, `after_etp`, `after_interception`, `after_snowevap`, `after_unsatzon`, …; `WASIM_DUMP_POINTS` restricts them; `identifications.txt` lists all grid names). This isolates one process even when other modules change the same grids in place. Inactive without the variables; golden tests unchanged. Used for the interception, snow evaporation, radiation correction and potential ET ports (experiments repo `gt4py/coupling.md`).

## Domain size and limits

See [docs/domain-size-limits.md](docs/domain-size-limits.md) for details.

- [x] netCDF I/O builds flat buffers sized `ncols*nrows` using `int` arithmetic (`src/modefuzr.cc:4063-4104`, `4358-4394`, `6014-6056`; `src/netcdfhlp.cc:880-960`). This overflows above 2,147,483,647 cells per grid (e.g. 46,341 × 46,341), which leads to a crash or heap corruption. Use `size_t`.
- [x] Zone and region code checks are off by one (`src/defltuzr.cc:1585-1616`, same pattern at `3783ff`): `zzarray` has `maxzonenwert` = 32767 entries, but code 32767 passes the check `zonenid > maxzonenwert`, and `zzone` has `maxzonenzahl` = 16000 entries, but `zonmax == 16000` passes `zonmax > maxzonenzahl`. Both write one element past the end. The comment there also still says 8000 zones.
- [x] `zonenid = (int)zn[i][j] + 0.5;` truncates instead of rounding (the cast happens before the `+ 0.5`). Same for regions.
- [x] The netCDF reader stores the `x`/`y` coordinate variables in `int` buffers (`src/netcdfhlp.cc:819`, `864`), so fractional coordinates are truncated. This gives a cellsize of 0 for geographic (degree) grids and wrong results for cells smaller than 1 m. It also reads `ybuf[1]` when `nc_rows == 1`.
- [x] The netCDF reader stores `add_offset` as `int` (`src/netcdfhlp.cc:785`), so fractional offsets of packed data (e.g. 273.15) are truncated.
  - The five items above are fixed (`size_t` sizes, one more array element for zone and region codes, rounding, `double` coordinates and offset, a guard for single-row grids). Golden tests unchanged; a Doubs run reading its precipitation as netCDF grids (method 9) gives identical results before and after. Grids with fractional coordinates, packed data and single rows are not covered by a test yet.
- [ ] `memory_handling::matrix()` computes sizes in `long`. On Windows `long` is 32 bit, so grids with more than about 268 million cells (2^31 bytes per grid, e.g. 16,384 × 16,384) overflow there. This is fine on Linux and macOS.
- [ ] MPI grid exchange uses the row number as message tag. The MPI standard only guarantees tags up to 32767; Open MPI allows much more, but other implementations may not. Check `MPI_TAG_UB` if domains have more than 32767 rows.
- [ ] Practical limits on writing out grids: is reading/writing handled in the main thread or delegated, can we tune a batch size for our system/memory availability...
  - Findings: all reading and writing happens synchronously in the main time loop (no I/O threads). Time series are read through a 1 MB buffer per file (`FileBuffer`), statistics files are buffered and written in blocks, grids are written when their write code says so. In the Rhone profile reading and writing take about 5% of CPU time and wall time equals CPU time, so asynchronous or batched I/O would gain at most a few percent there. It matters only when many grids are written every time step (write codes 1 and 2).

## Found while setting up the golden tests and benchmarks

- [ ] Relative `DefaultOutputDirectory` with multi-layer land use: statistics files of vegetation layers are created as `<dir>/<dir>/name_Layer_n.dat` and the run stops. `zonalstatistic` prepends the directory to the file name in place, and the layer file names are then derived from the already prefixed name (`etpetr.cc:878ff`, `interzep.cc:408ff`). Absolute paths are not prefixed again, so they work.
- [ ] A symbol followed by a single `/` (e.g. `$set $x = $outroot/`) makes the symbol expansion loop forever at "compiling symbols in control file" instead of stopping with an error (manual 4.3.3 says it stops).
- [ ] `$set $a = $b//` drops the trailing `/` of the value, so `$b//` cannot be used to end a path.
- [ ] Results depend on the number of OpenMP threads and on the schedule, beyond rounding: for the Rhone sample (10 days) about 30 of 46,000 values of `sb1_`, `uFCmm` and `qifl` differ by up to 1.0 between 1 and 4 threads with a dynamic schedule (static: up to 0.003). Suspects: the row-wise parallel groundwater iteration (`gwflow.cc:1392ff`, rows read neighbouring rows while they are updated) and shared state in the soil loop. Single-threaded runs are bit-identical.
- [ ] `-march=native` and `-march=haswell` (legacy Linux makefiles) change results through fused multiply-add, e.g. by up to 100 in a reservoir balance file of the Rhone sample; `-ffp-contract=off` would keep them reproducible across machines.
- [x] Overlapping `sprintf` (`sprintf(x, "%s...", x, ...)`, undefined behaviour) in 19 places. With g++ at -O3 the netCDF file names became just `.nc`, so every netCDF grid output failed; clang builds happened to work. Fixed by appending with `sprintf(x + strlen(x), ...)`; g++ and clang builds now give results identical to the reference. `-Wrestrict` (g++) finds new cases.
- [ ] MPI runs differ from serial runs, e.g. for Doubs on 2 ranks 0 instead of nodata in some grid cells and different balance grids (already before the changes above). Check whether the MPI output gathering writes nodata correctly.
- [ ] MPI decomposition depends on a timing benchmark at start-up (`defltuzr.cc:8391ff`), so it can change between runs; `MpiSpeedFactor.dat` in the working directory fixes it. Document or replace with a deterministic default.

## Found while porting to wasim-gt

wasim-gt reproduces each of these bit for bit, so a fix changes results and WaSiM and wasim-gt have to change together. The workarounds wasim-gt needs for that are listed in wasim-gt `TODO.md`. Sources: wasim-gt#19 (snow), #20 (meteo), #21 (permafrost), #22 (groundwater), #23 and #24 (glacier), #25 (radiation correction, ET, phenology), and the experiments repository's `gt4py/coupling.md`.

### Throughout

- [ ] Two values of pi: `teilmodellklasse::pi = 3.141592653` (`teilmod.cc`) and the macro `Pi = 3.14159265359` (`modconst.h:45`), both used in the same functions (e.g. `getsuncoordinates`); results depend on them and on how degree conversions are grouped (about 1e-10). Use `M_PI` (or `std::numbers::pi`).
- [ ] Calendar handling is inconsistent:
  - two leap-year rules, both wrong for 2000: `LUSEtabletype::GetJulDay` uses `year%4==0 && year%400!=0` (`defltuzr.cc:10021`; 2100 counts as a leap year), defaultini uses `year%4==0 && year%100!=0` on the run's end year (`defltuzr.cc:873`);
  - `teilmodellklasse` (sun geometry, time equation) always counts February with 28 days;
  - the step at 24:00 already carries the next day's date (wasim-gt#23).

### Meteo interpolation and precipitation correction

- [ ] Meteo: the regional regression is switched off in the constructor (`if ((1==2)&&(rarray))`, `meteoint.cc:3416`), so region and weight grids are read but ignored and one regression is fitted to all stations.
- [ ] Meteo: the regression fit counts stations with `value > -missing`, IDW with `value >= -missing` (`meteoint.cc:1189` vs `1551`, `1649`).
- [ ] Meteo, method 11: when all stations are missing, the old value is kept but `NewGrid` is not cleared, so scaling and clipping can be applied again (`meteoint.cc:3030-3065`; method 10 clears it at `2950`).
- [ ] Meteo: the precipitation profile refit changes the vertical coefficients and the inversion boundary but keeps the raw residual-plane coefficients (`meteoint.cc:4496-4660`).
- [ ] Meteo: a station at zero distance from a cell gets a squared distance of 1 m², so the other stations still count (`meteoint.cc:1556`, `1654`).

### Radiation correction and potential evapotranspiration

- [ ] Potential ET (`etpetr.cc`, `verdunstungsmodellklasse::potentiell`, line 2042): `LAIkum`, the LAI accumulated above a vegetation layer, adds `LAI[lay-1]` whenever layer `lay` has a land use, even where layer `lay-1` has none, i.e. adds -9999. The aerodynamic resistance `ra = ra_local + LAIscale*LAIkum` of the lower layer then becomes negative. Affects only cells whose top layer is empty and a lower one is not; Doubs and Rhone have none. Found while porting to GT4Py (wasim-gt#2).
- [ ] Potential ET corrects the global radiation in place, `RG = RGEX = RG*(1+(1-cr0)*ssd*(EXKORR-1))` (`etpetr.cc:1956`), but caps only `RGEX` at 1260 Wh/m² per bright hour (line 1973). Modules running later (snow) can see a corrected global radiation above the cap.
- [ ] Minor, in the same function: the constant `c64_0` holds 25.0 (`etpetr.cc:1349`, ratio of measurement height and roughness length for z0 > 2 m); `SSDTemp = 0.5` (line 1854) is unreachable, because `ssdFlag == 1` implies `ssdCnt >= 1`; the daylight test uses `Anteil_hell > 0.05` for the sunshine bookkeeping (line 1815) but `Anteil_hell < 0.05` for `nostrl` (`teilmod.cc:495`), so at exactly 0.05 the step is dark for the bookkeeping but gets radiation.

### Vegetation parameters and phenology

- [ ] Phenology, initialisation and daily update differ (`defltuzr.cc` around 2701 and `int_tvar_params`, Richards version): rsc and rse are clamped at 0 only at initialisation, z0 is raised to 0.01 below 0.01 at initialisation but only at or below 0 in the update, and the initialisation uses the elevation model before the lake model lowers lake cells (lake initialisation method 2), the update the lowered one.
- [ ] Phenology: before the first fix point of the year, values are interpolated from the last fix point minus 366 days, but without its altitude shift (`AltitudeCorrection`: `daystemp[0] = days[n] - 366`), so for cells away from 400 m the parameters jump on the day after the (shifted) last fix point.
- [ ] `LUSEtabletype::GetRsc` (`defltuzr.cc:10214`) tests `GetMethod()==4` twice where the other getters test 3, 4 and 6. Harmless today: DP4 (method 6) always takes the MultipleStartTimes branch before it.
- [ ] Dynamic phenology: the state grids (forcing and chilling units, F*, start of growth) are created filled with zeros (and -9999 where a layer has no land use), so the initialisation rules for unset values (`< 0`: start of growth from `maxStartJDforDP1`, F* from the table, DP1 fallback of DP2/DP4 while no chilling grid is read; `defltuzr.cc` around 2675) never apply to a fresh run. The start of growth stays 0 and F* stays 0 until the first `JDReset_TStart`, and DP2/DP4 switch off their DP1 fallback at once, unlike the manual (2.23.2, 2.23.4) says.
- [ ] Dynamic phenology: the day number `tagnr` uses defaultini's month table, whose February has 29 days in every year of a run whose end year is a leap year (`defltuzr.cc:873`), and 28 otherwise; the parameter interpolation uses `GetJulDay` with another rule.
- [ ] Dynamic phenology: the MultipleStartTimes are initialised with the altitude shift of the fix points but reset at `JDReset_TStart` (and for DP3 on day 1) without it; the initialisation also uses the elevation model before the lake model lowers lake cells. The start of leaf fall and its cold-day counter (`StartDayOfLeaveFall`, level 0 of `TemperatureMOVAVG7`) are one grid per cell, shared by all vegetation layers.

### Snow, snow evaporation and permafrost

- [ ] Snow, daily steps only (`InV > 720`), relevant for daily runs of the Swiss setup:
  - `run1` calls `calc_Accu_and_Melt` twice, for the night and the day fraction. Each call applies the wind redistribution `P += DEPOINDEX * P_SNOW` to the input grid in place, so the second call redistributes the already redistributed precipitation (`snowmod.cc:2863`, `4282-4286`). For DEPOINDEX = 0.3 the factor approaches 1.3 * 1.3 instead of 1.3, also for the water passed on to the soil.
  - In the sublimation case of the melt branch, `Mass_lat` is multiplied by `Anteil` twice (`snowmod.cc:3341`, `3410`), so sublimation shrinks with the square of the day or night fraction.
- [ ] Snow, melt branch: only negative `Mass_lat` (sublimation) is taken from the snow; condensation on melting snow (positive `Mass_lat`) is dropped (`snowmod.cc:3423`). Intended?
- [ ] Snow, refreezing branch: after the surface-temperature iteration the surface temperature can be limited (air - 13 K, or 0 degC), and only the latent heat flux is recomputed; refreezing still uses the energy balance of the unlimited temperature (`snowmod.cc:3087-3102`, `3186`).
- [ ] Snow evaporation adds each cell's EVRS to the zonal statistics twice (`etpetr.cc`, `verdunstungsmodellklasse::snow`, once in each of its two cell loops), and reduces the potential transpiration per vegetation layer but not the layer total `pot_transpiration0` (interception reduces both).
- [ ] Permafrost: WaSiM stops with error 33 unless `SnowCoverDaysGrid` is given as a standard grid. When it creates the counter grids itself, the second check tests `_SnowFreeDaysGrid_` again instead of `_SnowCoverDaysGrid_` (`snowmod.cc:2548`), so after creating the snow-free grid it looks up a snow-cover grid that does not exist. The `[permafrost]` sections of the Doubs and Rhone setups and of `docs/sample_control_file.wasim` also have 3 entries, while snow and soil model read 5 (statistics file and reduced k_sat, `snowmod.cc:2587-2589`, `unsatzon.cc:1406-1410`).

### Glacier

- [ ] Glacier: `Aice` and `Asnow` count the day of the year as `tag + (month-1) * mon[month]` (the loop adds `mon[(int)mt]` instead of `mon[i]`, `snowmod.cc:2660-2674`), so the seasonal sine of the radiation coefficients is shifted by up to a few days (1 July: day 187 instead of 182).
- [ ] Glacier dynamics: `CreateNewGlacierCells` returns at once for every existing glacier (`if (glacierindexmap.find(id) != glacierindexmap.end()) return;`, `snowmod.cc:6548`), so glaciers never grow into new cells; for an unknown id the next line dereferences `end()`.
- [ ] Glacier dynamics: `UpdateGlacierStatisticsZones` passes `isIce = false` to `addCell` (`snowmod.cc:6290`), so each glacier's elevation range is the elevation of its first cell and the elevation bands used for growing and shrinking collapse. Every added cell also extends the old area and volume again (`:6293-6306`).
- [ ] Glacier dynamics: `elevationbelt::addChangesToMarginCells` subtracts the negative value of a cell that drops below 0 from the remaining change (`remainingChanges -= temp`, `snowmod.cc:7208`) instead of the part of the change it absorbed.

### Soil hydraulic functions

- [ ] `SOILtabletype::get_theta` reads `h_tab[layer][0]` for suctions <= 0, but `LookUpTableGeneration` never sets index 0 (uninitialized memory; the result is θs as long as it happens to be 0). `k_calc` searches the k table with `anz_h` instead of `anz_k` (harmless while both are 30).

### Unsaturated zone (wasim-gt#7)

- [x] `hydraulic_heads_iterative` computes the old-state suction `h_temp1[k]` only in the first iteration (`z == 0`) and reuses it afterwards; the buffer is allocated per grid row. When the active depth grows in a later iteration, levels below the first iteration's depth are read from whatever the previous column of the row (or the row's initial zeros) left there (`unsatzon.cc:5239`, `5247`). In the Rhone 2021 year run this changes hydraulic heads in 11 cells on 14 days and 195 of 1213 output files (discharge at gauge 1 by up to 0.04 m³/s); reversing the column order within rows changes 214 output files, so the result depends on the processing order. Doubs 2021 and the golden runs never reach it. Fixed in this branch: levels the cell has not yet initialized in the current split are computed from its own `theta_sl_old`. The fix leaves the 10-day goldens unchanged and is order-independent (reversed column order: identical outputs).
- [ ] A soil-call water budget built from the terms of the debug balance leaves a remainder equal to minus the reported `GWEX` in about 7,800 calls of the Rhone 2021 run (largest 0.115 mm), with and without the fix above. Not yet clear whether this is an accounting boundary of the budget or a real loss; needs a coupled soil–groundwater–river budget.

### Groundwater (wasim-gt#22)

- [ ] The thickness of the first aquifer is overwritten with the total soil thickness at start-up (`gwflow.cc:568-584`).
- [ ] The soil model changes the storage coefficient S0 in every sub step (`unsatzon.cc:1629-1633`, `10888-10900`), so S0 is a state shared by two modules rather than a parameter.
- [ ] `gw_head_old_old_1` is registered with the pointer of `gw_head_old_1` (`gwflow.cc:789-791`), and the iteration grid reuses that name (`815-818`), so the second-to-last head cannot be written out or inspected.
- [ ] Asymmetries: the x transmissivity tests `> upper + 0.01` for the cell itself but `>=` for its neighbour and in y (`gwflow.cc:1172`, `1194`, `1272`); the head update mixes neighbours in north, east, south, west order but sums the transmissivities as north, south, east, west (`1713-1725`).
- [ ] Convergence is measured before the SOR step and divided by 5 in confined conditions (`gwflow.cc:1742-1782`).
- [ ] The balance adds the fluxes to the old water content and subtracts it again instead of summing the fluxes, and leaves out fixed fluxes (`gwflow.cc:2483-2495`).

## Setup

The setups we have inherited for the Swiss domain (hydrological Switzerland) are distributed across 17 basins, necessitating either a harness to run the catchments succesively and in concert, or a singular model that runs across the full domain.
