# Spectral Nudging for NEMO

A tool developed for NEMO ocean model to implement spectral nudging in nested simulations.

## Description

This tool implements spectral nudging for regional ocean simulations using NEMO v3.6.
Spectral nudging constrains the large-scales of the inner domain toward large-scale structures of a parent model or reanalysis, while allowing the smaller-scale dynamics to evolve freely.

The tool is applied on the fly during the simulation, and is included as a new key feature with additional requirements in the namelist.

The nudging is applied to temperature and salinity fields at pre defined time-steps in the model, following these steps:

1. Read T&S fields from a parent model or reanalysis
2. Interpolate them onto the model grid
3. Compute the difference between parent and nested model fields
4. Apply a 2D horizontal spatial filter at every vertical level to retain only large-scale differences
5. Apply the correction factor to the model fields weighted by a nudging coefficient

The spatial filter is a 2nd order Butterworth filter applied iteratively
using a gather → filter → scatter approach to handle the MPI parallel decomposition
of NEMO. The filter is currntly applied in serial (see Known Limitations).

The nudging coefficient is composed of a constant component (gamma0), and a spatially varying component provided as a 3d mask (Gamma map).

## Authors

- Renata Tatsch Eidt (CMCC Foundation)
- Anna Katavouta (NOC – National Oceanography Centre)

### Filter parameters

The filter is currently configured with the following hardcoded parameters in
`MY_SRC/sn_simple.F90`:

- **Filter type**: 2nd order zero-phase Butterworth (Wn = 2/7)
- **Number of iterations**: 5
- **Nudging intensity coefficient**: 0.2

These parameters can only be modified directly in the Fortran code at this stage.
Future releases are planned to expose them as namelist parameters for easier tuning.
The script `tools/filter_response.m` can be used to help defining the cutoff
wavelength for a given grid resolution, with the possibility to apply in multiple interations for a stronger smoothing.


## Requirements

- NEMO v3.6

This tool is currently designed and tested with NEMO v3.6.
It can be adapted however to be implemented in more recent versions of NEMO such as v4.2, which is expected to be future release.

## Compilation and Running

### 1. Copy MY_SRC files

Copy all files from `MY_SRC/` into your NEMO configuration `MY_SRC/` directory.

`sn_simple.F90` is the core module implementing the spectral nudging technique. All other files in `MY_SRC/` are standard NEMO v3.6 routines that were minimally modified to integrate the spectral nudging, primarily by adding the `#if defined key_sn_simple` condition to call the relevant routines in the model workflow.


### 2. Add CPP key

Add the following CPP key to your configuration's `cpp_MY_CONFIG.fcm` file:

```
key_sn_simple
```

### 3. Compile NEMO

Compile your NEMO configuration with the new CPP key and MY_SRC files:

```bash
./makenemo -n MY_CONFIG -m MY_ARCH -j 8
```

### 3. Add namelist block

Add the namelist block from `namelists/namelist_sn` to your `namelist_cfg`,
filling in the paths to your input files:

```fortran
!-----------------------------------------------------------------------
&namsn  !   spectral nudging                          ("key_sn_simple")
!-----------------------------------------------------------------------
!              !  file name          ! frequency (hours) ! variable  ! time interp. !  clim  ! 'yearly'/ ! weights               ! rotation ! land/sea mask  !
!              !                     !  (if <0  months)  !   name    !   (logical)  !  (T/F) ! 'monthly' ! filename              ! pairing  ! filename       !
   sr_sal  = 'SALINITY_FILE'         ,        24          , 'so'      ,  .true.      , .false. , 'daily'  , 'WEIGHTS_FILE.nc'    ,   ''     , 'LSM_FILE.nc'
   sr_tem  = 'TEMPERATURE_FILE'      ,        24          , 'thetao'  ,  .true.      , .false. , 'daily'  , 'WEIGHTS_FILE.nc'    ,   ''     , 'LSM_FILE.nc'
   sr_dep  = 'DEPTH_FILE'            ,       -12          , 'gdept_4D',  .false.     , .true.  , 'yearly' , 'WEIGHTS_FILE.nc'    ,   ''     , ''
   sr_msk  = 'MASK_FILE'             ,       -12          , 'tmask'   ,  .false.     , .true.  , 'yearly' , 'WEIGHTS_FILE.nc'    ,   ''     , ''
   !
   cn_gamma_file = 'Gamma_map_3d.nc'   ! 3D nudging coefficient file
   cn_dir        = './'                ! root directory for input files
/
```

### 4. Prepare input files

The following input files are required:

- **Parent model or reanalysis T&S files** — temperature and salinity on the original
  parent model grid and domain. The tool handles the interpolation onto the nested model grid
  on the fly using NEMO's standard `fld_read` function and the weights file
- **Parent model depth file** — depth levels of the parent model grid
- **Parent model mask file** — ocean mask of the parent model grid
- **Gamma map** — 3D NetCDF file containing the spatially-varying nudging
  coefficients on the model grid. It specifies how the strength of the nudging is expected to vary troughout the domain. See `examples/Adriatic/masks/` for an example
- **Weights file** — interpolation weights for regridding input files onto the
  model grid. Can be generated using the WEIGHT tool available in the NEMO
  tools package

## Known Limitations

- **Serial filtering**: the spatial filter gathers the entire global domain into
  1 processor, filters in serial, then redistributes. This can be a computational
  bottleneck for large domains or high grid discretization.
- **NEMO v3.6**: the code has been currently developed and tested with NEMO v3.6 only.
  Future releases are expected to be implemented in NEMO v4.2.
- **Hardcoded filter parameters**: the filter coefficients, constant nudging intensity
  and number of iterations (5) are currently hardcoded in `MY_SRC/sn_simple.F90`
  and require direct code modification to change. Future releases can expose
  these as namelist parameters.

## Example Configuration

The `examples/Adriatic/` directory contains a working configuration
for a regional model of the Adriatic Sea at 1/16° resolution (~6 km).


## Citation

Eidt, R. T., Katavouta, A., Santos da Costa, V., Navarra, A., Verri, G. Spectral nudging in nested ocean models: a marginal sea application for the Adriatic. *Under review.*
