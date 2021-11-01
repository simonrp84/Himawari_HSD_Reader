!******************************************************************************%
! *
! *    Copyright (C) 2016-2019 Simon Proud <simon.proud@physics.ox.ac.uk>
! *    License: CC BY-NC-ND 4.0
! *
! ******************************************************************************/


!*******************************************************************************
! Module containing the Fortran interface to himawari_native_util.
!*******************************************************************************

module himawari

    implicit none

    private

    public :: himawari_t_data, &
              himawari_t_navdata, &
              himawari_t_info, &
              himawari_t_extent, &
              himawari_t_struct

    ! Type kind value
    integer, parameter, public :: ahi_byte  = 1
    integer, parameter, public :: ahi_sint  = 2
    integer, parameter, public :: ahi_lint  = 4
    integer, parameter, public :: ahi_sreal = 4
    integer, parameter, public :: ahi_dreal = 8

    ! Fill values, equivalent to those in ORAC
    integer(kind=ahi_byte),  parameter, public :: him_byte_fill_value  = -127
    integer(kind=ahi_sint),  parameter, public :: him_sint_fill_value  = -32767
    integer(kind=ahi_lint),  parameter, public :: him_lint_fill_value  = -32767
    real(kind=ahi_sreal),    parameter, public :: him_sreal_fill_value = -999.0
    real(kind=ahi_dreal),    parameter, public :: him_dreal_fill_value = -999.0

    !V arious useful parameters for Himawari processing
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_BOUNDS_FULL_DISK = 0
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_BOUNDS_ACTUAL_IMAGE= 1
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_BOUNDS_LINE_COLUMN = 2
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_BOUNDS_LAT_LON  = 3

    integer(kind=ahi_sint), parameter, public :: HIMAWARI_UNIT_CNT = 0
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_UNIT_RAD = 1
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_UNIT_RBT = 2
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_UNIT_BRT = 3

    integer(kind=ahi_sint), parameter, public :: HIMAWARI_MAX_CHANS = 16
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_CHARLEN = 512

    integer(kind=ahi_sint), parameter, public :: HIMAWARI_SUCCESS = 0
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_FAILURE = -1

    integer(kind=ahi_sint), parameter, public :: HIMAWARI_TRUE   = 1
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_FALSE   = 0
    real(kind=ahi_sreal),    parameter, public :: HIMAWARI_PI  = 4 * atan (1.0_8)

    ! Image sizes, assumed to be full disk
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_IR_NLINES = 5500
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_IR_NCOLS = 5500

    integer(kind=ahi_sint), parameter, public :: HIMAWARI_VIS_NLINES = 11000
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_VIS_NCOLS = 11000

    integer(kind=ahi_sint), parameter, public :: HIMAWARI_HVI_NLINES = 22000
    integer(kind=ahi_sint), parameter, public :: HIMAWARI_HVI_NCOLS = 22000

    integer(kind=ahi_sint), parameter, public :: HIMAWARI_NCHANS  = 16

    ! Definitions of physical constants. Used for calculating EBT from RAD
    ! In most cases we will use the values from the header instead
    real(kind=ahi_sreal), parameter, public :: HIMAWARI_lgtspd = 299792458
    real(kind=ahi_sreal), parameter, public :: HIMAWARI_plkcons = 6.62606957E-034
    real(kind=ahi_sreal), parameter, public :: HIMAWARI_bltcons = 1.3806488E-023


    ! Definitions of navigation parameters. Used for transforming line/col to lat/lon
    real(kind=ahi_dreal), parameter, public :: HIMAWARI_DEGTORAD = (4*atan (1.0))/180.0
    real(kind=ahi_dreal), parameter, public :: HIMAWARI_RADTODEG = 180.0/(4*atan (1.0))
    real(kind=ahi_dreal), parameter, public :: HIMAWARI_SCLUNIT = 1.525878906250000e-05


    ! Main data arrays for the image, geoloc etc
    type :: himawari_t_data
        integer(kind=ahi_sint) :: memory_alloc_d
        integer(kind=ahi_sint) :: n_bands
        integer(kind=ahi_sint) :: n_lines
        integer(kind=ahi_sint) :: n_cols
        real(kind=ahi_dreal),    DIMENSION(:,:), ALLOCATABLE :: time(:, :)
        real(kind=ahi_sreal),    DIMENSION(:,:), ALLOCATABLE :: lat(:, :)
        real(kind=ahi_sreal),    DIMENSION(:,:), ALLOCATABLE :: lon(:, :)
        real(kind=ahi_sreal),    DIMENSION(:,:), ALLOCATABLE :: sza(:, :)
        real(kind=ahi_sreal),    DIMENSION(:,:), ALLOCATABLE :: saa(:, :)
        real(kind=ahi_sreal),    DIMENSION(:,:), ALLOCATABLE :: vza(:, :)
        real(kind=ahi_sreal),    DIMENSION(:,:), ALLOCATABLE :: vaa(:, :)
        real(kind=ahi_sreal),    DIMENSION(:,:,:), ALLOCATABLE:: indata(:, :, :)
        real(kind=ahi_sreal),    DIMENSION(:,:,:), ALLOCATABLE:: tmpdata(:, :, :)
        real(kind=ahi_sreal),    DIMENSION(:,:,:), ALLOCATABLE:: soldata(:, :, :)
        real(kind=ahi_sreal),    DIMENSION(:,:), ALLOCATABLE :: cal_slope(:)
    end type himawari_t_data

    ! Navigation parameters derived from the header file
    type :: himawari_t_navdata
        real(kind=ahi_dreal) :: subLon
        real(kind=ahi_dreal) :: cfac
        real(kind=ahi_dreal) :: lfac
        real(kind=ahi_dreal) :: coff
        real(kind=ahi_dreal) :: loff
        real(kind=ahi_dreal) :: satDis
        real(kind=ahi_dreal) :: eqtrRadius
        real(kind=ahi_dreal) :: polrRadius
        real(kind=ahi_dreal) :: projParam1
        real(kind=ahi_dreal) :: projParam2
        real (kind=ahi_dreal) :: projParam3
        real(kind=ahi_dreal) :: projParamSd
    end type himawari_t_navdata

    ! Useful bits and bobs for determining the segment filenames
    type :: himawari_t_info
        character(len=HIMAWARI_CHARLEN) :: indir
        character(len=HIMAWARI_CHARLEN) :: timeslot
        integer :: satnum !101 = H08, 102 = H08
        character(len=HIMAWARI_CHARLEN) :: region
        character(len=HIMAWARI_CHARLEN) :: regseg
    end type himawari_t_info

    ! Image extents
    type :: himawari_t_extent
        integer(kind=ahi_sint) :: x_min
        integer(kind=ahi_sint) :: y_min
        integer(kind=ahi_sint) :: x_max
        integer(kind=ahi_sint) :: y_max
        integer(kind=ahi_sint) :: x_size
        integer(kind=ahi_sint) :: y_size
        logical,dimension(10) :: procseg

        integer :: segdel_ir = 550
        integer :: segdel_vi = 1100
        integer :: segdel_hv = 2200

        integer, dimension(10) :: segpos_ir = (/1, 551,1101,1651,2201,2751, 3301, 3851, 4401, 4951/)
        integer, dimension(10) :: segpos_vi = (/1,1101,2201,3301,4401,5501, 6601, 7701, 8801, 9901/)
        integer, dimension(10) :: segpos_hv = (/1,2201,4401,6601,8801,11001,13201,15401,17601,19801/)

        integer,dimension(3) :: startpos
        integer,dimension(3) :: endpos
    end type himawari_t_extent

    ! Main struct
    type :: himawari_t_struct
        type(himawari_t_data) :: ahi_data
        type(himawari_t_navdata) :: ahi_navdata
        type(himawari_t_info) :: ahi_info
        type(himawari_t_extent) :: ahi_extent
        integer(kind=ahi_sint) :: inchans(HIMAWARI_MAX_CHANS)
        integer(kind=ahi_sint) :: convert(HIMAWARI_MAX_CHANS)
        logical :: do_solar
        logical :: do_solar_angles
        logical :: vis_res
        logical :: archive_struct
    end type himawari_t_struct

end module himawari
