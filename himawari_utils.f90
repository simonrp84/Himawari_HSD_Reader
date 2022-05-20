!******************************************************************************%
! *
! *    Copyright (C) 2016-2019 Simon Proud <simon.proud@physics.ox.ac.uk>
! *    License: CC BY-NC-ND 4.0
! *
! ******************************************************************************/


module himawari_utils

	use himawari
	implicit none

	public :: AHI_get_timeslot, &
				AHI_get_indir, &
				AHI_free_vals, &
				AHI_alloc_vals, &
				AHI_alloc_vals_data, &
				AHI_get_file_name, &
				AHI_file_exists

contains

integer function AHI_get_timeslot(filename,timeslot) result(status)

	character(len=*), intent(in) :: filename
	character(len=*), intent(out) :: timeslot

	integer pos

	pos	 = index(trim(filename),".DAT")
	if (pos .le. 0) then
		status = HIMAWARI_FAILURE
		return
	endif

	timeslot(1:8) = filename(pos-32:pos-32+8)
	timeslot(9:12) = filename(pos-23:pos-23+4)

	status = HIMAWARI_SUCCESS
	return

end function AHI_get_timeslot


integer function AHI_get_indir(filename,indir) result(status)

	character(len=*), intent(in) :: filename
	character(len=*), intent(out) :: indir

	integer pos

	pos	 = index(trim(filename),"HS_H08_")
	if (pos .le. 0) then
		pos	 = index(trim(filename),"HS_H09_")
	endif
	if (pos .le. 0) then
		status = HIMAWARI_FAILURE
		return
	endif

	indir = filename(1:pos-1)

	status = HIMAWARI_SUCCESS
	return

end function AHI_get_indir

integer function AHI_free_vals(ahi_main) result(status)

	type(himawari_t_struct), intent(inout) :: ahi_main

	if (allocated(ahi_main%ahi_data%lat)) deallocate(ahi_main%ahi_data%lat)
	if (allocated(ahi_main%ahi_data%lon)) deallocate(ahi_main%ahi_data%lon)
	if (allocated(ahi_main%ahi_data%vza)) deallocate(ahi_main%ahi_data%vza)
	if (allocated(ahi_main%ahi_data%vaa)) deallocate(ahi_main%ahi_data%vaa)
	if (allocated(ahi_main%ahi_data%sza)) deallocate(ahi_main%ahi_data%sza)
	if (allocated(ahi_main%ahi_data%saa)) deallocate(ahi_main%ahi_data%saa)
	if (allocated(ahi_main%ahi_data%time)) deallocate(ahi_main%ahi_data%time)
	if (allocated(ahi_main%ahi_data%indata)) deallocate(ahi_main%ahi_data%indata)
	if (allocated(ahi_main%ahi_data%tmpdata)) deallocate(ahi_main%ahi_data%tmpdata)
	if (allocated(ahi_main%ahi_data%soldata)) deallocate(ahi_main%ahi_data%soldata)
	if (allocated(ahi_main%ahi_data%cal_slope)) deallocate(ahi_main%ahi_data%cal_slope)

	status = HIMAWARI_SUCCESS
	return

end function AHI_free_vals

integer function AHI_free_vals_data(ahi_data,verbose) result(status)

	type(himawari_t_data), intent(inout) :: ahi_data
	logical, intent(in) :: verbose

	if (allocated(ahi_data%lat)) deallocate(ahi_data%lat)
	if (allocated(ahi_data%lon)) deallocate(ahi_data%lon)
	if (allocated(ahi_data%vza)) deallocate(ahi_data%vza)
	if (allocated(ahi_data%sza)) deallocate(ahi_data%vaa)
	if (allocated(ahi_data%sza)) deallocate(ahi_data%sza)
	if (allocated(ahi_data%saa)) deallocate(ahi_data%saa)
	if (allocated(ahi_data%time)) deallocate(ahi_data%time)
	if (allocated(ahi_data%indata))	deallocate(ahi_data%indata)
	if (allocated(ahi_data%soldata)) deallocate(ahi_data%soldata)
	if (allocated(ahi_data%tmpdata)) deallocate(ahi_data%tmpdata)
	if (allocated(ahi_data%cal_slope)) deallocate(ahi_data%cal_slope)

	status = HIMAWARI_SUCCESS
	return

end function AHI_free_vals_data

integer function AHI_alloc_vals_data(ahi_data,ahi_extent,nchans,do_solar,verbose) result(status)

	type(himawari_t_data), intent(inout) :: ahi_data
	type(himawari_t_extent), intent(inout) :: ahi_extent
	integer, intent(in) :: nchans
	logical, intent(in) :: do_solar
	logical, intent(in) :: verbose
	integer :: i,length

	if (verbose) write(*,*)"Number of bands to read:",nchans

	if (allocated(ahi_data%lat) .neqv. .true.)		allocate(ahi_data%lat(ahi_extent%x_size,ahi_extent%y_size))
	ahi_data%lat(:,:) = him_sreal_fill_value
	if (allocated(ahi_data%lon) .neqv. .true.)		allocate(ahi_data%lon(ahi_extent%x_size,ahi_extent%y_size))
	ahi_data%lon(:,:) = him_sreal_fill_value
	if (allocated(ahi_data%vza) .neqv. .true.)		allocate(ahi_data%vza(ahi_extent%x_size,ahi_extent%y_size))
	ahi_data%vza(:,:) = him_sreal_fill_value
	if (allocated(ahi_data%vaa) .neqv. .true.)		allocate(ahi_data%vaa(ahi_extent%x_size,ahi_extent%y_size))
	ahi_data%vaa(:,:) = him_sreal_fill_value
	if (allocated(ahi_data%sza) .neqv. .true.)		allocate(ahi_data%sza(ahi_extent%x_size,ahi_extent%y_size))
	ahi_data%sza(:,:) = him_sreal_fill_value
	if (allocated(ahi_data%saa) .neqv. .true.)		allocate(ahi_data%saa(ahi_extent%x_size,ahi_extent%y_size))
	ahi_data%saa(:,:) = him_sreal_fill_value
	if (allocated(ahi_data%time) .neqv. .true.)		allocate(ahi_data%time(ahi_extent%x_size,ahi_extent%y_size))
	ahi_data%time(:,:) = him_sreal_fill_value

	if (allocated(ahi_data%indata) .neqv. .true.)		allocate(ahi_data%indata(ahi_extent%x_size,ahi_extent%y_size,nchans))
	ahi_data%indata(:,:,:) = him_sreal_fill_value

	if (allocated(ahi_data%cal_slope) .neqv. .true.)		allocate(ahi_data%cal_slope(nchans))
	ahi_data%cal_slope(:) = him_sreal_fill_value

	if (do_solar .eqv. .true.) then
		if (allocated(ahi_data%soldata) .neqv. .true.)		allocate(ahi_data%soldata(ahi_extent%x_size,ahi_extent%y_size,3))
		ahi_data%soldata(:,:,:) = him_sint_fill_value
		if (allocated(ahi_data%tmpdata) .neqv. .true.)		allocate(ahi_data%tmpdata(ahi_extent%x_size,ahi_extent%y_size,3))
		ahi_data%tmpdata(:,:,:) = him_sreal_fill_value
	endif

	status = HIMAWARI_SUCCESS
	return

end function AHI_alloc_vals_data

integer function AHI_alloc_vals(ahi_main,ahi_extent,verbose) result(status)

	type(himawari_t_struct), intent(inout) :: ahi_main
	type(himawari_t_extent), intent(inout) :: ahi_extent
	logical, intent(in) :: verbose
	integer :: i,length

!	ahi_main%ahi_data%n_bands = 0
!	do i=1,16
!		if (ahi_main%inchans(i) .eq. 1) then
!			ahi_main%ahi_data%n_bands = ahi_main%ahi_data%n_bands+1
!		endif
!	end do

	if (ahi_main%ahi_data%n_bands <= 0) then
		write(*,*)"No bands are selected!"
		status = HIMAWARI_FAILURE
		return
	endif

	if (verbose) then
		write(*,*)"Number of bands to read:",ahi_main%ahi_data%n_bands
	endif

	if (allocated(ahi_main%ahi_data%lat) .neqv. .true.) allocate(ahi_main%ahi_data%lat(ahi_extent%x_size,ahi_extent%y_size))
	ahi_main%ahi_data%lat(:,:) = him_sreal_fill_value
	if (allocated(ahi_main%ahi_data%lon) .neqv. .true.) allocate(ahi_main%ahi_data%lon(ahi_extent%x_size,ahi_extent%y_size))
	ahi_main%ahi_data%lon(:,:) = him_sreal_fill_value
	if (allocated(ahi_main%ahi_data%vza) .neqv. .true.) allocate(ahi_main%ahi_data%vza(ahi_extent%x_size,ahi_extent%y_size))
	ahi_main%ahi_data%vza(:,:) = him_sreal_fill_value
	if (allocated(ahi_main%ahi_data%vaa) .neqv. .true.) allocate(ahi_main%ahi_data%vaa(ahi_extent%x_size,ahi_extent%y_size))
	ahi_main%ahi_data%vaa(:,:) = him_sreal_fill_value
	if (allocated(ahi_main%ahi_data%sza) .neqv. .true.) allocate(ahi_main%ahi_data%sza(ahi_extent%x_size,ahi_extent%y_size))
	ahi_main%ahi_data%sza(:,:) = him_sreal_fill_value
	if (allocated(ahi_main%ahi_data%saa) .neqv. .true.) allocate(ahi_main%ahi_data%saa(ahi_extent%x_size,ahi_extent%y_size))
	ahi_main%ahi_data%saa(:,:) = him_sreal_fill_value
	if (allocated(ahi_main%ahi_data%time) .neqv. .true.) allocate(ahi_main%ahi_data%time(ahi_extent%x_size,ahi_extent%y_size))
	ahi_main%ahi_data%time(:,:) = him_sreal_fill_value

	if (allocated(ahi_main%ahi_data%indata) .neqv. .true.) allocate(ahi_main%ahi_data%indata(ahi_extent%x_size,ahi_extent%y_size,ahi_main%ahi_data%n_bands))
	ahi_main%ahi_data%indata(:,:,:) = him_sreal_fill_value

	if (allocated(ahi_main%ahi_data%cal_slope) .neqv. .true.) allocate(ahi_main%ahi_data%cal_slope(ahi_main%ahi_data%n_bands))
	ahi_main%ahi_data%cal_slope(:) = him_sreal_fill_value

	if (ahi_main%do_solar .eqv. .true.) then
		if (allocated(ahi_main%ahi_data%soldata) .neqv. .true.) allocate(ahi_main%ahi_data%soldata(ahi_extent%x_size,ahi_extent%y_size,3))
		ahi_main%ahi_data%soldata(:,:,:) = him_sint_fill_value
		if (allocated(ahi_main%ahi_data%tmpdata) .neqv. .true.) allocate(ahi_main%ahi_data%tmpdata(ahi_extent%x_size,ahi_extent%y_size,3))
		ahi_main%ahi_data%tmpdata(:,:,:) = him_sreal_fill_value
	endif

	status = HIMAWARI_SUCCESS
	return

end function AHI_alloc_vals

integer function AHI_get_file_name(cnum, timeslot, satnum, indir, outfile,verbose) result(status)

	integer, intent(in) :: cnum
	character(len=*), intent(in) :: timeslot
	integer, intent(in) :: satnum
	character(len=*), intent(in) :: indir
	character(len=*), intent(out) :: outfile
	logical, intent(in) :: verbose

	if (satnum .eq. 101) then
		outfile = trim(outfile)//trim("himawari8")
	elseif (satnum .eq. 102) then
		outfile = trim(outfile)//trim("himawari9")
	else
		write(*,*) "ERROR: Unsupported platform: ",satnum
		stop
	endif
	outfile=indir
	outfile = trim(outfile)//trim("_ahi_le1b_")
	select case(cnum)
		case(1)
			outfile = trim(outfile)//trim("b01")
		case(2)
			outfile = trim(outfile)//trim("b02")
		case(3)
			outfile = trim(outfile)//trim("b03")
		case(4)
			outfile = trim(outfile)//trim("b04")
		case(5)
			outfile = trim(outfile)//trim("b05")
		case(6)
			outfile = trim(outfile)//trim("b06")
		case(7)
			outfile = trim(outfile)//trim("b07")
		case(8)
			outfile = trim(outfile)//trim("b08")
		case(9)
			outfile = trim(outfile)//trim("b09")
		case(10)
			outfile = trim(outfile)//trim("b10")
		case(11)
			outfile = trim(outfile)//trim("b11")
		case(12)
			outfile = trim(outfile)//trim("b12")
		case(13)
			outfile = trim(outfile)//trim("b13")
		case(14)
			outfile = trim(outfile)//trim("b14")
		case(15)
			outfile = trim(outfile)//trim("b15")
		case(16)
			outfile = trim(outfile)//trim("b16")
		case default
			status = HIMAWARI_FAILURE
			return
	end select
	outfile = trim(outfile)//trim("_org_f_")
	outfile = trim(outfile)//trim(timeslot)
	outfile = trim(outfile)//trim(".bin")

	status = HIMAWARI_SUCCESS
	return

end function AHI_get_file_name


integer function AHI_get_file_name_seg(cnum, seg, timeslot, satnum, indir, outfile, single_seg, verbose) result(status)

	integer, intent(in) :: cnum
	integer, intent(in) :: seg
	character(len=*), intent(in) :: timeslot
	integer, intent(in) :: satnum
	character(len=*), intent(in) :: indir

	character(len=*), intent(out) :: outfile
	logical, intent(in) :: single_seg
	logical, intent(in) :: verbose
	character(len=3) :: tstr
	character(len=32) :: tmpstr
	logical :: file_exists
	integer :: idx
	
	write(tmpstr,'(A,I2.2)') 'B',cnum

	outfile=indir
	if (satnum .eq. 101) then
		outfile = trim(outfile)//trim("HS_H08_")
	elseif (satnum .eq. 102) then
		outfile = trim(outfile)//trim("HS_H09_")
	elseif (satnum .eq. 100) then
	    outfile = trim(outfile)//trim(tmpstr)//"/"//trim("HS_H08_")
	else	
		write(*,*) "ERROR: Unsupported platform: ",satnum
		stop
	endif
    outfile = trim(outfile)//trim(timeslot(1:8))
    outfile = trim(outfile)//trim("_")
    outfile = trim(outfile)//trim(timeslot(9:12))
    outfile = trim(outfile)//trim("_")

	if (cnum .lt. 1 .or. cnum .gt. 16) then
		status = HIMAWARI_FAILURE
		return
	endif
	outfile = trim(outfile)//trim(tmpstr)
	outfile = trim(outfile)//trim("_FLDK_R")
	if (cnum.eq.1.or.cnum.eq.2.or.cnum.eq.4) then
		outfile = trim(outfile)//trim("10_S")
	else if (cnum.eq.3) then
		outfile = trim(outfile)//trim("05_S")
	else
		outfile = trim(outfile)//trim("20_S")
	endif
	write(tstr,"(I2.2)")seg
	outfile = trim(outfile)//trim(tstr)
	if (single_seg .neqv. .true.) then
    	outfile = trim(outfile)//trim("10.DAT")
    else
    	outfile = trim(outfile)//trim("01.DAT")
    endif
	
	INQUIRE(FILE=trim(outfile), EXIST=file_exists)
	if (file_exists .neqv. .true.) then
	    if (satnum .ne. 100) then
	        write(*,*) "ERROR: Input data file not found: ", trim(outfile)
		    stop
		else
		    idx = index(trim(outfile),"HS_H08_")
		    outfile = outfile(1:idx-1)//"HS_H09_"//outfile(idx+7:len(outfile))
		    INQUIRE(FILE=trim(outfile), EXIST=file_exists)
		    if (file_exists .neqv. .true.) then
	            write(*,*) "ERROR: Input data file not found: ", trim(outfile)
		        stop
		    endif
		endif
	endif
	

	status = HIMAWARI_SUCCESS
	return

end function AHI_get_file_name_seg

integer function AHI_file_exists(filename,verbose) result(status)

	character(len=*), intent(in) :: filename
	logical, intent(in) :: verbose
	logical exists

	inquire(file=filename, exist=exists)
	if (exists .eqv. .true.) then
  		status = HIMAWARI_SUCCESS;
	else
  		status = HIMAWARI_FAILURE;
  	endif
  	return
end function AHI_file_exists



integer function AHI_Create_TrueColour(ahi_data,ahi_extent,verbose) result(status)

	use omp_lib

	type(himawari_t_data), intent(inout) :: ahi_data
	type(himawari_t_extent), intent(in)    :: ahi_extent
	logical, intent(in) :: verbose

	integer :: x,y,n_threads
	real(kind=ahi_sreal) :: fval,rmin,rmax
	integer(kind=ahi_sint),dimension(4) :: scales

	fval = 0.13
	rmin = 0.00
	rmax = 1.25
#ifdef _OPENMP
	if (verbose) then
		n_threads = omp_get_max_threads()
		write(*,*) "Creating true-colour image using ",n_threads,'threads'
	endif
!$omp parallel DO PRIVATE(x,y)
#else
	if (verbose)write(*,*) 'Creating true-colour image without threading'
#endif

	do x=1,ahi_extent%x_max-ahi_extent%x_min
		do y=1,ahi_extent%y_max-ahi_extent%y_min
			ahi_data%tmpdata(x,y,1) = ahi_data%indata(x,y,1)
			ahi_data%tmpdata(x,y,2) = (1.-fval)*ahi_data%indata(x,y,2) + fval*ahi_data%indata(x,y,4)
			ahi_data%tmpdata(x,y,3) = ahi_data%indata(x,y,3)

		end do
	end do
#ifdef _OPENMP
!$omp end parallel do
#endif

	where(ahi_data%tmpdata(:,:,:) .gt. rmax) ahi_data%tmpdata(:,:,:) = rmax
	where(ahi_data%tmpdata(:,:,:) .lt. rmin) ahi_data%tmpdata(:,:,:) = rmin

	ahi_data%tmpdata(:,:,1) = ((ahi_data%tmpdata(:,:,1)-rmin)*255)/(rmax-rmin)
	ahi_data%tmpdata(:,:,2) = ((ahi_data%tmpdata(:,:,2)-rmin)*255)/(rmax-rmin)
	ahi_data%tmpdata(:,:,3) = ((ahi_data%tmpdata(:,:,3)-rmin)*255)/(rmax-rmin)

	ahi_data%tmpdata = sqrt(ahi_data%tmpdata)
	where(ahi_data%tmpdata(:,:,:) .gt. 180) ahi_data%tmpdata=sqrt(ahi_data%tmpdata)

	ahi_data%indata(:,:,1)=ahi_data%tmpdata(:,:,1)
	ahi_data%indata(:,:,2)=ahi_data%tmpdata(:,:,2)
	ahi_data%indata(:,:,3)=ahi_data%tmpdata(:,:,3)

	status = HIMAWARI_SUCCESS


end function AHI_Create_TrueColour

end module himawari_utils
