!******************************************************************************%
! *
! *    Copyright (C) 2016-2019 Simon Proud <simon.proud@physics.ox.ac.uk>
! *    License: CC BY-NC-ND 4.0
! * 
! 2020/05/18 CP made change to open call as later version of gfortran does not like 'newunit'  als small change to inquire call for same reason
!
! ******************************************************************************/


module himawari_readwrite
	use himawari
	use himawari_utils
	use himawari_navigation
	use himawari_headerinfo
	use himawari_readheader
	use netcdf
	use omp_lib

	implicit none

	public	::	AHI_Main_Read, &
				AHI_Setup_Read_Chans, &
				AHI_readchan, &
				AHI_resample_hres, &
				AHI_SavetoNCDF, &
				AHI_NCDF_check

contains


integer function AHI_Main_Read(filename, geofile, ahi_data2, &
                               ahi_extent, n_bands, band_ids,&
                               do_not_alloc, do_geo, predef_geo, &
                               do_solar, vis_res, satposstr, &
                               do_solar_angles, verbose) result(status)

	character(len=*), intent(in)				::	filename
	character(len=*), intent(in)				::	geofile
	type(himawari_t_data), intent(inout)	::	ahi_data2
	type(himawari_t_extent), intent(inout)	::	ahi_extent
	integer,intent(in)							::	n_bands
	integer,intent(in),dimension(16)			::	band_ids
	integer,intent(in)							::	do_not_alloc
	integer,intent(in)							::	do_geo
	logical,intent(in)							:: predef_geo
	logical,intent(in)							:: do_solar
	logical,intent(in)							:: vis_res
	character(len=128), intent(inout)		::	satposstr
	logical,intent(in)							:: do_solar_angles
	logical,intent(in)							:: verbose

	integer											::	satnum
	character(len=HIMAWARI_CHARLEN)			::	timeslot
	character(len=HIMAWARI_CHARLEN)			::	indir
	integer,dimension(HIMAWARI_NCHANS)		::	inbands
	type(himawari_t_struct)						::	ahi_main
	character(len=HIMAWARI_CHARLEN)			::	satlat,satlon,sathei,eqrrad,polrad

	integer		::	i,retval,pos


	satnum	=	101
	pos		=	index(trim(filename),"HS_H08_")
	if (pos .le. 0) then
		pos		=	index(trim(filename),"HS_H09_")
		if (pos .le. 0) then
			status	=	HIMAWARI_FAILURE
			return
		endif
		satnum	=	102
	endif

	ahi_main%ahi_data%memory_alloc_d	=	do_not_alloc
	ahi_main%ahi_data%n_bands		=	n_bands
	ahi_main%ahi_extent = ahi_extent

	retval	=	AHI_get_timeslot(filename,timeslot)
	retval	=	AHI_get_indir(filename,indir)

	ahi_main%ahi_info%indir		=	indir
	ahi_main%ahi_info%timeslot	=	timeslot
	ahi_main%ahi_info%satnum	=	satnum

	ahi_main%vis_res	=	vis_res
	ahi_main%do_solar	=	do_solar
	ahi_main%do_solar_angles = do_solar_angles

	if (verbose) write(*,*)"Reading AHI data for ",trim(ahi_main%ahi_info%timeslot(1:12))

	ahi_main%inchans(:)	=	0

	do i=1,n_bands
		if(band_ids(i)>0 .and. band_ids(i).le.HIMAWARI_NCHANS)ahi_main%inchans(band_ids(i))	=	1
	enddo

	ahi_main%convert=HIMAWARI_UNIT_RBT

	if (verbose) then
		write(*,*)"	-	Will process bands: ",ahi_main%inchans
		write(*,*)"	-	Will read data for region: "
		write(*,*)"			X: ",ahi_main%ahi_extent%x_min," to ", ahi_main%ahi_extent%x_max
		write(*,*)"			Y: ",ahi_main%ahi_extent%y_min," to ", ahi_main%ahi_extent%y_max
		write(*,*)"	-	Region size is: ",ahi_main%ahi_extent%x_size,"x",ahi_main%ahi_extent%y_size

	endif

	if (ahi_main%ahi_data%memory_alloc_d.ne.1) then
		retval	=	AHI_alloc_vals(ahi_main,ahi_extent,verbose)
		if (retval.ne.HIMAWARI_SUCCESS) then
			status	=	HIMAWARI_FAILURE
			return
		endif
	else
		ahi_main%ahi_data = ahi_data2
	endif

	retval	=	AHI_Setup_Read_Chans(ahi_main,verbose)
	if (retval.ne.HIMAWARI_SUCCESS) then
		status	=	HIMAWARI_FAILURE
		return
	endif
	if (do_geo.eq.1) then
		if (.not. predef_geo) then
			if(verbose)write(*,*)"Computing lat/lon and satellite angles"
			retval	=	AHI_Pix2Geo(ahi_main,verbose)
			if (retval.ne.HIMAWARI_SUCCESS) then
				status	=	HIMAWARI_FAILURE
				return
			endif
			retval	=	AHI_calc_satangs(ahi_main,verbose)
			if (retval.ne.HIMAWARI_SUCCESS) then
				status	=	HIMAWARI_FAILURE
				return
			endif
		else
			if(verbose)write(*,*)"Retrieving lat/lon and satellite angles from file"
			retval	=	AHI_Retrieve_Predef_Geo(ahi_main,geofile,verbose)
			if (retval.ne.HIMAWARI_SUCCESS) then
				status	=	HIMAWARI_FAILURE
				return
			endif
		endif
		retval	=	AHI_Calctime(ahi_main,verbose)
		if (retval.ne.HIMAWARI_SUCCESS) then
			status	=	HIMAWARI_FAILURE
			return
		endif
	endif

	! Define the satellite position string, used for parallax correction

   write(satlat,'(f10.7)')0.0
   write(satlon,'(f10.5)')ahi_main%ahi_navdata%sublon
   write(sathei,'(f10.2)')ahi_main%ahi_navdata%satDis
   write(eqrrad,'(f10.2)')ahi_main%ahi_navdata%eqtrRadius
   write(polrad,'(f10.2)')ahi_main%ahi_navdata%polrRadius
   satposstr=trim(satlat)//","//trim(satlon)//","//trim(sathei)//','//trim(eqrrad)//","//trim(polrad)

	ahi_data2	=	ahi_main%ahi_data


	if (ahi_main%ahi_data%memory_alloc_d.ne.1) then
		retval	=	AHI_free_vals(ahi_main)
		if (retval.ne.HIMAWARI_SUCCESS) then
			status	=	HIMAWARI_FAILURE
			return
		endif
	endif

	status	=	HIMAWARI_SUCCESS
	return

end function AHI_Main_Read

integer function AHI_Retrieve_Predef_Geo(ahi_main,geofile,verbose) result(status)
	use netcdf
	implicit none

	type(himawari_t_struct), intent(inout)		::	ahi_main
	character(len=*), intent(in)			      ::	geofile
	logical,intent(in)								:: verbose

	integer						::	ncid, varid
	integer,dimension(2)		::	start,countval

	start(1)	=	ahi_main%ahi_extent%x_min
	start(2)	=	ahi_main%ahi_extent%y_min

	countval(1)	=	ahi_main%ahi_extent%x_max - ahi_main%ahi_extent%x_min + 1
	countval(2)	=	ahi_main%ahi_extent%y_max - ahi_main%ahi_extent%y_min + 1

	call AHI_NCDF_check( nf90_open(geofile, NF90_NOWRITE, ncid) )
	call AHI_NCDF_check( nf90_inq_varid(ncid, "Lat", varid) )
	call AHI_NCDF_check( nf90_get_var(ncid, varid, ahi_main%ahi_data%lat, start = start, count = countval) )
	call AHI_NCDF_check( nf90_inq_varid(ncid, "Lon", varid) )
	call AHI_NCDF_check( nf90_get_var(ncid, varid, ahi_main%ahi_data%lon, start = start, count = countval) )
	call AHI_NCDF_check( nf90_inq_varid(ncid, "VZA", varid) )
	call AHI_NCDF_check( nf90_get_var(ncid, varid, ahi_main%ahi_data%vza, start = start, count = countval) )
	call AHI_NCDF_check( nf90_inq_varid(ncid, "VAA", varid) )
	call AHI_NCDF_check( nf90_get_var(ncid, varid, ahi_main%ahi_data%vaa, start = start, count = countval) )

	status	=	HIMAWARI_SUCCESS

end function AHI_Retrieve_Predef_Geo

integer function AHI_Get_SegStart_Point(ahi_main,curseg,verbose) result(status)

	type(himawari_t_struct), intent(inout)	::	ahi_main
	integer											::	curseg
	logical,intent(in)							:: verbose

	integer 	:: startline

	ahi_main%ahi_extent%startpos(1)	=	ahi_main%ahi_extent%y_min - ahi_main%ahi_extent%segpos_ir(curseg)
	ahi_main%ahi_extent%startpos(2)	=	ahi_main%ahi_extent%y_min*2 - ahi_main%ahi_extent%segpos_vi(curseg)
	ahi_main%ahi_extent%startpos(3)	=	ahi_main%ahi_extent%y_min*4 - ahi_main%ahi_extent%segpos_hv(curseg)

	if (ahi_main%ahi_extent%startpos(1) .le. 0.) then
		ahi_main%ahi_extent%startpos(1) = 1
		ahi_main%ahi_extent%startpos(2) = 1
		ahi_main%ahi_extent%startpos(3) = 1
	endif

	if (ahi_main%ahi_extent%startpos(1) .gt. ahi_main%ahi_extent%segdel_ir) then
		status = HIMAWARI_FAILURE
		return
	endif
	if (ahi_main%ahi_extent%startpos(2) .gt. ahi_main%ahi_extent%segdel_vi) then
		status = HIMAWARI_FAILURE
		return
	endif
	if (ahi_main%ahi_extent%startpos(3) .gt. ahi_main%ahi_extent%segdel_hv) then
		status = HIMAWARI_FAILURE
		return
	endif

	status = HIMAWARI_SUCCESS

	return

end function AHI_Get_SegStart_Point

integer function AHI_Get_SegEnd_Point(ahi_main,curseg,verbose) result(status)

	type(himawari_t_struct), intent(inout)	::	ahi_main
	integer											::	curseg
	logical,intent(in)							:: verbose

	integer 	:: startline

	ahi_main%ahi_extent%endpos(1)	=	ahi_main%ahi_extent%y_max - ahi_main%ahi_extent%segpos_ir(curseg)
	ahi_main%ahi_extent%endpos(2)	=	ahi_main%ahi_extent%y_max*2 - ahi_main%ahi_extent%segpos_vi(curseg)
	ahi_main%ahi_extent%endpos(3)	=	ahi_main%ahi_extent%y_max*4 - ahi_main%ahi_extent%segpos_hv(curseg)

	if (ahi_main%ahi_extent%endpos(1) .gt. ahi_main%ahi_extent%segdel_ir) ahi_main%ahi_extent%endpos(1) = ahi_main%ahi_extent%segdel_ir
	if (ahi_main%ahi_extent%endpos(2) .gt. ahi_main%ahi_extent%segdel_vi) ahi_main%ahi_extent%endpos(2) = ahi_main%ahi_extent%segdel_vi
	if (ahi_main%ahi_extent%endpos(3) .gt. ahi_main%ahi_extent%segdel_hv) ahi_main%ahi_extent%endpos(3) = ahi_main%ahi_extent%segdel_hv

	if (ahi_main%ahi_extent%endpos(1) .le. 0.) then
		status = HIMAWARI_FAILURE
		return
	endif
	if (ahi_main%ahi_extent%endpos(2) .le. 0.) then
		status = HIMAWARI_FAILURE
		return
	endif
	if (ahi_main%ahi_extent%endpos(3) .le. 0.) then
		status = HIMAWARI_FAILURE
		return
	endif

	status = HIMAWARI_SUCCESS

	return

end function AHI_Get_SegEnd_Point

integer function AHI_Setup_Segments(ahi_main,verbose) result(status)

	type(himawari_t_struct), intent(inout)	::	ahi_main
	logical,intent(in)							:: verbose

	integer,dimension(10)						::	seg_del_start
	integer,dimension(10)						::	seg_del_end

	logical,dimension(10)						::	proc_st
	logical,dimension(10)						::	proc

	seg_del_start 	=	ahi_main%ahi_extent%segpos_ir - ahi_main%ahi_extent%y_min
	seg_del_end 	=	ahi_main%ahi_extent%segpos_ir - ahi_main%ahi_extent%y_max

	proc_st(:)	=	.false.
	proc(:)	=	.false.

	ahi_main%ahi_extent%procseg(:)=.false.
	where (seg_del_start .gt. (-1. * ahi_main%ahi_extent%segdel_ir)) proc_st = .true.
	where (seg_del_end .lt. 0. .and. proc_st .eqv. .true.) ahi_main%ahi_extent%procseg=.true.

	status	=	HIMAWARI_SUCCESS

	return

end function AHI_Setup_Segments

integer function AHI_Setup_Read_Chans(ahi_main,verbose) result(status)

	type(himawari_t_struct), intent(inout)		::	ahi_main
	logical,intent(in)								:: verbose

	character(HIMAWARI_CHARLEN)				::	fname
	real(kind=ahi_sreal), DIMENSION(:,:), ALLOCATABLE::	tdata2
	real(kind=ahi_sreal), DIMENSION(:,:), ALLOCATABLE::	tseg

	integer,dimension(10)					::	segpos

	integer	::	i,j,minseg,maxseg
	integer	::	retval,bandpos
	integer	::	segdel
	integer	::	startl,endl
	integer	:: xsize,ysize
	integer	::	xmin,xmax
	integer	::	cur_y,cur_y_s,oc
	integer	::	indvar
	integer	::	y_start,y_end

	if (verbose) then
		write(*,*)"Reading image data"
	endif

	retval	=	AHI_Setup_Segments(ahi_main,verbose)
	if (retval .ne. HIMAWARI_SUCCESS) then
		write(*,*)"Error in AHI_Setup_Segments"
		stop
	endif

	bandpos	=	1
	minseg	=	1
	maxseg	=	10

	do i=1,HIMAWARI_NCHANS

		! Index of y position between segments
		cur_y		=	1
		! Index of how much data we're reading from a given segment
		cur_y_s	=	ahi_main%ahi_extent%segdel_ir

		if (ahi_main%inchans(i)==1) then
			if (verbose) then
				write(*,*)"Reading data for channel",i
			endif
			if (i.eq.1.or.i.eq.2.or.i.eq.4) then
				allocate(tseg(HIMAWARI_VIS_NCOLS,ahi_main%ahi_extent%segdel_vi))
				tseg(:,:) = him_sreal_fill_value
				if (ahi_main%vis_res .neqv. .true.) then
					allocate(tdata2(ahi_main%ahi_extent%x_size*2,ahi_main%ahi_extent%y_size*2))
					xsize	=	ahi_main%ahi_extent%x_size*2
					ysize	=	ahi_main%ahi_extent%y_size*2
					xmin	=	ahi_main%ahi_extent%x_min*2 - 1
					xmax	=	ahi_main%ahi_extent%x_max*2
				else
					allocate(tdata2(ahi_main%ahi_extent%x_size,ahi_main%ahi_extent%y_size))
					xsize	=	ahi_main%ahi_extent%x_size
					ysize	=	ahi_main%ahi_extent%y_size
					xmin	=	ahi_main%ahi_extent%x_min
					xmax	=	ahi_main%ahi_extent%x_max
				endif

				tdata2(:,:) = him_sreal_fill_value
				segdel	=	ahi_main%ahi_extent%segdel_vi
				segpos	=	ahi_main%ahi_extent%segpos_vi

				indvar	=	2
			else if (i.eq.3) then
				allocate(tseg(HIMAWARI_HVI_NCOLS,ahi_main%ahi_extent%segdel_hv))
				if (ahi_main%vis_res .neqv. .true.) then
					allocate(tdata2(ahi_main%ahi_extent%x_size*4,ahi_main%ahi_extent%y_size*4))
					xsize	=	ahi_main%ahi_extent%x_size*4
					ysize	=	ahi_main%ahi_extent%y_size*4
					xmin	=	ahi_main%ahi_extent%x_min*4 - 3
					xmax	=	ahi_main%ahi_extent%x_max*4
				else
					allocate(tdata2(ahi_main%ahi_extent%x_size*2,ahi_main%ahi_extent%y_size*2))
					xsize	=	ahi_main%ahi_extent%x_size*2
					ysize	=	ahi_main%ahi_extent%y_size*2
					xmin	=	ahi_main%ahi_extent%x_min*2 - 1
					xmax	=	ahi_main%ahi_extent%x_max*2
				endif
				tdata2(:,:) = him_sreal_fill_value
				tseg(:,:) = him_sreal_fill_value
				segdel	=	ahi_main%ahi_extent%segdel_hv
				segpos	=	ahi_main%ahi_extent%segpos_hv
				indvar	=	3
			else
				allocate(tseg(HIMAWARI_IR_NCOLS,ahi_main%ahi_extent%segdel_ir))
				tseg(:,:) = him_sreal_fill_value
				if (ahi_main%vis_res .neqv. .true.) then
					allocate(tdata2(ahi_main%ahi_extent%x_size,ahi_main%ahi_extent%y_size))
					xsize	=	ahi_main%ahi_extent%x_size
					ysize	=	ahi_main%ahi_extent%y_size
				else
					write(*,*)"Cannot process at visible resolution as you have selected a thermal channel."
					stop
				endif
				xmin		=	ahi_main%ahi_extent%x_min
				xmax		=	ahi_main%ahi_extent%x_max
				tdata2(:,:) = him_sreal_fill_value
				segdel	=	ahi_main%ahi_extent%segdel_ir
				segpos	=	ahi_main%ahi_extent%segpos_ir
				indvar	=	1
			endif

			do j=minseg,maxseg

				if (ahi_main%ahi_extent%procseg(j) .eqv. .true.) then
					if (verbose) then
						write(*,*)"	-	Reading segment number",j
					endif

					retval	=	AHI_Get_SegStart_Point(ahi_main,j,verbose)
					if (retval .ne. HIMAWARI_SUCCESS) then
						write(*,*) "Cannot get segment starting points."
						status	=	HIMAWARI_FAILURE
						return
					endif
					retval	=	AHI_Get_SegEnd_Point(ahi_main,j,verbose)
					if (retval .ne. HIMAWARI_SUCCESS) then
						write(*,*) "Cannot get segment ending points."
						status	=	HIMAWARI_FAILURE
						return
					endif

					if (j.eq.0) then
						retval	=	AHI_get_file_name(i,ahi_main%ahi_info%timeslot,ahi_main%ahi_info%satnum,ahi_main%ahi_info%indir,fname,verbose)
					else
						retval	=	AHI_get_file_name_seg(i,j,ahi_main%ahi_info%timeslot,ahi_main%ahi_info%satnum,ahi_main%ahi_info%indir,fname,verbose)
					endif
					if (retval.ne.HIMAWARI_SUCCESS) then
						write(*,*)"Cannot get filename for band: ",i
						status	=	HIMAWARI_FAILURE
						return
					endif
					retval	=	ahi_file_exists(fname,verbose)
					if (retval.ne.HIMAWARI_SUCCESS) then
						write(*,*)"File does not exist: ",fname
						status	=	HIMAWARI_FAILURE
						return
					endif
					startl	=	segpos(j)
					endl		=	segpos(j)+segdel-1
					retval	=	AHI_readchan(fname,tseg,i,ahi_main%convert(i),ahi_main%ahi_navdata,verbose)

					y_start	=	cur_y
					y_end		=	cur_y + ahi_main%ahi_extent%endpos(indvar) - ahi_main%ahi_extent%startpos(indvar)

					tdata2(:,y_start:y_end) = tseg(xmin:xmax,ahi_main%ahi_extent%startpos(indvar):ahi_main%ahi_extent%endpos(indvar))

					cur_y		=	y_end + 1

					if (retval.ne.HIMAWARI_SUCCESS) then
						write(*,*)"Failed to read the channel: ",i
						status	=	HIMAWARI_FAILURE
						return
					endif
				endif
			enddo

			retval	=	AHI_resample_hres(tdata2, ahi_main%ahi_data%indata(:,:,bandpos),ahi_main%ahi_extent,xsize,ysize,verbose)
			if (retval.ne.HIMAWARI_SUCCESS) then
				write(*,*)"Cannot resample data for channel: ",i
				deallocate(tdata2)
				deallocate(tseg)
				status	=	HIMAWARI_FAILURE
				return
			endif
			deallocate(tdata2)
			deallocate(tseg)

			bandpos	=	bandpos + 1
		endif
	enddo

	status	=	HIMAWARI_SUCCESS

	return

end function AHI_Setup_Read_Chans

integer function AHI_readchan(fname, indata,band,convert,ahi_nav,verbose)result(status)

	character(len=*), intent(in)									::	fname
	real(kind=ahi_sreal), DIMENSION(:,:), intent(inout)	::	indata
	integer, intent(in)												::	band
	integer(kind=ahi_sint), intent(in)							::	convert
	type(himawari_t_navdata), intent(inout)					::	ahi_nav
	logical,intent(in)												:: verbose

	integer(2), DIMENSION(:,:), ALLOCATABLE				::	tdata
	real(kind=ahi_sreal), DIMENSION(:,:), ALLOCATABLE	::	tdata2
	real(kind=ahi_sreal), DIMENSION(:,:), ALLOCATABLE	::	temp3

	type(himawari_t_VIS_Header)				::	ahi_hdrvis
	type(himawari_t_IR_Header)					::	ahi_hdrir

	integer					::	reclen,xsize,ysize,i,x,y,bval,retval
	real 						::	temp,temp2
	integer					::	arrxs,arrys,filelun,flen
	real(kind=ahi_dreal)	::	gain,offset,c0,c1,c2,lspd,plnk,bolz,clamb
        logical:: ef
 
	bval	=	band

	if (band.eq.1.or.band.eq.2.or.band.eq.4) then
		open(unit=filelun, file=fname,form='unformatted',action='read',status='old',access='stream',convert='little_endian')
		arrxs	=	HIMAWARI_VIS_NLINES/10
		arrys	=	HIMAWARI_VIS_NCOLS
	else if (band.eq.3) then
		open(unit=filelun, file=fname,form='unformatted',action='read',status='old',access='stream',convert='little_endian')
		arrxs	=	HIMAWARI_HVI_NLINES/10
		arrys	=	HIMAWARI_HVI_NCOLS
	else
		open(unit=filelun, file=fname,form='unformatted',action='read',status='old',access='stream',convert='little_endian')
		arrxs	=	HIMAWARI_IR_NLINES/10
		arrys	=	HIMAWARI_IR_NCOLS
	endif

	allocate(tdata(arrys,arrxs))
	allocate(tdata2(arrys,arrxs))
	allocate(temp3(arrys,arrxs))

	if (band.lt.7) then
		retval					=	AHI_readhdr_VIS(filelun,ahi_hdrvis,verbose)
		gain						=	ahi_hdrvis%him_calib%gain_cnt2rad
		offset					=	ahi_hdrvis%him_calib%cnst_cnt2rad
		clamb						=	ahi_hdrvis%him_calib%waveLen
		c0							=	ahi_hdrvis%him_chan_calib%rad2albedo
		ahi_nav%subLon			=	ahi_hdrvis%him_proj%subLon
		ahi_nav%cfac			=	ahi_hdrvis%him_proj%cfac
		ahi_nav%lfac			=	ahi_hdrvis%him_proj%lfac
		ahi_nav%coff			=	ahi_hdrvis%him_proj%coff
		ahi_nav%loff			=	ahi_hdrvis%him_proj%loff

		ahi_nav%satdis			=	ahi_hdrvis%him_proj%satdis
		ahi_nav%eqtrRadius	=	ahi_hdrvis%him_proj%eqtrRadius
		ahi_nav%polrRadius	=	ahi_hdrvis%him_proj%polrRadius
		ahi_nav%projParam1	=	ahi_hdrvis%him_proj%projParam1
		ahi_nav%projParam2	=	ahi_hdrvis%him_proj%projParam2
		ahi_nav%projParam3	=	ahi_hdrvis%him_proj%projParam3
		ahi_nav%projParamSd	=	ahi_hdrvis%him_proj%projParamSd

		ahi_nav%cfac			=	ceiling((ahi_nav%cfac)/2)
		ahi_nav%lfac			=	ceiling((ahi_nav%lfac)/2)
		ahi_nav%coff			=	(ahi_nav%coff-0.5)/2
		ahi_nav%loff			=	(ahi_nav%loff-0.5)/2
		if (band.eq.3) then
			ahi_nav%cfac		=	ceiling(ahi_nav%cfac/2)
			ahi_nav%lfac		=	ceiling(ahi_nav%lfac/2)
			ahi_nav%coff		=	ahi_nav%coff/2
			ahi_nav%loff		=	ahi_nav%loff/2
		endif

		ahi_nav%coff			=	ahi_nav%coff+0.5
		ahi_nav%loff			=	ahi_nav%loff+0.5

	else
		retval					=	AHI_readhdr_IR(filelun,ahi_hdrir,verbose)
		gain						=	ahi_hdrir%him_calib%gain_cnt2rad
		offset					=	ahi_hdrir%him_calib%cnst_cnt2rad
		clamb						=	ahi_hdrir%him_calib%waveLen
		c0							=	ahi_hdrir%him_chan_calib%btp2rad_c0
		c1							=	ahi_hdrir%him_chan_calib%btp2rad_c1
		c2							=	ahi_hdrir%him_chan_calib%btp2rad_c2
		lspd						=	ahi_hdrir%him_chan_calib%lightSpeed
		plnk						=	ahi_hdrir%him_chan_calib%planckConst
		bolz						=	ahi_hdrir%him_chan_calib%bolzConst

		ahi_nav%subLon			=	ahi_hdrir%him_proj%subLon
		ahi_nav%cfac			=	ahi_hdrir%him_proj%cfac
		ahi_nav%lfac			=	ahi_hdrir%him_proj%lfac
		ahi_nav%coff			=	ahi_hdrir%him_proj%coff
		ahi_nav%loff			=	ahi_hdrir%him_proj%loff
		ahi_nav%satdis			=	ahi_hdrir%him_proj%satdis
		ahi_nav%eqtrRadius	=	ahi_hdrir%him_proj%eqtrRadius
		ahi_nav%polrRadius	=	ahi_hdrir%him_proj%polrRadius
		ahi_nav%projParam1	=	ahi_hdrir%him_proj%projParam1
		ahi_nav%projParam2	=	ahi_hdrir%him_proj%projParam2
		ahi_nav%projParam3	=	ahi_hdrir%him_proj%projParam3
		ahi_nav%projParamSd	=	ahi_hdrir%him_proj%projParamSd
	endif

	ahi_nav%cfac				=	ahi_nav%cfac/HIMAWARI_DEGTORAD
	ahi_nav%lfac				=	ahi_nav%lfac/HIMAWARI_DEGTORAD

	INQUIRE(FILE=trim(fname), SIZE=flen, EXIST=ef)
	flen	=	flen-(arrxs*arrys*2)
	call	fseek(filelun,flen,0,retval)
	read(filelun)tdata(:,:)

	if (convert.eq.HIMAWARI_UNIT_RAD .or. convert.eq.HIMAWARI_UNIT_RBT) then
		tdata2	=	float(tdata)*gain+offset

		if (convert == HIMAWARI_UNIT_RBT) then
			if (band<7) then
				tdata2	=	tdata2*c0
				where(tdata2.lt. -10.0) tdata2=him_sreal_fill_value
				where(tdata2.gt. 10.0) tdata2=him_sreal_fill_value
			else
				temp		=	(lspd*plnk)/(bolz*clamb/1e6);
				temp2	=	(2.0*lspd*lspd*plnk)/((clamb/1e6)**5)
				temp3	=	log(temp2/(tdata2*1e6)+1)

				tdata2	=	temp / temp3;
				tdata2	=	c0 + c1*tdata2 + c2*tdata2*tdata2
				where(tdata2.le. 140.0) tdata2=him_sreal_fill_value
				where(tdata2.gt. 400.0) tdata2=him_sreal_fill_value
			endif
		endif
	else
		tdata2	=	float(tdata)
	endif
	where(tdata.le.0)tdata2=him_sreal_fill_value
	indata	=	tdata2

	close(filelun)
	deallocate(tdata)
	deallocate(tdata2)
	deallocate(temp3)

	status	=	HIMAWARI_SUCCESS
	return

end function AHI_readchan


integer function AHI_resample_hres(indata, outdata,ahi_extent,xsize,ysize,verbose) result(status)

	use omp_lib

	real(kind=ahi_sreal), DIMENSION(:,:), intent(in)	::	indata
	real(kind=ahi_sreal), DIMENSION(:,:), intent(out)	::	outdata

	type(himawari_t_extent), intent(in)	::	ahi_extent

	integer, intent(in)							::	xsize
	integer, intent(in)							::	ysize
	logical,intent(in)							:: verbose

	real,dimension(ahi_extent%x_size,ahi_extent%y_size) :: temparr

	integer	::	x,y
	integer	::	outx,outy
	integer	::	i,j
	integer	::	inposvar
	integer	::	outposvar
	integer	::	sizerx,sizery
	real 		::	val
	integer 	::	inpix
	integer	::	n_threads

	if (ahi_extent%x_size .eq. xsize .or. ahi_extent%y_size .eq. ysize) then
		outdata	=	indata
		status	=	HIMAWARI_SUCCESS
		return
	endif
	temparr(:,:)=0

	sizerx	=	xsize / ahi_extent%x_size
	sizery	=	ysize / ahi_extent%y_size

#ifdef __PGI 
	if (verbose)write(*,*) 'Resampling VIS grid to IR grid using PGI threading'
#else
#ifdef _OPENMP
	if (verbose) then
		n_threads	=	omp_get_max_threads()
		write(*,*) 'Resampling VIS grid to IR grid using',n_threads,'threads'
	endif
!$omp parallel DO PRIVATE(x,y,outx,outy,val,inpix)
#else
	if (verbose)write(*,*) 'Resampling VIS grid to IR grid without threading'
#endif
#endif
!$acc kernels
!$acc loop collapse(2) independent private(x,y,i,j,outx,outy,val,inpix)
	do x=1,xsize-sizerx
		do y=1,ysize-sizery
			outx	=	int(x/sizerx)+1
			outy	=	int(y/sizery)+1
			val	=	0
			inpix=	0
			!$acc loop collapse(2)
			do i=1,sizerx
				do j=1,sizery
					if (indata(x+i,y+j).gt.-100) then
						val	=	val + indata(x+i,y+j)
						inpix=	inpix + 1
					endif
				enddo
			enddo
			val	=	val/inpix
			if (outx .le. 0 .or. outx .ge. HIMAWARI_IR_NLINES .or. &
			    outy .le. 0 .or. outy .ge. HIMAWARI_IR_NCOLS) then
			    	continue
			endif

			if (temparr(outx,outy).le. 0) then
				temparr(outx,outy)=val
			endif
		enddo

	enddo

#ifdef _OPENMP
!$omp end parallel do
#endif
#ifdef __PGI
!$acc end kernels 
#endif
	outdata	=	temparr
	status	=	HIMAWARI_SUCCESS

	return

end function AHI_resample_hres

integer function AHI_SavetoNCDF(outdata,ahi_extent,fname,bname,newfile,verbose) result(status)

	character(len=*), intent(in)					::	fname
	character(len=*), intent(in)					::	bname
	real(kind=ahi_sreal), DIMENSION(:,:), intent(in)		::	outdata
	type(himawari_t_extent) , intent(in)		:: ahi_extent
	integer, intent(in)								::	newfile
	logical,intent(in)								:: verbose

	integer, parameter		::	NDIMS = 2
	integer						::	nx , ny
	integer						::	ncid, varid, dimids(NDIMS)
	integer						::	x_dimid, y_dimid
	integer						::	x, y

	nx	=	ahi_extent%x_size
	ny	=	ahi_extent%y_size

	if (newfile.eq.1) then
		call AHI_NCDF_check(nf90_create(fname, NF90_CLOBBER, ncid))
		call AHI_NCDF_check(nf90_def_dim(ncid, "x", nx, x_dimid))
		call AHI_NCDF_check(nf90_def_dim(ncid, "y", ny, y_dimid))
	else
		call AHI_NCDF_check(nf90_open(fname, NF90_WRITE, ncid))
		call AHI_NCDF_check(nf90_redef(ncid))
		call AHI_NCDF_check(nf90_inq_dimid(ncid, 'x', x_dimid))
		call AHI_NCDF_check(nf90_inq_dimid(ncid, 'y', y_dimid))
	endif

	dimids =  (/ x_dimid, y_dimid /)
	call AHI_NCDF_check(nf90_def_var(ncid,bname, NF90_FLOAT, dimids, varid))
	call AHI_NCDF_check(nf90_enddef(ncid))
	call AHI_NCDF_check(nf90_put_var(ncid, varid, outdata))
	call AHI_NCDF_check(nf90_close(ncid))
	status	= HIMAWARI_SUCCESS
	return

end function AHI_SavetoNCDF


integer function AHI_SavetoNCDF_int(outdata,ahi_extent,fname,bname,newfile,verbose) result(status)

	character(len=*), intent(in)					::	fname
	character(len=*), intent(in)					::	bname
	integer(kind=ahi_sint), DIMENSION(:,:), intent(in)		::	outdata
	type(himawari_t_extent) , intent(in)		:: ahi_extent
	integer, intent(in)								::	newfile
	logical,intent(in)								:: verbose

	integer, parameter		::	NDIMS = 2
	integer						::	nx , ny
	integer						::	ncid, varid, dimids(NDIMS)
	integer						::	x_dimid, y_dimid
	integer						::	x, y

	nx	=	ahi_extent%x_size
	ny	=	ahi_extent%y_size

	if (newfile.eq.1) then
		call AHI_NCDF_check(nf90_create(fname, NF90_CLOBBER, ncid))
		call AHI_NCDF_check(nf90_def_dim(ncid, "x", NX, x_dimid))
		call AHI_NCDF_check(nf90_def_dim(ncid, "y", NY, y_dimid))
	else
		call AHI_NCDF_check(nf90_open(fname, NF90_WRITE, ncid))
		call AHI_NCDF_check(nf90_redef(ncid))
	endif

	dimids =  (/ x_dimid, y_dimid /)
	call AHI_NCDF_check(nf90_def_var(ncid,bname, NF90_SHORT, dimids, varid))
	call AHI_NCDF_check(nf90_enddef(ncid))
	call AHI_NCDF_check(nf90_put_var(ncid, varid, outdata))
	call AHI_NCDF_check(nf90_close(ncid))
	status	= HIMAWARI_SUCCESS
	return

end function AHI_SavetoNCDF_int

subroutine AHI_NCDF_check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      write(*,*)trim(nf90_strerror(status))
      stop "Stopped"
    end if
end subroutine AHI_NCDF_check

end module himawari_readwrite
