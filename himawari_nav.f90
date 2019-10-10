!******************************************************************************%
! *
! *    Copyright (C) 2016-2019 Simon Proud <simon.proud@physics.ox.ac.uk>
! *    License: CC BY-NC-ND 4.0
! *
! ******************************************************************************/


module himawari_navigation
	use himawari
	use iso_c_binding
	use himawari_headerinfo

	implicit none

	public	::	AHI_Pix2Geo, &
			AHI_DefaultNav, &
			AHI_Calctime, &
			AHI_Solpos

	interface

		integer(c_int) function get_sza_saa(year,month,day,hour,minute,lat,lon,sza,saa) bind(C, name = 'get_sza_saa')

			use iso_c_binding

			implicit none

			integer(c_int),	intent(in), value	:: year
			integer(c_int),	intent(in), value	:: month
			integer(c_int),	intent(in), value	:: day
			integer(c_int),	intent(in), value	:: hour
			integer(c_int),	intent(in), value	:: minute
			real(c_float),	intent(in), value	:: lat
			real(c_float),	intent(in), value	:: lon
			real(c_float),	intent(out)		:: sza
			real(c_float),	intent(out)		:: saa

		end function get_sza_saa
	end interface
contains

integer function AHI_Solpos(year,month,day,hour,minute,lat,lon,sza,saa) result(status)

	implicit none

	integer,	intent(in), value	:: year
	integer,	intent(in), value	:: month
	integer,	intent(in), value	:: day
	integer,	intent(in), value	:: hour
	integer,	intent(in), value	:: minute
	real(kind=ahi_sreal),	intent(in), value	:: lat
	real(kind=ahi_sreal),	intent(in), value	:: lon
	real(kind=ahi_sreal),	intent(out)	:: sza
	real(kind=ahi_sreal),	intent(out)	:: saa

	integer	::	retval

	saa 		=	0.
	sza 		=	0.
	retval	=	0

    retval = get_sza_saa(year,month,day,hour,minute,lat,lon,sza,saa)
    if (retval .ne. 0) then
        write(6, *) 'ERROR: get_sza_saa()'
        return
    end if
	if (saa .gt. 360.0) then
			saa	=	him_sreal_fill_value
	endif
	if (saa .lt. 0.0) then
		saa	=	him_sreal_fill_value
	endif
	sza	=	abs(sza)
	if (sza .gt. 180.0) then
		sza	=	him_sreal_fill_value
	endif
	if (sza .lt. -180.0) then
		sza	=	him_sreal_fill_value
	endif

	status	=	HIMAWARI_SUCCESS
	return

end function AHI_Solpos

integer function AHI_DefaultNav(ahi_navi,him_nav,verbose) result(status)

	type(himawari_t_navdata), intent(inout)	::	ahi_navi
	type(himawari_t_Proj_Info), intent(in)	 	::	him_nav
	logical, intent(in)								:: verbose

	ahi_navi%subLon		=	him_nav%subLon
	ahi_navi%cfac			=	him_nav%cfac	/	HIMAWARI_DEGTORAD
	ahi_navi%lfac			=	him_nav%lfac	/	HIMAWARI_DEGTORAD
	ahi_navi%coff			=	him_nav%coff
	ahi_navi%loff			=	him_nav%loff
	ahi_navi%satDis		=	him_nav%satDis
	ahi_navi%eqtrRadius	=	him_nav%eqtrRadius
	ahi_navi%polrRadius	=	him_nav%polrRadius
	ahi_navi%projParam1	=	him_nav%projParam1
	ahi_navi%projParam2	=	him_nav%projParam2
	ahi_navi%projParam3	=	him_nav%projParam3
	ahi_navi%projParamSd	=	him_nav%projParamSd

	status	=	HIMAWARI_SUCCESS
	return

end function AHI_DefaultNav

integer function AHI_Pix2Geo(ahi_main,verbose) result(status)

	type(himawari_t_struct), intent(inout)	::	ahi_main
	logical, intent(in)							:: verbose

!	real		::	Sd,Sn,S1,S2,S3,Sxy
	integer	::	c1,l1
	integer	::	xsize
	integer	::	ysize

	real(kind=ahi_dreal),dimension(:,:),allocatable	::	c,l,x,y
	real(kind=ahi_dreal),dimension(:,:),allocatable	::	Sd,Sn,S1,S2,S3,Sxy

	xsize	=	HIMAWARI_IR_NLINES
	ysize	=	HIMAWARI_IR_NCOLS

	allocate(c(xsize,ysize))
	allocate(l(xsize,ysize))
	allocate(x(xsize,ysize))
	allocate(y(xsize,ysize))
	allocate(Sd(xsize,ysize))
	allocate(sn(xsize,ysize))
	allocate(s1(xsize,ysize))
	allocate(s2(xsize,ysize))
	allocate(s3(xsize,ysize))
	allocate(sxy(xsize,ysize))

	ahi_main%ahi_data%lat	=	him_sreal_fill_value
	ahi_main%ahi_data%lon	=	him_sreal_fill_value
	do c1=1,HIMAWARI_IR_NLINES
		do l1=1,HIMAWARI_IR_NCOLS
			c(c1,l1)=dble(c1)
			l(c1,l1)=dble(l1)
		enddo
	enddo

	x	=	( c - ahi_main%ahi_navdata%coff) / ( HIMAWARI_SCLUNIT * ahi_main%ahi_navdata%cfac)
	y	=	( l - ahi_main%ahi_navdata%loff) / ( HIMAWARI_SCLUNIT * ahi_main%ahi_navdata%lfac)

	sd	=	(ahi_main%ahi_navdata%satDis * cos(x) * cos(y)) * (ahi_main%ahi_navdata%satDis * cos(x) * cos(y)) -&
			 (cos(y) * cos(y) + ahi_main%ahi_navdata%projParam3 * sin(y) * sin(y)) *&
			  ahi_main%ahi_navdata%projParamSd
 	sd = sqrt(sd)

	sn = (ahi_main%ahi_navdata%satDis * cos(x) * cos(y) -sd) / (cos(y) * cos(y) +&
		 ahi_main%ahi_navdata%projParam3 * sin(y) * sin(y))
    	s1 = ahi_main%ahi_navdata%satDis - (sn * cos(x) * cos(y))
	s2 = sn * sin(x) * cos(y)
	s3 =-sn * sin(y)
	sxy=sqrt( s1 * s1 + s2 * s2)

	ahi_main%ahi_data%lon = sngl(HIMAWARI_RADTODEG * atan(s2/s1) + ahi_main%ahi_navdata%subLon)

	ahi_main%ahi_data%lat = sngl(atan(ahi_main%ahi_navdata%projParam3 * s3 / sxy) * HIMAWARI_RADTODEG)

	where (ahi_main%ahi_data%lon .gt. 180.0)
		ahi_main%ahi_data%lon	=	ahi_main%ahi_data%lon-360.0
	end where
	where (ahi_main%ahi_data%lon .lt. -180.0)
		ahi_main%ahi_data%lon	=	ahi_main%ahi_data%lon+360.0
	end where
!	ahi_main%ahi_data%lat	=	ahi_main%ahi_data%lat*-1


!	stop


	where (ahi_main%ahi_data%lon .gt. 180.0)
		ahi_main%ahi_data%lon	=	him_sreal_fill_value
	end where
	where (ahi_main%ahi_data%lon .lt. -180.0)
		ahi_main%ahi_data%lon	=	him_sreal_fill_value
	end where

	where (ahi_main%ahi_data%lat .gt. 90.0)
		ahi_main%ahi_data%lat	=	him_sreal_fill_value
	end where
	where (ahi_main%ahi_data%lat .lt. -90.0)
		ahi_main%ahi_data%lat	=	him_sreal_fill_value
	end where


	deallocate(x)
	deallocate(y)
	deallocate(c)
	deallocate(l)

	deallocate(sn)
	deallocate(sd)
	deallocate(s1)
	deallocate(s2)
	deallocate(s3)
	deallocate(sxy)

	status	=	HIMAWARI_SUCCESS
	return


end function AHI_Pix2Geo

integer function AHI_Calctime(ahi_main,verbose) result(status)

	use omp_lib

	type(himawari_t_struct), intent(inout)	::	ahi_main
	logical, intent(in)							:: verbose

	integer 			::	year,month,day,hour,minu,yearp,monthp,retval
	real(kind=ahi_dreal) 	::	a,b,c,d,jd,tfact
	integer			::	x,y
	real(kind=ahi_dreal)	::	start_jd,end_jd
	integer 			::	iye,mon,idy,ihr,min,ifail
	integer 			::	iyyy,jy,jm,igreg,ja,ijul
	integer			::	idint2,n_threads,tnr,i,t
	integer        :: xmin,ymin,xmax,ymax
	real(kind=ahi_dreal)	::	julian
	real(kind=ahi_sreal)	::	sza,saa
	real(kind=ahi_sreal)	::	sec
	real (kind=ahi_sreal)	::	doy

	real(kind=ahi_sreal):: bob1,bob2
	parameter (igreg=15+31*(10+12*1582))

	sec=0
	read(ahi_main%ahi_info%timeslot(1:4),'(i10)')iye
	read(ahi_main%ahi_info%timeslot(5:6),'(i10)')mon
	read(ahi_main%ahi_info%timeslot(7:8),'(i10)')idy
	read(ahi_main%ahi_info%timeslot(9:10),'(i10)')ihr
	read(ahi_main%ahi_info%timeslot(11:12),'(i10)')min

	if(iye.eq.0.or. iye.lt.-4713) then
		ifail=1
		return
	endif
	if(iye.lt.0) then
		iyyy=iye+1
	else
		iyyy=iye
	endif
	if(mon.gt.2) then
		jy=iyyy
		jm=mon+1
	else
		jy=iyyy-1
		jm=mon+13
	endif
	ijul=idint(365.25d0*dble(jy))+idint(30.6001d0*dble(jm))+idy+1720995
	if(idy+31*(mon+12*iyyy).ge.igreg) then
		ja=idint(0.01d0*dble(jy))
		ijul=ijul+2-ja+idint(0.25d0*dble(ja))
	endif
	julian=dble(ijul)+dble(ihr)/24.d0+dble(min)/1440.d0+dble(sec)/86400.d0-0.5d0

	tfact	=	10.0/(24.0*60.0)

	do y=ahi_main%ahi_extent%y_min,ahi_main%ahi_extent%y_max
			ahi_main%ahi_data%time(:,y - ahi_main%ahi_extent%y_min+1)=julian+tfact*(y/dble(HIMAWARI_IR_NLINES))
	enddo

	xmin = 1
	ymin = 1

	xmax = ahi_main%ahi_extent%x_max - ahi_main%ahi_extent%x_min
	ymax = ahi_main%ahi_extent%y_max - ahi_main%ahi_extent%y_min

#ifdef _OPENMP
!$omp parallel DO PRIVATE(i,x,y,sza,saa,tnr)
#endif
	do y=ymin,ymax
		do x=xmin,xmax
			retval	=	AHI_Solpos(iye,mon,idy,ihr,min,ahi_main%ahi_data%lat(x,y),ahi_main%ahi_data%lon(x,y),sza,saa)
			ahi_main%ahi_data%sza(x,y)=sza
			ahi_main%ahi_data%saa(x,y)=saa
		enddo
	enddo
#ifdef _OPENMP
!$omp end parallel do
#endif

	status	=	HIMAWARI_SUCCESS
	return


end function AHI_Calctime

integer function AHI_calc_satangs(ahi_main,verbose) result(status)

	type(himawari_t_struct), intent(inout)	::	ahi_main
	logical, intent(in)							:: verbose


	real(kind=ahi_dreal),dimension(:,:),allocatable	::	N
	real(kind=ahi_dreal),dimension(:,:),allocatable	::	x,y,z
	real(kind=ahi_dreal),dimension(:,:),allocatable	::	cos_lat,sin_lat,cos_lon,sin_lon
	real(kind=ahi_dreal),dimension(:,:),allocatable	::	qv1,qv2,qv3
	real(kind=ahi_dreal),dimension(:,:),allocatable	::	u1,u2,u3

	real(kind=ahi_dreal)		::	e2,a,b

	allocate(N(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))

	allocate(x(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))
	allocate(y(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))
	allocate(z(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))

	allocate(qv1(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))
	allocate(qv2(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))
	allocate(qv3(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))

	allocate(u1(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))
	allocate(u2(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))
	allocate(u3(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))

	allocate(cos_lat(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))
	allocate(sin_lat(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))
	allocate(cos_lon(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))
	allocate(sin_lon(HIMAWARI_IR_NLINES,HIMAWARI_IR_NCOLS))

	a	=	ahi_main%ahi_navdata%eqtrRadius
	b	=	ahi_main%ahi_navdata%polrRadius


   cos_lat = cos(ahi_main%ahi_data%lat * HIMAWARI_DEGTORAD)
   sin_lat = sin(ahi_main%ahi_data%lat * HIMAWARI_DEGTORAD)
   cos_lon = cos(ahi_main%ahi_data%lon * HIMAWARI_DEGTORAD)
   sin_lon = sin(ahi_main%ahi_data%lon * HIMAWARI_DEGTORAD)

   e2	=	1. - (b * b) / (a * a)

   N	=	a / sqrt(1. - e2 * sin_lat * sin_lat)

   x	=	N * cos_lat * cos_lon;
   y	=	N * cos_lat * sin_lon;
   z	=	((b * b) / (a * a) * N) * sin_lat;

   qv1	=	ahi_main%ahi_navdata%satdis * cos(ahi_main%ahi_navdata%subLon * HIMAWARI_DEGTORAD)
   qv2	=	ahi_main%ahi_navdata%satdis * sin(ahi_main%ahi_navdata%subLon * HIMAWARI_DEGTORAD)
   qv3	=	0

	qv1 = qv1 - x
	qv2 = qv2 - y
	qv3 = qv3 - z
	!s,e,z
   u1 = (-sin_lat * cos_lon * qv1) + (-sin_lat * sin_lon * qv2) + (cos_lat * qv3)
   u2 = (-sin_lon *           qv1) + ( cos_lon           * qv2)
   u3 = ( cos_lat * cos_lon * qv1) + (-cos_lat * sin_lon * qv2) + (sin_lat * qv3)

	ahi_main%ahi_data%vza = acos(u3 / sqrt(u1*u1 + u2*u2 + u3*u3)) * HIMAWARI_RADTODEG

	ahi_main%ahi_data%vaa = atan2(-u2, u1) * HIMAWARI_RADTODEG

	where (ahi_main%ahi_data%vaa .lt. 0)
		ahi_main%ahi_data%vaa	=	ahi_main%ahi_data%vaa+360.0
	end where

	where (ahi_main%ahi_data%vaa .gt. 360)
		ahi_main%ahi_data%vaa	=	him_sreal_fill_value
	end where
	where (ahi_main%ahi_data%vaa .lt. 0)
		ahi_main%ahi_data%vaa	=	him_sreal_fill_value
	end where
	where (ahi_main%ahi_data%vza .lt. 0)
		ahi_main%ahi_data%vza	=	him_sreal_fill_value
	end where
	where (ahi_main%ahi_data%vza .gt. 180)
		ahi_main%ahi_data%vza	=	him_sreal_fill_value
	end where

	deallocate(N)

	deallocate(x)
	deallocate(y)
	deallocate(z)

	deallocate(qv1)
	deallocate(qv2)
	deallocate(qv3)

	deallocate(u1)
	deallocate(u2)
	deallocate(u3)

	deallocate(cos_lat)
	deallocate(sin_lat)
	deallocate(cos_lon)
	deallocate(sin_lon)

	status	=	HIMAWARI_SUCCESS
     return

end function AHI_calc_satangs


end module himawari_navigation

