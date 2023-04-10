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
	use omp_lib

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
    if (retval /= 0) then
        write(6, *) 'ERROR: get_sza_saa()'
        return
    end if
	if (saa > 360.0) then
			saa	=	him_sreal_fill_value
	endif
	if (saa < 0.0) then
		saa	=	him_sreal_fill_value
	endif
	sza	=	abs(sza)
	if (sza > 180.0) then
		sza	=	him_sreal_fill_value
	endif
	if (sza < -180.0) then
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
	integer	::	c1, l1
	integer	::	xsize, ysize
	integer	::	cur_x, cur_y

	real(kind=ahi_dreal),dimension(:,:),allocatable	::	c,l,x,y
	real(kind=ahi_dreal),dimension(:,:),allocatable	::	Sd,Sn,S1,S2,S3,Sxy

	xsize	=	HIMAWARI_IR_NLINES
	ysize	=	HIMAWARI_IR_NCOLS

	allocate(c(ahi_main%ahi_extent%x_size, ahi_main%ahi_extent%y_size))
	allocate(l(ahi_main%ahi_extent%x_size, ahi_main%ahi_extent%y_size))
	allocate(x(ahi_main%ahi_extent%x_size, ahi_main%ahi_extent%y_size))
	allocate(y(ahi_main%ahi_extent%x_size, ahi_main%ahi_extent%y_size))
	allocate(Sd(ahi_main%ahi_extent%x_size, ahi_main%ahi_extent%y_size))
	allocate(sn(ahi_main%ahi_extent%x_size, ahi_main%ahi_extent%y_size))
	allocate(s1(ahi_main%ahi_extent%x_size, ahi_main%ahi_extent%y_size))
	allocate(s2(ahi_main%ahi_extent%x_size, ahi_main%ahi_extent%y_size))
	allocate(s3(ahi_main%ahi_extent%x_size, ahi_main%ahi_extent%y_size))
	allocate(sxy(ahi_main%ahi_extent%x_size, ahi_main%ahi_extent%y_size))

	ahi_main%ahi_data%lat	=	him_sreal_fill_value
	ahi_main%ahi_data%lon	=	him_sreal_fill_value
	do c1=ahi_main%ahi_extent%x_min, ahi_main%ahi_extent%x_max
		cur_x = 1 + (c1 - ahi_main%ahi_extent%x_min)
		do l1=ahi_main%ahi_extent%y_min, ahi_main%ahi_extent%y_max
			cur_y = 1 + (l1 - ahi_main%ahi_extent%y_min)
			c(cur_x,cur_y) = dble(c1)
			l(cur_x,cur_y) = dble(l1)
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

	where (ahi_main%ahi_data%lon > 180.0)
		ahi_main%ahi_data%lon	=	ahi_main%ahi_data%lon-360.0
	end where
	where (ahi_main%ahi_data%lon < -180.0)
		ahi_main%ahi_data%lon	=	ahi_main%ahi_data%lon+360.0
	end where
!	ahi_main%ahi_data%lat	=	ahi_main%ahi_data%lat*-1


!	stop

	where (ahi_main%ahi_data%lon > 180.0)
		ahi_main%ahi_data%lon	=	him_sreal_fill_value
	end where
	where (ahi_main%ahi_data%lon < -180.0)
		ahi_main%ahi_data%lon	=	him_sreal_fill_value
	end where

	where (ahi_main%ahi_data%lat > 90.0)
		ahi_main%ahi_data%lat	=	him_sreal_fill_value
	end where
	where (ahi_main%ahi_data%lat < -90.0)
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
	real(kind=ahi_sreal)	::	doy
	
	real(kind=ahi_sreal), allocatable :: tmparr(:,:)

	real(kind=ahi_sreal):: bob1,bob2
	parameter (igreg=15+31*(10+12*1582))

	sec=0
	read(ahi_main%ahi_info%timeslot(1:4),'(i10)')iye
	read(ahi_main%ahi_info%timeslot(5:6),'(i10)')mon
	read(ahi_main%ahi_info%timeslot(7:8),'(i10)')idy
	read(ahi_main%ahi_info%timeslot(9:10),'(i10)')ihr
	read(ahi_main%ahi_info%timeslot(11:12),'(i10)')min

	if(iye==0.or. iye<-4713) then
		ifail=1
		return
	endif
	if(iye<0) then
		iyyy=iye+1
	else
		iyyy=iye
	endif
	if(mon>2) then
		jy=iyyy
		jm=mon+1
	else
		jy=iyyy-1
		jm=mon+13
	endif
	ijul=idint(365.25d0*dble(jy))+idint(30.6001d0*dble(jm))+idy+1720995
	if(idy+31*(mon+12*iyyy)>=igreg) then
		ja=idint(0.01d0*dble(jy))
		ijul=ijul+2-ja+idint(0.25d0*dble(ja))
	endif
	julian=dble(ijul)+dble(ihr)/24.d0+dble(min)/1440.d0+dble(sec)/86400.d0-0.5d0

	tfact	=	10.0/(24.0*60.0)

#ifdef _OPENMP
!$omp parallel DO PRIVATE(y)
#endif
	do y=ahi_main%ahi_extent%y_min,ahi_main%ahi_extent%y_max
			ahi_main%ahi_data%time(:,y - ahi_main%ahi_extent%y_min+1)=julian+tfact*(y/dble(HIMAWARI_IR_NLINES))
	enddo
#ifdef _OPENMP
!$omp end parallel do
#endif

	xmin = 1
	ymin = 1

	xmax = ahi_main%ahi_extent%x_max - ahi_main%ahi_extent%x_min + 1 
	ymax = ahi_main%ahi_extent%y_max - ahi_main%ahi_extent%y_min + 1
    if(ahi_main%do_solar_angles .eqv. .true.) then
#ifdef _OPENMP
	if (verbose) then
		n_threads = omp_get_max_threads()
		write(*,*) 'Processing solar geometry using',n_threads,'threads'
	endif
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
    else
        ahi_main%ahi_data%sza(:,:) = him_sreal_fill_value
        ahi_main%ahi_data%saa(:,:) = him_sreal_fill_value
    endif
	status	=	HIMAWARI_SUCCESS
	return


end function AHI_Calctime

integer function AHI_calc_satangs(ahi_main,verbose) result(status)
    ! This function calculates viewing angles for AHI.
    ! Loosely based on the calculations described here:
    ! http://celestrak.com/columns/

    type(himawari_t_struct), intent(inout):: ahi_main
    logical, intent(in)	                  :: verbose

    real(kind=ahi_dreal) :: sin_o_lat, sin_o_lon, cos_o_lat, cos_o_lon
    real(kind=ahi_dreal) :: obs_x, obs_y, obs_z, obs_alt
    real(kind=ahi_dreal) :: sat_x, sat_y, sat_z, sat_alt
    real(kind=ahi_dreal) :: del_x, del_y, del_z
    real(kind=ahi_dreal) :: a, b, r, satlon, satlat, flatten, ut1
    real(kind=ahi_dreal) :: theta_o, theta_s
    real(kind=ahi_dreal) :: top_s, top_e, top_z, az_, rg_, el_
    
    real(kind=ahi_dreal) :: c, sq, achcp
	integer              :: xmin, ymin, xmax, ymax, xpos, ypos
	integer              :: n_threads
    
    ! Navigational data from the Himawari dataset
    a = ahi_main%ahi_navdata%eqtrRadius
    b = ahi_main%ahi_navdata%polrRadius
    r = ahi_main%ahi_navdata%satdis - ahi_main%ahi_navdata%eqtrRadius
    satlon = ahi_main%ahi_navdata%subLon
    satlat = 0 ! We set this to zero, in reality it varies by <1 degree.
    
    flatten = (a - b) / b
    
    ! Here we set the observer to be on the geoid surface
    ! This could be improved using a DEM or similar
    obs_alt = 0
    
	xmin = 1
	ymin = 1

	xmax = ahi_main%ahi_extent%x_max - ahi_main%ahi_extent%x_min + 1 
	ymax = ahi_main%ahi_extent%y_max - ahi_main%ahi_extent%y_min + 1
#ifdef _OPENMP
	if (verbose) then
		n_threads = omp_get_max_threads()
		write(*,*) 'Processing viewing geometry using',n_threads,'threads'
	endif
!$omp parallel DO PRIVATE(cos_o_lat, sin_o_lat, cos_o_lon, sin_o_lon, ut1, theta_o, c, sq, achcp, sat_x, sat_y, sat_z, obs_x, obs_y, obs_z, del_x, del_y, del_z, top_s, top_e, top_z, az_, rg_, el_)
#endif
    do ypos=ymin,ymax
	    do xpos=xmin,xmax

            cos_o_lat = cos(ahi_main%ahi_data%lat(xpos, ypos) * HIMAWARI_DEGTORAD)
            sin_o_lat = sin(ahi_main%ahi_data%lat(xpos, ypos) * HIMAWARI_DEGTORAD)
            cos_o_lon = cos(ahi_main%ahi_data%lon(xpos, ypos) * HIMAWARI_DEGTORAD)
            sin_o_lon = sin(ahi_main%ahi_data%lon(xpos, ypos) * HIMAWARI_DEGTORAD)
            
            ut1 = ahi_main%ahi_data%time(xpos, ypos) / 36525.0
            
            ! This calculates observer position
            theta_o = 67310.54841 + ut1 * (876600.0 * 3600.0 + 8640184.812866 + ut1 * (0.093104 - ut1 * 6.2 * 10e-6))
            theta_o = mod(theta_o / 240.0, 2 * HIMAWARI_PI)
            theta_o = mod(theta_o + ahi_main%ahi_data%lon(xpos, ypos) * HIMAWARI_DEGTORAD, 2 * HIMAWARI_PI)
            c = 1 / sqrt(1 + flatten * (flatten - 2) * sin_o_lat**2)
            sq = c * (1 - flatten)**2 
            achcp = (a * c + obs_alt) * cos_o_lat
            obs_x = achcp * cos(theta_o)
            obs_y = achcp * sin(theta_o)
            obs_z = (a * sq + obs_alt) * sin_o_lat
            
            ! This calculates satellite position
            theta_s = 67310.54841 + ut1 * (876600.0 * 3600.0 + 8640184.812866 + ut1 * (0.093104 - ut1 * 6.2 * 10e-6))
            theta_s = mod(theta_s / 240.0, 2 * HIMAWARI_PI)
            theta_s = mod(theta_s + satlon * HIMAWARI_DEGTORAD, 2 * HIMAWARI_PI)
            achcp = ahi_main%ahi_navdata%satdis * cos(satlat)
            sat_x = achcp * cos(theta_s)
            sat_y = achcp * sin(theta_s)
            sat_z = ahi_main%ahi_navdata%satdis * sin(satlat)
            
            del_x = sat_x - obs_x
            del_y = sat_y - obs_y
            del_z = sat_z - obs_z
            
            top_s =  sin_o_lat * cos(theta_o) * del_x +  sin_o_lat * sin(theta_o) * del_y - cos_o_lat * del_z
            top_e = -sin(theta_o) * del_x + cos(theta_o) * del_y
            top_z = cos_o_lat * cos(theta_o) * del_x + cos_o_lat * sin(theta_o) * del_y +  sin_o_lat * del_z

            az_ = abs(atan2(-top_e, top_s) + HIMAWARI_PI)
            rg_ = sqrt(del_x * del_x + del_y * del_y + del_z * del_z)
            el_ = asin(top_z / rg_)
            
            ahi_main%ahi_data%vza(xpos, ypos) = 90. - el_ * HIMAWARI_RADTODEG
            ahi_main%ahi_data%vaa(xpos, ypos) = az_ * HIMAWARI_RADTODEG
            
            
        enddo
    enddo

	status	=	HIMAWARI_SUCCESS
     return

end function AHI_calc_satangs


end module himawari_navigation

