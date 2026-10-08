MODULE sn_simple
   !!======================================================================
   !!                       ***  MODULE  sn_simple  ***
   !! Ocean dynamics : spectral nudging of temperature and salinity
   !!                 for nested regional ocean simulations
   !!======================================================================
   !! History : 1.0  ! 2023  (Renata Tatsch Eidt, Anna Katavouta)  Original code
   !!                !       CMCC Foundation & National Oceanography Centre
   !!           1.0  ! 2026  Public release
   !!======================================================================
#if defined key_sn_simple
   !!----------------------------------------------------------------------
   !!   key_sn_simple                                   spectral nudging
   !!----------------------------------------------------------------------
  
   !! * Modules used
   USE oce
   USE dom_oce         ! ocean space and time domain variables
   USE lib_mpp         ! distribued memory computing
   USE in_out_manager  ! I/O manager
   USE iom
   USE restart         ! only for lrst_oce
   USE phycst          ! Define parameters for the routines
!   USE dtatsd
   USE lbclnk
   USE daymod          ! calendar                         (day routine)
   USE eosbn2
   USE prtctl          ! Print control
   USE fldread
   USE timing
   USE wrk_nemo

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   PRIVATE
   SAVE

!  * Routine accessibility
   PUBLIC dta_sn_init, dta_sn, tra_sn_simple, gather_map, scatter_map, space_filt_2d

   INTEGER :: numsdt1, numsdt2, numGo, & ! logical unit for data
              nsal1,  &                  ! record useid
              Lorder, &                  ! to be used in space_filt_2d
              inum_dta, id, jpk_init     ! used in dta_sn_init

   CHARACTER(len=256) :: cn_gamma_file   ! path to 3D nudging coefficient file

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_sn   ! structure of input SST
                                                   ! (file informations, fields read)
   INTEGER             ::   zdim(4)                ! how many files read in namelist
   INTEGER , PARAMETER ::   jp_salin  = 1          ! index of sal at T-point
   INTEGER , PARAMETER ::   jp_tempe  = 2          ! index of temp at T-point
   INTEGER , PARAMETER ::   jp_depor  = 3          ! index of swpth at T-point
   INTEGER , PARAMETER ::   jp_maskor = 4          ! index of ocean mask

CONTAINS

!----------------------------------------------------------------------------------------
   SUBROUTINE dta_sn_init
!----------------------------------------------------------------------------------------
      INTEGER  :: ios, ierr0, ierr1, ierr2, & ! local integers
                      ierr3, ierr4, ierr5
      INTEGER  ::   ji, jj, jk, jm            ! dummy loop indices
      CHARACTER(len=100)        ::  cn_dir    ! Root directory for location of ssr files
      TYPE(FLD_N), DIMENSION(4) ::   slf_i    ! array of namelist informations on the fields to read
      TYPE(FLD_N) :: sr_sal, sr_tem, sr_dep, sr_msk
!----------------------------------------------------------------------------------------

      !!
      NAMELIST/namsn/   cn_dir, sr_sal, sr_tem, sr_dep, sr_msk, cn_gamma_file
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dta_sn_init')
      !
      !  Initialisation
      ierr0 = 0  ;  ierr1 = 0  ;  ierr2 = 0  ;  ierr3 = 0  ; ierr4 = 0  ;  ierr5 = 0
      !
      REWIND( numnam_ref )      ! Namelist namtsd in reference namelist :
      REWIND(numnam_cfg)

      READ  ( numnam_ref, namsn, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtsd in reference namelist', lwp )

      ! Namelist namtsd in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, namsn, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtsd in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namsn )

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'rean_ini : Temperature & Salinity data '
         WRITE(numout,*) '~~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namsn'
         WRITE(numout,*) '   Initialisation of reanalysis ocean T & S with T &S input data   ln_tsd_init   = '
      ENDIF

!     allocate the arrays (if necessary)

      ALLOCATE( sf_sn(4), STAT=ierr0 ) !
      IF( ierr0 > 0 ) THEN
         CALL ctl_stop( 'init_rean: unable to allocate sf_sn structure' ); RETURN
      ENDIF

      CALL iom_open ( trim(cn_dir) // trim(sr_dep%clname), inum_dta )
      id = iom_varid( inum_dta, sr_dep%clvar, zdim )
      jpk_init = zdim(3)
!      IF(lwp) WRITE(numout,*) 'Dimension of vertical coordinate in reanalysis: ', jpk_init
!      IF(lwp) WRITE(numout,*) 'Dimension of x coordinate in reanalysis: ', zdim(1)
!      IF(lwp) WRITE(numout,*) 'Dimension of y coordinate in reanalysis: ', zdim(2)
      CALL iom_close( inum_dta )   ! Close the input file
      !
      !                         ! fill sf_rean with sn_tem & sn_sal and control
      !                         print
      slf_i(jp_salin) = sr_sal
      slf_i(jp_tempe) = sr_tem
      slf_i(jp_depor) = sr_dep
      slf_i(jp_maskor) = sr_msk

      ALLOCATE( sf_sn(jp_salin)%fnow(jpi,jpj,jpk_init  ) , STAT=ierr0 )
      IF( slf_i(jp_salin)%ln_tint )      ALLOCATE( sf_sn(jp_salin)%fdta(jpi,jpj,jpk_init,2) , STAT=ierr1 )
      ALLOCATE( sf_sn(jp_tempe)%fnow(jpi,jpj,jpk_init  ) , STAT=ierr2 )
      IF( slf_i(jp_tempe)%ln_tint )      ALLOCATE( sf_sn(jp_tempe)%fdta(jpi,jpj,jpk_init,2) , STAT=ierr3 )
      ALLOCATE( sf_sn(jp_depor)%fnow(jpi,jpj,jpk_init  ) , STAT=ierr4 )
      ALLOCATE( sf_sn(jp_maskor)%fnow(jpi,jpj,jpk_init  ) , STAT=ierr5 )
      WRITE(numout,*) 'ln_tint: ', slf_i(jp_salin)%ln_tint
         !
      IF( ierr0 + ierr1 + ierr2 + ierr3 + ierr4 + ierr5 > 0 ) THEN
         CALL ctl_stop( 'nam_sn : unable to allocate T & S data arrays' ) ; RETURN
      ENDIF
      !
      CALL fld_fill( sf_sn, slf_i, cn_dir, 'dta_sn', 'Temperature & Salinity data', 'namsn')!, no_print )
         !

!      WRITE(numout,*) 'teste size jp_maskor(1)=', size(sf_sn(jp_maskor)%fnow,DIM =1)
!      WRITE(numout,*) 'teste size jp_maskor(2)=', size(sf_sn(jp_maskor)%fnow,DIM =2)
!      WRITE(numout,*) 'teste size jp_maskor(3)=', size(sf_sn(jp_maskor)%fnow,DIM =3)
      IF( nn_timing == 1 )  CALL timing_stop('dta_sn_init')


   END SUBROUTINE dta_sn_init

!----------------------------------------------------------------------------------------
   SUBROUTINE dta_sn( kt, ptsal, pttem)
!----------------------------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt     ! ocean time-step
      !
      INTEGER ::   ji, jj, jk, jl, jk_init   ! dummy loop indicies
      INTEGER ::   ik, il0, il1, ii0, ii1, ij0, ij1        ! local integers
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out) ::   ptsal   ! T & S data
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out) ::   pttem
      REAL(wp)::   zl, zi
!----------------------------------------------------------------------------------------

      IF( nn_timing == 1 )  CALL timing_start('dta_sn')
      !
      CALL fld_read( kt, 1, sf_sn )      !==   read T & S data at kt time step ==!
      !
!      WRITE(numout,*)'spectral_SIZE1', size(ptsal,DIM =1)
!      WRITE(numout,*)'spectral_SIZE2', size(ptsal,DIM =2)
!      WRITE(numout,*)'spectral_SIZE3', size(ptsal,DIM =3)

!      WRITE(numout,*)'SF_SALIN number', sf_sn(jp_salin)%fnow(:,:,50)

! test:
!         sntem=(sf_sn(jp_tempe)%fnow(:,:,1:120))

      IF( kt == nit000 .AND. lwp )THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dta_sn: interpolates T & S data onto current mesh'
      ENDIF
      !
      !
! teste
!         sntem=(sf_sn(jp_maskor)%fnow(:,:,1:120))

         DO jk = 1, jpk      ! determines the intepolated T-Snprofiles at each (i,j) points
            DO jj= 1, jpj
               DO ji= 1, jpi
                  zl = gdept_0(ji,jj,jk)
                  !WRITE(numout,*)'zl', zl
                  IF( zl < sf_sn(jp_depor)%fnow(ji,jj,1) ) THEN
! above the first level of data
                     ptsal(ji,jj,jk) = sf_sn(jp_salin)%fnow(ji,jj,1)
                     pttem(ji,jj,jk) = sf_sn(jp_tempe)%fnow(ji,jj,1)
                  ELSEIF( zl > sf_sn(jp_depor)%fnow(ji,jj,jpk_init) ) THEN
! below the last level of data
                     ptsal(ji,jj,jk) = sf_sn(jp_salin)%fnow(ji,jj,jpk_init)
                     pttem(ji,jj,jk) = sf_sn(jp_tempe)%fnow(ji,jj,jpk_init)
                  ELSE
! inbetween : vertical interpolation between jk_init & jk_init+1
                     DO jk_init = 1, jpk_init-1                     ! when gdept(jk_init) < zl < gdept(jk_init+1)
                        !WRITE(numout,*)'teste:jk_init=', jk_init
!                        IF( sf_sn(jp_maskor)%fnow(ji,jj,jk_init+1) == 1) THEN
 !if there is no data fill down
!                           sf_sn(jp_salin)%fnow(ji,jj,jk_init+1) = sf_sn(jp_salin)%fnow(ji,jj,jk_init)
!                           sf_sn(jp_tempe)%fnow(ji,jj,jk_init+1) = sf_sn(jp_tempe)%fnow(ji,jj,jk_init)
!                        ENDIF
                        IF( (zl-sf_sn(jp_depor)%fnow(ji,jj,jk_init)) * (zl-sf_sn(jp_depor)%fnow(ji,jj,jk_init+1)) <= 0._wp ) THEN
                           zi = ( zl - sf_sn(jp_depor)%fnow(ji,jj,jk_init) ) / &
                        & (sf_sn(jp_depor)%fnow(ji,jj,jk_init+1)-sf_sn(jp_depor)%fnow(ji,jj,jk_init))

                           ptsal(ji,jj,jk) = sf_sn(jp_salin)%fnow(ji,jj,jk_init) + &
                        & (sf_sn(jp_salin)%fnow(ji,jj,jk_init+1)-sf_sn(jp_salin)%fnow(ji,jj,jk_init)) * zi
                           pttem(ji,jj,jk) = sf_sn(jp_tempe)%fnow(ji,jj,jk_init) + &
                        & (sf_sn(jp_tempe)%fnow(ji,jj,jk_init+1)-sf_sn(jp_tempe)%fnow(ji,jj,jk_init)) * zi
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         !ptsal(:,:,:,jp_tem) = ptsal(:,:,:,jp_tem) *tmask(:,:,:)
         !ptsal(:,:,:) = sf_sn(jp_salin)%fnow *tmask(:,:,:)
! test:
! test:
!         sntem=(sf_sn(jp_tempe)%fnow(:,:,1:120))*tmask(:,:,:)
!          sntem=pttem(:,:,:) *tmask(:,:,:)

!         pttem(:,:,:)=(sf_sn(jp_tempe)%fnow(:,:,1:120))
!         sntem=(sf_sn(jp_tempe)%fnow(:,:,1:120))
! uncoment here:
         ptsal(:,:,:) = ptsal(:,:,:) *tmask(:,:,:)
         pttem(:,:,:) = pttem(:,:,:) *tmask(:,:,:)
    !  WRITE(numout,*)'SF_ptsal', ptsal(50,50,:)
    !  WRITE(numout,*)'SF_SALIN number', sf_sn(jp_salin)%fnow(50,50,:)
    !  WRITE(numout,*)'spectral_SIZE1', size(ptsal,DIM =1)
    !  WRITE(numout,*)'spectral_SIZE2', size(ptsal,DIM =2)
!      WRITE(numout,*)'spectral_SIZE3', size(ptsal,DIM =3)
      !   WRITE(numout,*)'spectral_SIZE4', size(ptsal,DIM =4)

      IF( nn_timing == 1 )  CALL timing_stop('dta_sn')
      !
   END SUBROUTINE dta_sn


!-----------------------------------------------------------------------------------------------
   SUBROUTINE tra_sn_simple(kt)
!-----------------------------------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt                                ! ocean time step
      INTEGER             :: kk                                ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: temdta2, saldta2, &  ! model output at specified timestep
                                          diff_tem, diff_sal, &! difference between child and parent
                                          temdta1, saldta1, &  ! external file temp and sal
                                          GammaO, &            ! Gamma map for spec nudg in 3d
                                          SN_T, SN_S           ! nudging correction
      REAL(wp), DIMENSION(jpi,jpj)           :: tem_par,sal_par   ! 2d field after filter
      REAL(wp), DIMENSION(jpiglo,jpjglo,jpk) :: tem_map,sal_map   ! nodes gathered in global domain
      REAL(wp), DIMENSION(jpiglo,jpjglo,jpk) :: mask_map          ! tmask gathered in global domain
      REAL(wp), DIMENSION(jpiglo,jpjglo,jpk) :: tem_filt,sal_filt ! filtered field global

      REAL(wp), DIMENSION(jpi,jpj,jpk) :: rean_tem, rean_sal

!-----------------------------------------------------------------------------------------------
    
         WRITE(numout,*)'Spectral'

!        Time step for the nudging:
         IF (nsec_day-nhre*3600==100) THEN  !every 1 hour

! GET REANALYSIS
            CALL dta_sn( kt, rean_sal, rean_tem)

!            WRITE(numout,*)'un_ptsal', rean_tem(50,50,10)
!            WRITE(numout,*)'spectral_SIZE1', size(rean_sal,DIM =1)
!            WRITE(numout,*)'spectral_SIZE2', size(rean_sal,DIM =2)
!            WRITE(numout,*)'spectral_SIZE3', size(rean_sal,DIM =3)
!            WRITE(numout,*)'spectral_nudging'

! OPEN GAMMA MAP
            IF( kt == nit000 ) THEN
               CALL iom_open (TRIM(cn_gamma_file), numGo)
            ENDIF
            CALL iom_get (numGo,jpdom_global,'Go',GammaO(:,:,:))

            temdta1(:,:,:)=rean_tem(:,:,:) ! reanalysis fields interpolated
            saldta1(:,:,:)=rean_sal(:,:,:) !  in model domain

            temdta2(:,:,:) = tsn(:,:,:,jp_tem)  ! temdta2= model temperature at time step
            saldta2(:,:,:) = tsn(:,:,:,jp_sal)  ! temdta2= model temperature at time step
! DIFFERENCE            
            diff_tem = (temdta1(:,:,:)-temdta2(:,:,:))*tmask(:,:,:) ! Parent minus child
            diff_sal = (saldta1(:,:,:)-saldta2(:,:,:))*tmask(:,:,:) ! Parent minus child

            DO kk=1,jpk  ! loop in each vertical level
! GATHER
               CALL gather_map(diff_tem(:,:,kk),tem_map)
               CALL gather_map(diff_sal(:,:,kk),sal_map)
               CALL gather_map(tmask(:,:,kk),mask_map)
! FILTER
               CALL space_filt_2d(tem_map,mask_map,5,tem_filt)
               CALL space_filt_2d(sal_map,mask_map,5,sal_filt)
! SCATTER
               CALL scatter_map(tem_filt,tem_par)
               CALL scatter_map(sal_filt,sal_par)
! NUDGING
!              nudging and nudge coefficient
               SN_T(:,:,kk) = tem_par(:,:)*GammaO(:,:,kk)*0.2*tmask(:,:,kk)
               SN_S(:,:,kk) = sal_par(:,:)*GammaO(:,:,kk)*0.2*tmask(:,:,kk)
!              update nemo field
               tsn(:,:,kk,jp_tem) = tsn(:,:,kk,jp_tem) + SN_T(:,:,kk)*tmask(:,:,kk)
               tsn(:,:,kk,jp_sal) = tsn(:,:,kk,jp_sal) + SN_S(:,:,kk)*tmask(:,:,kk)

            ENDDO

!           Output variable as a check:
            sntem = SN_T
            snsal = SN_S
!            sntem(:,:,1) = sf_sn(jp_salin)%fnow(:,:,1)
!             sntem = rean_tem

      ENDIF

      IF( kt == nitend ) THEN
         CALL iom_close (numGo)
      ENDIF

   END SUBROUTINE tra_sn_simple


!----------------------------------------------------------------
   SUBROUTINE gather_map(fldin,pngglo)
!----------------------------------------------------------------
      INTEGER  ::   ji, jj, jn                          ! dummy loop indices

      REAL(wp), DIMENSION(jpi,jpj)        :: fldin      ! input field in parallel
      REAL(wp), DIMENSION(jpi,jpj,jpnij)  :: png        ! parallel nodes gathered extra dimension
      REAL(wp), DIMENSION(jpiglo,jpjglo)  :: pngglo     ! nodes gathered in global map
!----------------------------------------------------------------

!??      CALL lbc_lnk( temdta1, 'T', 1. ) ! Ensure all haloes are filled in pn

      CALL mppgather (fldin(:,:),0,png(:,:,:))

      IF (nproc==0) THEN
         DO jn=1,jpnij
            DO jj=nldjt(jn),nlejt(jn)
               DO ji=nldit(jn),nleit(jn)
                  pngglo(ji+nimppt(jn)-1,jj+njmppt(jn)-1)=png(ji,jj,jn)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

   END SUBROUTINE gather_map

!-------------------------------------------------------------------------------------------------
   SUBROUTINE scatter_map(in_glo,out_par)
!--------------------------------------------------------------------------------------------------
      INTEGER  ::   ji, jj, jn                        ! dummy loop indices

      REAL(wp), DIMENSION(jpiglo,jpjglo)  :: in_glo   ! input field gathered in global domain
      REAL(wp), DIMENSION(jpi,jpj,jpnij)  :: png2     ! parallel nodes gathered in extra dimension
      REAL(wp), DIMENSION(jpi,jpj)        :: out_par  ! outut field in parallel
!---------------------------------------------------------------------------------------------------

!     Divide global domain in subdomains again, all together as an extra dimension
      IF (nproc==0) THEN
         DO jn=1,jpnij
            DO jj=nldjt(jn),nlejt(jn)
               DO ji=nldit(jn),nleit(jn)
                  png2(ji,jj,jn)=in_glo(ji+nimppt(jn)-1,jj+njmppt(jn)-1)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!      Scatter the nodes into parallel
      CALL mppsync ! like a "wait" function
         CALL mppscatter(png2,0,out_par)
      CALL mppsync

   END SUBROUTINE scatter_map


!---------------------------------------------------------------------------------------------------
   SUBROUTINE space_filt_2d(fldin,maskin,loop,fldout)
!---------------------------------------------------------------------------------------------------
      INTEGER            :: kk, jj, ji                    ! dummy loop indices
      INTEGER            :: iter, acount, bcount, cpass   ! dummy loop counters
      INTEGER            :: lx

      INTEGER            ::  loop                         ! number of iterations for the filter
      REAL(wp)           ::  alpha1, alpha2, alpha3, &    ! filter coefficients
                             beta1, beta2, beta3          !  (need to be provided manually!!!)
      REAL(wp), DIMENSION(2) ::  Inc                      ! dummy parameters for initialization
      REAL(wp), DIMENSION(jpiglo,jpjglo)  :: fldin        ! input map (2D field)
      REAL(wp), DIMENSION(jpiglo,jpjglo)  :: fldout       ! output filtered map (2D)
      REAL(wp), DIMENSION(jpiglo,jpjglo)  :: maskin       ! 2D mask for that specific vertical level
      REAL(wp), DIMENSION(jpiglo,jpjglo)  :: mask1        ! new 2D mask without borders

      REAL(wp), DIMENSION(jpiglo+(Lorder*2)) :: Tlargex   ! dummy parameters for filtering
      REAL(wp), DIMENSION(jpjglo+(Lorder*2)) :: Tlargey   ! the dimention is (filter order)*2
      REAL(wp), DIMENSION(jpiglo+(Lorder*2)) :: Tfx, Trx 
      REAL(wp), DIMENSION(jpjglo+(Lorder*2)) :: Tfy, Try 

      REAL(wp), DIMENSION(jpiglo,jpjglo) ::  TLout, Tgiv, TLh
!---------------------------------------------------------------------------------------------------

      cpass = loop ! number of iterations for the filter

      ! Filter coefficients
!      bx=(/0.1212,0.2425,0.1212/)
!      ax=(/1.0,-0.8030,0.2880/)

      alpha1=1.0
      alpha2=-0.8030
      alpha3=0.2880

      beta1=0.1212
      beta2=0.2425
      beta3=0.1212

      lx=3 ! size(alpha)
      Lorder=3*(lx-1)

!     fix the problem with initial x(1) and x(end) in the filter
!     based on octave and matlab flitflit
      Inc(1)=beta3-((beta1+beta2+beta3)/(alpha1+alpha2+alpha3))*alpha3 &
             +beta2-((beta1+beta2+beta3)/(alpha1+alpha2+alpha3))*alpha2
      Inc(2)=beta3-((beta1+beta2+beta3)/(alpha1+alpha2+alpha3))*alpha3

      ! Account for =0 border in nemo outputs
      fldin(1,:)=fldin(2,:)
      fldin(jpiglo,:)=fldin(jpiglo-1,:)
      fldin(:,1)=fldin(:,2)
      fldin(:,jpjglo)=fldin(:,jpjglo-1)

      mask1(:,:)=maskin(:,:)
      mask1(1,:)=mask1(2,:)
      mask1(jpiglo,:)=mask1(jpiglo-1,:)
      mask1(:,1)=mask1(:,2)
      mask1(:,jpjglo)=mask1(:,jpjglo-1)

      Tgiv(:,:) = fldin(:,:)*mask1(:,:) 

         DO iter=1,cpass ! iteration loop
                
            DO jj=1,jpjglo   !filtering in first dimension

               bcount=0 ! initialization
               DO acount= Lorder+1, 2, -1
                  bcount=bcount+1
                  Tlargex(bcount)=2*Tgiv(1,jj)-Tgiv(acount,jj)
               ENDDO
               Tlargex((bcount+1):(jpiglo+bcount))=Tgiv(:,jj)
               bcount=jpiglo+bcount
               DO acount= jpiglo-1,jpiglo-Lorder, -1
                  bcount=bcount+1
                  Tlargex(bcount)=2*Tgiv(jpiglo,jj)-Tgiv(acount,jj)
               ENDDO
               !Forward filtering
               Tfx(1)=beta1*Tlargex(1)+Inc(1)*Tlargex(1) !initial conditions
               Tfx(2)=beta1*Tlargex(2)+Inc(1)*Tlargex(1) ! "         "
               DO ji= 3, jpiglo+(Lorder*2)
                  Tfx(ji)=beta1*Tlargex(ji)+beta2*Tlargex(ji-1)+beta3*Tlargex(ji-2) &
                          -alpha2*Tfx(ji-1)-alpha3*Tfx(ji-2)
               ENDDO
               !Reverse filtering
               Tfx=Tfx(jpiglo+(Lorder*2):1:-1)
               Trx(1)=beta1*Tfx(1)+Inc(1)*Tfx(1) !initial conditions
               Trx(2)=beta1*Tfx(2)+Inc(1)*Tfx(1) !  "          "
               DO ji=3, jpiglo+(Lorder*2)
                  Trx(ji)=beta1*Tfx(ji)+beta2*Tfx(ji-1)+beta3*Tfx(ji-2) &
                          -alpha2*Trx(ji-1)-alpha3*Trx(ji-2)
              ENDDO
              Trx(:)=Trx(jpiglo+(Lorder*2):1:-1)
              TLh(:,jj)=Trx(1+Lorder:jpiglo+Lorder)
           ENDDO
           TLh(:,:)=TLh(:,:)*mask1(:,:)
 
           DO ji=1,jpiglo  !filtering in second dimension
              !initialization
              bcount=0
              DO acount= Lorder+1, 2, -1
                 bcount=bcount+1
                 Tlargey(bcount)=2*TLh(ji,1)-TLh(ji,acount)
              ENDDO
              Tlargey((bcount+1):(jpjglo+bcount))=TLh(ji,:)
              bcount=jpjglo+bcount
              DO acount= jpjglo-1,jpjglo-Lorder, -1
                 bcount=bcount+1
                 Tlargey(bcount)=2*TLh(ji,jpjglo)-TLh(ji,acount)
              ENDDO

             !Forward filtering
             Tfy(1)=beta1*Tlargey(1)+Inc(1)*Tlargey(1) !initial conditions
             Tfy(2)=beta1*Tlargey(2)+Inc(1)*Tlargey(1) ! "         "
             DO jj= 3, jpjglo+(Lorder*2)
                Tfy(jj)=beta1*Tlargey(jj)+beta2*Tlargey(jj-1)+beta3*Tlargey(jj-2) &
                        -alpha2*Tfy(jj-1)-alpha3*Tfy(jj-2)
             ENDDO
             !Reverse filtering
             Tfy=Tfy(jpjglo+(Lorder*2):1:-1)
             Try(1)=beta1*Tfy(1)+Inc(1)*Tfy(1) !initial conditions
             Try(2)=beta1*Tfy(2)+Inc(1)*Tfy(1) !  "
             DO jj= 3, jpjglo+(Lorder*2)
                Try(jj)=beta1*Tfy(jj)+beta2*Tfy(jj-1)+beta3*Tfy(jj-2) &
                        -alpha2*Try(jj-1)-alpha3*Try(jj-2)
             ENDDO
             Try(:)=Try(jpjglo+(Lorder*2):1:-1)
             TLout(ji,:)=Try(1+Lorder:jpjglo+Lorder)
          
           ENDDO
           Tgiv(:,:)=TLout(:,:)*mask1(:,:)
        ENDDO

        fldout = Tgiv

   END SUBROUTINE space_filt_2d


#endif
   !!======================================================================
END MODULE sn_simple
