! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras 

      use star_lib
      use star_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      
      implicit none

      integer :: time0, time1, clock_rate
      real(dp), parameter :: expected_runtime = 95 ! h
      logical :: reached_ZAMS = .false.      

      real(dp) :: Delta_Ma = 0.05d0
      real(dp) :: current_Ma_to_save 
      ! character(len=1024) :: Ma_profile_prefix
      character(len=1024) :: str1,Ma_profile_prefix
      integer istr,ls1,ls2
      real(dp) :: Ssurf0, Tsurf0, M0



      ! these routines are called by the standard run_star check_model
      contains
      
      ! include 'standard_run_star_extras.inc'
            subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         
         ! Uncomment these lines if you wish to use the functions in this file,
         ! otherwise we use a null_ version which does nothing.
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items


         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          s% job% warn_run_star_extras =.false.       
            
      end subroutine extras_controls
      
      ! None of the following functions are called unless you set their
      ! function point in extras_control.
      
      
      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: restart_time, prev_time_used
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0
         ! call system_clock(time0,clock_rate)
         if (.not. restart) then
            call system_clock(time0,clock_rate)
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
            call system_clock(restart_time,clock_rate)
            prev_time_used = time1 - time0
            time1 = restart_time
            time0 = time1 - prev_time_used
         end if

         Ssurf0 = s%entropy(1)
         Tsurf0 = s%T(1)
         M0 = s% star_mass
         current_Ma_to_save = s% star_mass - Delta_Ma

      end function extras_startup
      

      integer function extras_start_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s

         real(dp) :: dt

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         ! if (.false. .and. s% star_mass_h1 < 0.35d0) then
         !    ! stop when star hydrogen mass drops to specified level
         !    extras_check_model = terminate
         !    write(*, *) 'have reached desired hydrogen mass'
         !    return
         ! end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

        call system_clock(time1,clock_rate)
        dt = dble(time1-time0)/clock_rate / 60/60 ! hours



        if (dt > expected_runtime) then
          write(*,'(/,a30,2f18.6,a,/)') '>>>>>>> EXCESSIVE runtime', &
             dt, expected_runtime, '   <<<<<<<<<  ERROR'
          extras_check_model = terminate
          s% termination_code = t_xtra1
          termination_code_str(t_xtra1) = 'maximum runtime exceeded'
          return
        endif


         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 5   ! KDT: changed from 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
          ! logical :: sad_zone
         integer :: k
         real(dp) :: tau_th, mdot_crit_sad, mdot_crit_notsad, mout_sad, mout_notsad
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
        mdot_crit_sad = -15.0d0
        mdot_crit_notsad = -15.0d0

        mout_sad = 0.0d0
        mout_notsad = 0.0d0


        do k=1, s%nz
           if (s% grad_superad_actual(k) >= 0.1* s% actual_gradT(k)) then
              tau_th = sum(s% dm(1:k)*s% cp(1:k)*s% T(1:k))/s% L(1)
              tau_th = sum(s% dm(1:k))/tau_th/(Msun/secyer)
              if (safe_log10_cr(tau_th) > mdot_crit_sad) then
                 mdot_crit_sad = safe_log10_cr(tau_th)
                 mout_sad = sum(s% dm(1:k))/Msun
              endif
           else
              tau_th = sum(s% dm(1:k)*s% cp(1:k)*s% T(1:k))/s% L(1)
              tau_th = sum(s% dm(1:k))/tau_th/(Msun/secyer)
              ! mdot_crit_notsad = max(safe_log10_cr(tau_th), mdot_crit_notsad)
              if (safe_log10_cr(tau_th) > mdot_crit_notsad) then
                 mdot_crit_notsad = safe_log10_cr(tau_th)
                 mout_notsad = sum(s% dm(1:k))/Msun
              endif
           endif

           if (sum(s% dm(1:k))/Msun >= 0.1) then
              exit
           endif

        enddo

        names(1) = "log_mdot_crit_sad"
        vals(1) = mdot_crit_sad  ! in solar masses per year
        names(2) = "log_mdot_crit_notsad"
        vals(2) = mdot_crit_notsad  ! in solar masses per year

        names(3) = "ext_mass_sad"
        vals(3) = mout_sad  ! in solar masses per year
        names(4) = "ext_mass_notsad"
        vals(4) = mout_notsad  ! in solar masses per year
        
        ! s% xtra(2) = safe_log10_cr(s% r(1) /Rsun)
        ! s% xtra(3) = safe_log10_cr(s% star_mass /Msun)

        s% xtra(2) = s% lnR(1)
        s% xtra(3) = s% star_mass*ln10    ! ln(x) = log(x)*ln(10)

        names(5) = 'zeta'
        vals(5) = (s% xtra(2) - s% xtra_old(2)) / (s% xtra(3) - s% xtra_old(3))

        ! vals(5) = (s% xtra(2) - s% xtra_old(2))/&
        ! (s% xtra(3) - s% xtra_old(3))

        ! names(6) = ''
        ! (s% lnR(1) - s% lnR_start(1)) / &
        !              ((s% mstar - s% mstar_old)/(0.5d0*(s% mstar + s% mstar_old)))
         
         !note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.
         ! column 1
         

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns

      subroutine how_many_extra_history_header_items(id, id_extra, num_cols)
      integer, intent(in) :: id, id_extra
      integer, intent(out) :: num_cols
      num_cols=0
      end subroutine how_many_extra_history_header_items
      
      subroutine data_for_extra_history_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
      integer, intent(in) :: id, id_extra, num_extra_header_items
      character (len=*), pointer :: extra_header_item_names(:)
      real(dp), pointer :: extra_header_item_vals(:)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      !here is an example for adding an extra history header item
      !set num_cols=1 in how_many_extra_history_header_items and then unccomment these lines
      !extra_header_item_names(1) = 'mixing_length_alpha'
      !extra_header_item_vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_history_header_items


      subroutine how_many_extra_profile_header_items(id, id_extra, num_cols)
      integer, intent(in) :: id, id_extra
      integer, intent(out) :: num_cols
      num_cols = 0
      end subroutine how_many_extra_profile_header_items
      
      subroutine data_for_extra_profile_header_items( &
                  id, id_extra, num_extra_header_items, &
                  extra_header_item_names, extra_header_item_vals, ierr)
      integer, intent(in) :: id, id_extra, num_extra_header_items
      character (len=*), pointer :: extra_header_item_names(:)
      real(dp), pointer :: extra_header_item_vals(:)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      !here is an example for adding an extra profile header item
      !set num_cols=1 in how_many_extra_profile_header_items and then unccomment these lines
      !extra_header_item_names(1) = 'mixing_length_alpha'
      !extra_header_item_vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
         real(dp), parameter :: expected_runtime = 95 ! h
         ! integer :: time0=0, time1, clock_rate
         ! real(dp) :: dt

         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s

         integer :: f

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         extras_finish_step = keep_going
         call store_extra_info(s)


         ! MESA provides a number of variables that make it easy to get user input.
         ! these are part of the star_info structure and are named
         ! x_character_ctrl, x_integer_ctrl, x_logical_ctrl, and x_ctrl.
         ! by default there are num_x_ctrls, which defaults to 100, of each.
         ! they can be specified in the controls section of your inlist.

!          f = s% x_integer_ctrl(1)

         ! MESA also provides a number variables that are useful for implementing
         ! algorithms which require a state. if you just use these variables
         ! restarts, retries, and backups will work without doing anything special.
         ! they are named xtra1 .. xtra30, ixtra1 .. ixtra30, and lxtra1 .. lxtra30.
         ! they are automatically versioned, that is if you set s% xtra1, then
         ! s% xtra1_old will contains the value of s% xtra1 from the previous step
         ! and s% xtra1_older contains the one from two steps ago.

         s% xtra(1) = s% star_mass

        ! if ( ((s% L_nuc_burn_total / s% L_surf) > 0.9) .and. (.not. reached_ZAMS) ) then
        !     ! s% delta_HR_limit = 0.001d0
        !     ! extras_finish_step = terminate
        !     reached_ZAMS = .true.

        !     ! s% profile_interval = 3
        !     ! s% max_num_profile_models = 2000

        !     ! s% varcontrol_target = 1.0d-4
        !     ! s% delta_HR_limit = 1.0d-3

        !     ! write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        !     ! write(*,*) 'Reached (near) the ZAMS!'
        !     ! write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         
        ! endif

        ! if (s% c_core_mass > 0.1) then    ! on the AGB now!
        !     s% use_gold_tolerances = .false.
        !     ! s% do_conv_premix = .true.
        !     ! s% conv_premix_avoid_increase = .true. 

        !     if (s% initial_mass < 1.75) then
        !       s% convergence_ignore_equL_residuals = .false.
        !     endif

        !     if (s% initial_mass > 7.5) then
        !       ! s% varcontrol_target = 1.0d-3
        !       s% delta_HR_limit = 1.0d-2
        !     endif


        ! endif

        ! if ( (s% initial_mass < 1.75) .and. (s% power_he_burn > 1d-1) ) then  ! He burning is starting
        !   s% convergence_ignore_equL_residuals = .true.
        ! endif

        ! if (s%log_surface_radius .gt. log10(1.25)) then 
        !  ! write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        !  ! write(*,*) s%Teff, s%entropy(1)
        !  ! write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        !  ! s% need_to_save_profiles_now = .true.
        !  ! s% profile_data_prefix = 'profile_accretion_onset'

        !  s% use_other_surface_PT = .true.
        !  s% mass_change = 1d-5
        !  ! extras_finish_step = terminate
        ! endif

         if ((s% xtra(1) <= current_Ma_to_save) .and. (s% xtra_old(1) >= current_Ma_to_save) ) then
            write (Ma_profile_prefix, *) "Md_",current_Ma_to_save, ".mod"
            call star_write_model(id,Ma_profile_prefix, ierr)

            s% need_to_save_profiles_now = .true.
            write (Ma_profile_prefix, *) "profile_Md_",current_Ma_to_save, "_"
            s% profile_data_prefix = Ma_profile_prefix
            s% save_profiles_model_priority = 10

            current_Ma_to_save = current_Ma_to_save - Delta_Ma


            
         ! write(str1 , *) "profile_Md_",current_Ma_to_save, "_"
         ! ls1 = len_trim(str1)
         ! ls2 = 0
         ! do istr = 1,ls1
         !   if(str1(istr:istr).ne.' ') then
         !      ls2 = ls2 + 1
         !      Ma_profile_prefix(ls2:ls2) = str1(istr:istr)
         !   endif
         ! enddo

         ! ! s% need_to_save_profiles_now = .true.
         ! s% profile_data_prefix = str1
         ! ! s% save_profiles_model_priority = 10

         ! current_Ma_to_save = current_Ma_to_save + Delta_Ma
         
        endif

     ! if (current_Ma_to_save .lt. (M0 - 2.5d0)) then 
     !     extras_finish_step = terminate
     !  endif


         ! write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         ! write(*,*) current_Ma_to_save
         ! write(*,*) s% xtra_old(1)
         ! write(*,*) s% xtra(1)
         ! write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         ! extras_finish_step = terminate

         ! this expression will evaluate to true if f times the log center density
         ! has crossed an integer during the last step.  If f = 5, then we will get
         ! output at log center density = {... 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 ... }
         !if ((floor(f * s% xtra1_old) - floor(f * s% xtra1) .ne. 0)) then

            ! save a profile & update the history
          !  s% need_to_update_history_now = .true.
           ! s% need_to_save_profiles_now = .true.

            ! by default the priority is 1; you can change that if you'd like
            ! s% save_profiles_model_priority = ?

         !endif

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         ! call system_clock(time1,clock_rate)
         ! dt = dble(time1-time0)/clock_rate / 60/60 ! hours
         ! ! write(*,'(/,a30,2f18.6,a,/)') '>>>>>>> runtime', &
         ! !       dt*60*60, expected_runtime*60*60, '   <<<<<<<<<'

        
         ! if (dt > expected_runtime) then
         !    write(*,'(/,a30,2f18.6,a,/)') '>>>>>>> EXCESSIVE runtime', &
         !       dt, expected_runtime, '   <<<<<<<<<  ERROR'
         !    ! stop
         !    extras_finish_step = terminate
         !    return
         ! endif

         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring extra data so can do restarts
         
         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3
      
      
      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg    
         num_ints = i
         
         i = 0
         ! call move_dbl       
         
         num_dbls = i
         
         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      
      end subroutine move_extra_info

      
      end module run_star_extras
      
