!-----------------------------
!
!  weinercross.f90
!
! This program is a simple simulation of brownian diffusion with
! the intent of measuring the zero-crossing intervals. Also will
! be able to measure trapped and free intervals for other
! boundaries or derived quantities like energy.
!
! How will it work? The first version will simply be x=x+dx where
! dx = dt * N(0,sigma^2) where N is a a sample from a normal
! distribution. Each time step we check to see if it has
! changed state, and if so, store a value of stuff. Then, 
!
!-----------------------------

! num_steps pretty much gives the accuracy of a given walk
! num_walks gives the accuracy 

program weinercross

  use globals
  use random_pl
  use class_histogram
  use utilities
  
  implicit none
  
  ! declare variables to use below
  real(wp)               :: t = 0.0_wp
  real(wp)               :: dt, dx
  real(wp)               :: x = 0.0_wp
  real(wp)               :: area = 0.0_wp
  real(wp)               :: sigma = 1.0_wp
  real(wp)               :: dconst
  integer(kind=int16)    :: i_step,i_walk, num_steps, num_walks, j, k
  
  real(wp)               :: trap_timer = 0
  integer                :: trap_status = 0
  
  
  character(256)         :: format1, buff
  integer(kind=int16)    :: percent_complete = 0
  
  integer(rpk) :: seed             !      -> Seed for RNG
  integer, parameter     :: rands_num = 10**7
  integer                :: rands_index = 0
  real(wp), dimension(rands_num) :: rands_container
  integer                :: discrete_steps = 0
  
  ! crossing variables
  integer                :: crossings_max
  integer                :: crossings_num = 0
  type(histogram)        :: crossings_hist
  integer                :: crossings_logarithmic_binning = 1
  real(wp)               :: crossings_bin_start
  real(wp), allocatable, dimension(:) :: crossings_list
  real(wp)               :: crossings_bin_end
  integer                :: crossings_bin_num
  
  
  ! declare variables for histogram bins
  ! things about this section: it's for a specific axis (radial direction used)
  type(histogram) :: hist
  integer         :: logarithmic_binning = 1
  real(wp)        :: bin_start
  real(wp)        :: bin_end
  integer         :: bin_num
  
  
  ! declare variables for histogram bins
  ! things about this section: it's for a specific axis (radial direction used)
  type(histogram) :: hist2
  integer         :: logarithmic_binning2 = 0
  real(wp)        :: bin_start2 = -10.0_wp
  real(wp)        :: bin_end2   = 10.0_wp
  integer         :: bin_num2   = 100

  ! printing x
  integer         :: print_xt_pairs = 0
    
  ! Input related variables
  character(len=100) :: label, ctrl_file
  integer :: pos
  integer, parameter :: fh = 15
  integer :: ios = 0
  integer :: line = 0
  integer :: min_arguments = 1 !
  
  
  if ( command_argument_count() .ne. min_arguments )  call usage()
  call get_command_argument(1, ctrl_file)
  
!!!!!! Get arguments from the command line
  open(fh, file=ctrl_file)
  
  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.
  
  do while (ios == 0)
     read(fh, '(A)', iostat=ios) buff
     if (ios == 0) then
        line = line + 1
        
        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buff, '    ')
        label = buff(1:pos)
        buff = buff(pos+1:)
        
        
        select case (label)
           
        case('num_steps')
           num_steps = s2i(buff)
           dt = 1.0_wp/real(num_steps, wp)
           sigma = sqrt(dt)
           write(0,*) 'Read ', label, ': ', num_steps
           
        case('num_walks')
           num_walks = s2i(buff)
           write(0,*) 'Read ', label, ': ', num_walks

        case('discrete_steps')
           discrete_steps = s2i(buff)
           write(0,*) 'Read ', label, ': ', discrete_steps

        case('print_xt_pairs')
           print_xt_pairs = s2i(buff)
           write(0,*) 'Read ', label, ': ', print_xt_pairs

!axe it           
        case('crossings_max')
           crossings_max = s2i(buff)
           allocate(crossings_list(crossings_max))
           write(0,*) 'Read ', label, ': ', crossings_max
           
        case('seed')
           seed = s2i(buff)
           write(0,*) 'Read ', label, ': ', seed
           call init_rand_pl(seed1=seed)
           
        case('logarithmic_binning')
           logarithmic_binning = s2i(buff)
           write(0,*) 'Read ', label, ': ', logarithmic_binning
           
        case('bin_start')
           bin_start = s2r(buff)
           write(0,*) 'Read ', label, ': ', bin_start
           
        case('bin_end')
           bin_end = s2r(buff)
           write(0,*) 'Read ', label, ': ', bin_end
           
        case('bin_num')
           bin_num = s2i(buff)
           write(0,*) 'Read ', label, ': ', bin_num
       
    !axe these
        case('crossings_bin_start')
           crossings_bin_start = s2r(buff)
           write(0,*) 'Read ', label, ': ', crossings_bin_start
           
        case('crossings_bin_end')
           crossings_bin_end = s2r(buff)
           write(0,*) 'Read ', label, ': ', crossings_bin_end
           
        case('crossings_bin_num')
           crossings_bin_num = s2i(buff)
           write(0,*) 'Read ', label, ': ', crossings_bin_num
           
        case('crossings_logarithmic_binning')
           crossings_logarithmic_binning = s2i(buff)
           write(0,*) 'Read ', label, ': ', crossings_logarithmic_binning
           
        case default
           write(0,*) 'Unknown label, Read ', label, ': ', buff
           
        end select
     end if
  end do
  
  write (0,*) "sigma: ", sigma
  write (0,*) "dt: ", dt
  
  
  
  
  ! set up hist
  call hist%ctor(logarithmic_binning, bin_start, bin_end, bin_num)
  call hist2%ctor(logarithmic_binning2, bin_start2, bin_end2, bin_num2)
  call crossings_hist%ctor(crossings_logarithmic_binning, crossings_bin_start, crossings_bin_end, crossings_bin_num)
  
  do i_walk = 1_int16, num_walks

     ! percent complete printout
     if( (100_int16*i_walk/num_walks) .gt. percent_complete ) then
        percent_complete = 100_int16*i_walk/num_walks
        write(0,*) "percent complete: ", percent_complete, " i_walk/num_walks: ", i_walk, "/", num_walks
        ! the flush command forces the stderr to be unbuffered
        flush(0)
        ! call fsync(fnum(10)) ! this turned out to not be nessessary, I think.
     end if

     do i_step = 1_int16, num_steps
             
        ! handle random normals
        rands_index = rands_index + 1
        ! reset random counter if index is one greater than the max
        if ( rands_index .eq. rands_num + 1 ) then
           rands_index = 1
        end if
        ! refill container if we just reset the index to 1 (also covers first loop)
        if ( rands_index .eq. 1 ) then
           rands_container = nrand_pl_vec(rands_num)
        end if
        
        
        ! step forward
        dx = sigma * rands_container(rands_index)
        if( discrete_steps .eq. 1 ) then
           if( dx .ge. 0_wp ) then
              dx = 1.0_wp
           else 
              dx = -1.0_wp
           end if
        end if
        t = t + dt
        x = x + dx
        area = area + x*dt
                
        
        ! handle recording crossing intervals 
        if ( ((trap_status .eq. 0) .and. ( x .gt. 0.0_wp )) .or. &
             ((trap_status .eq. 1) .and. ( x .lt. 0.0_wp )) ) then
        
           ! switch trap state
           trap_status = mod(trap_status+1,2)
           
           ! store interval
           call hist%add_item(t-trap_timer)         
           
           ! printing xt pairs:
           if ( print_xt_pairs .eq. 1 ) then
              if ( i_step .eq. 1 ) then
                 write(*,*) ' '
                 write(*,*) ' '
                 write(*,*) '"xt_pairs"'
              end if
              write(*,*) t-trap_timer, abs(area)
           end if

           ! reset timer
           trap_timer = t
           area = 0

           !
           crossings_num = crossings_num + 1

        end if
           
     end do

     ! testing endpoint histogram:
     call hist2%add_item(x)      
     x = 0
     t = 0
     trap_timer = 0

  end do
  
  ! print simulation information
  
  write(0,*) 'Crossing Interval Stats:'
  write(0,*) 'histogram entries skipped (no matching bins): ', hist%bin_skipped
  write(0,*) 'histogram entries skipped randoms (no matching bins): ', hist2%bin_skipped
  write(0,*) 'histogram entries skipped crossings (no matching bins): ', crossings_hist%bin_skipped
  write(0,*) 'num_crossings: ', crossings_num
  

  write(*,*) ' '
  write(*,*) ' '

  write(*,*) '"cross_interval"'
  call hist%print
  
  write(*,*) ' '
  write(*,*) ' '
  
  write(*,*) '"rand_test_final_position"'
  call hist2%print
  
  write(*,*) ' '
  write(*,*) ' '
  
contains
  
  subroutine usage()
    ! write (0,*) means write to standard error (0), using default
    !   formatting (*)
    write(0,*) ''
    write(0,*) ''
    write(0,*) 'Usage: '
    write(0,*) ' Rabi        ->  scaled Rabi frequency'
    write(0,*) ' Delta_base  ->  base detning = omega-omega0'
    stop
  end subroutine usage
  
end program weinercross

