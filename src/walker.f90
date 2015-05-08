
!-----------------------------------------------------------------------
! 
!  walker.f90
!
!  Does basic random walks with arbitrary presision
!
!-----------------------------------------------------------------------


! okay... so number of steps
! start with 2 steps, x0=0 x0.5=gauss x1=x0.5+gauss
! now, with 4 steps I x0.25 = (x0+x0.5)/2 + gauss*sigma/sqrt(2)

! num_steps or depth

! tfinal
! sigma
! dconst

! generate
! sort
! find_crosses (should find the number of crosses, then make space, then return the size and array)

! so I want this to be consistent acruss runs, but multiple walkers might exist?
! or I might only need one walker at a time...

! the easiest thing to do is to take a seed,




module class_walker

  use globals
  use random_pl

  implicit none

  private
  integer                                   :: index


  type, public :: walker
     ! user defined variables
     integer                                :: num_points = 0
     integer                                :: depth = 0
     real(wp)                               :: dconst
     real(wp)                               :: sigma
     integer                                :: seed   

     ! managed internally
     real(wp), allocatable, dimension(:)    :: time
     real(wp), allocatable, dimension(:)    :: pos
     integer                                :: num_crossings = 0

   contains
     ! prototypes? sorta.
     procedure :: ctor => walker_ctor
     procedure :: allocate => walker_allocate
     procedure :: deallocate => walker_deallocate
     procedure :: generate => walker_generate
     procedure :: print => walker_print
  end type walker

contains

  ! walker_ctor constructor( logarithmic_binning, bin_start, bin_end, bin_number )
  subroutine walker_ctor(this, v1, v2, v3, v4) ! (depth, sigma, dconst)
    class(walker), intent(inout)  :: this
    integer, intent(in)              :: v1, v4
    real(wp), intent(in)             :: v2, v3

    depth = v1
    sigma = v2
    dconst = v3
    seed = v4

    if ( depth .ge. 1 ) then
       num_points = 2**depth
    else

    call this%allocate
    call this%init       

  end subroutine walker_ctor
  
  ! allocate
  subroutine walker_allocate(this)
    class(walker), intent(inout) :: this
    allocate(this%time(this%num_points))
    allocate(this%pos(this%num_points))

    ! initialize arrays
    this%time = 0
    this%pos = 0

  end subroutine walker_allocate

  ! deallocate
  subroutine walker_deallocate(this)
    class(walker), intent(inout) :: this
    deallocate(this%time)
    deallocate(this%pos)
  end subroutine walker_deallocate

  ! init
  subroutine walker_generate(this)
    class(walker), intent(inout) :: this

  end subroutine walker_generate

  ! walker_print
  subroutine walker_print(this)
    class(walker), intent(in) :: this

    write(*,*) '# bin# binstart binend count normalized_count'
        
    do bin_i = 1,this%bin_number
       write(*,*) bin_i, this%bin_boundaries(bin_i), this%bin_boundaries(bin_i+1), &
            this%bin_count(bin_i), real(this%bin_count(bin_i))/real(this%hist_count)
    end do
    
  end subroutine walker_print


end module class_walker


