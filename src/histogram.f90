
!-----------------------------------------------------------------------
! 
!  histogram.f90
!
!  Does basic histograms.
!
!-----------------------------------------------------------------------


module class_histogram
  use globals

  implicit none

  private
  integer                                   :: bin_i 
  real(wp)                                  :: base                  ! spacing between bins (for both logarithmic and otherwise)
  integer                                   :: found_bin             ! marker for whether a bin was found for an entry

  type, public :: histogram
     ! user defined variables
     integer                                :: logarithmic_binning
     real(wp)                               :: bin_start
     real(wp)                               :: bin_end
     integer                                :: bin_number
   
     ! managed internally
     real(wp), allocatable, dimension(:)    :: bin_boundaries
     integer, allocatable, dimension(:)     :: bin_count             ! number of entries in bin
     integer                                :: hist_count            ! number of entries total in histogram
     integer                                :: bin_skipped = 0       ! number of skipped entries

   contains
     ! prototypes? sorta.
     procedure :: ctor => histogram_ctor
     procedure :: allocate => histogram_allocate
     procedure :: deallocate => histogram_deallocate
     procedure :: init => histogram_init
     procedure :: add_item => histogram_add_item
     procedure :: print => histogram_print
  end type histogram

contains

  ! histogram_ctor constructor( logarithmic_binning, bin_start, bin_end, bin_number )
  subroutine histogram_ctor(this, v1, v2, v3, v4)
    class(histogram), intent(inout)  :: this
    integer, intent(in)              :: v1, v4
    real(wp), intent(in)                 :: v2, v3
    
    this%logarithmic_binning=v1
    this%bin_start=v2
    this%bin_end=v3
    this%bin_number=v4
    
    call this%allocate
    call this%init       

  end subroutine histogram_ctor
  
  ! allocate
  subroutine histogram_allocate(this)
    class(histogram), intent(inout) :: this
    allocate(this%bin_boundaries(this%bin_number+1))
    allocate(this%bin_count(this%bin_number))
    this%bin_count = 0
  end subroutine histogram_allocate

  ! deallocate
  subroutine histogram_deallocate(this)
    class(histogram), intent(inout) :: this
    deallocate(this%bin_boundaries)
    deallocate(this%bin_count)
  end subroutine histogram_deallocate

  ! init
  subroutine histogram_init(this)
    class(histogram), intent(inout) :: this

    ! calculate bin spacing
    if ( this%logarithmic_binning .eq. 1 ) then 
       base = exp( log( this%bin_end/this%bin_start ) / real(this%bin_number) )
    else
       base = (this%bin_end - this%bin_start) / real(this%bin_number)
    end if

    ! define bin boundaries
    do bin_i = 1,this%bin_number+1
       if ( this%logarithmic_binning .eq. 1 ) then
          this%bin_boundaries(bin_i) = this%bin_start*base**(real(bin_i-1))
       else
          this%bin_boundaries(bin_i) = this%bin_start + base*real(bin_i-1)
       end if
    end do
    
  end subroutine histogram_init

  ! add_item
  subroutine histogram_add_item(this,item)
    class(histogram), intent(inout) :: this
    real(wp), intent(in)         :: item

    found_bin = 0
    
    ! check each bin dumbly
    do bin_i = 1,this%bin_number
       if (( item .ge. this%bin_boundaries(bin_i) ) .and. ( item .lt. this%bin_boundaries(bin_i+1) )) then
          this%bin_count(bin_i) = this%bin_count(bin_i) + 1
          this%hist_count = this%hist_count + 1
          found_bin = 1
          exit
       endif
    end do
  
    if (found_bin .eq. 0) then
       this%bin_skipped = this%bin_skipped + 1
    endif
  end subroutine histogram_add_item

  ! histogram_print
  subroutine histogram_print(this)
    class(histogram), intent(in) :: this

    write(*,*) '# bin# binstart binend count normalized_count'
        
    do bin_i = 1,this%bin_number
       write(*,*) bin_i, this%bin_boundaries(bin_i), this%bin_boundaries(bin_i+1), &
            this%bin_count(bin_i), real(this%bin_count(bin_i))/real(this%hist_count)
    end do
    
  end subroutine histogram_print


end module class_histogram


