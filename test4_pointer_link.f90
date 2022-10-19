program test4
    ! linked_list
    implicit none
    type :: real_value
        real :: value
        type(real_value), pointer :: p
    end type
    type(real_value), pointer :: head
    character(len = 20) :: filename
    integer :: nvals = 0
    type(real_value), pointer :: ptr
    type(real_value), pointer :: tail
    integer :: istat
    real :: temp

    write(*,*) "enter the file containing the data to read:"
    read(*, '(A20)') filename
    open(unit=9, file=filename, status='old', action='read', iostat=istat)
    fileopen: if(istat == 0)then
        input: do
        read
end program test4