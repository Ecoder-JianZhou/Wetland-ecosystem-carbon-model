
module testMod
    implicit none
    real x, y, a, b
    common a, b
    ! b = 1
    contains
    subroutine mod_a()
        implicit none
        x = 3.
        y = 4.
        ! call mod_c(a, b)
    end subroutine mod_a

    subroutine mod_b()
        implicit none
        ! x = 1.
        ! y = 2.
        write(*,*) x+y
    end subroutine mod_b

    subroutine mod_c()
        implicit none
        a = 2.
        b = 3.
    end subroutine mod_c
end module testMod

program test1
    Use testMod
    implicit none
    ! real :: testA
    ! integer :: a = 1
    ! type :: abc
    !     real gg
    ! end type abc
    ! namelist /b/ testA, gg

    ! call test(a, b)
    ! write(*, nml = b)

    ! contains
    ! subroutine test(a, b)
    !     integer a
    !     if(a == 1) then
    !         namelist /b/ testA
    !     endif
    ! return
    ! end subroutine test
    write(*,*)"res1--x, y:", x, y
    call mod_a
     write(*,*)"res2--x, y:", x, y
     x = 10.
     y = 11.
    call mod_b
     write(*,*)"res3--x, y:", x, y
end program test1