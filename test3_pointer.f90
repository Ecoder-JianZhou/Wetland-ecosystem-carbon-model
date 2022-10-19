program test3
    implicit none
    ! call sub_pointer_ptr4()
    ! call sub_pointer_test1()
    ! call sub_pointer_test0()
    call sub_pointer_ptr5()
end 

subroutine sub_pointer()
    implicit none
    integer, pointer :: p1, p2, p3, p4
    integer, target :: t1 = 1
    integer :: t2 = 2, t3 =3
    target :: t2
    p1 => t1
    p2 => t2
    ! p3 => t3 ! target is neither target nor pointer ...
    p4 => p1
    print *, p1, p2, p3, p4
    ! p1 = 10 ! not change ?
    p1 => t2
    print *, p1, p2, p3, p4
end subroutine sub_pointer

subroutine sub_pointer_ptr4()
    implicit none
    real, pointer :: p1, p2, p3
    real, target :: a = 11., b=12.5, c
    Nullify(p1, p2, p3) ! nullify pointer p1 = null()
    p1 => a
    p2 => b
    p3 => c
    p3 = p1 + p2
    write(*,*) "p3 = ", p3
    p2 => p1
    p3 = p1 + p2
    write(*,*) "p3 = ", p3
    p3 = p1
    p3 => p1
    write(*,*) "p3 = ", p3
    write(*,*) "a, b, c = ", a, b, c
end subroutine sub_pointer_ptr4

subroutine sub_pointer_test0()
    implicit none
    real, pointer :: p1
    real, target :: t1 = 1
    p1 => t1
    write(*,*) "test1:", t1
    p1 = 2
    write(*,*) "test2:", t1
end subroutine sub_pointer_test0

subroutine sub_pointer_test1() ! just change the value of pointers
    implicit none
    real, pointer :: p1, p2, temp
    real, target :: t1, t2
    t1 = 1.5
    t2 = 2.5
    p1 => t1
    p2 => t2
    temp => p1
    write(*,*) "test1: t1, t2, p1, p2, temp =", t1, t2, p1, p2, temp
    p1 => p2
    p2 => temp
    write(*,*) "test2: t1, t2, p1, p2, temp =", t1, t2, p1, p2, temp
end subroutine sub_pointer_test1

subroutine sub_pointer_ptr5()
    implicit none
    integer :: i
    integer, dimension(16), target :: info = (/(i, i=1,16)/)
    integer, dimension(:), pointer :: p1, p2, p3, p4, p5
    p1 => info
    p2 => p1(2::2)
    p3 => p2(2::2)
    p4 => p3(2::2)
    p5 => p4(2::2)
    write(*,'(A,16I3)') "p1 = ", p1
    write(*,'(A,16I3)') "p2 = ", p2
    write(*,'(A,16I3)') "p3 = ", p3
    write(*,'(A,16I3)') "p4 = ", p4
    write(*,'(A,16I3)') "p5 = ", p5
end subroutine sub_pointer_ptr5