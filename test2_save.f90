! save: save the value 
program test2
    implicit none
    call sub1()
    call sub1()
    call sub2()
    call sub2()
    
end 
subroutine sub1()
        implicit none
        integer :: a = 1
        integer :: b = 2
        integer c
        c=a+b
        a=c+1
        b=c+2
        write(*,*)"sub1:", a, b 
    end subroutine sub1

    subroutine sub2()
        implicit none
        integer a, b
        ! integer :: b = 2
        integer c
        save a, b
        a = 1
        b = 2
        c=a+b
        a=c+1
        b=c+2
        write(*,*)"sub2:", a, b 
    end subroutine sub2