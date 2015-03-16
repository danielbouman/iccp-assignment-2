subroutine hello(a,b)
    implicit none
    
    integer, intent(in) :: a
    integer, intent(out) :: b
    
    b = a*2
end subroutine