subroutine func(a,size,b)
    implicit none
   
    integer, intent(in) :: size
    real(8), dimension(size,2), intent(in) :: a
    real(8), dimension(size), intent(out) :: b
    
    b = a(:,2)**(-1)
    
end subroutine