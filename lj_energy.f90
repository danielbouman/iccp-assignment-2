! lj_energy.f90
subroutine lj_energy(possible_pos,N_possible_pos,pos,N_beads,epsilon,sigma_squared, energies)
    implicit none

    ! input variables
    integer, intent(in) :: N_possible_pos, N_beads
    real(8), intent(in) :: epsilon, sigma_squared
    real(8), intent(in), dimension(N_possible_pos,2) :: possible_pos
    real(8), intent(in), dimension(N_beads,2):: pos
    
    ! local variables
    real(8) :: abs_distance_squared, V
    integer :: i, j
    
    ! output variables
    real(8), intent(out), dimension(N_possible_pos) :: energies
    
    ! calculate energies
    do i=0,N_possible_pos
        do j=0,N_beads
            abs_distance_squared = ( sum( (possible_pos(i,:) - pos(j,:))**2 ) )/sigma_squared
            V = 4*epsilon*( (abs_distance_squared)**(-6) - (abs_distance_squared)**(-3))
            energies(i) = energies(i)+V
        end do
    end do
     
end subroutine lj_energy