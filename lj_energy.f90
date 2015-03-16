subroutine lj_energy(possible_positions,N_beads,epsilon,sigma)
    implicit none
    
    ! input variables
    integer, intent(in) :: N_beads, possible_positions
    real(8), intent(in) :: epsilon, sigma

    ! output variables
    real(8), intent(out), dimension(size(N_possible_positions)) :: energies
    
    ! local variables
    real(8) :: abs_distance_squared
    
    ! calculate energies
    do i=0,size(N_possible_positions)-1
        do j=0,size(N_beads)-1
            abs_distance_squared 
        end do
    end do
    
end subroutine lj_energy