program lj_energy
    implicit none
    
    ! input variables
    integer :: pos, possible_pos
    real(8) :: epsilon, sigma
    
    ! local variables
    real(8) :: abs_distance_squared
    integer :: i,j, N_possible_positions

    ! output variables
    real(8), dimension(size(N_possible_positions)) :: energies
    
    ! calculate energies
    do i=0,size(N_possible_positions)-1
        do j=0,size(N_beads)-1
            abs_distance_squared = 
        end do
    end do
    
end subroutine lj_energ