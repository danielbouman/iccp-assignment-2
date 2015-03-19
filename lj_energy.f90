subroutine func(pos,candidate_pos,N_candidates,N_existing,sigma_squared,epsilon,energies)
    implicit none
    
    integer, intent(in) :: N_candidates, N_existing
    real(8), intent(in) :: sigma_squared, epsilon
    
    real(8), dimension(N_existing,2), intent(in) :: pos
    real(8), dimension(N_candidates,2), intent(in) :: candidate_pos
    real(8), dimension(N_candidates), intent(out) :: energies
    real(8) :: V
    
    real(8) :: abs_distance_squared
    integer :: i,j
    
    do i=1,N_candidates
        do j=1,N_existing
            abs_distance_squared = sum((candidate_pos(i,:) - pos(j,:))**2) / sigma_squared
            V = 4*epsilon*( (abs_distance_squared)**(-6) - (abs_distance_squared)**(-3))
            energies(i) = energies(i)+V
        end do
    end do
    
end subroutine