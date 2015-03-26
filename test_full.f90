program test
    implicit none
    
    real(8), dimension(3,2) :: pos
    real(8), dimension(5,2) :: candidate_pos
    real(8), dimension(5) :: energies, possib_angles
    real(8) :: sigma_squared, epsilon, bend_energy, last_angle
    integer :: N_candidates, N_existing
    
    sigma_squared = 0.64
    epsilon = 0.25
    N_candidates = 5
    N_existing = 3
    pos = reshape((/ 0, 1, 2, 3, 4, 5 /), shape(pos))
    candidate_pos = reshape((/ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 /), shape(candidate_pos))
    possib_angles = reshape((/ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 /), shape(possib_angles))
    bend_energy = 10
    last_angle = 0
    
    
    call func(pos,candidate_pos,N_candidates,N_existing,sigma_squared,epsilon,bend_energy,last_angle,possib_angles,energies)
    
    print *, energies
end program test

subroutine func(pos,candidate_pos,N_candidates,N_existing,sigma_squared,epsilon,bend_energy,last_angle,possib_angles,energies)
    implicit none
    
    integer, intent(in) :: N_candidates, N_existing
    real(8), intent(in) :: sigma_squared, epsilon, last_angle, bend_energy
    
    real(8), dimension(N_existing,2), intent(in) :: pos
    real(8), dimension(N_candidates,2), intent(in) :: candidate_pos
    real(8), dimension(N_candidates), intent(in) :: possib_angles
    real(8), dimension(N_candidates), intent(out) :: energies
    real(8) :: V
    
    real(8) :: abs_distance_squared
    integer :: i,j
    
    do i=1,N_candidates
        do j=1,N_existing-1
            abs_distance_squared = sum((candidate_pos(i,:) - pos(j,:))**2) / sigma_squared
            V = 4*epsilon*( (abs_distance_squared)**(-6) - (abs_distance_squared)**(-3))
            energies(i) = energies(i)+V
        end do
        energies(i) = energies(i)+bend_energy*cos(0.5*(possib_angles(i)-last_angle))
    end do
    
end subroutine