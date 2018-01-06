subroutine calculation(node, element, P)

    use types
    use parameters
    use invert
    implicit none

    type(fem_node), intent(in) :: node
    type(fem_element), intent(in) :: element
    
    real :: l12(num_elem), &
            l23(num_elem), &
            l31(num_elem), &
            A(num_elem), C(num_elem, 2:3, 3), &
            K(num_nodes, num_nodes), &
            f(num_nodes), fc, &
            K_free(count(node%stat == 1), count(node%stat == 1)), &
            f_free(count(node%stat == 1)), &
            P_free(count(node%stat == 1)), P(num_nodes)

    integer :: isfixed(count(node%stat == 0)), &
               isfree(count(node%stat == 1))


    ! Generate triangles' sides' lengths matrices
    l12 = 0
    l12 = sqrt((node%x(element%node(:,2))-node%x(element%node(:,1)))**2 + &
               (node%y(element%node(:,2))-node%y(element%node(:,1)))**2)
    l23 = 0
    l23 = sqrt((node%x(element%node(:,3))-node%x(element%node(:,2)))**2 + &
               (node%y(element%node(:,3))-node%y(element%node(:,2)))**2)
    l31 = 0
    l31 = sqrt((node%x(element%node(:,1))-node%x(element%node(:,3)))**2 + &
               (node%y(element%node(:,1))-node%y(element%node(:,3)))**2)

    ! Generate the 2nd and 3rd row of C matrix for each element
    C(:, 2, 1) = node%y(element%node(:, 2)) - node%y(element%node(:, 3))
    C(:, 2, 2) = node%y(element%node(:, 3)) - node%y(element%node(:, 1))
    C(:, 2, 3) = node%y(element%node(:, 1)) - node%y(element%node(:, 2))
    C(:, 3, 1) = node%x(element%node(:, 3)) - node%x(element%node(:, 2))
    C(:, 3, 2) = node%x(element%node(:, 1)) - node%x(element%node(:, 3))
    C(:, 3, 3) = node%x(element%node(:, 2)) - node%x(element%node(:, 1))

    ! Calculate element area
    ! If the triangles are orthogonal then C is singular
    ! and the area is calculated by width*height/2
    A = (C(:, 2, 1)*C(:, 3, 2) - C(:, 2, 2)*C(:, 3, 1))/2

    where (A == 0)
        A = (l12*l31)/2
    end where


    ! Generate stiffness matrix
    K = 0
    do n = 1, num_elem
        print "(16F7.1)", K
        print *,
        do j = 1, 3
            do i = 1, 3
                ni = element%node(n, i)
                nj = element%node(n, j)
                K(ni, nj) = K(ni, nj) + A(n)*(hx*C(n, 2, i)*C(n, 2, j) + &
                                              hy*C(n, 3, i)*C(n, 3, j))
            end do
        end do
    end do


    ! Generate f matrix
    f = 0
    fc = 0
    do n = 1, num_elem
        n1 = element%node(n, 1)
        n2 = element%node(n, 2)
        n3 = element%node(n, 3)
        fc = element%c(n)*A(n)/3
        f(n1) = f(n1) - fc + (l12(n)*element%q(n, 1))/2 + &
                             (l31(n)*element%q(n, 3))/2
        f(n2) = f(n2) - fc + (l12(n)*element%q(n, 1))/2 + &
                             (l23(n)*element%q(n, 2))/2
        f(n3) = f(n3) - fc + (l23(n)*element%q(n, 2))/2 + &
                             (l31(n)*element%q(n, 3))/2
    end do


    ! Generate isfree and isfixed matrices
    n1 = 0
    n2 = 0
    do i = 1, size(node%stat)
        if (node%stat(i) == 0) then
            n1 = n1 + 1
            isfixed(n1) = i
        else if (node%stat(i) == 1) then
            n2 = n2 + 1
            isfree(n2) = i
        end if
    end do


    ! Delete fixed degrees of freedom rows and move the columns
    ! to the right side of the equation
    K_free = 0
    f_free = 0
    do i = 1, size(isfree)
        do j = 1, size(isfree)
            K_free(i, j) = K(isfree(i), isfree(j))
        end do

        f_free(i) = f(isfree(i))
        do j = 1, size(isfixed)
            f_free(i) = f_free(i) - K(i, isfixed(j))*node%P(isfixed(j))
        end do
    end do


    ! Generate P matrix
    P_free = matmul(inv(K_free), f_free)
    P = node%P
    do i = 1, size(P_free)
        P(isfree(i)) = P_free(i)
    end do


    ! DEBUG
    print *, "** DEBUG **"
    print *,
    print *, "STIFFNESS MATRIX"
    print "(16F7.1)", K
    print *,
    print *, "INVERTED (FREE) STIFFNESS MATRIX"
    print "(12F7.1)", inv(K_free)
    print *,
    print *, "f MATRIX"
    print "(F8.4)", f
    print *,
    print *, "P MATRIX"
    print "(F8.4)", P
    print *,
    print *, "Free:", isfree
    print *, "Fixed:", isfixed

end subroutine calculation
