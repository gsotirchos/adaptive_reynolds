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
            K_full(num_nodes,num_nodes), &
            f_full(num_nodes), fc, &
            K(num_nodes,num_nodes), &
            f(num_nodes), &
            invK(num_nodes,num_nodes), &
            P(num_nodes), &
    ! array to map "reduced" matrices' rows to originals'
    integer :: isfixed(count(node%stat == 0)) 

    ! Generate triangles' sides' lengths matrices
    l12 = sqrt((node%x(element%node(:,2))-node%x(element%node(:,1)))**2 + &
               (node%y(element%node(:,2))-node%y(element%node(:,1)))**2)
    l23 = sqrt((node%x(element%node(:,3))-node%x(element%node(:,2)))**2 + &
               (node%y(element%node(:,3))-node%y(element%node(:,2)))**2)
    l31 = sqrt((node%x(element%node(:,1))-node%x(element%node(:,3)))**2 + &
               (node%y(element%node(:,1))-node%y(element%node(:,3)))**2)

    ! Generate the 2nd and 3rd row of C matrix for each element
    C(:, 2, 1) = node%y(element%node(:, 2)) - node%y(element%node(:, 3))
    C(:, 2, 2) = node%y(element%node(:, 3)) - node%y(element%node(:, 1))
    C(:, 2, 3) = node%y(element%node(:, 1)) - node%y(element%node(:, 2))
    C(:, 3, 1) = node%y(element%node(:, 3)) - node%y(element%node(:, 2))
    C(:, 3, 2) = node%y(element%node(:, 1)) - node%y(element%node(:, 3))
    C(:, 3, 3) = node%y(element%node(:, 2)) - node%y(element%node(:, 1))

    ! Calculate element area
    ! If the triangles are orthogonal then C is singular
    ! and the area is calculated by width*height/2
    A = (C(:, 2, 1)*C(:, 3, 2) - C(:, 2, 2)*C(:, 3, 1))/2

    do n = 1, num_elem
        if (A(n) == 0) A(n) = (l12(n)*l31(n))/2
        C(n, :, :) = C(n, :, :)/(2*A(n))
    end do


    K = 0
    K_full = 0

    ! Generate full stiffness matrix
    do n = 1, num_elem
        do j = 1, 3
            do i = 1, 3
                ni = element%node(n, i)
                nj = element%node(n, j)
                K_full(ni, nj) = K_full(ni, nj) + &
                                 A(n)*(hx*C(n, 2, i)*C(n, 2, j) + &
                                       hy*C(n, 3, i)*C(n, 3, j))
            end do
        end do
    end do


    f = 0
    f_full = 0
    fc = 0

    ! Generate full f matrix
    do n = 1, num_elem
        n1 = element%node(n, 1)
        n2 = element%node(n, 2)
        n3 = element%node(n, 3)
        fc = element%c(n)*A(n)/3
        f_full(n1) = f_full(n1) - fc + (l12(n)*element%q(n, 1) + &
                                        l31(n)*element%q(n, 3))/2
        f_full(n2) = f_full(n2) - fc + (l12(n)*element%q(n, 1) + &
                                        l23(n)*element%q(n, 2))/2
        f_full(n3) = f_full(n3) - fc + (l23(n)*element%q(n, 2) + &
                                        l31(n)*element%q(n, 3))/2
    end do


    ! Generate map matrix
    n = 0
    do i = 1, size(node%stat)
        if (node%stat(i) == 0) then
            n = n + 1
            isfixed(n) = i
        end if
    end do

    !Re
    do i = 1, size(isfixed)
        do j = 1, size(isfixed)
            K(i, j) = K_full(isfixed(i), isfixed(j))
        end do

        do j = 1, 
        f(i) = 
    end do



    ! Generate P matrix
    invK = inv(K)
    P = matmul(invK, f)

    ! DEBUG
    print *, "** DEBUG **"
    print *,
    print *, "STIFFNESS MATRIX"
    do i = 1, num_nodes
        print "(16F4.1)", K_full(i, :)
    end do
    print *,
    print *, "INVERTED STIFFNESS MATRIX"
    do i = 1, num_nodes
        print "(16F4.1)", invK(i, :)
    end do
    print *,
    print *, "f MATRIX"
    print "(F4.1)", f_full
    print *,
    print *, "P MATRIX"
    print "(F4.1)", P

    print *, map
    print *, node%stat

end subroutine calculation
