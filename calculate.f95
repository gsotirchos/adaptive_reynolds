subroutine calculate(node, element, P)

    use types
    use parameters
    use invert

    implicit none

    type(fem_node),    intent(in)  :: node
    type(fem_element), intent(in)  :: element
    real,              intent(out) :: P(size(node%x))
    
    real :: l12, l23, l31, &
            C(2:3, 3), A, &
            KM(size(node%x), size(node%x)), m, &
            fc, f(size(node%x)), &
            KM_free(count(node%stat == 0), count(node%stat == 0)), &
            f_free(count(node%stat == 0)), &
            P_free(count(node%stat == 0))

    integer :: isfixed(count(node%stat == 1)), &
               isfree(count(node%stat == 0))


    KM = 0
    fc = 0
    f = 0

    do n = 1, size(element%node, dim = 1)

        ! Calculate triangles' sides' lengths
        l12 = dist([node%x(element%node(n, 1)),  &
                    node%y(element%node(n, 1))], &
                   [node%x(element%node(n, 2)),  &
                    node%y(element%node(n, 2))])

        l23 = dist([node%x(element%node(n, 2)),  &
                    node%y(element%node(n, 2))], &
                   [node%x(element%node(n, 3)),  &
                    node%y(element%node(n, 3))])

        l31 = dist([node%x(element%node(n, 3)),  &
                    node%y(element%node(n, 3))], &
                   [node%x(element%node(n, 1)),  &
                    node%y(element%node(n, 1))])

        ! Generate the 2nd and 3rd row of C matrix
        C(2, 1) = node%y(element%node(n, 2)) - node%y(element%node(n, 3))
        C(2, 2) = node%y(element%node(n, 3)) - node%y(element%node(n, 1))
        C(2, 3) = node%y(element%node(n, 1)) - node%y(element%node(n, 2))
        C(3, 1) = node%x(element%node(n, 3)) - node%x(element%node(n, 2))
        C(3, 2) = node%x(element%node(n, 1)) - node%x(element%node(n, 3))
        C(3, 3) = node%x(element%node(n, 2)) - node%x(element%node(n, 1))

        ! Calculate element area
        ! If the triangles are orthogonal then C is singular
        ! and the area is calculated by width*height/2
        A = (C(2, 1)*C(3, 2) - C(2, 2)*C(3, 1))/2
        !if (A == 0) A = (l12*l31)/2

        ! Calculate mass matrix element m
        m = lambda*A/12


        ! Generate stiffness matrix
        do j = 1, 3
            nj = element%node(n, j)
            do i = 1, 3
                ni = element%node(n, i)

                KM(ni, nj) = KM(ni, nj) + A*(  hx*C(2, i)*C(2, j)   &
                                             + hy*C(3, i)*C(3, j)  )

                KM(ni, nj) = KM(ni, nj) - m
            end do

            KM(nj, nj) = KM(nj, nj) - m
        end do


        ! Generate f matrix
        n1 = element%node(n, 1)
        n2 = element%node(n, 2)
        n3 = element%node(n, 3)
        fc = element%c(n)*A/3

        f(n1) = f(n1) - fc + (l12*element%q(n, 1))/2 &
                           + (l31*element%q(n, 3))/2
        f(n2) = f(n2) - fc + (l12*element%q(n, 1))/2 &
                           + (l23*element%q(n, 2))/2
        f(n3) = f(n3) - fc + (l23*element%q(n, 2))/2 &
                           + (l31*element%q(n, 3))/2

    end do


    ! Generate isfree and isfixed matrices
    n1 = 0
    n2 = 0

    do i = 1, size(node%stat)
        if (node%stat(i) == 1) then
            n1 = n1 + 1
            isfixed(n1) = i
        else if (node%stat(i) == 0) then
            n2 = n2 + 1
            isfree(n2) = i
        end if
    end do


    ! Delete fixed degrees of freedom rows and move the columns
    ! to the right side of the equation
    KM_free = 0
    f_free = 0

    do i = 1, size(isfree)
        do j = 1, size(isfree)
            KM_free(i, j) = KM(isfree(i), isfree(j))
        end do

        f_free(i) = f(isfree(i))
        do j = 1, size(isfixed)
            f_free(i) = f_free(i) - KM(i, isfixed(j))*node%P(isfixed(j))
        end do
    end do


    ! Generate P matrix
    P_free = matmul(inv(KM_free), f_free)
    P = node%P

    do i = 1, size(P_free)
        P(isfree(i)) = P_free(i)
    end do


    !! DEBUG
    !print *, "** DEBUG **"
    !print *,
    !print *, "STIFFNESS MATRIX"
    !print "(16F7.1)", KM
    !print *,
    !print *, "INVERTED (FREE) STIFFNESS MATRIX"
    !print "(12F7.1)", inv(KM_free)
    !print *,
    !print *, "f MATRIX"
    !print "(F8.4)", f
    !print *,
    print *, "p matrix"
    print "(f8.4)", p
    print *,
    !print *, "Free:", isfree
    !print *, "Fixed:", isfixed

end subroutine calculate
