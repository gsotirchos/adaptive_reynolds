subroutine evaluation(node, element, P, split)

    use types
    use parameters

    implicit none

    type(fem_node),       intent(in)  :: node
    type(fem_element),    intent(in)  :: element
    real,                 intent(in)  :: P(size(node%x))
    integer,              intent(out) :: split(size(element%node, dim = 1))

    real    :: C(2:3, 3), &
               l12, l23, l31, &
               A(size(element%node, dim = 1)), &
               elem_dP(size(element%node, dim = 1), 2), &
               node_dP(size(node%x), 2), &
               estm_dP(size(element%node, dim = 1), 2), &
               err_norm(size(element%node, dim = 1)), &
               mean_err, var_err, stdev_err

    integer :: cnt(size(node%x))


    elem_dP = 0
    cnt = 0
    node_dP = 0

    do n = 1, size(elem_dP, dim = 1)

        ! Generate the 2nd and 3rd row of C matrix
        C(2, 1) = node%y(element%node(n, 2)) - node%y(element%node(n, 3))
        C(2, 2) = node%y(element%node(n, 3)) - node%y(element%node(n, 1))
        C(2, 3) = node%y(element%node(n, 1)) - node%y(element%node(n, 2))
        C(3, 1) = node%x(element%node(n, 3)) - node%x(element%node(n, 2))
        C(3, 2) = node%x(element%node(n, 1)) - node%x(element%node(n, 3))
        C(3, 3) = node%x(element%node(n, 2)) - node%x(element%node(n, 1))

        ! Calculate element's Pressure 1st Derivatives
        do i = 1, size(element%node, dim = 2)
            elem_dP(n, 1) = elem_dP(n, 1) + C(2, i)*P(element%node(n, i))
            elem_dP(n, 2) = elem_dP(n, 2) + C(3, i)*P(element%node(n, i))
        end do

        ! Apply element's contribution on its nodes' derivatives' estimation
        do i = 1, size(element%node, dim = 2)
            node_dP(element%node(n, i), 1) = node_dP(element%node(n, i), 1)&
                                             + elem_dP(n, 1)
            node_dP(element%node(n, i), 2) = node_dP(element%node(n, i), 2)&
                                             + elem_dP(n, 2)

            cnt(element%node(n, i)) = cnt(element%node(n, i)) + 1
        end do

        ! Calculate triangles' sides' lengths
        l12 = sqrt(  (  node%x(element%node(n, 2))         &
                      - node%x(element%node(n, 1))  )**2   &
                   + (  node%y(element%node(n, 2))         &
                      - node%y(element%node(n, 1))  )**2  )

        l23 = sqrt(  (  node%x(element%node(n, 3))         &
                      - node%x(element%node(n, 2))  )**2   &
                   + (  node%y(element%node(n, 3))         &
                      - node%y(element%node(n, 2))  )**2  )

        l31 = sqrt(  (  node%x(element%node(n, 1))         &
                      - node%x(element%node(n, 3))  )**2   &
                   + (  node%y(element%node(n, 1))         &
                      - node%y(element%node(n, 3))  )**2  )

        ! Calculate element's area
        A(n) = (C(2, 1)*C(3, 2) - C(2, 2)*C(3, 1))/2
        if (A(n) == 0) A(n) = (l12*l31)/2

    end do

    node_dP(:, 1) = node_dP(:, 1)/cnt
    node_dP(:, 2) = node_dP(:, 2)/cnt

    estm_dP = 0
    err_norm = 0

    do n = 1, size(elem_dP, dim = 1)
        do i = 1, size(element%node, dim = 2)
            estm_dP(n, 1) = estm_dP(n, 1) + node_dP(element%node(n, i), 1)/3
            estm_dP(n, 2) = estm_dP(n, 2) + node_dP(element%node(n, i), 2)/3
        end do

        err_norm(n) = sqrt((  (estm_dP(n, 1) - elem_dP(n, 1))**2   &
                            + (estm_dP(n, 2) - elem_dP(n, 2))**2  )*A(n))
    end do


    ! Calculate Mean Error Norm
    mean_err = sum(err_norm)/size(err_norm)

    ! Calculate Error Norm Variance and Standard Deviation
    var_err = sum(err_norm**2)
    var_err = (var_err - size(err_norm)*mean_err**2)/(size(err_norm) - 1)

    stdev_err = sqrt(var_err)


    ! Assign splits
    split = 0

    do n = 1, size(element%node, dim = 1)
        do i = 1, 7
            if (      err_norm(n) < mean_err + i*0.5*stdev_err &
                .and. err_norm(n) > mean_err + (i-1)*0.5*stdev_err) then
            split(n) = i
            exit
            end if
        end do
    end do


    ! DEBUG
    !print *, "** DEBUG **"
    print *,
    print *, "Element dP"
    print *, "/dx"
    print "(6F6.2)", elem_dP(:, 1)
    print *, "/dy"
    print "(6F6.2)", elem_dP(:, 2)
    print *,
    print *, "Node dP"
    print *, "/dx"
    print "(4F6.2)", node_dP(:, 1)
    print *, "/dy"
    print "(4F6.2)", node_dP(:, 2)
    print *,
    print *, "Element Estimated dP"
    print *, "/dx"
    print "(6F6.2)", estm_dP(:, 1)
    print *, "/dy"
    print "(6F6.2)", estm_dP(:, 2)
    print *,
    print *, "Error Norm"
    print "(6F6.2)", err_norm
    print *,
    print *, "Mean:", mean_err
    print *, "Variance", var_err
    print *, "Standard Deviation", stdev_err
    print *,
    print *, "Splits"
    print "(I1)", split

end subroutine evaluation
