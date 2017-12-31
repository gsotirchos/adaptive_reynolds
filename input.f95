subroutine input(node, element)

    use types
    use parameters
    implicit none

    type(fem_node), intent(out) :: node
    type(fem_element), intent(out) :: element

    ! Allocate arrays
    allocate(node%x(num_nodes))
    allocate(node%y(num_nodes))
    allocate(node%stat(num_nodes))
    allocate(node%P(num_nodes))
    allocate(node%H(num_nodes))

    allocate(element%node(num_elem, elem_nodes))
    allocate(element%q(num_elem, elem_nodes))
    allocate(element%He(num_elem))
    allocate(element%c(num_elem))

    ! Generate node positions
    n = 0

    do j = 0, ysub
        do i = 0, xsub
            n = n + 1
            node%x(n) = i*dL
            node%y(n) = j*dB
        end do
    end do


    ! Generate elements' nodes
    n = -1

    do j = 1, ysub
        do i = 1, xsub
        n = n + 2
        element%node(n, 1) = (j - 1)*(xsub + 1) + i
        element%node(n, 2) = element%node(n, 1) + 1
        element%node(n, 3) = element%node(n, 2) + xsub

        element%node(n + 1, 1) = element%node(n, 3) + 1
        element%node(n + 1, 2) = element%node(n, 3)
        element%node(n + 1, 3) = element%node(n, 2)
        end do
    end do


    ! Generate boundary values
    ! Here: at nodes where x = 0 or L or where y = B, p is prescribed
    where(node%y == B .or. node%x == 0 .or. node%x == L)
        node%stat = 0
        node%P = P_bound
    elsewhere
        node%stat = 1
        node%p = 0
    end where


    ! Generate boundary flow
    ! Element pairs are defined as follows
    !   3     2_____1 
    !   |\     \    |
    !   | \     \(2)|
    !   |  \     \  |    y
    !   |(1)\     \ |   |
    !  1|____\2    \|3  |___x
    ! horizontal border with y = C -> sides "1-2"
    !   vertical border with x = C -> sides "3-1"
    ! Here: on elements' sides with y = 0, q is nonzero
    ! Because dp/dy = 0 -> q = dp/dx (on y = 0)

    where(node%y(element%node(:, 1)) == 0 .and. &
          node%y(element%node(:, 2)) == 0)
        element%q(:, 1) = q
        element%q(:, 2) = 0
        element%q(:, 3) = 0
    elsewhere
        element%q(:, 1) = 0
        element%q(:, 2) = 0
        element%q(:, 3) = 0
    end where


    ! Generate nodal H
    ! H is described by an equation of the form:
    !     H = Ha*x + Hb
    ! where Ha = dH/dx = const.
    node%H = Ha*node%x + Hb


    ! Generate elements' He
    element%He = (node%H(element%node(:, 1)) + &
                  node%H(element%node(:, 2)) + &
                  node%H(element%node(:, 3)))/3


    ! Generate elements' c
    element%c = 6*Ha/(element%He**3)

    ! DEBUG
    !print *, "** DEBUG **"
    !print *, "NODES"
    !print "(A5, 5A5)", " ", "x", "y", "stat", "P", "H"
    !do i = 0, (num_nodes + 1)
    !    print "(I4, A1, 2F5.2, I5, 2F5.2)", i, ": ", &
    !        node%x(i), node%y(i), node%stat(i), node%P(i), node%H(i)
    !end do
    !print *, " "

    !print *, "ELEMENTS"
    !print "(A5, 3A4, 5A6)", " ", "n1", "n2", "n3", &
    !                        "q12", "q23", "q31", &
    !                        "He", "c"
    !do i = 0, (num_elem + 1)
    !    print "(I4, A1, 3I4, 5F6.2)", i, ": ", &
    !        element%node(i, 1), element%node(i, 2), element%node(i, 3), &
    !        element%q(i, 1), element%q(i, 2), element%q(i, 3), &
    !        element%He(i), element%c(i)
    !end do

end subroutine input
