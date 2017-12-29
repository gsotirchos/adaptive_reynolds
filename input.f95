subroutine input(node, element)

    use types
    implicit none

    real :: hx, hy, P_bound, q, &
            Ha, Hb, &
            L, B, &
            dL, dB

    integer :: xsub, ysub, &
               num_nodes, num_elem, &
               tot_df, &
               i, j, cnt

    integer, parameter :: node_df = 1, ele_nodes = 3

    type(el_node), intent(out) :: node
    type(tri_element), intent(out) :: element

    ! Reynolds Parameters
    hx = 1.
    hy = 1.2
    p_bound = 0.
    q = 5.
    Ha = 0.5
    Hb = 4.

    ! Dimensions
    L = 3.
    B = 3.
    xsub = 3
    ysub = 3
    dL = L/xsub
    dB = B/ysub
    num_nodes = (xsub + 1)*(ysub + 1)
    num_elem = 2*xsub*ysub
    tot_df = node_df*num_nodes

    ! Allocate Arrays
    allocate(node%x(num_nodes))
    allocate(node%y(num_nodes))
    allocate(node%stat(num_nodes))
    allocate(node%P(num_nodes))
    allocate(node%H(num_nodes))

    allocate(element%n1(num_elem))
    allocate(element%n2(num_elem))
    allocate(element%n3(num_elem))
    allocate(element%q12(num_elem))
    allocate(element%q23(num_elem))
    allocate(element%q31(num_elem))
    allocate(element%He(num_elem))
    allocate(element%c(num_elem))

    ! Generate Node Positions
    cnt = 0

    do j = 0, ysub
        do i = 0, xsub
            cnt = cnt + 1
            node%x(cnt) = i*dL
            node%y(cnt) = j*dB
        end do
    end do


    ! Generate Elements' Nodes
    cnt = -1

    do j = 1, ysub
        do i = 1, xsub
        cnt = cnt + 2
        element%n1(cnt) = (j - 1)*(xsub + 1) + i
        element%n2(cnt) = element%n1(cnt) + 1
        element%n3(cnt) = element%n2(cnt) + xsub

        element%n1(cnt + 1) = element%n2(cnt)
        element%n2(cnt + 1) = element%n3(cnt) + 1
        element%n3(cnt + 1) = element%n3(cnt)
        end do
    end do


    ! Generate Boundary Values
    ! Here: at nodes where x = 0 or L or where y = B, p is prescribed.
    where(node%y == B .or. node%x == 0 .or. node%x == L)
        node%stat = 0
        node%P = P_bound
    elsewhere
        node%stat = 1
        node%p = 0
    end where


    ! Generate Boundary Flow
    ! Element pairs are defined as follows
    !   3     3_____2 
    !   |\     \    |
    !   | \     \(2)|
    !   |  \     \  |    y
    !   |(1)\     \ |   |
    !  1|____\2    \|1  |___x
    ! horizontal side with y = 0 -> sides "1-2" of of odd numbered elements
    ! horizontal side with y = B -> sides "2-3" of of even numbered elements
    !   vertical side with x = 0 -> sides "3-1" of of odd numbered elements
    !   vertical side with x = L -> sides "1-2" of of even numbered elements
    ! Here: on elements' sides with y = 0, q is nonzero.
    ! Because dp/dy = 0 -> q = dp/dx (on y = 0).

    where(node%y(element%n1) == 0 .and. node%y(element%n2) == 0)
        element%q12 = q
        element%q23 = 0
        element%q31 = 0
    elsewhere
        element%q12 = 0
        element%q23 = 0
        element%q31 = 0
    end where


    ! Generate Nodal H
    ! H is described by an equation of the form:
    !     H = Ha*x + Hb
    ! where Ha = dH/dx = const.
    node%H = Ha*node%x + Hb


    ! Generate Elements' He
    element%He = (node%H(element%n1) + &
                  node%H(element%n2) + &
                  node%H(element%n3))/3


    ! Generate Elements' c
    element%c = 6*Ha/(element%He**3)

    ! DEBUG
    ! print nodes
    print *, "** DEBUG **"
    print *, "NODES"
    print "(A5, 5A5)", " ", "x", "y", "stat", "P", "H"
    do i = 0, (num_nodes + 1)
        print "(I4, A1, 2F5.2, I5, 2F5.2)", i, ": ", &
            node%x(i), node%y(i), node%stat(i), node%P(i), node%H(i)
    end do
    print *, " "

    ! print elements
    print *, "ELEMENTS"
    print "(A5, 3A4, 5A6)", " ", "n1", "n2", "n3", &
                            "q12", "q23", "q31", &
                            "He", "c"
    do i = 0, (num_elem + 1)
        print "(I4, A1, 3I4, 5F6.3)", i, ": ", &
            element%n1(i), element%n2(i), element%n3(i), &
            element%q12(i), element%q23(i), element%q31(i), &
            element%He(i), element%c(i)
    end do

end subroutine input
