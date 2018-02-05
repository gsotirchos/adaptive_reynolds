subroutine gen_nodes(node)

    use types
    use parameters

    implicit none

    type(fem_node), intent(out) :: node

    call alloc_node(node, num_nodes)

    ! Generate node positions
    ! Element pairs are defined as follows
    !   3     2_____1 
    !   |\     \    |
    !   | \     \(2)|
    !   |  \     \  |    y
    !   |(1)\     \ |   |
    !  1|____\2    \|3  |___x
    n = 0

    do j = 0, ysub
        do i = 0, xsub
            n = n + 1

            node%x(n) = i*dL
            node%y(n) = j*dB
        end do
    end do


    !! DEBUG
    !print *, "** DEBUG **"
    !print *,
    !print *, "NODES"
    !print "(A5, 5A5)", " ", "x", "y", "stat", "P", "H"
    !do i = 0, (size(node%x) + 1)
    !    print "(I4, A1, 2F5.2, I5, 2F5.2)", i, ": ", &
    !        node%x(i), node%y(i), node%stat(i), node%P(i), node%H(i)
    !end do
    !print *,

end subroutine gen_nodes
