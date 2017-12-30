program reynolds_main

    use types
    implicit none

    type(fem_node) :: node
    type(fem_element) :: element

    integer :: num_nodes, num_elem

    real :: hx, hy, Ha, Hb
    real, allocatable :: P(:)

    ! input data
    call input(num_nodes, num_elem, node, element, hx, hy, Ha, Hb)

    allocate(P(num_nodes))

    ! calculation
    call calculation(num_nodes, num_elem, node, element, hx, hy, Ha, Hb, P)

end program reynolds_main
