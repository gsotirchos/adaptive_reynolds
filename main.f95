program reynolds_main

    use types
    implicit none

    type(fem_node) :: node
    type(fem_element) :: element

    real :: hx, hy, Ha, Hb

    integer :: num_nodes, num_elem

    ! input data
    call input(num_nodes, num_elem, node, element, hx, hy, Ha, Hb)

    ! calculation
    call calculation(num_nodes, num_elem, node, element, hx, hy, Ha, Hb)

end program reynolds_main
