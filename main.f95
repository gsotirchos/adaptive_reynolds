program reynolds_main

    use types
    use parameters

    implicit none

    type(fem_node)    :: node
    type(fem_element) :: element
    real, allocatable :: P(:)

    allocate(P(num_nodes))

    ! input data
    call input(node, element)

    ! calculation
    call calculation(node, element, P)

end program reynolds_main
