program reynolds_main

    use types

    implicit none

    type(fem_node)       :: node
    type(fem_element)    :: element
    real,    allocatable :: P(:)
    integer, allocatable :: split(:)


    ! input data
    call input(node, element)

    allocate(P(size(node%x)))

    ! calculation
    call calculation(node, element, P)

    allocate(split(size(element%node, dim = 1)))

    call evaluation(node, element, P, split)

end program reynolds_main
