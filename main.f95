program reynolds_main

    use types
    implicit none
    type(el_node) :: node
    type(tri_element) :: element

    ! input data
    call input(node, element)

end program reynolds_main
