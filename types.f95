module types

    ! Define Data Types for Nodes and Elements
    
    type fem_node
        real, allocatable :: x(:)
        real, allocatable :: y(:)
        integer, allocatable :: stat(:) ! 0: prescribed, 1: free
        real, allocatable :: P(:) ! nodal prescribed dimensionless pressure
        real, allocatable :: H(:) ! nodal dimensionless height
    end type fem_node

    type fem_element
        integer, allocatable :: node(:, :) ! nodes of element
        real, allocatable :: q(:, :) ! flow on side
        real, allocatable :: He(:) ! average H of element
        real, allocatable :: c(:) ! parameter c of element
    end type fem_element

contains

    ! subroutine to allocate nodes
    subroutine alloc_node(node, n)

        implicit none

        integer, intent(in) :: n
        type(fem_node)      :: node

        allocate(node%x(n))
        allocate(node%y(n))
        allocate(node%stat(n))
        allocate(node%P(n))
        allocate(node%H(n))

    end subroutine alloc_node

    ! subroutine to allocate element
    subroutine alloc_elem(element, n, elem_nodes)

        implicit none

        integer, intent(in) :: n, elem_nodes
        type(fem_element)   :: element

        allocate(element%node(n, elem_nodes))
        allocate(element%q(n, elem_nodes))
        allocate(element%He(n))
        allocate(element%c(n))

    end subroutine alloc_elem

end module types
