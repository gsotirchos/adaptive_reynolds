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

    ! subroutine to deallocate nodes
    subroutine dealloc_node(node)

        implicit none

        type(fem_node) :: node

        deallocate(node%x)
        deallocate(node%y)
        deallocate(node%stat)
        deallocate(node%P)
        deallocate(node%H)

    end subroutine dealloc_node

    ! subroutine to allocate elements
    subroutine alloc_elem(element, n, elem_nodes)

        implicit none

        integer, intent(in) :: n, elem_nodes
        type(fem_element)   :: element

        allocate(element%node(n, elem_nodes))
        allocate(element%q(n, elem_nodes))
        allocate(element%He(n))
        allocate(element%c(n))

    end subroutine alloc_elem

    ! subroutine to deallocate nodes
    subroutine dealloc_elem(element)

        implicit none

        type(fem_element) :: element

        deallocate(element%node)
        deallocate(element%q)
        deallocate(element%He)
        deallocate(element%c)

    end subroutine dealloc_elem


    ! Function to calculate distance between two points
    pure function dist(point1, point2)

        implicit none

        real, intent(in)  :: point1(2), point2(2)
        real :: dist

        dist = sqrt(  (point2(1) - point1(1))**2   &
                    + (point2(2) - point1(2))**2  )

    end function dist

end module types
