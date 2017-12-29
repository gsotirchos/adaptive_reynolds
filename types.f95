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

end module types
