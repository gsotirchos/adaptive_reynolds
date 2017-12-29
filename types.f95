module types
    
    type el_node
        real, allocatable :: x(:)
        real, allocatable :: y(:)
        integer, allocatable :: stat(:) ! 0: prescribed, 1: free
        real, allocatable :: P(:) ! nodal prescribed dimensionless pressure
        real, allocatable :: H(:) ! nodal dimensionless height
    end type el_node

    type tri_element
        integer, allocatable :: n1(:)
        integer, allocatable :: n2(:)
        integer, allocatable :: n3(:)
        real, allocatable :: q12(:) ! flow on 1-2 side
        real, allocatable :: q23(:) ! flow on 2-3 side
        real, allocatable :: q31(:) ! flow on 3-1 side
        real, allocatable :: He(:) ! average H of element
        real, allocatable :: c(:) ! parameter c of element
    end type tri_element

end module types
