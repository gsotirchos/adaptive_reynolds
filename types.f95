module types

    ! Define Data Types for Nodes and Elements
    
    type fem_node
        real, allocatable :: x(:)
        real, allocatable :: y(:)
        integer, allocatable :: stat(:) ! 0: free, 1: prescribed
        real, allocatable :: P(:) ! nodal prescribed dimensionless pressure
        real, allocatable :: H(:) ! nodal dimensionless height
    end type fem_node

    type fem_element
        integer, allocatable :: node(:, :) ! nodes of element
        real, allocatable :: q(:, :) ! flow on side
        real, allocatable :: He(:) ! average H of element
        real, allocatable :: c(:) ! parameter c of element
    end type fem_element

    ! Interface the two types of "extend" subroutines
    interface extend
        module procedure ExtendArray, ExtendMatrix
    end interface

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

    ! Subroutine to add "len" number of rows after "loc" location
    ! in a matrix
    subroutine ExtendMatrix(A, loc, len)

        implicit none

        real, allocatable   :: A(:, :), temp(:, :)
        integer             :: n, m
        integer, intent(in) :: loc, len

        n = size(A, dim = 1)
        m = size(A, dim = 2)
        allocate(temp(n, m))


        if (loc > n .or. loc < 0) then
            print *, "Error: can't extend matrix on given location"
            return
        end if

        if (len + loc < 0) then
            print *, "Error: can't extend matrix with given length"
            return
        end if


        temp = A

        deallocate(A)
        allocate(A((n + len), m))

        A(:loc, :) = temp(:loc, :)
        A((loc + len + 1):, :) = temp((loc + 1):, :)
        A((loc + 1):(loc + len), :) = 0

        deallocate(temp)

    end subroutine ExtendMatrix

    ! Subroutine to add "len" number of rows after "loc" location
    ! in an array
    subroutine ExtendArray(A, loc, len)

        implicit none

        real, allocatable   :: A(:), temp(:)
        integer             :: n
        integer, intent(in) :: loc, len

        n = size(A)
        allocate(temp(n))


        if (loc > n .or. loc < 0) then
            print *, "Error: can't extend matrix on given location"
            return
        end if

        if (len + loc < 0) then
            print *, "Error: can't extend matrix with given length"
            return
        end if


        temp = A

        deallocate(A)
        allocate(A(n + len))

        A(:loc) = temp(:loc)
        A((loc + len + 1):) = temp((loc + 1):)
        A((loc + 1):(loc + len)) = 0

        deallocate(temp)

    end subroutine ExtendArray

end module types
