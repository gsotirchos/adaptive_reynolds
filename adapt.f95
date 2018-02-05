subroutine adapt(node, element, splits)

    use types
    use parameters

    implicit none

    type(fem_node)      :: node
    type(fem_element)   :: element
    integer, intent(in) :: splits(size(element%c))

    integer :: remain

    ! Split elements in reverse order because every time a new element
    ! is added each one after it is shifted by 1
    do n = size(splits), 1, -1
        if (splits(n) == 1) then

            ! Add three new nodes in the middles of n-th element's sides
            j = ubound(node%x, dim = 1)
            call extend(node%x, j, 3)
            call extend(node%y, j, 3)

            do n1 = 1, 3
                n2 = mod(n1 + 3, 3) + 1

                node%x(j + n1) = (  node%x(element%node(n, n1)) &
                                  + node%x(element%node(n, n2)))/2
                node%y(j + n1) = (  node%y(element%node(n, n1)) &
                                  + node%y(element%node(n, n2)))/2
            end do

        end if
    end do

    ! Reallocate the rest of the nodes' properties
    deallocate(node%stat)
    deallocate(node%P)
    deallocate(node%H)

    allocate(node%stat(size(node%x)))
    allocate(node%P(size(node%x)))
    allocate(node%H(size(node%x)))

    ! Deallocate old elements
    call dealloc_elem(element)

end subroutine adapt
