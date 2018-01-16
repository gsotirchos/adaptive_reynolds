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
        remain = splits(n)
        j = -1

        do while (remain > 0)
            j = j + 1

            do i = 3**j - 1, 0, -1
                remain = remain - 1

                ! Split (n + i)-th element
                call tr_split(node, element, n + i)

                if (remain == 0) exit
            end do

        end do

    end do


    ! DEBUG
    print *, "** DEBUG **"
    print *,
    print *, "NODES"
    print "(A5, 5A5)", " ", "x", "y", "stat", "P", "H"
    do i = 0, (size(node%x) + 1)
        print "(I4, A1, 2F5.2, I5, 2F5.2)", i, ": ", &
            node%x(i), node%y(i), node%stat(i), node%P(i), node%H(i)
    end do
    print *,
    print *, "ELEMENTS"
    print "(A5, 3A4, 5A6)", " ", "n1", "n2", "n3", &
                            "q12", "q23", "q31", &
                            "He", "c"
    do i = 0, (size(element%c) + 1)
        print "(I4, A1, 3I4, 5F6.2)", i, ": ", &
            element%node(i, 1), element%node(i, 2), element%node(i, 3), &
            element%q(i, 1), element%q(i, 2), element%q(i, 3), &
            element%He(i), element%c(i)
    end do
    print *,


contains

    subroutine tr_split(nodes, elements, n)

        implicit none

        integer, intent(in) :: n

        type(fem_node)    :: nodes, ntemp
        type(fem_element) :: elements, etemp
        integer           :: i, n1, n2, new_node


        ! Extend current nodes
         call alloc_node(ntemp, size(nodes%x))

        ntemp%x    = nodes%x
        ntemp%y    = nodes%y
        ntemp%stat = nodes%stat
        ntemp%P    = nodes%P
        ntemp%H    = nodes%H

        call dealloc_node(nodes)
        call alloc_node(nodes, size(ntemp%x) + 1)

        nodes%x(:size(ntemp%x))       = ntemp%x
        nodes%y(:size(ntemp%y))       = ntemp%y
        nodes%stat(:size(ntemp%stat)) = ntemp%stat
        nodes%P(:size(ntemp%P))       = ntemp%P
        nodes%H(:size(ntemp%H))       = ntemp%H

        ! Generate new node attributes
        nodes%x(ubound(nodes%x)) = sum(nodes%x(elements%node(n, :)))/3
        nodes%y(ubound(nodes%y)) = sum(nodes%y(elements%node(n, :)))/3
        nodes%stat(ubound(nodes%stat)) = 0
        nodes%P(ubound(nodes%P)) = sum(nodes%P(elements%node(n, :)))/3
        nodes%H(ubound(nodes%H)) = sum(nodes%H(elements%node(n, :)))/3

        call dealloc_node(ntemp)


        ! Extend current elements
        call alloc_elem(etemp, size(elements%c), 3)

        etemp%node(:size(elements%c), :) = elements%node
        etemp%q(:size(elements%c), :)    = elements%q
        etemp%He(:size(elements%c))      = elements%He
        etemp%c(:size(elements%c))       = elements%c

        call dealloc_elem(elements)
        call alloc_elem(elements, size(etemp%c) + 2, 3)

        elements%node(1:(n - 1), :) = etemp%node(1:(n - 1), :)
        elements%q(1:(n - 1), :) = etemp%q(1:(n - 1), :)
        elements%He(1:(n - 1)) = etemp%He(1:(n - 1))
        elements%c(1:(n - 1)) = etemp%c(1:(n - 1))

        elements%node((n + 3):, :) = etemp%node((n + 1):, :)
        elements%q((n + 3):, :) = etemp%q((n + 1):, :)
        elements%He((n + 3):) = etemp%He((n + 1):)
        elements%c((n + 3):) = etemp%c((n + 1):)

        new_node = size(nodes%x)

        ! Generate new elements' attributes
        ! the new elements are the n-th, (n+1)-th and (n+2)-th
        do n1 = 1, 3
            n2 = mod(n1 + 3, 3) + 1

            elements%node(n+(n1-1), 1) = etemp%node(n, n1)
            elements%node(n+(n1-1), 2) = etemp%node(n, n2)
            elements%node(n+(n1-1), 3) = new_node

            elements%q(n+(n1-1), 1) = etemp%q(n, n1)
            elements%q(n+(n1-1), 2) = 0
            elements%q(n+(n1-1), 3) = 0

            elements%He(n+(n1-1))=sum(nodes%H(elements%node(n+(n1-1), :)))/3

            elements%c(n+(n1-1)) = 6*Ha/(elements%He(n+(n1-1))**3)
        end do

        call dealloc_elem(etemp)

    end subroutine tr_split

end subroutine adapt
