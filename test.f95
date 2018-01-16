program test

    use types

    implicit none

    type(fem_node) :: node
    type(fem_element) :: element


    call tr_split(node, element, 2)

contains

    subroutine tr_split(nodes, elements, n)

        implicit none

        integer, intent(in) :: n

        type(fem_node)    :: nodes, ntemp
        type(fem_element) :: elements, etemp
        real              :: la, lb, lc
        integer           :: n1, n2, n3


        do n1 = 1, 3
            n2 = mod(n1 + 3, 3) + 1
            n3 = mod(n2 + 3, 3) + 1
            
            ! Find the longest side
            la = dist([nodes%x(elements%node(n, n1)),  &
                       nodes%y(elements%node(n, n1))], &
                      [nodes%x(elements%node(n, n2)),  &
                       nodes%y(elements%node(n, n2))])

            lb = dist([nodes%x(elements%node(n, n2)),  &
                       nodes%y(elements%node(n, n2))], &
                      [nodes%x(elements%node(n, n3)),  &
                       nodes%y(elements%node(n, n3))])

            lc = dist([nodes%x(elements%node(n, n3)),  &
                       nodes%y(elements%node(n, n3))], &
                      [nodes%x(elements%node(n, n1)),  &
                       nodes%y(elements%node(n, n1))])


            if (la > lb .and. la > lc) then
 
                ! Generate new node
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
                nodes%x(ubound(nodes%x)) = (nodes%x(elements%node(n, n1))  &
                                          + nodes%x(elements%node(n, n2)))/2

                nodes%y(ubound(nodes%y)) = (nodes%y(elements%node(n, n1))  &
                                          + nodes%y(elements%node(n, n2)))/2

                nodes%stat(ubound(nodes%stat)) = &
                                          nodes%stat(elements%node(n, n1)) &
                                         *nodes%stat(elements%node(n, n2))

                nodes%P(ubound(nodes%P)) = (nodes%P(elements%node(n, n1))  &
                                          + nodes%P(elements%node(n, n2)))/2

                nodes%H(ubound(nodes%H)) = (nodes%H(elements%node(n, n1))  &
                                          + nodes%H(elements%node(n, n2)))/2

                call dealloc_node(ntemp)


                ! Generate new elements
                call alloc_elem(etemp, size(elements%c), 3)

                etemp%node(:size(elements%c), :) = elements%node
                etemp%q(:size(elements%c), :)    = elements%q
                etemp%He(:size(elements%c))      = elements%He
                etemp%c(:size(elements%c))       = elements%c

                call dealloc_elem(elements)
                call alloc_elem(elements, size(etemp%c) + 1, 3)

                elements%node(1:(n - 1), :) = etemp%node(1:(n - 1), :)
                elements%q(1:(n - 1), :) = etemp%q(1:(n - 1), :)
                elements%He(1:(n - 1)) = etemp%He(1:(n - 1))
                elements%c(1:(n - 1)) = etemp%c(1:(n - 1))

                elements%node((n + 2):, :) = &
                          etemp%node((n + 1):, :)

                elements%q((n + 2):, :) = &
                          etemp%q((n + 1):, :)

                elements%He((n + 2):) = &
                          etemp%He((n + 1):)

                elements%c((n + 2):) = &
                          etemp%c((n + 1):)

                ! Generate new elements' attributes
                elements%node(n, 1)     = etemp%node(n, n3)
                elements%node(n, 2)     = etemp%node(n, n1)
                elements%node(n, 3)     = ubound(nodes%x, dim = 1)
                elements%node(n + 1, 1) = etemp%node(n, n3)
                elements%node(n + 1, 2) = ubound(nodes%x, dim = 1)
                elements%node(n + 1, 3) = etemp%node(n, n2)

                elements%q(n, 1)     = etemp%q(n, n3)
                elements%q(n, 2)     = etemp%q(n, n1)
                elements%q(n, 3)     = 0
                elements%q(n + 1, 1) = 0
                elements%q(n + 1, 2) = etemp%q(n, n1)
                elements%q(n + 1, 3) = etemp%q(n, n2)

                elements%He(n) = etemp%He(n)
                elements%He(n + 1) = etemp%He(n)

                elements%c(n) = etemp%c(n)
                elements%c(n + 1) = etemp%c(n)

                call dealloc_elem(etemp)

                exit

            end if

        end do

    end subroutine tr_split

end program test
