program main

    use types

    implicit none

    type(fem_node)       :: node
    type(fem_element)    :: element
    real,    allocatable :: P(:)
    integer, allocatable :: splits(:)

    integer :: tries_left, i

    ! Data input
    call input(node, element)

    allocate(P(size(node%x)))

    ! Calculation
    call calculate(node, element, P)

    allocate(splits(size(element%c)))

    ! Evaluation
    call evaluate(node, element, P, splits)
    
    ! If there is any split necessary then adapt
    tries_left = 15

    !do while (any(splits > 0))
        tries_left = tries_left - 1

        ! Adaptation
        print *, "adapting...", tries_left, "tries left"
        call adapt(node, element, splits)

        ! Re-calculation
        call calculate(node, element, P)

        ! Re-evaluation
        call evaluate(node, element, P, splits)

        ! If no tries are left exit
        if (tries_left == 0) then
            print *, "Aborted, ran out of tries"
    !        exit
        end if
    !end do

end program main
