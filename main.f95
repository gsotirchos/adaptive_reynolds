program main

    use types

    implicit none

    type(fem_node)       :: node
    type(fem_element)    :: element
    real,    allocatable :: P(:)
    integer, allocatable :: splits(:)

    integer :: tries_left

    ! Data input
    call input(node, element)

    allocate(P(size(node%x)))

    ! Calculation
    call calculation(node, element, P)

    allocate(splits(size(element%c)))

    ! Evaluation
    call evaluation(node, element, P, splits)
    
    ! If there is any split necessary then adapt
    tries_left = 20

    do while (any(splits > 0))
        tries_left = tries_left - 1

        ! Adaptation
        print *, "adapting...", tries_left, "tries left"
        call adaptation(node, element, splits)

        ! Re-calculation
        call calculation(node, element, P)

        ! Re-evaluation
        call evaluation(node, element, P, splits)

        ! If no tries are left exit
        if (tries_left == 0) then
            print *, "Aborted at 20 tries"
            exit
        end if
    end do

end program main
