subroutine gen_elements(node, element)

    use types
    use parameters
    use linear_algebra

    implicit none

    type(fem_node),    intent(in)  :: node
    type(fem_element), intent(out) :: element
    real, allocatable              :: angle(:)
    integer, allocatable           :: parent(:, :), neighbor(:, :), peer(:)
    real    :: p(2), p1(2), p2(2), p3(2), Dln(3, 3), D(3, 3)
    integer :: corner(4), cnt


    ! Determine corners of space
    do n = 1, size(node%x)
        if (node%x(n) == 0 .and. node%y(n) == 0) then
            corner(1) = n
        else if (node%x(n) == L .and. node%y(n) == 0) then
            corner(2) = n
        else if (node%x(n) == L .and. node%y(n) == B) then
            corner(3) = n
        else if (node%x(n) == 0 .and. node%y(n) == B) then
            corner(4) = n
        end if
    end do

    ! Generate first two elements
    allocate(element%node(2, 3))

    element%node(1, :) = corner([1, 2, 3])
    element%node(2, :) = corner([1, 3, 4])


    ! Insert nodes one by one
    do n = 1, size(node%x)
        ! Skip corners
        if (any(corner == n)) cycle

        allocate(parent(0, 3))
        parent = 0
        j = 0

        ! Find parents
        do i = size(element%node, dim = 1), 1, -1
            n1 = element%node(i, 1)
            n2 = element%node(i, 2)
            n3 = element%node(i, 3)

            p  = [node%x(n), node%y(n)]
            p1 = [node%x(n1), node%y(n1)]
            p2 = [node%x(n2), node%y(n2)]
            p3 = [node%x(n3), node%y(n3)]

            ! Check if point falls in triangle
            if (PointInTriangle(p, p1, p2, p3)) then
                ! Add parent
                j = j + 1
                call extend(parent, j - 1, 1)
                parent(j, :) = element%node(i, :)

                ! Remove the element
                call extend(element%node, i, -1)
            end if
        end do

        allocate(neighbor(0, 3))
        neighbor = 0
        j = 0

        ! Find neighbors
        do i = size(element%node, dim = 1), 1, -1
            do k = 1, size(parent, dim = 1)

                ! Check how many common nodes between i-th element
                ! and k-th parent
                cnt = 0
                do n1 = 1, 3
                    if (any(element%node(i, :) == parent(k, n1))) then
                        cnt = cnt + 1
                    end if
                end do

                if (cnt >= 2) then
                    n1 = element%node(i, 1)
                    n2 = element%node(i, 2)
                    n3 = element%node(i, 3)

                    p  = [node%x(n), node%y(n)]
                    p1 = [node%x(n1), node%y(n1)]
                    p2 = [node%x(n2), node%y(n2)]
                    p3 = [node%x(n3), node%y(n3)]

                    ! Construct matrix for Delaunay criterion
                    Dln(1, :) = [p1(1)-p(1), p1(2)-p(2), (p1(1)-p(1))**2 + (p1(2)-p(2))**2]
                    Dln(2, :) = [p2(1)-p(1), p2(2)-p(2), (p2(1)-p(1))**2 + (p2(2)-p(2))**2]
                    Dln(3, :) = [p3(1)-p(1), p3(2)-p(2), (p3(1)-p(1))**2 + (p3(2)-p(2))**2]

                    ! Check for Delaunay criterion
                    if (det(Dln) > 0) then
                        ! Add neighbor
                        j = j + 1
                        call extend(neighbor, j - 1, 1)
                        neighbor(j, :) = element%node(i, :)

                        ! Remove the element
                        call extend(element%node, i, -1)
                    
                        ! Stop parent check and move to the next element
                        exit
                    end if
                end if

            end do
        end do


        allocate(peer(0))
        peer = 0
        j = 0

        ! Store parents' nodes
        do i = 1, size(parent, dim = 1)
            do n1 = 1, 3
                if (any(peer == parent(i, n1))) cycle

                ! Add the node
                j = j + 1
                call extend(peer, j - 1, 1)
                peer(j) = (parent(i, n1))
            end do
        end do

        ! Store neighbors' nodes
        do i = 1, size(neighbor, dim = 1)
            do n1 = 1, 3
                if (any(peer == neighbor(i, n1))) cycle

                ! Add the node
                j = j + 1
                call extend(peer, j - 1, 1)
                peer(j) = (neighbor(i, n1))
            end do
        end do


        ! Calculate angle of every peer
        allocate(angle(size(peer)))
        angle = mod(atan(node%y(peer)-node%y(n), node%x(peer)-node%x(n)) &
                         + 2*pi, 2*pi)

        ! Sort peer nodes clockwise
        peer = peer(order(angle))


        ! Add new elements
        do i = 1, size(peer)
            n1 = n
            n2 = peer(i)
            n3 = peer(mod(i + size(peer), size(peer)) + 1)

            D(1, :) = [node%x(n1), node%y(n1), 1.0]
            D(2, :) = [node%x(n2), node%y(n2), 1.0]
            D(3, :) = [node%x(n3), node%y(n3), 1.0]

            if (det(D) <= 1.0E-05*dL*dB/2 .or. isnan(det(D))) cycle

            j = ubound(element%node, dim = 1) + 1
            call extend(element%node, j - 1, 1)

            element%node(j, 1) = n1
            element%node(j, 2) = n2
            element%node(j, 3) = n3
        end do


        deallocate(parent)
        deallocate(neighbor)
        deallocate(peer)
        deallocate(angle)

        !! DEBUG
        !print *, "ELEMENTS:"
        !do i = 1, size(element%node, 1)
        !    print *, element%node(i, :)
        !end do
        !print *, "==================================="
        !print *,

    end do
    
    ! Allocate the rest of the elements' properties matrices
    n = size(element%node, dim = 1)
    allocate(element%q(n, 3))
    allocate(element%He(n))
    allocate(element%c(n))


    !! DEBUG
    !print *, "** DEBUG **"
    !print *,
    !print *, "ELEMENTS"
    !print "(A5, 3A4, 5A6)", " ", "n1", "n2", "n3", &
    !                        "q12", "q23", "q31", &
    !                        "He", "c"
    !do i = 1, size(element%c)
    !    print "(I4, A1, 3I4, 5F6.2)", i, ": ", &
    !        element%node(i, 1), element%node(i, 2), element%node(i, 3), &
    !        element%q(i, 1), element%q(i, 2), element%q(i, 3), &
    !        element%He(i), element%c(i)
    !end do
    !print *,

end subroutine gen_elements
