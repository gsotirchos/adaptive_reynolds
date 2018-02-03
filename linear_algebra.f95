module linear_algebra

contains

    ! Function to calculate distance between two points p1 and p2
    pure function dist(p1, p2)

        implicit none

        real, intent(in)  :: p1(2), p2(2)
        real :: dist

        dist = sqrt(  (p2(1) - p1(1))**2   &
                    + (p2(2) - p1(2))**2  )

    end function dist


    ! Simple matrix inversion function
    ! The calculation is done using Gaussian Elimination
    pure function inv(AA) result(invA)

        implicit none
    
        ! Declaration of variables
        real, intent(in) :: AA(:, :)
        real :: pivot, xnum
    
        real, allocatable :: A(:, :), B(:, :), invA(:, :)
    
        integer :: n, i, j, k

    
        ! Allocation of matrices
        n = size(AA, dim = 1)
        allocate(A(n, n))
        allocate(B(n, 2*n))
        allocate(invA(n, n))
    

        ! Initialization
        A = AA
        invA = 0.0
        

        ! Generate augmented matrix
        do j = 1, n
            do i = 1, n
                B(i, j) = 0.0
                B(i, j + n) = 0.0

                B(i, j) = A(i, j)
                if (i == j) B(i, j + n) = 1.0
            end do
        end do


        do i = 1, n
            ! Chose the leftmost non-zero element as pivot
            do j = 1, n
                if (abs(B(i, j)) > 0) then
                    pivot = B(i, j)
                    exit
                end if
            end do

            ! Step 1: change the chosen pivot into 1 by dividing
            ! the pivot's row by the pivot number
            B(i, :) = B(i, :)/pivot

            ! update pivot value
            pivot = B(i, i)


            ! Step 2: change the remainder of the pivot's column into
            ! 0s by adding to each row a suitable multiple of the pivot row
            do k = 1, n
                if (k /= i) then
                    ! same column with the curent pivot
                    xnum = B(k, i)/pivot
                    col: do j = 1, 2*n
                        B(k, j) = B(k, j) - xnum*B(i, j)
                    end do col
                end if
            end do
        end do

            ! Prepare the final inverted matrix
            do j = 1, n
                do i = 1, n
                    invA(i, j) = B(i, j + n)
                end do
            end do

    end function inv

    ! Function which calculates sign of dot-product
    ! between vectors (p1, p) and (p1, p2)
    real function xProdSign(p, p1, p2)
            
        implicit none

        real :: p(2), p1(2), p2(2)

        xProdSign =   (p(1) - p2(1))*(p1(2) - p2(2)) &
                    - (p1(1) - p2(1))*(p(2) - p2(2))

    end function xProdSign


    ! Fuction which tells if a point p is located inside
    ! a triangle (p1, p2, p3)
    logical function PointInTriangle(p, p1, p2, p3)
      
        implicit none

        real    :: p(2), p1(2), p2(2), p3(2)
        logical :: b1, b2, b3

        b1 = xProdSign(p, p1, p2) < 0
        b2 = xProdSign(p, p2, p3) < 0
        b3 = xProdSign(p, p3, p1) < 0

        PointInTriangle = (b1 .eqv. b2) .and. (b2 .eqv. b3)

    end function PointInTriangle


    ! Function which calculates the upper triangular of matrix  A
    pure function upper(A)

        implicit none
        
        real, intent(in)  :: A(:, :)
        real, allocatable :: upper(:, :)
        real              :: E
        integer           :: n, m, i, j

        n = size(A, dim = 1)
        m = size(A, dim = 2)
        allocate(upper(n, m))

        upper = A

        do j = 1, m - 1
            do i = j + 1, n
                E = upper(i, j)/upper(j, j)
                upper(i, j:n) = upper(i, j:n) - upper(j, j:n)*E
            end do
        end do

    end function upper


    ! Function which calculates the determinant of square matrix A
    pure function det(A) result(detA)

        implicit none

        real, intent(in)  :: A(:, :)
        real, allocatable :: upA(:, :)
        real              :: detA
        integer           :: n, i

        n = size(A, dim = 1)
        allocate(upA(n, n))

        upA = upper(A)

        detA = 1

        do i = 1, n
            detA = detA*upA(i, i)
        end do

    end function det

end module linear_algebra
