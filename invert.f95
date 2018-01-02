module invert

contains

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

end module invert
