module parameters

    integer, parameter :: xsub = 3, ysub = 3, &
                          num_nodes = (xsub + 1)*(ysub + 1), &
                          ! Triangular elements
                          elem_nodes = 3, num_elem = 2*xsub*ysub, &
                          ! Only dof is node pressure
                          node_dof = 1, tot_dof = node_dof*num_nodes

    real, parameter    :: L = 10.0, B = 10.0, &
                          dL = L/xsub, dB = B/xsub, &
                          ! Reynolds equation parameters
                          hx = 1.0, hy = 1.2, &
                          P_bound = 0.0, q = 0.1, &
                          ! H = Ha*x + Hb
                          Ha = -0.1, Hb = 2.0

    integer            :: i, j, n, ni, nj, n1, n2, n3

end module parameters
