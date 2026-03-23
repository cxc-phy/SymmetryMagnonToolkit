!---------------------------------------------------------------------------------------------------------
module Tool

contains
    subroutine get_inverse ( matrix )
        implicit none
        complex ( 8 ) , intent ( inout ) :: matrix ( : , : )
        ! --- variables for general matrix inverse ---
        integer ( 4 ) :: N , M , LDA , INFO , LWORK
        integer ( 4 ) , allocatable :: IPIV ( : )
        complex ( 8 ) , allocatable :: WORK ( : )

        M = size ( matrix , dim = 1 )
        N = size ( matrix , dim = 2 )

        LDA = max ( 1 , N )
        LWORK = max ( 1 , N )
        allocate ( IPIV ( N ) , WORK ( LWORK ) )
        call ZGETRF ( M , N , matrix , LDA , IPIV , INFO )
        if ( INFO /= 0 ) then
            write ( * , * ) 'ZGETRF went wrong! aborting ...'
            stop
        end if

        call ZGETRI ( N , matrix , LDA , IPIV , WORK , LWORK , INFO )
        if ( INFO /= 0 ) then
            write ( * , * ) 'ZGETRI went wrong! aborting ...'
            stop
        end if

        deallocate ( IPIV , WORK )
        return
    end subroutine

end module Tool
!---------------------------------------------------------------------------------------------------------
module rotation_matrix
    use var
    use Tool
    complex ( 8 ) , allocatable :: rot_diag ( : , : )

contains
    subroutine get_rot ( RM , n )
        ! get rotation matrix that rotate global z axis to local z axis
        implicit none
        integer ( 4 ) , intent ( in ) :: n
        real ( 8 ) , intent ( inout ) :: RM ( 3 , 3 )
        ! --- local variable ---
        ! r^2 = sx^2 + sy^2 + sz^2
        ! d^2 = sx^2 + sy^2
        ! theta : rotation angle along y axis
        ! phi : rotation angle along z axis
        integer ( 4 ) :: i
        real ( 8 ) :: r , d , theta , phi

        if ( size ( s_coord , dim = 2 ) == 3 ) then
            r = sqrt ( s_coord ( n , 1 ) ** 2 + s_coord ( n , 2 ) ** 2 + s_coord ( n , 3 ) ** 2 )
            d = sqrt ( s_coord ( n , 1 ) ** 2 + s_coord ( n , 2 ) ** 2 )
            theta = acos ( s_coord ( n , 3 ) / r )
            if ( d > 1.d-5 ) then
                phi = acos ( s_coord ( n , 1 ) / d )
            else
                phi = 0.d0
            end if
            if ( s_coord ( n , 2 ) < 0.d0 ) then
                phi = 2.d0 * pi - phi
            end if
        else
            theta = s_coord ( n , 1 ) / dble ( 360 ) * 2.d0 * pi
            phi   = s_coord ( n , 2 ) / dble ( 360 ) * 2.d0 * pi
        end if

        RM = 0.d0
        RM ( 1 , 1 ) = cos ( phi ) * cos ( theta )
        RM ( 1 , 2 ) = sin ( phi ) * ( -1.d0 )
        RM ( 1 , 3 ) = cos ( phi ) * sin ( theta )
        RM ( 2 , 1 ) = sin ( phi ) * cos ( theta )
        RM ( 2 , 2 ) = cos ( phi )
        RM ( 2 , 3 ) = sin ( phi ) * sin ( theta )
        RM ( 3 , 1 ) = sin ( theta ) * ( -1.d0 )
        RM ( 3 , 3 ) = cos ( theta )
        return
    end subroutine get_rot

    subroutine get_rot_tmul ( RM , n , neq , i , j )
        implicit none
        integer ( 4 ) , intent ( in ) :: n , neq , i , j
        real ( 8 ) , intent ( inout ) :: RM ( 3 , 3 )
        ! --- local variable ---
        integer ( 4 ) :: k , l , o
        real ( 8 ) , dimension ( 3 , 3 ) :: RMi , RMj
        real ( 8 ) :: JJ ( size ( Jvalue , dim = 3 ) )

        call get_rot ( RMi , i )
        call get_rot ( RMj , j )

        JJ = Jvalue ( n , neq , : )

        RM = 0.d0
        do k = 1 , 3
            do l = 1 , 3
                do o = 1 , 3
                    RM ( k , l ) = RM ( k , l ) + JJ ( o ) * RMi ( o , k ) * RMj ( o , l )
                end do
            end do
        end do

        return
    end subroutine get_rot_tmul

    subroutine get_rot_coeff_Jex ( rot_coeff_Jex , n , neq , i , j , label )
        implicit none
        integer ( 4 ) , intent ( in ) :: n , neq , i , j , label
        complex ( 8 ) , intent ( inout ) :: rot_coeff_Jex
        ! --- local variable ---
        real ( 8 ) :: RM ( 3 , 3 )

        rot_coeff_Jex = v0
        call get_rot_tmul ( RM , n , neq , i , j )
        if ( label == 1 ) then
            rot_coeff_Jex = v1 * RM ( 1 , 1 ) &
                         &- vi * RM ( 1 , 2 ) &
                         &+ vi * RM ( 2 , 1 ) &
                         &+ v1 * RM ( 2 , 2 )
        else if ( label == 2 ) then
            rot_coeff_Jex = v1 * RM ( 1 , 1 ) &
                         &+ vi * RM ( 1 , 2 ) &
                         &+ vi * RM ( 2 , 1 ) &
                         &- v1 * RM ( 2 , 2 )
        else if ( label == 3 ) then
            rot_coeff_Jex = v1 * RM ( 1 , 1 ) &
                         &- vi * RM ( 1 , 2 ) &
                         &- vi * RM ( 2 , 1 ) &
                         &- v1 * RM ( 2 , 2 )
        else if ( label == 4 ) then
            rot_coeff_Jex = v1 * RM ( 1 , 1 ) &
                         &+ vi * RM ( 1 , 2 ) &
                         &- vi * RM ( 2 , 1 ) &
                         &+ v1 * RM ( 2 , 2 )
        else if ( label == 5 ) then
            rot_coeff_Jex = v1 * RM ( 3 , 3 )
        else
            write ( * , * ) 'error occurs at get_rot_coeff_Jex'
        end if

        return
    end subroutine get_rot_coeff_Jex

    subroutine get_rot_coeff_DMI ( rot_coeff_DMI , i , j , label )
        implicit none
        integer ( 4 ) , intent ( in ) :: i , j , label
        complex ( 8 ) , intent ( inout ) :: rot_coeff_DMI ( 3 , 3 )
        ! --- local variable ---
        integer ( 4 ) :: u , v
        real ( 8 ) :: RMi ( 3 , 3 ) , RMj ( 3 , 3 )

        call get_rot ( RMi , i )
        call get_rot ( RMj , j )
        if ( label == 1 ) then
            do u =  1 , 3 , 1
                do v = 1 , 3 , 1
                    rot_coeff_DMI ( u , v ) = v1 * RMi ( u , 1 ) * RMj ( v , 1 ) &
                                           &- vi * RMi ( u , 1 ) * RMj ( v , 2 ) &
                                           &+ vi * RMi ( u , 2 ) * RMj ( v , 1 ) &
                                           &+ v1 * RMi ( u , 2 ) * RMj ( v , 2 )
                end do
            end do
        else if ( label == 2 ) then
            do u =  1 , 3 , 1
                do v = 1 , 3 , 1
                    rot_coeff_DMI ( u , v ) = v1 * RMi ( u , 1 ) * RMj ( v , 1 ) &
                                           &+ vi * RMi ( u , 1 ) * RMj ( v , 2 ) &
                                           &+ vi * RMi ( u , 2 ) * RMj ( v , 1 ) &
                                           &- v1 * RMi ( u , 2 ) * RMj ( v , 2 )
                end do
            end do
        else if ( label == 3 ) then
            do u =  1 , 3 , 1
                do v = 1 , 3 , 1
                    rot_coeff_DMI ( u , v ) = v1 * RMi ( u , 1 ) * RMj ( v , 1 ) &
                                           &- vi * RMi ( u , 1 ) * RMj ( v , 2 ) &
                                           &- vi * RMi ( u , 2 ) * RMj ( v , 1 ) &
                                           &- v1 * RMi ( u , 2 ) * RMj ( v , 2 )
                end do
            end do
        else if ( label == 4 ) then
            do u =  1 , 3 , 1
                do v = 1 , 3 , 1
                    rot_coeff_DMI ( u , v ) = v1 * RMi ( u , 1 ) * RMj ( v , 1 ) &
                                           &+ vi * RMi ( u , 1 ) * RMj ( v , 2 ) &
                                           &- vi * RMi ( u , 2 ) * RMj ( v , 1 ) &
                                           &+ v1 * RMi ( u , 2 ) * RMj ( v , 2 )
                end do
            end do
        else if ( label == 5 ) then
            do u =  1 , 3 , 1
                do v = 1 , 3 , 1
                    rot_coeff_DMI ( u , v ) = v1 * RMi ( u , 3 ) * RMj ( v , 3 )
                end do
            end do
        else
            write ( * , * ) 'error occurs at get_rot_coeff_DMI'
        end if

        return
    end subroutine get_rot_coeff_DMI

    subroutine get_rot_diag ( dim , label )
        implicit none
        integer ( 4 ) , intent ( in ) :: dim
        character ( 1 ) , intent ( in ) :: label
        !complex ( 8 ) , allocatable , intent ( inout ) :: rot_diag ( : , : )
        ! --- local variable ---
        integer ( 4 ) :: i
        real ( 8 ) :: RM ( 3 , 3 )

        if ( FMConfig .eqv. .False. ) then
            allocate ( rot_diag ( 2 * n_ineq , 2 * n_ineq ) )
        else
            allocate ( rot_diag ( n_ineq , n_ineq ) )
        end if

        rot_diag = v0
        if ( label == '+' ) then
            do i = 1 , dim , 1
                call get_rot ( RM , i )
                rot_diag ( i , i ) = v1 * RM ( 1 , 1 ) &
                                  &- vi * RM ( 1 , 2 ) &
                                  &+ vi * RM ( 2 , 1 ) &
                                  &+ v1 * RM ( 2 , 2 )
                if ( rot_diag ( i , i ) == 0.d0 ) rot_diag ( i , i ) = 1.d-14
                if ( FMConfig .eqv. .False. ) rot_diag ( dim + i , dim + i ) = cmplx ( rot_diag ( i , i ) )
            end do
            call get_inverse ( rot_diag )
        else if ( label == '-' ) then
            do i = 1 , dim , 1
                call get_rot ( RM , i )
                rot_diag ( i , i ) = v1 * RM ( 1 , 1 ) &
                                  &- vi * RM ( 1 , 2 ) &
                                  &- vi * RM ( 2 , 1 ) &
                                  &- v1 * RM ( 2 , 2 )
                if ( rot_diag ( i , i ) == 0.d0 ) rot_diag ( i , i ) = 1.d-14
                if ( FMConfig .eqv. .False. ) rot_diag ( dim + i , dim + i ) = cmplx ( rot_diag ( i , i ) )
            end do
            call get_inverse ( rot_diag )
        else
            write ( * , * ) 'error occurs at get_rot_diag : wrong label'
            stop
        end if

        return
    end subroutine get_rot_diag

end module rotation_matrix
!--------------------------------------------------------------------------------------------------------------------
module interaction_coeff
    use var
    use rotation_matrix

contains
    ! --- get exchange coefficient J_k \times S_i \cdot S_j ---
    subroutine get_Jk ( Jk , kp , layer_mark , label )
        implicit none
        real ( 8 ) , intent ( in ) :: kp ( : )
        integer ( 4 ) , intent ( in ) :: label
        complex ( 8 ) , intent ( inout ) :: Jk ( : , : )
        character ( 2 ) , intent ( in ) :: layer_mark
        ! --- local variables ---
        ! a : lattice vector
        ! as : lattice vector on stack direction
        ! as = 0 : only consider contribution from neighbors in same layers
        ! as = 1/-1 : only consider contribution from neighbors in neighbor layers in two directions
        real ( 8 ) , dimension ( size ( kp ) ) :: a
        integer ( 4 ) :: i , j , k
        integer ( 4 ) :: count , n , neq , source , goal
        integer ( 4 ) :: as
        complex ( 8 ) :: rot_coeff_Jex

        Jk = v0
        if ( layer_mark == 'bk' ) then
            ! --- achieve bulk hamiltonian ---
            do count = 1 , size ( neighbor_list , dim = 1 ) , 1
                n = neighbor_list ( count , 1 )
                neq = neighbor_list ( count , 2 )
                source = neighbor_list ( count , 3 )
                goal   = neighbor_list ( count , 4 )
                a = neighbor_list ( count , 5 : 4 + size ( kp ) )
                call get_rot_coeff_Jex ( rot_coeff_Jex , n , neq , source , goal , label )
                Jk ( source , goal ) = Jk ( source , goal ) &
                                    &+ rot_coeff_Jex * exp( vi * dot_product ( kp , a ) )
            end do
        else
            ! --- achieve intra-layer / inter-layer hamiltonian ---
            if ( ( layer_mark == '0u' ) .or. ( layer_mark == '0b' ) .or. ( layer_mark == '0l' )  .or. ( layer_mark == '00' ) ) then
                as = 0
            else if ( layer_mark == '10' ) then
                as = -1
            else if ( layer_mark == '01' ) then
                as = 1
            else
                write(*,*) 'error occurs while getting Jk'
                stop
            end if

            if ( dir == 1 ) then ! stacking along x direction
                i = 6
                j = 7
                k = 5
            else if ( dir == 2 ) then ! stacking along y direction
                i = 5
                j = 7
                k = 6
            else if ( dir == 3 ) then ! stacking along z direction
                i = 5
                j = 6
                k = 7
            else
                write ( * , * ) 'error occurs while getting Jk'
                stop
            end if
            do count = 1 , size ( neighbor_list , dim = 1 ) , 1
                if ( neighbor_list ( count , k ) /= as ) cycle
                n  = neighbor_list ( count , 1 )
                neq = neighbor_list ( count , 2 )
                source = neighbor_list ( count , 3 )
                goal   = neighbor_list ( count , 4 )
                if ( any ( (/ slab2D , surf2D /) ) ) then
                    a ( 1 ) = neighbor_list ( count , i )
                    a ( 2 ) = neighbor_list ( count , j )
                else if ( any ( (/ slab1D , surf1D /) ) ) then
                    a = neighbor_list ( count , i )
                else
                    write ( * , * ) 'lattice vector not read in while getting Dk'
                    stop
                end if

                call get_rot_coeff_Jex ( rot_coeff_Jex , n , neq , source , goal , label )
                Jk ( source , goal ) = Jk ( source , goal ) &
                                    &+ rot_coeff_Jex * exp ( vi * dot_product ( kp , a ) )
            end do
        end if
        return
    end subroutine get_Jk

    ! --- get DMI coefficient D_k \cdot ( S_i \cross S_j ) ---
    subroutine get_Dk ( Dk , kp , layer_mark )
        implicit none
        real ( 8 ) , intent ( in ) :: kp ( : )
        character ( 2 ) , intent ( in ) :: layer_mark
        complex ( 8 ) , allocatable , intent ( inout ) :: Dk ( : , : , : )
        ! --- local variable ---
        real ( 8 ) :: a ( size ( kp ) )
        integer ( 4 ) :: i , j , k
        integer ( 4 ) :: count , n , neq , source , goal
        integer ( 4 ) :: as

        Dk = v0
        if ( layer_mark == 'bk' ) then
            ! --- achieve bulk hamiltonian ---
            do count = 1 , size ( Dvalue , dim = 1 ) , 1
                source = Dvalue ( count , 1 )
                goal   = Dvalue ( count , 2 )
                a = Dvalue ( count , 3 : 2 + size ( kp ) )
                Dk ( source , goal , : ) = Dk ( source , goal , : ) &
                                        &+ Dvalue ( count , 6 : ) * exp ( vi * dot_product ( kp , a ) )
            end do
        else
            ! --- achieve intra-layer / inter-layer hamiltonian ---
            if ( ( layer_mark == '0u' ) .or. ( layer_mark == '0b' ) .or. ( layer_mark == '0l' ) .or. ( layer_mark == '00' ) ) then
                as = 0
            else if ( layer_mark == '10' ) then
                as = -1
            else if ( layer_mark == '01' ) then
                as = 1
            else
                write(*,*) 'error occurs while getting Dk'
                stop
            end if

            if ( dir == 1 ) then ! stacking along x direction
                i = 4
                j = 5
                k = 3
            else if ( dir == 2 ) then ! stacking along y direction
                i = 3
                j = 5
                k = 4
            else if ( dir == 3 ) then ! stacking along z direction
                i = 3
                j = 4
                k = 5
            else
                write ( * , * ) 'error occurs while getting Dk'
                stop
            end if

            do count = 1 , size ( Dvalue , dim = 1 ) , 1
                if ( Dvalue ( count , k ) /= as ) cycle
                source = Dvalue ( count , 1 )
                goal   = Dvalue ( count , 2 )
                if ( any ( (/ slab2D , surf2D /) ) ) then
                    a ( 1 ) = Dvalue ( count , i )
                    a ( 2 ) = Dvalue ( count , j )
                else if ( any ( (/ slab1D , surf1D /) ) ) then
                    a = Dvalue ( count , i )
                else
                    write ( * , * ) 'lattice vector not read in while getting Dk'
                    stop
                end if

                Dk ( source , goal , : ) = Dk ( source , goal , : ) &
                                        &+ Dvalue ( count , 6 : ) * exp ( vi * dot_product ( kp , a ) )
            end do
        end if

        return
    end subroutine get_Dk

    subroutine get_rot_Dk_cp ( RD , i , j , label , kp , layer_mark )
        implicit none
        complex ( 8 ) , intent ( inout ) :: RD
        integer ( 4 ) , intent ( in ) :: i , j , label
        real    ( 8 ) , intent ( in ) :: kp ( : )
        character ( 2 ) , intent ( in ) :: layer_mark
        ! --- local variable ---
        complex ( 8 ) :: rot_coeff_DMI ( 3 , 3 )
        complex ( 8 ) , allocatable :: Dk ( : , : , : )

        allocate ( Dk ( n_ineq , n_ineq , 3 ) )
        call get_Dk ( Dk , kp , layer_mark )
        call get_rot_coeff_DMI ( rot_coeff_DMI , i , j , label )
        RD = rot_coeff_DMI ( 1 , 2 ) * Dk ( i , j , 3 ) - rot_coeff_DMI ( 2 , 1 ) * Dk ( i , j , 3 ) &
          &+ rot_coeff_DMI ( 3 , 1 ) * Dk ( i , j , 2 ) - rot_coeff_DMI ( 1 , 3 ) * Dk ( i , j , 2 ) &
          &+ rot_coeff_DMI ( 2 , 3 ) * Dk ( i , j , 1 ) - rot_coeff_DMI ( 3 , 2 ) * Dk ( i , j , 1 )

        deallocate ( Dk )

        return
    end subroutine get_rot_Dk_cp

end module interaction_coeff
!--------------------------------------------------------------------------------------------------------------------
module Hamiltonian
    use var
    use interaction_coeff

contains
    subroutine get_H ( Hk , kp , layer_mark )
        implicit none
        real ( 8 ) , intent ( in ) :: kp ( : )
        complex ( 8 ) , intent ( inout ) :: Hk ( : , : )
        character ( 2 ) , intent ( in ) :: layer_mark
        ! --- local variable ---
        integer ( 4 ) :: i , j
        real    ( 8 ) :: k0 ( size ( kp ) )
        complex ( 8 ) :: RD
        complex ( 8 ) , allocatable :: D_0_rot ( : ) , Jk ( : , : )

        ! --- achieve Jex part ---
        allocate ( D_0_rot ( 2 * n_ineq ) , Jk ( n_ineq , n_ineq ) )
        Hk = v0
        call get_Jk ( Jk , kp , layer_mark , 1 )
        Hk ( : n_ineq , : n_ineq ) = 0.5 * S * Jk
        call get_Jk ( Jk , kp , layer_mark , 2 )
        Hk ( : n_ineq , n_ineq + 1 : ) = 0.5 * S * Jk
        call get_Jk ( Jk , kp , layer_mark , 3 )
        Hk ( n_ineq + 1 : , : n_ineq ) = 0.5 * S * Jk
        call get_Jk ( Jk , kp , layer_mark , 4 )
        Hk ( n_ineq + 1 : , n_ineq + 1 : ) = 0.5 * S * Jk

        k0 = 0.0
        if ( layer_mark == 'bk' ) then
            call get_Jk ( Jk , k0 , layer_mark , 5 )
            do i = 1 , n_ineq
                Hk ( i , i ) = Hk ( i , i ) - S * sum ( Jk ( i , : ) )
                Hk ( i + n_ineq , i + n_ineq ) = Hk ( i + n_ineq , i + n_ineq ) - S * sum ( Jk ( : , i ) )
            end do
        else
            if ( any ( layer_mark == (/ '0u' , '0b' , '0l' , '00' /) ) ) then
                call get_Jk ( Jk , k0 , layer_mark , 5 )
                do i = 1 , n_ineq
                    Hk ( i , i ) = Hk ( i , i ) - S * sum ( Jk ( i , : ) )
                    Hk ( i + n_ineq , i + n_ineq ) = Hk ( i + n_ineq , i + n_ineq ) - S * sum ( Jk ( : , i ) )
                end do
            end if
            if ( any ( layer_mark == (/ '0u' , '0b' /) ) ) then
                call get_Jk ( Jk , k0 , '10' , 5 )
                do i = 1 , n_ineq
                    Hk ( i , i ) = Hk ( i , i ) - S * sum ( Jk ( i , : ) )
                end do
                call get_Jk ( Jk , k0 , '01' , 5 )
                do i = 1 , n_ineq
                    Hk ( i + n_ineq , i + n_ineq ) = Hk ( i + n_ineq , i + n_ineq ) - S * sum ( Jk ( : , i ) )
                end do
            end if
            if ( any ( layer_mark == (/ '0b' , '0l' /) ) ) then
                call get_Jk ( Jk , k0 , '01' , 5 )
                do i = 1 , n_ineq
                    Hk ( i , i ) = Hk ( i , i ) - S * sum ( Jk ( i , : ) )
                end do
                call get_Jk ( Jk , k0 , '10' , 5 )
                do i = 1 , n_ineq
                    Hk ( i + n_ineq , i + n_ineq ) = Hk ( i + n_ineq , i + n_ineq ) - S * sum ( Jk ( : , i ) )
                end do
            end if
        end if

        ! --- achieve DMI part ---
        if ( DMI .eqv. .True. ) then ! include DMI term
            D_0_rot = v0
            if ( layer_mark == 'bk' ) then
                do i = 1 , n_ineq , 1
                    do j = 1 , n_ineq , 1
                        call get_rot_Dk_cp ( RD , i , j , 5 , k0 , layer_mark )
                        D_0_rot ( i ) = D_0_rot ( i ) + RD
                        call get_rot_Dk_cp ( RD , j , i , 5 , k0 , layer_mark )
                        D_0_rot ( i + n_ineq ) = D_0_rot ( i + n_ineq ) + RD
                    end do
                end do
            else
                if ( any ( layer_mark == (/ '0u' , '0b' , '0l' , '00' /) ) ) then
                    do i = 1 , n_ineq , 1
                        do j = 1 , n_ineq , 1
                            call get_rot_Dk_cp ( RD , i , j , 5 , k0 , layer_mark )
                            D_0_rot ( i ) = D_0_rot ( i ) + RD
                            call get_rot_Dk_cp ( RD , j , i , 5 , k0 , layer_mark )
                            D_0_rot ( i + n_ineq ) = D_0_rot ( i + n_ineq ) + RD
                        end do
                    end do
                end if
                if ( any ( layer_mark == (/ '0u' , '0b' /) ) ) then
                    do i = 1 , n_ineq , 1
                        do j = 1 , n_ineq , 1
                            call get_rot_Dk_cp ( RD , i , j , 5 , k0 , '10' )
                            D_0_rot ( i ) = D_0_rot ( i ) + RD
                            call get_rot_Dk_cp ( RD , j , i , 5 , k0 , '01' )
                            D_0_rot ( i + n_ineq ) = D_0_rot ( i + n_ineq ) + RD
                        end do
                    end do
                end if
                if ( any ( layer_mark == (/ '0b' , '0l' /) ) ) then
                    do i = 1 , n_ineq , 1
                        do j = 1 , n_ineq , 1
                            call get_rot_Dk_cp ( RD , i , j , 5 , k0 , '01' )
                            D_0_rot ( i ) = D_0_rot ( i ) + RD
                            call get_rot_Dk_cp ( RD , j , i , 5 , k0 , '10' )
                            D_0_rot ( i + n_ineq ) = D_0_rot ( i + n_ineq ) + RD
                        end do
                    end do
                end if
            end if

            do i = 1 , n_ineq , 1
                do j = 1 , n_ineq , 1
                    call get_rot_Dk_cp ( RD , i , j , 1 , kp , layer_mark )
                    Hk ( i , j ) = Hk ( i , j ) + 0.5 * S * RD
                    call get_rot_Dk_cp ( RD , i , j , 2 , kp , layer_mark )
                    Hk ( i , j + n_ineq ) = Hk ( i , j + n_ineq ) + 0.5 * S * RD
                    call get_rot_Dk_cp ( RD , i , j , 3 , kp , layer_mark )
                    Hk ( i + n_ineq , j ) = Hk ( i + n_ineq , j ) + 0.5 * S * RD
                    call get_rot_Dk_cp ( RD , i , j , 4 , kp , layer_mark )
                    Hk ( i + n_ineq , j + n_ineq ) = Hk ( i + n_ineq , j + n_ineq ) + 0.5 * S * RD
                end do
                Hk ( i , i ) = Hk ( i , i ) - S * D_0_rot ( i )
                Hk ( i + n_ineq , i + n_ineq ) = Hk ( i + n_ineq , i + n_ineq ) - S * D_0_rot ( i + n_ineq )
            end do
        end if

        return
        Hk(1,2) = Hk(1,2) + vi * 1d-2 * exp(vi*(kp(2)-kp(1))/2.d0)
        !Hk(1,2) = vi * 1d-8 * exp(vi*(kp(2)-kp(1))/2.d0)
        Hk(2,1) = conjg(Hk(1,2))
        Hk(4,5) = conjg(Hk(1,2))
        Hk(5,4) = Hk(1,2)

        Hk(7,8) = Hk(1,2)
        Hk(8,7) = conjg(Hk(1,2))
        Hk(10,11) = conjg(Hk(1,2))
        Hk(11,10) = Hk(1,2)

        if (abs(Hk(1,2)*exp(vi*(kp(1)-kp(2))/2.d0)-conjg(Hk(1,2)*exp(vi*(kp(1)-kp(2))/2.d0))) > 1d-8 ) then
            write(*,*) "not conjugate"
        end if
        if (abs(Hk(1,3)*exp(vi*(kp(1)-kp(3))/2.d0)-conjg(Hk(1,3)*exp(vi*(kp(1)-kp(3))/2.d0))) > 1d-8 ) then
            write(*,*) "oops,not conjugate"
            stop
        end if
        ! --- achieve single ion anisotropy part ---
        if ( SIA .eqv. .True. ) then
            do i = 1 , 2 * n_ineq
                Hk ( i , i ) = Hk ( i , i ) - cmplx ( 2 , 0 ) * S * Avalue
            end do
        end if

        return
    end subroutine get_H

end module Hamiltonian
!--------------------------------------------------------------------------------------------------------------------
module Diagonalization
    use Tool
    use var
    use Hamiltonian

contains
    ! --- for diagonalizing positive semi-definite BdG Hamiltonian ---
    ! colpa 1978
    ! BdG = KM^+ * KM
    ! VR2^{-1} * [ KM * tz * KM^+ ] * VR2 = tz * EV1
    ! VR1 = KM^{-1} * VR2 * sqrt ( EV1 )
    ! VR1^+ * BdG * VR1 = EV1
    subroutine get_diag ( matrix , EV , VR )
        implicit none
        real    ( 8 ) , intent ( inout ) :: EV ( : )
        complex ( 8 ) , intent ( inout ) :: VR ( : , : )
        complex ( 8 ) , intent ( in ) :: matrix ( : , : )
        ! --- local variable ---
        integer ( 4 ) :: i , j , dim
        logical :: label
        complex ( 8 ) , allocatable :: VR1 ( : , : ) , VR2 ( : , : )
        ! --- variables for cholesky decomposition and
        ! for complex hermitian matrix diagonalization ---
        character :: JOBZ = 'V' , UPLO = 'U'
        integer ( 4 ) :: INFO , LWORK , LDA
        real    ( 8 ) , dimension ( : ) , allocatable :: RWORK , EV1
        complex ( 8 ) , dimension ( : ) , allocatable :: WORK
        complex ( 8 ) , dimension ( : , : ) , allocatable :: INPUT , KM

        dim = size ( matrix , dim = 1 )
        LDA = max ( 1 , dim )
        LWORK = max ( 1 , 2 * dim - 1 )
        allocate ( INPUT ( LDA , dim ) )
        allocate ( RWORK ( 3 * dim - 2 ) , WORK ( LWORK ) , EV1 ( dim ) )
        allocate ( VR1  ( dim , dim ) , VR2 ( dim , dim ) )

        ! check hermicity of hamiltonian
        INPUT = matrix - transpose ( conjg ( matrix ) )
        if ( any ( abs ( INPUT ) > 1d-8 ) ) then
            write ( * , * ) 'hamiltonian not hermitian, please check'
            write ( * , * ) INPUT
            stop
        end if

        ! decompose BdG into KM^+ * KM
        if ( FMConfig .eqv. .False. ) then
            allocate ( KM ( LDA , dim ) )
            label = .True.
            ! --- perform cholesky decomposition ---
            ! only suitable for postive definite case
            if ( label .eqv. .False. ) then
                INPUT = matrix
                call ZPOTRF (UPLO,dim,INPUT,LDA,INFO)
                if ( INFO /= 0 ) then
                    write ( * , * ) 'Cholesky decomposition went wrong! aborting ...'
                    stop
                end if

                KM = v0
                do i = 1 , dim
                    do j = i , dim
                        KM ( i , j ) = INPUT ( i , j )
                    end do
                end do
            ! --- general decomposition ---
            ! suitable for positive semi-definite case
            else if ( label .eqv. .True. ) then
                INPUT = matrix
                call ZHEEV (JOBZ,UPLO,dim,INPUT,LDA,EV1,WORK,LWORK,RWORK,INFO)
                if ( INFO /= 0 ) then
                    write ( * , * ) 'Similarity diagonalization of BdG went wrong! aborting ...'
                    write ( * , * ) 'INFO =' , INFO
                    stop
                end if
                if ( any ( EV1 < -1d-5 ) ) then
                    if ( mode == 0 ) then
                        err = -1
                    else
                        write ( * ,* ) 'oops, BdG hamiltonian not positive semi-definite' , EV1
                        !stop
                    end if
                end if
                do i = 1 , dim
                    KM ( : , i ) = sqrt ( abs ( EV1 ( i ) ) ) * INPUT ( : , i )
                end do
                KM = matmul ( KM , transpose ( conjg ( INPUT ) ) )
            end if
        end if

        ! diagonalization
        if ( FMConfig .eqv. .False. ) then
            ! diagonalize complex hermitian matrix KM * tz * KM^+
            ! EV1 is also eigenvalues of dynamical matrix with form :
            ! { w1 , ... , wn , -wn , ... , -w1 }

            INPUT = matmul ( KM , matmul ( tz , transpose ( conjg ( KM ) ) ) )
        else if ( FMConfig .eqv. .True. ) then
            ! diagonalize hamiltonian matrix
            INPUT = matrix
        end if

        call ZHEEV (JOBZ,UPLO,dim,INPUT,LDA,EV1,WORK,LWORK,RWORK,INFO)
        if ( INFO /= 0 ) then
            write ( * , * ) 'Diagonalization went wrong! aborting ...'
            write ( * , * ) 'INFO =' , INFO
            stop
        end if
        EV = EV1! ( dim / 2 + 1 : )
        !if ( mode == 1 ) return

        ! get transformation matrix in descending order of eigenvalues
        ! for FM case : VR1 ^ {-1} * H * VR1 = E
        ! for other cases : VR1 ^ {+} * H * VR1 = E ; VR1^{-1} * BdG * VR1 = tz * E
        VR2 = v0
        do i = 1 , dim
            VR2 ( : , i ) = INPUT ( : , dim + 1 - i )
        end do

        if ( FMConfig .eqv. .False. ) then
            ! get eigenvalues of BdG Hamiltonian with form :
            ! { w1 , ... , wn , wn , ... , w1 }
            INPUT = v0
            do i = 1 , dim
                INPUT ( i , i ) = tz ( i , i ) * EV1 ( dim + 1 - i )
            end do
            ! get transformation matrix
            VR1 = v0
            if ( FMConfig .eqv. .False. ) then
                call get_inverse ( KM )
                VR1 = matmul( KM , matmul ( VR2 , sqrt ( INPUT ) ) )
            end if
        else if ( FMConfig .eqv. .True. ) then
            VR1 = VR2
        end if

        if ( mode /= 4 ) return
        ! save bandn eigenvectors
        VR ( : , 1 : bandn ) = VR1 ( : , dim / 2 + bandi : dim / 2 + bandf )
        !if ( mode == 4 ) VR ( : , 1 : bandn ) = VR1 ( : , dim / 2 - bandn + 1 : dim / 2 )
        !if ( mode == 4 ) VR ( : , 1 : bandn ) = VR1 ( : , dim - bandn + 1 : dim )

        ! verify property of transformation matrix VR :
        ! T^{+} * tz * T = T * tz * T^{+} = tz
        VR2 = matmul ( matmul ( VR1 , tz ) , transpose ( conjg ( VR1 ) ) ) - tz
        VR1 = matmul ( transpose ( conjg ( VR1 ) ) , matmul ( tz , VR1 ) ) - tz

        if ( any ( abs ( VR1 ) > 1d-8 ) .or. any ( abs ( VR2 ) > 1d-8 ) ) then
            write ( * , * ) 'error : commutation relation not kept'
            stop
        end if
        return
    end subroutine get_diag

end module Diagonalization
