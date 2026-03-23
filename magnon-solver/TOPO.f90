!---------------------------------------------------------------------------------------------------------
module TopoProp
    use Tool
    use var
    use Hamiltonian
    use Diagonalization
contains
    ! --- obtain km in effective BZ ---
    ! km is the kpoint mesh inside effective bz
    ! Re ( km ( i , j ) ) = kx_i , Im ( km ( i , j ) ) = ky_j
    ! if time reversal symmetric
    !   consider half BZ ,  kx : -0.5 - 0.5 , ky : -0.5 - 0 ( in unit of 2pi )
    !   kmesh has no periodicity on ky
    ! if not
    !   consider whole BZ , kx : -0.5 - 0.5 , ky : -0.5 - 0.5 ( in unit of 2pi )
    subroutine get_2dBZ_discrete ( km , n1 , n2 )
        implicit none
        integer ( 4 ) , intent ( inout ) :: n1
        integer ( 4 ) , intent ( in ) :: n2
        real    ( 8 ) , allocatable , intent ( inout ) :: km ( : , : , : )
        ! --- local variables ---
        integer ( 4 ) :: i , j , n1m
        real    ( 8 ) :: dk1 , dk2

        if ( trs .eqv. .True. ) then
            if ( mod ( n1 , 2 ) == 0 ) n1 = n1 + 1
            dk1 = dble ( 1 ) / dble ( n1 - 1 )
            dk2 = dble ( 1 ) / dble ( 2 ) / dble ( n2 - 1 )
            n1m = ( n1 + 1 ) / 2
        else
            dk1 = dble ( 1 ) / dble ( n1 )
            dk2 = dble ( 1 ) / dble ( n2 )
        end if
        allocate ( km ( n1 , n2 , 2 ) )

        km = -1.0
        km (  1  , 1  , : ) = (/ -0.5 , -0.5 /)
        km (  1  , n2 , : ) = (/ -0.5 ,  0.0 /)
        km ( n1m , 1  , : ) = (/  0.0 , -0.5 /)
        km ( n1m , n2 , : ) = (/  0.0 ,  0.0 /)

        do i = 1 , n1m , 1
            do j = 1 , n2 , 1
                if ( any ( km ( i , j , : ) /= -1.0 )  ) cycle
                if ( ( j == 1 ) .or. ( j == n2 ) ) then
                    km ( i , j , 1 ) = km ( i - 1 , j , 1 ) + dk1
                    km ( i , j , 2 ) = km ( i - 1 , j , 2 )
                else
                    km ( i , j , 1 ) = km ( i , j - 1 , 1 )
                    km ( i , j , 2 ) = km ( i , j - 1 , 2 ) + dk2
                end if
            end do
        end do

        do i = n1m + 1 , n1 , 1
            do j = 1 , n2 , 1
                if ( i == n1 ) then
                    km ( i , j , : ) = km ( 2 * n1m - i , j , : )
                else
                    km ( i , j , 1 ) = -km ( 2 * n1m - i , j , 1 )
                    km ( i , j , 2 ) =  km ( 2 * n1m - i , j , 2 )
                end if
            end do
        end do

        return
    end subroutine get_2dBZ_discrete

    subroutine get_U1 ( D1 , D2 , u1 )
        implicit none
        complex ( 8 ) , intent ( inout ) :: u1
        complex ( 8 ) , intent ( in ) :: D1 ( : , : ) , D2 ( : , : )
        ! --- local variable ---
        integer ( 4 ) :: sgn , i
        complex ( 8 ) :: det
        complex ( 8 ) :: D3 ( size ( D1 , dim = 2 ) , size ( D2 , dim = 2 ) )
        ! --- variables for LU factorization ---
        integer ( 4 ) :: dim , LDA , INFO
        integer ( 4 ) , allocatable :: IPIV ( : )

        ! u1 = det ( D1^+ * tz * D2 ) / | det ( D1^+ * tz * D2 ) |
        D3  = matmul ( conjg ( transpose ( D1 ) ) , matmul ( tz , D2 ) )
        dim = size ( D3 , dim = 1 )
        LDA = max ( 1 , dim )
        allocate ( IPIV ( dim ) )

        IPIV = 0
        call ZGETRF ( dim , dim , D3 , LDA , IPIV , INFO )
        if ( INFO /= 0 ) then
            write ( * , * ) 'LU factorization went wrong! aborting ...'
            write ( * , * ) INFO
            stop
        end if

        det = v1
        sgn = 1
        do i = 1 , dim
            det = det * D3 ( i , i )
            if ( IPIV ( i ) /= i ) sgn = -1 * sgn
        end do
        det = cmplx ( sgn , 0 ) * det
        u1 = det / abs ( det )

        if ( abs ( abs ( u1 ) - 1.d0 ) > 1d-10 ) then
            write ( * , * ) 'u1 not unitary' , u1
            stop
        end if

        return
    end subroutine get_U1

    subroutine get_Bfam ( B , kp , mark )
        implicit none
        real ( 8 ) , intent ( in ) :: kp ( : , : )
        complex ( 8 ) , intent ( inout ) :: B
        character ( 1 ) , intent ( in ) :: mark
        ! --- local variable ---
        integer ( 4 ) :: dim , i , j
        logical , allocatable :: tri ( : , : )
        complex ( 8 ) :: u1
        complex ( 8 ) , dimension ( : ) , allocatable :: u2
        complex ( 8 ) , dimension ( : , : ) , allocatable :: Hk , H_bulk , D1 , D2
        complex ( 8 ) , dimension ( : , : , : ) , allocatable :: VR
        real    ( 8 ) , dimension ( : ) , allocatable :: EV , k2
        real    ( 8 ) , dimension ( : , : ) , allocatable :: k1

        if ( FMConfig .eqv. .False. ) then
            dim = 2 * n_ineq
        else
            dim = n_ineq
        end if

        if ( mark == 'F' ) then
            allocate ( VR ( 4 , dim , bandn ) , tri ( 4 , 2 ) , u2 ( 4 ) )
        else if ( mark == 'A' ) then
            allocate ( VR ( 2 , dim , bandn ) , tri ( 2 , 2 ) )
        end if

        allocate ( Hk ( 2 * n_ineq , 2 * n_ineq ) )
        !allocate ( H_bulk ( dim , dim ) , EV ( n_ineq ) )
        allocate ( H_bulk ( dim , dim ) , EV ( dim ) )
        allocate ( D1 ( dim , bandn ) , D2 ( dim , bandn ) )
        allocate ( k1 ( size ( kp , dim = 1 ) , size ( kp , dim = 2 ) ) , k2 ( sysdim ) )

        !get transformation matrix
        k1 = kp
        tri = .False.
        if ( trs .eqv. .True. ) then
            if ( insul2D .eqv. .True. ) then
                do i = 1 , size ( k1 , dim = 1 )
                    if ( ( ( k1 ( i , 1 ) == -0.5 ) .or. ( k1 ( i , 1 ) ==  0.0 ) ) .and. &
                        &( ( k1 ( i , 2 ) == -0.5 ) .or. ( k1 ( i , 2 ) ==  0.0 ) ) ) then
                        tri ( i , 1 ) = .True.
                    else if ( ( ( k1 ( i , 2 ) ==  0.0 ) .or.  ( k1 ( i , 2 ) == -0.5 ) ) .and. &
                               &( k1 ( i , 1 )  <  0.0 ) .and. ( k1 ( i , 1 )  > -0.5 ) ) then
                        tri ( i , 2 ) = .True.
                        k1 ( i , 1 ) = -k1 ( i , 1 )
                    end if
                end do
            else if ( semi1D .eqv. .True. ) then
                if ( sysdim == 3 ) then
                    do i = 1 , size ( k1 , dim = 1 )
                        !if ( ( ( k1 ( i , 1 ) == -0.5 ) .or. ( k1 ( i , 1 ) == 0.0 ) ) .and. &
                        !    &( ( k1 ( i , 2 ) == -0.5 ) .or. ( k1 ( i , 2 ) == 0.0 ) ) .and. &
                        !    &( ( k1 ( i , 3 ) == -0.5 ) .or. ( k1 ( i , 3 ) == 0.0 ) ) ) then
                        !    tri ( i , 1 ) = .True.
                        if ( ( k1 ( i , 1 ) <= 0.0 ) .and. ( k1 ( i , 1 ) >= -0.5 ) ) then
                            tri ( i , 2 ) = .True.
                            if ( .not. ( ( abs ( k1 ( i , 1 ) ) < 1d-8 ) .or. ( abs ( k1 ( i , 1 ) + 0.5 ) < 1d-8 ) ) ) &
                            &k1 ( i , 1 ) = -k1 ( i , 1 )
                            if ( .not. ( ( abs ( k1 ( i , 2 ) ) < 1d-8 ) .or. ( abs ( k1 ( i , 2 ) + 0.5 ) < 1d-8 ) ) ) &
                            &k1 ( i , 2 ) = -k1 ( i , 2 )
                            if ( .not. ( ( abs ( k1 ( i , 3 ) ) < 1d-8 ) .or. ( abs ( k1 ( i , 3 ) + 0.5 ) < 1d-8 ) ) ) &
                            &k1 ( i , 3 ) = -k1 ( i , 3 )
                        end if
                    end do
                else
                    write ( * , * ) 'error occurs at get_Bfam'
                    stop
                end if
            else
                write ( * , * ) 'error occurs at get_Bfam'
                stop
            end if
        end if

        VR = v0
        do i = 1 , size ( k1 , dim = 1 )
            D1 = v0
            H_bulk = v0
            k2 = k1 ( i , : )
            call get_H ( Hk ,  2.d0 * pi * k2 , 'bk' )
            H_bulk = Hk ( 1 : dim , 1 : dim )
            call get_diag ( H_bulk , EV , D1 )
            if ( trs .eqv. .True. ) then
                if ( tri ( i , 1 ) .eqv. .True. ) then
                    do j = 1 , bandn / 2
                        D1 ( : , 2 * j ) = matmul ( tro , conjg ( D1 ( : , 2 * j - 1 ) ) )
                    end do
                else if ( tri ( i , 2 ) .eqv. .True. ) then
                    D2 = D1
                    do j = 1 , bandn / 2
                        D1 ( : , 2 * j - 1 ) = matmul ( tro , conjg ( D2 ( : , 2 * j - 1 ) ) )
                        D1 ( : , 2 * j ) = matmul ( tro , conjg ( D2 ( : , 2 * j ) ) )
                    end do
                end if
            end if
            VR ( i , : , : ) = D1
        end do

        if ( mark == 'F' ) then
            u2 = v0
            D1 = VR ( 1 , : , : )
            D2 = VR ( 2 , : , : )
            call get_U1 ( D1 , D2 , u1 )
            u2 ( 1 ) = u1
            D1 = VR ( 2 , : , : )
            D2 = VR ( 4 , : , : )
            call get_U1 ( D1 , D2 , u1 )
            u2 ( 2 ) = u1
            D1 = VR ( 3 , : , : )
            D2 = VR ( 4 , : , : )
            call get_U1 ( D1 , D2 , u1 )
            u2 ( 3 ) = u1
            D1 = VR ( 1 , : , : )
            D2 = VR ( 3 , : , : )
            call get_U1 ( D1 , D2 , u1 )
            u2 ( 4 ) = u1

            B = v0
            B = log ( u2 ( 1 ) * u2 ( 2 ) / u2 ( 3 ) / u2 ( 4 ) )
            if ( trs .eqv. .True. ) B = log ( u2 ( 1 ) ) + log ( u2 ( 2 ) ) - log ( u2 ( 3 ) ) - log ( u2 ( 4 ) ) - B
            deallocate ( Hk , VR , u2 )
        else if ( mark == 'A' ) then
            D1 = VR ( 1 , : , : )
            D2 = VR ( 2 , : , : )
            call get_U1 ( D1 , D2 , u1 )

            B = v0
            B = log ( u1 )

            deallocate ( VR )
        else
            write ( * , * ) 'error occurs at get_Bfam'
            stop
        end if

        return
    end subroutine get_Bfam

    subroutine get_ti
        implicit none
        ! --- local variables ---
        ! Re ( km ( i , j ) ) = kx_i , Im ( km ( i , j ) ) = ky_j
        integer ( 4 ) :: i , j , n1 , n2
        real    ( 8 ) :: D
        real    ( 8 ) , allocatable :: km ( : , : , : ) , kp ( : , : )
        complex ( 8 ) :: F , ti
        character ( 50 ) :: term

        n1 = 20
        n2 = 30
        call get_2dBZ_discrete ( km , n1 , n2 )

        ! --- obtain the bz integral of F ---
        ti = v0
        allocate ( kp ( 4 , 2 ) )
        do i = 1 , n1 - 1
            do j = 1 , n2 - 1
                kp = 0.0
                kp ( 1 , : ) = km ( i , j , : )
                kp ( 2 , : ) = km ( i + 1 , j , : )
                kp ( 3 , : ) = km ( i , j + 1 , : )
                kp ( 4 , : ) = km ( i + 1 , j + 1 , : )
                call get_Bfam ( F , kp , 'F' )
                if ( abs ( dble ( F ) ) > 1d-10 ) then
                    write ( * , * ) 'sth wrong with F:' , F
                    stop
                end if
                ti = ti + F / cmplx ( 0 , 2.d0 * pi )
            end do
        end do
        deallocate ( kp )

        if ( trs .eqv. .True. ) then
            write ( term , '(F8.4)' ) dble( ti )
            read ( term , * ) D
            i = floor ( D / 2.d0 )
            write ( * , '(F8.5)' )  D - 2.d0 * dble ( i )
        else
            write ( * , '(F8.5)' )  dble( ti )
        end if

        return
    end subroutine get_ti

    subroutine get_ti_defect
        implicit none
        integer ( 4 ) :: i
        complex ( 8 ) :: A , ti
        character ( 50 ) :: term
        real ( 8 ) :: step , D , theta ( loopn )
        real ( 8 ) , allocatable :: kp ( : , : ) , kpath ( : , : )


        ! --- obtain kpoints on the loop ---
        !     for each k on k_i - k_j plane
        !     k_i = loopc_i + loopr * cos ( theta )
        !     k_j = loopc_j + loopr * sin ( theta )
        step = 2 * pi / dble ( loopn )
        do i = 1 , loopn
            theta ( i ) = dble ( i - 1 ) * step
        end do

        allocate ( kpath ( size ( theta ) + 1 , size ( loopc ) ) )
        kpath = 0.0
        if ( nordir == 1 ) then
            kpath ( : , 1 ) = loopc ( 1 )
            do i = 1 , size ( theta )
                kpath ( i , 2 ) = loopc ( 2 ) + loopr * cos ( theta ( i ) )
                kpath ( i , 3 ) = loopc ( 3 ) + loopr * sin ( theta ( i ) )
            end do
        else if ( nordir == 2 ) then
            kpath ( : , 2 ) = loopc ( 2 )
            do i = 1 , size ( theta )
                kpath ( i , 1 ) = loopc ( 1 ) + loopr * cos ( theta ( i ) )
                kpath ( i , 3 ) = loopc ( 3 ) + loopr * sin ( theta ( i ) )
            end do
        else if ( nordir == 3 ) then
            kpath ( : , 3 ) = loopc ( 3 )
            do i = 1 , size ( theta )
                kpath ( i , 1 ) = loopc ( 1 ) + loopr * cos ( theta ( i ) )
                kpath ( i , 2 ) = loopc ( 2 ) + loopr * sin ( theta ( i ) )
            end do
        end if
        kpath ( size ( theta ) + 1 , : ) = kpath ( 1 , : )

        allocate ( kp ( 2 , size ( loopc ) ) )
        ti = v0
        do i = 1 , loopn
            kp = kpath ( i : i + 1 , : )
            call get_Bfam ( A , kp , 'A' )
            if ( abs ( dble ( A ) ) > 1d-10 ) then
                write ( * , * ) 'sth wrong with A:' , A
                stop
            end if
            !write(*,*)A/(2.0*pi*vi)
            ti = ti + A
        end do

        if ( trs .eqv. .True. ) then
            write ( term , '(F8.4)' ) dble ( ti / ( pi * vi ) )
            read ( term , * ) D
            i = floor ( D / 2.d0 )
            write ( * , '(F8.5)' )  D - 2.d0 * dble ( i )
        else
            write ( * , '(F8.5)' )  dble( ti )
        end if

    end subroutine get_ti_defect

end module TopoProp
!---------------------------------------------------------------------------------------------------------
