!---------------------------------------------------------------------------------------------------------
module SPEC
    use Tool
    use var
    use Hamiltonian
    real ( 8 ) :: delta
    integer ( 4 ) :: n_iter

contains
    subroutine stop_check ( mat , stop_label )
        implicit none
        integer ( 4 ) , intent ( inout ) :: stop_label
        complex ( 8 ) , intent ( in ) :: mat ( : , : )
        integer ( 4 ) :: i , j , dim

        delta = 0.d0
        dim = size ( mat , dim = 1 )
        do i = 1 , dim , 1
            do j = 1 , dim , 1
                delta = delta + abs ( mat ( i , j ) )
            end do
        end do
        delta = delta / dim ** 2
        if ( delta < 1d-8 ) then
            stop_label = 0
        end if

        return
    end subroutine stop_check

    subroutine get_A_omega_k ( kp , count )
        implicit none
        real    ( 8 ) , intent ( in ) :: kp ( : )
        integer ( 4 ) , intent ( in ) :: count
        ! --- local variable ---
        real    ( 8 ) :: Ab , Au , Al
        real    ( 8 ) :: omega
        integer ( 4 ) :: dim , i , j
        integer ( 4 ) :: stop_label_01 , stop_label_10
        complex ( 8 ) :: trace
        complex ( 8 ) , dimension ( : , : ) , allocatable :: mat , wmat
        complex ( 8 ) , dimension ( : , : ) , allocatable :: Hk , H01 , H10 , H0l , H0u , Hnn
        complex ( 8 ) , dimension ( : , : ) , allocatable :: M01 , M10 , M00 , Mnn

        if ( FMConfig .eqv. .False. ) then
            dim = 2 * n_ineq
        else
            dim = n_ineq
        end if

        allocate ( Hk ( 2 * n_ineq , 2 * n_ineq ) )
        allocate ( H01 ( dim , dim ) )
        allocate ( H10 ( dim , dim ) )
        allocate ( H0l ( dim , dim ) )
        allocate ( H0u ( dim , dim ) )
        allocate ( Hnn ( dim , dim ) )
        allocate ( mat ( dim , dim ) )
        allocate ( wmat ( dim , dim ) )

        allocate ( M01 ( dim , dim ) )
        allocate ( M10 ( dim , dim ) )
        allocate ( M00 ( dim , dim ) )
        allocate ( Mnn ( dim , dim ) )

        call get_H ( Hk , kp , '0l' ) ! get H00_lower surface
        H0l = Hk ( 1 : dim , 1 : dim )
        call get_H ( Hk , kp , '0u' ) ! get H00_upper surface
        H0u = Hk ( 1 : dim , 1 : dim )
        call get_H ( Hk , kp , '0b' ) ! get H00
        Hnn = Hk ( 1 : dim , 1 : dim )
        call get_H ( Hk , kp , '01' ) ! get H_01
        H01 = Hk ( 1 : dim , 1 : dim )
        call get_H ( Hk , kp , '10' ) ! get H_10
        H10 = Hk ( 1 : dim , 1 : dim ) !transpose ( conjg ( H01 ) )


        H0l = matmul ( H0l , rot_diag )
        H0u = matmul ( H0u , rot_diag )
        Hnn = matmul ( Hnn , rot_diag )
        H01 = matmul ( H01 , rot_diag )
        H10 = matmul ( H10 , rot_diag )
        H0l = matmul ( conjg ( rot_diag ) , H0l )
        H0u = matmul ( conjg ( rot_diag ) , H0u )
        Hnn = matmul ( conjg ( rot_diag ) , Hnn )
        H01 = matmul ( conjg ( rot_diag ) , H01 )
        H10 = matmul ( conjg ( rot_diag ) , H10 )

        omega = wi
        do while ( omega <= wf )
            if ( FMConfig .eqv. .False. ) then
                wmat = omega * tz + vi * bf * Id
            else
                wmat = ( omega + vi * bf ) * Id
            end if
            wmat = matmul ( wmat , rot_diag )
            wmat = matmul ( conjg ( rot_diag ) , wmat )

            ! --- lower surface ---
            ! prepare initial value for lower surface
            M00 = H0l
            Mnn = Hnn
            M01 = H01
            M10 = H10

            n_iter = 1
            stop_label_01 = -1
            stop_label_10 = -1

            ! start iteration for lower surface
            do while ( ( stop_label_01 == -1 ) .or. ( stop_label_10 == -1 ) )
                mat = wmat - Mnn
                call get_inverse ( mat )
                Mnn = Mnn + matmul ( matmul ( M01 , mat ) , M10 ) &
                         &+ matmul ( matmul ( M10 , mat ) , M01 )
                M00 = M00 + matmul ( matmul ( M01 , mat ) , M10 )
                M01 = matmul ( matmul ( M01 , mat ) , M01 )
                M10 = matmul ( matmul ( M10 , mat ) , M10 )

                call stop_check ( M01 , stop_label_01 )
                call stop_check ( M10 , stop_label_10 )
                ! prevent endless loop
                n_iter = n_iter + 1
                if ( n_iter > 30 ) write ( * , * ) 'delta' , delta
                if ( n_iter > 50 ) then
                    write ( * , * ) 'endless iteration'
                    stop
                end if
            end do

            ! get lower surface layer spectra
            Al = 0.d0
            mat = wmat - M00
            call get_inverse ( mat )
            mat = 0.5 * S * mat
            mat = aimag ( mat )
            do i = 1 , dim , 1
                Al = Al - real ( mat ( i , i ) )
            end do
            if ( Al < 0.d0 ) write ( 98 , * ) 'negative spetral weights appears'

            ! --- upper surface ---
            ! prepare initial value for upper surface
            M00 = H0u
            Mnn = Hnn
            M01 = H01
            M10 = H10

            n_iter = 1
            stop_label_01 = -1
            stop_label_10 = -1

            ! start iteration for upper surface
            n_iter = 1
            do while ( ( stop_label_01 == -1 ) .or. ( stop_label_10 == -1 ) )
                mat = wmat - Mnn
                call get_inverse ( mat )
                Mnn = Mnn + matmul ( matmul ( M10 , mat ) , M01 ) &
                         &+ matmul ( matmul ( M01 , mat ) , M10 )
                M00 = M00 + matmul ( matmul ( M10 , mat ) , M01 )
                M01 = matmul ( matmul ( M01 , mat ) , M01 )
                M10 = matmul ( matmul ( M10 , mat ) , M10  )
                call stop_check ( M01 , stop_label_01 )
                call stop_check ( M10 , stop_label_10 )
                ! prevent endless loop
                n_iter = n_iter + 1
                if ( n_iter > 30 ) write ( * , * ) 'delta' , delta
                if ( n_iter > 50 ) then
                    write ( * , * ) 'endless iteration'
                    stop
                end if
            end do

            ! get upper surface layer spectra
            Au = 0.d0
            mat = wmat - M00
            call get_inverse ( mat )
            mat = 0.5 * S * mat
            mat = aimag ( mat )
            do i = 1 , dim , 1
                Au = Au - real ( mat ( i , i ) )
            end do
            if ( Au < 0.d0 ) write ( 98 , * ) 'negative spetral weights appears'

            ! get bulk layer spectra
            Ab = 0.d0
            mat = wmat - Mnn
            call get_inverse ( mat )
            mat = 0.5 * S * mat
            mat = aimag ( mat )
            do i = 1 , dim , 1
                Ab = Ab - real ( mat ( i , i ) )
            end do

            ! --- restore upper/lower/bulk spectra function ---
            write ( 99 , * ) count , omega , Au , Al , Ab
            omega = omega + ws
        end do

        deallocate ( mat , wmat , Hk , H01 , H10 , H0u , H0l , Hnn , M00 , Mnn , M01 , M10 )
        return
    end subroutine get_A_omega_k

end module SPEC
!---------------------------------------------------------------------------------------------------------
