!---------------------------------------------------------------------------------------------------------
module LSW
    use Tool
    use var
    use Hamiltonian
    use Diagonalization

contains
    subroutine get_bulk ( k , count )
        implicit none
        real ( 8 ) , dimension ( : ) , intent ( in ) :: k
        integer ( 4 ) , intent ( in ) :: count
        ! --- local variable ---
        integer ( 4 ) :: dim , i , j
        real ( 8 ) , dimension ( : ) , allocatable :: EV
        complex ( 8 ) , dimension ( : , : ) , allocatable :: Hbulk , Hk , VR

        if ( FMConfig .eqv. .False. ) then
            dim = 2 * n_ineq
        else
            dim = n_ineq
        end if

        allocate ( Hk ( 2 * n_ineq , 2 * n_ineq ) , Hbulk ( dim , dim ) )
        !allocate ( EV ( n_ineq ) , VR ( dim , dim ) )
        allocate ( EV ( dim ) , VR ( dim , dim ) )

        call get_H ( Hk , k , 'bk' )
        Hbulk = Hk ( 1 : dim , 1 : dim )
        call get_diag ( Hbulk , EV , VR )

        if ( mode == 0 ) then
            if ( err == -1 ) write ( 97 , * ) k / 2.0 / pi
            err = 0
        else
            write ( 99 , * ) count , EV , k / 2.0 / pi
        end if

        deallocate ( Hbulk , Hk , EV , VR )
        return
    end subroutine get_bulk

end module LSW
!---------------------------------------------------------------------------------------------------------
