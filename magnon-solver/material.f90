!---------------------------------------------------------------------------------------------------------
include 'parameter.f90'
include 'hmatrix.f90'
include 'LSW.f90'
include 'SLAB.f90'
include 'SPEC.f90'
include 'TOPO.f90'
!---------------------------------------------------------------------------------------------------------
program main
    use var
    use rotation_matrix
    use Hamiltonian
    implicit none

    logical :: fstat
    real ( 8 ) , allocatable :: delta_k ( : ) , kp ( : )
    integer ( 4 ) :: i , j , count , n
    character ( 60 ) :: format_0
    format_0 = '(A1,F5.1,X,A1,$)'

    open ( unit = 97 , file = 'output.txt' , status = 'replace' )
    call get_input
    do n = 1 , size ( modelist )
        mode = modelist ( n )
        if ( mode == -1 ) exit
        call get_cmat ( mode )
        if ( mode == 0 ) then
            write ( 97 , fmt = "(/,'BZ stability checking...')" )
            write ( * , fmt = "('BZ stability checking...')" )
            block
                use LSW
                if ( bz2D .eqv. .True. ) then
                    allocate ( kp ( 2 ) )
                else
                    allocate ( kp ( 3 ) )
                end if
                kp = -0.5
                count = 0
                do while ( kp ( 1 ) < 0.5 )
                    kp ( 2 ) = -0.5
                    do while ( kp ( 2 ) < 0.5 )
                        if ( bz2D .eqv. .True. ) then
                            call get_bulk ( 2.d0 * pi * kp , count )
                        else
                            kp ( 3 ) = -0.5
                            do while ( kp ( 3 ) < 0.5 )
                                write(*,'(A1,F5.1,A1,$)') char ( 13 ) , dble ( count ) / 132650.d0 * 100.d0 , '%'
                                call get_bulk ( 2.d0 * pi * kp , count )
                                count = count + 1
                                kp ( 3 ) = kp ( 3 ) + 0.02
                            end do
                        end if
                        kp ( 2 ) = kp ( 2 ) + 0.02
                    end do
                    kp ( 1 ) = kp ( 1 ) + 0.02
                end do
            end block
        else if ( mode == 1 ) then
            write ( 97 , fmt = "(/,'bulk computing...')" )
            write ( * , fmt = "('bulk computing...')" )
            block
                use LSW
                open ( unit = 99 , file = 'dispersion.txt' , status = 'unknown' )
                allocate ( kp ( size ( bulkkpath , dim = 2 ) ) )
                allocate ( delta_k ( size ( kp ) ) )
                count = 0
                ks = 5000
                do i = 1 , size ( bulkkpath , dim = 1 ) - 1 , 1
                    write ( * , '(A1,F8.5,X,F8.5,X,F8.5)' ) char ( 13 ) , bulkkpath ( i , : ) / 2.0 / pi
                    delta_k = ( bulkkpath ( i + 1 , : ) - bulkkpath ( i , : ) ) / dble ( ks )
                    do j = 0 , ks - 1 , 1
                        write ( * , format_0 ) char ( 13 ) , dble ( 100 * j ) / dble ( ks ) , '%'
                        kp = bulkkpath ( i , : ) + dble ( j ) * delta_k
                        call get_bulk ( kp , count )
                        count = count + 1
                    end do
                end do
                write ( * , '(A1,F8.5,X,F8.5,X,F8.5)' ) char ( 13 ) , bulkkpath ( size ( bulkkpath , dim = 1 ) , : ) / 2.0 / pi
                close ( 99 )
            end block
        else if ( mode == 2 ) then
            write ( 97 , fmt = "(/,'slab computing...')" )
            write ( * , fmt = "('slab computing...')" )
            block
                use Slab
                open ( unit = 99 , file = 'slab.txt' , status = 'unknown' )
                count = 0
                ks = 300
                dir = slabdir

                allocate ( kp ( size ( slabkpath , dim = 2 ) ) )
                allocate ( delta_k ( size ( kp ) ) )
                if ( size ( slabkpath , dim = 1 ) == 1 ) then
                    kp = slabkpath ( 1 , : )
                    write ( * , * ) kp
                    call get_slab ( kp , count )
                else
                    do i = 1 , size ( slabkpath , dim = 1 ) - 1 , 1
                        write ( * , '(A1,F8.5,X,F8.5)' ) char ( 13 ) , slabkpath ( i , : ) / 2.0 / pi
                        delta_k = ( slabkpath ( i + 1 , : ) - slabkpath ( i , : ) ) / dble ( ks )
                        do j = 0 , ks - 1 , 1
                            write ( * , format_0 ) char ( 13 ) , dble ( 100 * j ) / dble ( ks ) , '%'
                            kp = slabkpath ( i , : ) + dble ( j ) * delta_k
                            call get_slab ( kp , count )
                            count = count + 1
                        end do
                    end do
                    write ( * , '(A1,F8.5,X,F8.5)' ) char ( 13 ) , slabkpath ( size ( slabkpath , dim = 1 ) , : ) / 2.0 / pi
                end if
                close ( 99 )
            end block
        else if ( mode == 3 ) then
            write ( 97 , fmt = "(/,'surface green function computing...')" )
            write ( * , fmt = "('surface green function computing...')" )
            block
                use SPEC
                open ( unit = 99 , file = 'spec.txt' , status = 'unknown' )
                call get_rot_diag ( n_ineq , '-' )
                count = 0
                ks = 2000
                dir = surfdir
                allocate ( kp ( size ( surfkpath , dim = 2 ) ) )
                allocate ( delta_k ( size ( kp ) ) )
                if ( size ( surfkpath , dim = 1 ) == 1 ) then
                    kp = surfkpath ( 1 , : )
                    call get_A_omega_k ( kp , count )
                else
                    do i = 1 , size ( surfkpath , dim = 1 ) - 1 , 1
                        write ( * , '(A1,F8.5,X,F8.5)' ) char ( 13 ) , surfkpath ( i , : ) / 2.0 / pi
                        delta_k = ( surfkpath ( i + 1 , : ) - surfkpath ( i , : ) ) / dble ( ks )
                        do j = 0 , ks - 1 , 1
                            write ( * , format_0 ) char ( 13 ) , dble ( 100 * j ) / dble ( ks ) , '%'
                            kp = surfkpath ( i , : ) + dble ( j ) * delta_k
                            call get_A_omega_k ( -kp , count )
                            count = count + 1
                        end do
                    end do
                    write ( * , '(A1,F8.5,X,F8.5)' ) char ( 13 ) , surfkpath ( size ( surfkpath , dim = 1 ) , : ) / 2.0 / pi
                    close ( 99 )
                end if
            end block
        else if ( mode == 4 ) then
            write ( 97 , fmt = "(/,'topological invariant computing...')" )
            write ( * , fmt = "('topological invariant computing...')" )
            block
                use TopoProp
                if ( insul2D .eqv. .True. ) then
                    call get_ti
                else if ( semi1D .eqv. .True. ) then
                    call get_ti_defect
                end if
            end block
        end if
        write ( 97 , fmt = "('finish')" )
        if ( allocated ( delta_k ) ) deallocate ( delta_k )
        if ( allocated ( kp ) ) deallocate ( kp )
        if ( allocated ( tz ) ) deallocate ( tz )
        if ( allocated ( Id ) ) deallocate ( Id )
    end do

    write ( * , '(A1,A8)' ) char(13) , 'finished'
    close ( 97 )
    deallocate ( Jvalue )
    if ( allocated ( kpath ) ) deallocate ( kpath )
end program
!---------------------------------------------------------------------------------------------------------
