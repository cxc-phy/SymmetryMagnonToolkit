!---------------------------------------------------------------------------------------------------------
module SLAB
    use Tool
    use var
    use Hamiltonian
    use Diagonalization

contains
    subroutine get_Slab ( kp , count )
        implicit none
        real ( 8 ) , intent ( in ) :: kp ( : )
        integer ( 4 ) , intent ( in ) :: count
        ! --- local variable ---
        integer ( 4 ) :: i , d1 , d2 , d3 , d4 , dim , d0
        real    ( 8 ) , dimension ( : ) , allocatable :: EV
        complex ( 8 ) , dimension ( : , : ) , allocatable :: Hk , VR , tx
        complex ( 8 ) , dimension ( : , : ) , allocatable :: H0u , H0b , H0l , H01 , H10 , Hslab , H00
        complex ( 8 ) , dimension ( : , : ) , allocatable :: test
        integer ( 4 ) :: c1,c2,c3,c4

        if ( FMConfig .eqv. .False. ) then
            dim = 2 * n_ineq
        else
            dim = n_ineq
        end if

        allocate ( Hk ( 2 * n_ineq , 2 * n_ineq ) )
        allocate ( H0u ( dim , dim ) ) ! layer of upper surface
        allocate ( H0b ( dim , dim ) ) ! bulk laryer
        allocate ( H00 ( dim , dim ) ) ! single layer
        allocate ( H0l ( dim , dim ) ) ! layer of lower surface
        allocate ( H01 ( dim , dim ) )
        allocate ( H10 ( dim , dim ) )
        allocate ( Hslab ( dim * nlayer , dim * nlayer ) )
        H0u = v0
        H0b = v0
        H00 = v0
        H0l = v0
        H01 = v0
        H10 = v0
        Hslab = v0

        ! --- get H_slab ---
        if ( nlayer == 1 ) then
            call get_H ( Hk , kp , '00' ) ! get H00 of a single layer
            H00 = Hk ( 1 : dim , 1 : dim )
            Hslab = H00
        else
            call get_H ( Hk , kp , '0l' ) ! get H00 of layer on lower surface
            H0l = Hk ( 1 : dim , 1 : dim )
            call get_H ( Hk , kp , '0u' ) ! get H00 of layer on upper surface
            H0u = Hk ( 1 : dim , 1 : dim )
            call get_H ( Hk , kp , '0b' ) ! get H00 of bulk layer
            H0b = Hk ( 1 : dim , 1 : dim )
            call get_H ( Hk , kp , '01' ) ! get H_01
            H01 = Hk ( 1 : dim , 1 : dim )
            call get_H ( Hk , kp , '10' ) ! get H_10
            H10 = Hk ( 1 : dim , 1 : dim )
            do i = 1 , nlayer , 1
                d2 = dim * ( i - 1 )
                d3 = dim * i
                if ( i == 1 ) then
                    Hslab ( d2 + 1 : d3 , d2 + 1 : d3 ) = H0l
                else if ( i == nlayer ) then
                    Hslab ( d2 + 1 : d3 , d2 + 1 : d3 ) = H0u
                else
                    Hslab ( d2 + 1 : d3 , d2 + 1 : d3 ) = H0b
                end if

                if ( i /= 1 ) then
                    d1 = dim * ( i - 2 )
                    Hslab ( d2 + 1 : d3 , d1 + 1 : d2 ) = H10
                end if
                if ( i /= nlayer ) then
                    d4 = dim * ( i + 1 )
                    Hslab ( d2 + 1 : d3 , d3 + 1 : d4 ) = H01
                end if
            end do
        end if


        !allocate ( EV ( n_ineq * nlayer ) , VR ( dim * nlayer , dim * nlayer ) )
        allocate ( EV ( dim * nlayer ) , VR ( dim * nlayer , dim * nlayer ) )
        call get_diag ( Hslab , EV , VR )
        write ( 99 , * ) count , EV! , kp / 2.0 / pi
        return

        !!!!!
        !allocate (test(size(Hslab,dim=1)-4,size(Hslab,dim=2)-4))
        !test = v0
        !c1 = dim*(nlayer-1)+4
        !c2 = dim*(nlayer-1)+7

        !test(:c1,:c1) = Hslab(:c1,:c1)
        !test(:c1,c1+1:) = Hslab(:c1,c2:c2+3)
        !test(c1+1:,:c1) = Hslab(c2:c2+3,:c1)
        !test(c1+1:,c1+1:) = Hslab(c2:c2+3,c2:c2+3)
        !allocate ( EV ( dim * nlayer - 4 ) , VR ( dim * nlayer - 4 , dim * nlayer - 4 ) )

        !allocate (test(size(Hslab,dim=1)-8,size(Hslab,dim=2)-8))
        !test = v0

        !test(3:,3:) = Hslab(11:,11:)
        !test(3:,1:2) = Hslab(11:,5:6)
        !test(1:2,1:2) = Hslab(5:6,5:6)
        !test(1:2,3:) = Hslab(5:6,11:)
        !allocate ( EV ( dim * nlayer - 8 ) , VR ( dim * nlayer - 8 , dim * nlayer - 8 ) )
        !!!!!





        call get_diag ( test , EV , VR )
        write ( 99 , * ) count , EV! , kp / 2.0 / pi

        deallocate ( Hk , H0u , H0b , H0l , H01 , H10 , Hslab , VR , EV )
        return
    end subroutine get_slab

end module SLAB
!---------------------------------------------------------------------------------------------------------
