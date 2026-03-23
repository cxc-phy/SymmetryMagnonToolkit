module var
    implicit none
    ! --- global variables ---
    real    ( 8 ) :: S = 1.0
    real    ( 8 ) , parameter :: pi = 3.1415926
    complex ( 8 ) , parameter :: v0 = cmplx ( 0 , 0 )
    complex ( 8 ) , parameter :: v1 = cmplx ( 1 , 0 )
    complex ( 8 ) , parameter :: vi = cmplx ( 0 , 1 )
    complex ( 8 ) , allocatable :: Id ( : , : ) , tz ( : , : )
    complex ( 8 ) , allocatable :: tro ( : , : )
    logical :: trs
    integer ( 4 ) :: err = 0

    ! --- calculation mode specification ---
    integer ( 4 ) :: mode = -1
    integer ( 4 ) :: sysdim = 0
    integer ( 4 ) , dimension ( 4 ) :: modelist = -1
    logical :: bz3D = .False. , bz2D = .False.
    logical :: bulk = .False.
    logical :: slab2D = .False. , slab1D = .False.
    logical :: surf2D = .False. , surf1D = .False. , surfkp = .False. , surfe = .False.
    logical :: topo = .False. , semi1D = .False. , insul2D = .False.

    ! --- kpath ---
    real    ( 8 ) , dimension ( : , : ) , allocatable :: kpath
    real    ( 8 ) , dimension ( : , : ) , allocatable :: bulkkpath
    real    ( 8 ) , dimension ( : , : ) , allocatable :: slabkpath
    real    ( 8 ) , dimension ( : , : ) , allocatable :: surfkpath
    integer ( 4 ) :: knum ! number of k vertices in kpath
    integer ( 4 ) :: ks ! kpath step

    ! --- slab / surface green function ---
    integer ( 4 ) :: dir = 0
    integer ( 4 ) :: slabdir = 0 , surfdir = 0 ! stacking direction
    integer ( 4 ) :: nlayer ! layer of slab
    real ( 8 ) :: bf = 1d-2 ! broadening factor in green function
    real ( 8 ) :: wi ! initial omega value
    real ( 8 ) :: wf ! final omega value
    real ( 8 ) :: ws ! omega step
    real ( 8 ) :: ecut ! energy cut in equal energy mode

    ! --- Z2 invariant with time reversal symmetry ---
    integer ( 4 ) :: nordir = 0 ! normal direction of the loop
    integer ( 4 ) :: loopn ! number of kpoints on the loop
    integer ( 4 ) :: bandn ! number of bands below the cross
    integer ( 4 ) :: bandi
    integer ( 4 ) :: bandf
    real ( 8 ) :: loopr = 1d-3 ! loop radius
    real ( 8 ) , allocatable :: loopc ( : ) ! loop center

    ! --- J value / DMI value / Single Ion Anisotropy value ---
    logical :: Jani = .False. ! T : J_ex has different values on x/y/z , F : constant value , default : F
    logical :: DMI = .False. ! T : turn on DMI , F : turn off DMI , default : F
    logical :: SIA = .False. ! T : turn on SIA , F : turn off SIA , default : F
    integer ( 4 ) :: Jorder ! what order of nearest neighbors to be considered
    integer ( 4 ) :: Dorder ! what order of nearest neighbors to be considered
    integer ( 4 ) , allocatable :: neighbor_list ( : , : ) ! store n_neighbor of atom
    real ( 8 ) :: Avalue ! Single Ion Anisotropy value
    real ( 8 ) , allocatable :: Jvalue ( : , : , : ) ! contains J value of different order
    real ( 8 ) , allocatable :: Dvalue ( : , : )

    ! --- spin direction deviation ---
    logical :: FMConfig
    integer ( 4 ) :: n_ineq ! inequivalent atoms in the cell
    real ( 8 ) , allocatable :: s_coord ( : , : ) ! spin component in cartesian coord
    real ( 8 ) , allocatable :: theta ( : ) ! rotate along y-axis
    real ( 8 ) , allocatable :: phi ( : ) ! rotate along z-axis

contains
    subroutine remove_comment ( term )! skip comment component
        implicit none
        character ( 100 ) , intent ( inout ) :: term
        if ( index ( term , '#' ) > 0 ) term = term ( : index ( term , '#' ) - 1 )
        return
    end subroutine

    subroutine read_info
        implicit none
        integer ( 4 ) :: i , j , m , n , io
        !logical :: logic
        real ( 8 ) , allocatable :: list ( : )
        character ( 100 ) :: term , values

        ! --- read from info ---
        ! --- achieve calculate mode, Jvalue, spin config and kpath ---
        open ( unit = 9 , file = 'INFO_FORT' , form = 'formatted' , access = 'sequential' , status = 'old' )
        do
            read ( 9 , '(A)' , iostat = io ) term
            if ( io > 0 ) then
                write ( * , * ) 'oops, error occurs while reading in info'
                stop
            else if ( io == 0 ) then
                call remove_comment ( term )
                values = term ( index ( term , '=' ) + 1 : )
                term = term ( : index ( term , '=' ) - 1 )
                term = trim ( term )
                ! --- global parameters ---
                if ( term == 'SysDim' ) then
                    read ( values , * , iostat = io ) sysdim
                    if ( ( io /= 0 ) .or. ( sysdim <= 0 ) .or. ( sysdim > 3 ) )then
                        write ( * , * ) 'error occurs at SysDim: bad input'
                        stop
                    end if
                else if ( term == 'S_Value' ) then
                    read ( values , * , iostat = io ) S
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at S_Value: bad input'
                        stop
                    end if
                ! --- spin configuration/direction ---
                else if ( term == 'Spin_Config' ) then
                    i = index ( values , '*' )
                    if ( i == 0 ) then
                        write ( * , * ) 'error occurs at Spin_Config: bad input'
                        stop
                    end if
                    read ( values ( : i - 1 ) , '(I4)' ) n_ineq
                    read ( values ( i + 1 : ) , '(I4)' ) n
                    if ( ( n /= 2 ) .and. ( n /= 3 ) ) then
                        write ( * , * ) 'error occurs at Spin_Config: bad input'
                        stop
                    end if
                    allocate ( s_coord ( n_ineq , n ) )
                    i = 1
                    do while ( i <= n_ineq )
                        read ( 9 , '(A)' ) values
                        call remove_comment ( values )
                        read ( values , * , iostat = io ) s_coord ( i , : )
                        if ( io == 0 ) i = i + 1
                    end do
                    ! check FM spin configuration or not
                    FMConfig = .True.
                    do i = 1 , n_ineq , 1
                        if ( any ( s_coord ( i , : ) /= s_coord ( 1 , : ) ) ) then
                            FMConfig = .False.
                            exit
                        end if
                    end do
                ! --- exchange Jvalue ---
                ! whether consider anisotrpy of Jvalue
                else if ( term == 'J_Ani' ) then
                    read ( values , * , iostat = io ) Jani
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at J_Ani: bad input'
                        stop
                    end if
                else if ( term == 'J_value' ) then
                    i = index ( values , '*' )
                    if ( i == 0 ) then
                        read ( values , '(I4)' ) Jorder
                        n = 1
                    else
                        read ( values ( : i - 1 ) , '(I4)' ) Jorder
                        read ( values ( i + 1 : ) , '(I4)' ) n
                    end if
                    if ( ( n <= 0 ) .or. ( Jorder <= 0 ) ) then
                        write ( * , * ) 'error occurs at J_value: bad input'
                        stop
                    end if
                    if ( Jani .eqv. .True. ) then
                        allocate ( list ( 3 * n ) )
                    else
                        allocate ( list ( n ) )
                    end if
                    allocate ( Jvalue ( Jorder , n , 3 ) )
                    i = 1
                    Jvalue = 0.0
                    do while ( i <= Jorder )
                        read ( 9 , '(A)' ) values
                        call remove_comment ( values )
                        read ( values , * , iostat = io ) list ( : )
                        if ( Jani .eqv. .True. ) then
                            m = size ( Jvalue , dim = 2 )
                            n = size ( Jvalue , dim = 3 )
                            Jvalue ( i , : , : ) = transpose ( reshape ( list , (/ n , m /) ) )
                        else
                            do j = 1 , n
                                Jvalue ( i , j , : ) = list ( j )
                            end do
                        end if
                        if ( io == 0 ) i = i + 1
                    end do
                ! --- DMI ---
                else if ( term == 'DMI' ) then
                    read ( values , * , iostat = io ) DMI
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at DMI: bad input'
                        stop
                    end if
                else if ( term == 'D_value' ) then
                    if ( DMI .eqv. .False. ) cycle
                    read ( values , '(I4)' ) i
                    allocate ( Dvalue ( i , 8 ) )
                    Dvalue = 0.d0
                    i = 1
                    do while ( i <= size ( Dvalue , dim = 1 ) )
                        read ( 9 , '(A)' ) values
                        call remove_comment ( values )
                        read ( values , * , iostat = io ) Dvalue ( i , : )
                        if ( io == 0 ) then
                            i = i + 1
                        end if
                    end do

                ! --- single ion anisotropy ---
                else if ( term == 'SIA' ) then
                    read ( values , * , iostat = io ) SIA
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at SIA: bad input'
                        stop
                    end if
                else if ( term == 'A_value' ) then
                    if ( SIA .eqv. .False. ) cycle
                    read ( values , * , iostat = io ) Avalue
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at A_value: bad input'
                        stop
                    end if
                ! --- checking system stability in whole BZ ---
                else if ( term == 'BZ3D' ) then
                    read ( values , * , iostat = io ) bz3D
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at BZ3D: bad input'
                        stop
                    end if
                    if ( bz3D .eqv. .True. ) then
                        if ( sysdim /= 3 ) then
                            write ( * , * ) 'error occurs at BZ3D: dimension not match'
                            stop
                        end if
                        i = size ( modelist )
                        modelist ( 2 : i ) = modelist ( 1 : i - 1 )
                        modelist ( 1 ) = 0
                    end if
                else if ( term == 'BZ2D' ) then
                    read ( values , * , iostat = io ) bz2D
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at BZ2D: bad input'
                        stop
                    end if
                    if ( bz2D .eqv. .True. ) then
                        if ( sysdim /= 2 ) then
                            write ( * , * ) 'error occurs at BZ2D: dimension not match'
                            stop
                        end if
                        i = size ( modelist )
                        modelist ( 2 : i ) = modelist ( 1 : i - 1 )
                        modelist ( 1 ) = 0
                    end if
                ! --- bulk magnon calculation ---
                else if ( term == 'Bulk' ) then
                    read ( values , * , iostat = io ) bulk
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at Bulk: bad input'
                        stop
                    end if
                    if ( bulk .eqv. .True. ) then
                        do i = 1 , size ( modelist )
                            if ( modelist ( i ) == -1 ) then
                                modelist ( i ) = 1
                                exit
                            end if
                        end do
                    end if
                ! read in kpath for bulk1D/2D/3D
                else if ( term == 'kpath_bulk' ) then
                    if ( .not. bulk ) cycle
                    read ( values , '(I4)' ) knum
                    allocate ( bulkkpath ( knum , sysdim ) )
                    i = 1
                    do while ( i <= knum )
                        read ( 9 , '(A)' ) values
                        call remove_comment ( values )
                        read ( values , * , iostat = io ) bulkkpath ( i , : )
                        if ( io > 0 ) then
                            write ( * , * ) 'error occurs at kpath_bulk: bad input'
                            stop
                        else if ( io == 0 ) then
                            i = i + 1
                        end if
                    end do
                    bulkkpath = 2.0 * pi * bulkkpath
                ! --- slab calculation ---
                ! 2D surface state
                else if ( term == 'Slab2D' ) then
                    read ( values , * , iostat = io ) slab2D
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at Slab2D: bad input'
                        stop
                    end if
                    if ( slab2D .eqv. .True. ) then
                        if ( sysdim <= 2 ) then
                            write ( * , * ) 'error occurs at Slab2D: dimension not match'
                            stop
                        end if
                        do i = 1 , size ( modelist )
                            if ( modelist ( i ) == -1 ) then
                                modelist ( i ) = 2
                                exit
                            end if
                        end do
                    end if
                ! 2D surface state
                else if ( term == 'Slab1D' ) then
                    read ( values , * , iostat = io ) slab1D
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at Slab1D: bad input'
                        stop
                    end if
                    if ( slab1D .eqv. .True. ) then
                        if ( sysdim <= 1 ) then
                            write ( * , * ) 'error occurs at Slab1D: dimension not match'
                            stop
                        end if
                        do i = 1 , size ( modelist )
                            if ( modelist ( i ) == -1 ) then
                                modelist ( i ) = 2
                                exit
                            end if
                        end do
                    end if
                ! read in stacking direction for slab2D/1D
                else if ( term == 'stackdir_slab' ) then
                    if ( .not. any ( (/ slab1D , slab2D /) ) ) cycle
                    i = verify ( 'x' , values )
                    i = i + verify ( 'y' , values )
                    i = i + verify ( 'z' , values )
                    if ( i /= 2 ) then
                        write ( * , * ) 'error occurs at stackdir_slab: bad input'
                        stop
                    else
                        if ( verify ( 'x' , values ) == 0 ) then
                            slabdir = 1
                        else if ( verify ( 'y' , values ) == 0 ) then
                            slabdir = 2
                        else if ( verify ( 'z' , values ) == 0 ) then
                            slabdir = 3
                        end if
                    end if
                ! read in number of stacking layers for slab2D/1D
                else if ( term == 'nlayer' ) then
                    if ( .not. any ( (/ slab1D , slab2D /) ) ) cycle
                    read ( values , '(I4)' ) nlayer
                ! read kpath for slab2D/1D
                else if ( term == 'kpath_slab' ) then
                    if ( .not. any ( (/ slab1D , slab2D /) ) ) cycle
                    read ( values , '(I4)' ) knum
                    if ( slab2D .eqv. .True. ) then
                        allocate ( slabkpath ( knum , 2 ) )
                    else if ( slab1D .eqv. .True. ) then
                        allocate ( slabkpath ( knum , 1 ) )
                    end if
                    i = 1
                    do while ( i <= knum )
                        read ( 9 , '(A)' ) values
                        call remove_comment ( values )
                        read ( values , * , iostat = io ) slabkpath ( i , : )
                        if ( io /= 0 ) then
                            write ( * , * ) 'error occurs at kpath_slab: bad input'
                            stop
                        else if ( io == 0 ) then
                            i = i + 1
                        end if
                    end do
                    slabkpath = 2.0 * pi * slabkpath
                ! --- surface green function calculation ---
                ! 2D surface
                else if ( term == 'Surf2D' ) then
                    read ( values , * , iostat = io ) surf2D
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at Surf2D: bad input'
                        stop
                    end if
                    if ( surf2D .eqv. .True. ) then
                        if ( sysdim <= 2 ) then
                            write ( * , * ) 'error occurs at Surf2D: dimension not match'
                            stop
                        end if
                        do i = 1 , size ( modelist )
                            if ( modelist ( i ) == -1 ) then
                                modelist ( i ) = 3
                                exit
                            end if
                        end do
                    end if
                ! 1D surface
                else if ( term == 'Surf1D' ) then
                    read ( values , * , iostat = io ) surf1D
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at Surf1D: bad input'
                        stop
                    end if
                    if ( surf1D .eqv. .True. ) then
                        if ( sysdim <= 1 ) then
                            write ( * , * ) 'error occurs at Surf1D: dimension not match'
                            stop
                        end if
                        do i = 1 , size ( modelist )
                            if ( modelist ( i ) == -1 ) then
                                modelist ( i ) = 3
                                exit
                            end if
                        end do
                    end if
                ! read in stacking direction for surface green function
                else if ( term == 'stackdir_surf' ) then
                    if ( .not. any ( (/ surf2D , surf1D /) ) ) cycle
                    i = verify ( 'x' , values )
                    i = i + verify ( 'y' , values )
                    i = i + verify ( 'z' , values )
                    if ( i /= 2 ) then
                        write ( * , * ) 'error occurs at stackdir_surf: bad input'
                        stop
                    else
                        if ( verify ( 'x' , values ) == 0 ) then
                            surfdir = 1
                        else if ( verify ( 'y' , values ) == 0 ) then
                            surfdir = 2
                        else if ( verify ( 'z' , values ) == 0 ) then
                            surfdir = 3
                        end if
                    end if
                ! read in broadening factor for surface green function
                else if ( term == 'broaden' ) then
                    if ( .not. any ( (/ surf2D , surf1D /) ) ) cycle
                    read ( values , * ) bf
                ! kpath mode for surface green function
                else if ( term == 'kp_mode' ) then
                    if ( .not. any ( (/ surf2D , surf1D /) ) ) cycle
                    read ( values , * , iostat = io ) surfkp
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at kp_mode: bad input'
                        stop
                    end if
                ! read in energy range and step for gfunc kpath mode
                else if ( term == 'E_min' ) then
                    if ( ( .not. any ( (/ surf2D , surf1D /) ) ) .or. ( .not. surfkp ) ) cycle
                    read ( values , * ) wi
                else if ( term == 'E_max' ) then
                    if ( ( .not. any ( (/ surf2D , surf1D /) ) ) .or. ( .not. surfkp ) ) cycle
                    read ( values , * ) wf
                else if ( term == 'E_step' ) then
                    if ( ( .not. any ( (/ surf2D , surf1D /) ) ) .or. ( .not. surfkp ) ) cycle
                    read ( values , * ) ws
                ! read in kpath for gfunc kpath mode
                else if ( term == 'kpath_surf' ) then
                    if ( ( .not. any ( (/ surf2D , surf1D /) ) ) .or. ( .not. surfkp ) ) cycle
                    read ( values , '(I4)' ) knum
                    if ( surf2D .eqv. .True. ) then
                        allocate ( surfkpath ( knum , 2 ) )
                    else if ( surf1D .eqv. .True. ) then
                        allocate ( surfkpath ( knum , 1 ) )
                    end if
                    i = 1
                    do while ( i <= knum )
                        read ( 9 , '(A)' ) values
                        call remove_comment ( values )
                        read ( values , * , iostat = io ) surfkpath ( i , : )
                        if ( io > 0 ) then
                            write ( * , * ) 'error occurs at kpath_surf: bad input'
                            stop
                        else if ( io == 0 ) then
                            i = i + 1
                        end if
                    end do
                    surfkpath = 2.0 * pi * surfkpath
                ! equal energy mode for surface green function
                else if ( term == 'energy_mode' ) then
                    if ( .not. any ( (/ surf2D , surf1D /) ) ) cycle
                    read ( values , * , iostat = io ) surfe
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at energy_mode: bad input'
                        stop
                    end if
                else if ( term == 'E_cut' ) then
                    if ( ( .not. any ( (/ surf2D , surf1D /) ) ) .or. ( .not. surfe ) ) cycle
                    read ( values , * ) ecut
                ! --- topological property of system ---
                else if ( term == 'Topo' ) then
                    read ( values , * , iostat = io ) topo
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at Topo: bad input'
                        stop
                    end if
                    if ( topo .eqv. .True. ) then
                        do i = 1 , size ( modelist )
                            if ( modelist ( i ) == -1 ) then
                                modelist ( i ) = 4
                                exit
                            end if
                        end do
                    end if
                ! time reversal symmetric
                else if ( term == 'TRS' ) then
                    if ( .not. topo ) cycle
                    read ( values , * , iostat = io ) trs
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at TRS: bad input'
                        stop
                    end if
                ! read in number of bands to be included
                else if ( term == 'bands_include' ) then
                    if ( .not. topo ) cycle
                    i = index ( values , '-' )
                    read ( values ( : i - 1 ) , * , iostat = io ) bandi
                    if ( ( io /= 0 ) .or. ( bandi <= 0 ) ) then
                        write ( * , * ) 'error occurs at bands_include: bad input'
                        stop
                    end if
                    read ( values ( i + 1 : ) , * , iostat = io ) bandf
                    if ( ( io /= 0 ) .or. ( bandf <= 0 ) ) then
                        write ( * , * ) 'error occurs at bands_include: bad input'
                        stop
                    end if
                    bandn = bandf - bandi + 1
                    if ( bandn <= 0 ) then
                        write ( * , * ) 'error occurs at bands_include: bad input'
                        stop
                    end if
                ! 2D insulating system
                else if ( term == 'INSUL2D' ) then
                    if ( .not. topo ) cycle
                    read ( values , * , iostat = io ) insul2D
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at INSUL2D: bad input'
                        stop
                    end if
                    if ( ( insul2D .eqv. .True. ) .and. ( sysdim /= 2 ) ) then
                        write ( * , * ) 'error occurs at INSUL2D: dimension not match'
                        stop
                    end if
                ! topology defect
                else if ( term == 'SEMI1D' ) then
                    if ( .not. topo ) cycle
                    read ( values , * , iostat = io ) semi1D
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at SEMI1D: bad input'
                        stop
                    end if
                    if ( ( semi1D .eqv. .True. ) .and. ( sysdim <= 1 ) ) then
                        write ( * , * ) 'error occurs at SEMI1D: dimension not match'
                        stop
                    end if
                ! read in loop center
                else if ( term == 'loop_center' ) then
                    if ( .not. all ( (/ topo , semi1D /) ) ) cycle
                    allocate ( loopc ( sysdim ) )
                    read ( values , * , iostat = io ) loopc
                    if ( io /= 0 ) then
                        write ( * , * ) 'error occurs at loop_center: bad input'
                        stop
                    end if
                    !loopc = 2.0 * pi * loopc
                ! read in normal direction of the loop
                else if ( term == 'nor_dir' ) then
                    if ( .not. all ( (/ topo , semi1D /) ) ) cycle
                    if ( verify ( 'x' , values ) == 0 ) then
                        nordir = 1
                    else if ( verify ( 'y' , values ) == 0 ) then
                        nordir = 2
                    else if ( verify ( 'z' , values ) == 0 ) then
                        nordir = 3
                    else
                        write ( * , * ) 'error occurs at nor_dir: bad input'
                        stop
                    end if
                ! read in number of kpoints on the loop
                else if ( term == 'num_loopkp' ) then
                    if ( .not. all ( (/ topo , semi1D /) ) ) cycle
                    read ( values , * , iostat = io ) loopn
                    if ( ( io /= 0 ) .or. loopn <= 0 ) then
                        write ( * , * ) 'error occurs at num_loopkp: bad input'
                        stop
                    end if
                end if
            else if ( io < 0 ) then
                exit
            end if
        end do
        close ( 9 )

        return
    end subroutine read_info

    subroutine get_struct
        implicit none

        character ( 50 ) :: term
        integer ( 4 ) :: i , j , k , l
        integer ( 4 ) :: io , count , lines
        integer ( 4 ) , allocatable :: ineq_bond_num ( : ) ! store inequal number of each distance order
        real    ( 8 ) , allocatable :: list ( : , : , : )

        ! --- read result from output file ---
        open ( unit = 10 , file = 'neighbor_info.txt' , form = 'formatted' , access = 'sequential' )
        read ( unit = 10 , fmt = '(20A)' , iostat = io ) term
        ! find number of lines of neighbor info
        lines = 0
        do
            read ( unit = 10 , fmt = * , iostat = io )
            if ( io < 0 ) then
                exit
            else if ( io > 0 ) then
                write(*,*) 'error occurs when find number of lines of neighbor info'
                stop
            end if
            lines = lines + 1
        end do

        ! --- read in neighbor_list ---
        rewind ( 10 )
        ! skip sym group line
        read ( unit = 10 , fmt = * , iostat = io )
        ! --- consider inequivalent neighbors ---
        allocate ( neighbor_list ( lines , 7 ) )
        allocate ( ineq_bond_num ( Jorder ) )
        ineq_bond_num = 1

        do count = 1 , lines , 1
            ! (var1) n nearest neighbor , (var2) label for inequivalent part
            ! (var3) source atom , (var4) target atom
            ! (var5-7) cell difference = [ a , b , c ]
            read ( unit = 10 , fmt = * ) neighbor_list ( count , : )
            i = neighbor_list ( count , 1 )
            j = neighbor_list ( count , 2 )
            k = ineq_bond_num ( i )
            if ( j >= k ) ineq_bond_num ( i ) = j
        end do

        allocate ( list ( Jorder , size ( Jvalue , dim = 2 ) , 3 ) )
        list = Jvalue
        deallocate ( Jvalue )

        ! keep number of input Jvalue of each neighbor up to the inequivalent number
        ! if the inequivalent number is larger than the number of input J value
        ! then the remaining unfilled value will be set the same as the last value of that neighbor
        allocate ( Jvalue ( Jorder , maxval ( ineq_bond_num ) , 3 ) )
        Jvalue = 0
        i = size ( list , dim = 2 )
        do j = 1 , Jorder
            do l = 1 , 3
                k = ineq_bond_num ( j )
                if ( k <= i ) then
                    Jvalue ( j , : , l ) = list ( j , : k , l )
                else
                    Jvalue ( j , : i , l ) = list( j , : , l )
                    Jvalue ( j , i : k , l ) = list ( j , i , l )
                end if
            end do
        end do

        write ( unit = 97 , fmt = "(/,'Jvalue :')" )
        do i = 1 , Jorder
            do j = 1 , size ( Jvalue , dim = 2 )
                write ( unit = 97 , fmt = "(X,'J',I1,'_',I1,X,3F8.3)" ) i , j , Jvalue ( i , j , : )
            end do
        end do

        close ( 10 )
        return
    end subroutine get_struct

    subroutine get_cmat ( mode )
        implicit none
        integer ( 4 ) , intent ( in ) :: mode
        integer ( 4 ) :: i , dim

        ! prepare constant matrix required for each calculation mode
        if ( ( mode == 0 ) .or. ( mode == 1 ) .or. ( mode == 4 ) ) then
            allocate ( tz ( 2 * n_ineq , 2 * n_ineq ) )
            tz = v0
            do i = 1 , 2 * n_ineq , 1
                if ( i <= n_ineq ) then
                    tz ( i , i ) =  v1
                else
                    tz ( i , i ) = -v1
                end if
            end do

            if ( mode /= 4 ) return
            if ( trs .eqv. .False. ) return
            dim = 2 * n_ineq
            allocate ( tro ( dim , dim ) )
            tro = v0
            do i = 1 , dim / 4
                tro ( i , i + dim / 4 ) =  v1
                tro ( i + dim / 4 , i ) = -v1
            end do
            tro ( dim / 2 + 1 : , dim / 2 + 1 : ) = -tro ( : dim / 2 , : dim / 2 )
        else if ( mode == 2 ) then
            block
                integer ( 4 ) :: d1 , d2 , d3
                if ( FMConfig .eqv. .False. ) then
                    dim = 2 * n_ineq
                    allocate ( tz ( dim * nlayer , dim * nlayer ) )
                    allocate ( Id ( dim / 2 , dim / 2 ) )
                else
                    dim = n_ineq
                end if

                if ( allocated ( Id ) ) then
                    Id = v0
                    do i = 1 , dim / 2 , 1
                        Id ( i , i ) = v1
                    end do
                end if

                if ( allocated ( tz ) ) then
                    tz = v0
                    do i = 1 , nlayer , 1
                        d1 = dim * ( i - 1 )
                        d3 = dim * i
                        d2 = ( d1 + d3 ) / 2
                        tz ( d1 + 1 : d2 , d1 + 1 : d2 ) = Id
                        tz ( d2 + 1 : d3 , d2 + 1 : d3 ) = -1.d0 * Id
                    end do
                end if

            end block
        else if ( mode == 3 ) then
            if ( FMConfig .eqv. .False. ) then
                dim = 2 * n_ineq
                allocate ( tz ( dim , dim ) )
            else
                dim = n_ineq
            end if

            allocate ( Id ( dim , dim ) )
            Id = v0
            do i = 1 , dim , 1
                Id ( i , i ) = v1
            end do

            if ( allocated ( tz ) ) then
                tz = v0
                do i = 1 , dim , 1
                    if ( i <= dim / 2 ) then
                        tz ( i , i ) = v1
                    else
                        tz ( i , i ) = -v1
                    end if
                end do
            end if

        else
            write ( * , * ) 'constant matrix not prepared for current mode'
            stop
        end if
        return
    end subroutine get_cmat

    subroutine get_input
        implicit none
        logical :: fstat1 , fstat2
        ! --- check requied file ---
        inquire ( file = 'INFO_FORT' , exist = fstat1 )
        if ( fstat1 .eqv. .false. ) then
            write ( * , * ) 'INFO_FORT does not exist'
            stop
        end if

        ! --- get input ---
        call read_info
        call get_struct


        return
    end subroutine get_input

end module var
