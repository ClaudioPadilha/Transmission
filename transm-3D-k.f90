module const
    implicit none
    integer, parameter :: n = 41 ! gridpoints on real space
    integer, parameter :: q = 10 ! factor multipliying the number of gridpoints after interpolation (linear so far)
    integer, parameter :: nv = 60 ! number of voltage points for i x v curve
    integer, parameter :: nk = 200 ! number of gridpoints for kz integration on the E sphere
    integer, parameter :: nomp = 16 ! number of threads for OMP execution
	real(8), parameter :: EF = 0.0 ! Fermi energy
	real(8), parameter :: dE = 1E-3 ! gridsize over energy
	real(8), parameter :: dV = 5E-2 ! voltage graph points for i x v curve
	real(8), parameter :: kBT = 2.5E-2 ! for temperture, units eV
	real(8), parameter :: l = 20.0 ! total device length, units nm
 	real(8), parameter :: me = 0.5685629E-29 ! electron mass, units eV * s^2 * nm^{-2}
	real(8), parameter :: hbar = 0.6582119E-15 ! Planck's constant, units eV * s
	real(8), parameter :: meff = 8.0 !8.51 ! electron's effective mass on TiO2 rutile
	real(8), parameter :: l1 = 0.0
	real(8), parameter :: l2 = 0.3
	real(8), parameter :: U = 0.0
	real(8), parameter :: VR = 0.0 ! right electrode voltage. Should be zero according to Hannes' model
	real(8), parameter :: qe = 1.60217733E-19 ! electron chage, Coulombs
	real(8), parameter :: hp = 4.135667516E-15 ! Planck's constant, eV * s
	real(8), parameter :: pi = 3.1415926
	real(8), parameter :: lx = 1.0D3 !electrode's x dimension in nm
	real(8), parameter :: ly = 1.0D3 !electrode's y dimension in nm
end module const

function pot(x)
    use const
    real(8), intent(in)  :: x
    real(8)              :: pot
    
    if ((x >= l1 * l) .and. (x <= l2 * l)) then
        pot = U
    else
        pot = 0.0
    endif
end function pot

subroutine interpolate (xa, ya, xi, yi)
    use const
    
    real(8), intent(in)     :: xa(n), ya(n)
    real(8), intent(out)    :: xi(n * q), yi(n * q)
    
    integer                 :: i, j
    real(8)                 :: dx
    
    do i = 1, n - 1
        dx = xa(i + 1) - xa(i)
        do j = 1, m
            xi((i - 1) * q + j) = (real(q + 1 - j) / real(q)) * xa(i) + (real(j - 1) / real(q)) * xa (i + 1) 
            yi((i - 1) * q + j) = (real(q + 1 - j) / real(q)) * ya(i) + (real(j - 1) / real(q)) * ya (i + 1)
        end do
    end do
end subroutine interpolate

function FD (E, V)
    use const
	real(8), intent(in)     :: E, V
	real(8)                 :: FD

    FD = 1.0D0 / (exp((E - (EF + V)) / kbT) + 1.0D0)
end function FD
 
function newk (E, V)
    use const
    real(8), intent(in)     :: E, V
    complex(8)              :: newk
    real(8)                 :: k
    
    k = 2 * me * (E - V) / hbar ** 2 ! 2me/hbar^2 units eV^{-1} * nm^{-2}
    
    if (k < 0) then
        newk = cmplx (0.0D0, sqrt(-k))
    else
        newk = cmplx (sqrt(k), 0.0D0)
    end if
end function newk

subroutine transmfunc(VL, E, P, T, R)
    use const
  	real(8), intent(in)                             :: VL, E, P(n)
  	real(8), intent(out)                            :: T, R

  	real(8)                                         :: m2hbar2 = 2 * me / hbar ** 2 ! 2me/hbar^2 units eV^{-1} * nm^{-2}
  	real(8)                                         :: dx, kLz, Tk
    integer                                         :: i, j, ikx, iky
    integer                                         :: nrhs = 1, lda = n, ldb = n, info ! lapack stuff
    integer, dimension(n)                           :: ipiv ! lapack stuff
    complex(8), dimension(n,n)                      :: H ! n x n matrix of system H x = b
    complex(8), dimension(n)                        :: b ! right hand side of system H x = b
    complex(8)                                      :: kL, kR, kRz, transm, reflec, newk, im
    
    dx = l / (n)
    
    im = cmplx (0.0D0, 1.0D0)
    
    kL = newk (E, VL) ! kL ** 2 = kLx ** 2 + kLy ** 2 + kLz ** 2 
    kR = newk (E, VR) ! same thing
    
    Tk = 0.0

    ! we only calculate all the stuff if we are sure that kL is real. Otherwise the incident wave decays exponetially.
    ! The wave should be able to survive the travel through the right electrode if its real part is big enough.
    ! Given that kR is either purelly real or purelly imaginary, this means that it will not decay exponetially
    ! as it goes through this electrode if this is the case, i. e. Re(kR) > 0, Im(kR) == 0.				
    if (abs(real(kL)) > 0.0 .and. abs(real(kR)) > 1.0E-9) then 
        do ikx = 1, floor(sqrt(2 * me * (lx ** 2) * (E - VL) / (pi * hbar) ** 2))
            do iky = 1, floor(sqrt(2 * me * (ly ** 2) * (E - VL) / (pi * hbar) ** 2 - (ikx * ly / lx) ** 2))
        
                ! kLz is just a part of kL. kLx and kLy are the quantum numbers of standing waves inside a waveguide
                ! of lx \timex ly dimensions
                kLz = sqrt(2 * me * (E - VL) / hbar ** 2 - (ikx * pi / lx) ** 2 - (iky * pi / ly) ** 2) 
                ! kR keeps the same kLx and kLy, as the momentum is conserved in these directions. So, we can
                ! get kRz the same way we did for kLz, just being carefull to use VR instead of VL
                kRz = sqrt(2 * me * (E - VR) / hbar ** 2 - (ikx * pi / lx) ** 2 - (iky * pi / ly) ** 2)
                
                ! Thou shall not divide by zero! (See where we get T below...)
                ! kRz should be strong (!?!) enough to survive the travel through the right electrode!
                if (kLz > 1.0E-9 .and. abs(real(kRz)) > 1.0E-9) then             
                    ! inicializar a matriz H e o vetor b zerados 
                    H = cmplx (0.0D0, 0.0D0)
                    b = cmplx (0.0D0, 0.0D0)
                    
                    ! Define scattering matrix */
                    do i = 2, n - 1
                    H(i,i) = cmplx (2.0 / dx ** 2 + m2hbar2 * meff * (P(i) - E), 0.0)
                        H(i,i - 1) = cmplx (-1.0 / dx ** 2, 0.0)
                        H(i,i + 1) = cmplx (-1.0 / dx ** 2, 0.0)
                    end do
                    ! borders
                    H(1,1) = -1.0 /dx + im * kLz
                    H(1,2) = cmplx (1.0 / dx, 0.0)
                    H(n,n) = 1.0 / dx - im * kRz
                    H(n,n-1) = cmplx (-1.0 / dx, 0.0)
                    b(1) = 2 * im * kLz
                    
                    call zgesv (n, nrhs, H, lda, ipiv, b, ldb, info)
                    
                    transm = b(n) * exp (-im * kRz * l)
                    reflec = b(1) - 1.0
                    
                    T = zabs(kRz)/abs(kLz)*zabs(transm)**2
                    R = zabs(reflec)**2
                else ! the incoming or the transmitted wave is fully absorbed at the electrodes... =(((((
                    T = 0.0
                    R = 0.0
                end if
                
                Tk = Tk + T ! integrate Tk = T(E, kLz): T(E) = \int T(E,kLz) dkLz
                
                !if (Tk > 0.0 ) then
                !	write (*, '(E12.4xE12.4xE12.4)') Tk, kLz, E
                !end if
            end do
        end do
    else
        T = 0.0
        R = 0.0
    end if

    T = Tk ! * dkLz <-- figure out how to calculate this thing!

end subroutine transmfunc

program transm
    use const
    
    integer                                 :: i, j, m
    real(8)                                 :: dx, pot, V, Vmin, Vmax, E, E0, En, FD, curr, lixo, VL
    real(8), dimension(:), allocatable      :: P, Pint, x, xint, T, R
    character (50)                          :: entrada, temp
    
    call get_command_argument (1, temp)
    call get_command_argument (2, entrada)
    
    read (temp, *) VL
   
    write (*, '(F5.2)') VL

    dx = l / (n)
    
    allocate (P(n))
    !allocate (Pint(n * q))
    allocate (x(n))
    !allocate (xint(n * q))
    
    !open(99,file='pot.dat')
    !do i = 1, n
    !	P(i) = pot (i * dx)
    !	write(99, '(F12.8xF12.8)') i * dx, P(i)
    !end do
    !close (99)
    
    open (99, file=trim(adjustl(entrada)), status='old', action='read')
    do i = 1, n
        read (99, *) x(i), P(i), lixo, lixo
    end do
    close (99)
    
    !call interpolate (x, P, xint, Pint)
    
    VL = -VL
    
    open(98, file='ixv.dat')
    
    if (VL < VR) then
        E0 = EF + VL - 25 * kBT
        En = EF + VR + 25 * kBT
    else
        E0 = EF + VR - 25 * kBT
        En = EF + VL + 25 * kBT
    end if
    
    m = (En - E0) / dE
    write (*, '(A,I8,A,F12.4,A,F12.4)') 'm = ', m, 'E0 = ', E0, 'En = ', En
    
    allocate (T(m))
    allocate (R(m))
    
    !open(99, file='transm.dat')
    
    E = E0
    curr = 0.0
    
    !call OMP_set_num_threads(nomp)
    !call OMP_set_thread_num(nomp)
    
    !$OMP DO
    do i = 1, m
        call transmfunc(VL, E, P, T(i), R(i))
        !write (*, *) OMP_get_thread_num()
        write (*, '(E12.4,E12.4)') E, T(i) !, R(i), FD(E,VL)-FD(E,VR)
        E = E + dE
        curr = curr + T(i) * (FD(E,VL)-FD(E,VR))
    end do
    !$OMP END DO
    !close (99)
    
    write (98, '(F7.4xE12.4)') VL, curr * dE  * 2 * qe / hp
    
    close (98)
    
    deallocate (P)
    !deallocate (Pint)
    deallocate (x)
    !deallocate (xint)
    deallocate (T)
    deallocate (R)    
end program transm

