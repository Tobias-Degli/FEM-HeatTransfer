!gfortran -I/usr/local/include/fgsl programa.f90 -lfgsl -lm ; ./a.out ; gnuplot plot.plt
program diferencas_finitas
  use fgsl
  implicit none
  integer(fgsl_size_t):: n,i,j,nx,ctd
  integer(fgsl_int) :: status, signum
  type(fgsl_matrix) :: M_a
  type(fgsl_vector) :: V_b, V_x
  real(fgsl_double), dimension (:), allocatable::Vet_T,Vet_B
  real(fgsl_double), dimension (:,:), allocatable::Mat_T
  type(fgsl_permutation) :: p
  real(8)::dx,dz,Lx,Lz,T_direita,T_inf,k,q,h,t0,t1,t2
  open(2,file="dados.txt",status='unknown')

  nx=25
  n=nx**2
  allocate(Mat_T(n,n),Vet_T(n),Vet_B(n))

  !           DEFININDO PARÂMETROS
  Lx=6
  Lz=3
  dx=Lx/(nx-1)
  dz=Lz/(nx-1)
  k=60
  q=20
  h=350
  T_direita=300
  T_inf=310

  call cpu_time(t0)
!        DISCRETIZANDO O DIMÍNIO - MONTANDO A MATRIZ A E VETOR B
! ROTINA A (Lateral direita)
  do i=1,nx
    !print*,'i=',i*nx
    Mat_T(i*nx,i*nx)  = 1
    Vet_B(i*nx)  = T_direita
  end do

! ROTINA B (Linha inferior exceto as extremidades)
  do i=nx**2-nx+2,nx*nx-1
    !print*,'i=',i
    Mat_T(i,i)  = -2/dx**2  - 1/dz**2
    Mat_T(i,i-1)  = 1/dx**2
    Mat_T(i,i+1)  = 1/dx**2
    Mat_T(i,i-nx) = 1/dz**2
    Vet_B(i)  = 0
  end do

!VÉRTICE C (vertice superior esquerdo)
  i=1
  !print*,'i=',i
  Mat_T(i,i)  =  -1/dx**2 - 1/dz**2 !-(h*dz+1)
  Mat_T(i,i+1)  = 1/dx**2
  Mat_T(i,i+nx) = 1/dz**2
  Vet_B(i)  = -  q/k !- h*dz*T_inf

!ROTINA D (Lateral esquerda, exceto os vértices)
  do i=2,nx-1
    !print*,'i=',(i-1)*nx+1
    Mat_T((i-1)*nx+1,(i-1)*nx+1)  = -1/dx**2 - 2/dz**2
    Mat_T((i-1)*nx+1,(i-1)*nx+2)  = 1/dx**2
    Mat_T((i-1)*nx+1,(i-1)*nx+1-nx) = 1/dz**2
    Mat_T((i-1)*nx+1,(i-1)*nx+1+nx) = 1/dz**2
    Vet_B((i-1)*nx+1)  = -q/k
  end do

! VÉRTICE E (Vértice inferior esquerdo)
  i=nx**2-nx+1
  !print*,'i=',i
  Mat_T(i,i)  = -1/dx**2  - 1/dz**2
  Mat_T(i,i-nx) = 1/dz**2
  Mat_T(i,i+1)  = 1/dx**2
  Vet_B(i)  = -q/k

! ROTINA F (Linha superior, exceto os vértices)
  do i=2,nx-1
    !print*,'i=',i
    Mat_T(i,i)  = -2/dx**2  -  1/dz**2 -(h*dz+1)
    Mat_T(i,i-1)  = 1/dx**2
    Mat_T(i,i+1)  = 1/dx**2
    Mat_T(i,i+nx) = 1/dz**2
    Vet_B(i)  =  - h*dz*T_inf
  end do

! ROTINA G (Interior da malha Toda a malha, exceto o contorno)
  i=nx+1
  do ctd=1,nx-2
      do i = ctd*nx+2,(ctd+1)*nx-1
        !print*,'i=',i
        Mat_T(i,i)  = -2/dx**2 -2/dz**2
        Mat_T(i,i-1)  = 1/dx**2
        Mat_T(i,i+1)  = 1/dx**2
        Mat_T(i,i+nx) = 1/dz**2
        Mat_T(i,i-nx) = 1/dz**2
        Vet_B(i)=0
      end do
  end do
! IMPRIMINDO A MATRIZ
  !write(*,*)'MATRIZ'
  !do i=1 ,n
  !  write(*,*)Mat_T(i,:),Vet_B(i)
  !end do
  call cpu_time(t1)
! RESOLVENDO O SISTEMA
  Mat_T = transpose(Mat_T)
  M_a = fgsl_matrix_init(type=1.0_fgsl_double)
  V_b = fgsl_vector_init(type=1.0_fgsl_double)
  V_x = fgsl_vector_init(type=1.0_fgsl_double)
  p = fgsl_permutation_alloc(n)
  status = fgsl_matrix_align(Mat_T, n, n, n, M_a)
  status = fgsl_vector_align(Vet_B, n, V_b, n, 0_fgsl_size_t, 1_fgsl_size_t)
  status = fgsl_vector_align(Vet_T, n, V_x, n, 0_fgsl_size_t, 1_fgsl_size_t)
  status = fgsl_linalg_LU_decomp (M_a, p, signum)
  status = fgsl_linalg_LU_solve (M_a, p, V_b, V_x)
  call fgsl_matrix_free(M_a)
  call fgsl_vector_free(V_b)
  call fgsl_vector_free(V_x)
  call fgsl_permutation_free(p)
  call cpu_time(t2)
  print*,'Tempo para definir a matriz',t1-t0
  print*,'Tempo para resolver o sistema',t2-t1
  !print*,'VETOR SOLUÇÃO'
  !print*,Vet_T

  !VINCULANDO O VETOR SOLUÇÃO COM A POSIÇÃO DOS EIXOS CARTEZIANOS
  ctd  = 0
  do j=1,nx
    write(2,*) !Linha em branco para o Gnuplot
    do i=1,nx
      ctd  = ctd+1
      write(2,*)(i-1)*dx,(nx-j)*dz,Vet_T(ctd)
      !write(*,*)(j-1)*dx,(nx-i)*dz,Vet_T(ctd)
    end do
  end do
  close(2)
  deallocate(Mat_T,Vet_T,Vet_B)
end program diferencas_finitas
