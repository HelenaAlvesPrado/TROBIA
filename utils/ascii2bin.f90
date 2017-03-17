program ascii2bin

    implicit none
    
    CHARACTER *30 file_in, file_out
    CHARACTER *60 BUFFER 
    
    integer :: nx, ny, i, j
    
    parameter (nx=580, ny=9)   ! mudar ny para 360 linhas * 12 meses
    real*4, dimension(nx,ny):: arr_in
    
    CALL GETARG(1,BUFFER)
    READ(BUFFER,*) file_in
    CALL GETARG(2,BUFFER)
    READ(BUFFER,*) file_out
    
    open (unit=11,file=file_in,status='old',form='formatted',access='sequential',&
          action='read')
          
          
    open (unit=21,file=file_out,status='unknown',&
          form='unformatted',access='direct',recl=nx*ny*4)
          
          
     do j = 1, ny ! for each line do
             read(11,*) (arr_in(i,j), i=1,nx) ! read all elements in line j (implicit looping)
             !write(*,*) arr_in(:,j) 
     end do
     
     write(21,rec=1) arr_in 
  
    close(11)
    close(21)
    
end program ascii2bin
    
