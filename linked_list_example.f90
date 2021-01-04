PROGRAM test

  IMPLICIT NONE


  TYPE :: NODE2d
    
    INTEGER :: val_1d
    INTEGER :: val_2d(2)
    
    TYPE(NODE2d), POINTER :: next
  END TYPE NODE2d



  TYPE(NODE2d), POINTER :: box_test, head, p1, p2, p, p_prev, q_new, p_new

  INTEGER, ALLOCATABLE :: matrix(:)

  INTEGER :: lens, i, ss
  INTEGER :: Nx, Ny


  lens = 10
  Ny = 2

  ALLOCATE(matrix(Ny))
  matrix = (/1,10/)

  ! create first node
  ALLOCATE(p1)
  p1%val_1d = -1
  p1%val_2d = matrix
  NULLIFY(p1%next)
        
  head => p1


  ! create other nodes
  DO i=2,lens,1
    ALLOCATE(p2)
    p2%val_1d = -1*i
    p2%val_2d = matrix*i
    NULLIFY(p2%next)

    p1%next => p2
    p1 => p1%next

  ENDDO

  box_test => head

  WRITE(*,*) 11



  ! delete other node if node val is 5
  p => head
  p_prev => p
  i = 1
  DO WHILE(ASSOCIATED(p%next))
    p_prev => p
    p => p%next

    i = i+1
    WRITE(*,*)i
    IF(SUM(p%val_2d)==22.or.SUM(p%val_2d)==55)THEN
      p => p_prev%next
      p_prev%next => p%next
      DEALLOCATE(p)
      p => p_prev%next
    ENDIF

  ENDDO

  ! delete the node when val is 5
  p => head
  IF(SUM(p%val_2d)==22.or.SUM(p%val_2d)==55)THEN
    ! delete head node if node val is 5
    p => head
    head => p%next
    DEALLOCATE(p)
  ENDIF



  ! print out the linked list
  p => head
  DO WHILE(ASSOCIATED(p))
    WRITE(*,*) 'First:', p%val_1d, '*', p%val_2d
    p => p%next
  ENDDO


  ! insert new node in the head of linked list
  ALLOCATE(q_new)
  q_new%val_1d = -100
  q_new%val_2d = 100

  q_new%next => head
  head => q_new

  ! print out the linked list
  p => head
  DO WHILE(ASSOCIATED(p))
    WRITE(*,*)'Second', p%val_1d, '*', p%val_2d
    p => p%next
  ENDDO


END PROGRAM test
