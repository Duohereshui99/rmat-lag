ccccccc 
!this subroutine is used to calculate the inverse of a matrix
!with the subroutine of zgetrf and zgetri
!which comes from the lapack library
!for an arbitrary complex square matrix
ccccccc
            subroutine mat_inv(a,lda,n)
            integer::n,lda,info,lwork
            complex*16::a(n,n)
            complex*16,allocatable::work(:)  
            integer,allocatable::ipiv(:)

            lwork=n*n
            allocate(work(lwork),ipiv(n))

ccccccc
            call zgetrf(n, n, a, lda, ipiv, info)
ccccccc
            if (info .ne. 0) then
                print *, 'ZGETRF failed, info = ', info
                stop
            end if
ccccccc
            call zgetri(n, a, lda, ipiv, work, lwork, info)
ccccccc
            if (info .ne. 0) then
                print *, 'ZGETRI failed, info = ', info
                stop
            end if
ccccccc
            deallocate(work,ipiv)
ccccccc
            end subroutine mat_inv