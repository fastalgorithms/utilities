        subroutine peter_heapsort_rea(ra,ia,n)
        implicit real *8 (a-h,o-z)
        DIMENSION ia(n),ra(1)
c
c        This subroutine uses the standard heapsort scheme to
c        sort a real array ra. The sort is returned in the 
c        form of the integer array ia, such that the values
c
c        ra(ia(1)), ra(ia(2)), ra(ia(3)), . . .  ra(ia(n))        (1)
c
c        are in increasing order.
c
c            Input parameters:
c
c  ra - the array to be sorted
c  n - the number of elements in array ra
c
c            Output parameters:
c
c  ia - integer array containing the permutation that sorts ra 
c        (see (1) above)
c
c
        do 1200 i=1,n
c
        ia(i)=i
 1200 continue
c
        do 1600 ii=n,2,-2
c
        ii0=ii/2
        call peter_heapit_rea(ra,ia,n,ii0)
 1600 continue
c
c       . . . sort
c
        do 2000 i=1,n-1
c
        m=n-i+1
c
        jj=ia(m)
        ia(m)=ia(1)
        ia(1)=jj
c
        m=m-1
        call peter_heapit_rea(ra,ia,m,1)
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine peter_heapit_rea(ra,ia,n,ii)
        implicit real *8 (a-h,o-z)
        DIMENSION ia(n),ra(1)
c
c       heapify
c
        i=ii
c        
        do 1400 ijk=1,100
c
        ison1=i*2
        if(ison1 .gt. n) return
c
        ison2=i*2+1
        if(ison2 .gt. n) then
c
            j0=ia(i)
            j1=ia(ison1)
            if(ra(j0) .ge. ra(j1)) return
c
            ia(i)=j1
            ia(ison1)=j0
            return
         endif
c
        j0=ia(i)
        j1=ia(ison1)
        j2=ia(ison2)
c
        if( (ra(j0) .ge. ra(j1)) .and. (ra(j0) .ge. ra(j2))) return
c
        if(ra(j1) .ge. ra(j2)) then
            ia(i)=j1
            ia(ison1)=j0
            i=ison1
            goto 1400
        endif
c
            ia(i)=j2
            ia(ison2)=j0
            i=ison2
c
 1400 continue
            return
            end
