
        implicit real *8 (a-h,o-z)
        integer nmax, m
        parameter (nterms=17)
c        
        real *8 x, y(0:nterms,0:nterms), d(0:nterms,0:nterms)
        dimension w(1 000 000)
c        
        lw = 1 000 000
        nmax = 10
c
c       SET ALL PARAMETERS
c       
        call prini(6,13)
c
c
c       Evaluate Legendre functions
c
        x = 0.0d0
        call prin2('x=*',x,1)
        call ylgndr(nterms, x, y)
        call prin2('after ylgndr, y=*',y,0)
        call prinm0(y,nterms)
c
        do i = 0,nterms
        write(14,*)(y(i,j),j=0,i)
        enddo
c
c
c       ... and their derivatives
c        
        x = 0.0d0
        call prin2('x=*',x,1)
        call ylgndr2(nterms, x, y, d)
        call prin2('after ylgndr2, y=*',y,0)
        call prinm0(y,nterms)
        call prin2('after ylgndr2, d=*',d,0)
        call prinm0(d,nterms)
c
        do i = 0,nterms
        write(14,*)(y(i,j),j=0,i)
        write(14,*)(d(i,j),j=0,i)
        enddo
c
c       
c       ... remove singularity at the poles
c        
        x = 1.0d0
        call prin2('x=*',x,1)
        call ylgndr2(nterms, x, y, d)
        call prin2('after ylgndr2, y=*',y,0)
        call prinm0(y,nterms)
        call prin2('after ylgndr2, d=*',d,0)
        call prinm0(d,nterms)
c
        do i = 0,nterms
        write(14,*)(y(i,j),j=0,i)
        write(14,*)(d(i,j),j=0,i)
        enddo
c
c
c
        x = 1.0d0
        call prin2('x=*',x,1)
        call ylgndr2s(nterms, x, y, d)
        call prin2('after ylgndr2s, y=*',y,0)
        call prinm0(y,nterms)
        call prin2('after ylgndr2s, d=*',d,0)
        call prinm0(d,nterms)
c
        do i = 0,nterms
        write(14,*)(y(i,j),j=0,i)
        write(14,*)(d(i,j),j=0,i)
        enddo
c
c        
c       ... and the symmetries are
c        
        x = 0.4d0
        call prin2('x=*',x,1)
        call ylgndr2(nterms, x, y, d)
        call prin2('after ylgndr2, y=*',y,0)
        call prinm0(y,nterms)
        call prin2('after ylgndr2, d=*',d,0)
        call prinm0(d,nterms)
c
        do i = 0,nterms
        write(14,*)(y(i,j),j=0,i)
        write(14,*)(d(i,j),j=0,i)
        enddo
c
        x = -x
        call prin2('x=*',x,1)
        call ylgndr2(nterms, x, y, d)
        call prin2('after ylgndr2, y=*',y,0)
        call prinm0(y,nterms)
        call prin2('after ylgndr2, d=*',d,0)
        call prinm0(d,nterms)
c
        do i = 0,nterms
        write(14,*)(y(i,j),j=0,i)
        write(14,*)(d(i,j),j=0,i)
        enddo
c
        call ylgndr2pm(nterms,y,d)
        call prin2('after ylgndr2pm, y=*',y,0)
        call prinm0(y,nterms)
        call prin2('after ylgndr2pm, d=*',d,0)
        call prinm0(d,nterms)
c
c
c       ... fast evaluation of legendre functions 
c
        call ylgndrfwini(nmax, w, lw, lused)
c
        call prinf('nmax=*',nmax,1)
        call prinf('nterms=*',nterms,1)
        call prinf('lw=*',lw,1)
        call prinf('lused=*',lused,1)
c
        x = 0.0d0
        call prin2('x=*',x,1)
        call ylgndrfw(nterms, x, y, w, nmax)
        call prin2('after ylgndr2fw, y=*',y,0)
        call prinm0(y,nterms)
c
        do i = 0,nterms
        write(14,*)(y(i,j),j=0,i)
        enddo
c
c
c       ... and their derivatives
c        
        x = 0.0d0
        call prin2('x=*',x,1)
        call ylgndr2fw(nterms, x, y, d, w, nmax)
        call prin2('after ylgndr2fw, y=*',y,0)
        call prinm0(y,nterms)
        call prin2('after ylgndr2fw, d=*',d,0)
        call prinm0(d,nterms)
c
        do i = 0,nterms
        write(14,*)(y(i,j),j=0,i)
        write(14,*)(d(i,j),j=0,i)
        enddo
c
c       
c       ... remove singularity at the poles
c        
        x = 1.0d0
        call prin2('x=*',x,1)
        call ylgndr2fw(nterms, x, y, d, w, nmax)
        call prin2('after ylgndr2fw, y=*',y,0)
        call prinm0(y,nterms)
        call prin2('after ylgndr2fw, d=*',d,0)
        call prinm0(d,nterms)
c
        do i = 0,nterms
        write(14,*)(y(i,j),j=0,i)
        write(14,*)(d(i,j),j=0,i)
        enddo
c
        x = 1.0d0
        call prin2('x=*',x,1)
        call ylgndr2sfw(nterms, x, y, d, w, nmax)
        call prin2('after ylgndr2sfw, y=*',y,0)
        call prinm0(y,nterms)
        call prin2('after ylgndr2sfw, d=*',d,0)
        call prinm0(d,nterms)
c
        do i = 0,nterms
        write(14,*)(y(i,j),j=0,i)
        write(14,*)(d(i,j),j=0,i)
        enddo
c
c
c

        stop
        end
c
c
c
c
c
