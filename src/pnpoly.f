cccc
c     Original Version: Randolph Franklin
c     Modifications (minor): Travis Askham
cccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Copyright notice from original version
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Copyright (c) 1970-2003, Wm. Randolph Franklin

c     Permission is hereby granted, free of charge, to any person 
c     obtaining a copy of this software and associated documentation 
c     files (the "Software"), to deal in the Software without restriction, 
c     including without limitation the rights to use, copy, modify, merge, 
c     publish, distribute, sublicense, and/or sell copies of the Software, 
c     and to permit persons to whom the Software is furnished to do so, 
c     subject to the following conditions:

c     Redistributions of source code must retain the above copyright 
c     notice, this list of conditions and the following disclaimers.

c     Redistributions in binary form must reproduce the above copyright 
c     notice in the documentation and/or other materials provided with 
c     the distribution.

c     The name of W. Randolph Franklin may not be used to endorse or 
c     promote products derived from this Software without specific prior
c     written permission.

c     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
c     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
c     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
c     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
c     BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
c     ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
c     CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
c     SOFTWARE.


C>>>PNP1                                                                
C                                                                       
C     ..................................................................
C                                                                       
C        SUBROUTINE PNPOLY                                              
C                                                                       
C        PURPOSE                                                        
C           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
C                                                                       
C        USAGE                                                          
C           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
C                                                                       
C        DESCRIPTION OF THE PARAMETERS                                  
C           PX      - X-COORDINATE OF POINT IN QUESTION.                
C           PY      - Y-COORDINATE OF POINT IN QUESTION.                
C           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
C                     VERTICES OF POLYGON.                              
C           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
C                     VERTICES OF POLYGON.                              
C           N       - NUMBER OF VERTICES IN THE POLYGON.                
C           INOUT   - THE SIGNAL RETURNED:                              
C                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
C                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
C                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
C                                                                       
C        REMARKS                                                        
C           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
C           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
C           OPTIONALLY BE INCREASED BY 1.                               
C           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
C           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
C           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
C           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
C           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
C           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
C           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
C                                                                       
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
C           NONE                                                        
C                                                                       
C        METHOD                                                         
C           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
C           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
C           POINT IS INSIDE OF THE POLYGON.                             
C                                                                       
C     ..................................................................
C                                                                       
      SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT)      
      IMPLICIT NONE
      REAL*8 XX(*),YY(*), PX, PY
      REAL*8, ALLOCATABLE :: X(:), Y(:)
      LOGICAL MX,MY,NX,NY                                               
      INTEGER O, I, J, N, INOUT                                                         
C      OUTPUT UNIT FOR PRINTED MESSAGES                                 
      DATA O/6/                                                         
c      MAXDIM=200                                                        
c      IF(N.LE.MAXDIM)GO TO 6                                            
c      WRITE(O,7)                                                        
c7     FORMAT('0WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY.     
c     1RESULTS INVALID')                                                 
c      RETURN

      ALLOCATE(X(N))
      ALLOCATE(Y(N))
                                                            
6     DO 1 I=1,N                                                        
      X(I)=XX(I)-PX                                                     
1     Y(I)=YY(I)-PY                                                     
      INOUT=-1                                                          
      DO 2 I=1,N                                                        
      J=1+MOD(I,N)                                                      
      MX=X(I).GE.0.0                                                    
      NX=X(J).GE.0.0                                                    
      MY=Y(I).GE.0.0                                                    
      NY=Y(J).GE.0.0                                                    
      IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2       
      IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3  
      INOUT=-INOUT                                                      
      GO TO 2                                                           
3     IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5                       
4     INOUT=0                                                           
      RETURN                                                            
5     INOUT=-INOUT                                                      
2     CONTINUE                                                          
      RETURN                                                            
      END                                                               



      SUBROUTINE PNPOLY2(PX,PY,XX,YY,X,Y,N,INOUT)      
      IMPLICIT NONE
      REAL*8 XX(*),YY(*), PX, PY, X(*), Y(*)
      LOGICAL MX,MY,NX,NY                                               
      INTEGER O, I, J, N, INOUT                                                         
C      OUTPUT UNIT FOR PRINTED MESSAGES                                 
      DATA O/6/                                                         
c      MAXDIM=200                                                        
c      IF(N.LE.MAXDIM)GO TO 6                                            
c      WRITE(O,7)                                                        
c7     FORMAT('0WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY.     
c     1RESULTS INVALID')                                                 
c      RETURN
                
ccC$OMP PARALLEL DO PRIVATE(I)                                        
      DO I=1,N                                                        
         X(I)=XX(I)-PX                                                     
         Y(I)=YY(I)-PY                                                     
      ENDDO
ccC$OMP END PARALLEL DO

      INOUT=1 

ccC$OMP PARALLEL DO PRIVATE(I,J,MX,NX,MY,NY), REDUCTION(*:INOUT)
      DO 2 I=1,N                                                        
      J=1+MOD(I,N)                                                      
      MX=X(I).GE.0.0                                                    
      NX=X(J).GE.0.0                                                    
      MY=Y(I).GE.0.0                                                    
      NY=Y(J).GE.0.0                                                    
      IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2       
      IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3  
      INOUT=-INOUT                                                      
      GO TO 2                                                           
3     IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5                       
4     INOUT=0                                                           
5     INOUT=-INOUT                                                      
2     CONTINUE 
ccC$OMP END PARALLEL DO   
      INOUT= -INOUT
      RETURN                                                            
      END                                                               
