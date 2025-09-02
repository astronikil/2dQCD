            
**************************************************************************
c SRANF: scalar random number generator
c
c adaptation of knuth's algorithm, with throw-away of 910 rnds, HN, Dec. 8, 2003
c*************************************************************************
       DOUBLE PRECISION FUNCTION sranf ()
       IMPLICIT DOUBLE PRECISION (A)
       DIMENSION A(1009)
       INTEGER counter
c       EXTERNAL RNFARR, RNFSTR
       EXTERNAL RNFARR
       common /vsranfom1/ a
       common /vsranfom2/ counter
       save /vsranfom1/ 
       save /vsranfom2/ 

c      initialization must set counter to 100

       if(counter.lt.100) then
       sranf=a(counter)
       else
       counter=counter-99
       CALL RNFARR(A,1009)
       sranf=a(counter)
       end if       
       counter=counter+1
       return
       end


c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------

      SUBROUTINE RNFARR(AA,N)
C       FORTRAN 77 version of "ranf_array"
C       from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
C       including the MODIFICATIONS made in the 9th printing (2002)
C       ********* see the book for explanations and caveats! *********
      IMPLICIT DOUBLE PRECISION (A,R,Y)
      DIMENSION AA(*)
      PARAMETER (KK=100)
      PARAMETER (LL=37)
      COMMON /RSTATE/ RANX(KK)
      SAVE /RSTATE/
      DO 1 J=1,KK
      AA(J)=RANX(J)
 1    CONTINUE
      DO 2 J=KK+1,N
         Y=AA(J-KK)+AA(J-LL)
         AA(J)=Y-IDINT(Y)
 2    CONTINUE
      DO 3 J=1,LL
         Y=AA(N+J-KK)+AA(N+J-LL)
         RANX(J)=Y-IDINT(Y)
 3    CONTINUE
      DO 4 J=LL+1,KK
         Y=AA(N+J-KK)+RANX(J-LL)
         RANX(J)=Y-IDINT(Y)
 4    CONTINUE
      END

c-----------------------------------------------------
C Initilaization subroutine.
c----------------------------------------------------
      SUBROUTINE setrn (SEED)
      IMPLICIT DOUBLE PRECISION (A,R,U,V)
      IMPLICIT INTEGER (T)
      PARAMETER (KK=100)
      PARAMETER (LL=37)
      PARAMETER (MM=2**30)
      PARAMETER (ULP=1D0/(2D0**52))
      PARAMETER (TT=70)
      PARAMETER (KKK=KK+KK-1)
      INTEGER SEED,S,SSEED
      INTEGER counter
      DOUBLE PRECISION SS
      DIMENSION U(KKK) 
      common /vsranfom2/ counter
      COMMON /RSTATE/ RANX(KK)
      SAVE /RSTATE/
      save /vsranfom2/
      counter=100
      IF (SEED .LT. 0) THEN
         SSEED=MM-1-MOD(-1-SEED,MM)
      ELSE
         SSEED=MOD(SEED,MM)
      END IF
      SS=2D0*ULP*(SSEED+2)
      DO 1 J=1,KK
         U(J)=SS
         SS=SS+SS
         IF (SS .GE. 1D0) SS=SS-1D0+2*ULP
 1    CONTINUE
      U(2)=U(2)+ULP
      S=SSEED
      T=TT-1
 10   DO 12 J=KK,2,-1
         U(J+J-1)=U(J)
         U(J+J-2)=0
 12   CONTINUE
      DO 14 J=KKK,KK+1,-1
         V=U(J-(KK-LL))+U(J)
         U(J-(KK-LL))=V-IDINT(V)
         V=U(J-KK)+U(J)
         U(J-KK)=V-IDINT(V)
 14   CONTINUE
      IF (MOD(S,2) .EQ. 1) THEN
         DO 16 J=KK,1,-1
            U(J+1)=U(J)
 16      CONTINUE
         U(1)=U(KK+1)
         V=U(LL+1)+U(KK+1)
         U(LL+1)=V-IDINT(V)
      END IF
      IF (S .NE. 0) THEN
         S=S/2
      ELSE
         T=T-1
      END IF
      IF (T .GT. 0) GO TO 10
      DO 20 J=1,LL
        RANX(J+KK-LL)=U(J)
 20   CONTINUE
      DO 21 J=LL+1,KK
         RANX(J-LL)=U(J)
 21   CONTINUE
      DO 22 J=1,10
         CALL RNFARR(U,KKK)
 22   CONTINUE
      END

      

