C     MARKO-SIGGIA DUMBBELLS WITH 
C     INTERNAL VISCOSITY
c     + EXCLUDED VOLUME INTERACTIONS 
C     + HYDRODYNAMIC INTERACTIONS
C     SECOND ORDER SEMI-IMPLICIT PREDICTOR-CORRECTOR ALGORITHM
c
c     NOTE : WHEN RUNNING EQUILIBRIUM SIMULATIONS, DIVERT PROGRAM
c     FLOW FROM ENTERING TEXTRA. WHEN SR=0.0, MATERIAL FUNCTIONS BECOME UNDEFINED.
C     THIS CREATES HAVOC WITH TEXTRA
C
c     3-APR-2018  : DISCRIMINANT OF CUBIC EQUATIONS FOR  
c     MARKO-SIGGIA CASE IS POSITIVE, CANNOT USE TRIG. SOLVER
c     HAVE IMPLEMENTED NEWTON-RAPHSON+POLISHING(AS DONE IN
C     THE SINGLE CHAIN CODE)     

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NDATM=50,NOUT=100) 
      PARAMETER (NBINS=50)
      REAL*8 FENFAC,TEMPB
      REAL*8 AQ2,MQ2,VQ2
      REAL*8 QMAG,QLEN,TEMP12,TEMP23
      REAL*8 AVQ2(NOUT),MEQ2,ERQ2(NOUT)
      REAL*8 DTARR(10) 
      REAL*8 XARR(NDATM),YARR(NDATM),SIGARR(NDATM) 
      REAL*8 TEMP1(NDATM, NOUT),TEMP2(NDATM, NOUT),TEMP3(NDATM, NOUT)
      REAL*8 OUTIME(NOUT)
      REAL*8 Q(3)
      REAL*8 EQBFAC
      INTEGER NTIARR(10) 
      INTEGER TIME (8)
      CHARACTER (LEN=12) CLK(3)
      REAL*8 KAPPA_TEMP(3,3),FLOW_TMP(3,3)
      COMMON /STEPBL/ THI,B,ZMU,RMU2,BAUXQ,BAUXR,DTH,DTQ,SRDT,
     &SRDTH,C1P,C2P,E,AMPL2,BETA,AB,A2,A4,AUX1,AUX2,
     &AUX3,AUX4,AUX5,S,QALPH,GEE1,GEE2,GEE3,GEE4,FKAPPA(3,3),COEFF(4) 
      COMMON /EXTRP/ NOPT,NDUOPT,XOPT,YOPT,SIGOPT,ALIOPT,VLIOPT 

c      INTEGER :: THI
C     Excluded-volume and FENE parameters 
c     Z = Solvent Quality, RMU = EV parameter, B= FENE parameter
c     SR = Shear Rate, E = IV PARAMETER,
C     H0 = HYDRODYNAMIC INTERACTION PARAMETER
C     THI = EXPRESSION TO 
C     BE USED OFR THE HYDRODYNAMIC INTERACTION TENSOR.
C     THI = 1 FOR REGULARIZED OSEEN-BURGER
C     THI = 2 FOR ROTNE-PRAGER-YAMAKAWA 
C     NFLOW DECIDES TYPE OF FLOW TENSOR
C    (1) FOR STEADY SHEAR
C    (2) FOR PLANAR ELONGATION
C    (3) FOR UNIAXIAL ELONGATION

      open (unit=1, file='inp.dat') 
      READ (1,*) Z,RMU,B,SR,E,H0,THI,NFLOW 
      open (unit=2, file='tstep.dat') 
      open (unit=7, file='q2.dat', STATUS='UNKNOWN')
      open (unit=113, file='eqbconfigs.dat',STATUS='UNKNOWN')
      TMAX=10.D0
      READ (2,*) NTIWID, NTRAJ
      DO 15 I = 1, NTIWID
          READ (2,*) NTIARR(I)
          DTARR(I)=TMAX/NTIARR(I)
15    CONTINUE

C 
      ZMU=Z/RMU**(5.D0) 
      RMU2=(2.D0)*RMU*RMU
      SRTEMP=SR
      NTHI=NINT(THI) 
C 
      SELECT CASE (NTHI)
          CASE (1)
              WRITE(*,*) "THI : ",THI
              WRITE(*,*) "TYPE OF HI : ROBT"
              IF(H0.EQ.0) THEN
                  WRITE(*,*) "H0=0; SETTING THI=3"
                  THI=3.D0
              ENDIF
          CASE (2)
              WRITE(*,*) "THI : ",THI
              WRITE(*,*) "TYPE OF HI : RPY"
              IF(H0.EQ.0) THEN
                  WRITE(*,*) "H0=0; SETTING THI=3"
                  THI=3.D0
              ENDIF
          CASE DEFAULT
              WRITE(*,*) "HI OPTION NOT VALID"
              WRITE(*,*) "SETTING H0=0.0"
              H0=0.D0
      END SELECT

      CALL CPU_TIME(STARTTIME)


C    DECIDES ON THE TYPE OF VELOCITY FIELD TO BE IMPOSED ON SYSTEM : 
C    (1) FOR STEADY SHEAR
C    (2) FOR PLANAR ELONGATION
C    (3) FOR UNIAXIAL ELONGATION


      DO 75 I=1,3
          DO 76 J=1,3
              FKAPPA(I,J)=0.D0
              FLOW_TMP(I,J)=0.D0
76        CONTINUE
75    CONTINUE

      SELECT CASE (NFLOW)
          CASE (1)
              WRITE(*,*) "NFLOW : ",NFLOW
              WRITE(*,*) "TYPE OF FLOW : STEADY SHEAR"
              FLOW_TMP(1,2)=1.D0
          CASE (2)
              WRITE(*,*) "NFLOW : ",NFLOW
              WRITE(*,*) "TYPE OF FLOW : PLANAR ELONGATION"
              FLOW_TMP(1,1)=1.D0
              FLOW_TMP(2,2)=-1.D0
          CASE (3)
              WRITE(*,*) "NFLOW : ",NFLOW
              WRITE(*,*) "TYPE OF FLOW : UNIAXIAL ELONGATION"
              FLOW_TMP(1,1)=1.D0
              FLOW_TMP(2,2)=-0.5D0
              FLOW_TMP(3,3)=-0.5D0
          CASE DEFAULT
              WRITE(*,*) "FLOW OPTION NOT VALID"
              STOP "USE NFLOW=1 FOR SHEAR FLOW, NFLOW=2 FOR PLANAR ELONG
     &ATION OR NFLOW=3 FOR UNIAXIAL ELONGATION"
      END SELECT

      PI=3.1415926535897931D0
      AB=SQRT(PI)*H0
      A2=AB*AB
      A4=A2*A2
45    FORMAT(F20.16)
      WRITE(*,45) E
      WRITE(*,45) H0

   
c      CALL CPU_TIME(STARTTIME)
      DO 1000 IDT=1,NTIWID
C        Auxiliary parameters 
c         OPEN(UNIT=114,file=FPATH)
         WRITE(*,*) "DOING TIMESTEP WIDTH : ",DTARR(IDT)
         DELTAT=DTARR(IDT) 
         TEMPDT=DELTAT
         NTIME=NTIARR(IDT)/NOUT
         DTH=0.5D0*DELTAT 
         DTQ=0.25D0*DELTAT 
         SQDT=SQRT(DELTAT) 
         C1P=14.14855378D0*SQDT 
         C2P=1.21569221D0*SQDT 
         

C     Initializing averages and errors... 
C     AV... ARE AVERAGES OF CORRESPONDING MATL. FNS. OVER
C     ALL THE BLOCKS; ER... ARE STANDARD ERRORS FOR THE CORR-
C     -ESPONDING MATL. FNS.
         DO 25 I=1,NOUT
             AVQ2(I)=0.D0
             ERQ2(I)=0.D0
25       CONTINUE

C        A FRESH SEED IS GIVEN FOR EACH TIME-STEP WIDTH
C        TO ENSURE NON-OCCURENCE OF PERIOD EXHAUSTION

c         CALL DATE_AND_TIME(CLK(1),CLK(2),CLK(3),TIME)
c         ISEED=TIME(8)*100000+TIME(7)*1000+TIME(6)*10+TIME(5)
c         ISEED=ISEED+13998         
         ISEED=2231158
         CALL RANILS(ISEED)

         DO 100 ITRAJ=1,NTRAJ 
        
             Q(1)=RANGLS()
             Q(2)=RANGLS()
             Q(3)=RANGLS()

             IF(MODULO(ITRAJ,1000).EQ.0)THEN
             WRITE(*,*) "STATUS : EQB.TIME-STEP WIDTH : ",DELTAT,
     &"TRAJ # ",ITRAJ
             ENDIF
c             WRITE(*,*) "EQUILIBRATION AT SR = ",SR 
C            Relaxation of initial conditions
             NEQB=20.D0/DELTAT
             DO 78 I=1,3
                 DO 79 J=1,3
                     FKAPPA(I,J)=0.D0
79               CONTINUE
78           CONTINUE

             DO 50 ITIME=1,NEQB
                CALL SEMIMP(Q)
50           CONTINUE
C 
             DO 95 I=1,3
                 DO 96 J=1,3
                     FKAPPA(I,J)=FLOW_TMP(I,J)*SR
96               CONTINUE
95           CONTINUE

c             WRITE(*,*) "PRODUCTION AT SR = ",SR 
             IWAIT=0
             IOUT=0
C
C           Time integration: semi-implicit predictor-corrector scheme 
             DO 10 ITIME=1,NTIARR(IDT) 
                 CALL SEMIMP(Q) 
                 IWAIT=IWAIT+1
                 IF (IWAIT.EQ.NTIME) THEN
                     IWAIT=0
                     IOUT=IOUT+1
                     QMAG =Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3)
                     AVQ2(IOUT)=AVQ2(IOUT)+(QMAG)
                     ERQ2(IOUT)=ERQ2(IOUT)+(QMAG*QMAG)
                 ENDIF
10           CONTINUE 

             IF(IDT.EQ.NTIWID)THEN
                 WRITE(113,8) ITRAJ,Q(1),Q(2),Q(3)
             ENDIF
8            FORMAT(I10,4X,F20.16,4X,F20.16,4X,F20.16)

100      CONTINUE 
C 
C        Averages, statistical errors 

         DO 35 I=1,NOUT
             
             AVQ2(I)=AVQ2(I)/NTRAJ
             ERQ2(I)=ERQ2(I)/NTRAJ
             ERQ2(I)=(ERQ2(I)-AVQ2(I)*AVQ2(I))/(NTRAJ-1)
             ERQ2(I)=SQRT(ERQ2(I))

             OUTIME1=NTIME*I*DELTAT
c        Output of results 
             WRITE(7,4)  DELTAT,OUTIME1,AVQ2(I),ERQ2(I)
35       CONTINUE
c1        FORMAT(F11.8,4X,F6.2,2X,F16.5,2X,F16.5)
4        FORMAT(F11.8,4X,F10.5,4X,F30.12,4X,F30.12) 
c         CLOSE(UNIT=114)
1000  CONTINUE 
C 

      CALL CPU_TIME(ENDTIME)

      WRITE(*,*) "SHEAR RATE TO BE WRITTEN : ",SR

      BACKSPACE (UNIT=1)
      WRITE (1,3) Z,RMU,B,SR,E,H0,THI,NFLOW,ENDTIME-STARTTIME
3     FORMAT(F5.1,2X,F4.1,4X,F8.1,6X,F8.2,6X,F5.2,
     &2X,F5.2,2X,F5.2,4X,I1,F10.1)
      CLOSE (UNIT=1)
      close (unit=2) 
      CLOSE (UNIT=7)
      CLOSE (UNIT=89)
      CLOSE (UNIT=113)
1100  STOP

C 
      open (unit=15, file='q2.dat')
      open (unit=2, file='tstep.dat')
      open (unit=1, file='inp.dat')

      open (unit=40, file='q2x.dat',status='UNKNOWN')
      READ (1,*) Z,RMU,B,SR,E,H0,THI 
      READ (2,*) NTIWID, NTRAJ

1     FORMAT(F11.8,4X,F8.4,4X,F16.12,4X,F16.12)
2     FORMAT(2X,F8.4,3X,F5.2,3X,F5.2,4X,F8.4,4X,F16.5,4X,F16.5,
     &4X,I3,3X,I3,A4,I3,3X,F10.5,A4,F10.5)



C     EXTRAPOLATION TO ZERO STEP SIZE 
      DELTAT=0.D0
      IFLAG=0
C 
 
c     Q2 VALUES

      DO 31 J=1,NTIWID
        DO 31 I=1,NOUT
c         write(*,*) I,J
         READ(15,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)
c         WRITE(*,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)

31    CONTINUE

      DO 33 I=1,NOUT
        DO 34 J=1,NTIWID
          XARR(J)=TEMP1(J,I)
          YARR(J)=TEMP2(J,I)
          SIGARR(J)=TEMP3(J,I)
c          write(*,*) XARR(J),YARR(J),SIGARR(J)
34    CONTINUE
C     IF NOPT=-1, TEXTRA HAS NOT BEEN ABLE TO EXTRAPOLATE.
C     IN SUCH A CASE, THE ERRLEV PARAMETER IS REDUCED 10 TIMES
C     AND THE EXTRAPOLATION IS REPEATED. IF, EVEN AFTER 6 TIMES
C     THE RESULT IS THE SAME, ALL VALUES ARE REPORTED AS 0.0
C
      NOPT = -1
      ERRLEV = 0.25D0
      DO 32 WHILE ((NOPT.EQ.-1).AND.(ERRLEV.GT.0.0000025D0))
        CALL TEXTRA(XARR,YARR,SIGARR,NTIWID,NDATM,IFLAG,ERRLEV)
        ERRLEV = ERRLEV/10.D0
32    CONTINUE
      IF(NOPT.EQ.-1) THEN
         YOPT=0.D0
         SIGOPT=0.D0
         NOPT=0
         NDUOPT=0
         XOPT=0
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      IF(IFLAG.EQ.0) THEN
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      WRITE (40,2) SR, E , H0, OUTIME(I), YOPT, SIGOPT, NOPT, NDUOPT,
     &' of ', NTIWID, ALIOPT,' +- ',VLIOPT
33    CONTINUE

c1100  STOP 



      close (unit=15)
      close (unit=16)
      close (unit=17)
      close (unit=18)
      close (unit=19)
      close (unit=20)
      close (unit=21)

      close (unit=40)
      close (unit=41)
      close (unit=42)
      close (unit=43)
      close (unit=44)
      close (unit=45)
      close (unit=46)

      close (unit=1)
      close (unit=2)




c1100  STOP 


      END 
C 





      SUBROUTINE HISET(Y)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /STEPBL/ THI,B,ZMU,RMU2,BAUXQ,BAUXR,DTH,DTQ,SRDT,
     &SRDTH,C1P,C2P,E,AMPL2,BETA,AB,A2,A4,AUX1,AUX2,
     &AUX3,AUX4,AUX5,S,QALPH,GEE1,GEE2,GEE3,GEE4,FKAPPA(3,3),COEFF(4)  
       
      DATA C23/0.6666666666666667D0/
      DATA C43/1.3333333333333333D0/
      DATA C143/4.6666666666666670D0/
      DATA C83/2.6666666666666665D0/   
    
      REAL*8 COMP,COMPD,COMPI,Y1,ALPHA
      REAL*8 AUXD1,AUXD2,AUXD3,RES
      Y1=SQRT(Y)
      ALPHA=0.75D0*AB
      NTHI=NINT(THI)

c     Depending on THI, variables for the 
c     hydrodynamic tensor will be defined
C
C     I am considering HI tensor
C     to be of the form : 
C     \zeta Omega = (\alpha)/(Q)[A\delta+BQQ/Q^2]
C     so A = AUX2, B=AUX3 for both cases
C     In Regularized Oseen Burger case,
C     AUX4=M,AUX5=N
C
      SELECT CASE(NTHI)
C     FOR REGULARIZED OSEEN BURGERS          
          CASE(1)
              AUX1=Y+C43*A2
              AUX4=Y**(3.D0)+C143*A2*Y*Y+8.D0*A4*Y
              AUX5=Y**(3.D0)+2.D0*A2*Y*Y-C83*A4*Y
              AUX2=AUX4/(AUX1**(3.D0))
              AUX3=AUX5/(AUX1**(3.D00))
              AUXD1=((C43*A2)+(7.D0*Y))/(AUX1)
              AUXD2=(12.D0*(Y**(3.D0))+10.D0*C83*A2*Y*Y+4.D0*C83*A4*Y)
     &/(AUX1**(3.D0))
              AUXD3=(AUX2+AUX3)*AUXD1
              S=0.5D0*(AUXD3-AUXD2)
c     FOR RPY CASE      
          CASE(2)
              COMP=Y1/(2.D0*AB)
              COMPD=2.D0*COMP
              COMPI=1.D0/(COMPD) 
              IF(COMP.GE.1)THEN
                  AUX2=1.D0+(C23*COMPI*COMPI)
                  AUX3=1.D0-(2.D0*COMPI*COMPI) 
              ELSE
                  AUX2=(C43*COMPD)-(0.375D0*COMPD*COMPD)
                  AUX3=(0.125D0*COMPD*COMPD)
              ENDIF  
              S=AUX3
          CASE DEFAULT
              AUX1=Y
              AUX2=1.D0 
              AUX3=1.D0 
              S=1.D0
      END SELECT
      RES=S-AUX3
C     RES=S-B=A-r
      BETA=1.D0-(((ALPHA)/(Y1))*(AUX2+AUX3)) 
      AMPL2=(BETA)/((E*BETA)+1.D0)
      QALPH=(Y1)-(ALPHA*AUX2)
      GEE1=AMPL2*(((ALPHA*AUX3)/(BETA*QALPH))+E)
      GEE2=2.D0*(((AMPL2*AMPL2*ALPHA*AUX3)/(BETA*BETA*Y))
     &-(QALPH/Y)*GEE1 
     &-(E*(AMPL2*AMPL2*ALPHA*(RES)/(BETA*BETA*Y)))
     &-(E*ALPHA*(RES)*AMPL2/Y))
      GEE3=AMPL2/BETA
      GEE4=((2.D0*AMPL2*AMPL2*ALPHA*S)/(BETA*BETA*Y1))+(3.D0*AMPL2)
 
      RETURN
      END

      SUBROUTINE SEMIMP(Q)
C     Time step: semi-implicit predictor-corrector scheme for FENE+IV model
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (FOURPI=12.5663706143591725D0)
      COMMON /STEPBL/ THI,B,ZMU,RMU2,BAUXQ,BAUXR,DTH,DTQ,SRDT,
     &SRDTH,C1P,C2P,E,AMPL2,BETA,AB,A2,A4,AUX1,AUX2,
     &AUX3,AUX4,AUX5,S,QALPH,GEE1,GEE2,GEE3,GEE4,FKAPPA(3,3),COEFF(4) 

      DATA C23/0.6666666666666667D0/
      DATA C43/1.3333333333333333D0/
      DATA C143/4.6666666666666670D0/
      DATA C83/2.6666666666666665D0/

      REAL*8 Q(3),QAUX(3),KDQ(3),KDQAUX(3)
      REAL*8 B_D(3,3)
      REAL*8 S1,S2,S3,IVH,KH,IVAUX,KAUX
c     DEFINING VARIABLES FOR MARKO-SIGGIA
      REAL*8 SQB,QF,YFAC,AI,DEN
      REAL*8 G
45    FORMAT(F20.16)
      SQB=SQRT(B)
      KH=0.D0
      KAUX=0.D0
      NTHI=NINT(THI)
      DO 179 I=1,3
          KDQ(I)=0.D0
          KDQAUX(I)=0.D0
179   CONTINUE




c      WRITE(*,*) "REACHED SI"
C     B_D refers to the square root of the diffusion tensor
c     To be constructed at the beginning of the timestep
c     based on initial values of Q
C     S refers to dot product of diffusion tensor with 
c     the weiner numbers, W

C     Auxiliary parameters 
      QL=Q(1)*Q(1)+Q(2)*Q(2)+Q(3)*Q(3)
      SQQL=SQRT(QL)
      QF=SQQL/SQB
      CALL HISET(QL)
      YFAC=((0.5D0)/((1.D0-QF)**2.D0))-0.5D0+(2.D0*QF)

      HF=((1.D0)/(3.D0))*(1.D0/QF)*YFAC
      GT=-(GEE2*DTH/SQQL)
C     GENERAL FORMULA FOR KAPPA:QQ/Q^2

      DO 22 I=1,3
          DO 23 J=1,3
              KH=KH+(FKAPPA(I,J)*Q(J)*Q(I))
23        CONTINUE
22    CONTINUE

      KH=(E*AMPL2*KH*2.D0*DTH/QL)
24    FORMAT(F20.16)

      HE=-ZMU*EXP(-QL/RMU2)
      T0=DTH*AMPL2*(HF+HE) + GT + KH
      T1=1.D0-T0
C     Construction of suitable random numbers 
      W1=RANULS()-0.5D0
      W2=RANULS()-0.5D0
      W3=RANULS()-0.5D0
      W1=W1*(C1P*W1*W1+C2P)
      W2=W2*(C1P*W2*W2+C2P)
      W3=W3*(C1P*W3*W3+C2P)
C 

C     Construction of B_D : square root of Diffusion Tensor
C     See Eqn (31) in Prabhakar's 2002 paper in J.Rheology

      GEE = (QALPH/SQQL)
      GEETIL = -(QALPH/SQQL)*GEE1
c      GEETIL=-(E*AMPL2*BETA)-((AUX3*0.75*AB)/SQQL)
c      BOTH DEFINITIONS OF GEETIL ARE ESSENTIALLY THE SAME
c      JUST MAKING IT MATCH THE EXPRESSION IN MY DERIVATION

      AUXF = SQRT(GEE)
      AUXS = (SQRT(GEE+GEETIL)-SQRT(GEE))/QL

      B_D(1,1) = AUXF + (AUXS*Q(1)*Q(1))
      B_D(1,2) = AUXS*Q(1)*Q(2)
      B_D(1,3) = AUXS*Q(1)*Q(3)

      B_D(2,1) = AUXS*Q(2)*Q(1)
      B_D(2,2) = AUXF + (AUXS*Q(2)*Q(2))
      B_D(2,3) = AUXS*Q(2)*Q(3)

      B_D(3,1) = AUXS*Q(3)*Q(1)
      B_D(3,2) = AUXS*Q(3)*Q(2)
      B_D(3,3) = AUXF + (AUXS*Q(3)*Q(3))

c     S = (B_D).W

      S1 = B_D(1,1)*W1 + B_D(1,2)*W2 + B_D(1,3)*W3
      S2 = B_D(2,1)*W1 + B_D(2,2)*W2 + B_D(2,3)*W3
      S3 = B_D(3,1)*W1 + B_D(3,2)*W2 + B_D(3,3)*W3

c     KDQ = (FKAPPA).Q

      DO 11 I=1,3
          DO 12 J=1,3
              KDQ(I)=KDQ(I)+(FKAPPA(I,J)*Q(J)*2.D0*DTH)
12        CONTINUE
11    CONTINUE

C     Predictor step 

      QAUX(1)=T1*Q(1)+KDQ(1)+S1
      QAUX(2)=T1*Q(2)+KDQ(2)+S2
      QAUX(3)=T1*Q(3)+KDQ(3)+S3
      QAUXL=QAUX(1)*QAUX(1)+QAUX(2)*QAUX(2)+QAUX(3)*QAUX(3)


c     All the auxilary functions have to be
c     calculated again. Bit of a pain.
c     Better to do this through function calls.

      SQAUXQL=SQRT(QAUXL)
      CALL HISET(QAUXL)

      AI=DTQ*AMPL2
      DEN=((2.D0*AI)+(3.D0))
      GTAUX=-(GEE2*DTQ/SQAUXQL)

C     GENERAL FORMULA FOR KAPPA:QQ/Q^2

      DO 32 I=1,3
          DO 33 J=1,3
              KAUX=KAUX+(FKAPPA(I,J)*QAUX(J)*QAUX(I))
33        CONTINUE
32    CONTINUE
      KAUX=(E*AMPL2*KAUX*DTH)/QAUXL
      HEAUX=-DTQ*AMPL2*ZMU*EXP(-QAUXL/RMU2)
      T4=GTAUX+KAUX+HEAUX
      T3=1.D0-T4

c     KDQAUX = (FKAPPA).QAUX

      DO 111 I=1,3
          DO 112 J=1,3
              KDQAUX(I)=KDQAUX(I)+(FKAPPA(I,J)*QAUX(J)*2.D0*DTH)
112        CONTINUE
111    CONTINUE

C     Corrector step 
      Q(1)=T3*QAUX(1)+0.5D0*(KDQAUX(1)-KDQ(1))+0.5D0*T0*Q(1)
      Q(2)=T3*QAUX(2)+0.5D0*(KDQAUX(2)-KDQ(2))+0.5D0*T0*Q(2)
      Q(3)=T3*QAUX(3)+0.5D0*(KDQAUX(3)-KDQ(3))+0.5D0*T0*Q(3)
C     Exact solution of the implicit equation for the 
C     length of the connector vector 
      QL=Q(1)*Q(1)+Q(2)*Q(2)+Q(3)*Q(3)
      SQQL=SQRT(QL)
      G=SQQL/SQB

      COEFF(4)=1.D0
      COEFF(3)=-(1.5D0*(4.D0+(3.D0*AI)+(2.D0*G)))/DEN
      COEFF(2)=(3.D0*(1.D0+AI+(2.D0*G)))/DEN
      COEFF(1)=-(3.D0*G)/DEN
c     INITIAL GUESS
      R=(1.D0)-SQRT(AI/(6.D0*G))
      IF(G < 1.D0) THEN
          R=G/(1.D0+AI)
      ENDIF
 
      IF(G>100.D0) THEN
c     DON'T POLISH ROOT IF G>100
          R=R
      ELSE
          CALL POLISH_POLY_ROOT(COEFF,R,1.D-06)
      ENDIF      

c     YADA
c     Rescaling to obtain the proper length 
      RED=R/G
      Q(1)=RED*Q(1)
      Q(2)=RED*Q(2)
      Q(3)=RED*Q(3)
      RETURN
      END


      SUBROUTINE RANILS(ISEED) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     Choice of ISEED: 0 <= ISEED <= 2000000000 (2E+9); 
C     ISEED can, for example, be formed from current time and date: 
C     2 digits each for seconds, minutes, hours, day, month 
C                    or minutes, hours, day, month, year 
      PARAMETER (IN=2147483563,IK=40014,IQ=53668,IR=12211,NTAB=32) 
      INTEGER IV(NTAB) 
      COMMON /RANBLS/ IDUM,IDUM2,IY,IV 
C     Initial seeds for two random number generators 
      IDUM=ISEED+123456789 
      IDUM2=IDUM 
C     Load the shuffle table (after 8 warm-ups) 
      DO 10 J=NTAB+8,1,-1 
         K=IDUM/IQ 
         IDUM=IK*(IDUM-K*IQ)-K*IR 
         IF(IDUM.LT.0) IDUM=IDUM+IN 
         IF(J.LE.NTAB) IV(J)=IDUM 
10    CONTINUE 
      IY=IV(1) 
      RETURN 
      END 
C 
      FUNCTION RANULS() 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER (IN1=2147483563,IK1=40014,IQ1=53668,IR1=12211, 
     &           IN2=2147483399,IK2=40692,IQ2=52774,IR2=3791, 
     &           NTAB=32,AN=1.D0/IN1,INM1=IN1-1,NDIV=1+INM1/NTAB) 
      INTEGER IV(NTAB) 
      COMMON /RANBLS/ IDUM,IDUM2,IY,IV 
C     Linear congruential generator 1 
      K=IDUM/IQ1 
      IDUM=IK1*(IDUM-K*IQ1)-K*IR1 
      IF(IDUM.LT.0) IDUM=IDUM+IN1 
C     Linear congruential generator 2 
      K=IDUM2/IQ2 
      IDUM2=IK2*(IDUM2-K*IQ2)-K*IR2 
      IF(IDUM2.LT.0) IDUM2=IDUM2+IN2 
C     Shuffling and subtracting 
      J=1+IY/NDIV 
      IY=IV(J)-IDUM2 
      IV(J)=IDUM 
      IF(IY.LT.1) IY=IY+INM1 
      RANULS=AN*IY 
      RETURN 
      END 

C 
      FUNCTION RANGLS()
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE IFLAG,GAUSS2
      DATA IFLAG/0/
      IF(IFLAG.EQ.0) THEN
10       CONTINUE
C        Pair of uniform random numbers in [-1,1]x[-1,1] 
         X1=2.D0*RANULS()-1.D0
         X2=2.D0*RANULS()-1.D0
C        If not in the unit circle, try again 
         XSQ=X1*X1+X2*X2
         IF(XSQ.GE.1.D0) GOTO 10
C        Pair of Gaussian random numbers; return one and 
C        save the other for next time 
         AUX=SQRT(-2.D0*LOG(XSQ)/XSQ)
         RANGLS=X1*AUX
         GAUSS2=X2*AUX
         IFLAG=1
      ELSE
         RANGLS=GAUSS2
         IFLAG=0
      ENDIF
      RETURN
      END
C 


C 
      SUBROUTINE TEXTRA(XARR,YARR,SIGARR,NDAT,NDATM,IFLAG,ERRLEV) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NDATP=50) 
      REAL*8 XARR(NDATM),YARR(NDATM),SIGARR(NDATM)   
      REAL*8 A(NDATP),COVAR(NDATP,NDATP) 
      INTEGER LFLAG(NDATP) 
      COMMON /EXTRP/ NOPT,NDUOPT,XOPT,YOPT,SIGOPT,ALIOPT,VLIOPT 
      EXTERNAL fpoly 
C 
C     Sorting the data such that XARR(1).LE.XARR(2).LE.XARR(3) ... 
C     (by straight insertion) 
      DO 10 J=2,NDAT 
         X=XARR(J) 
         Y=YARR(J) 
C         IF((Y+1.D0).EQ.Y) THEN
C             WRITE(*,*) "ERROR"
C         ENDIF 
         SIG=SIGARR(J) 
         DO 20 I=J-1,1,-1 
            IF(XARR(I).LE.X) GOTO 30 
            XARR(I+1)=XARR(I) 
            YARR(I+1)=YARR(I) 
            SIGARR(I+1)=SIGARR(I) 
20       CONTINUE 
         I=0 
30       XARR(I+1)=X 
         YARR(I+1)=Y 
         SIGARR(I+1)=SIG 
10    CONTINUE 
C 
      NOPT=-1 
      SIGOPT=6.022D23 
C     Fitting polynomials of various degrees N.LE.NMAX 
      NMAX=NDAT-2 
      IF(IFLAG.EQ.0) NMAX=NMAX+1 
      DO 1000 N=0,NMAX 
         NDATMI=N+2 
         IF(IFLAG.EQ.0.AND.N.GE.1) NDATMI=NDATMI-1 
C        Discarding data with large XARR 
         DO 500 NDATU=NDAT,NDATMI,-1 
            NDF=NDATU-NDATMI+1 
C           Least squares fit 
            DO 40 I=1,N+1 
               A(I)=0.D0 
               LFLAG(I)=1 
40          CONTINUE 
            DO 50 I=N+2,NDATP 
               A(I)=0.D0 
               LFLAG(I)=0 
50          CONTINUE 
            IF(N.GT.0) LFLAG(2)=IFLAG 
            CALL lfit(XARR,YARR,SIGARR,NDATU,A,LFLAG,NDATP, 
     &                                     COVAR,NDATP,TEST,fpoly) 
C           Chi-squared test; smaller statistical error bars? 
            IF(gammq(0.5D0*NDF,0.5D0*TEST).GT.ERRLEV) THEN 
               IF(SQRT(COVAR(1,1)).LT.SIGOPT) THEN 
                  YOPT=A(1) 
                  SIGOPT=SQRT(COVAR(1,1)) 
                  NOPT=N 
                  NDUOPT=NDATU 
                  XOPT=XARR(NDATU) 
                  ALIOPT=A(2) 
                  VLIOPT=SQRT(COVAR(2,2)) 
               ENDIF 
            ENDIF 
500      CONTINUE 
1000  CONTINUE 
C 
      RETURN 
      END   
C  
      SUBROUTINE lfit(x,y,sig,ndat,a,ia,ma,covar,npc,chisq,funcs) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER ma,ia(ma),npc,ndat,MMAX 
      REAL*8 chisq,a(ma),covar(npc,npc),sig(ndat),x(ndat),y(ndat) 
      EXTERNAL funcs 
      PARAMETER (MMAX=50) 
CU    USES covsrt,gaussj 
      INTEGER i,j,k,l,m,mfit 
      REAL*8 sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX) 
      mfit=0 
      do 11 j=1,ma 
        if(ia(j).ne.0) mfit=mfit+1 
11    continue 
      if(mfit.eq.0) pause 'lfit: no parameters to be fitted' 
      do 13 j=1,mfit 
        do 12 k=1,mfit 
          covar(j,k)=0.D0 
12      continue 
        beta(j)=0.D0 
13    continue 
      do 17 i=1,ndat 
        call funcs(x(i),afunc,ma) 
        ym=y(i) 
        if(mfit.lt.ma) then 
          do 14 j=1,ma 
            if(ia(j).eq.0) ym=ym-a(j)*afunc(j) 
14        continue 
        endif 
        sig2i=1.D0/sig(i)**(2.D0) 
        j=0 
        do 16 l=1,ma 
          if (ia(l).ne.0) then 
            j=j+1 
            wt=afunc(l)*sig2i 
            k=0 
            do 15 m=1,l 
              if (ia(m).ne.0) then 
                k=k+1 
                covar(j,k)=covar(j,k)+wt*afunc(m) 
              endif 
15          continue 
            beta(j)=beta(j)+ym*wt 
          endif 
16      continue 
17    continue 
      do 19 j=2,mfit 
        do 18 k=1,j-1 
          covar(k,j)=covar(j,k) 
18      continue 
19    continue 
      call gaussj(covar,mfit,npc,beta,1,1) 
      j=0 
      do 21 l=1,ma 
        if(ia(l).ne.0) then 
          j=j+1 
          a(l)=beta(j) 
        endif 
21    continue 
      chisq=0.D0 
      do 23 i=1,ndat 
        call funcs(x(i),afunc,ma) 
        sum=0. 
        do 22 j=1,ma 
          sum=sum+a(j)*afunc(j) 
22      continue 
        chisq=chisq+((y(i)-sum)/sig(i))**2 
23    continue 
      call covsrt(covar,npc,ma,ia,mfit) 
      return 
      END 
C   
      SUBROUTINE fpoly(x,p,np) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER np 
      REAL*8 x,p(np) 
      INTEGER j 
      p(1)=1.D0 
      do 11 j=2,np 
        p(j)=p(j-1)*x 
11    continue 
      return 
      END 
C   
      SUBROUTINE covsrt(covar,npc,ma,ia,mfit) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER ma,mfit,npc,ia(ma) 
      REAL*8 covar(npc,npc) 
      INTEGER i,j,k 
      REAL*8 swap 
      do 12 i=mfit+1,ma 
        do 11 j=1,i 
          covar(i,j)=0.D0 
          covar(j,i)=0.D0 
11      continue 
12    continue 
      k=mfit 
      do 15 j=ma,1,-1 
        if(ia(j).ne.0)then 
          do 13 i=1,ma 
            swap=covar(i,k) 
            covar(i,k)=covar(i,j) 
            covar(i,j)=swap 
13        continue 
          do 14 i=1,ma 
            swap=covar(k,i) 
            covar(k,i)=covar(j,i) 
            covar(j,i)=swap 
14        continue 
          k=k-1 
        endif 
15    continue 
      return 
      END 
C   
      SUBROUTINE gaussj(a,n,np,b,m,mp) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER m,mp,n,np,NMAX 
      REAL*8 a(np,np),b(np,mp) 
      PARAMETER (NMAX=50) 
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX) 
      REAL*8 big,dum,pivinv 
      do 11 j=1,n 
        ipiv(j)=0 
11    continue 
      do 22 i=1,n 
        big=0.D0 
        do 13 j=1,n 
          if(ipiv(j).ne.1)then 
            do 12 k=1,n 
              if (ipiv(k).eq.0) then 
                if (abs(a(j,k)).ge.big)then 
                  big=abs(a(j,k)) 
                  irow=j 
                  icol=k 
                endif 
              else if (ipiv(k).gt.1) then 
                pause 'singular matrix in gaussj' 
              endif 
12          continue 
          endif 
13      continue 
        ipiv(icol)=ipiv(icol)+1 
        if (irow.ne.icol) then 
          do 14 l=1,n 
            dum=a(irow,l) 
            a(irow,l)=a(icol,l) 
            a(icol,l)=dum 
14        continue 
          do 15 l=1,m 
            dum=b(irow,l) 
            b(irow,l)=b(icol,l) 
            b(icol,l)=dum 
15        continue 
        endif 
        indxr(i)=irow 
        indxc(i)=icol 
        if (a(icol,icol).eq.0.D0) pause 'singular matrix in gaussj' 
        pivinv=1.D0/a(icol,icol) 
        a(icol,icol)=1.D0 
        do 16 l=1,n 
          a(icol,l)=a(icol,l)*pivinv 
16      continue 
        do 17 l=1,m 
          b(icol,l)=b(icol,l)*pivinv 
17      continue 
        do 21 ll=1,n 
          if(ll.ne.icol)then 
            dum=a(ll,icol) 
            a(ll,icol)=0.D0 
            do 18 l=1,n 
              a(ll,l)=a(ll,l)-a(icol,l)*dum 
18          continue 
            do 19 l=1,m 
              b(ll,l)=b(ll,l)-b(icol,l)*dum 
19          continue 
          endif 
21      continue 
22    continue 
      do 24 l=n,1,-1 
        if(indxr(l).ne.indxc(l))then 
          do 23 k=1,n 
            dum=a(k,indxr(l)) 
            a(k,indxr(l))=a(k,indxc(l)) 
            a(k,indxc(l))=dum 
23        continue 
        endif 
24    continue 
      return 
      END 
C   
      FUNCTION gammq(a,x) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 a,gammq,x 
CU    USES gcf,gser 
      REAL*8 gammcf,gamser,gln 
      if(x.lt.0.D0.or.a.le.0.D0)pause 'bad arguments in gammq' 
      if(x.lt.a+1.D0)then 
        call gser(gamser,a,x,gln) 
        gammq=1.D0-gamser 
      else 
        call gcf(gammcf,a,x,gln) 
        gammq=gammcf 
      endif 
      return 
      END 
C  
      SUBROUTINE gcf(gammcf,a,x,gln) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER ITMAX 
      REAL*8 a,gammcf,gln,x,EPS,FPMIN 
      PARAMETER (ITMAX=100,EPS=3.D-07,FPMIN=1.D-30) 
CU    USES gammln 
      INTEGER i 
      REAL*8 an,b,c,d,del,h,gammln 
      gln=gammln(a) 
      b=x+1.D0-a 
      c=1.D0/FPMIN 
      d=1.D0/b 
      h=d 
      do 11 i=1,ITMAX 
        an=-i*(i-a) 
        b=b+2.D0 
        d=an*d+b 
        if(abs(d).lt.FPMIN)d=FPMIN 
        c=b+an/c 
        if(abs(c).lt.FPMIN)c=FPMIN 
        d=1.D0/d 
        del=d*c 
        h=h*del 
        if(abs(del-1.D0).lt.EPS)goto 1 
11    continue 
      pause 'a too large, ITMAX too small in gcf' 
1     gammcf=exp(-x+a*log(x)-gln)*h 
      return 
      END 
C   
      SUBROUTINE gser(gamser,a,x,gln) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER ITMAX 
      REAL*8 a,gamser,gln,x,EPS 
      PARAMETER (ITMAX=100,EPS=3.D-07) 
CU    USES gammln 
      INTEGER n 
      REAL*8 ap,del,sum,gammln 
      gln=gammln(a) 
      if(x.le.0.D0)then 
        if(x.lt.0.D0)pause 'x < 0 in gser' 
        gamser=0.D0 
        return 
      endif 
      ap=a 
      sum=1.D0/a 
      del=sum 
      do 11 n=1,ITMAX 
        ap=ap+1.D0 
        del=del*x/ap 
        sum=sum+del 
        if(abs(del).lt.abs(sum)*EPS)goto 1 
11    continue 
      pause 'a too large, ITMAX too small in gser' 
1     gamser=sum*exp(-x+a*log(x)-gln) 
      return 
      END 
C   
      FUNCTION gammln(xx) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 gammln,xx 
      INTEGER j 
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6) 
      SAVE cof,stp 
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, 
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, 
     *-.5395239384953d-5,2.5066282746310005d0/ 
      x=xx 
      y=x 
      tmp=x+5.5d0 
      tmp=(x+0.5d0)*log(tmp)-tmp 
      ser=1.000000000190015d0 
      do 11 j=1,6 
        y=y+1.d0 
        ser=ser+cof(j)/y 
11    continue 
      gammln=tmp+log(stp*ser/x) 
      return 
      END 
C  

      FUNCTION betafn(x,y)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 betafnans,x,y
C     USES gammln
c     Returns  the  value  of  the  beta  function B(x,y).
      REAL*8 gammln
      betafnans=exp(gammln(x)+gammln(y)-gammln(x+y))
      return
      END

      SUBROUTINE POLISH_POLY_ROOT(C,XIN,ATOL)
C     PAGE 470 OF NUMERICAL RECIPES,3RD EDITION
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 C(4)
      PARAMETER(IMAX=30)
      PARAMETER(N=4)
      X=XIN
      DO 21 ITER=1,IMAX
          X0 = X
          P = C(N)*X + C(N-1)
          P1 = C(N)
          DO 20 I=(N-2),1,-1
              P1 = P + P1*X
              P  = C(I) + P*X
20        CONTINUE
          IF (ABS(P1) .EQ. 0.D0) THEN
              P1 = 6.D0 * C(4)
              P1 = 6.D0*(P/P1) 
c             SIGN(A,B) = SGN(B)*MOD(A) 
              X  = X - SIGN(1.D0,P1)*(ABS(P1)**(1.D0/3.D0))
          ELSE
              X = X - (P/P1)
          ENDIF
          IF (ABS(P).LE.ATOL) EXIT
21    CONTINUE
      IF (ITER>IMAX) THEN
          WRITE(*,*) "POLY-ROOT : LOOP EXCEEDED"
      ENDIF
      XIN = X
C      WRITE(*,*) "YAY",X
      RETURN
      END
