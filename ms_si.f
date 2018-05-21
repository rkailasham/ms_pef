c     MARKO-SIGGIA DUMBBELLS WITH 
C     INTERNAL VISCOSITY
c     + EXCLUDED VOLUME INTERACTIONS 
C     + HYDRODYNAMIC INTERACTIONS
C     SECOND ORDER SEMI-IMPLICIT PREDICTOR-CORRECTOR ALGORITHM
c
c     NOTE : WHEN RUNNING EQUILIBRIUM SIMULATIONS, DIVERT PROGRAM
c     FLOW FROM ENTERING TEXTRA. WHEN SR=0.0, MATERIAL FUNCTIONS BECOME
c     UNDEFINED.
C     THIS CREATES HAVOC WITH TEXTRA
C
c     2-APR-2018  : FINDING THE RIGHT SOLUTION FROM THE
C     THREE ROOTS GIVEN BY THE CUBIC EQUATION SOLVER
C
c     3-APR-2018  : DISCRIMINANT OF CUBIC EQUATIONS FOR  
c     MARKO-SIGGIA CASE IS POSITIVE, CANNOT USE TRIG. SOLVER
c     HAVE IMPLEMENTED NEWTON-RAPHSON+POLISHING    
c
c     SEE SUPPLEMENTARY INFORMATION OF SCHROEDER'S PAPER
C     ["Nonequilibrium thermodynamics of dilute polymer solutions in flow"
C     JCP 141, 174903 (2014)]
C     TO UNDERSTAND REGIME 1,2,3
C
c     USED TO CALCULATE DISSIPATED WORK, ALONG WITH OTHER INFORMATION SUCH
c     AS THE AVERAGE WORK AND THE FREE-ENERGY DIFFERENCE
C
c     CHANGE TMAX TO CHANGE THE RAMP RATE. INITIAL AND FINAL SHEAR 
C     RATES ARE SPECIFIED IN THE INPUT FILE
C



      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NDATM=50,NOUT1=50,NOUT=50,NOUT3=2000) 
      PARAMETER (NBINS=50)
      PARAMETER (NDB=10000000)
      PARAMETER (PI = 3.1415926535897931D0)
      REAL*8 FENFAC,TEMPB
      REAL*8 META1,MPSI2
      REAL*8 AQ2,MQ2,VQ2
      REAL*8 QMAG,QLEN,TEMP12,TEMP23
      REAL*8 METE,METD,MET
      REAL*8 META
      REAL*8 AVQ1(NOUT1),MEQ1,ERQ1(NOUT1)
      REAL*8 AVQ2(NOUT),MEQ2,ERQ2(NOUT)
      REAL*8 AVQ3(NOUT3),MEQ3,ERQ3(NOUT3)
      REAL*8 AVTE(NOUT),ERTE(NOUT)
      REAL*8 AVTD(NOUT),ERTD(NOUT)
      REAL*8 AVT(NOUT),ERT(NOUT)
      REAL*8 AVETA(NOUT),ERETA(NOUT)
      REAL*8 AVETA1(NOUT),ERETA1(NOUT)
      REAL*8 AVPSI2(NOUT),ERPSI2(NOUT)
      REAL*8 DTARR(10) 
      REAL*8 XARR(NDATM),YARR(NDATM),SIGARR(NDATM) 
      REAL*8 TEMP1(NDATM, NOUT),TEMP2(NDATM, NOUT),TEMP3(NDATM, NOUT)
      REAL*8 OUTIME(NOUT)
      REAL*8 ESR(NOUT1),TSR(NOUT),RSR(NOUT3)
      REAL*8 Q(3)
      REAL*8 DB(NDB,3)
      REAL*8 EQBFAC
      REAL*8 RAMP,CF1,CF2
      REAL*8 QAV1,QER1,QAV3,QER3,FWAV1,FWER1,FWAV3,FWER3
      REAL*8 AVW,ERW,MEW,AW,VW
C      REAL*8 KAPPA(3,3)
      REAL*8 KAPPA_TEMP(3,3),FLOW_TMP(3,3)
      REAL*8 WORK2(NDB)
      REAL*8 DELF,ER_DELF
      INTEGER NTIARR(10) 
      INTEGER TIME (8)
      CHARACTER*1 TAB  
      CHARACTER (LEN=12) CLK(3)
      CHARACTER (LEN=512) FPATH
      CHARACTER (LEN=512) CPATH
      PARAMETER(FPATH="/home/rkai0001/wm73/rkai0001/marko_siggia/ms_data
     &base/b477/eqbconfigs.dat")
      PARAMETER(CPATH="finconfigs.dat")
      COMMON /STEPBL/ THI,B,ZMU,RMU2,BAUXQ,BAUXR,DTH,DTQ,SRDT,
     &SRDTH,C1P,C2P,E,AMPL2,BETA,AB,A2,A4,AUX1,AUX2,
     &AUX3,AUX4,AUX5,S,QALPH,GEE1,GEE2,GEE3,GEE4,FKAPPA(3,3),COEFF(4) 
      COMMON /EXTRP/ NOPT,NDUOPT,XOPT,YOPT,SIGOPT,ALIOPT,VLIOPT 
      LOGICAL THERE
      LOGICAL CTHERE
c      INTEGER :: THI
C     Excluded-volume and FENE parameters 
c     Z = Solvent Quality, RMU = EV parameter, B= FENE parameter
c     SR1 = INITIAL PECLET NUMBER, SR2=FINAL PECLET NUMBER 
C     E = IV PARAMETER,
C     H0 = HYDRODYNAMIC INTERACTION PARAMETER
c     B_RAD = BEAD_RADIUS, IN DIMENSIONLESS UNITS
C     THI = EXPRESSION TO 
C     BE USED OFR THE HYDRODYNAMIC INTERACTION TENSOR.
C     THI = 1 FOR REGULARIZED OSEEN-BURGER
C     THI = 2 FOR ROTNE-PRAGER-YAMAKAWA 
      open (unit=1, file='inp.dat') 
      READ (1,*) Z,RMU,B,SR1,SR2,CIV,ETA_S,B_RAD,THI,INPAR,NFLOW 
      open (unit=2, file='tstep.dat') 
c      open (unit=3, file='eta1.dat',STATUS='UNKNOWN')
      open (unit=91, file='q2_r1.dat', STATUS='UNKNOWN')
      open (unit=7, file='q2_r2.dat', STATUS='UNKNOWN')
      open (unit=8, file='q2_r3.dat', STATUS='UNKNOWN')
      open (unit=10, file='output.dat', STATUS='UNKNOWN')
      open (unit=11, file='work_traj.dat', STATUS='UNKNOWN')

      INQUIRE (FILE=FPATH,EXIST=THERE)
      INQUIRE (FILE=CPATH,EXIST=CTHERE)

C     REGIME 1
      TEQB=5.D0
C     REGIME 2
      TMAX=5.D0
C     REGIME 3
      TREL=200.D0

c

      READ (2,*) NTIWID, NTRAJ
      DO 15 I = 1, NTIWID
          READ (2,*) NTIARR(I)
          DTARR(I)=TMAX/NTIARR(I)
15    CONTINUE

      TAB=CHAR(9)
      E=CIV/(3.D0*PI*ETA_S*B_RAD)

C 
      ZMU=Z/RMU**(5.D0) 
      RMU2=(2.D0)*RMU*RMU
      NTHI=NINT(THI) 

c     WRITE FORMATS


45    FORMAT(F20.16)
8     FORMAT(I10,4X,F20.16,4X,F20.16,4X,F20.16)
9     FORMAT(F11.8,4X,F30.16,4X,F30.16,4X,F30.16)
4     FORMAT(F11.8,4X,F10.5,4X,F30.12,4X,F30.12) 
23    FORMAT(I12)
24    FORMAT(F24.16) 
29    FORMAT(I8,4X,F12.5)
81    FORMAT(1A,F12.5,2X,1A,F12.5)

C     DECIDES TYPE OF HI TENSOR TO BE USED : (1) REGULARIZED OSEEN-BURGERS (ROB) 
C                                            (2) ROTNE-PRAGER-YAMAKAWA (RPY)
      SELECT CASE (NTHI)
          CASE (1)
              WRITE(*,*) "THI : ",THI
              WRITE(*,*) "TYPE OF HI : ROBT"
              AB=B_RAD
          CASE (2)
              WRITE(*,*) "THI : ",THI
              WRITE(*,*) "TYPE OF HI : RPY"
              AB=B_RAD 
         CASE DEFAULT
              WRITE(*,*) "HI OPTION NOT VALID"
              WRITE(*,*) "SETTING H0=0.0"
              AB=0.D0
      END SELECT

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

C     STORING THE FLOW TENSOR. HELPS IN SWITCHING FROM
C     EQUILIBRATION TO PRODUCTION RUNS
      DO 112 I=1,3
          DO 113 J=1,3
              KAPPA_TEMP(I,J)=FLOW_TMP(I,J)*SR1
113        CONTINUE
112    CONTINUE 


      CALL CPU_TIME(STARTTIME)

C     DECIDES IF INPUT CONFIGS ARE TAKEN FROM (1) EQUILIBRATED DATABASE OR
C     (2) SAVED CONFIGS FROM PREVIOUS RUN 


      SELECT CASE (INPAR)
          CASE (1)
              WRITE(*,*) "INPAR : ",INPAR
              IF(THERE)THEN
                  OPEN(UNIT=114,file=FPATH)
                  WRITE(*,*) "LOADING FROM EQUILIBRATED DATABASE.."
                  DO 12 I=1,NDB
                      READ(114,8) K,DB(I,1),DB(I,2),DB(I,3)
12                CONTINUE
                  CLOSE(UNIT=114)
              ELSE
                  STOP "eqbconfigs.dat file not found. Execution 
     &terminated"
              ENDIF  
              ISEED=20171113              
              EQBFAC=1.D0
C             When loading from equilibrated database, regime 1 
c             is forced to be the equilibration regime
c              SR1=0.D0
              WRITE(*,*) "LOADED DATABASE"
          CASE (2)
              WRITE(*,*) "INPAR : ",INPAR
              IF(CTHERE)THEN
                  OPEN(UNIT=115,file=CPATH)
                  WRITE(*,*) "READING FROM FINAL CONFIGS OF PREVIOUS 
     &RUN.."
                  DO 13 I=1,(NTIWID*NTRAJ)
                      READ(115,9) TIMEI,DB(I,1),DB(I,2),DB(I,3)
13                CONTINUE
                  READ(115,*) ISEED
                  CLOSE(UNIT=115)
              ELSE
                  STOP "finconfigs.dat file not found. Execution 
     &terminated"
              ENDIF
              EQBFAC=0.D0
          CASE DEFAULT
              WRITE(*,*) "READ_FROM_INPUT OPTION NOT VALID"
              STOP "USE INPAR=1 FOR EQB DBASE, INPAR=2 FOR RESUMING 
     &FROM PREVIOUS RUN"
      END SELECT

      H0=AB/SQRT(PI)
      A2=AB*AB
      A4=A2*A2
      WRITE(*,45) E
      WRITE(*,45) AB


      RAMP=(SR2-SR1)/(TMAX)

      open (unit=112,file='finconfigs.dat',STATUS='UNKNOWN')

C     Loop for different time step widths 

      DO 1000 IDT=1,NTIWID
C        Auxiliary parameters 
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

         DO 86 I=1,NOUT1
             AVQ1(I)=0.D0
             ERQ1(I)=0.D0
86       CONTINUE

         DO 25 I=1,NOUT
             AVQ2(I)=0.D0
             ERQ2(I)=0.D0
             AVETA1(I)=0.D0
             ERETA1(I)=0.D0
25       CONTINUE

         DO 85 I=1,NOUT3
             AVQ3(I)=0.D0
             ERQ3(I)=0.D0
85       CONTINUE

         QAV1=0.D0
         QER1=0.D0
         QAV3=0.D0
         QER3=0.D0
         FWAV1=0.D0
         FWER1=0.D0
         FWAV3=0.D0
         FWER3=0.D0




C        A FRESH SEED IS GIVEN FOR EACH TIME-STEP WIDTH
C        TO ENSURE NON-OCCURENCE OF PERIOD EXHAUSTION

         CALL DATE_AND_TIME(CLK(1),CLK(2),CLK(3),TIME)
         ISEED=TIME(8)*100000+TIME(7)*1000+TIME(6)*10+TIME(5)
         ISEED=ISEED+111555         
         CALL RANILS(ISEED)
         CALL SRAND(ISEED)


         DO 100 ITRAJ=1,NTRAJ 
        
             WORK2(ITRAJ)=0.D0
             PICK=RAND()
             NSEED=NSEED+1
             CHANGE=NDB*PICK
             IF(INPAR.EQ.1)THEN
                 NCHOOSE=NINT(CHANGE)
             ELSE
                 NCHOOSE=((IDT-1)*NTRAJ)+(ITRAJ)
             ENDIF 
             Q(1)=DB(NCHOOSE,1)
             Q(2)=DB(NCHOOSE,2)
             Q(3)=DB(NCHOOSE,3)



             IF(MODULO(ITRAJ,100).EQ.0)THEN
             WRITE(*,*) "STATUS : EQB.TIME-STEP WIDTH : ",DELTAT,
     &"TRAJ # ",ITRAJ
             ENDIF

             NFIRST=EQBFAC/DELTAT
             DO 40 ITIME=1,NFIRST
                CALL SEMIMP(Q)
40           CONTINUE

             DO 78 I=1,3
                 DO 79 J=1,3
                     FKAPPA(I,J)=FLOW_TMP(I,J)*SR1
79               CONTINUE
78           CONTINUE 
 

             NSEC=TREL/DELTAT
             DO 55 ITIME=1,NSEC
                CALL SEMIMP(Q)
55           CONTINUE

c            REGIME 1 : RELAXATION OF INITIAL CONDITION
C            Relaxation of initial condition for 1 dimensionless time
c            irrespective of whether starting from equilibrated database
c            or finconfigsdat

             DO 95 I=1,3
                 DO 96 J=1,3
                     FKAPPA(I,J)=FLOW_TMP(I,J)*SR1
96               CONTINUE
95           CONTINUE 
             CF1=FKAPPA(1,1)
             CF2=FKAPPA(2,2)

             NEQB=TEQB/DELTAT
             IWAIT=0
             IOUT=0

             DO 50 ITIME=1,NEQB
                CALL SEMIMP(Q)
                IWAIT=IWAIT+1
                IF (IWAIT.EQ.NTIME) THEN
                    IWAIT=0
                    IOUT=IOUT+1
                    ESR(IOUT)=FKAPPA(1,1)
                    QMAG=Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3)
                    MEQ1=SQRT(QMAG)
                    AVQ1(IOUT)=AVQ1(IOUT)+(MEQ1)
                    ERQ1(IOUT)=ERQ1(IOUT)+(MEQ1*MEQ1)
                ENDIF
50           CONTINUE
             QMAG=Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3)
             MEQ1=SQRT(QMAG)

             QAV1=QAV1+(MEQ1)
             QER1=QER1+(MEQ1*MEQ1)
               
             TEMP12=SR1*(Q(1)*Q(1) - Q(2)*Q(2))
             FWAV1=FWAV1+(TEMP12)
             FWER1=FWER1+(TEMP12*TEMP12)
c            REGIME 2 : GRADUALLY RAMPING UP PECLET NUMBER
C            [Pe IS SAME AS DIMENSIONLESS SHEAR/ELONGATION RATE]


             IWAIT=0
             IOUT=0
C

C           Time integration: semi-implicit predictor-corrector scheme 
             DO 10 ITIME=1,NTIARR(IDT) 
                 DO 115 I=1,3
                     DO 116 J=1,3
                         FKAPPA(I,J)=FKAPPA(I,J)+(DELTAT*RAMP*
     &FLOW_TMP(I,J))
116                  CONTINUE
115              CONTINUE 

                 CALL SEMIMP(Q) 
                 IWAIT=IWAIT+1
                 QMAG =Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3)
                 TEMP12=(Q(1)*Q(1) - Q(2)*Q(2))
                 TEMP23=(Q(1)*Q(1)-0.5D0*Q(2)*Q(2)-0.5D0*Q(3)*Q(3))
                 IF (IWAIT.EQ.NTIME) THEN
                     IWAIT=0
                     IOUT=IOUT+1
                     TSR(IOUT)=FKAPPA(1,1)
                     MEQ2=SQRT(QMAG)
 
                     AVQ2(IOUT)=AVQ2(IOUT)+(MEQ2)
                     ERQ2(IOUT)=ERQ2(IOUT)+(MEQ2*MEQ2)

                 ENDIF
                 WORK2(ITRAJ)=WORK2(ITRAJ)+(DELTAT*RAMP*TEMP12)
10           CONTINUE 
C 

c            REGIME 3 : RELAXATION AT FINAL CONDITION
             DO 119 I=1,3
                 DO 120 J=1,3
                     FKAPPA(I,J)=FLOW_TMP(I,J)*SR2
120              CONTINUE
119          CONTINUE 

             NFIN=TREL/DELTAT
             IWAIT=0
             IOUT=0

C           Time integration: semi-implicit predictor-corrector scheme 
             DO 20 ITIME=1,NFIN
                 CALL SEMIMP(Q) 
                 IWAIT=IWAIT+1 
                 IF (IWAIT.EQ.NTIME) THEN
                     IWAIT=0
                     IOUT=IOUT+1
                     RSR(IOUT)=FKAPPA(1,1)
                     QMAG=Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3)
                     MEQ3=SQRT(QMAG)
                     AVQ3(IOUT)=AVQ3(IOUT)+(MEQ3)
                     ERQ3(IOUT)=ERQ3(IOUT)+(MEQ3*MEQ3)
                 ENDIF
20           CONTINUE 
             WRITE(112,9) DELTAT,Q(1),Q(2),Q(3)

             QMAG=Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3)
             MEQ3=SQRT(QMAG)

             QAV3=QAV3+(MEQ3)
             QER3=QER3+(MEQ3*MEQ3)
      
             TEMP12=SR2*(Q(1)*Q(1) - Q(2)*Q(2))
             FWAV3=FWAV3+(TEMP12)
             FWER3=FWER3+(TEMP12*TEMP12)

100      CONTINUE 

C        Averages, statistical errors 

         DO 94 I=1,NOUT1
             AVQ1(I)=AVQ1(I)/NTRAJ
             ERQ1(I)=ERQ1(I)/NTRAJ
             ERQ1(I)=(ERQ1(I)-AVQ1(I)*AVQ1(I))/(NTRAJ-1)
             ERQ1(I)=SQRT(ERQ1(I))
             OUTIME1=(NTIME*I*DELTAT)
c        Output of results 
             WRITE(91,4) ESR(I),OUTIME1,AVQ1(I),ERQ1(I)
94       CONTINUE

         DO 35 I=1,NOUT
             AVQ2(I)=AVQ2(I)/NTRAJ
             ERQ2(I)=ERQ2(I)/NTRAJ
             ERQ2(I)=(ERQ2(I)-AVQ2(I)*AVQ2(I))/(NTRAJ-1)
             ERQ2(I)=SQRT(ERQ2(I))

             OUTIME2=OUTIME1+(NTIME*I*DELTAT)
c        Output of results 
             WRITE(7,4) TSR(I),OUTIME2,AVQ2(I),ERQ2(I)
35       CONTINUE


         DO 36 I=1,NOUT3
             AVQ3(I)=AVQ3(I)/NTRAJ
             ERQ3(I)=ERQ3(I)/NTRAJ
             ERQ3(I)=(ERQ3(I)-AVQ3(I)*AVQ3(I))/(NTRAJ-1)
             ERQ3(I)=SQRT(ERQ3(I))
             OUTIME3=OUTIME2+(NTIME*I*DELTAT)
c        Output of results 
             WRITE(8,4) RSR(I),OUTIME3,AVQ3(I),ERQ3(I)
36       CONTINUE

         QAV1=QAV1/NTRAJ
         QER1=QER1/NTRAJ
         QER1=(QER1-QAV1*QAV1)/(NTRAJ-1)
         QER1=SQRT(QER1)
 
         QAV3=QAV3/NTRAJ
         QER3=QER3/NTRAJ
         QER3=(QER3-QAV3*QAV3)/(NTRAJ-1)
         QER3=SQRT(QER3)

         FWAV1=FWAV1/NTRAJ
         FWER1=FWER1/NTRAJ
         FWER1=(FWER1-FWAV1*FWAV1)/(NTRAJ-1)
         FWER1=SQRT(FWER1)
        
         FWAV3=FWAV3/NTRAJ
         FWER3=FWER3/NTRAJ
         FWER3=(FWER3-FWAV3*FWAV3)/(NTRAJ-1)
         FWER3=SQRT(FWER3)
        

1000  CONTINUE 

      AVW=0.D0
      ERW=0.D0
      AW=0.D0
      VW=0.D0

      DO 92 ITRAJ=1,NTRAJ
          WORK2(ITRAJ)=FWAV3-FWAV1-WORK2(ITRAJ)
          WRITE(11,29) ITRAJ,WORK2(ITRAJ)
          MEW=EXP(-WORK2(ITRAJ))
          AVW=AVW+(MEW)
          ERW=ERW+(MEW*MEW)
          AW=AW+(WORK2(ITRAJ))
          VW=VW+(WORK2(ITRAJ)*WORK2(ITRAJ))

92    CONTINUE

      AVW=AVW/NTRAJ
      ERW=ERW/NTRAJ
      ERW=(ERW-AVW*AVW)/(NTRAJ-1)
      ERW=SQRT(ERW)

      AW=AW/NTRAJ
      VW=VW/NTRAJ
      VW=(VW-AW*AW)/(NTRAJ-1)
      VW=SQRT(VW)


      DELF=-LOG(AVW)
      ER_DELF=ERW/AVW

      WDIS=AW-DELF
      ER_WDIS=(VW*VW)+(ER_DELF*ER_DELF)
      ER_WDIS=SQRT(ER_WDIS)

      CALL CPU_TIME(ENDTIME)

      WRITE(*,*) "SHEAR RATE TO BE WRITTEN : ",SR2
      WRITE(112,23) ISEED

      BACKSPACE (UNIT=1)
      WRITE (1,3) Z,RMU,B,SR1,SR2,CIV,ETA_S,B_RAD,THI,INPAR,NFLOW,
     &(ENDTIME-STARTTIME)
3     FORMAT(F5.1,2X,F4.1,4X,F8.1,6X,F8.2,6X,F8.2,6X,F5.2,
     &4X,F5.2,4X,F5.2,4X,F5.2,4X,I1,4X,I1,4X,F10.1)
      WRITE(1,*) "     "
      WRITE(1,'(23A)')'## Z',TAB,'RMU',TAB,'B',TAB,'SR1',TAB,'SR2',
     &TAB,'KIV',TAB,'ETA_S',TAB,'B_RAD',TAB,'THI',TAB,'INPAR',
     &TAB,'NFLOW',TAB,'EXEC.TIME ##'
      WRITE(1,'(1A,F5.2)')'## IV PARAMETER, \epsilon     : ',E
      WRITE(1,'(1A,F5.2)')'## HI PARAMETER, h*           : ',H0

      WRITE(10,'(1A,I10)') 'NO. OF TRAJECTORIES   : ',NTRAJ
      WRITE(10,'(1A,F7.5)')'DELTAT                :       ',DELTAT
      WRITE(10,'(1A,F7.5)')'RAMP RATE             :       ',RAMP
      WRITE(10,'(1A,F7.2)')'RAMP TIME             :       ',TMAX
      WRITE(10,'(1A,F8.2)')'Pe_A                  :       ',SR1
      WRITE(10,'(1A,F8.2)')'RELAXATION IN REG.1   :       ',TEQB
      WRITE(10,81) '<|Q|> IN REGIME 1     : ',QAV1,'+/-',QER1
      WRITE(10,81) '<\Chi> IN REGIME 1    : ',FWAV1,'+/-',FWER1
      WRITE(10,'(1A,F8.2)')'Pe_B                  :       ',SR2
      WRITE(10,'(1A,F8.2)')'RELAXATION IN REG.3   :       ',TREL
      WRITE(10,81) '<|Q|> IN REGIME 3     : ',QAV3,'+/-',QER3
      WRITE(10,81) '<\Chi> IN REGIME 3    : ',FWAV3,'+/-',FWER3
      WRITE(10,*) "          "
      WRITE(10,81) '\Delta F              : ',DELF,'+/-',ER_DELF
      WRITE(10,81) '<W>                   : ',AW,'+/-',VW
      WRITE(10,81) '<W_dis>               : ',WDIS,'+/-',ER_WDIS

22    FORMAT(F5.2,2X,F5.2,2X,F5.2,2X,F8.2,2X,F8.2,6X,F7.5,
     &4X,F12.5,4X,F12.5,4X,F12.5,4X,F12.5,4X,F12.5,4X,F12.5)

      WRITE(10,*) "          "
      WRITE(10,'(23A)')'## E',TAB,'KIV',TAB,'ETA_S',TAB,'SR1',TAB,'SR2',
     &TAB,'RAMP',TAB,'W_DIS',TAB,'ER_WDIS',TAB,'DEL_F',TAB,'ER_DEL_F',
     &TAB,'AV_WORK',TAB,'ER_AV_WORK ##'
      WRITE(10,*) "          "
      WRITE(10,22) E,CIV,ETA_S,SR1,SR2,RAMP,WDIS,ER_WDIS,DELF,ER_DELF,
     &AW,VW
    

      CLOSE (UNIT=1)
      close (unit=2) 
c      close (unit=3)
      close (UNIT=91)
      CLOSE (UNIT=7)
      close (UNIT=8)
      CLOSE (UNIT=112)
      CLOSE (UNIT=10)
      CLOSE (UNIT=11) 

1100  STOP
C 
      open (unit=15, file='q2.dat')
      open (unit=2, file='tstep.dat')
      open (unit=1, file='inp.dat')

      open (unit=40, file='q2x.dat',status='UNKNOWN')
      READ (1,*) Z,RMU,B,SR,E,B_RAD,THI 
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
         READ(15,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)

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
      WRITE (40,2) SR, E , B_RAD, OUTIME(I), YOPT, SIGOPT, NOPT, NDUOPT,
     &' of ', NTIWID, ALIOPT,' +- ',VLIOPT
33    CONTINUE

      close (unit=15)
      close (unit=40)
      close (unit=1)
      close (unit=2)

C1100  STOP 


      END 
 
c
c

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
c
      SUBROUTINE SEMIMP(Q)
C     Time step: semi-implicit predictor-corrector scheme for FENE+IV
C     model
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
      RED=R/G
      Q(1)=RED*Q(1)
      Q(2)=RED*Q(2)
      Q(3)=RED*Q(3)
      RETURN
      END

c
C 
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

