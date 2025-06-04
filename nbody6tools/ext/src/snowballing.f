      PROGRAM main
      IMPLICIT NONE

      END PROGRAM

      SUBROUTINE SNOWBALLING_METHOD(NSTARS,NCL,STARS,CLUSTERS
     & ,RCORES,SCREEN,CLSTARS,G)
      IMPLICIT NONE
      !INPUTS
      INTEGER,intent(in):: NSTARS,NCL !Number of stars, Number of clumps
      DOUBLE PRECISION,intent(in) :: RCORES(NCL),STARS(8,NSTARS)
      DOUBLE PRECISION,intent(in) :: CLUSTERS(3,NCL)
      LOGICAL,intent(in):: SCREEN  !verbose
      DOUBLE PRECISION,intent(in) :: G
      !OUTPUTS
      INTEGER,intent(OUT):: CLSTARS(NCL,NSTARS) !cluster membership in flags
      !INTERNAL
      !not anoutput for python

      INTEGER I,CL,K,NIN,NBOUND,TOL
      DOUBLE PRECISION M,MI,R,CM(3),CV(3),VSIGMAS(NCL),VCORES(3,NCL)
      DOUBLE PRECISION RAD(NSTARS),RMAX,RTOCL(NCL,NSTARS)
      TOL=0
      RMAX=1.E20
      !K iterates dimension 
      !I iterates particle
      !CL iterates cluster
      IF (SCREEN) THEN
              PRINT*, "UPDATED 6/4/23"
      END IF
      !Preparing data to the subroutines
      DO CL=1,NCL
            IF (SCREEN) THEN 
               PRINT*, "Cluster # ",CL
               PRINT*, "Position # ",(CLUSTERS(K,CL),K=1,3)
               PRINT *,"Initial guess :", RCORES(CL)
            END IF
            MI=0
            CM(1:3)=(/0.0,0.0,0.0/)
            CV(1:3)=(/0.0,0.0,0.0/)
            NIN=0
            VSIGMAS(CL) = 0.0
            DO I=1,NSTARS
               R=0.0
               M=STARS(7,I)
               DO K=1,3
                  R=R+(STARS(K,I)-CLUSTERS(K,CL))**2
               END DO
               R=SQRT(R)
               RAD(I)=R
               RTOCL(CL,I)=R
               IF (R.LE.RCORES(CL)) THEN
                       CLSTARS(CL,I)=1
                       DO K=1,3
                          CM(K)=CM(K)+STARS(K,I)*M
                          CV(K)=CV(K)+STARS(K+3,I)*M
                       END DO
                       MI=MI+M
                       NIN = NIN + 1
               ELSE
                       CLSTARS(CL,I)=0
               END IF
            END DO
            DO K=1,3
               CM(K)=CM(K)/MI
               CV(K)=CV(K)/MI
            END DO
            !PRINT * ,"CM", (CM(K),k=1,3)
            !PRINT * ,"CV", (CV(K),k=1,3)

            IF (SCREEN) THEN
                    PRINT*, "Mass inside Rcore :", MI
                    PRINT*, "Part. inside Rcore :", NIN
                    PRINT*, "mean mass :", MI/NIN, STARS(7,1)
            END IF

            CALL ELIMINATION_STEP(NSTARS,STARS,RCORES(CL),CLSTARS(CL,:),
     &                 VSIGMAS(CL),VCORES(:,CL),NIN,NBOUND,RAD,SCREEN,G)

            CALL SB(NSTARS,STARS,RMAX,CLSTARS(CL,:),VSIGMAS(CL),
     &               VCORES(:,CL),RAD,SCREEN,TOL,G)

      END DO
      !CALL SB_SIME(NSTARS,NCL,STARS,RMAX,CLSTARS,VSIGMAS,VCORES,RTOCL,
      !&                                                   SCREEN,TOL,G)


      END

      SUBROUTINE ELIMINATION_STEP(NPART,STARS,RCORE,FLAG,VSIGMA,VCORE,
     & NIN,NBOUND,RAD,SCREEN,G)
      IMPLICIT NONE
      ! First step of the "Snowballing Method"
      ! This subroutine modifies the FLAG array.
      ! FLAG is an array with 1 and 0 for a selected set of particles.
      ! For a given radius RCORE it calculates if each particle is bound
      ! to the others and if not it flag the particle as 0 (unbound).
      ! This subroutine just eliminate particles (correcting the average
      ! velocity of the core every time) it DO NOT add particles to the
      ! core.
      !
      ! NPART          : Number of particles
      ! STARS(8,NPART) : Array containing position, velocity and mass of
      !                  the particles in format: x,y,z,vx,vy,vz,m,epot
      !                  positions in parsec, velocities in km/s,
      !                  masses in MSun, epot in Msun*kms**2
      ! NOTE: Positions are supposed to be corrected in order to have
      ! the center of the cluster at (0,0,0) 
      !
      ! updates: 
      ! 2/06/2018 : parallelized loops with Openmp
      ! 25/04/2020 : It requires a value for G 


      INTEGER NPART,FLAG(NPART),NIN,NBOUND
      DOUBLE PRECISION STARS(8,NPART),RCORE,RAD(NPART),G
      LOGICAL SCREEN

      INTEGER ITERFLAG,I,J,K,NOLD,ITER,TOL,ITERMAX,NB(NPART),CNT
      INTEGER DN,N,FORBIDEN(NPART),NRM,NREMOVE, IMAX(10),NRMMAX
      DOUBLE PRECISION MASSIN,KE,POT(NPART),IPOT(NPART),DX,DY,DZ,GR
      DOUBLE PRECISION SEP,V2,TEMP,MBOUNDI,COREPOT,COREKIN,MBOUND
      DOUBLE PRECISION VCORE(1:3),FLAG0(NPART),EBMAX(10),IKIN(NPART)
      DOUBLE PRECISION VPART(NPART),VSIGMA,VMEAN
      !DOUBLE PRECISION W(NPART),WSUM 
      !G=4.3024e-3
      DO I=1,NPART 
         FORBIDEN(I) = 0 !when excluding loose particles
         POT(I) = STARS(8,I)
         IF (RAD(I).GT.RCORE) THEN 
                 FLAG(I) = 0
         ELSE
                 FLAG(I) = 1
         ENDIF
      ENDDO

      NRM = 0
      NREMOVE = 1
      NRMMAX = 100
      TOL=1
      ITERFLAG=TOL+10
      ITER=0
      NBOUND=0
      ITERMAX=20 + NRMMAX
      DN=1
      ITER = 0
      N=0
10    CONTINUE
      IF (ABS(DN).GT.0) THEN
        ITER = ITER+1
        NOLD = N
        VMEAN = 1
        MBOUNDI=1.0
        vcore=(/0,0,0/)
        ! First we need a good frame velocity to start with. 
        ! note that the potential is oly calculated once in this
        ! subroutine, then POT(I) is only started once
        DO I=1,NPART 
           IF  (FLAG(I).EQ.0) CYCLE
           DO K = 1,3
              VCORE(K)=VCORE(K)+STARS(3+K,I)*STARS(7,I)!*W(I)
           END DO
           MBOUNDI=MBOUNDI+STARS(7,I)
        END DO
        IF (MBOUNDI.GT.0) THEN
                TEMP=1./MBOUNDI
        ELSE
                TEMP=0
        END IF

        DO K = 1,3
           VCORE(K)=VCORE(K)*TEMP
           VMEAN = VMEAN + VCORE(K)*VCORE(K)
        END DO
        VMEAN = SQRT(VMEAN)

        ! initialize  zero velocity and calculate standard deviation
        CNT = 0
        VSIGMA = 0.0
        DO I=1,NPART
           VPART(I) =(STARS(4,I)-VCORE(1))**2 + (STARS(5,I)-VCORE(2))**2
     &          +    (STARS(6,I)-VCORE(3))**2
           VPART(I) = SQRT(VPART(I))
           IF (FLAG(I).EQ.0) CYCLE
           VSIGMA  = VSIGMA + (VPART(I)-VMEAN)**2
           CNT=CNT+1
        END DO
        VSIGMA = SQRT(VSIGMA/(CNT-1))
        N=0
        DO I=1,NPART
                IF (VPART(I).GT.10*VSIGMA ) THEN
                        FLAG(I) = 0
                ELSE        
                        N = N + 1
                END IF
        ENDDO
        DN = ABS(NOLD-N)
        GOTO 10
      END IF

      IF (SCREEN)  print*,"VMEAN:VSIGMA", VMEAN,VSIGMA,CNT
      IF (SCREEN)  print*,"STARTING VELOCITY", VCORE(:)

      ! Start the elimination step
      ITER=0
80    CONTINUE
      IF ((ABS(ITERFLAG).GT.TOL).AND.(ITER.LT.ITERMAX)) THEN
         ITER=ITER+1
         !ITERFLAG=NOLD-NBOUND
         NOLD=NBOUND
         NIN=0
         NBOUND=0
         MASSIN=0
         MBOUND=0
         !WSUM=0.0
         MBOUNDI=0.0
         DO I=1,NPART
            FLAG0(I) = FLAG(I)
            NB(I) = 0
         END DO

         IF (ITER.EQ.1) THEN !only do it once
         DO I=1,NPART-1
            IF (RAD(I).GT.RCORE) CYCLE
            NIN=NIN+1
            MASSIN=MASSIN+STARS(7,I)

            DO J=I+1,NPART
               IF (I.EQ.J) CYCLE
               IF (RAD(J).GT.RCORE) CYCLE
               IF (FLAG0(J).eq.0) CYCLE

               DX=(STARS(1,I)-STARS(1,J))
               DY=(STARS(2,I)-STARS(2,J))
               DZ=(STARS(3,I)-STARS(3,J))
               SEP=SQRT(DX*DX + DY*DY + DZ*DZ)
               GR=G/SEP
               POT(I)=POT(I)-GR*STARS(7,J)
               POT(J)=POT(J)-GR*STARS(7,I)

            END DO
         END DO
         END IF

         DO I=1,NPART
            V2=0
            DO K = 1,3
                V2 = V2 + (STARS(3+K,I) - VCORE(K))**2
            END DO
            VPART(I) = SQRT(V2)
            
            IKIN(I) = 0.5*STARS(7,I)*V2
            IF (IKIN(I)+POT(I)*STARS(7,I).GE.0) THEN
               FLAG(I) = 0
            ELSE
               IF (FORBIDEN(I).EQ.1) THEN
                       FLAG(I) = 0 
               ELSE
                       NBOUND=NBOUND+1
                       FLAG(I)=1
               END IF
            END IF
         ENDDO
         mbound=0
         vcore=(/0,0,0/)
         TEMP = 0
         COREPOT = 0.0
         COREKIN = 0.0
         DO K=1,10 !10 maximum removals at each step 
                 EBMAX(K) = -1.0e30
                 IMAX(K) = NPART*10!value used later to check if changed
         END DO
         DO I=1,NPART
           IF (FLAG(I).EQ.1) THEN
            IPOT(I) = 0
            DO J=1,NPART
               IF (I.EQ.J) CYCLE
               IF (RAD(J).GT.RCORE) CYCLE
               IF (FLAG(J).EQ.0) CYCLE
               DX=(STARS(1,I)-STARS(1,J))
               DY=(STARS(2,I)-STARS(2,J))
               DZ=(STARS(3,I)-STARS(3,J))
               SEP=SQRT(DX*DX + DY*DY + DZ*DZ)
               GR=G/SEP
               IPOT(I)=IPOT(I)-GR*STARS(7,J)*STARS(7,I)
            END DO
           END IF
           
           IF (FORBIDEN(I).EQ.0.AND.FLAG(I).EQ.1) THEN
                   DO K=1,NREMOVE
                   IF (IPOT(I)+IKIN(I).GT.EBMAX(K)) THEN
                           IMAX(K) = I
                           EBMAX(K) = IPOT(I) + IKIN(I)
                           EXIT
                   END IF
                   END DO
           END IF


         IF (FLAG(I).EQ.1) THEN
               COREPOT = COREPOT + IPOT(I)
               COREKIN = COREKIN + IKIN(I)
               MBOUND=MBOUND+STARS(7,I)

               IF (VPART(I).LE.10*VSIGMA) THEN
                 TEMP = TEMP + STARS(7,I)
                 DO K=1,3
                    VCORE(K)=VCORE(K)+STARS(3+K,I)*STARS(7,I)
                 END DO
               ELSE
                 FLAG(I) = 0
               END IF
         END IF
         END DO


         IF (TEMP.GT.0) THEN 
                TEMP=1./TEMP
         ELSE
                TEMP=0
         END IF
         DO K=1,3
               VCORE(K)=VCORE(K)*TEMP
         END DO
         
         !stop
         ITERFLAG=NOLD-NBOUND
         IF (SCREEN) THEN
            !PRINT*,"(ES)nbound,nin",NBOUND,NIN
            !PRINT*,"-----ITER:",ITER,"-------"
            !PRINT*,"(ES)Mass inside Rcore:",MASSIN
            !PRINT*,"(ES)Mass bound in it.:",MBOUND
            !PRINT*,"-------------------------"
           PRINT*, "VCORE:", VCORE(:),"KIN,POT,Q",COREKIN,COREPOT,
     &              abs(COREKIN/COREPOT)
           PRINT*, ITER, NBOUND , MBOUND , NIN, MASSIN," ES",ITERFLAG
         END IF
         GOTO 80
      END IF

      IF ( ABS(COREKIN/COREPOT).GT.1.0 ) THEN
              NRM = NRM +1
             IF (SCREEN) THEN
                 print*,'Q core > 1, removing a loose particle, NRM',NRM
                 print*,COREKIN/COREPOT,NBOUND,IMAX(:),EBMAX(:)
             END IF
             DO K=1,NREMOVE
                     IF (IMAX(K).LE.NPART+1) THEN
                !print*, 'forbidding: k',K,IMAX(K)!,FLAG(IMAX(K))
                       FORBIDEN(IMAX(K)) = 1
                       FLAG(IMAX(K)) = 0 
                     END IF
             END DO
             ITERFLAG = 2*TOL
             !IF (SCREEN) THEN
             !        DO I=1,NPART 
             !        IF (FORBIDEN(I).EQ.1) THEN
             !                print*,'FB',I,EBMAX(:)
             !        END IF
             !        END DO
             !END IF
             if (NRM.LT.NRMMAX) THEN !max iterations
                !ITER = 2
                GO TO 80
             END IF
      END IF
      END

      SUBROUTINE SB_SIME_OLD(NPART,NCL,STARS,RMAX,CFLAGS,VSIGMA,VCORE,
     & R,SCREEN,TOL,G)
      IMPLICIT NONE

      INTEGER NPART,NCL,CFLAGS(NCL,NPART),INCORE(NPART)
      DOUBLE PRECISION STARS(8,NPART),RCORE,RMAX,MBOUND(NCL)
      DOUBLE PRECISION R(NCL,NPART),G,VSIGMA(NCL),VCORE(3,NCL)
      LOGICAL SCREEN

      INTEGER I,J,K,BMIN,L,ITERFLAG,NBOUND(NCL),NBOUNDOLD(NCL),CL
      INTEGER ITERMAX,TOL, ITEMP
      DOUBLE PRECISION SEP,KE,POT,IVEL(4,NPART), TEMP, MCOUNT(NCL)
      DOUBLE PRECISION PTOT(NCL),KTOT(NCL)
      DOUBLE PRECISION AVGVX(NCL),AVGVY(NCL),AVGVZ(NCL),FR,BE(NCL,NPART)
      DOUBLE PRECISION DX,DY,DZ, VPART(NPART), VK
      LOGICAL ALLPARTS

      !G=4.3024e-3
      !G=4317.15
      !G=4302.008
      ITERMAX=40
      RCORE=RMAX
      FR=1
      L=0
      ITERFLAG=TOL+1000
      ITEMP = 0
      DO I=1,NPART
         INCORE(I)=0
         DO CL=1,NCL
            IF (CFLAGS(CL,I).EQ.1) THEN
               INCORE(I)=1
            END IF
         END DO
      ENDDO
      !! DEPRECATING ALLPARTS, allways false now
      IF (NCL.EQ.1) THEN
              ALLPARTS=.TRUE.
      ELSE
              ALLPARTS=.FALSE.
      END IF

90    CONTINUE
      IF (ABS(ITERFLAG).GT.TOL.AND.L.LE.ITERMAX) THEN
         L=L+1
         IF (L.GT.1.OR.NCL.EQ.1) THEN
            RCORE=RMAX
         END IF
         DO CL=1,NCL
            NBOUNDOLD(CL)=NBOUND(CL)
         END DO
         ! Obtain mean velocity with currently bound particles
         DO CL=1,NCL
            MBOUND(CL)= 0.
            MCOUNT(CL)= 0.
            PTOT(CL) = 0.
            KTOT(CL) = 0.
            AVGVX(CL) = 0.
            AVGVY(CL) = 0.
            AVGVZ(CL) = 0.
            DO I=1,NPART
               IF (CFLAGS(CL,I).EQ.1) THEN 
                  MBOUND(CL)=MBOUND(CL)+STARS(7,I)
                  VPART(I) = 0.0
                  DO K=1,3 
                        VK = STARS(3+K,I) - VCORE(K,CL)
                        VPART(I) = VPART(I) + VK*VK
                  END DO
                  VPART(I) = DSQRT(VPART(I))
                  IF ( VPART(I).LT.(2*VSIGMA(CL)) ) THEN
                          MCOUNT(CL) = MCOUNT(CL) + STARS(7,I)
                          AVGVX(CL) = AVGVX(CL)+STARS(4,I)*STARS(7,I)
                          AVGVY(CL) = AVGVY(CL)+STARS(5,I)*STARS(7,I)
                          AVGVZ(CL) = AVGVZ(CL)+STARS(6,I)*STARS(7,I)
                  END IF
               END IF
            END DO
            IF (MCOUNT(CL).GT.0) THEN
               VCORE(1,CL)=AVGVX(CL)/MCOUNT(CL)
               VCORE(2,CL)=AVGVY(CL)/MCOUNT(CL)
               VCORE(3,CL)=AVGVZ(CL)/MCOUNT(CL)
            END IF
         END DO
      !Now obtain the binding energy for each particle to each core

        DO CL=1,NCL
            DO I=1,NPART
               !print*,CL,I, R(CL,I),R(CL,I)/RMAX
               IF (R(CL,I).GT.RMAX) CYCLE !OPTIMIZATION
               IF (.NOT.ALLPARTS.AND.INCORE(I).EQ.1) CYCLE
               IVEL(4,I) = 0.0
               DO K=1,3 
                       IVEL(K,I)=STARS(3+K,I)-VCORE(K,CL)
                       IVEL(4,I) = IVEL(4,I) + IVEL(K,I)*IVEL(K,I) 
               END DO
               POT=0.0
               DO J=1,NPART
                  IF (J.EQ.I) CYCLE
                  !IF (CFLAGS(CL,J).EQ.0) CYCLE
                  IF (R(CL,I).GT.RCORE) CYCLE
                  DX=STARS(1,I)-STARS(1,J)
                  DY=STARS(2,I)-STARS(2,J)
                  DZ=STARS(3,I)-STARS(3,J)
                  SEP=SQRT(DX*DX+DY*DY+DZ*DZ)
                  POT=POT-STARS(7,J)/SEP
               END DO
               KE=0.5*STARS(7,I)*IVEL(4,I)
               BE(CL,I)=KE + (G*POT)*STARS(7,I) + STARS(8,I)
               PTOT(CL) = PTOT(CL) + BE(CL,I)
               KTOT(CL) = KTOT(CL) + KE
            END DO
        END DO

        !Now, Check which core a particle is more bound
          DO I=1,NPART
             IF (.NOT.ALLPARTS.AND.INCORE(I).EQ.1) CYCLE
             BMIN=1 !index of the cluster where more bound
             DO CL=1,NCL
                IF (BE(CL,I).LT.BE(BMIN,I)) THEN
                   BMIN=CL
                END IF
             END DO

             ! mark the flag of the proper clump and unmark others
             IF (BE(BMIN,I).LT.0) THEN
                CFLAGS(BMIN,I)=1 
                DO CL=1,NCL
                   IF (CL.EQ.BMIN) CYCLE
                   CFLAGS(CL,I)=0
                END DO
             ELSE ! BE > 0 not bound to any cluster
                DO CL=1,NCL
                   CFLAGS(CL,I)=0
                END DO
             END IF
          END DO !parts
          ITERFLAG=0
          DO CL=1,NCL
             NBOUND(CL)=0
             DO I=1,NPART
                NBOUND(CL)=NBOUND(CL)+CFLAGS(CL,I)
             END DO
             ITERFLAG=ITERFLAG + NBOUND(CL) - NBOUNDOLD(CL)
          END DO
          IF (SCREEN) THEN
             DO CL = 1,NCL
                     PRINT *,L,NBOUND(CL)," SB (Rmax:",RMAX,")"
                     PRINT*,'VCORE', VCORE(:,CL)
                     PRINT*,'POT,KIN,Q', PTOT(CL),KTOT(CL),
     &                                   KTOT(CL)/ABS(PTOT(CL))
             END DO
          ENDIF
c          IF (ABS(ITERFLAG).LE.TOL.AND.(.NOT.ALLPARTS))
c     &    THEN !If no no changes then consider all particles now
c !ALLPARTS=.TRUE. !Deprecating
c             ALLPARTS=.FALSE.
c             ITERFLAG=TOL+1000
c          ENDIF

          GOTO 90
      ENDIF
      END

      SUBROUTINE SB(NPART,STARS,RMAX,FLAGS,VSIGMA,VCORE,
     & R,SCREEN,TOL,G)
      IMPLICIT NONE

      INTEGER NPART,FLAGS(NPART),INCORE(NPART)
      DOUBLE PRECISION STARS(8,NPART),RCORE,RMAX,MBOUND
      DOUBLE PRECISION R(NPART),G,VSIGMA,VCORE(3)
      LOGICAL SCREEN

      INTEGER I,J,K,BMIN,L,ITERFLAG,NBOUND,NBOUNDOLD,CL
      INTEGER ITERMAX,TOL, ITEMP
      DOUBLE PRECISION SEP,KE,POT,IVEL(4,NPART), TEMP, MCOUNT
      DOUBLE PRECISION PTOT,KTOT,MINSEP
      DOUBLE PRECISION AVGVX,AVGVY,AVGVZ,FR,BE(NPART)
      DOUBLE PRECISION DX,DY,DZ, VPART(NPART), VK
      LOGICAL ALLPARTS

      !G=4.3024e-3
      !G=4317.15
      !G=4302.008
      MINSEP = 0.0
      ITERMAX=40
      RCORE=RMAX
      FR=1
      L=0
      ITERFLAG=TOL+1000
      ITEMP = 0
      DO I=1,NPART
         INCORE(I)=0
         IF (FLAGS(I).EQ.1) THEN
            INCORE(I)=1
         END IF
      ENDDO

90    CONTINUE
      IF (ABS(ITERFLAG).GT.TOL.AND.L.LE.ITERMAX) THEN
         L=L+1
         NBOUNDOLD=NBOUND
         ! Obtain mean velocity with currently bound particles
         MBOUND= 0.
         MCOUNT= 0.
         PTOT = 0.
         KTOT = 0.
         AVGVX = 0.
         AVGVY = 0.
         AVGVZ = 0.
         DO I=1,NPART
            IF (FLAGS(I).EQ.1) THEN 
               MBOUND=MBOUND+STARS(7,I)
               VPART(I) = 0.0
               DO K=1,3 
                     VK = STARS(3+K,I) - VCORE(K)
                     VPART(I) = VPART(I) + VK*VK
               END DO
               VPART(I) = DSQRT(VPART(I))
               !IF ( VPART(I).LT.(10*VSIGMA) ) THEN
               MCOUNT = MCOUNT + STARS(7,I)
               AVGVX = AVGVX+STARS(4,I)*STARS(7,I)
               AVGVY = AVGVY+STARS(5,I)*STARS(7,I)
               AVGVZ = AVGVZ+STARS(6,I)*STARS(7,I)
               !END IF
            END IF
         END DO
         IF (MCOUNT.GT.0) THEN
            VCORE(1)=AVGVX/MCOUNT
            VCORE(2)=AVGVY/MCOUNT
            VCORE(3)=AVGVZ/MCOUNT
         END IF
      !print*, "avx avy avz mcount: ", AVGVX(CL), AVGVY(CL), AVGVZ(CL),
      !&   MCOUNT(CL)
      !Now obtain the binding energy using only bound stars

        
         DO I=1,NPART
            IF (R(I).GT.RMAX) CYCLE !OPTIMIZATION
            !IF (.NOT.ALLPARTS.AND.INCORE(I).EQ.1) CYCLE
            IVEL(4,I) = 0.0
            DO K=1,3 
                    IVEL(K,I)=STARS(3+K,I)-VCORE(K)
                    IVEL(4,I) = IVEL(4,I) + IVEL(K,I)*IVEL(K,I) 
            END DO
            POT=0.0
            DO J=1,NPART
               IF (J.EQ.I) CYCLE
               IF (FLAGS(J).EQ.0) CYCLE
               IF (R(I).GT.RCORE) CYCLE
               DX=STARS(1,I)-STARS(1,J)
               DY=STARS(2,I)-STARS(2,J)
               DZ=STARS(3,I)-STARS(3,J)
               SEP=SQRT(DX*DX+DY*DY+DZ*DZ)
               POT = POT - STARS(7,J)/SEP
               CL = CL+1
            END DO
            KE = 0.5*STARS(7,I)*IVEL(4,I)
            POT = G*POT*STARS(7,I) + STARS(8,I)
            BE(I) = KE + POT
            !print*, I,FLAGS(I),KE,POT,BE(I)
            IF (BE(I).LT.0.OR.INCORE(I).EQ.1) THEN
                    FLAGS(I) = 1
                    PTOT = PTOT + POT
                    KTOT = KTOT + KE
            ELSE
                    FLAGS(I) = 0
            END IF
         END DO

         ITERFLAG=0
         NBOUND=0
         DO I=1,NPART
            NBOUND=NBOUND+FLAGS(I)
         END DO
         ITERFLAG=ITERFLAG + NBOUND - NBOUNDOLD
         IF (SCREEN) THEN
            PRINT *,L,NBOUND," SB (Rmax:",RMAX,")"
            PRINT*,'VCORE', VCORE(:)
            PRINT*,'POT,KIN,Q', PTOT,KTOT,
     &                          KTOT/ABS(PTOT)
         ENDIF

         GOTO 90
      END IF
      END


