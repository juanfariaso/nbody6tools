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

      END DO
      CALL SB_SIME(NSTARS,NCL,STARS,RMAX,CLSTARS,VSIGMAS,VCORES,RTOCL,
     &                                                     SCREEN,TOL,G)

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
      DOUBLE PRECISION MASSIN,KE,POT(NPART),DX,DY,DZ,GR,MBOUND
      DOUBLE PRECISION SEP,V2,TEMP,MBOUNDI
      DOUBLE PRECISION VCORE(1:3),FLAG0(NPART)
      DOUBLE PRECISION VPART(NPART),VSIGMA,VMEAN
      !DOUBLE PRECISION W(NPART),WSUM 
      !G=4.3024e-3

      TOL=1
      ITERFLAG=TOL+10
      ITER=0
      NBOUND=0
      ITERMAX=20
      VMEAN = 0
      !WSUM=0
      MBOUNDI=0.0
      vcore=(/0,0,0/)
      ! First we need a good frame velocity to start with. 
      ! note that the potential is oly calculated once in this
      ! subroutine, then POT(I) is only started once
      DO I=1,NPART 
         POT(I) = STARS(8,I)
         !IF (RAD(I).LT.0.5*RCORE) THEN
         !        W(I) = 1.0
         !ELSE
         !        W(I) = 1.0 - ((RAD(I)-0.1*RCORE)/RAD(I))
         !END IF

         IF (RAD(I).GT.RCORE) CYCLE
         DO K = 1,3
            VCORE(K)=VCORE(K)+STARS(3+K,I)*STARS(7,I)!*W(I)
         END DO
         MBOUNDI=MBOUNDI+STARS(7,I)
         !WSUM=WSUM+W(I)
      END DO
      IF (MBOUNDI.GT.0) THEN
              !TEMP=1./MBOUNDI/WSUM
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
         !DO K = 1,3
         !    VPART(I) = VPART(I) + (STARS(3+K,I) - VCORE(K))**2
         !END DO
         VPART(I) = (STARS(4,I)-VCORE(1))**2 + (STARS(5,I)-VCORE(2))**2
     &        +     (STARS(6,I)-VCORE(3))**2
         VPART(I) = SQRT(VPART(I))
         IF (RAD(I).GT.RCORE) CYCLE
         !VSIGMA  = VSIGMA + (VPART(I)-VMEAN)**2/W(I)
         VSIGMA  = VSIGMA + (VPART(I)-VMEAN)**2
         CNT=CNT+1
         !print*, I , VSIGMA, VPART(I), CNT
         !if (I.GT.14) STOP
      END DO
      !VSIGMA = SQRT(CNT*VSIGMA/((CNT-1)/WSUM))
      VSIGMA = SQRT(VSIGMA/(CNT-1))
      IF (SCREEN)  print*,"VMEAN:VSIGMA", VMEAN,VSIGMA,CNT
      IF (SCREEN)  print*,"STARTING VELOCITY", VCORE(:)

      ! Start the elimination step
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
         DO I=1,NPART
            IF (RAD(I).GT.RCORE) CYCLE
            NIN=NIN+1
            MASSIN=MASSIN+STARS(7,I)
            DO J=1,NPART
               IF (I.EQ.J) CYCLE
               IF (RAD(J).GT.RCORE) CYCLE
               IF (FLAG0(J).eq.0) CYCLE
               IF (J.GT.I) THEN 
               !IF (J.GT.I) THEN ! do it allways
               !if (i.eq.1.and.j.eq.2) then
               !   print *, "measuring"
               !end if
               DX=(STARS(1,I)-STARS(1,J))
               DY=(STARS(2,I)-STARS(2,J))
               DZ=(STARS(3,I)-STARS(3,J))
               SEP=SQRT(DX*DX + DY*DY + DZ*DZ)
               GR=G/SEP
               POT(I)=POT(I)-GR*STARS(7,J)
               POT(J)=POT(J)-GR*STARS(7,I)
               END IF
            END DO
         END DO
         END IF

         DO I=1,NPART
            V2=0
            DO K = 1,3
                V2 = V2 + (STARS(3+K,I) - VCORE(K))**2
            END DO
            VPART(I) = SQRT(V2)
            KE=0.5*STARS(7,I)*V2
            IF (KE+POT(I)*STARS(7,I).GT.0) THEN
               FLAG(I)=0
            ELSE
               NBOUND=NBOUND+1
               !MBOUND=MBOUND+STARS(7,I)
               FLAG(I)=1
            END IF
         ENDDO
         mbound=0
         vcore=(/0,0,0/)
         TEMP = 0
         DO I=1,NPART
         IF (FLAG(I).EQ.1) THEN
               MBOUND=MBOUND+STARS(7,I)
               IF (VPART(I).LE.2*VSIGMA) THEN
                 TEMP = TEMP + STARS(7,I)
                 DO K=1,3
                    VCORE(K)=VCORE(K)+STARS(3+K,I)*STARS(7,I)
                 END DO
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
           PRINT*, "VCORE:", VCORE(:)
           PRINT*, ITER, NBOUND , MBOUND , NIN, MASSIN," ES",ITERFLAG
         END IF
         GOTO 80
      END IF
      !DO I=1,NPART
      !      if (flag(I).eq.1) then
      !            write(6,"(I5)",advance="no") I
      !      end if
      !end do
      !print *,""
      END

      SUBROUTINE SB_SIME(NPART,NCL,STARS,RMAX,CFLAGS,VSIGMA,VCORE,
     & R,SCREEN,TOL,G)
      IMPLICIT NONE

      INTEGER NPART,NCL,CFLAGS(NCL,NPART),INCORE(NPART)
      DOUBLE PRECISION STARS(8,NPART),RCORE,RMAX,MBOUND(NCL)
      DOUBLE PRECISION R(NCL,NPART),G,VSIGMA(NCL),VCORE(3,NCL)
      LOGICAL SCREEN

      INTEGER I,J,K,BMIN,L,ITERFLAG,NBOUND(NCL),NBOUNDOLD(NCL),CL
      INTEGER ITERMAX,TOL, ITEMP
      DOUBLE PRECISION SEP,KE,POT,IVEL(4,NPART), TEMP, MCOUNT(NCL)
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
      !print*, "avx avy avz mcount: ", AVGVX(CL), AVGVY(CL), AVGVZ(CL),
      !&   MCOUNT(CL)
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
                     PRINT *,L,NBOUND(CL) ," SB"
                     PRINT*,'VCORE', VCORE(:,CL)
             END DO
          ENDIF
          IF (ABS(ITERFLAG).LE.TOL.AND.(.NOT.ALLPARTS))
     &    THEN !If no no changes then consider all particles now
             ALLPARTS=.TRUE.
             ITERFLAG=TOL+1000
          ENDIF

          GOTO 90
      ENDIF
      END
