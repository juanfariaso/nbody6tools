      ! JP: This whole module I obtained it at the beggining of my research life.
      ! Not sure who wrote it but I have some suspects.
      ! TODO: Find out and put the proper credits here.
      SUBROUTINE qparameter(x,y,z,nstr,qpar,sigma,average,zeroaxis,rmax)
      IMPLICIT none

      !EXTERNAL mst
      integer, intent(in) :: nstr
      real, intent(in) :: x(nstr),y(nstr),z(nstr),rmax !maximun size of the cluster (for normalization)
      logical, intent(in) :: average
      integer,intent(in) :: zeroaxis !if average=False zeroaxis is the axis
                                      !that will be taked as zero

      real,intent(out) :: qpar, sigma !if average true, qpar is an average of
                                      ! the projections xy,xz,yz and sigma the
                                      ! stdv

      !internals
      real cut,mbar,sbar,mtot,stot,pi,Q(1:3),sbar2
      real xdumy(nstr),ydumy(nstr),zdumy(nstr),length(1:nstr-1)
      integer node(1:nstr-1,1:2),treemem(1:nstr),ncount
      integer i,j,cnt,d,ax,za


      pi=atan(1.0)*4.0
      cut=1.E20 !no cut implemented
      za=zeroaxis
     
 
      if (average) then
           d=3
           za=1
      else
           d=1
      end if 
      Qpar=0.0
      do ax=1,d
          do i=1,nstr
             if (za==1) then
                xdumy(i)=0.0
                ydumy(i)=y(i)
                zdumy(i)=z(i)
             else if (za==2) then
                xdumy(i)=x(i)
                ydumy(i)=0
                zdumy(i)=z(i)
             else if (za==3) then
                xdumy(i)=x(i)
                ydumy(i)=y(i)
                zdumy(i)=0
             else 
                 PRINT*, "Wrong zeroaxis"
                 STOP
             end if
          enddo

          
          call mst(nstr,xdumy,ydumy,zdumy,cut,node,length,treemem,ncount)
          mbar=0
          do i=1,nstr-1
             mbar=mbar+length(i)
          end do
          mbar=mbar
          !mbar is normalized
          mbar=mbar/SQRT(nstr*pi*rmax**2)
          sbar=0
          DO i=1,(nstr-1)
             DO j=i,nstr
                   sbar=sbar+ sqrt( (xdumy(i)-xdumy(j))**2 + (ydumy(i)-ydumy(j))**2 +  (zdumy(i)-zdumy(j))**2)
             ENDDO
          ENDDO
          sbar= 2.0*sbar/(nstr**2 - nstr)/rmax


          Qpar=mbar/sbar
          if (average) then
             Q(ax)=Qpar
          end if
          za=za+1
      end do
          sigma=0

      if (average) then
           Qpar=(Q(1)+Q(2)+Q(3))/3.0
           sigma=SQRT( ((Q(1)-Qpar)**2+(Q(2)-Qpar)**2+(Q(3)-Qpar)**2 )/3.0 )
      end if
      end subroutine
! ======================================================================
! ======================================================================
!
   SUBROUTINE mst(n,x,y,z,cut,node,length,treemem,ncount)
! generates a minimum spanning tree 
! if 2d one of the x,y,z arrays must be a zero array
   IMPLICIT NONE
! inputs
   INTEGER, INTENT(in) :: n
   REAL, INTENT(in) :: x(1:n),y(1:n),z(1:n),cut
! outputs
   INTEGER, INTENT(out) :: node(1:n-1,1:2) ! node connections
   REAL, INTENT(out) :: length(1:n-1)      ! connection lengths
   INTEGER, INTENT(out) :: treemem(1:n),ncount
! internals
   INTEGER, DIMENSION(:), ALLOCATABLE :: list
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: idents
!   INTEGER ::i,j,memnum,told,ncount,imin,jmin,nlist
   INTEGER ::i,j,memnum,told,imin,jmin,nlist

   REAL, DIMENSION(:), ALLOCATABLE :: sep
!
! how long is the list?
   nlist=0
   DO i=1,n-1
     DO j=i+1,n
       nlist=nlist + 1
     END DO
   END DO
!
   ALLOCATE(sep(1:nlist))
   ALLOCATE(list(1:nlist))
   ALLOCATE(idents(1:2,1:nlist))
!
! get distances between *all* points and put in a list
! I'm sure this could be clever and work-out what i and j are 
! from the position in the array, but I can't be bothered
   nlist=0
   DO i=1,n-1
     DO j=i+1,n
       nlist=nlist + 1
       sep(nlist)=SQRT((x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2)
       list(nlist)=nlist
       idents(1,nlist)=i
       idents(2,nlist)=j
     END DO
   END DO
!
! heapsort separations
   CALL heapsort(nlist,sep,list)
!
! loop to allocate nodes
   node=0
   length=0.
   ncount=1
   treemem=0
   memnum=1
   nlist=0
   DO
! smallest distance
     nlist=nlist + 1
     length(ncount)=sep(list(nlist))
     imin=idents(1,list(nlist))
     jmin=idents(2,list(nlist))
!
     node(ncount,1)=imin
     node(ncount,2)=jmin
!
     IF (treemem(imin)/=0 .AND. treemem(imin)==treemem(jmin)) GOTO 23
!
! check if neither ptcl is already in a tree -------------------
     IF (treemem(imin)==0 .AND. treemem(jmin)==0) THEN
! not already part of a tree, so give a new treenum
       treemem(imin)=memnum
       treemem(jmin)=memnum
! advance memnum and ncount
       memnum=memnum + 1
       ncount=ncount + 1
       GOTO 23
     END IF
!
! check if one or the other is ---------------------------------
     IF (treemem(imin)==0 .AND. treemem(jmin)/=0) THEN
       treemem(imin)=treemem(jmin)
       ncount=ncount + 1
       GOTO 23
     END IF
!
     IF (treemem(imin)/=0 .AND. treemem(jmin)==0) THEN
       treemem(jmin)=treemem(imin)
       ncount=ncount + 1
       GOTO 23
     END IF
!
! are they in different trees? ---------------------------------
     IF (treemem(imin)/=0 .AND. treemem(jmin)/=treemem(imin)) THEN
! they'll all be part of this tree now
       told=treemem(jmin)
! loop over all current nodes
       DO i=1,n
         IF (treemem(i)==told) treemem(i)=treemem(imin)
       END DO
       ncount=ncount + 1
       GOTO 23
     END IF
!
23  IF (ncount==n.OR.sep(list(nlist+1))>cut) EXIT
!
   END DO
!
   RETURN
!
   END SUBROUTINE mst


   SUBROUTINE heapsort(psort,measureof,pwhichhas)
! does a heapsort (by APW)
   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: psort              ! number of values to be sorted.
   REAL,    INTENT(IN)  :: measureof(1:psort) ! values to be sorted.
   INTEGER, INTENT(OUT) :: pwhichhas(1:psort) ! identifier of value.
   INTEGER              :: rank               ! rank of value.
   INTEGER              :: ranknow            ! dummy rank.
   INTEGER              :: ranktest           ! dummy rank.
!
   DO rank=2,psort                ! THIS DO-LOOP BUILDS THE BINARY HEAP
     ranknow=rank
1    IF (ranknow==1) CYCLE
     ranktest=ranknow/2
     IF (measureof(pwhichhas(ranktest))>=measureof(pwhichhas(ranknow))) CYCLE
     CALL swapi(pwhichhas(ranknow),pwhichhas(ranktest))
     ranknow=ranktest
     GOTO 1
   END DO
!
   DO rank=psort,2,-1             ! AND THIS DO-LOOP INVERTS THE BINARY HEAP
     CALL swapi(pwhichhas(rank),pwhichhas(1))
     ranknow=1
2    ranktest=2*ranknow
     IF (ranktest>=rank) CYCLE
     IF ((measureof(pwhichhas(ranktest+1))>measureof(pwhichhas(ranktest))) &
                            & .AND.(ranktest+1<rank)) ranktest=ranktest+1
     IF (measureof(pwhichhas(ranktest))<=measureof(pwhichhas(ranknow))) CYCLE
     CALL swapi(pwhichhas(ranknow),pwhichhas(ranktest))
     ranknow=ranktest
     GOTO 2
   END DO
!
   RETURN
   END SUBROUTINE heapsort
!
! ===========================================================================
! ===========================================================================
!
   SUBROUTINE swapi(item1,item2)
!
   IMPLICIT NONE
   INTEGER :: item0,item1,item2
!
   item0=item1; item1=item2; item2=item0
!
   RETURN
   END SUBROUTINE swapi
