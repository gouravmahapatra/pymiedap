* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE hdf5deriv2ascii(foufile,rfou,outputname)
      USE hdf5
Cf2py intent(in) foufile,rfou,outputname
*----------------------------------------------------------------------------
*     Read and re-write Fourier coefficients file.
*
*     Author: Ashwyn Groot
*     Date: November 2018
*----------------------------------------------------------------------------
*----------------------------------------------------------------------------
*     Open and read and re-write the Fourier coefficients file:
*----------------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER iunf,iunfou
      PARAMETER (iunf=23)
      DOUBLE PRECISION eps,t1,t2,t3,t4
      PARAMETER (eps=1.D-6)

      INTEGER m,i,j,ibase,nfou,nmat,nmugs

      DOUBLE PRECISION xmu(nmuMAX),rfou(nmatMAX*nmuMAX,nmuMAX,0:nfouMAX)
      DOUBLE PRECISION wmu(nmuMAX)
      REAL*8, DIMENSION(:), ALLOCATABLE :: xmu_out
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: rfou_out
      CHARACTER*200 fst
      CHARACTER(len=200) :: foufile
      CHARACTER(len=200) :: outputname
! Names (file and HDF5 objects)
      CHARACTER(LEN=5), PARAMETER :: groupname1 = "props" ! Sub-Group 1 name
      CHARACTER(LEN=4), PARAMETER :: groupname2 = "rfou" ! Sub-Group 1 name
! Dataset 1 name
      CHARACTER(LEN=12), PARAMETER :: dsetname1 = "Array-counts"
! Dataset 2 name
      CHARACTER(LEN=9), PARAMETER :: dsetname2 = "XMU-array"
! Dataset 3 name
      CHARACTER(LEN=10), PARAMETER :: dsetname3 = "RFOU-array"

! Identifiers
      INTEGER(HID_T) :: file_id = 0      ! File identifier
      INTEGER(HID_T) :: group1_id = 1    ! Group 1 identifier
      INTEGER(HID_T) :: group2_id = 2    ! Group 2 identifier
      INTEGER(HID_T) :: dset1_id = 6     ! Dataset 1 identifier
      INTEGER(HID_T) :: dset2_id = 7     ! Dataset 2 identifier
      INTEGER(HID_T) :: dset3_id = 8     ! Dataset 3 identifier


! Dimension array (nfou,nmat,nmugs)
      INTEGER :: rank                 ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims1
      INTEGER, DIMENSION(3) :: counts  ! Data buffers

! xmu array
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims2

! Rfou array
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims3

! Misc variables (e.g. loop counters)
      INTEGER :: error ! Error flag
      rfou=0.D0

      CALL h5open_f(error)
      CALL h5fopen_f (foufile, H5F_ACC_RDONLY_F, file_id, error)

      CALL h5gopen_f(file_id, groupname1, group1_id, error)
      rank=1
C      CALL h5sopen_f(rank, dims1, dspace1_id, error)
      CALL h5dopen_f(file_id, dsetname1, dset1_id, error)
      CALL h5dread_f(dset1_id, H5T_NATIVE_INTEGER, counts, data_dims1,
     .                  error)
      CALL h5dclose_f(dset1_id, error)
C      CALL h5sclose_f(dspace1_id, error)
      nfou=counts(1)
      nmat=counts(2)
      nmugs=counts(3)
      ALLOCATE(xmu_out(nmugs),rfou_out(nmugs*nmat,nmugs,nfou+1))
C      CALL h5sopen_f(rank, dims2, dspace2_id, error)
      CALL h5dopen_f(file_id, dsetname2, dset2_id, error)
      CALL h5dread_f(dset2_id, H5T_NATIVE_DOUBLE, xmu_out, data_dims2,
     .                  error)
      CALL h5dclose_f(dset2_id, error)
C      CALL h5sclose_f(dspace2_id, error)
      CALL h5gclose_f(group1_id, error)

      xmu(:nmugs)=xmu_out

      CALL h5gopen_f(file_id, groupname2, group2_id, error)
      rank=3
C      CALL h5sopen_f(rank, dims3, dspace3_id, error)
      CALL h5dopen_f(file_id, dsetname3, dset3_id, error)
      CALL h5dread_f(dset3_id, H5T_NATIVE_DOUBLE, rfou_out, data_dims3,
     .                  error)
      CALL h5dclose_f(dset3_id, error)
C      CALL h5sclose_f(dspace3_id, error)
      CALL h5gclose_f(group2_id, error)

      rfou(:nmat*nmugs,:nmugs,0:nfou)=rfou_out

      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error)

*----------------------------------------------------------------------
*     Re-write the Fourier-coefficients to ASCII file:
*----------------------------------------------------------------------
      iunfou= 123
*       Value set arbitrarily

      OPEN(unit=iunfou,file=outputname)
      WRITE(iunfou,330)
      WRITE(iunfou,331)
      WRITE(iunfou,332)
      WRITE(iunfou,333)

      WRITE(iunfou,'(I3)') nmat
      WRITE(iunfou,'(I3)') nmugs
      CALL gauleg(nmuMAX,nmugs-1,0.D0,1.D0,xmu,wmu)
      wmu(nmugs)= 0.5D0
      DO i=1,nmugs
         WRITE(iunfou,'(E16.8,3X,E16.8)') xmu(i),wmu(i)
      ENDDO

      DO m=0,nfou
         fst = '(I4,2X,2(I3,2X),'
         DO i=1,nmugs
            ibase= (i-1)*nmat
            DO j=1,nmugs
               fst = '(I4,2X,2(I3,2X),'
               IF (nmat.EQ.1) THEN
                  t1= rfou(ibase+1,j,m)
                  IF (DABS(t1).LT.eps) THEN
                     t1= 0.D0
                     fst = trim(fst) // 'F2.0)'
                  ELSE
                     fst = trim(fst) // 'E16.8)'
                  ENDIF

                  WRITE(iunfou,fmt=fst) m,i,j,t1

               ELSEIF (nmat.EQ.3) THEN
                  t1= rfou(ibase+1,j,m)
                  t2= rfou(ibase+2,j,m)
                  t3= rfou(ibase+3,j,m)
                  IF (DABS(t1).LT.eps) THEN
                     t1= 0.D0
                     fst = trim(fst) // 'F2.0,1X,'
                  ELSE
                     fst = trim(fst) // 'E16.8,1X,'
                  ENDIF

                  IF (DABS(t2).LT.eps) THEN
                     t2= 0.D0
                     fst = trim(fst) // 'F2.0,1X,'
                  ELSE
                     fst = trim(fst) // 'E16.8,1X,'
                  ENDIF

                  IF (DABS(t3).LT.eps) THEN
                     t3= 0.D0
                     fst = trim(fst) // 'F2.0)'
                  ELSE
                     fst = trim(fst) // 'E16.8)'
                  ENDIF

                  WRITE(iunfou,fmt=fst) m,i,j,t1,t2,t3

               ELSE
                  t1= rfou(ibase+1,j,m)
                  t2= rfou(ibase+2,j,m)
                  t3= rfou(ibase+3,j,m)
                  t4= rfou(ibase+4,j,m)
                  IF (DABS(t1).LT.eps) THEN
                     t1= 0.D0
                     fst = trim(fst) // 'F2.0,1X,'
                  ELSE
                     fst = trim(fst) // 'E16.8,1X,'
                  ENDIF

                  IF (DABS(t2).LT.eps) THEN
                     t2= 0.D0
                     fst = trim(fst) // 'F2.0,1X,'
                  ELSE
                     fst = trim(fst) // 'E16.8,1X,'
                  ENDIF

                  IF (DABS(t3).LT.eps) THEN
                     t3= 0.D0
                     fst = trim(fst) // 'F2.0,1X,'
                  ELSE
                     fst = trim(fst) // 'E16.8,1X,'
                  ENDIF

                  IF (DABS(t4).LT.eps) THEN
                     t4= 0.D0
                     fst = trim(fst) // 'F2.0)'
                  ELSE
                     fst = trim(fst) // 'E16.8)'
                  ENDIF

                  WRITE(iunfou,fmt=fst) m,i,j,t1,t2,t3,t4

               ENDIF

            ENDDO
         ENDDO

      ENDDO

      CLOSE(iunfou)
      GOTO 1000
*-----------------------------------------------------------------------
*     Formats:
*-----------------------------------------------------------------------
330   FORMAT('# For the correct use of this datafile, see Stam and',
     .       ' de Haan 2013')
331   FORMAT('# Please refer to Stam and de Haan 2013 when using',
     .       ' this data')
332   FORMAT('# The wavelength for the calculations:',F12.6)
333   FORMAT('#')
*-----------------------------------------------------------------------
1000  RETURN
      END
