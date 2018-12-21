* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE ascii2hdf5(foufile,rfou,outputname)
      USE hdf5
Cf2py intent(in) foufile,rfou,outputname
*----------------------------------------------------------------------------
*     Read and re-write Fourier coefficients file.
*
*     Author: Ashwyn Groot
*     Date: November 2018
!     NOTE: When using this routine outside the python rewrite function,
!           the model class with information is not rewritten in the new file
*----------------------------------------------------------------------------
*----------------------------------------------------------------------------
*     Open and read and re-write the Fourier coefficients file:
*----------------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER iunf 
      PARAMETER (iunf=23)

      INTEGER m,i1,i2,i,j,ibase,nfou,nmat,nmugs,k

      DOUBLE PRECISION xmu(nmuMAX),rfou(nmatMAX*nmuMAX,nmuMAX,0:nfouMAX)
      CHARACTER ch*1
      CHARACTER(len=200) :: foufile
      CHARACTER(len=200) :: outputname
! Names (file and HDF5 objects)
      CHARACTER(LEN=5), PARAMETER :: groupname1 = "props"
! Sub-Group 1 name
      CHARACTER(LEN=4), PARAMETER :: groupname2 = "rfou"
! Sub-Group 2 name

! Dataset 1 name
      CHARACTER(LEN=12), PARAMETER :: dsetname1 = "Array-counts"
! Dataset 2 name
      CHARACTER(LEN=9), PARAMETER :: dsetname2 = "XMU-array"
! Dataset 3 name
      CHARACTER(LEN=10), PARAMETER :: dsetname3 = "RFOU-array"

! Identifiers
      INTEGER(HID_T) :: file_id = 0
! File identifier
      INTEGER(HID_T) :: group1_id = 1
! Group 1 identifier
      INTEGER(HID_T) :: group2_id = 2
! Group 2 identifier
      INTEGER(HID_T) :: dspace1_id = 4
! Dataspace 1 identifier
      INTEGER(HID_T) :: dspace2_id = 5
! Dataspace 2 identifier
      INTEGER(HID_T) :: dspace3_id = 6
! Dataspace 3 identifier
      INTEGER(HID_T) :: dset1_id = 8
! Dataset 1 identifier
      INTEGER(HID_T) :: dset2_id = 9
! Dataset 2 identifier
      INTEGER(HID_T) :: dset3_id = 10
! Dataset 3 identifier


! Dimension array (nfou,nmat,nmugs)
      INTEGER :: rank
! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(1) :: dims1 = (/3/)
! Dataset dimensions
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims1
      INTEGER, DIMENSION(3) :: dset_data1
! Data buffers

! xmu array
      INTEGER(HSIZE_T), DIMENSION(1) :: dims2
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims2
      REAL*8, DIMENSION(:), ALLOCATABLE :: dset_data2

! Rfou array
      INTEGER(HSIZE_T), DIMENSION(3) :: dims3
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims3
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: dset_data3

! Misc variables (e.g. loop counters)
      INTEGER :: error
! Error flag
      rfou=0.D0
      OPEN(unit=iunf,file=foufile,status='old',err=999)

*----------------------------------------------------------------------------
*     Read the header (the dap.in file contents):
*----------------------------------------------------------------------------
2     READ(iunf,'(a1)',err=997,end=998) ch
      IF (ch.EQ.'#') GOTO 2
      BACKSPACE(iunf)

*----------------------------------------------------------------------------
*     Read the accuracy, the matrix size, the number of Gauss points and the
*     Gaussian integration points:
*----------------------------------------------------------------------------
      READ(iunf,*) nmat
      READ(iunf,*) nmugs
      DO i=1,nmugs
         READ(iunf,*) xmu(i)
      ENDDO

      nfou=0
20    DO i=1,nmugs
         ibase= (i-1)*nmat
         DO j=1,nmugs
            READ(iunf,*,END=21) m,i1,i2,
     .                        (rfou(ibase+k,j,nfou),k=1,nmat)
         ENDDO
      ENDDO
      nfou= nfou+1
      GOTO 20
21    CLOSE(iunf)

! Size of nfou's:
      nfou= nfou+1
*----------------------------------------------------------------------
*     Re-write the Fourier-coefficients to hdf5 file:
*----------------------------------------------------------------------
! Initialize Fortran interface
      CALL h5open_f(error)
! Create a new file
C      CALL h5fcreate_f(outputname, H5F_ACC_TRUNC_F, file_id, error)
      CALL h5fopen_f (outputname, H5F_ACC_RDWR_F, file_id, error)

*-----------------------------------------------------------------------
* Put nfou,nmat,nmugs,xmu in group 1:
*-----------------------------------------------------------------------
! Create a group in the HDF5 file: properties group
      CALL h5gcreate_f(file_id, groupname1, group1_id, error)

*********************************************
! Create dataspace: counts (the dataset is next)
*********************************************
      data_dims1(1) = 3
      rank = 1
      dset_data1(:) = (/nfou,nmat,nmugs/)
! Create dataspace 1
      CALL h5screate_simple_f(rank, dims1, dspace1_id, error)
! Create dataset 1 with default properties
      CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_INTEGER,
     .            dspace1_id,dset1_id, error)
! Write dataset 1
      CALL h5dwrite_f(dset1_id, H5T_NATIVE_INTEGER, dset_data1,
     .            data_dims1,error)
! Close access to dataset 1
      CALL h5dclose_f(dset1_id, error)
! Close access to data space 1
      CALL h5sclose_f(dspace1_id, error)

*********************************************
! Create dataspace: xmu (the dataset is next)
*********************************************
      data_dims2(1) = nmugs
      dims2(1)=nmugs
      rank = 1
      ALLOCATE(dset_data2(nmugs))
      dset_data2 = xmu(:nmugs)
! Create dataspace 2
      CALL h5screate_simple_f(rank, dims2, dspace2_id, error)
! Create dataset 2 with default properties
      CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE,
     .            dspace2_id,dset2_id, error)
! Write dataset 2
      CALL h5dwrite_f(dset2_id, H5T_NATIVE_DOUBLE, dset_data2,
     .            data_dims2,error)
! Close access to dataset 2
      CALL h5dclose_f(dset2_id, error)
! Close access to data space 2
      CALL h5sclose_f(dspace2_id, error)
      DEALLOCATE(dset_data2)
*********************************************

! Close the group
      CALL h5gclose_f(group1_id, error)
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
* Put rfou in group 2:
*-----------------------------------------------------------------------
! Create a group in the HDF5 file: rfou group
      CALL h5gcreate_f(file_id, groupname2, group2_id, error)

*********************************************
! Create dataspace: rfou (the dataset is next)
*********************************************
      data_dims3(1) = nmugs*nmat
      data_dims3(2) = nmugs
      data_dims3(3) = nfou
      dims3(1)=nmugs*nmat
      dims3(2)=nmugs
      dims3(3)=nfou
      rank = 3
      ALLOCATE(dset_data3(nmugs*nmat,nmugs,0:nfou-1))
      dset_data3 = rfou(:nmat*nmugs,:nmugs,0:nfou-1)
! Create dataspace 1
      CALL h5screate_simple_f(rank, dims3, dspace3_id, error)
! Create dataset 1 with default properties
      CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE,
     .            dspace3_id,dset3_id, error)
! Write dataset 1
      CALL h5dwrite_f(dset3_id, H5T_NATIVE_DOUBLE, dset_data3,
     .            data_dims3,error)
! Close access to dataset 1
      CALL h5dclose_f(dset3_id, error)
! Close access to data space 1
      CALL h5sclose_f(dspace3_id, error)
      DEALLOCATE(dset_data3)
*********************************************

! Close the group
      CALL h5gclose_f(group2_id, error)
*-----------------------------------------------------------------------
! Close the file
      CALL h5fclose_f(file_id, error)
! Close FORTRAN interface
      CALL h5close_f(error)

      GOTO 1000


*----------------------------------------------------------------------------
999   WRITE(*,*) 'Error opening Fourier file!'
      STOP
998   WRITE(*,*) 'Error: unexpected end of the Fourier file!'
      STOP
997   WRITE(*,*) 'Error reading the Fourier file!'
      STOP

*----------------------------------------------------------------------------
1000  RETURN
      END
