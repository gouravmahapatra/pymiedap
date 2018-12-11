* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE newfou(nmat,xmu,rfou,outputname,lg,nfou,nmu)
      USE HDF5
Cf2py intent(in) nmat,xmu,rfou,outputname,lg,nfou,nmu
      IMPLICIT NONE

      INCLUDE 'max_incl'
      INTEGER nfou,nmat,nmu,lg
      DOUBLE PRECISION xmu(nmu),rfou(nmu*nmat,nmu,0:nfou)
      CHARACTER(len=200) :: outputname

*----------------------------------------------------------------------
*     Write the Fourier-coefficients to file:
*----------------------------------------------------------------------
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
      INTEGER(HID_T) :: dspace1_id = 3   ! Dataspace 1 identifier
      INTEGER(HID_T) :: dspace2_id = 4   ! Dataspace 2 identifier
      INTEGER(HID_T) :: dspace3_id = 5   ! Dataspace 2 identifier
      INTEGER(HID_T) :: dset1_id = 6     ! Dataset 1 identifier
      INTEGER(HID_T) :: dset2_id = 7     ! Dataset 2 identifier
      INTEGER(HID_T) :: dset3_id = 8     ! Dataset 3 identifier


! Dimension array (nfou,nmat,nmugs)
      INTEGER :: rank                 ! Dataset rank
      INTEGER(HSIZE_T), DIMENSION(1) :: dims1 = (/3/) ! Dataset dimensions
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims1
      INTEGER, DIMENSION(3) :: dset_data1   ! Data buffers

! xmu array
      INTEGER(HSIZE_T), DIMENSION(1) :: dims2
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims2
      REAL*8, DIMENSION(nmu) :: dset_data2

! Rfou array
      INTEGER(HSIZE_T), DIMENSION(3) :: dims3
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims3
      REAL*8, DIMENSION(nmu*nmat,nmu,0:nfou) :: dset_data3


! Misc variables (e.g. loop counters)
      INTEGER :: error ! Error flag
! =====================================================================

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
      dset_data1(:) = (/nfou+1,nmat,nmu/)
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
      data_dims2(1) = nmu
      dims2(1)=nmu
      rank = 1
      dset_data2 = xmu
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
      data_dims3(1) = nmu*nmat
      data_dims3(2) = nmu
      data_dims3(3) = nfou+1
      dims3(1)=nmu*nmat
      dims3(2)=nmu
      dims3(3)=nfou+1
      rank = 3
      dset_data3 = rfou
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
*********************************************

! Close the group
      CALL h5gclose_f(group2_id, error)
*-----------------------------------------------------------------------
! Close the file
      CALL h5fclose_f(file_id, error)
! Close FORTRAN interface
      CALL h5close_f(error)

*-----------------------------------------------------------------------
      RETURN
      END
