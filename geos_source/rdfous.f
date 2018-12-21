* This file is part of PyMieDAP, released under GNU General Public License.
* See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

      SUBROUTINE rdfous(foufile,rfou,nfou,nmat,nmugs,xmu)
      USE hdf5
Cf2py intent(in) foufile
Cf2py intent(out) nfou,nmat,nmugs,xmu
Cf2py intent(in,out) rfou
*----------------------------------------------------------------------------
*     Read Fourier coefficients from file.
*
*     Author: Ashwyn Groot, based on rdfous.f routine of Daphne M. Stam
*     Date: November 2018
*----------------------------------------------------------------------------
*----------------------------------------------------------------------------
*     Open and read the Fourier coefficients file:
*----------------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER iunf 
      PARAMETER (iunf=23)

      INTEGER nfou,nmat,nmugs

      DOUBLE PRECISION xmu(nmuMAX),rfou(nmatMAX*nmuMAX,nmuMAX,0:nfouMAX)

      CHARACTER(len=200) :: foufile
! Names (HDF5 objects)
      CHARACTER(LEN=5), PARAMETER :: groupname1 = "props"
! Sub-Group 1 name
      CHARACTER(LEN=4), PARAMETER :: groupname2 = "rfou"
! Sub-Group 1 name
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
      INTEGER(HID_T) :: dset1_id = 6
! Dataset 1 identifier
      INTEGER(HID_T) :: dset2_id = 7
! Dataset 2 identifier
      INTEGER(HID_T) :: dset3_id = 8
! Dataset 3 identifier

! Dimension array (nfou,nmat,nmugs)
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims1
      INTEGER, DIMENSION(3) :: counts
! Data buffers

! xmu array
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims2

! Rfou array
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims3

! Misc variables (e.g. loop counters)
      INTEGER :: error
! Error flag

*----------------------------------------------------------------------------
*     Open the Fourier coefficients file:
*----------------------------------------------------------------------------
      CALL h5open_f(error)
      CALL h5fopen_f (foufile, H5F_ACC_RDONLY_F, file_id, error)

      CALL h5gopen_f(file_id, groupname1, group1_id, error)
      CALL h5dopen_f(file_id, dsetname1, dset1_id, error)
      data_dims1(1)=3
      CALL h5dread_f(dset1_id, H5T_NATIVE_INTEGER, counts, data_dims1,
     .                  error)
      CALL h5dclose_f(dset1_id, error)
      nfou=counts(1)
      nmat=counts(2)
      nmugs=counts(3)
      data_dims2(1) = nmugs
      CALL h5dopen_f(file_id, dsetname2, dset2_id, error)
      CALL h5dread_f(dset2_id, H5T_NATIVE_DOUBLE, xmu(:nmugs),
     .                  data_dims2,error)
      CALL h5dclose_f(dset2_id, error)
      CALL h5gclose_f(group1_id, error)

      CALL h5gopen_f(file_id, groupname2, group2_id, error)
      data_dims3(1) = nmugs*nmat
      data_dims3(2) = nmugs
      data_dims3(3) = nfou
      CALL h5dopen_f(file_id, dsetname3, dset3_id, error)
      CALL h5dread_f(dset3_id, H5T_NATIVE_DOUBLE,
     .              rfou(:nmat*nmugs,:nmugs,0:nfou-1), data_dims3,
     .              error)
      CALL h5dclose_f(dset3_id, error)
      CALL h5gclose_f(group2_id, error)

      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error)

      GOTO 1000

*----------------------------------------------------------------------------
C996   WRITE(*,*) 'Error reading dimensions!'
C      STOP

C997   WRITE(*,*) 'Error reading mu values!'
C      STOP

C998   WRITE(*,*) 'Error reading rfou matrix!'
C      STOP

C999   WRITE(*,*) 'Error reading the Fourier file!'
C      STOP

*----------------------------------------------------------------------------
1000  RETURN
      END
