      SUBROUTINE writsc(outfile,idis,nsubr,ngaur,coefs,ncoefs,csth,miec,
     .                  wav,nr,ni,rmin,rmax,par1,par2,par3)

*-----------------------------------------------------------------------
*  On entry :                                                        
*      outfile     name of the coefficient file             
*      coefs    array containing the coefficients, according to : 
*                                                               
*                         ( alpha1(l) beta2(l)     0         0      ) 
*      coefs(i,j,l)   =  ( beta2(l)  alpha2(l)    0         0      )
*                         (   0         0       alpha3(l) beta2(l)  )
*                         (   0         0      -beta2(l)  alpha4(l) )
*                                                                 
*      ncoefs      maximum index l for which coefficients are non-zero 
*      a         albeDO for single scattering                   
*      csth      asymmetry parameter <cos th>                
*      Qext      extinction efficiency                     
*      Qsca      scattering efficiency                    
*      G         average geometrical cross section       
*      V         average volume                       
*      s1,...,s6 strings with description of scattering matrix   
*  On exit :                                                
*      The coefficients and other inFORMATion is written on the file 
*      with unit number iunit in the standard FORMAT.               
*                                                 VLD January 9, 1989 
*-----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER l,iunit,ncoefs,idis,nsubr,ngaur

      DOUBLE PRECISION a,csth,Qext,Qsca,G,V,wav,ni,nr,
     .                 rmin,rmax,par1,par2,par3

      DOUBLE PRECISION coefs(nmatMAX,nmatMAX,0:ncoefsMAX),
     .                 miec(13)

      DOUBLE PRECISION eps
      PARAMETER (eps=1.D-6)

      CHARACTER*60 s1,s2,s3,s4,s5,s6
      CHARACTER outfile*16

Cf2py intent(in) outfile,idis,nsubr,ngaur,coefs,ncoefs,csth,miec
Cf2py intent(in) wav,nr,ni,rmin,rmax,par1,par2,par3

*-----------------------------------------------------------------------
*     Put miec parameters in own variables:
*-----------------------------------------------------------------------
      Qsca= miec(3)
      Qext= miec(4)
      a= miec(5)
      G= miec(6)
      V= miec(10)

*-----------------------------------------------------------------------
*     Check parameters:
*-----------------------------------------------------------------------
      IF (ncoefs.GT.ncoefsMAX) THEN
         WRITE(*,*) 'writsc: ncoefs= ',ncoefs,
     .             ' maximum ncoefs = ',ncoefsMAX
         STOP
      ENDIF
      IF (DABS(csth-coefs(1,1,1)/3.D0).GT.eps) THEN
         WRITE(*,*) 'writsc: warning:    asymmetry parameter = ',csth
         WRITE(*,*) '         1/3 first order Legendre coeff. = '
     .                                   ,coefs(1,1,1)/3.D0
      ENDIF
      IF ((a.LT.0.D0).OR.(a.GT.1.D0)) THEN
         WRITE(*,*) 'writsc: warning: single scattering albeDO = ',a
     .                                       ,' outside interval [0;1]'
      ENDIF
      IF (Qext.LT.0.D0) THEN
         WRITE(*,*) 'writsc: warning: extinction efficiency = ',Qext
     .                                                   ,' negative !'
      ENDIF
      IF (Qsca.LT.0.D0) THEN
         WRITE(*,*) 'writsc: warning: scattering efficiency = ',Qsca
     .                                                   ,' negative !'
      ENDIF
      IF (Qsca.GT.Qext) THEN
         WRITE(*,*) 'writsc: warning:    scattering efficiency = ',Qsca
         WRITE(*,*) '         larger than extinction efficiency = ',Qext
      ENDIF
      IF (G.LT.0.D0) THEN
         WRITE(*,*) 'writsc: warning:    geometrical cross section =',G
     .                                                   ,' negative !'
      ENDIF
      IF (V.LT.0.D0) THEN
         WRITE(*,*) 'writsc: warning:    average volume = ',V
     .                                                   ,' negative !'
      ENDIF
      IF (DABS(coefs(1,1,0)-1.D0).GT.eps) THEN
         WRITE(*,*) 'writsc: warning: zero Legendre coeff = ',
     .               coefs(1,1,0),'         not equal to one !'
      ENDIF

*-----------------------------------------------------------------------
*     Make the headlines:
*-----------------------------------------------------------------------
      s1='PRODUCED BY THE MEERHOFF MIE PROGRAM VERSION 3.1'
      WRITE(s2,1001) wav,nr,ni

      IF (idis.NE.0) THEN
          WRITE(s5,1007) rmin,rmax
          WRITE(s6,1009) nsubr,ngaur
      ENDIF
      IF (idis.EQ.0) THEN
         s5='                                                          '
         s6='                                                          '
         s3='No size distribution (single particle)                    '
         WRITE(s4,1003) par1
      ELSEIF (idis.EQ.1) THEN
         s3='Size distribution : two parameter gamma                   '
         WRITE(s4,1011) par1, par2
      ELSE IF (idis.EQ.2) THEN
         s3='Size distribution : two parameter gamma                   '
         WRITE(s4,1012) par1, par2
      ELSEIF (idis.EQ.3) THEN
         s3='Size distribution : bimodal gamma with equal mode weights '
         WRITE(s4,1013) par1, par2, par3
      ELSEIF (idis.EQ.4) THEN
         s3='Size distribution : log normal                            '
         WRITE(s4,1014) par1, par2
      ELSEIF (idis.EQ.5) THEN
         s3='Size distribution : log normal                            '
         WRITE(s4,1015) par1, par2
      ELSEIF (idis.EQ.6) THEN
         s3='Size distribution : power law                             '
         WRITE(s4,1016) par1, par2, par3
      ELSEIF (idis.EQ.7) THEN
         s3='Size distribution : modified gamma                        '
         WRITE(s4,1017) par1, par2, par3
      ELSEIF (idis.EQ.8) THEN
         s3='Size distribution : modified gamma                        '
         WRITE(s4,1018) par1, par2, par3
      ELSE
         WRITE(*,*) ' MAIN: illegal size distribution nr:',idis
      ENDIF

*-----------------------------------------------------------------------
*     Write some parameters to file:
*-----------------------------------------------------------------------
      OPEN (file=outfile,unit=iunit)
      WRITE(iunit,1)
      WRITE(iunit,2) s1
      WRITE(iunit,2) s2
      WRITE(iunit,2) s3
      WRITE(iunit,2) s4
      WRITE(iunit,2) s5
      WRITE(iunit,2) s6
      WRITE(iunit,3) a,csth,Qext,Qsca,G,V

*-----------------------------------------------------------------------
*     Write the expansion coefficients to file:
*-----------------------------------------------------------------------
      DO l=0,ncoefs
         WRITE(iunit,4) l,coefs(1,1,l),coefs(2,2,l),coefs(3,3,l)
     .                   ,coefs(4,4,l),coefs(1,2,l),coefs(3,4,l)
         IF ((DABS(coefs(1,3,l)).gt.eps).or.
     .       (DABS(coefs(1,4,l)).gt.eps).or.
     .       (DABS(coefs(2,3,l)).gt.eps).or.
     .       (DABS(coefs(2,4,l)).gt.eps).or.
     .       (DABS(coefs(3,1,l)).gt.eps).or.
     .       (DABS(coefs(3,2,l)).gt.eps).or.
     .       (DABS(coefs(4,1,l)).gt.eps).or.
     .       (DABS(coefs(4,2,l)).gt.eps))
     .       WRITE(*,*) 'writsc: warning: nonzero coeff. in upper'
     .                 ,'        right or lower left corner l = ',l
      ENDDO

      WRITE(iunit,4) -1,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0
      CLOSE(iunit)

*-----------------------------------------------------------------------
*     Formats:
*-----------------------------------------------------------------------
   1  FORMAT(' EXPANSION COEFFICIENTS SCATTERING MATRIX')
   2  FORMAT(' ',a60)
   3  FORMAT(' single scattering albedo  a  =',e22.14,/
     .      ,' asymmetry parameter <cos th> =',e22.14,/
     .      ,' extinction efficiency  Qext  =',e22.14,/
     .      ,' scattering efficiency  Qsca  =',e22.14,/
     .      ,' geometrical cross section G  =',e22.14,/
     .      ,' average volume            V  =',e22.14,/
     .      ,'   l'
     .      ,'       alpha1             alpha2      '
     .      ,'       alpha3             alpha4      '
     .      ,'       beta1              beta2       ')
   4  FORMAT(i4,6f19.14)

1001  FORMAT('lambda= ',f11.7,' Re(m) = ',f11.7,' Im(m) = ',f11.7)
1003  FORMAT('Radius = ',1f11.7
     .                       ,'                                       ')
1007  FORMAT('rmin  = ',f11.7,' rmax  = ',f11.7,'                     ')
1009  FORMAT('r integration with ',i3,' interval(s) of '
     .                                            ,i3,' Gauss points  ')
1011  FORMAT('alpha = ',f11.7,' b    = ',f11.7,'                    ')
1012  FORMAT('reff  = ',f11.7,' veff  = ',f11.7,'                   ')
1013  FORMAT('reff1 = ',f11.7,' reff2 = ',f11.7,' veff = ',f11.7)
1014  FORMAT('rg     = ',f11.7,' sigma = ',f11.7,'                    ')
1015  FORMAT('reff   = ',f11.7,' veff  = ',f11.7,'                    ')
1016  FORMAT('alpha = ',f11.7,' rmin  = ',f11.7,' rmax  = ',f11.7)
1017  FORMAT('alpha = ',f11.7,' rc    = ',f11.7,' gamma = ',f11.7)
1018  FORMAT('alpha = ',f11.7,' b     = ',f11.7,' gamma = ',f11.7)
*-----------------------------------------------------------------------
      RETURN
      END

