c   diese version 'd20'  hat reinen h2o (CH4)  - trend
c   temperatur von 1878 (1. sol min)
c   ijahr_fixtemp = 1878  siehe eingabe
c   jahr 1976 : lesen der 1h bin files
c   ijahr_fixdyn = 1976   siehe eingabe



c  diese version basiert auf 'l5_eiszp_lima_lc25_d16.f'
c  aber hier ... siehe c c, filenamen mit endung:  'd18'
c
c  hier kein CH4 Anstieg beruecksichtigt (wie reiner transport in, Ly-alpha an, files 'hs02ch4_n'
c                                         l5_h2otrans_zpress_lc25_2trend_ch4.f)
c           aber alle Jahre CH4 wert von 1871 (5 Jahre Shift)
c  hier NEU
c         if (float(i_jahr)-5.eq.a1)  ch4_akt = a2    ! shift 5 jahre: age of air
cc         if (1881.-5.eq.a1)  ch4_akt = a2    ! shift 5 jahre: age of air
c ende NEU
c--------------------------------------------------------------------------
c CH4: hier wird h2o erniedrigt WICHTIG !!!
c ch4 werte in tabelle : CH4_1861_2008
c   ...  in sub_init_global    ->  scale_ch4 (in prozent)



c
c
c  neu zusaetzlich staub auf mesopause 20%  : wie d13, jetzt aber d14
c     + j_zbound_ob um 800 m hoeher
c     d14:  staeube werden auf das Druckniveau z_p = 89.1 km gelegt,
c           dies entspricht 86 km geohoehe fuer 2008 1-20 Juli wie SOFIE
c
c  neu fuer era interim  beta wr zu gross
c  1) wsedi  erhoehter Wert der Formel 
c  2) staub rinit   xhisto 0.2 kleiner
c  3) partikel temp: worka_tab  hoeherer wert
c  4)  volufak  (stellvertreter)  60 -> 80   und 0.6 (0.7)
c  5)  hminit3d um 10 % veringert
c  6)  volufak nochmal um 10 % veringert
c  7)  staub rinit nochmal um 0.2 verkleinert  
c
c
c   diese version hat stellvertreter= mal 1/2  (siehe volufak)
c
c
c    neu:    hergestellt aus version erainterim_c7
c              1)  Staub fix 85.5 - 86.5
c                                         in SUBROUTINE sub_init_staub
c                                         + 2 mal in hauptschleife

c              2) partikel temp von 200% -> 100 %

c              3)  volufak Stellvertreter       volufak/100.
c                                     Zeile 2761:   change_n(n) =  backgr*volufak*0.7*8./20.   ! neu 
c                                     Zeile 2804:   change_n(n) =  backgr*volufak*0.7*8./20.   ! neu 

c              4)  Kzz: turbkzz  in SUBROUTINE sub_init_global
c                  Zeile 2411     <= 83 km

c              5)   Zeile 2849   wrando = 0.5*(....)

c              6)  staub rinit   data xhisto +0.1 groesser  SUBROUTINE sub_init_staub
c  
c              7)  diese version berechnet iwc density (Variable densice) u. schreibr raus (steht auf Feld pmseorg


c                                            -> filenamen betd20_n
c                                                         hd20_n, 40md20_n




c  beta, h2o files alle 6 h,  40mio files alle 20 tage
c  version l5 mit H2O const auf druckniveaus -> filenamen beta1_n
c                                                         hm1_n, 40mi1_n
c
c   h2o schmelz  auf 50 % veringert
c
c  vorsicht lyman alpha tag einstellen
c
c    vorsicht ellipsefak (itag) fuer sued ist noch gespiegelt, -> aendern!
c
c    wie test 9 , d.h lyman = solar mean, aber volufak_tab (= stellver 
c                      fuer h2o background aus, kommtschelzone besser raus?
c                      (aber nicht zu gross)
c     +subroutine zenit:  secchi(j,i), zgeo_86, sunris_86   
c      sub zenit jeden 2ten zeitschritt ( shift geloescht)
c
c    ohne volufak_tab (Eis <-> hintergrundwasser)
c    geht es nicht weil sonst der Wasserdampf in der Schmelzregion
c    viel zu gross wird !!
c
c  vorsicht auskommentieren  bei single mode !!!

c ... DO PRIVATE (i,j,k,rwert0,tback,tp,rwert1,work1akt,itemp,
c     & vpsatt_tp,satt_tp,drdt,backgr,volufak,dback,dreferenz)

c  parallel: ifort -openmp -convert big_endian -save -fpp -O2 -o
c
c    ifort -convert big_endian -fpp -O2 -o
c

c----------------- Steuerparameter -------------
#define NSTART    0        ! 0 (null-start)  oder 1 (zwischen)
#define HEMIS     2        ! 1 sued          oder  2 nord
c#define JAHR      2008   ! wird jetzt eingelesen !!!!  
#define MONSTART  5
#define TAGSTART  10

#define MONENDE   8
#define TAGENDE   31

#define JAHRRESTART '1961' ! Jahr des Restarts   vorsicht  !!!!
#define RUNINFO   'run.out'

! Input files

#define CRINIT        '../trac_3d_ini_40mio_l3.dat'
#define H2OINIT      '../trach2o_1dstart_zpress.dat'
#define CTURBKZZ      '/home/oa034/knocluc3d_lima/luebkzza.dat'
#define CBACKSCATTER532 '/home/oa034/knocluc3d_lima/backscatter532.dat'
#define CELLIPSE_EARTH '/home/oa034/kmoda/ellipse_earth.dat'
#define CH4_1861_2135 '/home/oa034/atmos/ch4_annual_MIMAS_0001.tab'
#define UT1_1861_2100 '/home/oa034/atmos/ut1_10_may_MIMAS_MLFU.tab'

#define LIMA_PFAD    '/home/oa034/cop/LC33/LIMA-ICE/backgr/'

! Restart files
#define RESTART_PFAD'/home/optik/dmf/LC33/LIMA-ICE/40mio/'

! Output Pfad
#define OUTPUT_PFAD '/home/oa034/dmf/AV06/MIMAS_RunA_LC33_d19_'

! Statistik Pfad
#define STAT_PFAD  './stat/'

c-----------------------------------------------
c      nstart = 0     ! ##NSTART## 0 (null-start)  oder 1 (zwischen)
c      i_hemi  = 1    ! ##HEMIS##  1 sued       oder  2 nordc

c     c_jahrst='2001' wird nur bei nstart=1 gebraucht in suedhemi !!
c      c_jahrst='2001'           ! ##YEAR## 
                                ! wichtig bei jahreswende zu silvester
      !  auf vergangenes jaht zurueck gesetzt werden, ansonsten =ijahr
      !  checken ob lese file stimmt bei if-anweisung vor haupt schleife
c      i_jahr =  2000         ! ##NYEAR##  2003
c      i_monatstart =  11     ! ##MONSTART## monat 05   also mai   05
c      i_monatende  =  14     ! ##MONENDE## monat 08   also august 08
c      i_tagstart  = 15       ! ##TAGSTART## tag des monats 25 also 25mai
c      i_tagende   = 28       ! ##TAGENDE## tag des monats 25
                             ! (14-15 august letztes eis 2005)

c   f90 +U77 +DD64 +DSnative +Oopenmp +Oparallel +Ovectorize +allow_ivdep +Ofltacc +Ofaster +Ofastaccess -o eis3_20mio_para.exe eis3_20mio_para.f
c   f90 +U77 +DD64 +DSnative +Oall -o eis3_20mioa8.exe eis3_20mioa8.f

c
c   hauptpfad 20mio      alte saettigung schwach                 files 20 mi1,hm1,2001
c             20 mi2     neue mauersberger u. krankowsky stark   files 20 mi1,hm1,2001
c

c  volufak geaendert
c  photolyse raten gibt es nicht an transport boundaries
c            also nur innerhalb k=24,153   j=2,nbneu   i=1,igitneu
c  nord-sued transport neu ab j=2,nbneu  vorher j=1,nbneu
c  vertikaler transport neu ab  k=24,153
c  -->  tunen fuer eis:
c   mauersberger saettigungsformel (neu ist besser!)
c   staubinit
c   turbkzz   ( und lunten)
c   drdt  formel
c   wsedi formel
c   partikel temperatur  -> oberschranke fuer max eisradi

c______________________________________________________________
c    work1  * 1.6    mit 1.6 fuer 2003, 2005 NH, andere jahre (2001) ohne
c    worka  * 1.6    mit 1.6 fuer 2003, 2005 NH
c    1.6 ist aber falsch  !!!!
c______________________________________________________________
c    alter mauersberger
c    particle temp 100 %   (200 % alt)
c    staub *1. und nicht 0.5 doch 0.5 an

c   i_hemi  = 1      sued hemisp
c   i_hemi  = 2      nord hemisp
c default integer size is integer*4 (32 bit) --> bis zu 2*1E9


#define ntrac       40000000
#define kgitzpr     166  !  166 hoehen ab 72,0 72,2 72,4  ... 104,8 105,0 km als druckhoehen
c                           für breite alomar
#define kgitneu     163
#define kgitneu1    164
#define lunten      12        ! Achtung GB hat im code 22 entdeckt, oben war 12 gesetzt
#define igitneu     120
#define igitneu1    121
#define nbneu       53
#define nbneu1      54
#define tg          456
#define pi          3.14159265358979
#define dttrans     90.0                ! zeitschritt = 960 pro tag
#define dttrans_2   180.0
#define rearth      6366197.0
#define gmeso       9.55
#define xkbolz      1.3805e-23
#define xmh2o       2.99e-26
#define xmair       4.845e-26
#define rhoice      932.0
#define NMAX        500  ! Memory size in tridiag
#define TAB_RPTP    200
#define TAB_RMAX    300
#define TAB_TMAX    71
#define beta05      0.24  ! 12.0/50.0
#define t82km       13.0

c-----------------------------------------------------------------------
      MODULE common_h2o

      real*4, dimension(kgitzpr) :: zpress_init,zpress_init_h2o
      real*4, dimension(nbneu,igitneu,kgitneu) :: hm,hminit3d,hmdummy
      real*4, dimension(nbneu,igitneu) :: secchi
      real*4, dimension(nbneu,kgitneu) :: hminit
      real*4, dimension(365) :: ellipsefak,oblinoa,true_solar_time
      END MODULE common_h2o
c-----------------------------------------------------------------------
      MODULE common_const

      real*4, dimension(kgitneu) :: zgeo,turbkzz
      real*4, dimension(nbneu)   :: dx
      real*4, dimension(1500):: backradius,crosssec532, crosssec126
      real*4, dimension(1500):: crosssec200, crosssec1000, crosssec3000
      real*4, dimension(1500):: ext532, ext126, ext200, ext1000,ext3000
      real*4, parameter :: dz = 100.0
      real*4, parameter :: dy = 111111.1111
      END MODULE common_const
c-----------------------------------------------------------------------
      MODULE common_tab

      real*4, dimension(nbneu)    :: horiwicht_tab
      real*4, dimension(kgitneu)  :: wrichtig_tab,volufak_tab
      real*4, dimension(TAB_TMAX) :: vpsatt_tab,sqrt1t_tab,sqrtsedi_tab
      real*4, dimension(TAB_RMAX,TAB_TMAX)    :: faktkelvin_tab
      real*4, dimension(TAB_RPTP,kgitneu)     :: worka_tab
      END MODULE common_tab
c-----------------------------------------------------------------------
      MODULE common_uvwtd

      real*4, dimension(nbneu,igitneu,kgitneu)  :: um,tm,dm
      real*4, dimension(nbneu1,igitneu,kgitneu) :: vm
      real*4, dimension(nbneu,igitneu,kgitneu1) :: wm
      END MODULE common_uvwtd
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      MODULE common_uvwtd48

      real*4, dimension(nbneu,igitneu,kgitneu)  :: um1,um2,
     &                                         tm1,tm2,dm1,dm2,zm1,zm2
      real*4, dimension(nbneu1,igitneu,kgitneu) :: vm1,vm2
      real*4, dimension(nbneu,igitneu,kgitneu1) :: wm1,wm2
      END MODULE common_uvwtd48
c-----------------------------------------------------------------------
      MODULE common_eis

      real*4, dimension(ntrac) :: xfeld,yfeld,zfeld,
     &                              rfeld,rinit
      END MODULE common_eis
c-----------------------------------------------------------------------
      MODULE common_netcdf
      real*4, dimension(kgitneu,igitneu,nbneu,tg) :: aq,bq,cq,dq
      real*4, dimension(kgitneu,igitneu,nbneu,tg) :: eq,fq,gq,hq
      real*4, dimension(kgitneu,igitneu,nbneu,tg) :: iq,jq,kq,lq
      real*4, dimension(kgitneu,igitneu,nbneu,tg) :: mq,nq,oq,pq
      END MODULE common_netcdf
c_______________________________________________________________________



      PROGRAM walcek
c
c    neustart nur bei vollen tagen (0:00 UT) moeglich !!!!
c    siehe eingabe
c
      use OMP_lib
!       USE nrtype
      use ifport
      USE common_const
      USE common_h2o
      USE common_uvwtd
      USE common_uvwtd48
      USE common_eis

      implicit none
!       real*4, dimension(nbneu,igitneu,kgitneu) :: xnluft
      real*4, dimension(igitneu,60)  :: tfeld2d

      real*4, dimension(35)  :: bbreite,bsum   ! neu 
c-----------------------------------------------------------------------
      real*4, dimension(kgitneu,40) ::                 !  local felder
c                                 3 stationen 40 zeitschritte pro h
     &              verti_1_um,verti_1_vm,verti_1_wm,verti_1_tm,
     &              verti_1_hm,verti_1_dm,verti_1_sm,
     &              verti_2_um,verti_2_vm,verti_2_wm,verti_2_tm,
     &              verti_2_hm,verti_2_dm,verti_2_sm,
     &              verti_3_um,verti_3_vm,verti_3_wm,verti_3_tm,
     &              verti_3_hm,verti_3_dm,verti_3_sm
      real*4, dimension(2400000) ::        !maximal 60000x40 zeiten pro h
     &              verti_1_eisradius, verti_1_zfeld, 
     &              verti_2_eisradius, verti_2_zfeld, 
     &              verti_3_eisradius, verti_3_zfeld
      integer, dimension(2400000) ::       !maximal 60000x40 zeiten pro h
     &              nverti_1_eisnr,
     &              nverti_2_eisnr,
     &              nverti_3_eisnr

      real*4,  dimension(kgitneu) :: tmean,zpressmean
      real*4   xlyman_obs(366) 
      real*4   xlyman_const(366) 
      integer, dimension(nbneu)   :: j_mesopause,j_zpress86km

      integer, dimension(40) ::  
     &              nverti_1_anz,nverti_2_anz,nverti_3_anz
c-----------------------------------------------------------------------
      real*4, dimension(1000,960) ::            !  1000er gruppe
     &              grup_x,grup_y,grup_z,grup_r,
     &              grup_tm,grup_hm,grup_sm
c-----------------------------------------------------------------------
      integer record_number
c-----------------------------------------------------------------------
      integer, dimension(12) :: ntage_monat
      integer, dimension(nbneu) :: j_zbound_ob

      integer ijahr_fixdyn,ijahr_fixtemp
      integer :: i_jahr,i_monat,i_monat2,i_tag,i_std,i_hemi,nstart,
     &           i_tagstart,i_tagende,j1,i,j,jj,k,kk,n,n_offset,
     &           i_monatstart, i_monatende, i_jahrst, i_monatst,
     &           i_tagst, i_stdst, monanf,monend,ii,idummy,
     &           itag_im_jahr,ntg,ntg_pro_h,ntot1,ntot2,ntot3,
     &           nlauf,nsterb,nbabies,nloka1,nloka2,nloka3,
     &           nrando,intx,inty,ngrup,nsum,max_doy,kmeso,kk1,kk2
      integer :: ios,iter,i0,i1,i2,j0,j2,mnudge,JAHR
      real*4, dimension(nbneu,igitneu,kgitneu1) :: workw

      real*4 :: dsum, a1,a2,worka,vpsatta,bgrenze,
     &           tminwert, rando, rando1, sum, tmin,
     &           sump, scale_ch4   !!!,d_std, dttt
      
      character*100 charnamegrup,charname1,charname2
      character*100 charnameloca1
      character*100 charnameloca2
      character*100 charnameloca3
      character*2 c_monat(12),c_tag(31),c_stunde(24)
      character*1 c_hemi(2)
      character*4 c_jahr,c_jahrst
      character*20 c_readbuffer
c-----------------------------------------------------------------------
! changes for netcdf output
c-----------------------------------------------------------------------
      integer*4 :: tu
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      real*4  ypos(igitneu)
      real*4  bwidth(nbneu)
      real*4  times(tg)
      real*4   unix_time	   
      character*48 nc_outfile

c-----------------------------------------------------------------------
      tu=0     
c-----------------------------------------------------------------------
      data c_hemi   /'s','n'/
      data c_monat  /'01','02','03','04','05','06','07','08','09','10',
     &               '11','12'/
      data c_tag    /'01','02','03','04','05','06','07','08','09','10',
     &               '11','12','13','14','15','16','17','18','19','20',
     &               '21','22','23','24','25','26','27','28','29','30',
     &               '31'/
      data c_stunde /'01','02','03','04','05','06','07','08','09','10',
     &               '11','12','13','14','15','16','17','18','19','20',
     &               '21','22','23','24'/
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c Achtung ##IDENTIFIER## Zeilen werden durch parameter ersetzt
c-- !$OMP MASTER
      nstart = NSTART     ! 0 (null-start)  oder 1 (zwischen)
      i_hemi  = HEMIS    !  1 sued       oder  2 nord

c     c_jahrst='2001' wird nur bei nstart=1 gebraucht in suedhemi !!
c     c_jahrst=JAHRRESTART           !  
                                ! wichtig bei jahreswende zu silvester
      !  auf vergangenes jaht zurueck gesetzt werden, ansonsten =ijahr
      !  checken ob lese file stimmt bei if-anweisung vor haupt schleife
      i_jahr = 2014 !Grygalashvyly! JAHR         !  2003
      ijahr_fixdyn = 1976
c     ijahr_fixtemp = 1878   !  1. solare min !!
      i_monatstart =  MONSTART     ! monat 05   also mai   05
      i_monatende  =  MONENDE     ! monat 08   also august 08
      i_tagstart  = TAGSTART       ! tag des monats 25 also 25mai
      i_tagende   = TAGENDE       ! tag des monats 25
                             ! (14-15 august letztes eis 2005)
c-----------------------------------------------------------------------
c     ende eingabe !!!
c-----------------------------------------------------------------------
      write(c_jahr,'(I4)') i_jahr
      print*,'c_jahr ',c_jahr
c-----------------------------------------------------------------------
c /home/dmf/oa015/
c oder lokal ./dmf
      
c-----------------------------------------------------------------------
      call set_monate(ntage_monat,max_doy,i_jahr)

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
!       goto 8080
      print*,'lesen vor sub_init_global'
      call sub_init_global (i_jahr,xlyman_obs, scale_ch4, unix_time)
      print*,'lesen nach sub_init_global'
c-----------------------------------------------------------------------
      
      print*,'lesen vor sub_leseh2o_zprinit'
      call sub_leseh2o_zprinit
      print*,'lesen nach sub_leseh2o_zprinit'
      print*,''
      do k=1,kgitzpr
         zpress_init(k)=72.+float(k-1)*0.2
         print*,k,zpress_init(k),zpress_init_h2o(k)
      enddo

c-----------------------------------------------------------------------
      print*,'start nur zum vollen tag : lesen vor sub_lesedyn'
      i_jahrst = i_jahr
      i_monatst = i_monatstart 
      i_tagst   = i_tagstart-1  
      i_stdst   = 24
  
      if (i_tagstart.eq.1) then
       i_monatst=i_monatst-1
       if (i_monatst.eq.0) then
        i_monatst = 12
        i_jahrst = i_jahrst-1
       endif
       i_tagst   = ntage_monat(i_monatst)
      endif

      print*, 'start vor sublesedyn: ',
     &        i_hemi,i_jahrst,i_monatst,i_tagst,i_stdst
      call sub_lesedyn (i_hemi,ijahr_fixdyn,
     &   i_jahrst,i_monatst,i_tagst,i_stdst,um1,tm1,dm1,vm1,wm1,zm1)
c     & ijahr_fixtemp,i_monatst,i_tagst,i_stdst,um1,tm1,dm1,vm1,wm1,zm1)

c      print*,'z press min max: ',minval(zm1),maxval(zm1)

c-----------------------------------------------------------------------
      do iter = 1,5
      workw= wm1

      do k=2,kgitneu1-1
       do i=1,igitneu
       i1=i
       i0=i-1
       if (i.eq.1) i0=igitneu
       i2=i+1
       if (i.eq.igitneu) i2=1

       j=nbneu-4
       j1=j
       j0=j-1
       j2=j+1
       wm1(j,i,k)=(workw(j0,i0,k)+ workw(j1,i0,k)+workw(j2,i0,k)+
     &          workw(j0,i1,k)+ 27.*workw(j1,i1,k)+workw(j2,i1,k)+
     &          workw(j0,i2,k)+ workw(j1,i2,k)+workw(j2,i2,k))/35.

       j=nbneu-5
       j1=j
       j0=j-1
       j2=j+1
       wm1(j,i,k)=(workw(j0,i0,k)+ workw(j1,i0,k)+workw(j2,i0,k)+
     &          workw(j0,i1,k)+ 27.*workw(j1,i1,k)+workw(j2,i1,k)+
     &          workw(j0,i2,k)+ workw(j1,i2,k)+workw(j2,i2,k))/35.

       j=nbneu-6
       j1=j
       j0=j-1
       j2=j+1
       wm1(j,i,k)=(workw(j0,i0,k)+ workw(j1,i0,k)+workw(j2,i0,k)+
     &          workw(j0,i1,k)+ 27.*workw(j1,i1,k)+workw(j2,i1,k)+
     &          workw(j0,i2,k)+ workw(j1,i2,k)+workw(j2,i2,k))/35.

       j=nbneu-7
       j1=j
       j0=j-1
       j2=j+1
       wm1(j,i,k)=(workw(j0,i0,k)+ workw(j1,i0,k)+workw(j2,i0,k)+
     &          workw(j0,i1,k)+ 27.*workw(j1,i1,k)+workw(j2,i1,k)+
     &          workw(j0,i2,k)+ workw(j1,i2,k)+workw(j2,i2,k))/35.

       j=nbneu-8
       j1=j
       j0=j-1
       j2=j+1
       wm1(j,i,k)=(workw(j0,i0,k)+ workw(j1,i0,k)+workw(j2,i0,k)+
     &          workw(j0,i1,k)+ 27.*workw(j1,i1,k)+workw(j2,i1,k)+
     &          workw(j0,i2,k)+ workw(j1,i2,k)+workw(j2,i2,k))/35.

       j=nbneu-9
       j1=j
       j0=j-1
       j2=j+1
       wm1(j,i,k)=(workw(j0,i0,k)+ workw(j1,i0,k)+workw(j2,i0,k)+
     &          workw(j0,i1,k)+ 27.*workw(j1,i1,k)+workw(j2,i1,k)+
     &          workw(j0,i2,k)+ workw(j1,i2,k)+workw(j2,i2,k))/35.

      enddo  ! i
      enddo  ! k
      enddo ! iter
       
      do k=1,kgitneu1
       do i=1,igitneu 
        a1=0.4*wm1(49,i,k)
        a2=0.4*wm1(48,i,k)
        wm1(50,i,k)= 0.2*wm1(50,i,k)+a1+a2
        wm1(51,i,k)= 0.2*wm1(51,i,k)+a1+a2
        wm1(52,i,k)= 0.2*wm1(52,i,k)+a1+a2
        wm1(53,i,k)= 0.2*wm1(53,i,k)+a1+a2
       enddo
      enddo
      print*, 'nach sublesedyn'
      call sub_maxmin (um1,vm1,wm1)      
c_______________________________________________________________________
c    initialisierung staubpartikel
c_______________________________________________________________________
! 8080  continue

c
      do k=1,kgitneu
      zgeo(k) = float(k-1)*0.1+77.8
      enddo

c   mesopausenhoehe + j_zpress86km
      do j =1,nbneu      
      do k=1,kgitneu
       sum = 0.
       sump= 0.
       do i=1,igitneu 
       sum = sum + tm1(j,i,k)
       sump= sump+ zm1(j,i,k)
       enddo
       tmean(k) =      sum/float(igitneu)
       zpressmean(k) = sump/float(igitneu)
      enddo
      kmeso = -9999 
      tmin = minval (tmean)
      do k=1,kgitneu
      if (tmean(k).eq.tmin) kmeso = k 
      enddo
c  neu
      if (kmeso.gt.122) kmeso = 122
      j_mesopause(j) = kmeso
        do k=1,kgitneu-1
        if (zpressmean(k).le. 89.1 .and. zpressmean(k+1).ge. 89.1) 
     &                              j_zpress86km(j) = k
        enddo
c        print*,'hallo zgeo_zpress86km(j): ',j,zgeo(j_zpress86km(j))
      enddo
      
      call sub_init_staub (nstart,bbreite,bsum,j_mesopause,
     &                                         j_zpress86km)
                                      ! rinit,... bei nstart = 0
     
      if (nstart.eq.0)  then
       i_std = 0
 6601  print*,' vor dmfnam_ini : ',i_std,CRINIT

c       goto 6603   !   neu, um schreiben zu verhindern !!

       open (66,err=6602,status='unknown',form='unformatted',
     &      file=CRINIT)
       write (66,err=6602) i_std
       write (66,err=6602) xfeld,yfeld,zfeld,rfeld
       print*,' vor close (66) '
       close (66,err=6602)
       print*,' nach close (66) '
      
       goto 6603

 6602   close(66,iostat=ios)
        print*,'Error writing init: ',CRINIT
        call sleep (2)
        goto 6601

 6603  rinit = rfeld
      endif

c
c wichtig !!!    erster zeitschritt 15.mai
c     erste initialisierung auf null, nur wenn hm(1,1,1)=0.
c     dann wird start hm = hminit3d
c
      hm = 0. ! erste initialisierung auf null,


      call sub_tabelle
      
      itag_im_jahr = 0
      do i=1,i_monatstart
        monanf = 1
        monend= ntage_monat(i)
c        if (i.eq.i_monatstart) monend=i_tagende
c  GB/UB 2007.01.15
        if (i.eq.i_monatstart) monend=i_tagstart
c        print*,' monanf/monend : ',monanf,monend
        do ii=monanf,monend
         itag_im_jahr = itag_im_jahr + 1
        enddo
      enddo
      itag_im_jahr = itag_im_jahr - 1
c      print*,' itag_im_jahr-Start : ',itag_im_jahr 

      do 1002 i_monat2=i_monatstart,i_monatende

c Monatsschleife

      monanf = 1
      i_monat=i_monat2
      if(i_monat2.eq.13) i_jahr=i_jahr+1
      write(c_jahr,'(I4)') i_jahr

      if(i_monat.gt.12) i_monat=i_monat-12
      monend= ntage_monat(i_monat)

      if (i_monat.eq.i_monatstart) monanf=i_tagstart
      if (i_monat2.eq.i_monatende)  monend=i_tagende

      do 1001 i_tag= monanf,monend
c______________________________
c Tagesschleife
c
      itag_im_jahr = itag_im_jahr + 1
c      goto 2222
      if (itag_im_jahr.gt.max_doy) itag_im_jahr=1
      call set_monate(ntage_monat,max_doy,i_jahr)
      print *,' itag_im_jahr : ', itag_im_jahr ,i_jahr,'-',
     &     i_monat,'-',i_tag

c Zufallszahlgenerator-Initialisierung:
c Jeden Tag wird die seed auf die Tagesnummer im Jahr gesetzt, um
c einerseits die Qualität der Zufallszahlen nicht zu stören, und
c andererseits Wiederholbarkeit im Fall eines Restarts zu gewährleisten

	call srand(itag_im_jahr)

      grup_x = 0.
      grup_y = 0.
      grup_z = 0.
      grup_r = 0.
      grup_tm = 0.
      grup_hm = 0.
      grup_sm = 0.

      ntg = 0

      do 1024 i_std=1,24
CC              Grygalashvyly
CC     call sub_dmget (i_hemi,i_jahr,i_monat,i_tag)

c Stundenschleife
c-------------------------  def stündliches hminit --------------------
      um = um1 
      tm = tm1 
      dm = dm1 
      vm = vm1 
      wm = wm1

      call sub_h2oinit_zpr_to_zgeo (zm1) ! berechnet hminit3d

c   neu bei era interim runs !!!
      hminit3d = hminit3d *0.9   !  minus 10 prozent

c----------------------------------------------------------------------
c CH4: hier wird h2o erniedrigt WICHTIG !!!
c----------------------------------------------------------------------
      hminit3d = hminit3d * scale_ch4   ! scale_ch4 in prozent/100.

c
c wichtig !!! siehe oben wo start hm auf null steht
c
      if (hm(1,1,1).eq.0.) hm = hminit3d
   
      call sub_zenit (1,itag_im_jahr,i_hemi)
      

      hmdummy = hm
      hm = hminit3d*1.2  ! der faktor (15 proz) um solar min zu kompensieren
c
c    sub_photolyse rechnet mit hm, deshalb die kopie mit hmdummy
c    hier soll ja hminit3d berechnet werden,
c    deshalb auch mnudge=0, sonnendauer 100000. sec

      mnudge = 0
      xlyman_const = 4.46597 !4.75
c      call sub_photolyse (mnudge,100000.,itag_im_jahr,i_hemi,xlyman_obs)
      call sub_photolyse (mnudge,100000.,itag_im_jahr,i_hemi,
     &                    xlyman_const)

      hminit3d = hm
      hm = hmdummy

      do k=1,kgitneu
       do j=1,nbneu
        dsum=0.
        do i=1,igitneu
         dsum = dsum + hminit3d(j,i,k)
        enddo
        hminit(j,k) = dsum/float(igitneu)
       enddo
      enddo
      print*,' '
      do k=23,kgitneu,10
       write(6,1000) k,(hminit(j,k),j=1,53,5),k
      enddo
 1000 format(i4,2x,11(f6.3,1x),2x,i3)
      
      ntg_pro_h = 0
      ntot1 = 0
      ntot2 = 0
      ntot3 = 0

      verti_1_um = 0.
      verti_1_vm = 0.
      verti_1_wm = 0.
      verti_1_tm = 0.
      verti_1_hm = 0.
      verti_1_dm = 0.
      verti_1_sm = 0.
      verti_1_eisradius = 0.
      verti_1_zfeld = 0.
      nverti_1_anz = 0 
      nverti_1_eisnr = 0

      verti_2_um = 0.
      verti_2_vm = 0.
      verti_2_wm = 0.
      verti_2_tm = 0.
      verti_2_hm = 0.
      verti_2_dm = 0.
      verti_2_sm = 0.
      verti_2_eisradius = 0.
      verti_2_zfeld = 0.
      nverti_2_anz = 0 
      nverti_2_eisnr = 0

      verti_3_um = 0.
      verti_3_vm = 0.
      verti_3_wm = 0.
      verti_3_tm = 0.
      verti_3_hm = 0.
      verti_3_dm = 0.
      verti_3_sm = 0.
      verti_3_eisradius = 0.
      verti_3_zfeld = 0.
      nverti_3_anz = 0 
      nverti_3_eisnr = 0

c      print*,'jetzt lesedyn mit Datum: ',i_jahr,i_monat,i_tag,i_std
c_______________________________________________________________________
      call sub_lesedyn (i_hemi,ijahr_fixdyn,
     &       i_jahr,i_monat,i_tag,i_std,um2,tm2,dm2,vm2,wm2,zm2)
c     & ijahr_fixtemp,i_monat,i_tag,i_std,um2,tm2,dm2,vm2,wm2,zm2)

      do iter = 1,5
      workw= wm2

      do k=2,kgitneu1-1
       do i=1,igitneu
       i1=i
       i0=i-1
       if (i.eq.1) i0=igitneu
       i2=i+1
       if (i.eq.igitneu) i2=1

       j=nbneu-4
       j1=j
       j0=j-1
       j2=j+1
       wm2(j,i,k)=(workw(j0,i0,k)+ workw(j1,i0,k)+workw(j2,i0,k)+
     &          workw(j0,i1,k)+ 27.*workw(j1,i1,k)+workw(j2,i1,k)+
     &          workw(j0,i2,k)+ workw(j1,i2,k)+workw(j2,i2,k))/35.

       j=nbneu-5
       j1=j
       j0=j-1
       j2=j+1
       wm2(j,i,k)=(workw(j0,i0,k)+ workw(j1,i0,k)+workw(j2,i0,k)+
     &          workw(j0,i1,k)+ 27.*workw(j1,i1,k)+workw(j2,i1,k)+
     &          workw(j0,i2,k)+ workw(j1,i2,k)+workw(j2,i2,k))/35.

       j=nbneu-6
       j1=j
       j0=j-1
       j2=j+1
       wm2(j,i,k)=(workw(j0,i0,k)+ workw(j1,i0,k)+workw(j2,i0,k)+
     &          workw(j0,i1,k)+ 27.*workw(j1,i1,k)+workw(j2,i1,k)+
     &          workw(j0,i2,k)+ workw(j1,i2,k)+workw(j2,i2,k))/35.

       j=nbneu-7
       j1=j
       j0=j-1
       j2=j+1
       wm2(j,i,k)=(workw(j0,i0,k)+ workw(j1,i0,k)+workw(j2,i0,k)+
     &          workw(j0,i1,k)+ 27.*workw(j1,i1,k)+workw(j2,i1,k)+
     &          workw(j0,i2,k)+ workw(j1,i2,k)+workw(j2,i2,k))/35.

       j=nbneu-8
       j1=j
       j0=j-1
       j2=j+1
       wm2(j,i,k)=(workw(j0,i0,k)+ workw(j1,i0,k)+workw(j2,i0,k)+
     &          workw(j0,i1,k)+ 27.*workw(j1,i1,k)+workw(j2,i1,k)+
     &          workw(j0,i2,k)+ workw(j1,i2,k)+workw(j2,i2,k))/35.

       j=nbneu-9
       j1=j
       j0=j-1
       j2=j+1
       wm2(j,i,k)=(workw(j0,i0,k)+ workw(j1,i0,k)+workw(j2,i0,k)+
     &          workw(j0,i1,k)+ 27.*workw(j1,i1,k)+workw(j2,i1,k)+
     &          workw(j0,i2,k)+ workw(j1,i2,k)+workw(j2,i2,k))/35.

      enddo  ! i
      enddo  ! k
      enddo ! iter

      do k=1,kgitneu1
      do i=1,igitneu 
      a1=0.4*wm2(49,i,k)
      a2=0.4*wm2(48,i,k)
      wm2(50,i,k)= 0.2*wm2(50,i,k)+a1+a2
      wm2(51,i,k)= 0.2*wm2(51,i,k)+a1+a2
      wm2(52,i,k)= 0.2*wm2(52,i,k)+a1+a2
      wm2(53,i,k)= 0.2*wm2(53,i,k)+a1+a2
      enddo
      enddo

      call sub_maxmin (um2,vm2,wm2) 
c______________________________
      do 998 idummy=1, 20   !1,2    !1,5

c______________________________
      a1 = float(mod(ntg,20))   !  40 pro stunde
      a1 = 1. - a1/20.
      a2 = 1. - a1
      um = a1*um1 + a2*um2
      tm = a1*tm1 + a2*tm2
      dm = a1*dm1 + a2*dm2
      vm = a1*vm1 + a2*vm2
      wm = a1*wm1 + a2*wm2
c______________________________


c  erster zeitschritt
      ntg = ntg + 1
      ntg_pro_h = ntg_pro_h + 1
      hm = amax1(0.01,hm)
      call sub_transwalcek (dx,dy,dz,um,vm,wm,hm)
c
c  zweiter zeitschritt

      ntg = ntg + 1
      ntg_pro_h = ntg_pro_h + 1
      hm = amax1(0.01,hm)
      call sub_transwalcek (dx,dy,dz,um,vm,wm,hm)
c  
      call sub_diffu_h2o (hm,dz,turbkzz)
      call sub_zenit (ntg,itag_im_jahr,i_hemi)
      mnudge = 1
      call sub_photolyse (mnudge, dttrans_2, itag_im_jahr,
     &                    i_hemi,xlyman_obs)

c_______________________________________________________________________
c    eis mikrophsik und eis transport
c_______________________________________________________________________
      call sub_tracertransp     ! also jetzt dttrans_2 = dttrans * 2.
                                ! plus wrando geviertelt
c_________________________________________
c    geburtenkontrolle, zurueckwerfen der staubteilchen
c_________________________________________
      bgrenze = 16.   !  neu  (alt=3.5)
c      print*,' bgrenze: ',bgrenze

      nbabies = 0
      nloka1=0
      nloka2=0
      nloka3=0
      ngrup=0

c      j_zbound_ob = 92
      j_zbound_ob = 100  ! + 800 m
      do j= 1,10
      jj = 43 + j
c      j_zbound_ob (jj) = 133
      j_zbound_ob (jj) = 141     ! + 800 m
      jj = 33 + j
c      j_zbound_ob (jj) = 122 + j
      j_zbound_ob (jj) = 130 + j   ! + 800 m
      jj = 23 + j
c      j_zbound_ob (jj) = 112 + j
      j_zbound_ob (jj) = 120 + j   ! + 800 m
      jj = 13 + j
c      j_zbound_ob (jj) = 102 + j
      j_zbound_ob (jj) = 110 + j   ! + 800 m
      jj =  3 + j
c      j_zbound_ob (jj) =  92 + j
      j_zbound_ob (jj) = 100 + j   ! + 800 m
      enddo
c_______________________________________________________________________
c
c   mesopausenhoehe
      do j =1,nbneu
      do k=1,kgitneu
       sum = 0.
       do i=1,igitneu
       sum = sum + tm(j,i,k)
       enddo
       tmean(k) = sum/float(igitneu)
      enddo
      kmeso = -9999
      tmin = minval (tmean)
      do k=1,kgitneu
      if (tmean(k).eq.tmin) kmeso = k
      enddo
c  neu
      if (kmeso.gt.122) kmeso = 122
      j_mesopause(j) = kmeso
      enddo
c_______________________________________________________________________
c  alle 6 h (20*6 = 120) sind alle teilchen abgefragt worden
      n_offset = int(rand()*120)+1

      do 2001 n = n_offset,ntrac,120  !  ->  Abfrage zufaelliger pro 120

      if (rfeld(n).eq.rinit(n)) then
c     xfeld bleibt gleich

      rando = rand()*100.       ! breitenkreis zwischen 55 - 90 grad
      do i=1,35
          if (rando .le. bsum(i)) then
            rando1=rand()          !  0-1
            yfeld(n)= bbreite(i)+rando1
            if (yfeld(n).gt.53.9999) yfeld(n) = 53.999
            goto 3332
          endif
      enddo
c      print*,' rando - grenze FEHLER !!! ',rando
 3332    continue


      j = yfeld(n)
      rando = rand()  ! 77.8 km=1, 80 km=23, 83.5 km=58
c      zfeld(n) = 73. + 50. * rando !  neu: 85 - 90 km
c      zfeld(n)=23.+60. +  10.* rando !  neu: 86 - 87 km
c      zfeld(n)= float(j_mesopause(j)) -  10.* rando !  neu: 1 km unterhalb der mesopa
c      zfeld(n)= float(j_mesopause(j)) +10. -  20.* rando !  neu: +-1 km um mesopause
c   neu
c   neu
         if (rando.gt. 0.8) then
            rando = (rando - 0.8)*5.                      !  wieder zwischen 0:1
c            kk1 = j_mesopause(j) - (23+55+10)            ! mesop - 86.5 km
            kk1 = j_mesopause(j) - (5+j_zpress86km(j))    ! mesop - 86.5 km
            kk2 = 15
            if (kk1.ge.1) kk2 = kk2 + kk1
c            zfeld(n)= 23.+57. +  float(kk2)* rando       !  neu: 1 km unterhalb der mesopa
            zfeld(n)= float(j_zpress86km(j)) - 3. +  
     &                   float(kk2)* rando                !  neu: 1 km unterhalb der mesopa
         else
         rando = rando * 1./0.8                           ! wieder 0:1
c         zfeld(n)=23.+55. +  10.* rando                  !  neu: 85.5 - 86.5 km
         zfeld(n)= float(j_zpress86km(j)) - 5. +  10.* rando !  neu: 85.5 - 86.5 km
         endif
c  ende neu      


      endif
c            staub - geburte
 2001 continue
c_______________________________________________________________________

c   klappt noch nicht wegen nbabies !!!
c!$OMP PARALLEL DO PRIVATE (i,j,k,n,ii,nbabies,rando,rando1)

      do 2002 n=1,ntrac     !  ->  Abfrage aller 40 Millionen
c            staub - geburte
c_______________________________________________________________________
      if (rfeld(n).eq.rinit(n)) then !  neu

      k=int(zfeld(n))
      j=int(yfeld(n))
      if (yfeld(n) .le. bgrenze .or.
     &    k .le. 62 .or. k.gt. j_zbound_ob(j)) then   ! neu

c-----------------------------------------------------------------------

      nbabies=nbabies+1


      rfeld(n)=rinit(n)

 3335    continue
         rando = float(igitneu)*rand()
         xfeld(n)=1.+rando     ! 120 laemgen =2 mal 60
         ii=int(xfeld(n))
         if(ii.gt.igitneu.or.ii.lt.1) goto 3335


      rando = rand()*100.       ! breitenkreis zwischen 55 - 90 grad
      do i=1,35
          if (rando .le. bsum(i)) then
            rando1=rand()          !  0-1
            yfeld(n)= bbreite(i)+rando1
            if (yfeld(n).gt.53.9999) yfeld(n) = 53.999
            goto 3333
          endif
      enddo
c      print*,' rando - grenze FEHLER !!! ',rando
 3333    continue

      j = yfeld(n)
      rando = rand()  ! 77.8 km=1, 80 km=23, 83.5 km=58
c      zfeld(n) = 73. + 50. * rando !  neu: 85 - 90 km
c      zfeld(n)=23.+60. +  10.* rando !  neu: 86 - 87 km
c      zfeld(n)= float(j_mesopause(j)) -  10.* rando !  neu: 2 km unterhalb der mesopa
c      zfeld(n)= float(j_mesopause(j)) +15. -  30.* rando !  neu: +-1.5 km um mesopa
c   neu
c   neu
         if (rando.gt. 0.8) then
            rando = (rando - 0.8)*5.                      !  wieder zwischen 0:1
c            kk1 = j_mesopause(j) - (23+55+10)            ! mesop - 86.5 km
            kk1 = j_mesopause(j) - (5+j_zpress86km(j))    ! mesop - 86.5 km
            kk2 = 15
            if (kk1.ge.1) kk2 = kk2 + kk1
c            zfeld(n)= 23.+57. +  float(kk2)* rando       !  neu: 1 km unterhalb der mesopa
            zfeld(n)= float(j_zpress86km(j)) - 3. +  
     &                   float(kk2)* rando                !  neu: 1 km unterhalb der mesopa
         else
         rando = rando * 1./0.8                           ! wieder 0:1
c         zfeld(n)=23.+55. +  10.* rando !  neu: 85.5 - 86.5 km
         zfeld(n)= float(j_zpress86km(j)) - 5. +  10.* rando !  neu: 85.5 - 86.5 km
         endif
c  ende neu
        

      endif
      endif                          ! neu
 2002 continue
c!$OMP END PARALLEL DO

CC     print*,ntg,' geburtenkontrolle: ',nbabies,
CC     &           ' bgrenze: ',bgrenze



  998 continue    !  5er
c_______________________________________________________________________
c     rausschreiben    der wasserdamp alle 6h files

      if (i_std.eq.6 .or. i_std.eq.12 .or. i_std.eq.18 .or. i_std.eq.24)
     &    then
      tu = tu+1
      unix_time = unix_time +21600.
      times(tu)= unix_time
      call schreibe_beta0 (tu)

      endif

c      umbesetzen background
c_______________________________________________________________________
          um1 = um2
          tm1 = tm2
          dm1 = dm2
          vm1 = vm2
          wm1 = wm2

c-----------------------------------------------------------------------
 1024 continue   !  24 h in einem tag
c Ende der Tagesschleife
c      call sub_dmput (i_hemi,i_jahr,i_monat,i_tag)
c_______________________________________________________________________

 1001 continue   !  tage im monat
 1002 continue   !  monats schliefe
c_______________________________________________________________________

      do j=1,nbneu
        bwidth(j)=37.5+float(j-1)
      enddo

      do i=1,igitneu
        ypos(i) = float(i-1)*3.
      enddo
 		
	nc_outfile= OUTPUT_PFAD//c_jahr//'.nc'
 
! write netcdf file for year 
         call writegrid(nc_outfile,bwidth,ypos,zgeo,times)
c_______________________________________________________________________

 8888 format ('zeitschritt : ',4i8,'   cpu: ',3x,f6.1)
c-- !$OMP END MASTER
      stop
      END PROGRAM


!WRITEGRID - write a netCDF gridfile
!:=========================================================================
      SUBROUTINE writegrid(outfile,xpos,ypos,zpos,time)
      USE netcdf
      USE common_netcdf
      IMPLICIT NONE
      REAL(KIND=4), DIMENSION(nbneu), INTENT(IN) :: xpos
      REAL(KIND=4), DIMENSION(igitneu), INTENT(IN) :: ypos
      REAL(KIND=4), DIMENSION(kgitneu), INTENT(IN) :: zpos
      REAL(KIND=4), DIMENSION(tg), INTENT(IN) :: time
	  
      INTEGER :: t_dimid
      INTEGER :: t_varid
 	  
      INTEGER(KIND=4) :: ncid, x_dimid, y_dimid, z_dimid
      INTEGER(KIND=4) :: x_varid, y_varid, z_varid
      INTEGER(KIND=4) :: varid1, varid2, varid3, varid4
      INTEGER(KIND=4) :: varid5, varid6, varid7
      INTEGER(KIND=4) :: varid8, varid9, varid10, varid11
      INTEGER(KIND=4) :: varid12, varid13, varid14, varid15
      INTEGER(KIND=4) :: varid16

      INTEGER(KIND=4), DIMENSION(4) :: dimids1, dimids2
      INTEGER(KIND=4), DIMENSION(4) :: dimids3, dimids4
      INTEGER(KIND=4), DIMENSION(4) :: dimids5, dimids6
      INTEGER(KIND=4), DIMENSION(4) :: dimids7, dimids8
      INTEGER(KIND=4), DIMENSION(4) :: dimids9, dimids10
      INTEGER(KIND=4), DIMENSION(4) :: dimids11,dimids12
      INTEGER(KIND=4), DIMENSION(4) :: dimids13,dimids14
      INTEGER(KIND=4), DIMENSION(4) :: dimids15,dimids16 
	  
      CHARACTER(LEN=48), INTENT(IN) :: outfile

!Create the netCDF file.
      CALL check(nf90_create(outfile, NF90_NETCDF4, ncid))

!Define the dimensions.
      CALL check(nf90_def_dim(ncid, "height", kgitneu, z_dimid))
      CALL check(nf90_def_dim(ncid, "lon", igitneu, y_dimid))
      CALL check(nf90_def_dim(ncid, "lat", nbneu, x_dimid))
      CALL check(nf90_def_dim(ncid, "time", tg, t_dimid))

!Define coordinate variables
      CALL check(nf90_def_var(ncid, "height", NF90_REAL, z_dimid, 
     & z_varid))
      CALL check(nf90_def_var(ncid,"lon",NF90_REAL,y_dimid,y_varid))
      CALL check(nf90_def_var(ncid,"lat",NF90_REAL,x_dimid,x_varid))
      CALL check(nf90_def_var(ncid,"time",NF90_REAL,t_dimid,t_varid))	  


      dimids1 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids2 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids3 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids4 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids5 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids6 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids7 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids8 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids9 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids10 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids11 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids12 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids13= (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids14= (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids15 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)
      dimids16 = (/ z_dimid, y_dimid, x_dimid, t_dimid /)


!Define variable
      CALL check(nf90_def_var(ncid,"beta532",NF90_FLOAT,dimids1,
     &  varid1))
      CALL check(nf90_def_var_deflate(ncid, varid1,
     &  shuffle = 1, deflate = 1, deflate_level = 5))	 
	 
      CALL check(nf90_def_var(ncid,"densice",NF90_FLOAT,dimids2,
     &  varid2))
      CALL check(nf90_def_var_deflate(ncid, varid2,
     &  shuffle = 1, deflate = 1, deflate_level = 5))		 
	 
      CALL check(nf90_def_var(ncid,"r_ice",NF90_FLOAT,dimids3,
     &  varid3))
       CALL check(nf90_def_var_deflate(ncid, varid3,
     &  shuffle = 1, deflate = 1, deflate_level = 5))
	 
      CALL check(nf90_def_var(ncid,"n_ice",NF90_FLOAT,dimids4,
     &  varid4))
       CALL check(nf90_def_var_deflate(ncid, varid4,
     &  shuffle = 1, deflate = 1, deflate_level = 5))	 
	 
      CALL check(nf90_def_var(ncid,"n_dust",NF90_FLOAT,dimids5,
     &  varid5))    
       CALL check(nf90_def_var_deflate(ncid, varid5,
     &  shuffle = 1, deflate = 1, deflate_level = 5))	
	 
      CALL check(nf90_def_var(ncid,"r_ice_eff",NF90_FLOAT,dimids6,
     &  varid6))  
       CALL check(nf90_def_var_deflate(ncid, varid6,
     &  shuffle = 1, deflate = 1, deflate_level = 5))	

      CALL check(nf90_def_var(ncid,"H2O",NF90_FLOAT,dimids7,
     &  varid7))  
       CALL check(nf90_def_var_deflate(ncid, varid7,
     &  shuffle = 1, deflate = 1, deflate_level = 5))

      CALL check(nf90_def_var(ncid,"zon_win",NF90_FLOAT,dimids8,
     &  varid8))
       CALL check(nf90_def_var_deflate(ncid, varid8,
     &  shuffle = 1, deflate = 1, deflate_level = 5))	 
	 
      CALL check(nf90_def_var(ncid,"mer_win",NF90_FLOAT,dimids9,
     &  varid9))    
       CALL check(nf90_def_var_deflate(ncid, varid9,
     &  shuffle = 1, deflate = 1, deflate_level = 5))	
	 
      CALL check(nf90_def_var(ncid,"tem",NF90_FLOAT,dimids10,
     &  varid10))  
       CALL check(nf90_def_var_deflate(ncid, varid10,
     &  shuffle = 1, deflate = 1, deflate_level = 5))	

      CALL check(nf90_def_var(ncid,"press_alt",NF90_FLOAT,dimids11,
     &  varid11))  
       CALL check(nf90_def_var_deflate(ncid, varid11,
     &  shuffle = 1, deflate = 1, deflate_level = 5))

      CALL check(nf90_def_var(ncid,"ext_126",NF90_FLOAT,dimids12,
     &  varid12))  
       CALL check(nf90_def_var_deflate(ncid, varid12,
     &  shuffle = 1, deflate = 1, deflate_level = 5))

      CALL check(nf90_def_var(ncid,"ext_200",NF90_FLOAT,dimids13,
     &  varid13))  
       CALL check(nf90_def_var_deflate(ncid, varid13,
     &  shuffle = 1, deflate = 1, deflate_level = 5))

      CALL check(nf90_def_var(ncid,"ext_532",NF90_FLOAT,dimids14,
     &  varid14))  
       CALL check(nf90_def_var_deflate(ncid, varid14,
     &  shuffle = 1, deflate = 1, deflate_level = 5))

      CALL check(nf90_def_var(ncid,"ext_1000",NF90_FLOAT,dimids15,
     &  varid15))  
       CALL check(nf90_def_var_deflate(ncid, varid15,
     &  shuffle = 1, deflate = 1, deflate_level = 5))

      CALL check(nf90_def_var(ncid,"ext_3000",NF90_FLOAT,dimids16,
     &  varid16))  
       CALL check(nf90_def_var_deflate(ncid, varid16,
     &  shuffle = 1, deflate = 1, deflate_level = 5))	 

!Define unis,standard_name and axis
	  
      CALL check(nf90_put_att(ncid,x_varid,'units','degree_north'))
      CALL check(nf90_put_att(ncid,x_varid,'standard_name','latitude'))
      CALL check(nf90_put_att(ncid,x_varid,'axis','X'))	  
c-----------------------------------------------------------------------	  
      CALL check(nf90_put_att(ncid,y_varid,'units','degree_east'))
      CALL check(nf90_put_att(ncid,y_varid,'standard_name','longitude'))
      CALL check(nf90_put_att(ncid,y_varid,'axis','Y'))
	  
      CALL check(nf90_put_att(ncid, z_varid, 'units', 'km'))
	  
      CALL check(nf90_put_att(ncid, z_varid,
     &  'standard_name', 'geometric height'))  
      CALL check(nf90_put_att(ncid, z_varid, 'axis', 'Z'))	  
	  
	  
      CALL check(nf90_put_att(ncid, t_varid,
     &  'units','seconds since 1970-01-01 00:00:00'))  
      CALL check(nf90_put_att(ncid, t_varid, 'axis', 'T'))	      

      CALL check(nf90_put_att(ncid,varid1,'units','10E-10m^-1sr^-1'))
      CALL check(nf90_put_att(ncid,varid2,'units','g/km^3'))	  
      CALL check(nf90_put_att(ncid,varid3,'units','nm'))		
      CALL check(nf90_put_att(ncid,varid4,'units','cm^-3'))		
      CALL check(nf90_put_att(ncid,varid5,'units','cm^-3'))
      CALL check(nf90_put_att(ncid,varid6,'units','nm'))	 
      CALL check(nf90_put_att(ncid,varid7,'units','ppmv'))	  
      CALL check(nf90_put_att(ncid,varid8,'units','m/s'))		
      CALL check(nf90_put_att(ncid,varid9,'units','m/s'))
      CALL check(nf90_put_att(ncid,varid10,'units','K'))	 
      CALL check(nf90_put_att(ncid,varid11,'units','km'))
      CALL check(nf90_put_att(ncid,varid12,'units','m^-1'))
      CALL check(nf90_put_att(ncid,varid13,'units','m^-1'))
      CALL check(nf90_put_att(ncid,varid14,'units','m^-1'))
      CALL check(nf90_put_att(ncid,varid15,'units','m^-1'))
      CALL check(nf90_put_att(ncid,varid16,'units','m^-1'))

      CALL check(nf90_enddef(ncid)) !End Definitions

!Write Data
      CALL check(nf90_put_var(ncid, t_varid, time))
      CALL check(nf90_put_var(ncid, y_varid, ypos))
      CALL check(nf90_put_var(ncid, x_varid, xpos))
      CALL check(nf90_put_var(ncid, z_varid, zpos))	  
	  
      CALL check(nf90_put_var(ncid, varid1, aq))
      CALL check(nf90_put_var(ncid, varid2, bq))
      CALL check(nf90_put_var(ncid, varid3, cq))
      CALL check(nf90_put_var(ncid, varid4, dq))
      CALL check(nf90_put_var(ncid, varid5, eq))
      CALL check(nf90_put_var(ncid, varid6, fq))
      CALL check(nf90_put_var(ncid, varid7, gq))
      CALL check(nf90_put_var(ncid, varid8, hq))
      CALL check(nf90_put_var(ncid, varid9, iq))
      CALL check(nf90_put_var(ncid, varid10,jq))
      CALL check(nf90_put_var(ncid, varid11,kq))
      CALL check(nf90_put_var(ncid, varid12,lq))
      CALL check(nf90_put_var(ncid, varid13,mq))
      CALL check(nf90_put_var(ncid, varid14,nq))
      CALL check(nf90_put_var(ncid, varid15,oq))
      CALL check(nf90_put_var(ncid, varid16,pq))

      CALL check(nf90_close(ncid))

      END SUBROUTINE writegrid
!:=========================================================================


!Check (ever so slightly modified from www.unidata.ucar.edu)
!:======================================================================
      SUBROUTINE check(istatus)
      USE netcdf
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: istatus
      IF (istatus /= nf90_noerr) THEN
      write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
      END IF
      END SUBROUTINE check
!:======================================================================


c_______________________________________________________________________
c
      SUBROUTINE schreibe_beta0(tu)

      USE common_const
      USE common_eis
      USE common_h2o
      USE common_uvwtd
      USE common_uvwtd48
      USE common_netcdf

      implicit none
CC      integer*4 i_jahr,i_monat,i_tag,i_hemi
      integer*4 j,i,k,n,neis,nstaub,icross,ios,tu
CC      integer record_number, varnum, ngridp
      real*4 volufak,xpi,radi3,radi2,volume,rhoicecgs !,d_std

      real*4 xneis(nbneu,igitneu,kgitneu)
      real*4 xnstaub(nbneu,igitneu,kgitneu)
      real*4 eisradius(nbneu,igitneu,kgitneu)
      real*4 eisradiuseff(nbneu,igitneu,kgitneu)
      real*4 reff_zaehler(nbneu,igitneu,kgitneu) 
      real*4 reff_nenner(nbneu,igitneu,kgitneu)
      real*4 betaorg(nbneu,igitneu,kgitneu) 
      real*4 ext_coef126(nbneu,igitneu,kgitneu)
      real*4 ext_coef200(nbneu,igitneu,kgitneu)
      real*4 ext_coef532(nbneu,igitneu,kgitneu)
      real*4 ext_coef1000(nbneu,igitneu,kgitneu)
      real*4 ext_coef3000(nbneu,igitneu,kgitneu)

      real*4 pmseorg(nbneu,igitneu,kgitneu)
      real*4 densice(nbneu,igitneu,kgitneu)

      real*4 bbreite(nbneu),bgewicht(nbneu)

CC      character*4 c_jahr
CC      character*6 c_stunde,c_std
CC      character*2 c_monat,c_tag,c_hemi
CC      character*120 charname2
C----*-----------------------------------------------------------------*
CC     varnum = 11
CC      ngridp =  nbneu *igitneu *kgitneu
C----*-----------------------------------------------------------------*
      print*,' ---  vor geschrieben beta0 ----  '  !!,d_std
c_______________________________________________________________________
      do k=1,kgitneu
      zgeo(k) = float(k-1)*0.1+77.8
      enddo

      do j=1,nbneu
      bbreite(j)=37.5+float(j-1)
c      bbreite(j)=36.5+float(j-1)  !  bgewicht um eine Breite zurueck
                                  !  (um starken polwert zu vermeiden)
      bgewicht(j) = 1./cos(bbreite(j)*pi/180.)
c      print*,j,bbreite(j),bgewicht(j)
      enddo
c----------------------------------------------------------------------
      ext_coef126 = 0.
      ext_coef200 = 0.
      ext_coef532 = 0.
      ext_coef1000 = 0.
      ext_coef3000 = 0.
      betaorg = 0.
      pmseorg = 0.
      densice = 0.
      eisradius = 0.
      eisradiuseff = 0.
        reff_zaehler = 0.
        reff_nenner = 0.
      xneis = 0.
      xnstaub = 0.
      neis=0
      nstaub=0
c-----------------------------------------------------------------------
c  start
c-----------------------------------------------------------------------
       xpi     =     3.14159265358979
       rhoicecgs=932.0*1.e-3  !   in g / cm**3

      do n=1,ntrac
        if (rfeld(n).gt.rinit(n))  then

        neis=neis+1
        i = aint(xfeld(n))
        j = aint(yfeld(n))
        k = aint(zfeld(n))
        xneis(j,i,k)=xneis(j,i,k) + 1.
        eisradius(j,i,k) = eisradius(j,i,k) + rfeld(n)
          radi2 = rfeld(n)*rfeld(n)
          radi3 = rfeld(n) * radi2
          reff_zaehler(j,i,k) = reff_zaehler(j,i,k) + radi3
          reff_nenner(j,i,k)  = reff_nenner(j,i,k)  + radi2
          volume =  (4./3.)*xpi*radi3*1.e-21                         ! cgs system 
          densice(j,i,k) = densice(j,i,k) + volume*rhoicecgs*1.e15      !  in g/km**3 

        icross = int(rfeld(n)*10.+1.1)
        if (icross.gt.1500)  icross=1500
        betaorg(j,i,k) = betaorg(j,i,k) + 1.e16*crosssec532(icross)
        pmseorg(j,i,k) = pmseorg(j,i,k) + (rfeld(n) * rfeld(n))
        ext_coef126(j,i,k)=ext_coef126(j,i,k)+ 1.e6*ext126(icross)
        ext_coef200(j,i,k)=ext_coef200(j,i,k)+ 1.e6*ext200(icross)
        ext_coef532(j,i,k)=ext_coef532(j,i,k)+ 1.e6*ext532(icross)
        ext_coef1000(j,i,k)=ext_coef1000(j,i,k)+ 1.e6*ext1000(icross)
        ext_coef3000(j,i,k)=ext_coef3000(j,i,k)+ 1.e6*ext3000(icross)
        else
        nstaub = nstaub + 1
        i = aint(xfeld(n))
        j = aint(yfeld(n))
        k = aint(zfeld(n))
        xnstaub(j,i,k) = xnstaub(j,i,k) + 1.
        endif
      enddo
      print*,' anzahl staub teilchen: ',nstaub
      print*,' anzahl eis teilchen: ',neis

c eisradi mean  + umrechnen auf 1/cm**3
c      exp90km = exp(-88./4.0)   ! k=123 orginal code =90 km
c-----------------------------------------------------------------------
c     volufak :  vertikale gewichte
c    plus gewicht_kzz fuer pmse  (spaeter in plot routine)
c     82 km: 10 m**2/s    ->   89 km 185 m**2/s  (Raketen Luebken)
c-----------------------------------------------------------------------
      do k = 1,kgitneu
      volufak =  1. ! 1./exp((88.-zgeo(k))/4.0)   ! also ausgeschaltet
      volufak =  volufak * float(igitneu)

      do i = 1,igitneu
      do j = 1,nbneu
      if (xneis(j,i,k).ge.1.)
     &              eisradius(j,i,k) = eisradius(j,i,k)/xneis(j,i,k)
c
c  60 als skalierung ist so gewaehlt, dass dass fuer initstaub
c  hunten werte folgen  !!!!!
c      xneis(j,i,k)=xneis(j,i,k)* bgewicht(j)*volufak/60.
c      xnstaub(j,i,k)=xnstaub(j,i,k)* bgewicht(j)*volufak/60.
c      betaorg(j,i,k) = betaorg(j,i,k)* bgewicht(j)*volufak/60.
c      pmseorg(j,i,k) = pmseorg(j,i,k)* bgewicht(j)*volufak/60.

c      xneis(j,i,k)=xneis(j,i,k)* bgewicht(j)*volufak/80.
c      xnstaub(j,i,k)=xnstaub(j,i,k)* bgewicht(j)*volufak/80.
c      betaorg(j,i,k) = betaorg(j,i,k)* bgewicht(j)*volufak/80.
c      pmseorg(j,i,k) = pmseorg(j,i,k)* bgewicht(j)*volufak/80.

c      xneis(j,i,k)=xneis(j,i,k)* bgewicht(j)*volufak/100.
c      xnstaub(j,i,k)=xnstaub(j,i,k)* bgewicht(j)*volufak/100.
c      betaorg(j,i,k) = betaorg(j,i,k)* bgewicht(j)*volufak/100.
c      pmseorg(j,i,k) = pmseorg(j,i,k)* bgewicht(j)*volufak/100.

      xneis(j,i,k)=xneis(j,i,k)* bgewicht(j)*volufak/100.
      xnstaub(j,i,k)=xnstaub(j,i,k)* bgewicht(j)*volufak/100.
      betaorg(j,i,k) = betaorg(j,i,k)* bgewicht(j)*volufak/100.
      pmseorg(j,i,k) = pmseorg(j,i,k)* bgewicht(j)*volufak/100.
      ext_coef126(j,i,k) = ext_coef126(j,i,k)* bgewicht(j)*volufak/100.
      ext_coef200(j,i,k) = ext_coef200(j,i,k)* bgewicht(j)*volufak/100.
      ext_coef532(j,i,k) = ext_coef532(j,i,k)* bgewicht(j)*volufak/100.
      ext_coef1000(j,i,k)= ext_coef1000(j,i,k)*bgewicht(j)*volufak/100.
      ext_coef3000(j,i,k)= ext_coef3000(j,i,k)*bgewicht(j)*volufak/100.

      densice(j,i,k) = densice(j,i,k)* bgewicht(j)*volufak/100.   ! jetzt echte g/km**3

      if (reff_nenner(j,i,k).ge.1.)
     &  eisradiuseff(j,i,k) = reff_zaehler(j,i,k) / reff_nenner(j,i,k)

      enddo
      enddo
      enddo
c_______________________________________________________________________
CC      print*,'vor open von:  '!,charname2
CC 6651 open (66,file=charname2, access='direct', recl=ngridp)
cc     print*,'nach open'
CC      write (66, rec=record_number*varnum+1, err=6652) betaorg
CC      write (66, rec=record_number*varnum+2, err=6652) densice !g/km**3
CC      write (66, rec=record_number*varnum+3, err=6652) eisradius
CC      write (66, rec=record_number*varnum+4, err=6652) xneis
CC      write (66, rec=record_number*varnum+5, err=6652) xnstaub
CC      write (66, rec=record_number*varnum+6, err=6652) eisradiuseff
CC      write (66, rec=record_number*varnum+7, err=6652) hm
CC      write (66, rec=record_number*varnum+8, err=6652) um
CC      write (66, rec=record_number*varnum+9, err=6652) vmx
CC      write (66, rec=record_number*varnum+10, err=6652) tm
CC      write (66, rec=record_number*varnum+11, err=6652) zm1
CC      close (66)
CC      record_number = record_number +1
cc     print*,'nach close'
CC      goto 6653
CC 6652   close(66,iostat=ios)
CC        print*,'Error writing beta0 '
CC        call sleep (10)
CC        goto 6651
CC 6653 print*,' ---  beendet schreiben beta0 ----  ' !!!,d_std
c_______________________________________________________________________
      do j = 1, nbneu
	do i= 1, igitneu
	  do k = 1, kgitneu		
	    aq(k,i,j,tu) = betaorg(j,i,k)
	    bq(k,i,j,tu) = densice(j,i,k)
	    cq(k,i,j,tu) = eisradius(j,i,k)
	    dq(k,i,j,tu) = xneis(j,i,k)
	    eq(k,i,j,tu) = xnstaub(j,i,k)
	    fq(k,i,j,tu) = eisradiuseff(j,i,k)
	    gq(k,i,j,tu) = hm (j,i,k)
	    hq(k,i,j,tu) = um(j,i,k)
	    iq(k,i,j,tu) = vm(j,i,k)
	    jq(k,i,j,tu) = tm(j,i,k)
	    kq(k,i,j,tu) = zm1(j,i,k)
	    lq(k,i,j,tu) = ext_coef126(j,i,k)
	    mq(k,i,j,tu) = ext_coef200(j,i,k)
	    nq(k,i,j,tu) = ext_coef532(j,i,k)
	    oq(k,i,j,tu) = ext_coef1000(j,i,k)
	    pq(k,i,j,tu) = ext_coef3000(j,i,k)
          enddo
	enddo
      enddo
		

      END SUBROUTINE schreibe_beta0
c_______________________________________________________________________
c      
c_______________________________________________________________________
c
c
      SUBROUTINE set_monate (ntage_monat,max_doy,i_jahr)
      integer ntage_monat(12),max_doy,i_jahr,mschalt

      ntage_monat(1) = 31
      ntage_monat(2) = 28
      ntage_monat(3) = 31
      ntage_monat(4) = 30
      ntage_monat(5) = 31
      ntage_monat(6) = 30
      ntage_monat(7) = 31
      ntage_monat(8) = 31
      ntage_monat(9) = 30
      ntage_monat(10) = 31
      ntage_monat(11) = 30
      ntage_monat(12) = 31
      mschalt = mod(i_jahr,4)   ! falls modulo = 0  --> schaltjahr
c                                z.B 2000 ist schaltjahr

      if (mschalt.eq.0) ntage_monat(2) = 29  ! schaltjahr
      max_doy=365
      if(ntage_monat(2).eq.29) max_doy=366

      END SUBROUTINE set_monate


      SUBROUTINE sub_transwalcek (dx,dy,dz,um,vm,wm,hm)
c-----------------------------------------------------------------------
      real*4 q0i(0:igitneu+1),qni(igitneu), ui(0:igitneu),
     &        den0i(igitneu),den1i(igitneu),dd0i(0:igitneu)
      real*4 q0j(0:nbneu+1),qnj(nbneu), uj(0:nbneu),
     &        den0j(nbneu),den1j(nbneu),dd0j(0:nbneu)
      real*4 q0k(0:kgitneu-lunten),qnk(kgitneu-lunten-1),
     &        uk(0:kgitneu-lunten-1),
     &        den0k(kgitneu-lunten-1),den1k(kgitneu-lunten-1),
     &        dd0k(0:kgitneu-lunten-1)
      real*4 hm(nbneu,igitneu,kgitneu)
      real*4 um(nbneu,igitneu,kgitneu)
      real*4 vm(nbneu1,igitneu,kgitneu)
      real*4 wm(nbneu,igitneu,kgitneu1)
      real*4 dx(nbneu),d0,deltax,dy,dz
      integer i,j,k,ip1,im1,iperiod,kk

c Parallelisierung von sub_transwalcek hat das Programm gebremst
c-----------------------------------------------------------------------
      d0 = 1.
c-----------------------------------------------------------------------
!       den0 = -99999.
!       den1 = -99999.
!       dd0  = -99999.
!       q0   = -99999.
!       qn   = -99999.
c-----------------------------------------------------------------------
c-----------  ost - west transport
c-----------------------------------------------------------------------
!       print *,'OST-WEST'
c   dd0 (0:ix) :=: d0    dichte am Zellenrand    
! walcek empfiehlt upstream, nicht mittel
      dd0i  =  d0
      iperiod = 1
      do 100 k=22,kgitneu-2
       do 101 j=2,nbneu

        deltax = dx(j)    !  wichtig

!         ui = -99999.
!         do i=1,igitneu
!          q0i(i) = hm(j,i,k)
!         enddo
!         q0i(0) = q0i(igitneu)
!         q0i(igitneu+1)=q0i(1)
! 
!         do i=1,igitneu
!          ui(i-1) = um(j,i,k)
!         enddo
!         ui(igitneu)=ui(0)

        do i=1,igitneu
         q0i(i) = hm(j,i,k)
         ui(i-1) = um(j,i,k)
        enddo
        ui(igitneu)=ui(0)
        q0i(0) = q0i(igitneu)
        q0i(igitneu+1)=q0i(1)

        do i=1,igitneu
         den0i(i) = d0  !  dichte im Inneren der Zelle
         den1i(i) = den0i(i) - dttrans/dx(j)*(d0*ui(i) - d0*ui(i-1)) ! eq. 5
        enddo
c


c-----------------------------------------------------------------------

        call sub_advec1d(igitneu,igitneu1,iperiod,
     &                 deltax,ui,den0i,den1i,dd0i,q0i,qni)

        do i=1,igitneu
         hm(j,i,k) = qni(i)
        enddo
c-----------------------------------------------------------------------
 101   continue
 100  continue
c-----------------------------------------------------------------------
c-----------  nord - sued transport
c-----------------------------------------------------------------------
!       print *,'NORD-SUED'
      deltax = dy    !  wichtig
      iperiod = 0
c   dd0 (0:ix) :=: d0    dichte am Zellenrand    
! walcek empfiehlt upstream, nicht mittel
      dd0j  =  d0
      do 200 k=22,kgitneu-2
       do 201 i=1,igitneu
        ip1=i+1
        im1=i
        if (i.eq.igitneu) ip1=1

!         q0j=-9999.
!         uj = -9999.

!         do j=1,nbneu
!          q0j(j) = hm(j,i,k)
!         enddo
!         q0j(0) = q0j(1)
!         q0j(nbneu+1)=q0j(nbneu)
!         do j=1,nbneu1
!          uj(j-1) = vm(j,i,k)
!         enddo

        do j=1,nbneu
         q0j(j) = hm(j,i,k)
         uj(j-1) = vm(j,i,k)
        enddo
        uj(nbneu)=vm(nbneu+1,i,k)
        q0j(0) = q0j(1)
        q0j(nbneu+1)=q0j(nbneu)

        do j=1,nbneu
c      urechts=um(j,ip1,k)
c      ulinks=um(j,im1,k)
         den0j(j) = d0 -dttrans/dx(j)*(d0*um(j,ip1,k)-d0*um(j,im1,k))
         den1j(j) = den0j(j) - dttrans/dy*(d0*uj(j) - d0*uj(j-1))
        enddo
c

c-----------------------------------------------------------------------
        call sub_advec1d(nbneu,nbneu1,iperiod,
     &                 deltax,uj,den0j,den1j,dd0j,q0j,qnj)

        do j=2,nbneu
         hm(j,i,k) = qnj(j)
        enddo
c-----------------------------------------------------------------------
 201   continue
 200  continue
c-----------------------------------------------------------------------
c-----------   vertikal - transport
c-----------------------------------------------------------------------
      klevels = kgitneu-lunten-1     
!  von einschliesslich 22 - 161 = 140 werte
      klevels1 = klevels+1
!       print *,'Vertikal'
      deltax = dz    !  wichtig
      iperiod = 0
c   dd0 (0:ix) :=: d0    dichte am Zellenrand    
! walcek empfiehlt upstream, nicht mittel
      dd0k  =  d0
      do 300 i=1,igitneu
       ip1=i+1
       im1=i
       if (i.eq.igitneu) ip1=1
       do 301 j=2,nbneu

!         q0k=-9999.
!         uk = -9999.

!         do kk=0,klevels1
!          k=kk+lunten-1
!          q0k(kk) = hm(j,i,k)
!         enddo
! 
!         do k=0,klevels
!          uk(k) = wm(j,i,k+lunten)
!         enddo

        do kk=1,klevels1
         k=kk+lunten-1
         q0k(kk) = hm(j,i,k)
         uk(kk-1) = wm(j,i,k)
        enddo
        q0k(0)=hm(j,i,lunten-1)


        do kk=1,klevels
         k = kk+(lunten-1)

         den0k(kk) = d0 - dttrans/dx(j)*(d0*um(j,ip1,k)- d0*um(j,im1,k))
     &             - dttrans/dy   *(d0*vm(j+1,i,k) - d0*vm(j,i,k))
         den1k(kk) = den0k(kk) -dttrans/dz*(d0*uk(kk) -d0*uk(kk-1))
        enddo
c

c-----------------------------------------------------------------------
        call sub_advec1d(klevels,klevels1,iperiod,
     &                 deltax,uk,den0k,den1k,dd0k,q0k,qnk)

        do kk=1,klevels
         k=kk+(lunten-1)
         hm(j,i,k) = qnk(kk)
        enddo
c-----------------------------------------------------------------------
 301   continue
 300  continue
      END SUBROUTINE sub_transwalcek

      SUBROUTINE sub_advec1d(ix,ix1,iperiod,
     &                            dxx,u,den0,den1,dd0,q0,qn)
c     This subroutine calculates change in mixing ratio (q0) during time
c     step dt due to advection along a grid idim in length. Mixing ratios
c     from host code (c) are loaded into q0 array, which is updated to qn.
c     Velocities (u) and fluxes (flux) are specified at cell faces, having
c     dimensions 0:idim. u, q0, qn, dxx and flux indices defined here:
c     Fluid densities flowing across each face (dd0), & beginning and end
c     of each dimension step (den0, den1) are defined in HOST CODE
c i grid->   |   1   |   2   |  i-1  |   i   | .. .. |   idim   | <- host grid
c u-array-> u(0)    u(1)    u(2)   u(i-1)   u(i)             u(idim)
c flux ->  flux(0)               flux(i-1) flux(i)          flux(idim)
c dd0  ->   dd0(0)                dd0(i-1)  ddo(i)           dd0(idim)
c c-array->  |  c(1) |  c(2) | c(i-1)|  c(i) | .. .. |  c(idim) | mixing ratio
c dxx-arry-> |  dx1  |  dx2  | dxi-1 |  dxi  | .. .. |  dxidim  |
c density->  |  dd1  |  dd2  | ddi-1 |  ddi  | .. .. |  ddidim  |
c                   q0 defined along 0 - idim+1 cells:
c    |       |   qn  |   qn  |   qn  |   qn  |       |     qn   |        |
c    |   q0--|---q0--|---q0--|---q0--|---q0--| .. .. |----q0----|---q0   |
c    |   0   |   1   |   2   |  i-1  |   i   |       |   idim   |  idim+1|
c   lower BC |                <---  q0 grid  --->               | upper BC
c           Boundary conditions are stored in q0 cells 0 & idim+1
c  im orginal ist dxx(ix)   sogar variabel !!!
c-----------------------------------------------------------------------
c
c q0 (0: ix+1) input mixing ratio with two additional boundary values
c dd0(0: ix) input fluid density flowing between each grid cell
c            remains constant for all dimensions at the initial
c            fluid density of the 1st dimension of a 2-3 D calculation
c            one can use upstream density here (dd0(i) = rho0(i) if u>0
c            or dd0(i)=rho0(i+1) if u<0)  where rho0 is the initial fluid 
c            density at the beginning of the 1st dimensional advection step
c            of a 2 or 3 D advection calculation done one step at a time 
c-----------------------------------------------------------------------
      implicit none
      integer idim,i,ix,ix1,iperiod

      real*4 q0(0:ix1), qn(ix), u(0:ix), den0(ix),
     &        den1(ix),dd0(0:ix)
      real*4 dxx

      real*4 flux(0:ix), vcmax(ix), vcmin(ix)
      real*4 ck1,ck2,cf,cf1,zr0,x1,x1n
      logical imxmn(0:ix1)

      zr0=0.
      imxmn(0) = .FALSE. 
      idim = ix
      imxmn(idim +1) = .FALSE.

C     Identify local max and min, specify mixing ratio limits at new time

c-----------------------------------------------------------------------
      do 5 i=1,idim

        imxmn(i) = (q0(i-1)-q0(i))*(q0(i)-q0(i+1)).le.zr0
c  the statement above is about 10 % faster on many cpu's
c  on sun's the staement below is faster
c        imxmn(i) = q0(i).ge.max(q0(i-1), q0(i+1)) .or.      !=true if local
c     &             q0(i).le.min(q0(i-1), q0(i+1))           !   extreme at i

      ck1 = q0(i)
      ck2 = q0(i)
      if (u(i).lt.zr0) ck1 = q0(i+1)
      if (u(i-1).gt.zr0) ck2 = q0(i-1)
      vcmax(i) = max( q0(i), ck1, ck2 )                     ! Eq.7
      vcmin(i) = min( q0(i), ck1, ck2 )                     ! Eq.7
    5 continue
c-----------------------------------------------------------------------
c     Update mixing ratios and limit Fluxes going UP where u>0
c-----------------------------------------------------------------------

      if (u(0).ge.zr0) flux(0) = q0(0)*u(0)*dttrans*dd0(0)
              !upstream at boundary
c--- Diese Parallelisierung funktioniert nicht wegen flux(i-1)
C-- !$OMP PARALLEL DO PRIVATE (x1,x1n,cf,cf1)
      do 10 i=1, idim
      if (u(i).lt.zr0) goto 10
      if (u(i-1).lt.zr0) then            ! outflow-only cell: upstream
        flux(i) = q0(i)*u(i)*dttrans*dd0(i)
      else
        x1 = u(i)*dttrans/dxx                                 !Courant number
        x1n = (1.-x1)*(q0(i+1) -q0(i-1))/4.
        cf = q0(i) +x1n                                ! Eq.4a

        if (imxmn(i-1))cf=q0(i)+max(1.5,1.2+0.6*x1)*x1n ! Eq.10b
        if (imxmn(i+1)) cf=q0(i)+(1.75-0.45*x1)*x1n       ! Eq.10a
c       cf = q0(i) +5.*x1n   ! uncomment this line for 'full sharp' ?

        cf1 = min(max(cf,min(q0(i),q0(i+1))), max(q0(i), q0(i+1)) )
c-----------------------------------------------------------------------
c
c  DEBUG   DEBUG    DEBUG    DEBUG   DEBUG   DEBUG
c      cf1 = cf   ! this statement ignores monotonic limitations.
c                  you should uncomment this line and run the calculation
c                  with constant initial mixing ratios everywhere. If you have
c                  properly implemented this subroutine constant mixing ratios
c                  should be maintained !!
c  DEBUG   DEBUG    DEBUG    DEBUG   DEBUG   DEBUG
c-----------------------------------------------------------------------

        qn(i) = max( vcmin(i), min( vcmax(i),             ! Eq.3 & 8
     &      (q0(i)*den0(i) -x1*cf1*dd0(i) +flux(i-1)/dxx)/den1(i) ))

c-----------------------------------------------------------------------
c  DEBUG   DEBUG    DEBUG    DEBUG   DEBUG   DEBUG
c      qn(i)=(q0(i)*den0(i) -x1*cf1*dd0(i) +flux(i-1)/dxx)/den1(i)
c  DEBUG   DEBUG    DEBUG    DEBUG   DEBUG   DEBUG
c-----------------------------------------------------------------------

        flux(i) = dxx*(q0(i)*den0(i) -qn(i)*den1(i)) +flux(i-1)
      endif                                                     ! Eq.9a
   10 continue
C-- !$OMP END PARALLEL DO
c  if periodic boundary conditions are assumed, it is necessary to recalculate
c  the updated mixing ratio at cell 1 if there inflow to that cell from the
c  boundary between IDIM and 1. Set IPERIOD = 1 for periodic boundary conditions      
      if (iperiod.eq.1) then
         if (u(idim-1).ge.zr0.and.u(idim).ge.zr0)
     &      qn(1)=(q0(1)*den0(1)-flux(1)/dxx+flux(idim)/dxx)/den1(1)
         endif
   
c-----------------------------------------------------------------------
c     Update mixing ratios and limit Fluxes going DOWN where u<0
c-----------------------------------------------------------------------

      if (u(idim).lt.zr0) 
     &   flux(idim)=q0(idim+1)*u(idim)*dttrans*dd0(idim)
c--- Diese Parallelisierung funktioniert nicht wegen flux(i-1)
c--!$OMP PARALLEL DO PRIVATE (x1,x1n,cf,cf1)
      do 20 i=idim, 1, -1                                      ! Eq.13
      if (u(i-1).ge.zr0) then                  ! Inflow-only cell
      if (u(i).lt.zr0) qn(i) = max( vcmin(i), min( vcmax(i),
     &   (q0(i)*den0(i) -flux(i)/dxx +flux(i-1)/dxx)/den1(i)))
      else
        x1 = dttrans*abs(u(i-1))/dxx             ! Courant number
        x1n = (1.d0 -x1)*(q0(i-1) -q0(i+1))/4.d0
        cf = q0(i) +x1n                                          ! Eq.4b
        if (imxmn(i+1)) cf=q0(i)+max(1.5d0,1.2d0+.6d0*x1)*x1n  ! Eq.10b
        if (imxmn(i-1)) cf=q0(i)+(1.75d0-.45d0*x1)*x1n         ! eq.10a
        cf1 = min(max(cf,min(q0(i), q0(i-1)) ), max(q0(i), q0(i-1)) )

c-----------------------------------------------------------------------
c  DEBUG   DEBUG    DEBUG    DEBUG   DEBUG   DEBUG
c      cf1 = cf   ! this statement ignores monotonic limitations.
c  DEBUG   DEBUG    DEBUG    DEBUG   DEBUG   DEBUG
c-----------------------------------------------------------------------

        if (u(i).ge.zr0) cf1 = q0(i)        ! outflow-only cell upstream
        qn(i) = max( vcmin(i), min( vcmax(i),                 ! Eq.3 & 8
     &      (q0(i)*den0(i)-flux(i)/dxx-x1*cf1*dd0(i-1))/den1(i)))

c-----------------------------------------------------------------------
c  DEBUG   DEBUG    DEBUG    DEBUG   DEBUG   DEBUG
c      qn(i)=(q0(i)*den0(i)-flux(i)/dxx-x1*cf1*dd0(i-1))/den1(i)
c  DEBUG   DEBUG    DEBUG    DEBUG   DEBUG   DEBUG
c-----------------------------------------------------------------------

        flux(i-1) = dxx*(qn(i)*den1(i) -q0(i)*den0(i)) +flux(i)
      endif                                                      ! Eq.9b
   20 continue
c-- !$OMP END PARALLEL DO

c-----------------------------------------------------------------------

c  if periodic boundary conditions are assumed, it is necessary to recalculate
      if (iperiod.eq.1) then
        if (u(1).lt.zr0.and.u(idim).lt.zr0)
     &   qn(idim)=(q0(idim)*den0(idim)-flux(0)/dxx+flux(idim-1)/dxx)
     &               /den1(idim)
        endif
      
      END SUBROUTINE sub_advec1d
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
      SUBROUTINE sub_lesedyn (ihemi,ijahr_fixdyn,
     &                        ijahr,imonat,itag,istd,um,tm,dm,vm,wm,zm)
c-----------------------------------------------------------------------
!       USE nrtype

      integer ::ijahr,ijahr_fixdyn,imonat,itag,istd,ihemi,i,j,k

c  neu
      real*4 dm_out_1run_fix(31,kgitneu,nbneu)
      real*4 um_out_1run_fix(31,kgitneu,nbneu)
      real*4 vm_out_1run_fix(31,kgitneu,nbneu)
      real*4 wm_out_1run_fix(31,kgitneu,nbneu)
      real*4 tm_out_1run_fix(31,kgitneu,nbneu)
      real*4 zm_out_1run_fix(31,kgitneu,nbneu)

      real*4 dm_out_6run(31,kgitneu,nbneu)
      real*4 um_out_6run(31,kgitneu,nbneu)
      real*4 vm_out_6run(31,kgitneu,nbneu)
      real*4 wm_out_6run(31,kgitneu,nbneu)
      real*4 tm_out_6run(31,kgitneu,nbneu)
      real*4 zm_out_6run(31,kgitneu,nbneu)

      real*4 dm_out_6run_fix(31,kgitneu,nbneu)
      real*4 um_out_6run_fix(31,kgitneu,nbneu)
      real*4 vm_out_6run_fix(31,kgitneu,nbneu)
      real*4 wm_out_6run_fix(31,kgitneu,nbneu)
      real*4 tm_out_6run_fix(31,kgitneu,nbneu)
      real*4 zm_out_6run_fix(31,kgitneu,nbneu)
c----



      real*4 dichte_lima(nbneu,igitneu,kgitneu)
      real*4 u_lima(nbneu,igitneu,kgitneu)
      real*4 v_lima(nbneu,igitneu,kgitneu)
      real*4 w_lima(nbneu,igitneu,kgitneu)
      real*4 t_lima(nbneu,igitneu,kgitneu)
      real*4 v_shift(nbneu,igitneu,kgitneu)
      real*4 w_shift(nbneu,igitneu,kgitneu)
      real*4 dichte_shift(nbneu,igitneu,kgitneu)

      real*4, dimension(nbneu,igitneu,kgitneu)  :: um,tm,dm,zm
      real*4, dimension(nbneu1,igitneu,kgitneu) :: vm
      real*4, dimension(nbneu,igitneu,kgitneu1) :: wm
      
      real*4 rstern
      character*100 charname
      character*21 input_pfad_2
      character*26 input_pfad_2a
      character*27 input_pfad_2b
      character*14 input_pfad_3,input_pfad_4

      integer ios

c-----------------------------------------------------------------------
      

      if (ihemi.eq.1) input_pfad_3='.hem_s_icegrid'
      if (ihemi.eq.2) input_pfad_3='.hem_n_icegrid'

      if (ihemi.eq.1) input_pfad_4='.zpr_s_icegrid'
      if (ihemi.eq.2) input_pfad_4='.zpr_n_icegrid'

      write(input_pfad_2,125) ijahr_fixdyn,imonat,
     &                        ijahr_fixdyn,imonat,itag,istd
  125 format(I4.4,'/',I2.2,'/',I4.4,'-',I2.2,'-',I2.2,'.',I2.2)
c-----------------------------------------------------------------------
      charname= LIMA_PFAD // input_pfad_2 // input_pfad_4
      open (98,
     &         status='old',form='unformatted',file=charname)
      read(98,err=4004) zm  ! druckhoehen auf den geom eishohen !
      close (98)
c-----------------------------------------------------------------------
      charname= LIMA_PFAD // input_pfad_2 // input_pfad_3
c-----------------------------------------------------------------------

c     tag 8 : t_min = 147.3138 satt0 = 0.734578   -> 22.mai  : start
c     tag 9 : t_min = 145.1851 satt0 = 1.06721    -> 23.mai
      print*,charname
 4001 continue
      open (98,iostat=ios,err=4001,
     &         status='old',form='unformatted',file=charname)
      if (ios.ne.0)  then
        print*,'alarm open fuer lesen backgr : ',charname,ios
        call sleep (2)
        goto 4001
      endif
      read(98,err=4004) dichte_lima ! dichte kg/m**3  /frueher pdruck
      read(98,err=4004) u_lima
      read(98,err=4004) v_lima
      read(98,err=4004) w_lima
      read(98,err=4004) t_lima
      close (98)
      goto 4003

 4004 close(98)
      print*,'Error reading from backgr: ',charname
      call sleep (2)
      goto 4001

 4003 continue

c----------------------------------------------------------------------
c   neu
      write(input_pfad_2a,126) ijahr_fixdyn,imonat,
     &                            ijahr_fixdyn,imonat
  126 format(I4.4,'/',I2.2,'/uvwtd_mean_',I4.4,'_',I2.2)
      charname= LIMA_PFAD // input_pfad_2a //   '.bin'
      print*,charname
      open (22,status='old',form='unformatted',file=charname)
      read (22) dm_out_1run_fix
      read (22) um_out_1run_fix
      read (22) vm_out_1run_fix
      read (22) wm_out_1run_fix
      read (22) tm_out_1run_fix
      close (22)

      write(input_pfad_2a,126) ijahr_fixdyn,imonat,
     &                            ijahr_fixdyn,imonat
      charname= LIMA_PFAD // input_pfad_2a //   '.binzz'
      print*,charname
      open (22,status='old',form='unformatted',file=charname)
      read (22) dm_out_6run_fix
      read (22) um_out_6run_fix
      read (22) vm_out_6run_fix
      read (22) wm_out_6run_fix
      read (22) tm_out_6run_fix
      close (22)

      write(input_pfad_2a,126) ijahr,imonat,ijahr,imonat
      charname= LIMA_PFAD // input_pfad_2a //   '.binzz'
      print*,charname
      open (22,status='old',form='unformatted',file=charname)
      read (22) dm_out_6run
      read (22) um_out_6run
      read (22) vm_out_6run
      read (22) wm_out_6run
      read (22) tm_out_6run
      close (22)

c  jetzt druckhoehen  in km


      write(input_pfad_2b,127) ijahr_fixdyn,imonat,
     &                            ijahr_fixdyn,imonat
  127 format(I4.4,'/',I2.2,'/zpress_mean_',I4.4,'_',I2.2)
      charname= LIMA_PFAD // input_pfad_2b //   '.bin_zpr'
      print*,'zm_out_1run_fix: ',charname
      open (22,status='old',form='unformatted',file=charname)
      read (22) zm_out_1run_fix
      close (22)

      write(input_pfad_2b,127) ijahr_fixdyn,imonat,
     &                            ijahr_fixdyn,imonat
      charname= LIMA_PFAD // input_pfad_2b //   '.bin_zprzz'
      print*,'zm_out_6run_fix: ',charname
      open (22,status='old',form='unformatted',file=charname)
      read (22) zm_out_6run_fix
      close (22)

      write(input_pfad_2b,127) ijahr,imonat,ijahr,imonat
      charname= LIMA_PFAD // input_pfad_2b //   '.bin_zprzz'
      print*,'zm_out_6run: ',charname
      open (22,status='old',form='unformatted',file=charname)
      read (22) zm_out_6run
      close (22)

c---------
c  testausdruck
      mschreib_test = 0
      if (mschreib_test .eq. 1) then
      i=1
      j=33
      print*,'itag: ',itag
      print*,' k,t_lima, 1run_fix, 6run, 6run_fix:'
      do k=1,kgitneu
      print*,k,t_lima (j,i,k),tm_out_1run_fix(itag,k,j),
     &         tm_out_6run(itag,k,j),tm_out_6run_fix(itag,k,j)
      enddo
      print*,'t: alle 120 laengen bei k=100'
      k = 100
      sum1 = 0.
      do i=1,igitneu
         sum1 = sum1 + t_lima (j,i,k)
         print*,i,t_lima (j,i,k)
      enddo
         print*,'mean: ',sum1/float(igitneu)

      i=1
      j=33
      print*,' k,zm, 1run_fix, 6run, 6run_fix:'
      do k=1,kgitneu
      print*,k,zm(j,i,k),zm_out_1run_fix(itag,k,j),
     &         zm_out_6run(itag,k,j),zm_out_6run_fix(itag,k,j)
      enddo
      print*,'zm: alle 120 laengen bei k=100'
      k = 100
      sum1 = 0.
      do i=1,igitneu
         sum1 = sum1 + zm(j,i,k)
         print*,i,zm(j,i,k)
      enddo
         print*,'mean: ',sum1/float(igitneu)
      endif
c--------

      do k=1,kgitneu
      do i=1,igitneu
      do j=1,nbneu

c      u,v,w immer gleich des jahres ijahr_fixdyn (=1976)
c      aber mit 6run_fix mittel 'out_binzz'
       u_lima(j,i,k) = u_lima(j,i,k) - um_out_1run_fix(itag,k,j)  
     &                               + um_out_6run_fix(itag,k,j)
       v_lima(j,i,k) = v_lima(j,i,k) - vm_out_1run_fix(itag,k,j)  
     &                               + vm_out_6run_fix(itag,k,j)
       w_lima(j,i,k) = w_lima(j,i,k) - wm_out_1run_fix(itag,k,j)  
     &                               + wm_out_6run_fix(itag,k,j)

c  dichte,zm, T laufen mit 6run !!
       dichte_lima (j,i,k) = 
     &         dichte_lima (j,i,k)   - dm_out_1run_fix(itag,k,j) 
     &                               + dm_out_6run    (itag,k,j)
       t_lima(j,i,k) = t_lima(j,i,k) - tm_out_1run_fix(itag,k,j) 
     &                               + tm_out_6run    (itag,k,j)
       zm    (j,i,k) = zm    (j,i,k) - zm_out_1run_fix(itag,k,j)
     &                               + zm_out_6run    (itag,k,j) 
      enddo
      enddo
      enddo
c
c   ende neu !!!
c----------------------------------------------------------------------


c-----------------------------------------------------------------------
c
c  neu:  1981 wurde einmal eine Tem unter 100 K erreicht
c        deshalb setzte min Temp auf 102 K
c-----------------------------------------------------------------------
      do k=1,kgitneu
        do i=1,igitneu
          do j=1,nbneu
          if (t_lima(j,i,k).lt.102.) t_lima(j,i,k) = 102.
      enddo
      enddo
      enddo
c-----------------------------------------------------------------------
      if(ihemi.eq.1) v_lima = -1.*v_lima
      um = u_lima
C Hier parallelisieren. cpu-nr abfragen, dann schleife     
      do k=1,kgitneu
       do j=1,nbneu
        do i=1,igitneu-1
         w_shift(j,i,k)     =(w_lima(j,i,k)+w_lima(j,i+1,k))*0.5
         tm(j,i,k)          =(t_lima(j,i,k)+t_lima(j,i+1,k))*0.5
         v_shift(j,i,k)     =(v_lima(j,i,k)+v_lima(j,i+1,k))*0.5
         dichte_shift(j,i,k) =
     &               (dichte_lima(j,i,k)+dichte_lima(j,i+1,k))*0.5
        enddo
        w_shift(j,igitneu,k)  = (w_lima(j,igitneu,k)+w_lima(j,1,k))*0.5
        tm(j,igitneu,k)       = (t_lima(j,igitneu,k)+t_lima(j,1,k))*0.5
        v_shift(j,i,k)        = (v_lima(j,igitneu,k)+v_lima(j,1,k))*0.5
        dichte_shift(j,i,k)   =
     &                (dichte_lima(j,igitneu,k)+dichte_lima(j,1,k))*0.5
       enddo
      enddo      
      
c-----------------------------------------------------------------------
      do k=1,kgitneu
       do i=1,igitneu
        do j=2,nbneu
         vm(j,i,k) = (v_shift(j-1,i,k)+v_shift(j,i,k))*0.5
        enddo
        vm(1,i,k) = vm(2,i,k) ! wichtig
        vm(nbneu1,i,k) = 0.   ! wichtig
       enddo
      enddo      
c----*------------------------------------------------------------------
      do i=1,igitneu
       do j=1,nbneu
        do k=2,kgitneu
         wm(j,i,k) = (w_shift(j,i,k)+w_shift(j,i,k-1))*0.5
        enddo
        wm(j,i,1) = wm(j,i,2)
        wm(j,i,kgitneu1) = wm(j,i,kgitneu)
       enddo
      enddo      
c-----------------------------------------------------------------------
      rstern = 8.31432*1000.  ! N*m/(K*kmol) univers. gaskonst
      do k=1,kgitneu
       do i=1,igitneu
        do j=1,nbneu
         dm(j,i,k) = rstern * dichte_shift(j,i,k) * tm(j,i,k) / 28.96  ! druck
        enddo
       enddo
      enddo      

      END SUBROUTINE sub_lesedyn 
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
      SUBROUTINE sub_leseh2o_zprinit 
c-----------------------------------------------------------------------

!       USE nrtype
      USE common_h2o
      open (98,status='old',form='unformatted',file=H2OINIT)
c  
c     166 hoehen ab 72,0 72,2 72,4  ... 104,8 105,0 km als druckhoehen
c     für breite alomar 

      read(98) zpress_init
      read(98) zpress_init_h2o
      close (98)


      END SUBROUTINE sub_leseh2o_zprinit 
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
      SUBROUTINE sub_maxmin (um,vm,wm) 
c-----------------------------------------------------------------------
!       USE nrtype
      USE common_const
      real*4, dimension(nbneu,igitneu,kgitneu)  :: um
      real*4, dimension(nbneu1,igitneu,kgitneu) :: vm
      real*4, dimension(nbneu,igitneu,kgitneu1) :: wm
      real*4 umax,ucour_max,vmax,vcour_max,wcour_max,wmax
      integer i,j,k
c-----------------------------------------------------------------------
      umax=0.
      do k=1,kgitneu
      do i=1,igitneu
      do j=1,nbneu
         if (abs(um(j,i,k)/dx(j)).gt.umax) umax=abs(um(j,i,k)/dx(j))
         enddo
         enddo
         enddo
      ucour_max = dttrans*umax
c-----------------------------------------------------------------------

      vmax=0.
      do k=1,kgitneu
      do i=1,igitneu
      do j=1,nbneu1
         if (abs(vm(j,i,k)/dy).gt.vmax) vmax=abs(vm(j,i,k)/dy)
         enddo
         enddo
         enddo
      vcour_max = dttrans*vmax
c-----------------------------------------------------------------------

      wmax=0.
      do k=1,kgitneu1
      do i=1,igitneu
      do j=1,nbneu
c
c   neu
      wcour_max = dttrans*abs(wm(j,i,k))/dz
         if (wcour_max.gt.0.9) wm(j,i,k) =
     &                       abs(wm(j,i,k))*0.9*dz/(dttrans*wm(j,i,k))

         if (abs(wm(j,i,k)/dz).gt.wmax) wmax=abs(wm(j,i,k)/dz)
         enddo
         enddo
         enddo
      wcour_max = dttrans*wmax
c-----------------------------------------------------------------------
      print*,'ucour_max : ',ucour_max
      print*,'vcour_max : ',vcour_max
      print*,'wcour_max : ',wcour_max
      if (max(ucour_max,vcour_max,wcour_max) .gt.1.) then
         print*,'alarm courant > 1',ucour_max,vcour_max,wcour_max
c         stop
      endif
      END SUBROUTINE sub_maxmin 
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
      SUBROUTINE sub_zenit (ntg,itag1,ihemi)
c
c  ntg: zeitschritt waehrend des tages
c   vorsicht ntg <-> dttrans muss passen
c
!       USE nrtype
      USE common_const
      USE common_h2o
      implicit none

c-----------------------------------------------------------------------
      real*4 skh,za,zb,zd,zccc,zf,zg,day,dek,day16,day2,dayl,
     &       breite,ortzeit,rncom,stuwi,coschi,chi,x,
     &       y,erfy,arg,zeit,zgeo_86,sunris_86
      integer i,j,k,itag,itag1,ihemi,ntg

      data skh /7./,zgeo_86/86./
      data za/1.0606963/ ,zb/0.55643831/ ,zd/1.7245609/  ,
     *     zccc/1.0619896/,zf/0.56498823/,zg/0.06651874/
c-----------------------------------------------------------------------

       sunris_86 = acos(6373./(6373.+zgeo_86))*
     &           180./pi + 90.
c      print*,' zgeo_86, sunris_86 :',zgeo_86, sunris_86

c-----------------------------------------------------------------------
c     vorlauf = 120+15 tage   bis 15.mai
c-----------------------------------------------------------------------
      itag=itag1
      if (ihemi.eq.1) itag=mod(itag1+182,365)

c   alt
c      day =float(itag)
c          ! vorsicht gilt sowohl fuer nord als auch sued hemi
c      dek=sin(-(81.-day)*360./365.*pi/180.) *23.5 *pi/180.

      dek = oblinoa(itag)*pi/180.

      dayl=86400./dttrans     !  zeitschritte pro tag
      day16=dayl/float(igitneu)
      day2=dayl/2.
      rncom=float(ntg)
      zeit=amod(rncom,dayl)

      do i=1,igitneu
      do j=1,nbneu

      breite = (37.5 +  float(j-1)) *pi/180.
      ortzeit=zeit+float(i-1)*day16
      ortzeit=amod(ortzeit,dayl)
      stuwi=ortzeit/dayl*2.*pi
      coschi=-cos(breite)*cos(dek)*cos(stuwi)+sin(breite)*sin(dek)
      if (coschi.gt.1.) coschi=1.
      if (coschi.lt.-1.) coschi=-1.
      chi=acos(coschi)*180./pi
c
      secchi(j,i)=0.
      if (chi.le.75.) then

             secchi(j,i)= 1./coschi
      else if (chi.gt.75.and.chi.lt.90.)  then
             x=(6373.+zgeo_86)/skh
             y=sqrt(.5*x)*abs(cos(chi*pi/180.))
             erfy=(za+zb*y)/(zccc+zd*y+y*y)
             if (y.gt.8.)     erfy=zf/(zg+y)
             secchi(j,i)=sqrt(pi/2.*x)*erfy
      else if (chi.gt.sunris_86)         then
           secchi(j,i)=0.
      else
              x=(6373.+zgeo_86)/skh
              y=sqrt(.5*x)*abs(cos(chi*pi/180.))
              erfy=(za+zb*y)/(zccc+zd*y+y*y)
              if (y.gt.8.)     erfy=zf/(zg+y)
              arg=x*(1.-sin (chi*pi/180.))
              secchi(j,i)=sqrt(2.*pi*x)*(sqrt(sin(chi*pi/180.))
     *                    *exp(arg)-0.5*erfy)
      endif
c
      enddo
      enddo

      return
      END SUBROUTINE sub_zenit
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
      SUBROUTINE sub_photolyse 
     &                     (mnudge,delta_t,ihemi,itag1,xlyman_obs)

      USE common_h2o
      USE common_uvwtd
      implicit none
      
      real*4 gesamt_photo_h2o,photodiss_rate_h2o,photodiss_cross_h2o
      real*4 bfak(nbneu),xnudge,faknudge,xno2
      real*4 hfak(kgitneu),so2neu(nbneu,igitneu)
      real*4  xlyman_obs(366) 
c-----------------------------------------------------------------------
      integer k,i,j,ihemi,itag1,itag,mnudge
      real*4 delta_t,xnluft,xnh2o,
     &       b1_h2o,c1_h2o,b2_h2o,c2_h2o,b3_h2o,c3_h2o,
     &       d_disso_o2_ev,d_disso_o2_erg,fak_energie,
     &       epslyman,hplank,clicht,phitop_ly,xlamda_ly,
     &       exply1,exply2,exply3,reduct_h2o,wert

c      data photodiss_cross__h2o/1.4e-17/   ! 1.4e-17 cm**2 fuer Solomon u. Brasser buch S. 227
      data photodiss_cross_h2o/1.53e-17/   ! 1.53e-17 cm**2 aus Chabrillat u. Kockarts, 1997, GRL
c=======================================================================
c     lyman alpha constanten  u. sol max - min faktoren

      itag=itag1
      if (ihemi.eq.1) itag=mod(itag1+182,365)

      phitop_ly = xlyman_obs(itag) * 1.e11   ! anzahl photonen  ! cgs - system (cm^-2 * sec^-1)
c      phitop_ly = phitop_ly * ellipsefak(itag)

c  neu
c      phitop_ly = 4.5e11   ! anzahl photonen  ! cgs - system (cm^-2 * sec^-1)

c      phitop_ly = 6.5e11   ! anzahl photonen  ! cgs - system (cm^-2 * sec^-1)
c      phitop_ly = 3.5e11   ! anzahl photonen  ! cgs - system (cm^-2 * sec^-1)

      hplank = 6.6262e-34   ! (J * s) = 4.135667e-15 eV*s
      clicht = 2.998e8      ! m/s
c
c  c    nach hamonnia (Lean et al.)strobel 1978 Formel 4 gibt Photon - Disso energie an
c
c (1 erg = 1.e-7 J)
c  1 eV = 1.60218e-12 erg
c z.B Lyman_alpha

      xlamda_ly = 1215. ! wellenlaenge des photons in Angstroem
      d_disso_o2_ev = 5.12 ! disso von o2 in eV
      d_disso_o2_erg = d_disso_o2_ev * 1.60218e-12 ! disso von o2 in erg

      fak_energie = 1.e7*hplank * clicht / (xlamda_ly*1.e-10) ! erg
      epslyman = (fak_energie - d_disso_o2_erg) / fak_energie

c      print*,' lyman fak_energie: ',fak_energie,' epslyman: ',epslyman
c-----------------------
c  photo dissoziations rate fuer o2 durch Ly_alpha
c-----------------------
c
c    phitop_ly = solar irradiance at the top of atmos at 121.5 nm (ly)
c    reduct_o2 = reduction factor (Chabrillat and Kockarts (GRL, 1997)
c                ist funktion der Hoehe (eq. 8)  einheit: cm^2
c      d1_o2 = 6.007e-21     ! 6.0073e-21
c      e1_o2 = 8.217e-21     ! 8.21666e-21
c      d2_o2 = 4.28569e-21
c      e2_o2 = 1.633e-20     ! 1.63296e-20
c      d3_o2 = 1.2806e-20    !1.28059e-20
c      e3_o2 = 4.85e-17      !4.85121e-17

c    reduct_h2o = reduction factor (Chabrillat and Kockarts (GRL, 1997)
c                ist funktion der Hoehe (eq. 8)  einheit: cm^2
      b1_h2o = 0.68431 
      c1_h2o = 8.2214e-21
      b2_h2o = 0.229841
      c2_h2o = 1.77556e-20
      b3_h2o = 0.0865412
      c3_h2o = 8.22112e-21
c=======================================================================

      xnudge = 90./(1.*86400.)
      bfak=0.1

      bfak(1) = 1.5 
      bfak(2) = 1.5
      bfak(3) = 1.5 
      bfak(4) = 1.5
      bfak(5) = 1.5 
      bfak(6) = 1.4
      bfak(7) = 1.3
      bfak(8) = 1.2
      bfak(9) = 1.1
      bfak(10) = 1.0
      bfak(11) = 1.0
      bfak(12) = 1.0
      bfak(13) = 0.9
      bfak(14) = 0.8
      bfak(15) = 0.7
      bfak(16) = 0.6
      bfak(17) = 0.5
      bfak(18) = 0.4
      bfak(19) = 0.3
      bfak(20) = 0.2
           !  vertikaler transport von einschliesslich
           ! 17 - 160 = 144 werte
      do k=1,5
       hfak(k)=40.
      enddo
      hfak(6)=40.
      hfak(7)=40.
      hfak(8)=40.
      hfak(9)=40.
      hfak(10)=40.
      hfak(11)=38.
      hfak(12)=36.
      hfak(13)=34.
      hfak(14)=32.
      hfak(15)=30.
      hfak(16)=28.
      hfak(17)=26.
      hfak(18)=24.
      hfak(19)=22.
      hfak(20)=20.
      hfak(21)=18.
      hfak(22)=16.
      hfak(23)=14.
      hfak(24)=12.
      hfak(25)=10.
      hfak(26)= 8.
      hfak(27)= 6.
      hfak(28)= 4.
      hfak(29)= 2.

      do k=30,kgitneu
       hfak(k)=1.
      enddo

      do k=1,29
       hfak(kgitneu+1-k)=hfak(k)*0.5
      enddo

      so2neu = 0.  ! wert = saule o2 cm**-2
      do i=1,igitneu
      do j=13,nbneu
         xnluft= 1.e-6*dm(j,i,kgitneu)/(xkbolz*tm(j,i,kgitneu))
         xno2= 0.2 * xnluft
c     anzahl o2 molekuele:   xnluft*0.2   ! also 20 prozent nach Sonnemann !
         so2neu(j,i) = xno2*4.5*1.e5 ! 4.5 km (94 km)skalenhoehe in cm, 
c                                    !  muesste sauele sein, mit gerd abgesprochen   
      enddo
      enddo

      do k=kgitneu,3,-1
       do i=1,igitneu
        do j=1,nbneu
       
c  xnluft: anzahl molekuele pro cm**3
c  xnh2o : anzahl h2o-molekuele pro cm**3

         xnluft= 1.e-6*dm(j,i,k)/(xkbolz*tm(j,i,k))
         xnh2o = hm(j,i,k)*1.e-6*xnluft
         
         xno2 = 0.2 * xnluft
         so2neu(j,i) = so2neu(j,i) + xno2*1.e4  ! * 100 m in cm 
c
c__________  nur wenn sonne scheint: secchi 1 - sehr gross
        if (secchi(j,i).gt.0.1)  then  ! nur wenn sonne: secchi 1 - sehr gross
         wert = so2neu(j,i) * secchi(j,i)  ! wert= saule o2 cm**2

         exply1 = -c1_h2o*wert
           if (exply1.lt. -70.) exply1 = -70.
         exply2 = -c2_h2o*wert
           if (exply2.lt. -70.) exply2 = -70.
         exply3 = -c3_h2o*wert
           if (exply3.lt. -70.) exply3 = -70.
      reduct_h2o = b1_h2o * exp(exply1) +
     &             b2_h2o * exp(exply2) +
     &             b3_h2o * exp(exply3) 

         photodiss_rate_h2o = phitop_ly * photodiss_cross_h2o 
     &                                  * reduct_h2o 
         gesamt_photo_h2o = photodiss_rate_h2o*xnh2o*delta_t
c         print*,j,i,k,xnh2o,gesamt_photo_h2o,secchi(j,i)
         xnh2o=xnh2o - gesamt_photo_h2o
         hm(j,i,k)=xnh2o*1.e6/xnluft
        endif
        enddo
       enddo
      enddo

      if (mnudge.eq.1) then

      do k=1,kgitneu
       do j=1,nbneu
c        faknudge=0.5*amin1 (4.,bfak(j)*hfak(k))    ! nochmal halbieren
c        faknudge= 1.* amin1 (4.,2.*bfak(j)*hfak(k))    ! neu
c        faknudge= 1.* amin1 (4.,3.*bfak(j)*hfak(k))    ! neu
        faknudge= 1.* amin1 (4.,5.*bfak(j)*hfak(k))    ! neu2
        do i=1,igitneu
         hm(j,i,k)=hm(j,i,k) + xnudge*faknudge*(hminit(j,k)-hm(j,i,k))
        enddo
       enddo
      enddo

c      print*,' in photo'
c      do k=23,kgitneu,10
c       write(6,1000) k,(hminit(j,k),j=1,53,5),k
c      enddo
c 1000 format(i4,2x,11(f6.3,1x),2x,i3)

      endif

      return
      END SUBROUTINE sub_photolyse 
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
      SUBROUTINE sub_ausdruck(feld)

      real*4 feld(nbneu,igitneu,kgitneu)
      integer*4 i,j,k
c-----------------------------------------------------------------------

      i=1
      do k=1,kgitneu,2
       write(6,1000) k,(feld(j,i,k),j=1,53,5),k
      enddo
 1000 format(i4,2x,11(f6.3,1x),2x,i3)
      return
      END SUBROUTINE sub_ausdruck
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
      SUBROUTINE sub_diffu_h2o (hm,dz,turbkzz)
      implicit none
!       integer lunten,kgitneu,igitneu,nbneu
!       parameter (lunten=22,kgitneu=163,igitneu=120,nbneu=53)

      real*4 hm(nbneu,igitneu,kgitneu)
      real*4 turbkzz(kgitneu),a(kgitneu),b(kgitneu),c(kgitneu),
     &          r(kgitneu),u(kgitneu)

      real*4 dz,alpha
      integer i,j,k

c
c alle 45 sec zeitschritt !!!
c
c  crank-nicholson
C Hier parallelisieren. parallel do auf i mit private ...
      do i=1,igitneu
       do j=3,nbneu
        do k=1,kgitneu
         alpha = turbkzz(k)*dttrans_2/(dz*dz) ! alt
         a(k)= - alpha
         b(k)= 1.+2.*alpha
         c(k)= - alpha
         r(k)= hm(j,i,k)
        enddo
        a(1)=0.
        c(kgitneu)=0.

        call tridiag(a,b,c,r,u,kgitneu)

        do k=lunten,kgitneu-2
         hm(j,i,k)=u(k)
        enddo
 
       enddo
      enddo
      return
      END SUBROUTINE sub_diffu_h2o 
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
      SUBROUTINE tridiag (a,b,c,r,u,n)

      integer j,n
      real*4 a(n),b(n),c(n),r(n),u(n)
      real*4 gam(NMAX),bet

      if (b(1).eq.0.) pause 'tridag: rewrite equtios'
      bet = b(1)
      u(1)=r(1)/bet

      do j=2,n
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j)*gam(j)
         if (bet.eq.0.) pause 'tridag failed'
         u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
      do j=n-1,1,-1
         u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
      return
      end
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
      SUBROUTINE sub_init_global(i_jahr,xlyman_obs,scale_ch4,unix_time)

!       USE nrtype
      USE common_const
      USE common_h2o
      implicit none
c-----------------------------------------------------------------------
      integer*4, parameter :: kgitneu2d=286
      integer*4 :: i,ii,j,k,kk,k2,i_jahr,idim,mschalt, x1, x2
      real*4    :: xwert,zuntenkm,dzkm,faktkzz,wert,a1,a2,b1,b2,
     &             ch4_akt,ch4_2008,h2oanteil_2008,h2oanteil_akt,
     &             fr_alpha,basis_h2o,h2o_true,scale_ch4,unix_time

      real*4, dimension(kgitneu2d) :: turbkzz2d,wturb2d
      real*4  xlyman_obs(366)
      character*44    char_lyman
      character*4    char_year
      character*100  c_dummy

c      zkm = 0.1 km ist 77.8 km also k=1
c      bis         k=kgitneu = 163

      do j=1,nbneu
       xwert = 37.5 +  float(j-1)
       dx(j)=2.*pi/float(igitneu)*rearth*cos(xwert*pi/180.)
      enddo
      dzkm=0.1
      zuntenkm=77.8
      do k=1,kgitneu
         zgeo(k)=zuntenkm+float(k-1)*dzkm
c         print*,k,zgeo(k)
      enddo
c-----------------------------------------------------------------------
      open (98,status='old',form='unformatted',
     *     file=CTURBKZZ)
      read (98) turbkzz2d,wturb2d
      k2=0
      print*,' kzz - daten !!!!!!!!!!!!'
      do k=2,kgitneu2d,2
       k2=k2+1
       turbkzz(k2) = turbkzz2d(k)*0.25
c      print*,k2,turbkzz(k2)
      enddo
      kk=0
      do k=144,163
       kk=kk+1
       faktkzz=1.+float(kk)/20.
       turbkzz(k)=turbkzz(143)/faktkzz
c      print*,k,turbkzz(k)
      enddo
      close (98)
      print*,'luebkzza.dat gelesen !!'
c
c   neu
c      k2 = 23+30        !  werte  = Kzz bei 83 km
c      k2 = 23+40        !  werte  = Kzz bei 84 km
      k2 = 23+50        !  werte  = Kzz bei 85 km
      do k=1,k2
         turbkzz(k)=turbkzz(k2+1)
      enddo
c-----------------------------------------------------------------------

      open (98,status='old',form='unformatted', file=CBACKSCATTER532)
      read (98) backradius, crosssec532, ext532, crosssec126, ext126,
     &         crosssec200, ext200, crosssec1000, ext1000,
     &         crosssec3000, ext3000
      close (98)
      print*,'backscatter532.dat gelesen !!'
c-----------------------------------------------------------------------
      open (98,status='old',form='unformatted',
     *     file=CELLIPSE_EARTH)
      do i=1,365
      read(98) ii,ellipsefak(i),oblinoa(i),true_solar_time(i)
      print*,ii,ellipsefak(i),oblinoa(i),true_solar_time(i)
      enddo
      close (98)
c-----------------------------------------------------------------------


c ersetzen der if-abragen (jahreszahlen) DR
c konvertiert die int jahreszahl in einen string
      write(char_year,'(i4)')i_jahr

      char_lyman= '/home/oa034/klymanalpha/'//char_year//
     & '.lyman_alpha.txt'
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!       if (i_jahr.eq.2010)
!      &    char_lyman= '/home/oa015/klymanalpha/2010.lyman_alpha.bin'


      print*,char_lyman
ccc      open (98,status='old',form='unformatted',
ccc     *     file=char_lyman)
         open (98,status='old',form='formatted', file=char_lyman)

      idim = 365
      mschalt = mod(i_jahr,4)   ! falls modulo = 0  --> schaltjahr
      if (i_jahr.eq.1900) mschalt=1  ! ausnahme 1900 u. 2100
      if (mschalt.eq.0) idim = 366

      do i=1,idim
         read(98, 777) x1, wert, x2
ccc         read(98) wert
         xlyman_obs(i) = wert !   wichtig, faktor steht in photolyse* 1.e11
c         xlyman_obs(i) = 4.75  ! solar mean
         print*,i_jahr,i,xlyman_obs(i)
      enddo
      close (98)
c-----------------------------------------------------------------------
      open (98,status='old',form='formatted',file= UT1_1861_2100)
      read (98,*)  c_dummy
      print*, '-- Tabelle unix_time --'
      do i=1,240   ! 1861 - 2100  unix_time 10.05 0:0:0 GMT
         read (98,*) b1, b2

         if (float(i_jahr).eq.b1)  unix_time = b2    
      enddo
      close (98)
          print*,b1,b2, unix_time
      close (98)
c-----------------------------------------------------------------------
      open (98,status='old',form='formatted',file= CH4_1861_2135)
      read (98,*)  c_dummy
      print*, '-- Tabelle CH4 --'
      do i=1,275   ! 1861 - 2008  ch4 werte
         read (98,*) a1, a2
c  hier NEU
         if (float(i_jahr)-5.eq.a1)  ch4_akt = a2    ! shift 5 jahre: age of air
c         if (1881.-5.eq.a1)  ch4_akt = a2    ! shift 5 jahre: age of air
c ende NEU
         if (a1.eq.2008.)            ch4_2008 = a2    ! shift 5 jahre: age of air
c         print*,a1,a2
      enddo
      close (98)
         fr_alpha = 0.95                              ! fractional release factor alpha
         ch4_akt = ch4_akt / 1000.                    ! ppmv
         ch4_2008 = ch4_2008 / 1000.                  ! ppmv
         h2oanteil_2008 = fr_alpha * ch4_2008 * 2.    ! 1 methan molekuele -> 2 Wassermolek.
         h2oanteil_akt  = fr_alpha * ch4_akt * 2.    ! 1 methan molekuele -> 2 Wassermolek.
c
c   6 ppmv sind typisch bei 78 km  (z.B. AIM 2008)
         basis_h2o = 6. - h2oanteil_2008      ! soviel basis-wasser ohne ch4 beitrag
         h2o_true = basis_h2o + h2oanteil_akt
         print*,'6 ppmv sind typisch bei 78 km  (z.B. AIM 2008)'
         print*,'basis_h2o,h2o_true: ', basis_h2o, h2o_true

         scale_ch4 = 100.*h2o_true/6.        !    in %
         print*,' scale_ch4, h2o_true ?: ',scale_ch4,
     &            h2o_true,scale_ch4*0.01*6.

         print*,
     &  'ppmv -> ch4_akt, ch4_2008, h2oanteil_2008,h2oanteil_akt: ',
     &           ch4_akt, ch4_2008, h2oanteil_2008,h2oanteil_akt

      scale_ch4 = scale_ch4 * 0.01
  777 format (I7, F6.2, I3)
      return
      END SUBROUTINE sub_init_global
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
      SUBROUTINE sub_init_staub (nstart,bbreite,bsum,j_mesopause,
     &                                               j_zpress86km)
!       USE nrtype
      use ifport
      USE common_const
      USE common_eis

c
c    25 radi-klassen a 0.1 nm ab 1.0 nm
c
      implicit none
      integer, dimension(nbneu)   :: j_mesopause,j_zpress86km
      integer, dimension(10,25)   :: nhunt
      integer, dimension(5)       :: nhisto
      integer                     :: i,k1,ii,jj,kk,kk1,kk2
      integer nstart,nklasse,nsum,nzahl,nnn,nklasse_haupt,
     &        nklasse_zwi,nwert1,nwert2,nwert3,nwert4,
     &        nwert5,n,nwerta,nwertb

      real*4, dimension(25)   :: xhisto
      real*4, dimension(35)  :: bbreite,bsum   ! neu
      real*4 wert,fak,rando,xdiff,rando1,scalewert,sum,
     &       wert1,wert2,wert3,wert4,wert5,sfactor,invfactor
c-----------------------------------------------------------------------
c      data xhisto/1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
c     &            2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
c     &            3.0,3.1,3.2,3.3,3.4,3.5,3.6/

      data xhisto/0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,    !  vorher ab 0.8
     &            2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
     &            3.0,3.1,3.2,3.3/

c-------------------------------------------
c   richtige hunten vert. nach Klostermeyer

      data nhisto/864900,117000,15800,2100,200/

c-------------------------------------------
c    wichtig es nhisto an fuer gleiche Anzahl
c    in allen hoehen fuer gl. mischungsverhaeltnis
c    also nhunt aus  --> nhisto
c    d.h data nhunt wird ueberschrieben !!
c____________________________________________
c  exp. klassen initialisierung   !!!
c      goto 4711
      print*,' pi : ',pi
      print*,' breiten grenzen !!'
      fak = 100./10.36195
      bsum = 0.

      i = 1
      bsum(1)= 5.46616
      bbreite(1)=19.    !19. =  55.0 grad
c       print*,i,bbreite(i),'bsum: ',bsum(i)

      do i=2,35
       bbreite(i)=bbreite(i-1)+1.
       wert =fak* cos((54.5+float(i))*pi/180.)
       bsum(i) = bsum(i-1) + wert
c        print*,i,bbreite(i),'bsum: ',bsum(i)
      enddo
c-----------------------------------------------------------------------
c Hier werden der Sinus des Anfangsbreitenkreises (55°) und sein Kehrwert
c ausgerechnet, damit diese Kalkulation nicht in der 40-mio Schleife liegt
c-----------------------------------------------------------------------
        sfactor = sin(55.*pi/180.)
        invfactor = 1.-sfactor

c____________________________________________
      if (nstart.eq.1) goto 4711
c____________________________________________
c
c  neu    ueberschreiben von nhunt
c         zur erstellung von mischungsverhaeltnissen
c         in allen heoehen gleich
c-----------------------------------------------------------------------
      do k1=1,10
       do nklasse_haupt=1,5   ! jede Klasse wird linear in 5 segmente a 0.1 mn geteilt

        nwerta = nhisto(nklasse_haupt)
        nwertb = 1
        if (nklasse_haupt.lt.5) nwertb = nhisto(nklasse_haupt+1)

        wert1 = float(nwertb) + 5.*float(nwerta-nwertb)
        wert2 = float(nwertb) + 4.*float(nwerta-nwertb)
        wert3 = float(nwertb) + 3.*float(nwerta-nwertb)
        wert4 = float(nwertb) + 2.*float(nwerta-nwertb)
        wert5 = float(nwertb) + 1.*float(nwerta-nwertb)

        sum = wert1+wert2+wert3+wert4+wert5
        scalewert = float(nhisto(nklasse_haupt))/sum

        nwert1 = int(wert1*scalewert)
        nwert2 = int(wert2*scalewert)
        nwert3 = int(wert3*scalewert)
        nwert4 = int(wert4*scalewert)
        nwert5 = int(wert5*scalewert)
        nsum = nwert1+nwert2+nwert3+nwert4+nwert5

c   wichtig
        nwert1 = nwert1 +(nhisto(nklasse_haupt)-nsum)

        if (k1.eq.10) then
        print*,nwert1
        print*,nwert2
        print*,nwert3
        print*,nwert4
        print*,nwert5
        endif
        nsum = nwert1+nwert2+nwert3+nwert4+nwert5
        print*,nsum,nhisto(nklasse_haupt)
        print*,' '

        nhunt(k1,1+(nklasse_haupt-1)*5) = nwert1
        nhunt(k1,2+(nklasse_haupt-1)*5) = nwert2
        nhunt(k1,3+(nklasse_haupt-1)*5) = nwert3
        nhunt(k1,4+(nklasse_haupt-1)*5) = nwert4
        nhunt(k1,5+(nklasse_haupt-1)*5) = nwert5

       enddo
      enddo

      nsum=0

      do k1=1,10                       ! 10 hoehen
       do nklasse=1,25                   ! 25 klassen
        nzahl = nhunt(k1,nklasse) * 4   ! von 10 auf 40 mio
        print*,k1,nklasse,nzahl
        do nnn=1,nzahl
         nsum=nsum + 1
c        print*, nsum
         rando = rand()
         xdiff=0.1     ! klasse von delta 0.1 nm
         rinit(nsum)= xhisto(nklasse) + rando*xdiff

 3335    continue
c       print*, '***'
         rando = float(igitneu)*rand()
         xfeld(nsum)=1.+rando     ! 120 laemgen =2 mal 60
         ii=int(xfeld(nsum))
         if(ii.gt.igitneu.or.ii.lt.1) goto 3335
c       print*, '*'
c
c  nord - sued verteilung 1 : 86
c
 3336   continue
        rando = sfactor + invfactor*rand() ! erzeugt Zufallszahl sin(55°) < rand < 1
         ! Transformation der gleichförmigen Verteilung in am Pol gestreckte
        wert = asin(rando)*180./pi
        yfeld(nsum) = wert-36
        jj = int(yfeld(nsum))
c-----------------------------------------------------------------------
c Wiederholung, falls der y-Wert außerhalb der spezifizierten Schranke liegt
c-----------------------------------------------------------------------
        if (jj.gt.nbneu.or.jj.lt.1) then
                print*, sfactor, invfactor
c               print*,' rando - grenze FEHLER !!! ',rando, wert
                goto 3336
        endif

c       print*, '**'
        if (yfeld(nsum).gt.53.999) yfeld(nsum) = 53.999

 3334    continue
c       print*, '*4*'
         rando = rand()  ! 77.8 km=1, 80 km=23, 83.5 km=58
c         zfeld(nsum)=73. + float(k1-1)*5. +  5.* rando !  neu: 85 - 90 km
c         zfeld(nsum)=23.+60. +  10.* rando !  neu: 86 - 87 km

c   neu
         if (rando.gt. 0.8) then
            rando = (rando - 0.8)*5.                        !  wieder zwischen 0:1
c            kk1 = j_mesopause(jj) - (23+55+10)             ! mesop - 86.5 km
            kk1 = j_mesopause(jj) - (5+j_zpress86km(jj))    ! mesop - 86.5 km
            kk2 = 15
            if (kk1.ge.1) kk2 = kk2 + kk1 
c            zfeld(nsum)= 23.+57. +  float(kk2)* rando      !  neu: 1 km unterhalb der mesopa
            zfeld(nsum)= float(j_zpress86km(jj)) - 3. +  
     &                   float(kk2)* rando                  !  neu: 1 km unterhalb der mesopa
c            print*,'kk1,kk2,rando:', kk1,kk2,rando
         else
         rando = rando * 1./0.8  ! wieder 0:1
c         zfeld(nsum)=23. + 55. +  10.* rando !  neu: 85.5 - 86.5 km
         zfeld(nsum)= float(j_zpress86km(jj)) - 5. +  10.* rando !  neu: 85.5 - 86.5 km
         endif
c  ende neu

         kk=int(zfeld(nsum))
        if(kk.lt.1.or.kk.gt.kgitneu) goto 3334
c       print*, '*5*'
        enddo

       enddo
      enddo
      print*,' nsum = 40 mio = ',nsum

      rfeld = rinit
 4711 continue

c-----------------------------------------------------------------------
      return
      END SUBROUTINE sub_init_staub
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      SUBROUTINE sub_tracertransp
!       USE nrtype
      use ifport
      USE common_const
      USE common_tab
      USE common_eis
      USE common_uvwtd
      USE common_h2o

      implicit none
      integer*4     intel_no, intel_offset
      integer*4     i,j,k,n,itemp,nl,j0,j1,j2,i0,i1,i2,ii,jj,kk
      real*4, dimension(ntrac) :: change_n
      real*4, dimension(nbneu) :: ax
      real*4, dimension(nbneu,igitneu,kgitneu)  :: work1,change
      real*4 eps,aconst1,ay,az,rwert0,tback,tp,work1akt,volufak,
     &       drdt,rwert1,vpsatt_tp,backgr,dsum,wfall,wrando,
     &       satt_tp,dback,dreferenz

c-----------------------------------------------------------------------
c      print*,' rearth = 6366197.0 ',rearth
c      print*,' gmeso  = 9.55 ',gmeso
c      print*,' xkbolz = 1.3805e-23 ',xkbolz
c      print*,' xmh2o  = 2.99e-26 ',xmh2o  !  kg
c      print*,' xmair  = 4.845e-26 ',xmair
c      print*,' rhoice = 932.0 ',rhoice   !  kg/m**3
c-----------------------------------------------------------------------
c   xnluft  : anzahl luft molekuele   pro cm**3
c   xnh2o   : anzahl h2o-molekuele    pro cm**3
      change = 0.    ! anz wasser molekuele
      change_n = 0.    ! anz wasser molekuele
c-----------------------------------------------------------------------
c temp0 : T in 90 km auf laenge 1, k-level: 123
c druck0: p in 90 km auf laenge 1, k-level: 123

c-----------------------------------------------------------------------
      eps = 1.e-4
      aconst1=1.e-27*rhoice*4.*pi/(3.*xmh2o)
c    1/(sqrt(2.)*2.)  -->   1/4    veringert fall geschw.
c-----------------------------------------------------------------------

      ax = dttrans_2/dx
      ay = dttrans_2/dy
      az = dttrans_2/dz
c-----------------------------------------------------------------------

      work1 = 1.e-6 * hm *dm

c-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE (i,j,k,rwert0,tback,tp,rwert1,work1akt,itemp,
!$OMP+ vpsatt_tp,satt_tp,drdt,backgr,volufak,dback,dreferenz)

      do 1000 n=1,ntrac
       i = int(xfeld(n))
       j = int(yfeld(n))
       k = int(zfeld(n))
       rwert0 = rfeld(n)

       tback = tm(j,i,k)
       tp = tback

c     eis --> eis  waechst  oder schmilzt
c-----------------------------------------------------------------------
      if (rwert0.gt.rinit(n))   then
          work1akt= work1(j,i,k)

c
c    partikel temp wert1 von worka_tab(k) = 0.45 bis 0.1
c    die 0.1 damit auf keinen Fall, der index 0 auftritt,
c    falls staubradien bei 1.0 nm beginnen
c_________
c          tp = worka_tab(int(rwert0+0.1),k) * 2. + tp  ! 1.=200 %
          tp = worka_tab(int(rwert0+0.1),k) * 1. + tp  ! 1.=100 %

          if (tp.le.100.) tp = 100.1
          itemp = int(tp-99.)
          itemp = min0(itemp,71)
c
c    Kelvin formel
c_________
          vpsatt_tp= vpsatt_tab(itemp) *
     &               faktkelvin_tab(min0(300,int(10.*rwert0)),itemp)
          satt_tp = work1akt / vpsatt_tp    !
c-----------------------------------------------------------------------
c
c    eisradi  --> wachsen schmelzen
c_________
         drdt =  (satt_tp-1.)*vpsatt_tp*sqrt1t_tab(itemp)
       rwert1 = rwert0 + drdt*dttrans_2
       rwert1 = amax1(rinit(n),rwert1)
       rfeld(n) = rwert1
c-----------------------------------------------------------------------
c     backgr:  anzahl wasser molekuele in eisschale
       backgr = aconst1*(rwert1*rwert1*rwert1 - rwert0*rwert0*rwert0)
c       dback = dm(j,i,k)         !  alt
c       dreferenz = dm(j,i,123)   !  alt
c
c  volufak: umrechnen auf 1/cm**3
c
c       volufak = horiwicht_tab(j) *dreferenz/dback !  alt
       volufak = horiwicht_tab(j) !  * volufak_tab(k) ! immer an
c       change_n(n) =  backgr*volufak
c       change_n(n) =  backgr*volufak*0.7   ! neu : um 30 proz kleiner
c       change_n(n) =  backgr*volufak*0.7*8./10.   ! neu 
       change_n(n) =  backgr*volufak*0.7*8./10.   ! neu 
c-----------------------------------------------------------------------

      else

c     staub nukleation zu eis     ! immer unter 150 K
c-----------------------------------------------------------------------
c      if (tback.lt. 150.) then    ! alt
      if (tback.lt. 155.) then
          work1akt= work1(j,i,k)

          itemp = int(tp-99.)
          itemp = min0(itemp,71)
          vpsatt_tp = vpsatt_tab(itemp)

c
c    Kelvin formel
c_________
       vpsatt_tp= vpsatt_tp * faktkelvin_tab(int(10.*rwert0),itemp)
       satt_tp = work1akt / vpsatt_tp

c
c    eisbildung nur wenn ...
c_________
      if (satt_tp.gt.1.)   then

       drdt =  (satt_tp-1.)*vpsatt_tp*sqrt1t_tab(itemp)
       rwert1 = rwert0 + drdt*dttrans_2
       rwert1 = amax1(rinit(n),rwert1)
       rfeld(n) = rwert1
c-----------------------------------------------------------------------
c     backgr:  anzahl wasser molekuele in eisschale
       backgr = aconst1*(rwert1*rwert1*rwert1 - rwert0*rwert0*rwert0)
c         dback = dm(j,i,k)        ! alt
c         dreferenz = dm(j,i,123)  ! alt
c
c  volufak: umrechnen auf 1/cm**3
c
c       volufak = horiwicht_tab(j) *dreferenz/dback   ! alt
       volufak = horiwicht_tab(j)  ! * volufak_tab(k)    ! immer an
c       change_n(n) =  backgr*volufak
c       change_n(n) =  backgr*volufak*0.7   ! neu : um 30 proz kleiner
c       change_n(n) =  backgr*volufak*0.7*8./10.   ! neu 
       change_n(n) =  backgr*volufak*0.7*8./10.   ! neu 
         endif   ! ende eisbildung nur wenn ...
c-----------------------------------------------------------------------
      endif

      endif
 1000 continue

!$OMP END PARALLEL DO
c-----------------------------------------------------------------------

       intel_no=int(rand()*49.) + 1
       intel_offset=int(rand()*1.e6)

       if (intel_no.gt.50)  intel_no = 50
       if (intel_no.lt.1)   intel_no = 1

!$OMP PARALLEL DO PRIVATE (nl,i,j,k,itemp,wrando,wfall,ii,jj,kk)
      do 2000 n=1,ntrac
       i = int(xfeld(n))
       j = int(yfeld(n))
       k = int(zfeld(n))

c  wturb richtig !!!   w_turb = sqrt (D/a)
c  0.5 weil wrandos auf delta t = 45 sec def. waren
c                  nun          = 90 sec
c  also halbieren der zufallswege
c  wenn 1.0, dann sind nach einer Stunde nach oben bis zu 5 km erreicht, bzw.
c            nach unten mehr als 3 km

c         wrando =  0.5*(2.* rand() - 1.)* wrichtig_tab(k)

c        wrando =  0.5*(2.* ran(n*intel_no+intel_offset) - 1.)
c    &                  * wrichtig_tab(k)

c   neu
         wrando =  0.5*(2.* ran(n*intel_no+intel_offset) - 1.)
     &                  * wrichtig_tab(k)
       


c 02.04.2008 nutze 2.0*... vorher (Langzeitläufe 0.5*...)
!  sensitivitaet mit 4facher
c     wfall:  fallgeschwindigkeit [m/s]   pos. nach unten !!
c
c    /(2.*sqrt(2.)) --> 1./4. veringert fallgeschw.(klostermeyer-formel)

c      wfall = 1.e-9*(rhoice*gmeso*rwert/(4.*pdruck(i,k)))*
c     *         sqrt(pi*xkbolz*ti/xmair)

       itemp = int(tm(j,i,k)-99.)
       itemp = min0(itemp,71)
       wfall = rfeld(n)*sqrtsedi_tab(itemp)/dm(j,i,k)

       xfeld(n)=xfeld(n) + ax(j) * um(j,i,k)
       yfeld(n)=yfeld(n) + ay * vm(j,i,k)
       zfeld(n)=zfeld(n) + az * (wm(j,i,k)+wrando-wfall)

       jj=int(yfeld(n))
       kk=int(zfeld(n))
C -- Teste Feldgrenzenueberschreitung nach transport
       if (jj.gt.nbneu) then
        yfeld(n) = 2.*(nbneu.+1.) - yfeld(n)
        if (jj.gt.nbneu) then
         yfeld(n) = yfeld(n) -eps
         jj=int(yfeld(n))
!          print *,'jj left north (used eps) jj,xfeld(n)',jj,yfeld(n),n
        endif
        xfeld(n) = 1.+amod(xfeld(n)-1.+igitneu./2.,igitneu.)
       else
        if (jj.lt.1) yfeld(n) = 1.0001
       endif

       ii=int(xfeld(n))
       if (ii.gt.igitneu) then
!         print *,'ii left right ii,xfeld(n)',ii,xfeld(n),n
        xfeld(n)=xfeld(n) - igitneu.
        ii=int(xfeld(n))
        if (ii.lt.1) xfeld(n)=xfeld(n) + eps
!         print *,'xfeld neu(n),n',xfeld(n),n
       else
        if (ii.lt.1) then
         xfeld(n)=igitneu. + xfeld(n)
         ii=int(xfeld(n))
         if (ii.gt.igitneu) then
          xfeld(n)=xfeld(n)-eps
          ii=int(xfeld(n))
!           print *,'ii left out (used eps) ii,xfeld(n)',ii,xfeld(n),n
         endif
        endif
       endif
c
c  neu .ge. statt .gt.    und .le. statt .lt.
       if (kk.ge.kgitneu) then
        zfeld(n)=kgitneu.999
       else
        if (kk.le.1) zfeld(n)=1.001
       endif

! For Debug only    , ist nie aufgetreten !!! also auskommentieren
c       nl = int(zfeld(n))
c       if(nl.lt.1.or.nl.gt.kgitneu) then
c         print *,'k outof bounds i,j,k,nl',i,j,k,nl
c         print *,'n,z(n),az,wm(j,i,k),wrando,wfall',n,zfeld(n),az,
c     &         wm(j,i,k),wrando,wfall
c         print *,'rfeld(n),itemp,sqrtsedi_tab(itemp),dm(j,i,k)',
c     &            rfeld(n),itemp,sqrtsedi_tab(itemp),dm(j,i,k)
c       endif
! ----
 2000 continue
!$OMP END PARALLEL DO
C-- !$OMP FLUSH
C --- Wasserdampf aenderung ausrechnen
      do 3000 n=1,ntrac
        i = int(xfeld(n))
        j = int(yfeld(n))
        k = int(zfeld(n))
        change(j,i,k)=change(j,i,k) + change_n(n)
 3000 continue

C-- !$OMP PARALLEL DO PRIVATE(dsum)
      do k=17,kgitneu-3
      dsum=0.
       do i=1,igitneu
        dsum=dsum+change(nbneu,i,k)
       enddo
       dsum=dsum/float(igitneu)
       do i=1,igitneu
        change(nbneu,i,k)=dsum
       enddo
      enddo
C-- !$OMP END PARALLEL DO
c      a1=0.4*change(49,i,k)
c      a2=0.4*change(48,i,k)
c      change(50,i,k)= 0.2*change(50,i,k)+a1+a2
c      change(51,i,k)= 0.2*change(51,i,k)+a1+a2
c      change(52,i,k)= 0.2*change(52,i,k)+a1+a2
c      change(53,i,k)= 0.2*change(53,i,k)+a1+a2
c      enddo
c      enddo

      change = hm - change*xkbolz*tm/dm
      hm = amax1(0.01,change)

C-- !$OMP PARALLEL DO PRIVATE(dsum)
      do k=17,kgitneu-3
       dsum=0.
       do i=1,igitneu
        dsum=dsum+hm(nbneu,i,k)
       enddo
       dsum=dsum/float(igitneu)
       do i=1,igitneu
        hm(nbneu,i,k)=dsum
       enddo
      enddo
C-- !$OMP END PARALLEL DO

      work1=hm
C-- !$OMP PARALLEL DO PRIVATE(i,i0,i1,i2,j,j0,j1,j2)
      do k=17,kgitneu-3
       do i=1,igitneu
       i1=i
       i0=i-1
       if (i.eq.1) i0=igitneu
       i2=i+1
       if (i.eq.igitneu) i2=1
       j=nbneu
       j1=j
       j0=j-1
       j2=nbneu
       hm(j,i,k)=(work1(j0,i0,k)+ work1(j1,i0,k)+work1(j2,i0,k)+
     &          work1(j0,i1,k)+ work1(j1,i1,k)+work1(j2,i1,k)+
     &          work1(j0,i2,k)+ work1(j1,i2,k)+work1(j2,i2,k))/9.
       j=nbneu-1
       j1=j
       j0=j-1
       j2=j+1
       hm(j,i,k)=(work1(j0,i0,k)+ work1(j1,i0,k)+work1(j2,i0,k)+
     &          work1(j0,i1,k)+ 3.*work1(j1,i1,k)+work1(j2,i1,k)+
     &          work1(j0,i2,k)+ work1(j1,i2,k)+work1(j2,i2,k))/11.
       j=nbneu-2
       j1=j
       j0=j-1
       j2=j+1
       hm(j,i,k)=(work1(j0,i0,k)+ work1(j1,i0,k)+work1(j2,i0,k)+
     &          work1(j0,i1,k)+ 9.*work1(j1,i1,k)+work1(j2,i1,k)+
     &          work1(j0,i2,k)+ work1(j1,i2,k)+work1(j2,i2,k))/17.
       j=nbneu-3
       j1=j
       j0=j-1
       j2=j+1
       hm(j,i,k)=(work1(j0,i0,k)+ work1(j1,i0,k)+work1(j2,i0,k)+
     &          work1(j0,i1,k)+ 27.*work1(j1,i1,k)+work1(j2,i1,k)+
     &          work1(j0,i2,k)+ work1(j1,i2,k)+work1(j2,i2,k))/35.
       j=nbneu-4
       j1=j
       j0=j-1
       j2=j+1
       hm(j,i,k)=(work1(j0,i0,k)+ work1(j1,i0,k)+work1(j2,i0,k)+
     &          work1(j0,i1,k)+ 81.*work1(j1,i1,k)+work1(j2,i1,k)+
     &          work1(j0,i2,k)+ work1(j1,i2,k)+work1(j2,i2,k))/89.
       enddo
      enddo
C-- !$OMP END PARALLEL DO

      return
      END SUBROUTINE sub_tracertransp
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      SUBROUTINE sub_tabelle

!       USE nrtype
      USE common_const
      USE common_tab

      implicit none
      real*4  breite,rwert0,sigma,tp,wert,alt,v1,v2,v3,scale
      integer j,k,iradi,itemp
c-----------------------------------------------------------------------
      print*,' rearth = 6366197.0 ',rearth
      print*,' gmeso  = 9.55 ',gmeso
      print*,' xkbolz = 1.3805e-23 ',xkbolz
      print*,' xmh2o  = 2.99e-26 ',xmh2o
      print*,' xmair  = 4.845e-26 ',xmair
      print*,' rhoice = 932.0 ',rhoice
c-----------------------------------------------------------------------
c
c      wachstum
c
c      drdt =  1.e9*0.83*(satt_tp-1.)*vpsatt_tp*
c     &             sqrt(xmh2o/(2.*pi*xkbolz*tp)) / rhoice

c
c  kelvin formel:     faktkelvin_tab  0.1 - 30 nm = (90 radien, 71 temp)
c
c              sigma=0.12/(1.+2.*0.15/rwert0)
c              fakt=exp(1.e9*2.0*xmh2o*sigma/
c     &			  (rhoice*xkbolz*tp*rwert0))
c-----------------------------------------------------------------------
c   60 als skalierung ist so gewaehlt, dass dass fuer initstaub
c  hunten werte folgen  !!!!!
c 1.e12 folgen aus hm in ppmv u. xnluft mit 1.e-6
c   float(igitneu)  pro laengen segment
c!$OMP MASTER
      do j=1,nbneu
         breite = 37.5+float(j-1) 
                                  
c
c  horiwicht_tab u. volufak gibt kopplung an hintergrundwasserdampf an !
c
c         horiwicht_tab(j) =  (1./2.)*   ! 1./2. wegen 20 -> 40 Mio
         horiwicht_tab(j) =  (1.)*   ! nehme 1 dies steht auch in beta routinen
     &     1.e12*float(igitneu)/   !  die 1.e12 siehe 5 zeilen drueber
     &     (cos(breite*pi/180.)*60.)
      enddo
c-----------------------------------------------------------------------
c       vpsatt = exp(28.548 - 6077.4/ti)    ! alte gadsden u. Schroeder 1989
c       vpsatt = 10.**(14.88 - 3059./ti)      ! neue Mauersberger uu. Krankowsky, 2003
c       vpsatt = exp(9.550426 - 5723.265/tm(j,i,k) +   ! Murph a. Koop 2005 aus Rapp
c     &           3.53068*alog(tm(j,i,k)) - 0.00728332*tm(j,i,k))
c-----------------------------------------------------------------------
      print*,'tp,iradi,rwert0,faktkelvin_tab(ir,it)'
      vpsatt_tab = 0.
      sqrt1t_tab     = 0. 
      do itemp=1,TAB_TMAX       ! 100 - 170 K
        tp = float(itemp)+99.
c        vpsatt_tab(itemp) = sngl(dexp(dble(28.548 - 6077.4/tp)))
c        vpsatt_tab(itemp) = 10.**(dble(14.88 - 3059./tp))
       vpsatt_tab(itemp) = exp(9.550426 - 5723.265/tp +   ! Murph a. Koop 2005 aus Rapp
     &           3.53068*alog(tp) - 0.00728332*tp)
c        v1 = sngl(dexp(dble(28.548 - 6077.4/tp)))
c        v2 = 10.**(dble(14.88 - 3059./tp))
c        v3 = exp(9.550426 - 5723.265/tp +   ! Murph a. Koop 2005 aus Rapp
c     &           3.53068*alog(tp) - 0.00728332*tp)
c      print*,itemp,tp,v1,v2,v3

        sqrt1t_tab(itemp) = sngl(1.e9*0.83*
     &     dsqrt(dble(xmh2o/(2.*pi*xkbolz*tp)))/rhoice)
c      print*,itemp,vpsatt_tab(itemp)

c    1/(sqrt(2.)*2.)  -->   1/4    veringert fall geschw.     
c      sqrtsedi_tab(itemp) = sngl(1.e-9*rhoice*gmeso*
c     &     dsqrt(dble(pi*xkbolz*tp/xmair)/4.)) ! mittlerer wert fuer 140 K

      sqrtsedi_tab(itemp) = sngl(1.e-9*rhoice*gmeso*
     &     dsqrt(dble(pi*xkbolz*tp/xmair)/(sqrt(2.)*2.))) ! mittlerer wert fuer 140 K

       do iradi=1,TAB_RMAX     !  0.1 - 30 nm radius in 0.1 schritten 
        rwert0= float(iradi)*0.1

c	 sigma=0.12/(1.+2.*0.15/rwert0)
        sigma=sngl(dble(0.141 - 1.5e-4*tp)/dble(1.+2.*0.15/rwert0))    ! formel 9
        faktkelvin_tab(iradi,itemp) = sngl(
     &     dexp(dble(1.e9*2.0*xmh2o*sigma/
     &         (rhoice*xkbolz*tp*rwert0))))
       enddo
      enddo

c    hier wird der 30 nm wert abgezogen um stetigkeit zu gantieren zwioschen
c    wachsen und schmelzen !!   wichtig

      print*,' '
      print*,' neu faktkelvin_tab : '
      print*,' '
      do itemp=1,TAB_TMAX       ! 100 - 170 K
       do iradi=1,TAB_RMAX     !  0.1 - 30 nm radius in 0.1 schritten 
        wert = faktkelvin_tab(TAB_RMAX,itemp) - 1.
        alt = faktkelvin_tab(iradi,itemp)
        faktkelvin_tab(iradi,itemp) = alt - wert
c        if (iradi.eq.10)
c     & print*,iradi/10,itemp,faktkelvin_tab(iradi,itemp)
       enddo
c       print*,' kelvin: ',alt,itemp,faktkelvin_tab(300,itemp)
      enddo
c-----------------------------------------------------------------------
      print*,'k,zgeo(k),turbkzz(k),volufak_tab(k)'
      do k=1,kgitneu
       wrichtig_tab(k) = sqrt(turbkzz(k)/7.5)
       volufak_tab(k) = 1./exp((88.-zgeo(k))/4.0)
       print*,k,zgeo(k),turbkzz(k),volufak_tab(k)
      enddo
c-----------------------------------------------------------------------
      print*,'k,zgeo(k),worka_tab(k)'

c  neu  :  zwischen 1- 10 nm gibt es keine Partikeltemp 
c          zwischen 11 - 29 nm steigt linear gewicht von 0 auf 1

      do k=1,kgitneu
       do iradi=1,TAB_RPTP  !  (1 - 200 nm)
       scale = 1.
       if (iradi.le.10) scale = 0.

       if (iradi.eq.11) scale = 0.05
       if (iradi.eq.12) scale = 0.1
       if (iradi.eq.13) scale = 0.15
       if (iradi.eq.14) scale = 0.2
       if (iradi.eq.15) scale = 0.25
       if (iradi.eq.16) scale = 0.3
       if (iradi.eq.17) scale = 0.35
       if (iradi.eq.18) scale = 0.4
       if (iradi.eq.19) scale = 0.45
       if (iradi.eq.20) scale = 0.5
       if (iradi.eq.21) scale = 0.55
       if (iradi.eq.22) scale = 0.6
       if (iradi.eq.23) scale = 0.65
       if (iradi.eq.24) scale = 0.7
       if (iradi.eq.25) scale = 0.75
       if (iradi.eq.26) scale = 0.8
       if (iradi.eq.27) scale = 0.85
       if (iradi.eq.28) scale = 0.9
       if (iradi.eq.29) scale = 0.95
c
c   falls alter zustand gewuenscht, setze ->
c       scale = 1.  ! altes feld  ( hier mehr kleines Eis, wegen freeze dry)
       
       worka_tab(iradi,k) = scale*float(iradi) * 
     &    sngl(0.01*beta05 * t82km * dexp(dble((zgeo(k)-82.)/4.5))) 
       enddo
c      print*,zgeo(k),worka_tab(10,k),worka_tab(20,k),worka_tab(30,k)
      enddo

c!$OMP END MASTER
      return
      END SUBROUTINE sub_tabelle
c
c
c-----------------------------------------------------------------------
      SUBROUTINE sub_dmget (ihemi,ijahr,imonat,itag)
c-----------------------------------------------------------------------
c     LIMA_PFAD    './LC33/LIMA-ICE/backgr/'
     
      integer :: ihemi,ijahr,imonat,itag
      
      character*100 charname
      character*18 input_pfad_2
      character*18 input_pfad_3
      character*18 input_pfad_4

      if (ihemi.eq.1) input_pfad_3='.*.hem_s_icegrid &'
      if (ihemi.eq.2) input_pfad_3='.*.hem_n_icegrid &'
      if (ihemi.eq.1) input_pfad_4='.*.zpr_s_icegrid &'
      if (ihemi.eq.2) input_pfad_4='.*.zpr_n_icegrid &'

      write(input_pfad_2,125) ijahr,imonat,ijahr,imonat,itag
  125 format(I4.4,'/',I2.2,'/',I4.4,'-',I2.2,'-',I2.2)

      charname = LIMA_PFAD // input_pfad_2 // input_pfad_3
      print*,'dmget von ', charname
      res = SYSTEM('dmget '//charname)    !  res = result

      charname = LIMA_PFAD // input_pfad_2 // input_pfad_4
      print*,'dmget von ', charname
      res = SYSTEM('dmget '//charname)    !  res = result

      return
      END SUBROUTINE sub_dmget
c
c
c-----------------------------------------------------------------------
      SUBROUTINE sub_dmput (ihemi,ijahr,imonat,itag)
c-----------------------------------------------------------------------
c     LIMA_PFAD    './LC33/LIMA-ICE/backgr/'
     
      integer :: ihemi,ijahr,imonat,itag
      
      character*100 charname
      character*18 input_pfad_2
      character*18 input_pfad_3
      character*18 input_pfad_4

      if (ihemi.eq.1) input_pfad_3='.*.hem_s_icegrid &'
      if (ihemi.eq.2) input_pfad_3='.*.hem_n_icegrid &'
      if (ihemi.eq.1) input_pfad_4='.*.zpr_s_icegrid &'
      if (ihemi.eq.2) input_pfad_4='.*.zpr_n_icegrid &'

      write(input_pfad_2,125) ijahr,imonat,ijahr,imonat,itag
  125 format(I4.4,'/',I2.2,'/',I4.4,'-',I2.2,'-',I2.2)

      charname = LIMA_PFAD // input_pfad_2 // input_pfad_3
      print*,'dmput von ', charname
      res = SYSTEM('dmput -r '//charname)    !  res = result

      charname = LIMA_PFAD // input_pfad_2 // input_pfad_4
      print*,'dmput von ', charname
      res = SYSTEM('dmput -r '//charname)    !  res = result

      return
      END SUBROUTINE sub_dmput
c
c
c
c-----------------------------------------------------------------------
      SUBROUTINE sub_h2oinit_zpr_to_zgeo (zm)
c-----------------------------------------------------------------------
      USE common_h2o

      real*4 zm(nbneu,igitneu,kgitneu)
      real*4 zpressakt(kgitneu)
      real*4 fneu(kgitneu)
c-----------------------------------------------------------------------


c      do k=1,kgitzpr
c      print*,k,zpress_init(k),zpress_init_h2o(k)
c      enddo

      do j=1,nbneu
      do i=1,igitneu
      do k=1,kgitneu
         zpressakt(k)=zm(j,i,k)
      enddo

      call sub_intpol(kgitzpr,zpress_init,zpress_init_h2o,
     &                kgitneu,zpressakt,fneu)

      do k=1,kgitneu
         hminit3d(j,i,k)=fneu(k)
      enddo

      enddo
      enddo


      END SUBROUTINE sub_h2oinit_zpr_to_zgeo
c-----------------------------------------------------------------------
c
c
c
      subroutine sub_intpol(n,x,y,ni,xi,yi)
c=======================================================================
      implicit none
c-----------------------------------------------------------------------
      integer   n,ni,i
      real*4    x(n),y(n),xi(ni),yi(ni),
     &          p(380),q(380),y1(380),
     &          a(380),b(380),c(380),d(380)
      real*4 ps
c-----------------------------------------------------------------------
      ps=4.1
      do 10 i=1,n
            p(i)=ps
            q(i)=ps
  10        continue
c
c
      y1(1)=(y(2)-y(1))/(x(2)-x(1))
      y1(n)=(y(n)-y(n-1))/(x(n)-x(n-1))
      call sub_raspl1(n,x,y,p,q,y1,a,b,c,d)
c
      do 20 i=1,ni
      call sub_yspline(n,x,xi(i),a,b,c,d,yi(i),p,q)
 20   continue
c
      return
      end
c
c
c
      subroutine sub_yspline(n,x,xwert,a,b,c,d,ywert,p,q)
      real*4    x(380),a(380),b(380),c(380),d(380),
     *          p(380),q(380)

      i=1
 11   if(xwert.le.x(i+1) )         goto 10
      i=i+1
      if(i.eq.n)                  goto 8
                                  goto 11
 8    i=n-1
 10   t=(xwert-x(i))/(x(i+1)-x(i))
      u=1.-t
c
      ywert=a(i)*u+b(i)*t+c(i)*u*u*u/(p(i)*t+1.)+d(i)*t*t*t/(
     &q(i)*u+1.)
      return
      end
c
c
c
      subroutine sub_raspl1(n,x,y,p,q,y1,a,b,c,d)
c
      real*4    x(n),y(n),p(380),q(380),y1(380),
     &          a(380),b(380),
     &          c(380),d(380)
      n1=n-1
      c(1)=0.
      d(1)=0.
c
      do 2 k=1,n1
           j2=k+1
           pp=p(k)
           qq=q(k)
           pp2=pp*(pp+3.)+3.
           qq2=qq*(qq+3.)+3.
           p22=2.+pp
           q22=2.+qq
           a(k)=x(j2)-x(k)
           h=1./a(k)
           b(k)=1./(p22*q22-1.)
           h2=h*b(k)
           r2=h*h2*(y(j2)-y(k))
           if(k.eq.1)                               goto 1
           hq=h1*qq1
           hp=h2*pp2
           z=1./(hq*(p21-c(j1))+hp*q22)
           c(k)=z*hp
           h=r1*qq1*(1.+p21)+r2*pp2*(1.+q22)
           if(k.eq.2)   h=h-hq*y1(1)
           if(k.eq.n1)  h=h-hp*y1(n)
           d(k)=z*(h-hq*d(j1))
  1        j1=k
           p21=p22
           qq1=qq2
           h1=h2
           r1=r2
  2   continue
      y1(n1)=d(n1)
      if(n1.le.2)                                   goto 4
      n2=n1-1
      do 3 j1=2,n2
           k=n-j1
           y1(k)=d(k)-c(k)*y1(k+1)
  3   continue
  4   do 5 k=1,n1
           j2=k+1
           h=b(k)*(y(j2)-y(k))
           z=b(k)*a(k)
           p2=2.+p(k)
           q2=2.+q(k)
           c(k)= (1.+q2)*h-z*(y1(j2)+q2*y1(k))
           d(k)=-(1.+p2)*h+z*(p2*y1(j2)+y1(k))
           a(k)=y(k)-c(k)
           b(k)=y(j2)-d(k)
  5   continue
      return
      end

