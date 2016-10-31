  subroutine op_peachin

 !! opk_peach reads in the opacity tables for C and N from the 
 !! file peach_tables
 !! file is found in lte_codes/share/osdata/data/peach_tables
 !! Updated version of Sterne subroutine PEACH

 !! reads in total opacity for KIEL: TAB...
 !! and only free-free for OP: TFF...

    use op_peach
 
    implicit none

 !! unit number for file
      integer :: NP=10
      integer :: i,j

 !! open peach_tables and read into arrays
 
      open (unit=NP,file='PEACH_TABLES',status='old',action='read')

 !! read in the Carbon opacities

      read(NP,*) ((TAB1C1(i,j),i=1,6),j = 1,37)
      read(NP,*) ((TAB2C1(i,j),i=1,6),j = 1,37)
      read(NP,*) ((TAB3C1(i,j),i=1,6),j = 1,37)
      read(NP,*) ((TAB4C1(i,j),i=1,6),j = 1,37)
      read(NP,*) ((TAB5C1(i,j),i=1,6),j = 1,37)
      read(NP,*) ((TAB1C2(i,j),i=1,6),j = 1,33)
      read(NP,*) ((TAB2C2(i,j),i=1,6),j = 1,33)
      read(NP,*) ((TAB3C2(i,j),i=1,6),j = 1,33)
      read(NP,*) ((TAB4C2(i,j),i=1,6),j = 1,33)
      read(NP,*) ((TAB5C2(i,j),i=1,6),j = 1,33)
      read(NP,*) ((TAB1C3(i,j),i=1,6),j = 1,40)
      read(NP,*) ((TAB2C3(i,j),i=1,6),j = 1,40)
      read(NP,*) ((TAB3C3(i,j),i=1,6),j = 1,40)
      read(NP,*) ((TAB4C3(i,j),i=1,6),j = 1,40)
      read(NP,*) ((TAB5C3(i,j),i=1,6),j = 1,40)
      read(NP,*) ((TAB1C4(i,j),i=1,6),j = 1,48)
      read(NP,*) ((TAB2C4(i,j),i=1,6),j = 1,48)
      read(NP,*) ((TAB3C4(i,j),i=1,6),j = 1,48)
      read(NP,*) ((TAB4C4(i,j),i=1,6),j = 1,48)
      read(NP,*) ((TAB5C4(i,j),i=1,6),j = 1,48)

 !! read in the Nitrogen opacities

      read(NP,*) ((TAB1N1(i,j),i=1,6),j = 1,45)
      read(NP,*) ((TAB2N1(i,j),i=1,6),j = 1,45)
      read(NP,*) ((TAB3N1(i,j),i=1,6),j = 1,45)
      read(NP,*) ((TAB4N1(i,j),i=1,6),j = 1,45)
      read(NP,*) ((TAB5N1(i,j),i=1,6),j = 1,45)
      read(NP,*) ((TAB1N2(i,j),i=1,6),j = 1,40)
      read(NP,*) ((TAB2N2(i,j),i=1,6),j = 1,40)
      read(NP,*) ((TAB3N2(i,j),i=1,6),j = 1,40)
      read(NP,*) ((TAB4N2(i,j),i=1,6),j = 1,40)
      read(NP,*) ((TAB1N3(i,j),i=1,6),j = 1,27)
      read(NP,*) ((TAB2N3(i,j),i=1,6),j = 1,27)
      read(NP,*) ((TAB3N3(i,j),i=1,6),j = 1,27)

      close(NP)

  end subroutine op_peachin
