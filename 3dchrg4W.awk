#!/bin/bash -l
#==========================================================================
# 3dchrg4W.avg ### calculate avg charge in bins of volume for a moiety  ###
#              ### W is the water atoms                                 ###
# parameters
# $1         number of horizontal slices
# $2         file name, contains the charges for each atom in the .gro file
# $3         name of the output file
# $4         membrane channel to analyze
#   options  0 - top membrane, 1 - bottom membrane
# $5         number of extra water molecules in simulation
#   #! NOTE  there should have always been 11480 H2O molecules in each
#            simulation. Somewhere, somehow, another variant with 11484 H2O
#            molecules crept in. We have to know which we're dealing with.
# $6         number of ions in alpha bath
# $7         number of ions in beta bath
#__________________________________________________________________________

ARG1=${1:-'50'}
ARG2=${2:-'chrgs.xvg'}
ARG3=${3:-'chrg4W.xvg'}
ARG4=${4:-'0'}
ARG5=${5:-'0'}
ARG6=${6:-'11'}
ARG7=${7:-'11'}

if [ -s "${ARG3}" ]; then s ../jobs/roll. ${ARG3}; fi

declare -a boxdims
boxdims=($(cat ../boxdims.xvg))

awk	-v PREC="quad" \
	-v zd=${boxdims[2]} -v xd=${boxdims[0]} -v yd=${boxdims[1]} -v ch=${ARG4} \
	-v pags=${ARG1} -v cfyl=${ARG2} -v wa=${ARG5} -v io0=${ARG6} -v io1=${ARG7} \
    '
BEGIN	{ frms=2001
	  dlns=0; drcs=0
          dza=zd/pags
          cols=int(xd/dza)+1; rows=int(yd/dza)+1
          dxa=xd/cols; dya=yd/rows
          xs=int(2/dxa); ys=int(2/dya); zs=int(2/dza)
          iz=cols*rows; iy=rows
          ibxs=cols*rows*pags
          for ( i=1; i<=ibxs; i++ ) {
            chrg[i]=0;
        } }

BEGINFILES
        /^[^\n#@]/ \
        { if ( FILENAME != "energy.xvg" ) {
            split( $0, chrgs )
          }
          if ( FILENAME == "energy.xvg" ) {
            if ( drcs%10 == 0 ) {
              dlns++
              dims[dlns][""]
              split( $0, dims[dlns])
            }
            drcs++
        } }
ENDFILES

END	\
        { fr0=386+22080
	  fr1=fr0+386+22080+34440+wa*3+io0*2
	  for ( irecs=1; irecs<=frms; irecs++ ) {
            getline < "../coord.xvg" > 0
            FIRST_POSN=match( $0, /[\n@#]/ )
              while ( FIRST_POSN == 1 ) {
                getline < "../coord.xvg" > 0
                FIRST_POSN=match( $0, /[\n@#]/ )
              }
	    split( $0, posns );
            dx=dims[irecs][2]/cols
            dy=dims[irecs][3]/rows
            dz=dims[irecs][4]/pags
            k=fr0
            for ( i=1; i<=34440+wa*3; i++ ) {
              idx0=((k+i)-1)*3+2;
              x=int((posns[idx0]  +2)/dx)-xs
              y=int((posns[idx0+1]+2)/dy)-ys
              z=int((posns[idx0+2]+2)/dz)-zs
              if ( x>cols-1 ) { x-=cols }
              if ( y>rows-1 ) { y-=rows }
              if ( z>pags-1 ) { z-=pags }
              if ( x<0      ) { x+=cols }
              if ( y<0      ) { y+=rows }
              if ( z<0      ) { z+=pags }
              chrg[z*iz+x+y*iy+1]+=chrgs[i];
	    }
            k=fr1
            for ( i=1; i<=34440+wa*3; i++ ) {
              idx0=((k+i)-1)*3+2;
              x=int((posns[idx0]  +2)/dx)-xs
              y=int((posns[idx0+1]+2)/dy)-ys
              z=int((posns[idx0+2]+2)/dz)-zs
              if ( x>cols-1 ) { x-=cols }
              if ( y>rows-1 ) { y-=rows }
              if ( z>pags-1 ) { z-=pags }
              if ( x<0      ) { x+=cols }
              if ( y<0      ) { y+=rows }
              if ( z<0      ) { z+=pags }
              chrg[z*iz+x+y*iy+1]+=chrgs[i]
	  } }
	  adjc=0
	  if ( substr(cfyl,1,4) != "nden" ) {
	    totc=0; ipos=0; bxls=0
	    for ( i=0; i<pags; i++ ) {
	      for ( j=0; j<rows; j++ ) {
	        for ( k=1; k<=cols; k++ ) {
		  ipos++
                  if ( chrg[ipos] != 0 ) {
                    totc+=chrg[ipos]
                    bxls++
            } } } }
            adjc=totc/bxls/frms
          }
	  printf("\n\n")
          ipos=0; zero=0
	  for ( i=0; i<pags; i++ ) {
	    for ( j=0; j<rows; j++ ) {
	      for ( k=1; k<=cols; k++ ) {
		ipos++
                if ( chrg[ipos] != 0 ) {
                  printf(" %9.6g",chrg[ipos]/frms-adjc) }
                else {
                  printf(" %9.6g",zero)
              } }
	      printf("\n")
	    }
	    printf("\n\n")
	} }
    '	../../../${ARG2} energy.xvg \
	> ${ARG3}
