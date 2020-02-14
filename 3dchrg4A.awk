#!/bin/bash -l
#===================================================================================
# 3dchrg4A.awk ### calculate avg charge in bins of volume for a trajectory       ###
#              ###                                                               ###
# parameters
# $1         number of horizontal (x-y) slices to create along z axis (2-?)
#   options: the output of this algorithm can be used to calculate non-uniform
#            electric fields or to produce computational x-rays. since the
#            computational cost of calculating electric fields scales as N6, and
#            that algorithm must process all the slices created here, balancing
#            resolution and computational time of both steps is required subject
#            to the requirements and limitations of each system on which it is used.
#            higher resolutions can be produced here if the intention is to
#            produce computational x-rays since they incure negligible scaling
# $2         name of file containing charges of each atom in system
# $3         output filname containing bins of average charge for each system volume
#            element
#___________________________________________________________________________________

ARG1=${1:-'50'}
ARG2=${2:-'chrgs2.xvg'}
ARG3=${3:-'chrg4.xvg'}

if [ -s "${ARG3}" ]; then s ../jobs/roll. ${ARG3}; fi

declare -a boxdims
boxdims=($(cat ../boxdims.xvg))

awk     -v PREC="quad" \
        -v zd=${boxdims[2]} -v xd=${boxdims[0]} -v yd=${boxdims[1]} \
        -v pags=${ARG1} -v cfyl=${ARG2} \
    '
BEGIN   { frms=2001
          dlns=0; drcs=0
          dza=zd/pags
          cols=int(xd/dza)+1; rows=int(yd/dza)+1
          dxa=xd/cols; dya=yd/rows
          xs=int(2/dxa); ys=int(2/dya); zs=int(2/dza)
          iz=cols*rows; iy=cols
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

END     \
        { for ( irecs=1; irecs<=frms; irecs++ ) {
            getline < "../coord.xvg" > 0
            FIRST_POSN=match( $0, /[\n@#]/ )
              while ( FIRST_POSN == 1 ) {
                getline < "../coord.xvg" > 0
                FIRST_POSN=match( $0, /[\n@#]/ )
              }
            split( $0, posns )
            ipos=0
            dx=dims[irecs][2]/cols
            dy=dims[irecs][3]/rows
            dz=dims[irecs][4]/pags
            for ( i=2; i<=NF; i+=3 ) {
              idx0=(i-2)+2
              x=int((posns[idx0]  +2)/dx)-xs
              y=int((posns[idx0+1]+2)/dy)-ys
              z=int((posns[idx0+2]+2)/dz)-zs
              if ( x>cols-1) { x-=cols }
              if ( y>rows-1) { y-=rows }
              if ( z>pags-1) { z-=pags }
              if ( x<0     ) { x+=cols }
              if ( y<0     ) { y+=rows }
              if ( z<0     ) { z+=pags }
              ipos++
              chrg[z*iz+x+y*iy+1]+=chrgs[ipos]
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
    '   ../../../${ARG2} energy.xvg \
        > ${ARG3}
