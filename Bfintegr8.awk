#!/usr/bin bash -l
#==============================================================
# Bfintegr8.awk ### performs brute force calculation of     ###
#               ### electric field at location of each bin  ###
#               ### in simulation from charge in each bin   ###
# parameters                                                
# $1    number of bins in z direction                       
# $2    charge file to integrate                            
# $3    output file name - outputs 2 sections, 1st is field on
#                          just that particular bin      
#                          2nd is total field, i.e. accumulated
#                          field from pb to this posn 
#   #!NOTE - the computational cost of this script scales as N6
#            of parameter $1. it gets very expensive very quick
#______________________________________________________________

ARG1=${1:-'50'}

declare -a boxdims
boxdims=($(cat ../boxdims.xvg))

awk     -v PREC="quad" \
        -v xd=${boxdims[0]} -v yd=${boxdims[1]} -v zd=${boxdims[2]} \
        -v pags=${ARG1} \
    '
BEGIN   { lns=0; dza=zd/pags;
          cols=int(xd/dza)+1; rows=int(yd/dza)+1
          dxa=xd/cols; dya=yd/rows
          divf=cols*rows*pags-1; ibxs=divf+1
          hx=(cols-1)/2; hy=(rows-1)/2; hz=(pags-1)/2
        }

BEGINFILES
        /^ / \
        { lns++
          vals[lns][""]
          split( $0, vals[lns])
        }
ENDFILES

END     \
        { FACT=18.0951281698764944
          for ( l=1; l<=pags; l++ ) {
            for ( m=1; m<=rows; m++ ) {
              indx=(l-1)*rows+m
              ft[indx][""]
              for ( n=1; n<=cols; n++ ) {
                lf=(n-1)*4+1
                for ( i=1; i<=pags; i++ ) {
                  for ( j=1; j<=rows; j++ ) {
                    for ( k=1; k<=cols; k++ ) {
                      lndx=(i-1)*rows+j
                      dstx=(n-k)*dxa; dsty=(m-j)*dya; dstz=(l-i)*dza
                      if ( k-n > hx || k-n < -hx ) {
                        if ( k < n ) {
                          dstx-=xd }
                        else {
                          dstx+=xd
                      } }
                      if ( j-m > hy || j-m < -hy ) {
                        if ( j < m ) {
                          dsty-=yd }
                        else {
                          dsty+=yd
                      } }
                      if ( i-l > hz || i-l < -hz ) {
                        if ( i < l ) {
                          dstz-=zd }
                        else {
                          dstz+=zd
                      } }
                      dist2=dstx^2+dsty^2+dstz^2
                      if ( dist2 != 0 ) {
                        qxk=vals[lndx][k]*FACT
                        ftot=qxk/dist2
                        dist=sqrt(dist2)
                        if ( ftot != 0 ) {
                          ft[indx][lf+1]+=ftot*dstx/dist
                          ft[indx][lf+2]+=ftot*dsty/dist
                          ft[indx][lf+3]+=ftot*dstz/dist
                } } } } }
                ft[indx][lf+1]/=divf
                tx+=ft[indx][lf+1]
                ft[indx][lf+2]/=divf
                ty+=ft[indx][lf+2]
                ft[indx][lf+3]/=divf
                tz+=ft[indx][lf+3]
          } } }
          adjx=-tx/ibxs
          adjy=-ty/ibxs
          adjz=-tz/ibxs

          printf("\n\n")
          for ( i=1; i<=pags; i++ ) {
            for ( j=1; j<=rows; j++ ) {
              indx=(i-1)*rows+j
              for ( k=1; k<=cols; k++ ) {
                lf=(k-1)*4+1
                ft[indx][lf+1]+=adjx
                ft[indx][lf+2]+=adjy
                ft[indx][lf+3]+=adjz
                ft[indx][lf]=sqrt((ft[indx][lf+1])^2+(ft[indx][lf+2])^2+(ft[indx][lf+3])^2)
                printf(" %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g\n",k*dxa,j*dya,i*dza,ft[indx][lf],ft[indx][lf+1],ft[indx][lf+2],ft[indx][lf+3])
            } }
            print("\n\n")
          }

          tx=0; ty=0; tz=0
          for ( l=1; l<=pags; l++ ) {
            for ( m=1; m<=rows; m++ ) {
              indx=(l-1)*rows+m
              for ( n=1; n<=cols; n++ ) {
                lf=(n-1)*4+1
                if ( m == 1 ) {
                  ft[indx][lf+1]+=ft[indx][lf-3]
                  tx+=ft[indx][lf+1]
                  ft[indx][lf+2]+=ft[indx+rows-1][lf+2]
                  ty+=ft[indx][lf+2]
                  ft[indx][lf+3]+=ft[indx-rows][lf+3]
                  tz+=ft[indx][lf+3] }
                else {
                  if ( m == rows ) {
                    ft[indx][lf+1]+=ft[indx][lf-3]
                    tx+=ft[indx][lf+1]
                    ft[indx][lf+2]+=ft[indx-rows+1][lf+2]
                    ty+=ft[indx][lf+2]
                    ft[indx][lf+3]+=ft[indx-rows][lf+3]
                    tz+=ft[indx][lf+3] }
                  else {
                    ft[indx][lf+1]+=ft[indx][lf-3]
                    tx+=ft[indx][lf+1]
                    ft[indx][lf+2]+=ft[indx-1][lf+2]
                    ty+=ft[indx][lf+2]
                    ft[indx][lf+3]+=ft[indx-rows][lf+3]
                    tz+=ft[indx][lf+3]
          } } } } }
          adjx=-tx/ibxs
          adjy=-ty/ibxs
          adjz=-tz/ibxs

          printf("\n\n")
          for ( i=1; i<=pags; i++ ) {
            for ( j=1; j<=rows; j++ ) {
              indx=(i-1)*rows+j
              for ( k=1; k<=cols; k++ ) {
                lf=(k-1)*4+1
                ft[indx][lf+1]+=adjx
                ft[indx][lf+2]+=adjy
                ft[indx][lf+3]+=adjz
                ft[indx][lf]=sqrt((ft[indx][lf+1])^2+(ft[indx][lf+2])^2+(ft[indx][lf+3])^2)
                printf(" %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g\n",k*dxa,j*dya,i*dza,ft[indx][lf],ft[indx][lf+1],ft[indx][lf+2],ft[indx][lf+3])
            } }
            printf("\n\n")
        } }
    '   ${2} > ${3}
