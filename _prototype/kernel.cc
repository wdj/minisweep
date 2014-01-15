//-----------------------------------------------------------------------------

#include <sys/time.h>
#include <stdlib.h>
#include <iostream>

//-----------------------------------------------------------------------------
//---1. Original code.

void sub1( double* v,
           double* w,
           double* matrix1,
           double* matrix2,
           double* sX,
           int nE,
           int nU,
           int nM,
           int nA )
{

#define SX(iU) ( sX[iU] )
#define V(iE,iU,iM) ( v[(iE)+nE*((iU)+nU*(iM))] )
#define W(iE,iU,iM) ( w[(iE)+nE*((iU)+nU*(iM))] )
#define Matrix1(iA,iM) ( matrix1[(iA)+nA*(iM)] )
#define Matrix2(iA,iM) ( matrix2[(iA)+nA*(iM)] )

  for ( int iE=0; iE<nE; iE++ ) {
  for ( int iA=0; iA<nA; iA++ ) {

    //---First matvec.
    for ( int iU=0; iU<nU; iU++ ) {
      SX(iU) = 0.;
      for ( int iM=0; iM<nM; iM++ ) {
        SX(iU) += Matrix1(iA,iM) * V(iE,iU,iM);
      }
    }

    //---<<<Other computations here>>>

    //---Second matvec.
    for ( int iU=0; iU<nU; iU++ ) {
      for ( int iM=0; iM<nM; iM++ ) {
        if ( iA == 0 ) {
          W(iE,iU,iM) = 0.;
        }
        W(iE,iU,iM) += Matrix2(iA,iM) * SX(iU);
      }
    }

  }
  }

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//---2. Specify values at compile time: nU=4, nM=16.

void sub2( double* v,
           double* w,
           double* matrix1,
           double* matrix2,
           double* sX,
           int nE,
           int nU,
           int nM,
           int nA )
{

#define NU 4
#define NM 16
#define SX(iU) ( sX[iU] )
#define V(iE,iU,iM)     ( v[(iE)+nE*((iU)+nU*(iM))] )
#define W(iE,iU,iM)     ( w[(iE)+nE*((iU)+nU*(iM))] )
#define V_NU4(iE,iU,iM) ( v[(iE)+nE*((iU)+NU*(iM))] )
#define W_NU4(iE,iU,iM) ( w[(iE)+nE*((iU)+NU*(iM))] )
#define Matrix1(iA,iM) ( matrix1[(iA)+nA*(iM)] )
#define Matrix2(iA,iM) ( matrix2[(iA)+nA*(iM)] )

  if ( nU == 4 && nM == 16 ) {

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      for ( int iU=0; iU<NU; iU++ ) {
        SX(iU) = 0.;
        for ( int iM=0; iM<NM; iM++ ) {
          SX(iU) += Matrix1(iA,iM) * V_NU4(iE,iU,iM);
        }
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        for ( int iM=0; iM<nM; iM++ ) {
          if ( iA == 0 ) {
            W_NU4(iE,iU,iM) = 0.;
          }
          W_NU4(iE,iU,iM) += Matrix2(iA,iM) * SX(iU);
        }
      }

    }
    }

  } else { //---general case.

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        SX(iU) = 0.;
        for ( int iM=0; iM<nM; iM++ ) {
          SX(iU) += Matrix1(iA,iM) * V(iE,iU,iM);
        }
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        for ( int iM=0; iM<nM; iM++ ) {
          if ( iA == 0 ) {
            W(iE,iU,iM) = 0.;
          }
          W(iE,iU,iM) += Matrix2(iA,iM) * SX(iU);
        }
      }

    }
    }

  }

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//---3. Fully unroll iU loops.

void sub3( double* v,
           double* w,
           double* matrix1,
           double* matrix2,
           double* sX,
           int nE,
           int nU,
           int nM,
           int nA )
{

#define NU 4
#define NM 16
#define SX(iU) ( sX[iU] )
#define V(iE,iU,iM)     ( v[(iE)+nE*((iU)+nU*(iM))] )
#define W(iE,iU,iM)     ( w[(iE)+nE*((iU)+nU*(iM))] )
#define V_NU4(iE,iU,iM) ( v[(iE)+nE*((iU)+NU*(iM))] )
#define W_NU4(iE,iU,iM) ( w[(iE)+nE*((iU)+NU*(iM))] )
#define Matrix1(iA,iM) ( matrix1[(iA)+nA*(iM)] )
#define Matrix2(iA,iM) ( matrix2[(iA)+nA*(iM)] )

  if ( nU == 4 && nM == 16 ) {

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      SX(0) = 0.;
      for ( int iM=0; iM<NM; iM++ ) {
        SX(0) += Matrix1(iA,iM) * V_NU4(iE,0,iM);
      }
      SX(1) = 0.;
      for ( int iM=0; iM<NM; iM++ ) {
        SX(1) += Matrix1(iA,iM) * V_NU4(iE,1,iM);
      }
      SX(2) = 0.;
      for ( int iM=0; iM<NM; iM++ ) {
        SX(2) += Matrix1(iA,iM) * V_NU4(iE,2,iM);
      }
      SX(3) = 0.;
      for ( int iM=0; iM<NM; iM++ ) {
        SX(3) += Matrix1(iA,iM) * V_NU4(iE,3,iM);
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      for ( int iM=0; iM<nM; iM++ ) {
        if ( iA == 0 ) {
          W_NU4(iE,0,iM) = 0.;
        }
        W_NU4(iE,0,iM) += Matrix2(iA,iM) * SX(0);
      }
      for ( int iM=0; iM<nM; iM++ ) {
        if ( iA == 0 ) {
          W_NU4(iE,1,iM) = 0.;
        }
        W_NU4(iE,1,iM) += Matrix2(iA,iM) * SX(1);
      }
      for ( int iM=0; iM<nM; iM++ ) {
        if ( iA == 0 ) {
          W_NU4(iE,2,iM) = 0.;
        }
        W_NU4(iE,2,iM) += Matrix2(iA,iM) * SX(2);
      }
      for ( int iM=0; iM<nM; iM++ ) {
        if ( iA == 0 ) {
          W_NU4(iE,3,iM) = 0.;
        }
        W_NU4(iE,3,iM) += Matrix2(iA,iM) * SX(3);
      }

    }
    }

  } else { //---general case.

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        SX(iU) = 0.;
        for ( int iM=0; iM<nM; iM++ ) {
          SX(iU) += Matrix1(iA,iM) * V(iE,iU,iM);
        }
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        for ( int iM=0; iM<nM; iM++ ) {
          if ( iA == 0 ) {
            W(iE,iU,iM) = 0.;
          }
          W(iE,iU,iM) += Matrix2(iA,iM) * SX(iU);
        }
      }

    }
    }

  }

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//---4. Reorder statements, fuse loops, fuse if-blocks.

void sub4( double* v,
           double* w,
           double* matrix1,
           double* matrix2,
           double* sX,
           int nE,
           int nU,
           int nM,
           int nA )
{

#define NU 4
#define NM 16
#define SX(iU) ( sX[iU] )
#define V(iE,iU,iM)     ( v[(iE)+nE*((iU)+nU*(iM))] )
#define W(iE,iU,iM)     ( w[(iE)+nE*((iU)+nU*(iM))] )
#define V_NU4(iE,iU,iM) ( v[(iE)+nE*((iU)+NU*(iM))] )
#define W_NU4(iE,iU,iM) ( w[(iE)+nE*((iU)+NU*(iM))] )
#define Matrix1(iA,iM) ( matrix1[(iA)+nA*(iM)] )
#define Matrix2(iA,iM) ( matrix2[(iA)+nA*(iM)] )

  if ( nU == 4 && nM == 16 ) {

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      SX(0) = 0.;
      SX(1) = 0.;
      SX(2) = 0.;
      SX(3) = 0.;
      for ( int iM=0; iM<NM; iM++ ) {
        SX(0) += Matrix1(iA,iM) * V_NU4(iE,0,iM);
        SX(1) += Matrix1(iA,iM) * V_NU4(iE,1,iM);
        SX(2) += Matrix1(iA,iM) * V_NU4(iE,2,iM);
        SX(3) += Matrix1(iA,iM) * V_NU4(iE,3,iM);
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      for ( int iM=0; iM<nM; iM++ ) {
        if ( iA == 0 ) {
          W_NU4(iE,0,iM) = 0.;
          W_NU4(iE,1,iM) = 0.;
          W_NU4(iE,2,iM) = 0.;
          W_NU4(iE,3,iM) = 0.;
        }
        W_NU4(iE,0,iM) += Matrix2(iA,iM) * SX(0);
        W_NU4(iE,1,iM) += Matrix2(iA,iM) * SX(1);
        W_NU4(iE,2,iM) += Matrix2(iA,iM) * SX(2);
        W_NU4(iE,3,iM) += Matrix2(iA,iM) * SX(3);
      }

    }
    }

  } else { //---general case.

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        SX(iU) = 0.;
        for ( int iM=0; iM<nM; iM++ ) {
          SX(iU) += Matrix1(iA,iM) * V(iE,iU,iM);
        }
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        for ( int iM=0; iM<nM; iM++ ) {
          if ( iA == 0 ) {
            W(iE,iU,iM) = 0.;
          }
          W(iE,iU,iM) += Matrix2(iA,iM) * SX(iU);
        }
      }

    }
    }

  }

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//---5. Convert SX array to registers.

void sub5( double* v,
           double* w,
           double* matrix1,
           double* matrix2,
           double* sX,
           int nE,
           int nU,
           int nM,
           int nA )
{

#define NU 4
#define NM 16
#define SX(iU) ( sX[iU] )
#define V(iE,iU,iM)     ( v[(iE)+nE*((iU)+nU*(iM))] )
#define W(iE,iU,iM)     ( w[(iE)+nE*((iU)+nU*(iM))] )
#define V_NU4(iE,iU,iM) ( v[(iE)+nE*((iU)+NU*(iM))] )
#define W_NU4(iE,iU,iM) ( w[(iE)+nE*((iU)+NU*(iM))] )
#define Matrix1(iA,iM) ( matrix1[(iA)+nA*(iM)] )
#define Matrix2(iA,iM) ( matrix2[(iA)+nA*(iM)] )

  if ( nU == 4 && nM == 16 ) {

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      double x0 = 0.;
      double x1 = 0.;
      double x2 = 0.;
      double x3 = 0.;
      for ( int iM=0; iM<NM; iM++ ) {
        x0 += Matrix1(iA,iM) * V_NU4(iE,0,iM);
        x1 += Matrix1(iA,iM) * V_NU4(iE,1,iM);
        x2 += Matrix1(iA,iM) * V_NU4(iE,2,iM);
        x3 += Matrix1(iA,iM) * V_NU4(iE,3,iM);
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      for ( int iM=0; iM<nM; iM++ ) {
        if ( iA == 0 ) {
          W_NU4(iE,0,iM) = 0.;
          W_NU4(iE,1,iM) = 0.;
          W_NU4(iE,2,iM) = 0.;
          W_NU4(iE,3,iM) = 0.;
        }
        W_NU4(iE,0,iM) += Matrix2(iA,iM) * x0;
        W_NU4(iE,1,iM) += Matrix2(iA,iM) * x1;
        W_NU4(iE,2,iM) += Matrix2(iA,iM) * x2;
        W_NU4(iE,3,iM) += Matrix2(iA,iM) * x3;
      }

    }
    }

  } else { //---general case.

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        SX(iU) = 0.;
        for ( int iM=0; iM<nM; iM++ ) {
          SX(iU) += Matrix1(iA,iM) * V(iE,iU,iM);
        }
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        for ( int iM=0; iM<nM; iM++ ) {
          if ( iA == 0 ) {
            W(iE,iU,iM) = 0.;
          }
          W(iE,iU,iM) += Matrix2(iA,iM) * SX(iU);
        }
      }

    }
    }

  }

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//---6. Replicate iM loop, to move if-test outside loop.  Replace iA with 0 when possible.

void sub6( double* v,
           double* w,
           double* matrix1,
           double* matrix2,
           double* sX,
           int nE,
           int nU,
           int nM,
           int nA )
{

#define NU 4
#define NM 16
#define SX(iU) ( sX[iU] )
#define V(iE,iU,iM)     ( v[(iE)+nE*((iU)+nU*(iM))] )
#define W(iE,iU,iM)     ( w[(iE)+nE*((iU)+nU*(iM))] )
#define V_NU4(iE,iU,iM) ( v[(iE)+nE*((iU)+NU*(iM))] )
#define W_NU4(iE,iU,iM) ( w[(iE)+nE*((iU)+NU*(iM))] )
#define Matrix1(iA,iM) ( matrix1[(iA)+nA*(iM)] )
#define Matrix2(iA,iM) ( matrix2[(iA)+nA*(iM)] )

  if ( nU == 4 && nM == 16 ) {

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      double x0 = 0.;
      double x1 = 0.;
      double x2 = 0.;
      double x3 = 0.;
      for ( int iM=0; iM<NM; iM++ ) {
        x0 += Matrix1(iA,iM) * V_NU4(iE,0,iM);
        x1 += Matrix1(iA,iM) * V_NU4(iE,1,iM);
        x2 += Matrix1(iA,iM) * V_NU4(iE,2,iM);
        x3 += Matrix1(iA,iM) * V_NU4(iE,3,iM);
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      if ( iA == 0 ) {
        for ( int iM=0; iM<nM; iM++ ) {
          W_NU4(iE,0,iM)  = Matrix2(0,iM) * x0;
          W_NU4(iE,1,iM)  = Matrix2(0,iM) * x1;
          W_NU4(iE,2,iM)  = Matrix2(0,iM) * x2;
          W_NU4(iE,3,iM)  = Matrix2(0,iM) * x3;
        }
      } else {
        for ( int iM=0; iM<nM; iM++ ) {
          W_NU4(iE,0,iM) += Matrix2(iA,iM) * x0;
          W_NU4(iE,1,iM) += Matrix2(iA,iM) * x1;
          W_NU4(iE,2,iM) += Matrix2(iA,iM) * x2;
          W_NU4(iE,3,iM) += Matrix2(iA,iM) * x3;
        }
      }

    }
    }

  } else { //---general case.

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        SX(iU) = 0.;
        for ( int iM=0; iM<nM; iM++ ) {
          SX(iU) += Matrix1(iA,iM) * V(iE,iU,iM);
        }
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        for ( int iM=0; iM<nM; iM++ ) {
          if ( iA == 0 ) {
            W(iE,iU,iM) = 0.;
          }
          W(iE,iU,iM) += Matrix2(iA,iM) * SX(iU);
        }
      }

    }
    }

  }

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//---7. Change array axis order, for this code and caller.

void sub7( double* v,
           double* w,
           double* matrix1,
           double* matrix2,
           double* sX,
           int nE,
           int nU,
           int nM,
           int nA )
{

#define NU 4
#define NM 16
#define SX(iU) ( sX[iU] )
#define V(iE,iU,iM)     ( v[(iE)+nE*((iU)+nU*(iM))] )
#define W(iE,iU,iM)     ( w[(iE)+nE*((iU)+nU*(iM))] )
#define V_NU4_PERMUTED(iE,iU,iM) ( v[(iM)+NM*((iU)+NU*(iE))] )
#define W_NU4_PERMUTED(iE,iU,iM) ( w[(iM)+NM*((iU)+NU*(iE))] )
#define Matrix1(iA,iM) ( matrix1[(iA)+nA*(iM)] )
#define Matrix2(iA,iM) ( matrix2[(iA)+nA*(iM)] )

  if ( nU == 4 && nM == 16 ) {

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      double x0 = 0.;
      double x1 = 0.;
      double x2 = 0.;
      double x3 = 0.;
      for ( int iM=0; iM<NM; iM++ ) {
        x0 += Matrix1(iA,iM) * V_NU4_PERMUTED(iE,0,iM);
        x1 += Matrix1(iA,iM) * V_NU4_PERMUTED(iE,1,iM);
        x2 += Matrix1(iA,iM) * V_NU4_PERMUTED(iE,2,iM);
        x3 += Matrix1(iA,iM) * V_NU4_PERMUTED(iE,3,iM);
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      if ( iA == 0 ) {
        for ( int iM=0; iM<nM; iM++ ) {
          W_NU4_PERMUTED(iE,0,iM)  = Matrix2(0,iM) * x0;
          W_NU4_PERMUTED(iE,1,iM)  = Matrix2(0,iM) * x1;
          W_NU4_PERMUTED(iE,2,iM)  = Matrix2(0,iM) * x2;
          W_NU4_PERMUTED(iE,3,iM)  = Matrix2(0,iM) * x3;
        }
      } else {
        for ( int iM=0; iM<nM; iM++ ) {
          W_NU4_PERMUTED(iE,0,iM) += Matrix2(iA,iM) * x0;
          W_NU4_PERMUTED(iE,1,iM) += Matrix2(iA,iM) * x1;
          W_NU4_PERMUTED(iE,2,iM) += Matrix2(iA,iM) * x2;
          W_NU4_PERMUTED(iE,3,iM) += Matrix2(iA,iM) * x3;
        }
      }

    }
    }

  } else { //---general case.

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        SX(iU) = 0.;
        for ( int iM=0; iM<nM; iM++ ) {
          SX(iU) += Matrix1(iA,iM) * V(iE,iU,iM);
        }
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        for ( int iM=0; iM<nM; iM++ ) {
          if ( iA == 0 ) {
            W(iE,iU,iM) = 0.;
          }
          W(iE,iU,iM) += Matrix2(iA,iM) * SX(iU);
        }
      }

    }
    }

  }

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//---8. Extract V, U into temporary arrays, for more locality, reuse.

void sub8( double* v,
           double* w,
           double* matrix1,
           double* matrix2,
           double* sX,
           int nE,
           int nU,
           int nM,
           int nA )
{

#define NU 4
#define NM 16
#define SX(iU) ( sX[iU] )
#define V(iE,iU,iM)     ( v[(iE)+nE*((iU)+nU*(iM))] )
#define W(iE,iU,iM)     ( w[(iE)+nE*((iU)+nU*(iM))] )
#define V_NU4_PERMUTED(iE,iU,iM) ( v[(iM)+NM*((iU)+NU*(iE))] )
#define W_NU4_PERMUTED(iE,iU,iM) ( w[(iM)+NM*((iU)+NU*(iE))] )
#define Matrix1(iA,iM) ( matrix1[(iA)+nA*(iM)] )
#define Matrix2(iA,iM) ( matrix2[(iA)+nA*(iM)] )
#define SV(iU,iM) ( sV[(iU)+NU*(iM)] )
#define SW(iU,iM) ( sW[(iU)+NU*(iM)] )

  if ( nU == 4 && nM == 16 ) {

    double sV[ NU * NM ];
    double sW[ NU * NM ];

    for ( int iE=0; iE<nE; iE++ ) {

      //---Explicitly load V values, for reuse across iA.
      for ( int iM=0; iM<nM; iM++ ) {
        SV(0,iM) = V_NU4_PERMUTED(iE,0,iM);
        SV(1,iM) = V_NU4_PERMUTED(iE,1,iM);
        SV(2,iM) = V_NU4_PERMUTED(iE,2,iM);
        SV(3,iM) = V_NU4_PERMUTED(iE,3,iM);
      }

      for ( int iA=0; iA<nA; iA++ ) {

        //---First matvec.
        double x0 = 0.;
        double x1 = 0.;
        double x2 = 0.;
        double x3 = 0.;
        for ( int iM=0; iM<NM; iM++ ) {
          x0 += Matrix1(iA,iM) * SV(0,iM);
          x1 += Matrix1(iA,iM) * SV(1,iM);
          x2 += Matrix1(iA,iM) * SV(2,iM);
          x3 += Matrix1(iA,iM) * SV(3,iM);
        }

        //---<<<Other computations here>>>

        //---Second matvec.
        if ( iA == 0 ) {
          for ( int iM=0; iM<nM; iM++ ) {
            SW(0,iM)  = Matrix2(0,iM) * x0;
            SW(1,iM)  = Matrix2(0,iM) * x1;
            SW(2,iM)  = Matrix2(0,iM) * x2;
            SW(3,iM)  = Matrix2(0,iM) * x3;
          }
        } else {
          for ( int iM=0; iM<nM; iM++ ) {
            SW(0,iM) += Matrix2(iA,iM) * x0;
            SW(1,iM) += Matrix2(iA,iM) * x1;
            SW(2,iM) += Matrix2(iA,iM) * x2;
            SW(3,iM) += Matrix2(iA,iM) * x3;
          }
        }

      }

      //---Explicitly store W values, reused across iA.
      for ( int iM=0; iM<nM; iM++ ) {
        W_NU4_PERMUTED(iE,0,iM) = SW(0,iM);
        W_NU4_PERMUTED(iE,1,iM) = SW(1,iM);
        W_NU4_PERMUTED(iE,2,iM) = SW(2,iM);
        W_NU4_PERMUTED(iE,3,iM) = SW(3,iM);
      }

    }

  } else { //---general case.

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        SX(iU) = 0.;
        for ( int iM=0; iM<nM; iM++ ) {
          SX(iU) += Matrix1(iA,iM) * V(iE,iU,iM);
        }
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        for ( int iM=0; iM<nM; iM++ ) {
          if ( iA == 0 ) {
            W(iE,iU,iM) = 0.;
          }
          W(iE,iU,iM) += Matrix2(iA,iM) * SX(iU);
        }
      }

    }
    }

  }

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//---9. Unroll most expensive iM loops to depth 4 - note no cleanup code needed.

void sub9( double* v,
           double* w,
           double* matrix1,
           double* matrix2,
           double* sX,
           int nE,
           int nU,
           int nM,
           int nA )
{

#define NU 4
#define NM 16
#define SX(iU) ( sX[iU] )
#define V(iE,iU,iM)     ( v[(iE)+nE*((iU)+nU*(iM))] )
#define W(iE,iU,iM)     ( w[(iE)+nE*((iU)+nU*(iM))] )
#define V_NU4_PERMUTED(iE,iU,iM) ( v[(iM)+NM*((iU)+NU*(iE))] )
#define W_NU4_PERMUTED(iE,iU,iM) ( w[(iM)+NM*((iU)+NU*(iE))] )
#define Matrix1(iA,iM) ( matrix1[(iA)+nA*(iM)] )
#define Matrix2(iA,iM) ( matrix2[(iA)+nA*(iM)] )
#define SV(iU,iM) ( sV[(iU)+NU*(iM)] )
#define SW(iU,iM) ( sW[(iU)+NU*(iM)] )

  if ( nU == 4 && nM == 16 ) {

    double sV[ NU * NM ];
    double sW[ NU * NM ];

    for ( int iE=0; iE<nE; iE++ ) {

      //---Explicitly load V values, for reuse across iA.
      for ( int iM=0; iM<nM; iM++ ) {
        SV(0,iM) = V_NU4_PERMUTED(iE,0,iM);
        SV(1,iM) = V_NU4_PERMUTED(iE,1,iM);
        SV(2,iM) = V_NU4_PERMUTED(iE,2,iM);
        SV(3,iM) = V_NU4_PERMUTED(iE,3,iM);
      }

      for ( int iA=0; iA<nA; iA++ ) {

        //---First matvec.
        double x0 = 0.;
        double x1 = 0.;
        double x2 = 0.;
        double x3 = 0.;
        for ( int iM=0; iM<NM; iM+=4 ) {
          x0 += Matrix1(iA,iM  ) * SV(0,iM  )
              + Matrix1(iA,iM+1) * SV(0,iM+1)
              + Matrix1(iA,iM+2) * SV(0,iM+2)
              + Matrix1(iA,iM+3) * SV(0,iM+3);
          x1 += Matrix1(iA,iM  ) * SV(1,iM  )
              + Matrix1(iA,iM+1) * SV(1,iM+1)
              + Matrix1(iA,iM+2) * SV(1,iM+2)
              + Matrix1(iA,iM+3) * SV(1,iM+3);
          x2 += Matrix1(iA,iM  ) * SV(2,iM  )
              + Matrix1(iA,iM+1) * SV(2,iM+1)
              + Matrix1(iA,iM+2) * SV(2,iM+2)
              + Matrix1(iA,iM+3) * SV(2,iM+3);
          x3 += Matrix1(iA,iM  ) * SV(3,iM  )
              + Matrix1(iA,iM+1) * SV(3,iM+1)
              + Matrix1(iA,iM+2) * SV(3,iM+2)
              + Matrix1(iA,iM+3) * SV(3,iM+3);
        }

        //---<<<Other computations here>>>

        //---Second matvec.
        if ( iA == 0 ) {
          for ( int iM=0; iM<nM; iM++ ) {
            SW(0,iM)  = Matrix2(0,iM) * x0;
            SW(1,iM)  = Matrix2(0,iM) * x1;
            SW(2,iM)  = Matrix2(0,iM) * x2;
            SW(3,iM)  = Matrix2(0,iM) * x3;
          }
        } else {
          for ( int iM=0; iM<nM; iM+=4 ) {
            SW(0,iM  ) += Matrix2(iA,iM  ) * x0;
            SW(0,iM+1) += Matrix2(iA,iM+1) * x0;
            SW(0,iM+2) += Matrix2(iA,iM+2) * x0;
            SW(0,iM+3) += Matrix2(iA,iM+3) * x0;
            SW(1,iM  ) += Matrix2(iA,iM  ) * x1;
            SW(1,iM+1) += Matrix2(iA,iM+1) * x1;
            SW(1,iM+2) += Matrix2(iA,iM+2) * x1;
            SW(1,iM+3) += Matrix2(iA,iM+3) * x1;
            SW(2,iM  ) += Matrix2(iA,iM  ) * x2;
            SW(2,iM+1) += Matrix2(iA,iM+1) * x2;
            SW(2,iM+2) += Matrix2(iA,iM+2) * x2;
            SW(2,iM+3) += Matrix2(iA,iM+3) * x2;
            SW(3,iM  ) += Matrix2(iA,iM  ) * x3;
            SW(3,iM+1) += Matrix2(iA,iM+1) * x3;
            SW(3,iM+2) += Matrix2(iA,iM+2) * x3;
            SW(3,iM+3) += Matrix2(iA,iM+3) * x3;
          }
        }

      }

      //---Explicitly store W values, reused across iA.
      for ( int iM=0; iM<nM; iM++ ) {
        W_NU4_PERMUTED(iE,0,iM) = SW(0,iM);
        W_NU4_PERMUTED(iE,1,iM) = SW(1,iM);
        W_NU4_PERMUTED(iE,2,iM) = SW(2,iM);
        W_NU4_PERMUTED(iE,3,iM) = SW(3,iM);
      }

    }

  } else { //---general case.

    for ( int iE=0; iE<nE; iE++ ) {
    for ( int iA=0; iA<nA; iA++ ) {

      //---First matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        SX(iU) = 0.;
        for ( int iM=0; iM<nM; iM++ ) {
          SX(iU) += Matrix1(iA,iM) * V(iE,iU,iM);
        }
      }

      //---<<<Other computations here>>>

      //---Second matvec.
      for ( int iU=0; iU<nU; iU++ ) {
        for ( int iM=0; iM<nM; iM++ ) {
          if ( iA == 0 ) {
            W(iE,iU,iM) = 0.;
          }
          W(iE,iU,iM) += Matrix2(iA,iM) * SX(iU);
        }
      }

    }
    }

  }

}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

static double getWallclockTime()
{

  struct timeval tp;
  struct timezone tzp;
  int i = gettimeofday( &tp, &tzp );
  return ( (double) tp.tv_sec +
           (double) tp.tv_usec * 1.e-6 );

}


//-----------------------------------------------------------------------------

int main( int argc, char** argv )
{

  const int nE = argc > 1 ? atoi(argv[1]) : 1000 * 100;
  const int nU = argc > 2 ? atoi(argv[2]) : 4;
  const int nM = argc > 3 ? atoi(argv[3]) : 16;
  const int nA = argc > 4 ? atoi(argv[4]) : 32;

  double* v       = new double[ nE * nU * nM ];
  double* w       = new double[ nE * nU * nM ];
  double* matrix1 = new double[ nA * nM ];
  double* matrix2 = new double[ nA * nM ];
  double* sX      = new double[ nU ];

  for ( int i=0; i<nE*nU*nM; ++i )
  {
    v[i] = 1. * nE;
  }

  for ( int i=0; i<nA*nM; ++i )
  {
    matrix1[i] = 1.e-2 / ( 3. * nA * ( i + 1 ) );
    matrix2[i] = 1.e-2 / ( 4. * nM + i );
  }

  //---Run once to touch memory.

  sub1( v, w, matrix1, matrix2, sX, nE, nU, nM, nA );

  //--------------------

  {
    double t1 = getWallclockTime();
    sub1( v, w, matrix1, matrix2, sX, nE, nU, nM, nA );
    double t2 = getWallclockTime();

    double sum = 0.;
    for ( int i=0; i<nE*nU*nM; ++i )
    {
      sum += w[i] * w[i];
    }

    std::cout << "Kernel 1, time " << t2 - t1 << " sec, sum " << sum << ".\n";

  }

  //--------------------

  {
    double t1 = getWallclockTime();
    sub2( v, w, matrix1, matrix2, sX, nE, nU, nM, nA );
    double t2 = getWallclockTime();

    double sum = 0.;
    for ( int i=0; i<nE*nU*nM; ++i )
    {
      sum += w[i] * w[i];
    }

    std::cout << "Kernel 2, time " << t2 - t1 << " sec, sum " << sum << ".\n";

  }

  //--------------------

  {
    double t1 = getWallclockTime();
    sub3( v, w, matrix1, matrix2, sX, nE, nU, nM, nA );
    double t2 = getWallclockTime();

    double sum = 0.;
    for ( int i=0; i<nE*nU*nM; ++i )
    {
      sum += w[i] * w[i];
    }

    std::cout << "Kernel 3, time " << t2 - t1 << " sec, sum " << sum << ".\n";

  }

  //--------------------

  {
    double t1 = getWallclockTime();
    sub4( v, w, matrix1, matrix2, sX, nE, nU, nM, nA );
    double t2 = getWallclockTime();

    double sum = 0.;
    for ( int i=0; i<nE*nU*nM; ++i )
    {
      sum += w[i] * w[i];
    }

    std::cout << "Kernel 4, time " << t2 - t1 << " sec, sum " << sum << ".\n";

  }

  //--------------------

  {
    double t1 = getWallclockTime();
    sub5( v, w, matrix1, matrix2, sX, nE, nU, nM, nA );
    double t2 = getWallclockTime();

    double sum = 0.;
    for ( int i=0; i<nE*nU*nM; ++i )
    {
      sum += w[i] * w[i];
    }

    std::cout << "Kernel 5, time " << t2 - t1 << " sec, sum " << sum << ".\n";

  }

  //--------------------

  {
    double t1 = getWallclockTime();
    sub6( v, w, matrix1, matrix2, sX, nE, nU, nM, nA );
    double t2 = getWallclockTime();

    double sum = 0.;
    for ( int i=0; i<nE*nU*nM; ++i )
    {
      sum += w[i] * w[i];
    }

    std::cout << "Kernel 6, time " << t2 - t1 << " sec, sum " << sum << ".\n";

  }

  //--------------------

  {
    double t1 = getWallclockTime();
    sub7( v, w, matrix1, matrix2, sX, nE, nU, nM, nA );
    double t2 = getWallclockTime();

    double sum = 0.;
    for ( int i=0; i<nE*nU*nM; ++i )
    {
      sum += w[i] * w[i];
    }

    std::cout << "Kernel 7, time " << t2 - t1 << " sec, sum " << sum << ".\n";

  }

  //--------------------

  {
    double t1 = getWallclockTime();
    sub8( v, w, matrix1, matrix2, sX, nE, nU, nM, nA );
    double t2 = getWallclockTime();

    double sum = 0.;
    for ( int i=0; i<nE*nU*nM; ++i )
    {
      sum += w[i] * w[i];
    }

    std::cout << "Kernel 8, time " << t2 - t1 << " sec, sum " << sum << ".\n";

  }

  //--------------------

  {
    double t1 = getWallclockTime();
    sub9( v, w, matrix1, matrix2, sX, nE, nU, nM, nA );
    double t2 = getWallclockTime();

    double sum = 0.;
    for ( int i=0; i<nE*nU*nM; ++i )
    {
      sum += w[i] * w[i];
    }

    std::cout << "Kernel 9, time " << t2 - t1 << " sec, sum " << sum << ".\n";

  }

  //--------------------

  delete [] sX;
  delete [] matrix1;
  delete [] matrix2;
  delete [] w;
  delete [] v;

  return 0;

}

//-----------------------------------------------------------------------------
