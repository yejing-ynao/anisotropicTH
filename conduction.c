#include "../copyright.h"
/*============================================================================*/
/*! \file aniconduction.c
 *  \brief Adds semi-implicit thermal conduction term to the energy equation,
 *      dE/dt = Div(Q)
 *
 *   where
 *    - Q = spitzer model = heat flux
 *    - T = (P/d)*(mbar/k_B) = temperature
 *    - b = magnetic field unit vector
 *

 *
 * Note the kappa's are DIFFUSIVITIES, not CONDUCTIVITIES.  Also note this
 * version uses "dimensionless units" in that the factor (mbar/k_B) is not
 * included in calculating the temperature (instead, T=P/d is used).  For cgs
 * units, kappa must be entered in units of [cm^2/s], and the heat fluxes would
 * need to be multiplied by (k_B/mbar).
 *
 * The heat flux Q is calculated by calls to HeatFlux_* functions.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - conduction() - updates energy equation with thermal conduction
 * - conduction_init() - allocates memory needed
 * - conduction_destruct() - frees memory used */
/*============================================================================*/
#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef THERMAL_CONDUCTION

#ifdef BAROTROPIC
#error : Thermal conduction requires an adiabatic EOS
#endif
/*time step number*/
 int ts=0;
/* Arrays for the temperature and heat fluxes */
static Real ***Temp=NULL;//size: Nz*Ny*Nx
static Real **Temp1=NULL;//size: Ny*Nx
static Real **Temp2=NULL;//size: Ny*Nx
static Real ***Temp11=NULL;//size: 2*Ny*Nx for registering the temperature at n and n-1
static Real ***Temp22=NULL;//size: 2*Ny*Nx for registering the temperature at n and n-1
static Real ***Temp_=NULL;//size: Nz*Ny*Nx

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   HeatFlux_iso   - computes   isotropic heat flux
 *   HeatFlux_aniso - computes anisotropic heat flux
 *============================================================================*/


void tdma(Real x[], const int N, const Real a[], const Real b[], Real c[]);
static Real k_nl(const Real T);
static Real T_ex(Real ***T, const int Nt, const int j, const int i);

static Real limiter2(const Real A, const Real B);
static Real limiter4(const Real A, const Real B, const Real C, const Real D);
static Real vanleer (const Real A, const Real B);
static Real minmod  (const Real A, const Real B);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void conduction(DomainS *pD)
 *  \brief Semi-implicit thermal conduction
 */
void conduction(DomainS *pD)
{
    GridS *pG = (pD->Grid);
    int i, is = pG->is, ie = pG->ie;
    int j, jl, ju, js = pG->js, je = pG->je;
    int k, kl, ku, ks = pG->ks, ke = pG->ke;
    int ii, jj;
    int myL,myM,myN;
    int Nx1, Nx2 ,Ny1, Ny2;
    int co,co1;
    int np;
    int Nt=ts%2;
    Real dTdx1, dTdx2;
    Real dTdy1, dTdy2;
    Real x1[pG->Nx[0]*pD->NGrid[0]], a1[pG->Nx[0]*pD->NGrid[0]], b1[pG->Nx[0]*pD->NGrid[0]], c1[pG->Nx[0]*pD->NGrid[0]];
    Real x2[pG->Nx[1]*pD->NGrid[1]], a2[pG->Nx[1]*pD->NGrid[1]], b2[pG->Nx[1]*pD->NGrid[1]], c2[pG->Nx[1]*pD->NGrid[1]];
    Real T_cut[pG->Nx[1]/pD->NGrid[0]+2][pG->Nx[0]], T_cut1[pG->Nx[1]/pD->NGrid[0]+2][pG->Nx[0]];
    Real T_cut2[pG->Nx[1]/pD->NGrid[0]][pG->Nx[0]/pD->NGrid[1]+2], T_cut3[pG->Nx[1]/pD->NGrid[0]][pG->Nx[0]/pD->NGrid[1]+2];
    Real T_cut4[pG->Nx[1]][pG->Nx[0]/pD->NGrid[1]], T_cut5[pG->Nx[1]][pG->Nx[0]/pD->NGrid[1]];
    Real B01, B02;
    MPI_Request req1[pD->NGrid[0]],req2[pD->NGrid[1]], req3[pD->NGrid[1]], req4[pD->NGrid[1]];
    Real Lambda;
    Lambda = 1.67e-27/2./1.38e-23;
    Real dx1=pG->dx1, dx2=pG->dx2;
#ifdef STS
    Real my_dt = STS_dt;
#else
    Real my_dt = pG->dt;
#endif
    if (pG->Nx[1] > 1){
        jl = js - 1;
        ju = je + 1;
    } else {
        jl = js;
        ju = je;
    }
    if (pG->Nx[2] > 1){
        kl = ks - 1;
        ku = ke + 1;
    } else {
        kl = ks;
        ku = ke;
    }
    /* get (l,m,n) coordinates of Grid being updated on this processor */
    
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
  
    /*2-steps ADI method in two dimensions*/
    /*------------------------------update in x direction--------------------------------*/
    /*Domain decomposition along y-axis*/
    //rearrange Temp from other ID
    
    /* Zero heat heat flux array; compute temperature at cell centers.  Temperature
     * includes a factor [k_B/mbar].  For cgs units, the heat flux would have to
     * be multiplied by this factor.
     */
    
   /* for (k=kl; k<=ku; k++) {
        for (j=jl+myM*pG->Nx[1]/pD->NGrid[0]; j<=jl+(myM+1)*pG->Nx[1]/pD->NGrid[0]+1; j++) {
            for (i=is-Nx1*pG->Nx[0]; i<=ie+Nx2*pG->Nx[0]; i++) {
                
                Q[k][j][i].x1 = 0.0;
                Q[k][j][i].x2 = 0.0;
                Q[k][j][i].x3 = 0.0;
                
                Temp1[Nt][k][j][i] = pG->U[k][j][i].E - (0.5/pG->U[k][j][i].d)*
                (SQR(pG->U[k][j][i].M1) +SQR(pG->U[k][j][i].M2) +SQR(pG->U[k][j][i].M3));
#ifdef MHD
                Temp1[Nt][k][j][i] -= (0.5)*(SQR(pG->U[k][j][i].B1c) +
                                        SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
#endif
                Temp1[Nt][k][j][i] *= (Gamma_1/pG->U[k][j][i].d)*Lambda;
                //boundary conditions dT/dy=0
             if (myM==0 && myN==0)
                 Temp1[Nt][k][jl][i]=Temp1[Nt][k][js][i];
                if (myM==pD->NGrid[1]-1 && myN==pD->NGrid[0]-1)
                    Temp1[Nt][k][ju][i]=Temp1[Nt][k][je][i];
            }}}*/
    for (k=kl; k<=ku; k++) {
        for (j=jl; j<=ju; j++) {
            for (i=is-1; i<=ie+1; i++) {
                
                
                Temp[k][j][i] = pG->U[k][j][i].E - (0.5/pG->U[k][j][i].d)*
                (SQR(pG->U[k][j][i].M1) +SQR(pG->U[k][j][i].M2) +SQR(pG->U[k][j][i].M3));
#ifdef MHD
                Temp[k][j][i] -= (0.5)*(SQR(pG->U[k][j][i].B1c) +
                                        SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
#endif
                Temp[k][j][i] *= (Gamma_1/pG->U[k][j][i].d*Lambda);
                
            }}}
    //MPI arrangement
    //send data
    for (np=0; np<=pD->NGrid[0]-1; np++ ){
        
        for (k=kl; k<=ku; k++){
            co=0;
            for (j=jl+np*pG->Nx[1]/pD->NGrid[0]; j<=jl+(np+1)*pG->Nx[1]/pD->NGrid[0]+1; j++){
                 co1=0;
                for (i=is; i<=ie; i++){
                    T_cut[co][co1++]=Temp[k][j][i];
                }
                co++;
            }
        }
        MPI_Isend(&T_cut[0][0], (pG->Nx[1]/pD->NGrid[0]+2)*pG->Nx[0], MPI_DOUBLE, pD->GData[myL][myM][np].ID_Comm_Domain, myN, MPI_COMM_WORLD, &req1[np]);
    }
    
    //receive data
      for (np=0; np<=pD->NGrid[0]-1; np++ ){
           MPI_Irecv(&T_cut1[0][0], (pG->Nx[1]/pD->NGrid[0]+2)*pG->Nx[0], MPI_DOUBLE, myID_Comm_world, np, MPI_COMM_WORLD, &req2[np]);
       
       
            co=0;
            for (j=jl; j<=jl+pG->Nx[1]/pD->NGrid[0]+1; j++){
                co1=0;
                for (i=is; i<=ie; i++){
                    Temp1[j][i+np*pG->Nx[0]]=T_cut1[co][co1++];
                    if (ts==0){
                        Temp11[0][j][i+np*pG->Nx[0]]=Temp1[j][i+np*pG->Nx[0]];
                        Temp11[1][j][i+np*pG->Nx[0]]=Temp1[j][i+np*pG->Nx[0]];
                    }
                }
                co++;
            }
        
        
    }
   //tridiagonal matrix
    //thomas algorithm in x-direction
    for (k=ks; k<=ke; k++){
    for (j=js; j<=js+pG->Nx[1]/pD->NGrid[0]-1; j++){
        co=0;
        for (i=is; i<=is+pD->NGrid[0]*pG->Nx[0]-1; i++){
            B01 = SQR(pG->B1i[k][j][i])+SQR(0.5*(pG->U[k][j][i+1].B2c + pG->U[k][j][i].B2c));
            B01 = MAX(B01,TINY_NUMBER); /* limit in case B=0 */
            B02 = SQR(pG->B1i[k][j][i])+SQR(0.5*(pG->U[k][j][i-1].B2c + pG->U[k][j][i].B2c));
            B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */
            if (i==is){
                dTdy2 = limiter2(Temp1[j+1][i+1] - Temp1[j  ][i+1],
                                 Temp1[j  ][i+1] - Temp1[j-1][i+1]);
            a1[co]=0;
                
            }
            else{
                    dTdy2 = limiter4(Temp1[j+1][i  ] - Temp1[j  ][i  ],
                                      Temp1[j  ][i  ] - Temp1[j-1][i  ],
                                      Temp1[j+1][i-1] - Temp1[j  ][i-1],
                                     Temp1[j  ][i-1] - Temp1[j-1][i-1]);
                 a1[co]=-Lambda*Gamma_1/pG->U[k][j][i].d*pG->B1i[k][j][i]*pG->B1i[k][j][i]/B02*k_nl((T_ex(Temp11,Nt,j,i-1)+T_ex(Temp11,Nt,j,i))/2.)*my_dt/dx1/dx1;
            }
            if (i==is+pD->NGrid[0]*pG->Nx[0]-1){
                dTdy1 = limiter2(Temp1[j+1][i  ] - Temp1[j  ][i  ],
                                 Temp1[j  ][i  ] - Temp1[j-1][i  ]);
                 c1[co]=0;
            }
            else{
            dTdy1 = limiter4(Temp1[j+1][i+1] - Temp1[j  ][i+1],
                             Temp1[j  ][i+1] - Temp1[j-1][i+1],
                             Temp1[j+1][i  ] - Temp1[j  ][i  ],
                             Temp1[j  ][i  ] - Temp1[j-1][i  ]);
               c1[co]=-Lambda*Gamma_1/pG->U[k][j][i].d*(pG->B1i[k][j][i+1]*pG->B1i[k][j][i+1]/B01*k_nl((T_ex(Temp11,Nt,j,i)+T_ex(Temp11,Nt,j,i+1))/2.)*my_dt/dx1/dx1);
            }
            
            x1[co]=Temp1[j][i]+Lambda*Gamma_1/pG->U[k][j][i].d*(pG->B1i[k][j][i+1]*pG->B2i[k][j][i+1]/B01*k_nl((T_ex(Temp11,Nt,j,i)+T_ex(Temp11,Nt,j,i+1))/2.)*dTdy1*my_dt/dx1/dx2-pG->B1i[k][j][i]*pG->B2i[k][j][i]/B02*k_nl((T_ex(Temp11,Nt,j,i-1)+T_ex(Temp11,Nt,j,i))/2.)*dTdy2*my_dt/dx1/dx2);
            if (i==is)
           b1[co]=1+Lambda*Gamma_1/pG->U[k][j][i].d*pG->B1i[k][j][i+1]*pG->B1i[k][j][i+1]/B01*k_nl((T_ex(Temp11,Nt,j,i)+T_ex(Temp11,Nt,j,i+1))/2.)*my_dt/dx1/dx1;
            else if (i==is+pD->NGrid[0]*pG->Nx[0]-1)
           b1[co]=1+Lambda*Gamma_1/pG->U[k][j][i].d*pG->B1i[k][j][i]*pG->B1i[k][j][i]/B02*k_nl((T_ex(Temp11,Nt,j,i-1)+T_ex(Temp11,Nt,j,i))/2.)*my_dt/dx1/dx1;
            else
            b1[co]=1+Lambda*Gamma_1/pG->U[k][j][i].d*(pG->B1i[k][j][i]*pG->B1i[k][j][i]/B02*k_nl((T_ex(Temp11,Nt,j,i-1)+T_ex(Temp11,Nt,j,i))/2.)*my_dt/dx1/dx1+pG->B1i[k][j][i+1]*pG->B1i[k][j][i+1]/B01*k_nl((T_ex(Temp11,Nt,j,i)+T_ex(Temp11,Nt,j,i+1))/2.)*my_dt/dx1/dx1);
            co++;
        }
            tdma(x1,pG->Nx[0]*pD->NGrid[0],a1,b1,c1);
            co1=0;
           for (i=is; i<=is+pD->NGrid[0]*pG->Nx[0]-1; i++){
               Temp1[j][i]=x1[co1++];
               if (i==is)
                   Temp1[j][i-1]=Temp1[j][i];
               if (i==is+pD->NGrid[0]*pG->Nx[0]-1)
                   Temp1[j][i+1]=Temp1[j][i];
           }
                                                                       
                                                                       
        }
    }
  /*------------------------------update in y direction--------------------------------*/
 //send Temp1 to Temp2 for y-direction

for (np=0; np<=pD->NGrid[0]*pD->NGrid[1]-1; np++){
        co=0;
        for (j=js; j<=js+pG->Nx[1]/pD->NGrid[0]-1; j++){
            co1=0;
            for (i=is-1+pG->Nx[0]/pD->NGrid[1]*np; i<=is+pG->Nx[0]/pD->NGrid[1]*(np+1); i++){
                T_cut2[co][co1++]=Temp1[j][i];
            }
        }
        co++;
    MPI_Isend(&T_cut2[0][0], pG->Nx[1]/pD->NGrid[0]*(pG->Nx[0]/pD->NGrid[1]+2), MPI_DOUBLE, pD->GData[myL][np/pD->NGrid[0]][np%pD->NGrid[0]].ID_Comm_Domain, myN+myM*pD->NGrid[0], MPI_COMM_WORLD, &req3[np]);
      }
//receive data
    for (np=0; np<=pD->NGrid[0]*pD->NGrid[1]-1; np++){
        MPI_Irecv(&T_cut3[0][0], pG->Nx[1]/pD->NGrid[0]*(pG->Nx[0]/pD->NGrid[1]+2), MPI_DOUBLE, myID_Comm_world, np, MPI_COMM_WORLD, &req4[np]);
            co=0;
            for (j=js+np*pG->Nx[1]/pD->NGrid[0]; j<=js+(np+1)*pG->Nx[1]/pD->NGrid[0]-1; j++ ){
                co1=0;
                for (i=is-1; i<=is+pG->Nx[0]/pD->NGrid[1]; i++){
                    Temp2[j][i]=T_cut3[co][co1++];
                    if (ts==0){
                        Temp22[0][j][i]=Temp2[j][i];
                        Temp22[1][j][i]=Temp2[j][i];
                    }
                        
                }
                co++;
            }

     }
 //tridiagonal matrix
 //thomas algorithm in y-direction
      for (k=kl; k<=ku; k++){
          for (i=is-1; i<=is+pG->Nx[0]/pD->NGrid[1]; i++){
              co=0;
              for (j=js; j<=js+pD->NGrid[1]*pG->Nx[1]-1; j++){
                  B01 = SQR(pG->B2i[k][j][i])+SQR(0.5*(pG->U[k][j+1][i].B1c + pG->U[k][j][i].B1c));
                  B01 = MAX(B01,TINY_NUMBER); /* limit in case B=0 */
                  B02 = SQR(pG->B2i[k][j][i])+SQR(0.5*(pG->U[k][j-1][i].B1c + pG->U[k][j][i].B1c));
                  B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */
                  if (j==js){
                      dTdx2 = limiter2(Temp2[j  ][i+1] - Temp2[j  ][i ],
                                       Temp2[j  ][i  ] - Temp2[j  ][i ]);
                      a2[co]=0;
                      
                  }
                  else{
                      dTdx2 = limiter4(Temp2[j  ][i+1] - Temp2[j  ][i  ],
                                       Temp2[j  ][i  ] - Temp2[j  ][i-1],
                                       Temp2[j-1][i-1] - Temp2[j-1][i  ],
                                       Temp2[j-1][i  ] - Temp2[j-1][i-1]);
                      a2[co]=-Lambda*Gamma_1/pG->U[k][j][i].d*pG->B2i[k][j][i]*pG->B2i[k][j][i]/B02*k_nl((T_ex(Temp22,Nt,j-1,i)+T_ex(Temp22,Nt,j,i))/2.)*my_dt/dx2/dx2;
                  }
                  if (j==js+pD->NGrid[1]*pG->Nx[1]-1){
                      dTdx1 = limiter2(Temp2[j  ][i+1] - Temp2[j  ][i  ],
                                       Temp2[j  ][i  ] - Temp2[j  ][i-1]);
                      c1[co]=0;
                  }
                  else{
                      dTdx1 = limiter4(Temp2[j+1][i+1] - Temp2[j+1][i  ],
                                       Temp2[j+1][i  ] - Temp2[j+1][i-1],
                                       Temp2[j  ][i+1] - Temp2[j  ][i  ],
                                       Temp2[j  ][i  ] - Temp2[j  ][i-1]);
                      c2[co]=-Lambda*Gamma_1/pG->U[k][j][i].d*(pG->B2i[k][j+1][i]*pG->B2i[k][j+1][i]/B01*k_nl((T_ex(Temp22,Nt,j+1,i)+T_ex(Temp22,Nt,j,i))/2.)*my_dt/dx2/dx2);
                  }
                  
                  x2[co]=Temp2[j][i]+Lambda*Gamma_1/pG->U[k][j][i].d*(pG->B2i[k][j+1][i]*pG->B1i[k][j+1][i]/B01*k_nl((T_ex(Temp22,Nt,j+1,i)+T_ex(Temp22,Nt,j,i))/2.)*dTdx1*my_dt/dx1/dx2-pG->B2i[k][j][i]*pG->B1i[k][j][i]/B02*k_nl((T_ex(Temp22,Nt,j-1,i)+T_ex(Temp22,Nt,j,i))/2.)*dTdx2*my_dt/dx1/dx2);
            if (j==js)
            b2[co]=1+Lambda*Gamma_1/pG->U[k][j][i].d*pG->B1i[k][j+1][i]*pG->B1i[k][j+1][i]/B01*k_nl((T_ex(Temp22,Nt,j+1,i)+T_ex(Temp22,Nt,j,i))/2.)*my_dt/dx2/dx2;
            else if (j==js+pD->NGrid[1]*pG->Nx[1]-1)
            b2[co]=1+Lambda*Gamma_1/pG->U[k][j][i].d*pG->B2i[k][j][i]*pG->B2i[k][j][i]/B02*k_nl((T_ex(Temp22,Nt,j-1,i)+T_ex(Temp22,Nt,j,i))/2.)*my_dt/dx2/dx2;
            else
            b2[co]=1+Lambda*Gamma_1/pG->U[k][j][i].d*(pG->B2i[k][j][i]*pG->B2i[k][j][i]/B02*k_nl((T_ex(Temp22,Nt,j-1,i)+T_ex(Temp22,Nt,j,i))/2.)*my_dt/dx2/dx2+pG->B1i[k][j+1][i]*pG->B1i[k][j+1][i]/B01*k_nl((T_ex(Temp22,Nt,j+1,i)+T_ex(Temp22,Nt,j,i))/2.)*my_dt/dx2/dx2);
             co++;
                                                                             }
              
              tdma(x2,pG->Nx[1]*pD->NGrid[1],a2,b2,c2);
                co1=0;
              for (j=js; j<=js+pD->NGrid[1]*pG->Nx[1]-1; j++){
                  Temp2[j][i]=x2[co1++];
                  Temp22[Nt][j][i]=Temp2[j][i];
              }
              }
          }
/*ADI  done */
//calculate the flux to add on the energy density in the initial domain

for (np=0; np<=pD->NGrid[1]-1; np++){

             co=0;
             for (j=js+np*pG->Nx[1]; j<=js+(np+1)*pG->Nx[1]; j++){
                 co1=0;
                 for (i=is; i<=pG->Nx[0]/pD->NGrid[1]-1; i++){
                     T_cut4[co][co1++]=Temp2[j][i];
                 }
                 co++;
             }

         MPI_Isend(&T_cut4[0][0], pG->Nx[1]*pG->Nx[0]/pD->NGrid[1], MPI_DOUBLE, pD->GData[myL][np][myN/pD->NGrid[1]].ID_Comm_Domain, myM, MPI_COMM_WORLD, &req3[np]);
          }
      for (np=0; np<=pD->NGrid[1]-1; np++){
          MPI_Irecv(&T_cut5[0][0], pG->Nx[1]*pG->Nx[0]/pD->NGrid[1], MPI_DOUBLE, myID_Comm_world, np, MPI_COMM_WORLD, &req4[np]);
          for (k=kl; k<=ku; k++){
              co=0;
              for (j=js; j<=je; j++){
                  co1=0;
                  for (i=is+np*pG->Nx[0]/pD->NGrid[1]; i<=is+(np+1)*pG->Nx[0]/pD->NGrid[1]-1; i++){
                      Temp_[k][j][i]=T_cut5[co][co1++];
                      pG->U[k][j][i].E=pG->U[k][j][i].d/Lambda/Gamma_1*Temp_[k][j][i]+(0.5/pG->U[k][j][i].d)*
                      (SQR(pG->U[k][j][i].M1) +SQR(pG->U[k][j][i].M2) +SQR(pG->U[k][j][i].M3))+(0.5)*(SQR(pG->U[k][j][i].B1c)+SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));;
                  }
                  co++;
              }
          }
       }
    //MPI arrangement
    //send data
    for (np=0; np<=pD->NGrid[0]-1; np++ ){
        for (k=kl; k<=ku; k++){
            co=0;
            for (j=jl+np*pG->Nx[1]/pD->NGrid[0]; j<=jl+(np+1)*pG->Nx[1]/pD->NGrid[0]+1; j++){
                co1=0;
                for (i=is; i<=ie; i++){
                    T_cut[co][co1++]=Temp_[k][j][i];
                }
                co++;
            }
        MPI_Isend(&T_cut[0][0], (pG->Nx[1]/pD->NGrid[0]+2)*pG->Nx[0], MPI_DOUBLE, pD->GData[myL][myM][np].ID_Comm_Domain, myN, MPI_COMM_WORLD, &req1[np]);
    }
    }
    
    //receive data
    for (np=0; np<=pD->NGrid[0]-1; np++ ){
        MPI_Irecv(&T_cut1[0][0], (pG->Nx[1]/pD->NGrid[0]+2)*pG->Nx[0], MPI_DOUBLE, myID_Comm_world, np, MPI_COMM_WORLD, &req2[np]);
        co=0;
        for (j=jl; j<=jl+pG->Nx[1]/pD->NGrid[0]+1; j++){
            co1=0;
            for (i=is; i<=ie; i++){
                Temp11[Nt][j][i+np*pG->Nx[0]]=T_cut1[co][co1++];
            }
            co++;
        }
        
        
    }
    ++ts;
}
/*----------------------------------------------------------------------------*/
/* limiter2 and limiter4: call slope limiters to preserve monotonicity
 */

static Real limiter2(const Real A, const Real B)
{
    /* van Leer slope limiter */
    return vanleer(A,B);
    
    /* monotonized central (MC) limiter */
    /* return minmod(2.0*minmod(A,B),0.5*(A+B)); */
}

static Real limiter4(const Real A, const Real B, const Real C, const Real D)
{
    return limiter2(limiter2(A,B),limiter2(C,D));
}

/*----------------------------------------------------------------------------*/
/* vanleer: van Leer slope limiter
 */

static Real vanleer(const Real A, const Real B)
{
    if (A*B > 0) {
        return 2.0*A*B/(A+B);
    } else {
        return 0.0;
    }
}

/*----------------------------------------------------------------------------*/
/* minmod: minmod slope limiter
 */

static Real minmod(const Real A, const Real B)
{
    if (A*B > 0) {
        if (A > 0) {
            return MIN(A,B);
        } else {
            return MAX(A,B);
        }
    } else {
        return 0.0;
    }
}
/*---------------------------------------------------------------------------*/
/* tdma: thomas algorithm for tridiagonal matrix
 */
void tdma(Real x[], const int N, const Real a[], const Real b[], Real c[])
{
    int n;
    
    c[0] = c[0] / b[0];
    x[0] = x[0] / b[0];
    
    for (n = 1; n < N; n++) {
        Real m = 1.0f / (b[n] - a[n] * c[n - 1]);
        c[n] = c[n] * m;
        x[n] = (x[n] - a[n] * x[n - 1]) * m;
    }
    
    for (n = N - 1; n-- > 0; )
        x[n] = x[n] - c[n] * x[n + 1];
}
/*spitzer model*/
static Real k_nl(const Real T)
{
    Real kappa;
    kappa = 1.84e-10/30.*T*T*sqrt(T);
    return kappa;
}

static Real T_ex(Real ***T, const int Nt, const int j, const int i)
{
    Real T_ap;
    if (Nt == 0)
        T_ap = 2*T[0][j][i]-T[1][j][i];
    else
        T_ap = 2*T[1][j][i]-T[0][j][i];
    return T_ap;
}

/*----------------------------------------------------------------------------*/
/*! \fn void conduction_init(MeshS *pM)
 *  \brief Allocate temporary arrays
 */

void conduction_init(MeshS *pM)
{
    DomainS **pD = pM->Domain;
    int nl,nd,size1=1,size2=1,size3=1,Nx1,Nx2,Nx3;
    int ng1,ng2,ng3,size11,size22,size33,Nxx1,Nxx2,Nxx3;
    
    
    /* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
    for (nl=0; nl<(pM->NLevels); nl++){
        for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
            if (pM->Domain[nl][nd].Grid != NULL) {
                if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
                    ng1 = pD[nl][nd].NGrid[0];
                    size1 = pM->Domain[nl][nd].Grid->Nx[0];
                    size11 = ng1*size1;
                }
                if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
                    ng2 = pD[nl][nd].NGrid[1];
                    size2 = pM->Domain[nl][nd].Grid->Nx[1];
                    size22 = ng2*size2;
                }
                if (pM->Domain[nl][nd].Grid->Nx[2] > size3){
                    ng3 = pD[nl][nd].NGrid[1];
                    size3 = pM->Domain[nl][nd].Grid->Nx[2];
                    size33 = pD[nl][nd].NGrid[2]*size3;
                }
            }
        }
    }
    
    Nx1 = size1 + 2*nghost;
    Nxx1 = size11+2*nghost;
    
    if (pM->Nx[1] > 1){
        Nx2 = size2 + 2*nghost;
        Nxx2 = size22+2*nghost;
    } else {
        Nx2 = size2;
        Nxx2 = size22;
    }
    
    if (pM->Nx[2] > 1){
        Nx3 = size3 + 2*nghost;
        Nxx3 = size33+2*nghost;
    } else {
        Nx3 = size3;
        Nxx3 = size33;
    }
    if ((Temp = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
        goto on_error;
    if ((Temp_ = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
        goto on_error;
    if ((Temp1 = (Real**)calloc_2d_array(size2/ng1+2*nghost,Nxx1, sizeof(Real))) == NULL)
        goto on_error;
    if ((Temp2 = (Real**)calloc_2d_array(Nxx2,size1/ng2+2*nghost, sizeof(Real))) == NULL)
        goto on_error;
    if ((Temp11 = (Real***)calloc_3d_array(2,size2/ng1+2*nghost,Nxx1, sizeof(Real))) == NULL)
        goto on_error;
    if ((Temp22 = (Real***)calloc_3d_array(2,Nxx2,size1/ng2+2*nghost, sizeof(Real))) == NULL)
        goto on_error;

    return;
    
on_error:
    conduction_destruct();
    ath_error("[conduct_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void conduction_destruct(void)
 *  \brief Free temporary arrays
 */

void conduction_destruct(void)
{
    if (Temp != NULL) free_3d_array(Temp);
    if (Temp_ != NULL) free_3d_array(Temp_);
    if (Temp11 != NULL) free_3d_array(Temp11);
    if (Temp22 != NULL) free_3d_array(Temp22);
    if (Temp1 != NULL) free_2d_array(Temp1);
    if (Temp2 != NULL) free_2d_array(Temp2);
    
    return;
}
#endif
