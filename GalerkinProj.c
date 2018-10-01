/*
 * This function finds the mode combinations which lead to non-zero coupling constant
 */

#include "struct_def.h"

double GalerkinProj(gsl_vector* A, gsl_vector* B, gsl_vector* C, disc_const dc){
  int i, j;
  double H, delxxB, delxxC, delyyB, delyyC, delxpypB, delxpypC, delxmymB, delxmymC, delxpymB, delxpymC, delxmypB, delxmypC, LijBC ;
  
  H = 0.0;
  for (j = 1; j < dc.Nfy-1; j++){
    for (i = 1; i < dc.Nfx-1; i++){
      delxxB      = 1.0/dc.hx/dc.hx*( gsl_vector_get(B, (i+1)+(j)*dc.Nfx) + gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxxC      = 1.0/dc.hx/dc.hx*( gsl_vector_get(C, (i+1)+(j)*dc.Nfx) + gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      
      delyyC      = 1.0/dc.hy/dc.hy*( gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j-1)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      delyyB      = 1.0/dc.hy/dc.hy*( gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j-1)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      
      delxpypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i+1)+(j+1)*dc.Nfx) - gsl_vector_get(B, (i+1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i+1)+(j+1)*dc.Nfx) - gsl_vector_get(C, (i+1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i-1)+(j-1)*dc.Nfx) - gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i-1)+(j-1)*dc.Nfx) - gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxpymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i+1)+(j-1)*dc.Nfx) - gsl_vector_get(B, (i+1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i+1)+(j-1)*dc.Nfx) - gsl_vector_get(C, (i+1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i-1)+(j+1)*dc.Nfx) - gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i-1)+(j+1)*dc.Nfx) - gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      LijBC       = delxxB*delyyC + delyyB*delxxC + -1.0/2.0*( delxpypB*delxpypC + delxmymB*delxmymC + delxmypB*delxmypC + delxpymB*delxpymC);
      H 	  = H + dc.hx*dc.hy*gsl_vector_get(A, (i)+(j)*dc.Nfx)*LijBC;

    }
  }

  for (j = 0; j < 1; j++){
    for (i = 0; i < 1; i++){
      delxxB      = 1.0/dc.hx/dc.hx*( gsl_vector_get(B, (i+1)+(j)*dc.Nfx) + gsl_vector_get(B, (0)+(j)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxxC      = 1.0/dc.hx/dc.hx*( gsl_vector_get(C, (i+1)+(j)*dc.Nfx) + gsl_vector_get(C, (0)+(j)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      
      delyyC      = 1.0/dc.hy/dc.hy*( gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(0)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      delyyB      = 1.0/dc.hy/dc.hy*( gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(0)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      
      delxpypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i+1)+(j+1)*dc.Nfx) - gsl_vector_get(B, (i+1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i+1)+(j+1)*dc.Nfx) - gsl_vector_get(C, (i+1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (0)+(0)*dc.Nfx) - gsl_vector_get(B, (0)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(0)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (0)+(0)*dc.Nfx) - gsl_vector_get(C, (0)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(0)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxpymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i+1)+(0)*dc.Nfx) - gsl_vector_get(B, (i+1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(0)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i+1)+(0)*dc.Nfx) - gsl_vector_get(C, (i+1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(0)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (0)+(j+1)*dc.Nfx) - gsl_vector_get(B, (0)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (0)+(j+1)*dc.Nfx) - gsl_vector_get(C, (0)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      LijBC       = delxxB*delyyC + delyyB*delxxC + -1.0/2.0*( delxpypB*delxpypC + delxmymB*delxmymC + delxmypB*delxmypC + delxpymB*delxpymC);
      H 	  = H + dc.hx*dc.hy*gsl_vector_get(A, (i)+(j)*dc.Nfx)*LijBC;

    }
  }

  for (j = 0; j < 1; j++){
    for (i = 1; i < dc.Nfx-1; i++){
      delxxB      = 1.0/dc.hx/dc.hx*( gsl_vector_get(B, (i+1)+(j)*dc.Nfx) + gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxxC      = 1.0/dc.hx/dc.hx*( gsl_vector_get(C, (i+1)+(j)*dc.Nfx) + gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      
      delyyC      = 1.0/dc.hy/dc.hy*( gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(0)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      delyyB      = 1.0/dc.hy/dc.hy*( gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(0)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      
      delxpypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i+1)+(j+1)*dc.Nfx) - gsl_vector_get(B, (i+1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i+1)+(j+1)*dc.Nfx) - gsl_vector_get(C, (i+1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i-1)+(0)*dc.Nfx) - gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(0)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i-1)+(0)*dc.Nfx) - gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(0)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxpymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i+1)+(0)*dc.Nfx) - gsl_vector_get(B, (i+1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(0)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i+1)+(0)*dc.Nfx) - gsl_vector_get(C, (i+1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(0)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i-1)+(j+1)*dc.Nfx) - gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i-1)+(j+1)*dc.Nfx) - gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      LijBC       = delxxB*delyyC + delyyB*delxxC + -1.0/2.0*( delxpypB*delxpypC + delxmymB*delxmymC + delxmypB*delxmypC + delxpymB*delxpymC);
      H 	  = H + dc.hx*dc.hy*gsl_vector_get(A, (i)+(j)*dc.Nfx)*LijBC;

    }
  }

  for (j = 0; j < 1; j++){
    for (i = dc.Nfx-1; i < dc.Nfx; i++){
      delxxB      = 1.0/dc.hx/dc.hx*( gsl_vector_get(B, (dc.Nfx-1)+(j)*dc.Nfx) + gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxxC      = 1.0/dc.hx/dc.hx*( gsl_vector_get(C, (dc.Nfx-1)+(j)*dc.Nfx) + gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      
      delyyC      = 1.0/dc.hy/dc.hy*( gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(0)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      delyyB      = 1.0/dc.hy/dc.hy*( gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(0)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      
      delxpypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (dc.Nfx-1)+(j+1)*dc.Nfx) - gsl_vector_get(B, (dc.Nfx-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (dc.Nfx-1)+(j+1)*dc.Nfx) - gsl_vector_get(C, (dc.Nfx-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i-1)+(0)*dc.Nfx) - gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(0)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i-1)+(0)*dc.Nfx) - gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(0)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxpymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (dc.Nfx-1)+(0)*dc.Nfx) - gsl_vector_get(B, (dc.Nfx-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(0)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (dc.Nfx-1)+(0)*dc.Nfx) - gsl_vector_get(C, (dc.Nfx-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(0)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i-1)+(j+1)*dc.Nfx) - gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i-1)+(j+1)*dc.Nfx) - gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      LijBC       = delxxB*delyyC + delyyB*delxxC + -1.0/2.0*( delxpypB*delxpypC + delxmymB*delxmymC + delxmypB*delxmypC + delxpymB*delxpymC);
      H 	  = H + dc.hx*dc.hy*gsl_vector_get(A, (i)+(j)*dc.Nfx)*LijBC;

    }
  }

  for (j = 1; j < dc.Nfy-1; j++){
    for (i = 0; i < 1; i++){
      delxxB      = 1.0/dc.hx/dc.hx*( gsl_vector_get(B, (i+1)+(j)*dc.Nfx) + gsl_vector_get(B, (0)+(j)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxxC      = 1.0/dc.hx/dc.hx*( gsl_vector_get(C, (i+1)+(j)*dc.Nfx) + gsl_vector_get(C, (0)+(j)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      
      delyyC      = 1.0/dc.hy/dc.hy*( gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j-1)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      delyyB      = 1.0/dc.hy/dc.hy*( gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j-1)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      
      delxpypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i+1)+(j+1)*dc.Nfx) - gsl_vector_get(B, (i+1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i+1)+(j+1)*dc.Nfx) - gsl_vector_get(C, (i+1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (0)+(j-1)*dc.Nfx) - gsl_vector_get(B, (0)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (0)+(j-1)*dc.Nfx) - gsl_vector_get(C, (0)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxpymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i+1)+(j-1)*dc.Nfx) - gsl_vector_get(B, (i+1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i+1)+(j-1)*dc.Nfx) - gsl_vector_get(C, (i+1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (0)+(j+1)*dc.Nfx) - gsl_vector_get(B, (0)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (0)+(j+1)*dc.Nfx) - gsl_vector_get(C, (0)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      LijBC       = delxxB*delyyC + delyyB*delxxC + -1.0/2.0*( delxpypB*delxpypC + delxmymB*delxmymC + delxmypB*delxmypC + delxpymB*delxpymC);
      H 	  = H + dc.hx*dc.hy*gsl_vector_get(A, (i)+(j)*dc.Nfx)*LijBC;

    }
  }

  for (j = 1; j < dc.Nfy-1; j++){
    for (i = dc.Nfx-1; i < dc.Nfx; i++){
      delxxB      = 1.0/dc.hx/dc.hx*( gsl_vector_get(B, (dc.Nfx-1)+(j)*dc.Nfx) + gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxxC      = 1.0/dc.hx/dc.hx*( gsl_vector_get(C, (dc.Nfx-1)+(j)*dc.Nfx) + gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      
      delyyC      = 1.0/dc.hy/dc.hy*( gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j-1)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      delyyB      = 1.0/dc.hy/dc.hy*( gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j-1)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      
      delxpypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (dc.Nfx-1)+(j+1)*dc.Nfx) - gsl_vector_get(B, (dc.Nfx-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (dc.Nfx-1)+(j+1)*dc.Nfx) - gsl_vector_get(C, (dc.Nfx-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i-1)+(j-1)*dc.Nfx) - gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i-1)+(j-1)*dc.Nfx) - gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxpymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (dc.Nfx-1)+(j-1)*dc.Nfx) - gsl_vector_get(B, (dc.Nfx-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (dc.Nfx-1)+(j-1)*dc.Nfx) - gsl_vector_get(C, (dc.Nfx-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i-1)+(j+1)*dc.Nfx) - gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j+1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i-1)+(j+1)*dc.Nfx) - gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j+1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      LijBC       = delxxB*delyyC + delyyB*delxxC + -1.0/2.0*( delxpypB*delxpypC + delxmymB*delxmymC + delxmypB*delxmypC + delxpymB*delxpymC);
      H 	  = H + dc.hx*dc.hy*gsl_vector_get(A, (i)+(j)*dc.Nfx)*LijBC;

    }
  }

  for (j = dc.Nfy-1; j < dc.Nfy; j++){
    for (i = 0; i < 1; i++){
      delxxB      = 1.0/dc.hx/dc.hx*( gsl_vector_get(B, (i+1)+(j)*dc.Nfx) + gsl_vector_get(B, (0)+(j)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxxC      = 1.0/dc.hx/dc.hx*( gsl_vector_get(C, (i+1)+(j)*dc.Nfx) + gsl_vector_get(C, (0)+(j)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      
      delyyC      = 1.0/dc.hy/dc.hy*( gsl_vector_get(C, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j-1)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      delyyB      = 1.0/dc.hy/dc.hy*( gsl_vector_get(B, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j-1)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      
      delxpypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i+1)+(dc.Nfy-1)*dc.Nfx) - gsl_vector_get(B, (i+1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i+1)+(dc.Nfy-1)*dc.Nfx) - gsl_vector_get(C, (i+1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (0)+(j-1)*dc.Nfx) - gsl_vector_get(B, (0)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (0)+(j-1)*dc.Nfx) - gsl_vector_get(C, (0)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxpymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i+1)+(j-1)*dc.Nfx) - gsl_vector_get(B, (i+1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i+1)+(j-1)*dc.Nfx) - gsl_vector_get(C, (i+1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (0)+(dc.Nfy-1)*dc.Nfx) - gsl_vector_get(B, (0)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (0)+(dc.Nfy-1)*dc.Nfx) - gsl_vector_get(C, (0)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      LijBC       = delxxB*delyyC + delyyB*delxxC + -1.0/2.0*( delxpypB*delxpypC + delxmymB*delxmymC + delxmypB*delxmypC + delxpymB*delxpymC);
      H 	  = H + dc.hx*dc.hy*gsl_vector_get(A, (i)+(j)*dc.Nfx)*LijBC;

    }
  }
  
  for (j = dc.Nfy-1; j < dc.Nfy; j++){
    for (i = 1; i < dc.Nfx-1; i++){
      delxxB      = 1.0/dc.hx/dc.hx*( gsl_vector_get(B, (i+1)+(j)*dc.Nfx) + gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxxC      = 1.0/dc.hx/dc.hx*( gsl_vector_get(C, (i+1)+(j)*dc.Nfx) + gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      
      delyyC      = 1.0/dc.hy/dc.hy*( gsl_vector_get(C, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j-1)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      delyyB      = 1.0/dc.hy/dc.hy*( gsl_vector_get(B, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j-1)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      
      delxpypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i+1)+(dc.Nfy-1)*dc.Nfx) - gsl_vector_get(B, (i+1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i+1)+(dc.Nfy-1)*dc.Nfx) - gsl_vector_get(C, (i+1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i-1)+(j-1)*dc.Nfx) - gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i-1)+(j-1)*dc.Nfx) - gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxpymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i+1)+(j-1)*dc.Nfx) - gsl_vector_get(B, (i+1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i+1)+(j-1)*dc.Nfx) - gsl_vector_get(C, (i+1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i-1)+(dc.Nfy-1)*dc.Nfx) - gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i-1)+(dc.Nfy-1)*dc.Nfx) - gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      LijBC       = delxxB*delyyC + delyyB*delxxC + -1.0/2.0*( delxpypB*delxpypC + delxmymB*delxmymC + delxmypB*delxmypC + delxpymB*delxpymC);
      H 	  = H + dc.hx*dc.hy*gsl_vector_get(A, (i)+(j)*dc.Nfx)*LijBC;

    }
  }
  
  for (j = dc.Nfy-1; j < dc.Nfy; j++){
    for (i = dc.Nfx-1; i < dc.Nfx; i++){
      delxxB      = 1.0/dc.hx/dc.hx*( gsl_vector_get(B, (dc.Nfx-1)+(j)*dc.Nfx) + gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxxC      = 1.0/dc.hx/dc.hx*( gsl_vector_get(C, (dc.Nfx-1)+(j)*dc.Nfx) + gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      
      delyyC      = 1.0/dc.hy/dc.hy*( gsl_vector_get(C, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j-1)*dc.Nfx) - 2*gsl_vector_get(C, (i)+(j)*dc.Nfx) );  
      delyyB      = 1.0/dc.hy/dc.hy*( gsl_vector_get(B, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j-1)*dc.Nfx) - 2*gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      
      delxpypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (dc.Nfx-1)+(dc.Nfy-1)*dc.Nfx) - gsl_vector_get(B, (dc.Nfx-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (dc.Nfx-1)+(dc.Nfy-1)*dc.Nfx) - gsl_vector_get(C, (dc.Nfx-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i-1)+(j-1)*dc.Nfx) - gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i-1)+(j-1)*dc.Nfx) - gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxpymB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (dc.Nfx-1)+(j-1)*dc.Nfx) - gsl_vector_get(B, (dc.Nfx-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(j-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxpymC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (dc.Nfx-1)+(j-1)*dc.Nfx) - gsl_vector_get(C, (dc.Nfx-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(j-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      delxmypB    = 1.0/dc.hx/dc.hy*( gsl_vector_get(B, (i-1)+(dc.Nfy-1)*dc.Nfx) - gsl_vector_get(B, (i-1)+(j)*dc.Nfx) - gsl_vector_get(B, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(B, (i)+(j)*dc.Nfx) );  
      delxmypC    = 1.0/dc.hx/dc.hy*( gsl_vector_get(C, (i-1)+(dc.Nfy-1)*dc.Nfx) - gsl_vector_get(C, (i-1)+(j)*dc.Nfx) - gsl_vector_get(C, (i)+(dc.Nfy-1)*dc.Nfx) + gsl_vector_get(C, (i)+(j)*dc.Nfx) );
      
      LijBC       = delxxB*delyyC + delyyB*delxxC + -1.0/2.0*( delxpypB*delxpypC + delxmymB*delxmymC + delxmypB*delxmypC + delxpymB*delxpymC);
      H 	  = H + dc.hx*dc.hy*gsl_vector_get(A, (i)+(j)*dc.Nfx)*LijBC;

    }
  }

  return H;
}
// Matlab code
// delxxB      = 1/dc.hx/dc.hx*( Bmat(i+1, j) + Bmat(i-1, j) - 2*Bmat(i, j) );
// delxxC      = 1/dc.hx/dc.hx*( Cmat(i+1, j) + Cmat(i-1, j) - 2*Cmat(i, j) );

// delyyC      = 1/dc.hy/dc.hy*( Cmat(i, j+1) + Cmat(i, j-1) - 2*Cmat(i, j) );
// delyyB      = 1/dc.hy/dc.hy*( Bmat(i, j+1) + Bmat(i, j-1) - 2*Bmat(i, j) );
        
// delxpypB    = 1/dc.hx/dc.hy*(  Bmat(i+1, j+1) - Bmat(i+1, j) - Bmat(i, j+1) + Bmat(i, j) );
// delxpypC    = 1/dc.hx/dc.hy*(  Cmat(i+1, j+1) - Cmat(i+1, j) - Cmat(i, j+1) + Cmat(i, j) );
        
// delxmymB    = 1/dc.hx/dc.hy*(  Bmat(i-1, j-1) - Bmat(i-1, j) - Bmat(i, j-1) + Bmat(i, j) );
// delxmymC    = 1/dc.hx/dc.hy*(  Cmat(i-1, j-1) - Cmat(i-1, j) - Cmat(i, j-1) + Cmat(i, j) );

// delxpymB    = 1/dc.hx/dc.hy*(  Bmat(i+1, j-1) - Bmat(i+1, j) - Bmat(i, j-1) + Bmat(i, j) );
// delxpymC    = 1/dc.hx/dc.hy*(  Cmat(i+1, j-1) - Cmat(i+1, j) - Cmat(i, j-1) + Cmat(i, j) );

// delxmypB    = 1/dc.hx/dc.hy*(  Bmat(i-1, j+1) - Bmat(i-1, j) - Bmat(i, j+1) + Bmat(i, j) );
// delxmypC    = 1/dc.hx/dc.hy*(  Cmat(i-1, j+1) - Cmat(i-1, j) - Cmat(i, j+1) + Cmat(i, j) );
 
