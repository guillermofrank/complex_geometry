#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

// Social Force Model

#define N         1          // cantidad de particulas (1)
#define M         70.0       // masa de las particulas (70)
#define TAU       0.50       // tiempo de reaction de peaton (0.5)
#define SAMPLE    0.05       // intervalo de sampleo (0.05)
#define A         2000.0     // intensidad fuerza social
#define B         0.08       // escala espacial de fuerza social
#define KAPPA     25000.0    // coeficiente de friccion
#define RADIUS    0.2        // radio de los peatones
#define VD        1.0        // velocidad deseada (modulo)
#define STEP      0.0001     // paso de integracion tipico 0.0001

// Cell list

#define EMPTY     -1         // empty value -1 (no not change)
#define RCELL     1.2        // lado de las celdas cuadradas 
#define PEDTAB    5000       // table length for social force between pedestrians
#define PEDRC     1.0        // cut-off distance between pedestrians (5*RADIUS)      
#define PEDRI     0.04       // minimum distance between pedestrians (RADIUS/5)
#define OBSTAB    5000       // table length for social force with obstacles
#define OBSRC     1.0        // cut-off distance for obstacles (5*RADIUS)      
#define OBSRI     0.02       // minimum distance for obstacles (RADIUS/5)
#define BIG       10000000.0 // any big number

// Plot

#define WIDTH     6.5         // .eps width 
#define HEIGHT    4.017       // .eps height 0.618 x 6.5 = 4.017
#define UNITS     72          // .eps units (in 'pt' = 1/72 inch)
#define OFFSET    0.2         // .eps offset for the bounding box


double  table_ped[4*PEDTAB],table_obs[4*OBSTAB];

void    checkinput(char **argv,int argc);
int     size_of_myfile(char *myfile);
int     size_of_myrule(char *myrule);
int     size_of_cells(double *box,double rcell);
void    read_initial(double *x,double *v,double *f,int n);
void    read_myfile(char *myfile,double *xb);
void    read_myrule(char *myrule,int *ti,double *xp,double *tp,double *box,int np);
void    build_table_ped(double ri,double rf,int ntable);
void    build_table_obs(double ri,double rf,int ntable);
int     buildlist(double *xp,double *x,int *head,int *lscl,double *box,int np,int n,int nc);
int     computelist(double *xb,double *xp,double *x,double *v,double *f,int *head,int *lscl,double *box,int np,int n,int nc,int nb);
void    force_ped(double *xb,double *x,double *v,double *f,double *shift,int i,int j,int np,int nb);
void    force_obs(double *xb,double *x,double *v,double *f,int n,int nb);
void    force_vd(double *v,double *v0,double *f,int n);
int     intersection(double *xb,double *p1,double *p2,int nb);
void    collision(double *x,double *v,double nx,double ny,int ii,int jj);
void    positions(double *x,double *v,double *f,double *f0,double *box,int n);
void    velocities(double *v,double *f,double *f0,int n);
void    desire_test(int *mem,int *ti,double *xp,double *tt,double *xb,double *x,double *v0,int *nm,int np,int n,int nb);
int     set_invisibles(int *mem,double *tt,double *xb,double *p1,int *nm,int n,int nb,int i);
int     neighbor(int *ti,double *xp,double *p,int np);
void    move2target(double *v0,double *tt,double *x,double *xb,int i,int j,int np,int nb);
int     targets(int *mem,int *ti,double *tp,double *tt,double *x,double *xb,int *nm,int np,int n,int nb);
double  open_plot(char *myplot,char *myfile,double *xb,double *box,int i,int nb);
void    itoa(int n,char s[]);
double  boundaries(FILE *fp,double *xb,double scale,int nb);
double  plot_people(char *myplot,double *v0,double *x,double scale,int n);
void    close_plot(char *myplot);


int main(int argc, char *argv[])
{
  char    myfile[100],myrule[100],myplot[100];
  int     i,n,nb,np,nc,*head,*lscl,*ti,*mem,*nm;
  double  t,scale;
  double  *x,*v,*f,*v0,*f0,*xb,*xp,*tp,*tt,*box;
  clock_t c0,c1;

  c0=clock();
  
  n = N;

  checkinput(argv,argc);

  strcpy(myfile,argv[1]);
  strcpy(myrule,argv[2]);

  nb = size_of_myfile(myfile);
  np = size_of_myrule(myrule);

  x  = (double *)malloc(2*n*sizeof(double));
  v  = (double *)malloc(2*n*sizeof(double));
  f  = (double *)malloc(2*n*sizeof(double));
  v0 = (double *)malloc(2*n*sizeof(double));
  f0 = (double *)malloc(2*n*sizeof(double));
  xb = (double *)malloc(4*nb*sizeof(double));
  xp = (double *)malloc(2*np*sizeof(double));
  tp = (double *)malloc(2*np*sizeof(double));
  tt = (double *)malloc(2*np*sizeof(double));
  ti = (int *)malloc(np*sizeof(int));
  box = (double *)malloc(6*sizeof(double));
  mem = (int *)malloc(n*np*sizeof(int));
  nm  = (int *)malloc(n*sizeof(int));

  read_initial(x,v,f,n);
  read_myfile(myfile,xb);
  read_myrule(myrule,ti,xp,tp,box,np);

  build_table_ped(PEDRI,PEDRC,PEDTAB);
  build_table_obs(OBSRI,OBSRC,OBSTAB);

  nc = size_of_cells(box,RCELL);

  head = (int *)malloc(nc*nc*sizeof(int));
  lscl = (int *)malloc((n+np)*sizeof(int));

  buildlist(xp,x,head,lscl,box,np,n,nc); 
  computelist(xb,xp,x,v,f,head,lscl,box,np,n,nc,nb);

  targets(mem,ti,tp,tt,x,xb,nm,np,n,nb);

  desire_test(mem,ti,xp,tt,xb,x,v0,nm,np,n,nb);

  for(i=0;i<300;i++)
    {
      for(t=0;t<SAMPLE;t=t+STEP)
        {
          positions(x,v,f,f0,box,n);
          buildlist(xp,x,head,lscl,box,np,n,nc);
          computelist(xb,xp,x,v,f,head,lscl,box,np,n,nc,nb);
          force_vd(v,v0,f,n);
          velocities(v,f,f0,n);
        }

      desire_test(mem,ti,xp,tt,xb,x,v0,nm,np,n,nb);
      printf("%d x=%lf y=%lf vdx=%lf vdy=%lf vx=%lf vy=%lf fx=%lf fy=%lf\n",i,*(x+0),*(x+1),*(v0+0),*(v0+1),*(v+0),*(v+1),*(f+0),*(f+1));

      if (i%1==0)
        {
          scale = open_plot(myplot,myfile,xb,box,i,nb);
          scale = plot_people(myplot,v0,x,scale,n); 
          close_plot(myplot);
        }
    }

  free(head);
  free(lscl);

  free(x);
  free(v);
  free(f);
  free(v0);
  free(f0);
  free(xb);
  free(xp);
  free(tp);
  free(tt);
  free(ti);
  free(box);
  free(mem);
  free(nm);

  c1=clock();

  printf("\nCPU time (sec.): %f\n\n",(float)(c1-c0)/CLOCKS_PER_SEC);

  return 1;
}


/*************************************************************/
/*      AUXILIARY FUNCTIONS: DECISION MAKING                 */
/*************************************************************/

void desire_test(int *mem,int *ti,double *xp,double *tt,double *xb,double *x,double *v0,int *nm,int np,int n,int nb)
{
   int    h,i,j;
   double *p1;

   for(i=0;i<n;i++)
     {
       p1 =  x+2*i;

       set_invisibles(mem,tt,xb,p1,nm,n,nb,i);  // update de "mem" vector

       h = neighbor(ti,xp,p1,np);               // get closest invisible target

       j = *(mem+n*i+h-1);

       printf("target = %d (%lf,%lf)\n\n",j,*(tt+2*(j-1)+0),*(tt+2*(j-1)+1));

       if (j>0)  move2target(v0,tt,x,xb,i,j,np,nb);
       //else zombi();
     }


   return;
}


int neighbor(int *ti,double *xp,double *p1,int np)
{
  // This function gets the closest point to the
  // current pedestrian. It returns the corresponding
  // target number (say, "mem" position plus one).

  int    h,i,j;
  double xi,yi,xj,yj,dx,dy,dd,rr;

  h  = -1;
  j  = -1;
  dd = BIG*BIG;

  xi = *(p1+0);
  yi = *(p1+1);

  for(i=0;i<np;i++)
    {
      xj = *(xp + 2*i+0);
      yj = *(xp + 2*i+1);

      dx = xi - xj;
      dy = yi - yj;
      rr = dx*dx + dy*dy;

      if (rr<dd) 
        {
          dd = rr;
          h  = i;
        }
    }

  if(h>0) j = *(ti+h);

  return j;
}

int set_invisibles(int *mem,double *tt,double *xb,double *p1,int *nm,int n,int nb,int i)
{
  int    h,j,k;
  double *p2;

  k = 0;

  j = *(nm+i);

  printf("\n");

  for(h=0;h<j;h++)
    {        
      p2 = tt+2*h;

      // Warning: ONLY if the target becomes visible, it switches
      // to a visible target, but formerly visible targets DO NOT
      // switch to invisible ones (this is the consequence of memory). 
      //
      // 0 (visible), 1 (invisible)

      if(!intersection(xb,p1,p2,nb)) *(mem+n*i+h) = -(h+1); // visible 

      if ((*(mem+n*i+h))>0) k++;                            // only count invisibles

      printf("%d ",*(mem+n*i+h));
    }

  printf("\n\n");

  return k;
}

void move2target(double *v0,double *tt,double *x,double *xb,int i,int j,int np,int nb)
{

  double xi,yi,tx,ty,dx,dy,r,vdx,vdy;
  double p1[2],p2[2];
  
  xi = *(x+2*i+0);
  yi = *(x+2*i+1);

  tx = *(tt+2*(j-1)+0);
  ty = *(tt+2*(j-1)+1);

  //----PARALLAX ADJUSTMENT-----------

  dx = tx - xi;
  dy = ty - yi;
          
  r = sqrt(dx*dx + dy*dy);  

  *(p1+0) = xi;
  *(p1+1) = yi;

  *(p2+0) = tx-RADIUS*dy/r;
  *(p2+1) = ty+RADIUS*dx/r;

  if(intersection(xb,p1,p2,nb)) // 0 (visible), 1 (invisible)
    {
      dx = -dx;
      dy = -dy;    
    }
 
  tx = tx-RADIUS*dy/r;
  ty = ty+RADIUS*dx/r;

  dx = tx - xi;
  dy = ty - yi;

  r = sqrt(dx*dx + dy*dy);   
   
  //----------------------------------

  vdx = VD*dx/r;
  vdy = VD*dy/r;
 
  *(v0+2*i+0) = VD*vdx;
  *(v0+2*i+1) = VD*vdy;
        
  return;
}

int targets(int *mem,int *ti,double *tp,double *tt,double *x,double *xb,int *nm,int np,int n,int nb)
{
  int    h,i,j,nt;
  double xi,yi,*p1,*p2;

  nt = 0;

  for(i=0;i<np;i++)
    {
      xi = *(tp+2*i+0);
      yi = *(tp+2*i+1);

      j = (*(ti+i)) - 1;

      *(tt+2*j+0) = xi;
      *(tt+2*j+1) = yi;

      if (nt<j) nt = j;
    }
  
  nt++;
  
  for(h=0;h<n;h++)
    {
      *(nm+h) = nt;

      p1 =  x+2*h;

      for(i=0;i<nt;i++)
        {        
          p2 = tt+2*i;

          // 0 (visible), 1 (invisible)

          if(intersection(xb,p1,p2,nb)) j = i+1; // invisible 
          else j = -(i+1);                       // visible

          *(mem+h*n+i) = j; 
        }
    }

  return 1;
}


/*************************************************************/
/*      AUXILIARY FUNCTIONS: VERLET                          */
/*************************************************************/


void positions(double *x,double *v,double *f,double *f0,double *box,int n)
{
  int    i;
  double xmax,ymax,dim;
  double px,py,vx,vy,fx,fy,dt,dt2;
 
  dt  = STEP;
  dt2 = STEP*STEP/(2.0*M);
 
  xmax = *(box+1);
  ymax = *(box+3);

  if (xmax>ymax) dim = xmax;
  else dim = ymax;

  for(i=0;i<n;i++)
    {
      px = *(x+2*i+0);
      py = *(x+2*i+1);
      vx = *(v+2*i+0);
      vy = *(v+2*i+1);
      fx = *(f+2*i+0);
      fy = *(f+2*i+1);

      px = px + vx*dt + fx*dt2;
      py = py + vy*dt + fy*dt2;

      if (px < 0.0) px += dim;
      if (px > dim) px -= dim;
      if (py < 0.0) py += dim;
      if (py > dim) py -= dim;

      *(x+2*i+0)  = px;
      *(x+2*i+1)  = py;

      *(f0+2*i+0) = fx;
      *(f0+2*i+1) = fy;
      *(f+2*i+0)  = 0.0;
      *(f+2*i+1)  = 0.0;
    }

  return;
}

void velocities(double *v,double *f,double *f0,int n)
{
  int    i;
  double vx,vy,fx,fy,ffx,ffy,dt;

  dt = STEP/(2.0*M);

  for(i=0;i<n;i++)
    {
      vx  = *(v+2*i+0);
      vy  = *(v+2*i+1);
      fx  = *(f+2*i+0);
      fy  = *(f+2*i+1);
      ffx = *(f0+2*i+0);
      ffy = *(f0+2*i+1);

      vx += (fx+ffx)*dt;
      vy += (fy+ffy)*dt;

      *(v+2*i+0) = vx;
      *(v+2*i+1) = vy;
    }

  return;
}


/*************************************************************/
/*      AUXILIARY FUNCTIONS: FORCES                          */
/*************************************************************/

void force_obs(double *xb,double *x,double *v,double *f,int n,int nb)
{
  int    h,hh,i,j;
  double px,py,vx,vy,xi,yi,xj,yj,xs,ys,ri,rj;
  double nx,ny,tx,ty,vn,vt,dv;
  double dx,dy,dd,invdd,t,r,rmin,rmax,dr;
  double f1,f2,frac,fr,fx,fy,ffx,ffy;
  double p1[2],p2[2];
   
  for(i=0;i<n;i++)
    {
      px = *(x+2*i+0);
      py = *(x+2*i+1);
      vx = *(v+2*i+0);
      vy = *(v+2*i+1);

      for(j=0;j<nb;j++)
        {
          xi = *(xb+4*j+0);
          yi = *(xb+4*j+1);
          xj = *(xb+4*j+2);
          yj = *(xb+4*j+3);

          dx = xj - xi;
          dy = yj - yi;

          dd = dx*dx + dy*dy;
   
          invdd = 1.0/dd;
          
          // find closest point (xs,ys)
  
          t = (dx*(px-xi) + dy*(py-yi))*invdd; 
              
          if ((t>0.0) && (t<1.0))  
            {
              xs = xi + t*dx;
              ys = yi + t*dy;
              
              r = sqrt((px-xs)*(px-xs) + (py-ys)*(py-ys));
            }
          else
            {
             ri = sqrt((px-xi)*(px-xi) + (py-yi)*(py-yi));
             rj = sqrt((px-xi-dx)*(px-xi-dx) + (py-yi-dy)*(py-yi-dy));
            
             if (ri<rj) 
               {
                 xs = xi;
                 ys = yi;
                 r  = ri;
               }
             else 
               {
                 xs = xi+dx;
                 ys = yi+dy;
                 r  = rj;
               } 
            }

          rmin = *(table_obs+0);
          rmax = *(table_obs+4*(OBSTAB-1)+0);

          dr = (rmax-rmin)/(double)(OBSTAB-1);

          if (r<rmax)
            {
              *(p1+0) = px;
              *(p1+1) = py;
              *(p2+0) = xs;
              *(p2+1) = ys;

              hh = intersection(xb,p1,p2,nb);
     
              if (!hh)
                {
                  // recall the table format "r,r2,v,f" 
         
                  nx = (px-xs)/r;
                  ny = (py-ys)/r; 

                  tx = -ny;
                  ty =  nx;

                  if (r<rmin) 
                    {
                      vn =  vx*nx + vy*ny;
                      vt = -vx*ny + vy*nx;

                      if (vn<0.0) vn = -vn;

                      *(v+2*i+0) = vn*nx - vt*ny;
                      *(v+2*i+1) = vn*ny + vt*nx;

                      fr = *(table_obs+3); 
                      fx = nx*fr;
                      fy = ny*fr;
                    }
                  else
                    {
                       h = (int)((r-rmin)/dr);
                       frac = (r-(*(table_obs+4*h+0)))/dr;

                       f1 = *(table_obs+4*h+3); 
                       f2 = *(table_obs+4*(h+1)+3); 

                       fr = f1*(1.0-frac)+f2*frac;

                       fx = nx*fr;
                       fy = ny*fr;
                    }   
  
                  *(f+2*i+0) += fx;
                  *(f+2*i+1) += fy;

                  if (r<RADIUS)
                    {
                      dv = vx*tx + vy*ty;

                      ffx = KAPPA*dv*tx*(RADIUS-r);
                      ffy = KAPPA*dv*ty*(RADIUS-r);

                      *(f+2*i+0) -= ffx;
                      *(f+2*i+1) -= ffy;
                    }                     // close (r<2.0*RADIUS)
                }                         // close (!hh)
            }                             // close (r<rmax)
        }                                 // close (j<nb)
    }                                     // close (i<n)

  return;
}

int intersection(double *xb,double *p1,double *p2,int nb)
{
  int    h,i;
  double x1,y1,x2,y2;
  double x3,y3,x4,y4;
  double det,a,b,ta,tb;

  // credit: http://www.cs.swan.ac.uk/~cssimon/line_intersection.html

  h = 0; // 0 = do not intersect. 1 = intersect.

  x1 = *(p1+0);
  y1 = *(p1+1);
  x2 = *(p2+0);
  y2 = *(p2+1);

  for(i=0;i<nb;i++)
    {  
      x3 = *(xb+4*i+0);
      y3 = *(xb+4*i+1);
      x4 = *(xb+4*i+2);
      y4 = *(xb+4*i+3);

      det = (x4-x3) * (y1-y2) - (x1-x2) * (y4-y3);

      if (det)
        {
          a = (y3-y4) * (x1-x3) + (x4-x3) * (y1-y3);
          b = (y1-y2) * (x1-x3) + (x2-x1) * (y1-y3);

          ta = a/det;
          tb = b/det;
      
          if ((ta>0.0) && (ta<1.0) && (tb>0.0) && (tb<1.0)) h = 1;
        }
    }

  return h;
}

void collision(double *x,double *v,double nx,double ny,int ii,int jj)
{
  // compute the (perfectly) elastic collision of two similar atoms.
  // credit: https://stackoverflow.com/questions/35211114/2d-elastic-ball-collision-physics

  double vxi,vyi,vxj,vyj,dvx,dvy;

  vxi = *(v+2*ii+0);
  vyi = *(v+2*ii+1);
  vxj = *(v+2*jj+0);
  vyj = *(v+2*jj+1);

  dvx = vxi - vxj;
  dvy = vyi - vyj;

  vxi = vxi - (dvx*nx+dvy*ny)*nx;
  vyi = vyi - (dvx*nx+dvy*ny)*ny;

  vxj = vxj + (dvx*nx+dvy*ny)*nx;
  vyj = vyj + (dvx*nx+dvy*ny)*ny;

  *(v+2*ii+0) = vxi;
  *(v+2*ii+1) = vyi;
  *(v+2*jj+0) = vxj;
  *(v+2*jj+1) = vyj;

  return;
}

void force_ped(double *xb,double *x,double *v,double *f,double *shift,int i,int j,int np,int nb)
{
  int    h,hh,ii,jj;
  double nx,ny,tx,ty,fx,fy,ffx,ffy;
  double dx,dy,dvx,dvy,dr,dv,fr,r,rmin,rmax;
  double frac,f1,f2;
  double *pi,*pj;
         
  ii = i - np;
  jj = j - np;
                                
  if ((ii>=0) && (jj>=0))
    {
      dx  = (*(x+2*ii+0)) - (*(x+2*jj+0)) - (*(shift+0));
      dy  = (*(x+2*ii+1)) - (*(x+2*jj+1)) - (*(shift+1));

      dvx = (*(v+2*ii+0)) - (*(v+2*jj+0));
      dvy = (*(v+2*ii+1)) - (*(v+2*jj+1));
  
      r = sqrt(dx*dx+dy*dy);

      if (r>0.0)
        {
          nx = dx/r;
          ny = dy/r; 

          rmin = *(table_ped+0);
          rmax = *(table_ped+4*(PEDTAB-1)+0);

          dr = (rmax-rmin)/(double)(PEDTAB-1);

          if (r<rmax) 
            { 
              // check for obstables in between

              pi = x + 2*ii;
              pj = x + 2*jj;

              hh = intersection(xb,pi,pj,nb);
     
              if (!hh)
                {
                  // recall the table format "r,r2,v,f" 
         
                  if (r<rmin) 
                    {
                      collision(x,v,nx,ny,ii,jj);

                      fr = *(table_ped+3);
                    }
                  else
                    {
                      h = (int)((r-rmin)/dr);
                      frac = (r-(*(table_ped+4*h+0)))/dr;

                      f1 = *(table_ped+4*h+3); 
                      f2 = *(table_ped+4*(h+1)+3); 
                      fr = f1*(1.0-frac)+f2*frac;
                    }     

                  fx = dx*fr/r;
                  fy = dy*fr/r;
            
                  *(f+2*ii+0) += fx;
                  *(f+2*ii+1) += fy;
                  *(f+2*jj+0) -= fx;
                  *(f+2*jj+1) -= fy;

                  if (r<2.0*RADIUS)
                    {
                      tx = -ny;
                      ty =  nx;

                      dv = dvx*tx + dvy*ty;

                      ffx = KAPPA*dv*tx*(2.0*RADIUS-r);
                      ffy = KAPPA*dv*ty*(2.0*RADIUS-r);

                      *(f+2*ii+0) -= ffx;
                      *(f+2*ii+1) -= ffy;
                      *(f+2*jj+0) += ffx;
                      *(f+2*jj+1) += ffy;
                    }                     // close (r<2.0*RADIUS)
                }                         // close (!hh)
            }                             // close (r<rmax)
        }                                 // close (r>0.0)
    }                                     // close (ii>np) (jj>np)
 
  return;
}

void force_vd(double *v,double *v0,double *f,int n)
{
  int    i;
  double vx,vy,v0x,v0y,mtau;

  mtau = M/TAU;

  for(i=0;i<n;i++)
    {
      vx  = *(v+2*i+0);
      vy  = *(v+2*i+1);

      v0x  = *(v0+2*i+0);
      v0y  = *(v0+2*i+1);

      *(f+2*i+0) += mtau*(v0x-vx);
      *(f+2*i+1) += mtau*(v0y-vy);
    }

  return;
}


int computelist(double *xb,double *xp,double *x,double *v,double *f,int *head,int *lscl,double *box,int np,int n,int nc,int nb)
{
  int    i,j,c,cx,cy,cc,ccx,ccy;
  double xmax,ymax,dim;
  double shift[2];

  xmax = *(box+1);
  ymax = *(box+3);

  if (xmax>ymax) dim = xmax;
  else dim = ymax;


  for (cx=0; cx<nc; cx++)                         // scan inner cells 
    for (cy=0; cy<nc; cy++)
      {
        c = cx*nc+cy;                             // calculate a scalar cell index
          
          for (ccx=cx-1; ccx<=cx+1; ccx++)        // scan the neighbor cells (including itself)
            for (ccy=cy-1; ccy<=cy+1; ccy++)
              {

                // periodic boundary condition by shifting coordinates

                *(shift+0) = 0.0;
                *(shift+1) = 0.0;

                if (ccx<0) *(shift+0) = -dim;
                else if (ccx>=nc) *(shift+0) = dim;

                if (ccy<0) *(shift+1) = -dim;
                else if (ccy>=nc) *(shift+1) = dim;
                
                cc = ((ccx+nc)%nc)*nc+((ccy+nc)%nc);  
                  
                i = head[c];                             // scan atom i in cell c
                while (i != EMPTY) 
                    {
                      j = head[cc];                      // scan atom j in cell cc

                      while (j != EMPTY) 
                        {
                          if (i < j)                     // avoid double counting and correct image positions 
                            {       
                              force_ped(xb,x,v,f,shift,i,j,np,nb);  
                            }
                          j = lscl[j];
                        }
                      i = lscl[i];
                    }
                }
      }

  force_obs(xb,x,v,f,n,nb);
 
  return 1;
}

/*************************************************************/
/*      AUXILIARY FUNCTIONS: TABLES AND CELLS                */
/*************************************************************/

int buildlist(double *xp,double *x,int *head,int *lscl,double *box,int np,int n,int nc)
{
  int    i,ii,c,cx,cy,ncells;
  double xmax,ymax,rcell;
  
  ncells = nc*nc;

  xmax = *(box+1);
  ymax = *(box+3);

  if (xmax>ymax) rcell = xmax/(double)nc;
  else rcell = ymax/(double)nc;

  for (c=0; c<ncells; c++) head[c] = EMPTY;    // reset the headers, head 

  for (i=0; i<np; i++)                         // scan atoms to head, lscl 
    { 
      cx = (int)(*(xp+2*i+0)/rcell);
      cy = (int)(*(xp+2*i+1)/rcell);
      c  = cx*nc+cy;

      lscl[i] = head[c];                       // link to previous occupant (or EMPTY) 
      head[c] = i;                             // the last one goes to the header
    }

  for (i=np; i<np+n; i++)                      // do the same for the moving pedestrians 
    { 
      ii = i - np;
      cx = (int)(*(x+2*ii+0)/rcell);
      cy = (int)(*(x+2*ii+1)/rcell);
      c  = cx*nc+cy;

      lscl[i] = head[c];                    
      head[c] = i;                          
    }

  return 1;
}



void build_table_ped(double ri,double rf,int ntable)
{
  int     i;
  double  dr,r,r2,rexp,rcexp,v,f;

  dr = (rf-ri)/(double)(ntable-1);

  for(i=0 ; i<ntable ; i++)
    {
      r  = ri + (double)i*dr;
      r2 = r*r;
     
      rexp  = (2.0*RADIUS-r)/B;
      rcexp = (2.0*RADIUS-rf)/B;

      v = A * B * (exp(rexp) - exp(rcexp));
      f = A * exp(rexp);
 
      // store as "r,r2,v,f" (similar to Lammps format)

      *(table_ped+4*i+0) = r;
      *(table_ped+4*i+1) = r2;
      *(table_ped+4*i+2) = v;
      *(table_ped+4*i+3) = f;
    }

  return;
}


void build_table_obs(double ri,double rf,int ntable)
{
  int     i;
  double  dr,r,r2,rexp,rcexp,v,f;

  dr = (rf-ri)/(double)(ntable-1);

  for(i=0 ; i<ntable ; i++)
    {
      r  = ri + (double)i*dr;
      r2 = r*r;
     
      rexp  = (RADIUS-r)/B;
      rcexp = (RADIUS-rf)/B;

      v = A * B * (exp(rexp) - exp(rcexp));
      f = A * exp(rexp);
 
      // store as "r,r2,v,f" (similar to Lammps format)

      *(table_obs+4*i+0) = r;
      *(table_obs+4*i+1) = r2;
      *(table_obs+4*i+2) = v;
      *(table_obs+4*i+3) = f;
    }

  return;
}




/*************************************************************/
/*      AUXILIARY FUNCTIONS: FILE MANAGEMENT                 */
/*************************************************************/

void read_initial(double *x,double *v,double *f,int n)
{
  //*(x+0) = 3.2;
  //*(x+1) = 3.0;
  //*(v+0) = -0.5338;
  //*(v+1) = -0.84561;

  //*(x+0) = 1.5;
  //*(x+1) = 0.5;
  ////*(v+0) = 0.5338;
  *(v+1) = 0.84561;

  *(x+0) = 3.2;
  *(x+1) = 5.0;
  *(v+0) = -0.5338;
  *(v+1) = -0.84561;

  //*(x+0) = 0.5;
  //*(x+1) = 3.5;
  //*(v+0) = -0.5338;
  //*(v+1) = -0.84561;

  *(f+0) = 0.0;
  *(f+1) = 0.0;

  return;
}


void read_myrule(char *myrule,int *ti,double *xp,double *tp,double *box,int np)
{
  int    h,i,j,flag;
  char   trash1[40],trash2[40],trash3[40],trash4[40];
  double xmin,ymin,xmax,ymax,zmin,zmax;
  double xi,yi,zi,txi,tyi,tzi;
  FILE   *fp;

  fp=fopen(myrule,"r");

  flag = fscanf(fp,"%s %s\n",trash1,trash2);
  flag = fscanf(fp,"%d\n",&i);
  flag = fscanf(fp,"%s %s %s %s\n",trash1,trash2,trash3,trash4);
  flag = fscanf(fp,"%d\n",&j);
  flag = fscanf(fp,"%s %s %s",trash1,trash2,trash3);
  flag = fscanf(fp,"%s %s %s\n",trash1,trash2,trash3);
  flag = fscanf(fp,"%lf %lf\n",&xmin,&xmax);
  flag = fscanf(fp,"%lf %lf\n",&ymin,&ymax);
  flag = fscanf(fp,"%lf %lf\n",&zmin,&zmax);
  flag = fscanf(fp,"%s %s %s %s",trash1,trash2,trash3,trash4);
  flag = fscanf(fp,"%s %s %s %s",trash1,trash2,trash3,trash4);
  flag = fscanf(fp,"%s %s\n",trash1,trash2);

  *(box+0) = xmin;
  *(box+1) = xmax;
  *(box+2) = ymin;
  *(box+3) = ymax;
  *(box+4) = zmin;
  *(box+5) = zmax;

  for(h=0;h<np;h++)
    {
      flag = fscanf(fp,"%d %d ",&i,&j);
      flag = fscanf(fp,"%lf %lf %lf ",&xi,&yi,&zi);
      flag = fscanf(fp,"%lf %lf %lf\n",&txi,&tyi,&tzi);

      *(ti+h+0) = j;
      *(xp+2*h+0) = xi;
      *(xp+2*h+1) = yi;
      *(tp+2*h+0) = txi;
      *(tp+2*h+1) = tyi;
    }

  fclose(fp);

  flag++; // useless to avoid warnings

  return;
}


void read_myfile(char *myfile,double *xb)
{
  int    h,i,j,k,flag;
  char   line[100];
  double xi,yi,xj,yj;
  FILE   *fp;

  i = 0;
  j = 0;
  k = 0;

  fp=fopen(myfile,"r");

  while ((h=fscanf(fp,"%s",line))>0)
    {
      if ((*line)=='#') while(getc(fp)!='\n');
      else 
        { 
          i = atoi(line);
          flag = fscanf(fp,"%lf %lf %lf %lf %d",&xi,&yi,&xj,&yj,&j);
          if ((i>0) && (flag==5)) 
            {
              *(xb+4*k+0) = xi;
              *(xb+4*k+1) = yi;
              *(xb+4*k+2) = xj;
              *(xb+4*k+3) = yj;
              k++;
            }
        }
    }
  
  fclose(fp);


  return;
}


int size_of_myfile(char *myfile)
{
  int    h,i,j,k,flag;
  char   line[100];
  double xi,yi,xj,yj;
  FILE   *fp;

  i = 0;
  j = 0;
  k = 0;

  fp=fopen(myfile,"r");

  while ((h=fscanf(fp,"%s",line))>0)
    {
      if ((*line)=='#') while(getc(fp)!='\n');
      else 
        { 
          i = atoi(line);
          flag = fscanf(fp,"%lf %lf %lf %lf %d",&xi,&yi,&xj,&yj,&j);
          if ((i>0) && (flag==5)) k++;
        }
    }
  
  fclose(fp);

  return k;
}

int size_of_myrule(char *myrule)
{
  int    n,step,flag;
  char   trash1[40],trash2[40];
  FILE   *fp;

  fp=fopen(myrule,"r");

  flag = fscanf(fp,"%s %s\n",trash1,trash2);
  flag = fscanf(fp,"%d\n",&step);
  flag = fscanf(fp,"%s %s ",trash1,trash2);
  flag = fscanf(fp,"%s %s\n",trash1,trash2);
  flag = fscanf(fp,"%d\n",&n);

  fclose(fp);

  flag++; // useless to avoid warnings

  return n;
}

int size_of_cells(double *box,double rcell)
{
  int    nc;
  double xmin,xmax,ymin,ymax,dim;

  xmin = *(box+0);
  ymin = *(box+2);

  if ((xmin!=0.0) || (ymin!=0.0)) 
    printf("\nWarning: box out of bounds!\n\n");

  xmax = *(box+1);
  ymax = *(box+3);

  if (xmax>ymax) dim = xmax;
  else dim = ymax;

  if (!(nc=(int)floor(dim/RCELL))) nc++;

  if (nc*nc<9) printf("\nWarning: too few cells. I suggest to resize the original cells\n\n");

  return nc;
}

void checkinput(char **argv,int argc)
{
  char myfile[100],myrule[100];
  FILE *fp;

   if (argc!=3)  
     {
       printf("\nPlease, make sure that you typed:\n\n");
       printf("\t./move.e [myfile.dat] [myfile_rule.lammpstrj].\n\n");
       exit(0);
    }
  else 
    {
      strcpy(myfile,argv[1]);
      strcpy(myrule,argv[2]);

      fp=fopen(myfile,"r");

      if (fp==NULL) 
         {
           printf("\nWrong name for first file!\n\n");
           exit(0);
         }
      else fclose(fp);

      fp=fopen(myrule,"r");

      if (fp==NULL) 
         {
           printf("\nWrong name for second file!\n\n");
           exit(0);
         }
      else fclose(fp);
    }

  return;
}

/*************************************************************/
/*      AUXILIARY FUNCTIONS: PLOTS                           */
/*************************************************************/

double plot_people(char *myplot,double *v0,double *x,double scale,int n)
{
  int    h,i,j,ii,jj,r;
  double xi,yi,vx,vy;
  FILE   *fp;

  // UNITS converts inches to pt

  r  = (int)(scale*RADIUS);

  fp=fopen(myplot,"a");

  fprintf(fp,"newpath\n"); 

  for(h=0;h<n;h++)
    {
      vx = *(v0+2*h+0);
      vy = *(v0+2*h+1);

      xi = (*(x+2*h+0));
      yi = (*(x+2*h+1));
      

      ii = (int)(scale*RADIUS*vx/VD);
      jj = (int)(scale*RADIUS*vy/VD);

      i = (int)(scale*xi) + UNITS*OFFSET;
      j = (int)(scale*yi) + UNITS*OFFSET;

      fprintf(fp,"%d %d moveto\n",i,j);
      fprintf(fp,"%d %d lineto\n",i+ii,j+jj); 
      fprintf(fp,"%d %d moveto\n",i+r,j);
      fprintf(fp,"%d %d %d %d %d arc\n",i,j,r,0,360); 
    }

  fprintf(fp,"1 0.1 0.1 setrgbcolor\n"); 
  fprintf(fp,"stroke\n"); 

  fclose(fp);

  return scale;
}


void close_plot(char *myplot)
{
  FILE *fp;

  fp=fopen(myplot,"a");

  fprintf(fp,"showpage");
  fclose(fp);

  return;
}

double open_plot(char *myplot,char *myfile,double *xb,double *box,int i,int nb)
{
  int    h,j,k;
  char   frame[15];
  double xmin,xmax,ymin,ymax,sx,sy,scale;
  FILE   *fp;

  xmin = *(box+0);
  xmax = *(box+1);
  ymin = *(box+2);
  ymax = *(box+3);

  scale = 1.0;

  // UNITS converts inches to pt

  sx = UNITS*(WIDTH-2.0*OFFSET)/(xmax-xmin);  
  sy = UNITS*(HEIGHT-2.0*OFFSET)/(ymax-ymin);
       
  if (sx<sy) scale = sx;
  else       scale = sy;

  j = 0;

  while((myfile[j]!='.') && (myfile[j]!='\0')) 
     {
       myplot[j]  = myfile[j];
       j++;
     }

  myplot[j]  = '\0';

  itoa(i,frame);

  strcat(myplot,"_");
  strcat(myplot,frame);

  strcat(myplot,".eps");

  fp=fopen(myplot,"w");

  h = (int)(scale*(xmax-xmin)+UNITS*2.0*OFFSET);
  k = (int)(scale*(ymax-ymin)+UNITS*2.0*OFFSET);

  fprintf(fp,"%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(fp,"%%%%BoundingBox: 0 0 %d %d\n",h,k);  // A4: 468 290
  fprintf(fp,"%%%%EndComments\n");

  scale =  boundaries(fp,xb,scale,nb);

  fclose(fp);

  return scale;
}

double boundaries(FILE *fp,double *xb,double scale,int nb)
{
  int    h,i,j,ii,jj;
  double xi,yi,xj,yj;


  fprintf(fp,"newpath\n"); 

  for(h=0;h<nb;h++)
    {
      xi = (*(xb+4*h+0));
      yi = (*(xb+4*h+1));

      xj = (*(xb+4*h+2));
      yj = (*(xb+4*h+3));

      i = (int)(scale*xi) + UNITS*OFFSET;
      j = (int)(scale*yi) + UNITS*OFFSET;

      ii = (int)(scale*xj) + UNITS*OFFSET;
      jj = (int)(scale*yj) + UNITS*OFFSET;

      fprintf(fp,"%d %d moveto\n",i,j); 
      fprintf(fp,"%d %d lineto\n",ii,jj); 
    }

   fprintf(fp,"stroke\n"); 
  
  return scale;
}


/*************************************************************/
/*      AUXILIARY FUNCTIONS: MISCELLANEOUS                   */
/*************************************************************/


void itoa(int n,char s[])
 {
   char c; 
   int  i,j;
     
   i=0;
   do{s[i++]=n%10+'0';} while ((n/=10)>0);   
   s[i]='\0';

   for(i=0,j=strlen(s)-1;i<j;i++,j--) 
     {
       c = s[i];
       s[i]=s[j];
       s[j]=c;
     }
 }


