#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define MAXLINE   10000       // maximum number of segments (borders)
#define MAXPOINTS 100000      // maximum number of points 
#define BIG       10000000.0  // any big number
#define NUM       10          // default number of points (pedestrians)
#define NUMB      100         // default granularity of segments (borders)
#define STEPS     10          // numbers of grains in borders to analyse
#define DSTEP     20          // default granularity for deeper search (borders)
#define WIDTH     6.5         // .eps width 
#define HEIGHT    4.017       // .eps height 0.618 x 6.5 = 4.017
#define UNITS     72          // .eps units (in 'pt' = 1/72 inch)
#define OFFSET    0.2         // .eps offset for the bounding box
#define RADIUS    0.2         // .eps radius of points
#define PLOT      0           // 0 = no plot, 1 = plot


/* Este programa abre el archivo geometry.dat 
 * separa contornos visibles de invisibles
 *
 * Uso: gcc -Wall -o hope.e hope.c -lm
 * 
 * ./hope.e myfile.dat [number of points]
 *
 */

double boundaries(FILE *fp,double *x,double *xt,double *xm,int k);
int    points(FILE *fp,double *xc,double *xm,double scale,int n);
int    view(double *xb,int *vi,double x1,double y1,int k);
int    intersection(int *vi,double *xb,double *x,double *xt,double *xc,int h,int nc,int k,int nk);
int    write_header(FILE *fp,double xmin,double ymin,double xmax,double ymax,int h,int k,int step);
int    write_data(FILE *fp,double *xb,int *vi,double px,double py,int h,int k,int step);
int    borders(FILE *fp,double *xb,double *xm,int *vi,double px,double py,double scale,int n);
void   itoa(int n,char s[]);

int main(int argc,char *argv[])
{
  int    h,i,j,k,n,num,flag,*p,*vi;
  char   geometry[100],number[50],outfile[100];
  char   epsfile[100],viewfile[100],viewer[15];
  double xi,yi,xj,yj,xmin,xmax,ymin,ymax,scale;
  double px,py,dx,dy,nx,ny,dd,rx,ry,la,mu,xs,ys,ds,dmin,mumin;
  double *x,*xt,*xn,*xc,*xm,*pt,*xb;
  FILE   *fp,*fpp;

  h = 1;
  i = 0;
  j = 0;
  k = 0;
  n = 0;

  num  = NUM;
  flag = 0;

  xmin =  BIG;
  ymin =  BIG;
  xmax = -BIG;
  ymax = -BIG;

  if((argc==2) || (argc==3)) 
   {
     if (argc==3) sscanf(argv[2],"%d",&num);
     else printf("\nNo number of points requested. Switching to default n=%d\n\n",num);

     sscanf(argv[1],"%s",geometry);

     if((fp=fopen(geometry,"r"))!=NULL)
      {
        p   = (int *)malloc(MAXLINE*sizeof(int));
        xm  = (double *)malloc(4*sizeof(double));
        x   = (double *)malloc(2*MAXLINE*sizeof(double));
        xt  = (double *)malloc(2*MAXLINE*sizeof(double));
        xn  = (double *)malloc(2*MAXLINE*sizeof(double));
        xc  = (double *)malloc(2*MAXPOINTS*sizeof(double));
        pt  = (double *)malloc(2*MAXPOINTS*sizeof(double));

        while((h=fscanf(fp,"%s",number))>0)
          {
            if ((*number)=='#') while(getc(fp)!='\n');
            else 
              { 
                j = 0;
                i = atoi(number);
                flag = fscanf(fp,"%lf %lf %lf %lf %d",&xi,&yi,&xj,&yj,&j);

                if ((i>0) && (abs(j)==1)) 
                  {
                    if (xmin>xi) xmin = xi;
                    if (xmin>xj) xmin = xj;
                    if (xmax<xi) xmax = xi;
                    if (xmax<xj) xmax = xj;
                    if (ymin>yi) ymin = yi;
                    if (ymin>yj) ymin = yj;
                    if (ymax<yi) ymax = yi;
                    if (ymax<yj) ymax = yj;

                    *(p+k) = i;

                    *(x+2*k+0) = xi;
                    *(x+2*k+1) = yi;

                    *(xt+2*k+0) = xj-xi;
                    *(xt+2*k+1) = yj-yi;

                    *(xn+2*k+0) = (double)j*(yj-yi);
                    *(xn+2*k+1) = (double)j*(xi-xj);

                    k++;
                  }
                else 
                  {
                    printf("\nCorrupted data. Check your data file...\n\n");
                    break;
                  }
               
                if (k==MAXLINE) 
                  {
                    printf("\nData exceeds the maximun allowed segments. Please modify MAXLINE.\n\n");
                    break;
                  }
              }
          }

        // locate points in allowed regions

        h = 0;
        i = num;

        while((h<i) && (n<10*MAXPOINTS))
          {
            dmin = BIG;

            mumin = 0.0;

            px = (double)rand()/(double)RAND_MAX;
            py = (double)rand()/(double)RAND_MAX;

            px = xmin + px*(xmax-xmin);
            py = ymin + py*(ymax-ymin);

            // compute all distances to borders

            for(j=0;j<k;j++)
              {
                xi = *(x+2*j+0);
                yi = *(x+2*j+1);
                dx = *(xt+2*j+0);
                dy = *(xt+2*j+1);
                nx = *(xn+2*j+0);
                ny = *(xn+2*j+1);

                dd = dx*dx + dy*dy;

                rx = dx/dd;
                ry = dy/dd;
           
                la = rx*(px-xi) + ry*(py-yi);
                
                if (la < 0.0) la = 0.0;
                if (la > 1.0) la = 1.0;

                xs = xi + la*dx;
                ys = yi + la*dy;

                // accept closest border from inside

                ds = (xs-px)*(xs-px) + (ys-py)*(ys-py); // squared distance 

                // I keep comments on the following during the testing stage:
                // mu = 0.0;
                // if ((la>0.0) && (la<1.0)) mu = nx*(xs-px) + ny*(ys-py);

                mu = nx*(xs-px) + ny*(ys-py);  

                if(ds<dmin) 
                 {
                   dmin = ds;
                   mumin = mu;

                   *(xc+2*h+0) = px; // observers
                   *(xc+2*h+1) = py;
                   *(pt+2*h+0) = xs; // intersections
                   *(pt+2*h+1) = ys;
                 }

              } // end of for   j<n

            if ((dmin<BIG) && (mumin>0.0)) h++;

            n++;
          }     // end of while h<i

        fclose(fp);
    
        i = 0;

        while((geometry[i]!='.') && (geometry[i]!='\0')) 
          {
            epsfile[i]  = geometry[i];
            outfile[i]  = geometry[i];
            i++;
          }

        epsfile[i]  = '\0';
        outfile[i]  = '\0';

        strcat(epsfile,".eps");
        strcat(outfile,"_map.lammpstrj");

        *(xm+0) = xmin;
        *(xm+1) = ymin;
        *(xm+2) = xmax;
        *(xm+3) = ymax;

        fp=fopen(epsfile,"w");

        fprintf(fp,"%%!PS-Adobe-3.0 EPSF-3.0\n");
        fprintf(fp,"%%%%BoundingBox: 0 0 468 290\n");
        fprintf(fp,"%%%%EndComments\n");
        
        scale = boundaries(fp,x,xt,xm,k);
        points(fp,xc,xm,scale,h);

        fprintf(fp,"showpage");
        fclose(fp);

        // compute visibility of borders from each observer

        vi = (int *)malloc(NUMB*k*sizeof(int));
        xb = (double *)malloc(2*NUMB*k*sizeof(double));

        for(i=0;i<k;i++)
          {
            xi = *(x+2*i+0);
            yi = *(x+2*i+1);
            dx = *(xt+2*i+0);
            dy = *(xt+2*i+1);

            for(j=0;j<NUMB;j++)
              {
                 *(xb+2*NUMB*i+2*j+0) = xi + ((double)j/(double)NUMB)*dx;
                 *(xb+2*NUMB*i+2*j+1) = yi + ((double)j/(double)NUMB)*dy;
              }
          }
     
        if (n<10*MAXPOINTS)
         {
           fpp=fopen(outfile,"w"); 

           for(i=0;i<h;i++) 
             {
               px = *(xc+2*i+0);
               py = *(xc+2*i+1);

               //view(xb,vi,px,py,k);
               intersection(vi,xb,x,xt,xc,i,h,k,NUMB*k);

               // PLOT into .eps (optional)

               if (PLOT)
                 {
                   itoa(i,viewer);

                   j = 0;

                   while((geometry[j]!='.') && (geometry[j]!='\0')) 
                     {
                       viewfile[j] = geometry[j];
                       j++;
                     }

                   viewfile[j]  = '\0';

                   strcat(viewfile,"_");
                   strcat(viewfile,viewer);
                   strcat(viewfile,".eps");

                   fp=fopen(viewfile,"w");

                   fprintf(fp,"%%!PS-Adobe-3.0 EPSF-3.0\n");
                   fprintf(fp,"%%%%BoundingBox: 0 0 468 290\n");
                   fprintf(fp,"%%%%EndComments\n");
        
                   scale = boundaries(fp,x,xt,xm,k);
                   borders(fp,xb,xm,vi,px,py,scale,NUMB*k);

                   fprintf(fp,"showpage");
                   fclose(fp);
                 }

               // end of PLOT

               write_header(fpp,xmin,ymin,xmax,ymax,NUMB,k,i);
               write_data(fpp,xb,vi,px,py,NUMB,k,i);

              }

           fclose(fpp);
         }
        else printf("\nSorry, no points could be established.\n\n");


        free(p);
        free(x);
        free(xm);
        free(xt);
        free(xn);
        free(xc);
        free(pt);
        free(vi);
        free(xb);

      }
     else printf("\nWrong file name.\n\n");
   }
  else   printf("\nNo input file.\n\n");


  flag++; // useless to avoid warnings 

  return 1;
}
 

double boundaries(FILE *fp,double *x,double *xt,double *xm,int k)
{
  int    h,i,j,ii,jj,ic,jc;
  double xmin,xmax,ymin,ymax;
  double sx,sy,scale;
  double xi,yi,xj,yj;

  xmin = *(xm+0);
  ymin = *(xm+1);
  xmax = *(xm+2);
  ymax = *(xm+3);

  scale = 1.0;

  // UNITS converts inches to pt

  sx = UNITS*(WIDTH-2.0*OFFSET)/(xmax-xmin);  
  sy = UNITS*(HEIGHT-2.0*OFFSET)/(ymax-ymin);
       
  if (sx<sy) scale = sx;
  else       scale = sy;

  ic = (UNITS*(WIDTH-2.0*OFFSET)-(int)(scale*(xmax-xmin)))/2;
  jc = (UNITS*(HEIGHT-2.0*OFFSET)-(int)(scale*(ymax-ymin)))/2;

  fprintf(fp,"newpath\n"); 

  for(h=0;h<k;h++)
    {
      xi = (*(x+2*h+0));
      yi = (*(x+2*h+1));

      xj = (*(xt+2*h+0)) + xi;
      yj = (*(xt+2*h+1)) + yi;

      i = (int)(scale*xi) + UNITS*OFFSET + ic;
      j = (int)(scale*yi) + UNITS*OFFSET + jc;

      ii = (int)(scale*xj) + UNITS*OFFSET + ic;
      jj = (int)(scale*yj) + UNITS*OFFSET + jc;

      fprintf(fp,"%d %d moveto\n",i,j); 
      fprintf(fp,"%d %d lineto\n",ii,jj); 
    }

   fprintf(fp,"stroke\n"); 
  
  return scale;
}

int points(FILE *fp,double *xc,double *xm,double scale,int n)
{
  int    h,i,j,r,ic,jc;
  double xmin,xmax,ymin,ymax;
  double xi,yi;

  // UNITS converts inches to pt

  r  = (int)(scale*RADIUS);

  xmin = *(xm+0);
  ymin = *(xm+1);
  xmax = *(xm+2);
  ymax = *(xm+3);

  ic = (UNITS*(WIDTH-2.0*OFFSET)-(int)(scale*(xmax-xmin)))/2;
  jc = (UNITS*(HEIGHT-2.0*OFFSET)-(int)(scale*(ymax-ymin)))/2;

  fprintf(fp,"newpath\n"); 

  for(h=0;h<n;h++)
    {
      xi = (*(xc+2*h+0));
      yi = (*(xc+2*h+1));

      i = (int)(scale*xi) + UNITS*OFFSET + ic;
      j = (int)(scale*yi) + UNITS*OFFSET + jc;

      fprintf(fp,"%d %d moveto\n",i+r,j); 
      fprintf(fp,"%d %d %d %d %d arc\n",i,j,r,0,360); 
    }

   fprintf(fp,"1 0.1 0.1 setrgbcolor\n"); 
   fprintf(fp,"stroke\n"); 
  
  return 1;
}

int view(double *xb,int *vi,double x1,double y1,int k) 
{
  // This function classifies the border points as 'visible'
  // or 'invisible'. It returns the 'vi' vector with '1' if 
  // visible, or, 0 is invisible.

  int    h,i,j,hh,ii,jj,steps,visible;
  double xi,yi,xj,yj,dx,dy;
  double x2,y2,x3,y3,x,y,a,b,c;

  steps = STEPS;

  for(i=0;i<NUMB*k;i++) *(vi+i) = 1; 

  for(i=0;i<k;i++) 
    for(j=0;j<NUMB-steps;j++)
      {
        x2 = *(xb+2*NUMB*i+2*j+0);
        y2 = *(xb+2*NUMB*i+2*j+1);
        x3 = *(xb+2*NUMB*i+2*(j+steps)+0);
        y3 = *(xb+2*NUMB*i+2*(j+steps)+1);

        // check if (x,y) lies inside the triangle 
        // (x1,y1) , (x2,y2) , (x3,y3)

        for(ii=0;ii<k;ii++) 
          {
            xi = *(xb+2*NUMB*ii+0);
            yi = *(xb+2*NUMB*ii+1);
            xj = *(xb+2*NUMB*ii+2);
            yj = *(xb+2*NUMB*ii+3);

            dx = (xj - xi)/(double)DSTEP;
            dy = (yj - yi)/(double)DSTEP;

            for(jj=0;jj<NUMB-1;jj++)
               {
                 visible = 1;

                 hh = 0;

                 x = *(xb+2*NUMB*ii+2*jj+0);
                 y = *(xb+2*NUMB*ii+2*jj+1);

                 while (visible && (hh<DSTEP)) 
                   {
                     x += dx;
                     y += dy;
                  
                     a = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
                     b = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
                     c = 1.0 - a - b;

                     if((a>0.0) && (a<=1.0) && (b>0.0) && (b<=1.0) && (c>0.0) && (c<=1.0))
                       {
                         visible = 0;
                         for(h=0;h<steps+1;h++) *(vi+NUMB*i+j+h) = 0;
                       }

                     hh++;
                   }
               }
          } 
      }
     
  return 1;
}

int borders(FILE *fp,double *xb,double *xm,int *vi,double px,double py,double scale,int n)
{
  int    h,i,j,r,s,ic,jc;
  double xmin,xmax,ymin,ymax;
  double xi,yi;

  // UNITS converts inches to pt

  r  = (int)(scale*0.5*RADIUS);

  xmin = *(xm+0);
  ymin = *(xm+1);
  xmax = *(xm+2);
  ymax = *(xm+3);

  ic = (UNITS*(WIDTH-2.0*OFFSET)-(int)(scale*(xmax-xmin)))/2;
  jc = (UNITS*(HEIGHT-2.0*OFFSET)-(int)(scale*(ymax-ymin)))/2;

  fprintf(fp,"newpath\n"); 

  i = (int)(scale*px) + UNITS*OFFSET + ic;
  j = (int)(scale*py) + UNITS*OFFSET + jc;

  fprintf(fp,"%d %d moveto\n",i+2*r,j); 
  fprintf(fp,"%d %d %d %d %d arc\n",i,j,2*r,0,360); 

  fprintf(fp,"0.1 0.1 1 setrgbcolor\n"); 
  fprintf(fp,"stroke\n"); 


  fprintf(fp,"newpath\n"); 

  for(h=0;h<n;h++)
    {
      s  = *(vi+h);

      if (s==0)
        {
          xi = (*(xb+2*h+0));
          yi = (*(xb+2*h+1));

          i = (int)(scale*xi) + UNITS*OFFSET + ic;
          j = (int)(scale*yi) + UNITS*OFFSET + jc;

          fprintf(fp,"%d %d moveto\n",i+r,j); 
          fprintf(fp,"%d %d %d %d %d arc\n",i,j,r,0,360); 
        }
    }

   fprintf(fp,"1 0.1 0.1 setrgbcolor\n"); 
   fprintf(fp,"stroke\n"); 

   fprintf(fp,"newpath\n"); 

   for(h=0;h<n;h++)
    {
      s  = *(vi+h);

      if (s==1)
        {
          xi = (*(xb+2*h+0));
          yi = (*(xb+2*h+1));

          i = (int)(scale*xi) + UNITS*OFFSET + ic;
          j = (int)(scale*yi) + UNITS*OFFSET + jc;

          fprintf(fp,"%d %d moveto\n",i+r,j); 
          fprintf(fp,"%d %d %d %d %d arc\n",i,j,r,0,360); 
        }
    }

   fprintf(fp,"0.1 1 0.1 setrgbcolor\n"); 
   fprintf(fp,"stroke\n"); 
 
  return 1;
}

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

int write_header(FILE *fp,double xmin,double ymin,double xmax,double ymax,int h,int k,int step)
{
  // Write the header of the ouput file in a .lammpstrj format
  //
  // Caution:
  //
  // TIMESTEP corresponds to the id of the inner particle 
  // h = granularity (NUMB)
  // k = borders 

  fprintf(fp,"ITEM: TIMESTEP\n%d\n",step);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n%d\n",h*k+1);
  fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
  fprintf(fp,"%lf %lf\n%lf %lf\n%lf %lf\n",xmin,xmax,ymin,ymax,0.0,2.0*RADIUS);
  fprintf(fp,"ITEM: ATOMS id type x y z c_group\n");
  
  return 1;
}

int write_data(FILE *fp,double *xb,int *vi,double px,double py,int h,int k,int step)
{
  // Write data in .lammpstrj format, as follows:
  // 
  // the first NUMB*k lines corrrespond to the border atoms 
  // (with non-visible = type 1 and visible = type 2)
  // 
  // The last line corresponds to the inner atom (type 3) 
  //
  // The c_group means the border number (for type 1 and 2) 
  // or the inner number (for type 3)
  //

  int    i,j,s;
  double x,y;
  
  for(i=0;i<k;i++)
    for(j=0;j<h;j++)
      {
        s = *(vi+h*i+j);
        x = *(xb+2*h*i+2*j+0);
        y = *(xb+2*h*i+2*j+1);

        fprintf(fp,"%d %d %lf %lf %lf %d\n",h*i+j+1,s+1,x,y,RADIUS,i+1);
      }

  fprintf(fp,"%d %d %lf %lf %lf %d\n",h*k+1,3,px,py,RADIUS,step);

  return 1;
}

int intersection(int *vi,double *xb,double *x,double *xt,double *xc,int h,int nc,int k,int nk)
{
  int    i,j;
  double x1,y1,x2,y2;
  double x3,y3,x4,y4,dx,dy;
  double det,a,b,ta,tb;

  // credit: http://www.cs.swan.ac.uk/~cssimon/line_intersection.html

  // vi = vector of visibilities for xb points (size nk = NUMB*k)
  // xb = vector of border points (size nk = 2*NUMB*k)
  // x  = vector of initial points for borders (size 2*k)
  // xt = vector of tangential direction of borders (size 2*k)
  // xc = vector of pedestrian points (size 2*nc)
  // h  = position of pedestrian point

  x1 = *(xc+2*h+0);
  y1 = *(xc+2*h+1);

  for(i=0;i<nk;i++)
    {
      *(vi+i) = 1; // 1 = do not intersect. 0 = intersect.

      x2 = *(xb+2*i+0);
      y2 = *(xb+2*i+1);

      for(j=0;j<nc;j++)
        { 
          x3 = *(x+2*j+0);
          y3 = *(x+2*j+1);

          dx = *(xt+2*j+0);
          dy = *(xt+2*j+1);

          x4 = x3 + dx;
          y4 = y3 + dy;

          det = (x4-x3) * (y1-y2) - (x1-x2) * (y4-y3);

          if (det)
            {
              a = (y3-y4) * (x1-x3) + (x4-x3) * (y1-y3);
              b = (y1-y2) * (x1-x3) + (x2-x1) * (y1-y3);

              ta = a/det;
              tb = b/det;
      
              if ((ta>0.0) && (ta<1.0) && (tb>0.0) && (tb<1.0)) *(vi+i) = 0;
            }
        }
    }
  return h;
}


