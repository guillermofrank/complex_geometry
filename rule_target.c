#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define BIG    10000000.0  // any big number

/*
 * rule_target.c (compute the target position by the rule of the shortest distance)
 *
 * Notice: the output file resembles a .lammpstrj file. The "type" column corresponds
 *         to the labels of the target set. The null label means the NULL target
 *         (say, NO definite target) 
 * 
 * gcc -O3 -Wall -O3 -o rule_target.e rule_target.c -lm
 * 
 * Autor: Guillermo Frank   22 enero 2020
 *
 * ./rule_target.e myfile_map.lammpstrj
 */

int  spyfile(char *infile,char *outfile,int *steps,int *atoms);
int  read_header(FILE *fp,double *dim,int *atoms);
int  read_data(FILE *fp,int *id,int *ti,double *x,double *y,double *z,int *c,int *atoms);
int  write_header(FILE *fp,double *dim,int *atoms);
int  write_data(FILE *fp,double *xx,double *yy,double *tx,double *ty,double *dim,int *ii,int *atoms);
int  compute(double *x,double *y,double *xx,double *yy,double *tx,double *ty,int *ii,int *ti,int *c,int *atoms);
int  label(int *ii,int *atoms,int *steps);

int main(int argc,char *argv[])
{
  char    infile[100],outfile[100];
  int     i,*ii,*id,*ti,*c,steps[1],atoms[1];
  double  *x,*y,*z,*xx,*yy,*tx,*ty,dim[6];
  FILE    *fp; 

  if (argc==2)
    {
      sscanf(argv[1],"%s",infile);

      spyfile(infile,outfile,steps,atoms);

      // Get data. 

      if ((*steps))
        {
          xx = (double*)malloc((*steps) * sizeof(double));
          yy = (double*)malloc((*steps) * sizeof(double));
          tx = (double*)malloc((*steps) * sizeof(double));
          ty = (double*)malloc((*steps) * sizeof(double));
          ii = (int*)malloc((*steps) * sizeof(int));

          id = (int*)malloc((*atoms) * sizeof(int));
          ti = (int*)malloc((*atoms) * sizeof(int));
          x  = (double*)malloc((*atoms) * sizeof(double));
          y  = (double*)malloc((*atoms) * sizeof(double));
          z  = (double*)malloc((*atoms) * sizeof(double));
          c  = (int*)malloc((*atoms) * sizeof(int));

          fp=fopen(infile,"r");

          for(i=0;i<(*steps);i++)
            {
              read_header(fp,dim,atoms);
              read_data(fp,id,ti,x,y,z,c,atoms);
              compute(x,y,xx,yy,tx,ty,ii,ti,c,atoms);
            }
            
          fclose(fp);

          // re-label the targets to consecutive numbers 1...max

          label(ii,atoms,steps);

          fp=fopen(outfile,"w");

          write_header(fp,dim,steps);
          write_data(fp,xx,yy,tx,ty,dim,ii,steps);

          fclose(fp);


          free(id);
          free(ti);
          free(x);
          free(y);
          free(z);
          free(c);

          free(xx);
          free(yy);
          free(tx);
          free(ty);
          free(ii);
        }
      else printf("\nSorry, no matching atoms.\n\n");

    }
  else printf("\nPlease, indicate the name of the input file (i.e.: myfile_map.lammpstrj).\n\n");


  return 1;
}

int spyfile(char *infile,char *outfile,int *steps,int *atoms)
{
  int    i,flag;
  char   word[40],trash1[40],trash2[40];
  FILE   *fp;

  fp=fopen(infile,"r");

  if (fp!=NULL)
    { 
      while((flag = fscanf(fp,"%s",word))==1)
        {
          if(!strcmp(word,"TIMESTEP")) 
           {
             flag = fscanf(fp,"%d\n",steps);
             flag = fscanf(fp,"%s %s ",trash1,trash2);
             flag = fscanf(fp,"%s %s\n",trash1,trash2);
             flag = fscanf(fp,"%d\n",atoms);
           }
        }

      (*steps)++;

      strcpy(outfile,infile);
   
      i = 0;

      while (*(outfile+i)!='_') i++;

      *(outfile+i)='\0';

      strcat(outfile,"_rule.lammpstrj");
    }      
  else printf("\nWrong input file!\n\n");
 
  fclose(fp);


  return 1;
}

int read_header(FILE *fp,double *dim,int *atoms)
{
  int    n,step,flag;
  char   trash1[40],trash2[40],trash3[40],trash4[40];
  double xmin,ymin,xmax,ymax,zmin,zmax;

  flag = fscanf(fp,"%s %s\n",trash1,trash2);
  flag = fscanf(fp,"%d\n",&step);
  flag = fscanf(fp,"%s %s %s %s\n",trash1,trash2,trash3,trash4);
  flag = fscanf(fp,"%d\n",&n);
  flag = fscanf(fp,"%s %s %s",trash1,trash2,trash3);
  flag = fscanf(fp,"%s %s %s\n",trash1,trash2,trash3);
  flag = fscanf(fp,"%lf %lf\n",&xmin,&xmax);
  flag = fscanf(fp,"%lf %lf\n",&ymin,&ymax);
  flag = fscanf(fp,"%lf %lf\n",&zmin,&zmax);
  flag = fscanf(fp,"%s %s %s %s",trash1,trash2,trash3,trash4);
  flag = fscanf(fp,"%s %s %s %s\n",trash1,trash2,trash3,trash4);

  *(dim+0) = xmin;
  *(dim+1) = xmax;
  *(dim+2) = ymin;
  *(dim+3) = ymax;
  *(dim+4) = zmin;
  *(dim+5) = zmax;

  flag++; // useless to avoid warnings

  return 1;
}

int read_data(FILE *fp,int *id,int *ti,double *x,double *y,double *z,int *c,int *atoms)
{
  int    h,i,j,k,flag;
  double xi,yi,zi;

  for(h=0;h<*atoms;h++)
    {
      flag = fscanf(fp,"%d %d %lf %lf %lf %d\n",&i,&j,&xi,&yi,&zi,&k);

      *(id+h) = i;
      *(ti+h) = j;

      *(x+h)  = xi;
      *(y+h)  = yi;
      *(z+h)  = zi;
      *(c+h)  = k;
    }

  flag++; // useless to avoid warnings.

  return 1;
}

int write_header(FILE *fp,double *dim,int *atoms)
{
  // Write the header of the ouput file in a .lammpstrj format
  //
  // Caution:
  //
  // Only the inner atoms are included.  
  // The velocity corresponds to the 'e_d' unit vector.

  double xmin,xmax,ymin,ymax,zmin,zmax;

  xmin = (*(dim+0)); 
  xmax = (*(dim+1));
  ymin = (*(dim+2));  
  ymax = (*(dim+3));
  zmin = (*(dim+4));
  zmax = (*(dim+5));

  fprintf(fp,"ITEM: TIMESTEP\n%d\n",0);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n%d\n",*atoms);
  fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
  fprintf(fp,"%lf %lf\n%lf %lf\n%lf %lf\n",xmin,xmax,ymin,ymax,zmin,zmax);
  fprintf(fp,"ITEM: ATOMS id type x y z tx ty tz\n");
  
  return 1;
}

int write_data(FILE *fp,double *xx,double *yy,double *tx,double *ty,double *dim,int *ii,int *atoms)
{
  int    i,jj,flag;
  double xi,yi,zi,xj,yj,zj;

  zj = (*(dim+5))/2.0; // observer
  zi = zj;             // target

  for(i=0;i<*atoms;i++)
    {
      xj = (*(xx+i));  // observer
      yj = (*(yy+i));  // observer
      xi = (*(tx+i));   // target
      yi = (*(ty+i));   // target
      jj = (*(ii+i));   // target

      flag = fprintf(fp,"%d %d %lf %lf %lf %lf %lf %lf\n",i+1,jj,xj,yj,zj,xi,yi,zi);
    }

  flag++; // useless to avoid warnings

  return 1;
}



int compute(double *x,double *y,double *xx,double *yy,double *tx,double *ty,int *ii,int *ti,int *c,int *atoms)
{
  // This function computes the 'target' position
  // to the shortest invisible distance.

  int    h,i,imin,vi,n;
  double xi,yi,xj,yj,dx,dy,dd,ddmin;

  // "n" is the number of points (atoms) 
  // belonging to the borders.
  // Targets point to 0...n-1
  // The NULL target is n.

  n = (*atoms); 

  imin = n;

  h  = (*(c+n-1));
  xj = (*(x+n-1)); // observer position
  yj = (*(y+n-1)); // observer position
  
  *(ii+h) = imin;  // default target label (NULL)
  *(tx+h) = 0.0;   // default target position
  *(ty+h) = 0.0;   // default target position

  ddmin = BIG*BIG;

  for(i=0;i<n-2;i++)
    {
      vi = (*(ti+i))-1; // invisible = 0; visible = 1;

      if (!vi)
        {
          xi = *(x+i);
          yi = *(y+i);

          dx = xi - xj;
          dy = yi - yj;

          dd = dx*dx + dy*dy;

          if (dd<ddmin) 
            { 
              ddmin = dd;
              imin  = i;
            }
        }     
    }

  *(xx+h) = xj;        // observer position
  *(yy+h) = yj;        // observer position

  if (imin<n)
    {
      xi = *(x+imin);
      yi = *(y+imin);

      *(ii+h) = imin;  // target label 
      *(tx+h) = xi;    // target position
      *(ty+h) = yi;    // target position
    }

  return imin;
}

int label(int *ii,int *atoms,int *steps)
{
  int h,hh,j,jj,n,*c;

  j = 1;

  n = (*atoms);

  c = (int*)malloc((n+1) * sizeof(int));

  *(c+n) = 0;

  for(h=0;h<n;h++) *(c+h) = -1;

  for(h=0;h<(*steps);h++)
    {
      // hh corresponds to imin=0...n 
      // hh = n means the NULL target

      hh = *(ii+h);  
      jj = *(c+hh);

      if (jj==-1) 
        {
          *(c+hh) = j;
          j++;
        }
    }


  for(h=0;h<(*steps);h++)
    {
      hh = *(ii+h);
      jj = *(c+hh);

      *(ii+h) = jj;
    }

  free(c);

  return 1;
}
