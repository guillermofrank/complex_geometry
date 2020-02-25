#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

/*
 * get_map.c (get the visibility map of a single atom)
 *
 * gcc -O3 -Wall -O3 -o get_map.e get_map.c -lm
 * 
 * Autor: Guillermo Frank   13 enero 2020
 *
 * get_map myfile_map.lammpstrj [number of particle]
 */

int  read_header(FILE *fp,double *dim,int *atoms);
int  read_data(FILE *fp,int *id,int *ti,double *x,double *y,double *z,int *c,int *atoms);
int  write_header(FILE *fp,double *dim,int *atoms,int step);
int  write_data(FILE *fp,int *id,int *ti,double *x,double *y,double *z,int *c,int *atoms);
void itoa(int n,char s[]);


int main(int argc,char *argv[])
{
  char    infile[100],outfile[100],number[15];
  int     h,i,step,atoms[1],*id,*ti,*c;
  double  dim[6],*x,*y,*z;
  FILE    *fp,*fpp; 

  h = 0;
  i = 0;

  if (argc==3)
    {
      sscanf(argv[1],"%s",infile);
      fp=fopen(infile,"r");

      sscanf(argv[2],"%d",&h); 

      if (fp!=NULL)
      {        
        while (*(infile+i)!='.') i++;
        *(infile+i)='\0';

        itoa(h,number);
        strcpy(outfile,infile);
        strcat(outfile,"_");
        strcat(outfile,number);
        strcat(outfile,".lammpstrj");       

        while ((step=read_header(fp,dim,atoms))>=0)
          {
            id = (int*)malloc((*atoms) * sizeof(int));
            ti = (int*)malloc((*atoms) * sizeof(int));
            x  = (double*)malloc((*atoms) * sizeof(double));
            y  = (double*)malloc((*atoms) * sizeof(double));
            z  = (double*)malloc((*atoms) * sizeof(double));
            c  = (int*)malloc((*atoms) * sizeof(int));

            read_data(fp,id,ti,x,y,z,c,atoms);

            if (step==h)
              {
                fpp=fopen(outfile,"w");

                write_header(fpp,dim,atoms,step);
                write_data(fpp,id,ti,x,y,z,c,atoms);

                fclose(fpp);

                break;
              }

            free(id);
            free(ti);
            free(x);
            free(y);
            free(z);
            free(c);
          }
    
         if (step<0) printf("\nSorry, no matching atom.\n\n");
     }
     else printf("\nWrong input file!\n\n");

      
    }
  else printf("\nPlease, indicate the name of the input file and the selected atom.\n\n");


  return 1;
}


int read_header(FILE *fp,double *dim,int *atoms)
{
  int    n,step,flag;
  char   trash1[40],trash2[40],trash3[40],trash4[40];
  double xmin,ymin,xmax,ymax,zmin,zmax;

  flag = fscanf(fp,"%s %s\n",trash1,trash2);
  flag = fscanf(fp,"%d\n",&step);

  if (flag!=1) step = -1;
  else 
    {
      flag = fscanf(fp,"%s %s %s %s\n",trash1,trash2,trash3,trash4);
      flag = fscanf(fp,"%d\n",&n);
      flag = fscanf(fp,"%s %s %s",trash1,trash2,trash3);
      flag = fscanf(fp,"%s %s %s\n",trash1,trash2,trash3);
      flag = fscanf(fp,"%lf %lf\n",&xmin,&xmax);
      flag = fscanf(fp,"%lf %lf\n",&ymin,&ymax);
      flag = fscanf(fp,"%lf %lf\n",&zmin,&zmax);
      flag = fscanf(fp,"%s %s %s %s",trash1,trash2,trash3,trash4);
      flag = fscanf(fp,"%s %s %s %s\n",trash1,trash2,trash3,trash4);
    
      *atoms = n;
  
      *(dim+0) =  xmin;
      *(dim+1) =  ymin;
      *(dim+2) =  xmax;
      *(dim+3) =  ymax;
    }

  return step;
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

int write_header(FILE *fp,double *dim,int *atoms,int step)
{
  // Write the header of the ouput file in a .lammpstrj format
  //
  // Caution:
  //
  // TIMESTEP corresponds to the id of the inner particle 
  // h = granularity (NUMB)
  // k = borders 

  int    n;
  double xmin,xmax,ymin,ymax,zmin,zmax;

  n = *atoms;

  xmin = (*(dim+0)); 
  xmax = (*(dim+1));
  ymin = (*(dim+2));  
  ymax = (*(dim+3));
  zmin = (*(dim+4));
  zmax = (*(dim+6));

  fprintf(fp,"ITEM: TIMESTEP\n%d\n",step);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n%d\n",n);
  fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
  fprintf(fp,"%lf %lf\n%lf %lf\n%lf %lf\n",xmin,xmax,ymin,ymax,zmin,zmax);
  fprintf(fp,"ITEM: ATOMS id type x y z c_group\n");
  
  return 1;
}

int write_data(FILE *fp,int *id,int *ti,double *x,double *y,double *z,int *c,int *atoms)
{
  int    h,i,j,k,flag;
  double xi,yi,zi;

  for(h=0;h<*atoms;h++)
    {
      i  = (*(id+h));
      j  = (*(ti+h));
      xi = (*(x+h));
      yi = (*(y+h));
      zi = (*(z+h));
      k  = (*(c+h));

      flag = fprintf(fp,"%d %d %lf %lf %lf %d\n",i,j,xi,yi,zi,k);
    }

  flag++; // useless to avoid warnings

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


