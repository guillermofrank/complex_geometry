#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#define INFTY  999999999   // "infinite" cost for dijkstra algorithm
#define W      1000.0      // "scale" between nodes

/*
 * rule_triangle.c (compute the target position by the rule of the shortest distance)
 *
 * Notice: the output file resembles a .lammpstrj file. The "type" column corresponds
 *         to the labels of the target set. The null label means the NULL target
 *         (say, NO definite target) 
 *         Make sure you ran "triangle -pDYY myfila_map.poly
 * 
 * gcc -O3 -Wall -O3 -o rule_triangle.e rule_triangle.c -lm
 * 
 * Autor: Guillermo Frank   24 febrero 2020
 *
 * ./rule_triangle.e myfile_map.lammpstrj
 */

int spyfile(char *infile,char *outfile,int *steps,int *atoms);
int read_header_ele(char *infile,char *elefile,char *nodefile,int *nodes);
int read_data_ele_node(char *elefile,char *nodefile,int *cost,int *nodes);
int read_header(FILE *fp,double *dim,int *atoms);
int read_data(FILE *fp,int *id,int *ti,double *x,double *y,double *z,int *c,int *atoms);
int write_header(FILE *fp,double *dim,int *atoms);
int write_data(FILE *fp,double *xx,double *yy,double *tx,double *ty,double *dim,int *ii,int *atoms);
int compute(double *x,double *y,double *xx,double *yy,double *tx,double *ty,int *ii,int *ti,int *c,int *cost,int *atoms,int *nodes,int source);
int dijkstra(int *cost,int source,int target,int n);
int label(int *ii,int *atoms,int *steps);

int main(int argc,char *argv[])
{
  char    infile[100],outfile[100],elefile[100],nodefile[100];
  int     i,*ii,*id,*ti,*c,*cost,steps[1],atoms[1],nodes[1],source;
  double  *x,*y,*z,*xx,*yy,*tx,*ty,dim[6];
  FILE    *fp; 
  clock_t c0,c1;

  c0=clock();

  if (argc==2)
    {
      sscanf(argv[1],"%s",infile);

      spyfile(infile,outfile,steps,atoms);

      read_header_ele(infile,elefile,nodefile,nodes);

      // Get data. 

      if ((*steps) && (*nodes))
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

          cost = (int*)malloc(((*nodes)+1)*((*nodes)+1) * sizeof(int));

          read_data_ele_node(elefile,nodefile,cost,nodes);

          fp=fopen(infile,"r");

          printf("\nDecision making to closest invisible:\n\n");

          for(i=0;i<(*steps);i++)
            {
              read_header(fp,dim,atoms);
              read_data(fp,id,ti,x,y,z,c,atoms);

              source = (*atoms)-1+i+1;
              compute(x,y,xx,yy,tx,ty,ii,ti,c,cost,atoms,nodes,source);
            }
            
          fclose(fp);

          printf("\n");

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

          free(cost);
        }
      else printf("\nSorry, no matching atoms.\n\n");

    }
  else printf("\nPlease, indicate the name of the input file (i.e.: myfile_map.lammpstrj).\n\n");

  c1=clock();

  printf("\nCPU  time (sec.): %f\n",(float)(c1-c0)/CLOCKS_PER_SEC);

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

int read_header_ele(char *infile,char *elefile,char *nodefile,int *nodes)
{
  int    h,i,j,k,n,num,attr,flag,maximum;
  FILE   *fp;

  strcpy(elefile,infile);
  strcpy(nodefile,infile);

  i = 0;

  n = 0;

  * nodes = 0;

  maximum = 0;

  while (*(elefile+i)!='_') i++;

  *(elefile+i)='\0';

  strcat(elefile,"_map.1.ele");

  i = 0;

  while (*(nodefile+i)!='_') i++;

  *(nodefile+i)='\0';

  strcat(nodefile,"_map.1.node");


  fp=fopen(elefile,"r");

  if (fp!=NULL)
    { 
      flag = fscanf(fp,"%d %d %d\n",&n,&num,&attr);

      if (flag==3) 
        {
          h = 0;
       
          flag = 4;

          while ((h<=n) && (flag==4))
            {
              flag = fscanf(fp,"%d %d %d %d\n",&h,&i,&j,&k);
          
              if (i>maximum) maximum = i;
              if (j>maximum) maximum = j;
              if (k>maximum) maximum = k;
              
            }

          *nodes = maximum;
        }
      else printf("\nWarning: something is wrong with the header of %s\n\n",elefile);
    }      
  else 
    {
      printf("\nWrong element's file!\n\n");
      exit(0);
    }
 
  fclose(fp);

  return 1;
}


int read_data_ele_node(char *elefile,char *nodefile,int *cost,int *nodes)
{
  int    h,i,j,k,m,n,flag;
  double *x,xi,yi,xj,yj,xk,yk,dx,dy,d,dmax,attr;
  FILE   *fp;

  h = 0;

  n = (*nodes);

  x = (double*)malloc(2 * n * sizeof(double));
 
  // We first read the nodes positions

  fp=fopen(nodefile,"r");

  flag = fscanf(fp,"%d %d %d %d\n",&m,&h,&i,&j);
  
  flag = 5;

  h = 0;

  while ((h<=m) && (flag==5))
    {
      flag = fscanf(fp,"%d %lf %lf %lf %d\n",&h,&xi,&yi,&attr,&k);
    
      *(x+2*(h-1)+0) = xi;
      *(x+2*(h-1)+1) = yi;
    }

  fclose(fp);

  // We secondly open the elements file  and compute 
  // the distances between nodes

  h = 0;

  dmax = 0.0;

  for(i=0;i<n*n;i++) *(cost+i) = INFTY;

  fp=fopen(elefile,"r");

  flag = fscanf(fp,"%d %d %d\n",&m,&i,&j);
  
  flag = 4;

  while ((h<=m) && (flag==4))
    {
      flag = fscanf(fp,"%d %d %d %d\n",&h,&i,&j,&k);

      xi = *(x+2*(i-1)+0);
      yi = *(x+2*(i-1)+1);

      xj = *(x+2*(j-1)+0);
      yj = *(x+2*(j-1)+1);

      xk = *(x+2*(k-1)+0);
      yk = *(x+2*(k-1)+1);

      dx = xi - xj;
      dy = yi - yj;
     
      d  = sqrt(dx*dx+dy*dy);

      if (d>dmax) dmax = d;

      *(cost+n*(i-1)+(j-1)) = (int)(W*d);
      *(cost+n*(j-1)+(i-1)) = (int)(W*d);

      dx = xi - xk;
      dy = yi - yk;
     
      d  = sqrt(dx*dx+dy*dy);

      if (d>dmax) dmax = d;

      *(cost+n*(i-1)+(k-1)) = (int)(W*d);
      *(cost+n*(k-1)+(i-1)) = (int)(W*d);

      dx = xj - xk;
      dy = yj - yk;
     
      d  = sqrt(dx*dx+dy*dy);

      if (d>dmax) dmax = d;

      *(cost+n*(k-1)+(j-1)) = (int)(W*d);
      *(cost+n*(j-1)+(k-1)) = (int)(W*d); 
    }

  fclose(fp);

  free(x);

  if (((int)(W*dmax))>=INFTY) printf("\nWarning: distance between nodes exceeds INFTY value!\n\n"); 

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



int compute(double *x,double *y,double *xx,double *yy,double *tx,double *ty,int *ii,int *ti,int *c,int *cost,int *atoms,int *nodes,int source)
{
  // This function computes the 'target' position
  // to the shortest invisible distance.

  int    h,i,imin,vi,m,n,d,dmin;
  double xi,yi,xj,yj;

  // "n" is the number of points (atoms) 
  // belonging to the borders.
  // Targets point to 0...n-1
  // The NULL target is n.
  //
  // "m" is the size of the cost matrix

  m = (*nodes);
  n = (*atoms); 

  imin = n;

  h  = (*(c+n-1));
  xj = (*(x+n-1)); // observer position
  yj = (*(y+n-1)); // observer position
  
  *(ii+h) = imin;  // default target label (NULL)
  *(tx+h) = 0.0;   // default target position
  *(ty+h) = 0.0;   // default target position

  dmin = INFTY;

  for(i=0;i<n-1;i++)
    {
      vi = (*(ti+i))-1; // invisible = 0; visible = 1;

      if (!vi)
        {
          xi = *(x+i); // target position
          yi = *(y+i); // target position

          // compute the Dijkstra algorithm
          // source = (atoms - 1) + (step number) + 1;
          // target = i+1;
         
          d = dijkstra(cost,source,i+1,m); 

          if ((d<dmin) && (d>0))
            { 
              dmin = d;
              imin = i;  
            }
        }     
    }

  *(xx+h) = xj;        // observer position
  *(yy+h) = yj;        // observer position

  if (imin<n)
    {
      xi = *(x+imin);
      yi = *(y+imin);

      *(ii+h) = imin;  // target label (starting from 0) 
      *(tx+h) = xi;    // target position
      *(ty+h) = yi;    // target position

      printf("observer: %d (%lf,%lf) --> ",source-n+1,xj,yj);
      printf("target: %d (%lf,%lf)\t",imin+1,xi,yi);
      printf("distance: %d\n",d);
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

int dijkstra(int *cost,int source,int target,int n)
{
  // credit: https://www.codewithc.com/dijkstras-algorithm-in-c/
  //
  // This is the usual Dijkstra algoritm. 
  // The "cost" matrix is an n*n matrix of weights between nodes. 
  // "source" and "target" are numbered 1..n (although "cost" is 
  // is filled with data indices from 0 to n-1).
  // The function returns the minimum cumulative weights (distance)
  // between the "source" node and the "target" node.
  // If source=target, the function returns 0.

  int *dist,*prev,*selected,*path,c,i,m,min,s,t,d,j;

  dist = (int*)malloc(n * sizeof(int));
  prev = (int*)malloc(n * sizeof(int));
  selected = (int*)malloc(n * sizeof(int));
  path = (int*)malloc(n * sizeof(int));

  for(i=0;i<n;i++)
    {
      dist[i] = INFTY;
      prev[i] = -1;
      selected[i] = 0;
    }

    s = source-1;   // "s" is the starting node (the -1 shifts index to 0)
    t = target-1;   // "t" is the ending node (the -1 shifts index to 0)

    selected[s]=1;
    dist[s]=0;

     while(selected[t]==0)
        {
          m=n+1;;
          min=INFTY;
        
         for(i=0;i<n;i++)
           {
             d=dist[s]+cost[n*s+i]; 

             if (d<dist[i] && selected[i]==0)
               {
                 dist[i]=d;
                 prev[i]=s;
               }
             if (dist[i]<min && selected[i]==0)
               {
                 min=dist[i];
                 m=i;
               }
           }

         if (m<n+1) {s=m; selected[s]=1;} 
         else {s=-1; break;}
        }

  j = 0;

  while(s != -1)
    {
      *(path+j) = s+1;
      s = prev[s];
      j++;
    }
  
  c = dist[t];

  free(dist);
  free(prev);
  free(selected);
  free(path);

  return c;
}


