#include<stdio.h>    
#include<stdlib.h>  
#include<math.h>
#include<time.h>
#include <stdbool.h>
#define PI 3.141592654

int look2(int*vector,int**matrix,int vn,int mn,int j);
int look(int*vector,int**matrix,int row,int col);
int max_array(int *a, int num_elements);
int ** Init_Matrix_2D(int  nx, int  ny);
double gamratio(double x1,double x2,double y1,double y2);
void award(int*A,int*B,int**C,int**D,int*A2,int*B2,int**C2,int**D2,int n,int m);
int check(int **lpropc,int **gpropc,int *ldegree ,int *gdegree,int n);
void dohits(int **lpropc,int **gpropc,int *ldegree,int *gdegree,int start,int n,int *hit);
long int factorial(int a);
int sum(int *pa, int num_elements);
long int choose(int n,int k);
double randn();
double randomInteger01(void);
int rand_int(int a,int b);
int randomIntegern(int n);
void loglikelihood(double *result,double lamdaL,double lamdaG,int **gpropc,int **lpropc,int*ldegree,int*gdegree,int *locsize,int n,int pop);
void display_mat_deg(int ** c, int * d, int n);

int main()
{
    int n=71,pop=217,m=7; /*m is the maximum group size, n is the final size*/
    int i,ii,j,k,kk,*grsize,*orio,*dg,*dl,*lpropd,*gpropd,*linf;
    int Nllinks,lsum;
    int miko,row,col,gia,top,ll,mikro2=0;
    double propl,propg,numer[2],denom[2],mean,mikro1=0,u,up,acc;
    double lambdaL=0.001,lambdaG=0.001,cholprob=0.5,slamL=0.01,slamG=0.01,prratio,lmean=0,gmean=0;  /*lthresh, gthresh*/
    long double infper;
    int burn=20000,samgap=700,samsize=10000,sumdat,looplen,counter,p,run=0,A=0,B=0,C=0,D=0,A2=0,B2=0,C2=0,D2=0,A3=0,A4=0,B3=0,B4=0,C3=0,C4=0,D3=0,D4=0;
    bool done;
    int **dat, **gpropc,**lpropc,**cg,**cl,**who;
    FILE *fd, *fc;

    dat=Init_Matrix_2D(m, m);
    gpropc=Init_Matrix_2D(n, n);
    lpropc=Init_Matrix_2D( n, m);
    cg=Init_Matrix_2D(n, n);
    cl=Init_Matrix_2D(n, m);
    who=Init_Matrix_2D(n, m); 
    
    orio=(int*)malloc((m+1)*sizeof(int));
    grsize=(int*)malloc(n*sizeof(int)); 
    linf=(int*)malloc(n*sizeof(int));
    dg=(int*)calloc(n,sizeof(int));
    dl=(int*)calloc(n,sizeof(int));
    lpropd=(int*)malloc(n*sizeof(int));
    gpropd=(int*)malloc(n*sizeof(int));

    if((fd=fopen("cowling.dat","w"))==NULL)   
        printf("File could not opened\n"); 
 if((fc=fopen("mikro.dat","w"))==NULL)   
        printf("File could not opened\n"); 
 
    fprintf(fd,"sum sum max max lambdaL lambdaG  counter|\n");
    fprintf(fd," dl  dg  dl  dg                         |\n");
    fprintf(fd,"----------------------------------------|\n");
    infper = 4.1;

        clock_t start = clock();

dat[0][2]=20;     
dat[0][3]=17;
dat[0][4]=7;
dat[0][5]=1;

dat[1][2]=2;
dat[1][3]=5;
dat[1][4]=3;

dat[2][2]=1;
dat[2][3]=1;
	

   sumdat = 0;  /*Count the number of infected households*/
    for(i=0; i<m; i++)       
    {
        for(j=i; j<m; j++)
        {
	      sumdat=sumdat+dat[i][j]; 
		
	}
    }  
                /*Create the group size of each individual*/          
    k = 0;       
    for(i=0; i<m; i++)   
    {
         for(j=i; j<m; j++)
         { 
              kk=(i+1)* dat[i][j]; 
	      if(kk>0)
 	      {
	           for(ii=0; ii<kk; ii++)  
	           { 
			  
	                grsize[k] = j+1; 
                        printf("%d ",grsize[k]);
   		        k=k+1; 			
                   }
              }     
	 }
    }  
               
   orio[0]=0;
   for(i=1; i<m+1; i++)       
   {
        for(j=i-1; j<m; j++)  
        { 
             orio[i] = orio[i] + (i)*dat[i-1][j];
        }
        orio[i] = orio[i] + orio[i-1];            
   } 
                 
   printf("\norio:\n");					
   for (i=0; i<m+1; i++)				 
   {	
       printf("%d ",orio[i]); 
   }


 miko=1;
   row=0;
   for(i=0; i<m; i++)  
   {
        for(j=i; j<m; j++)
        {   
             kk=(i+1)* dat[i][j]; 
             if(kk>0)  
	     {       
                for(ll=0; ll<dat[i][j]; ll++)
                {

                  top=row+kk/dat[i][j];
	          for(gia=row; gia<top; gia++)
		  {         
                       for(col=0; col<kk/dat[i][j]; col++)
                       {
		            who[gia][col]=miko;
		            miko++;	
		       }
		       miko=miko-kk/dat[i][j];
		  }
		  row=row+kk/dat[i][j];	
		  miko=miko+kk/dat[i][j];	
                  }
	     }
	}
   }
	
   printf("\nthe who matrix\n");
   for(i=0;i<n; i++)
   {  
      for(j=0;j<m; j++)
      {
          printf("%d ",who[i][j]);
      }
      printf("\n");
   }


 printf("\nlinf:");
   for(i=0; i<n; i++)
   {
        lsum=0;
        for(j=0; j<m; j++)
	{
	     if(who[i][j]>0)
	     {
	          lsum=lsum+1;
	     }
	}
	linf[i]=lsum;  /*The final size in i's household*/
        printf("%d ",linf[i]); 
   }

   lsum=0;
   for(i=1; i<m; i++)  
   {
        for(j=i; j<m; j++)
	{
	     lsum=lsum+(i)*(i+1)*dat[i][j];
	}
   }
   Nllinks = lsum;
printf(" Nllinks: %d \n",Nllinks);

   for (i=0; i<n-1; i++)
   {
   	cg[i][0]=i+2;   //At the biggining global links only
        dg[i]=1;
   }



   //look(dl,cl,90,m);
   //look(dg,cg,50,30);
   srand((unsigned int)time(NULL));
   srand(57);

   burn = burn-(2*burn);
   looplen=samgap*samsize;
   miko=0;
   for(counter=burn; counter<=looplen; counter++)
   {    

        propl=fabs(randn(lambdaL,slamL));
	mikro1=mikro1+propl;
	mikro2=mikro2+1;
        propg=fabs(randn(lambdaG,slamG));
        loglikelihood(numer,propl,propg,cl,cg,dl,dg,grsize,n,pop);
        if(numer[1]==1)   
        {   
	     loglikelihood(denom,lambdaL,lambdaG,cl,cg,dl,dg,grsize,n,pop);
           
	     prratio= gamratio(propl,propg,lambdaL,lambdaG);
	     acc= prratio+numer[0]-denom[0];
             u=randomInteger01();
             if(log(u)<acc)
	     { 	 
                 miko++;	
	         lambdaL=propl; 
                 lambdaG=propg;	
             }
        }    
        /*update graph */
        p=randomIntegern(2); /*if 1->add link if 2->delete link  */
        if(p==1)             /*we pick to add either local or global*/
        {    
             up=randomInteger01(); 
             if((up<cholprob)&&(Nllinks-(sum(dl,n))>0)) /*Choose to add local*/  
	     {    
                  A2++;
    	          i=randomIntegern(Nllinks-(sum(dl,n)));              //  1  //
	          j=orio[1]+1;  
	          k=0;

	          do
                  {    
	               k=k+(linf[j-1]-1)-dl[j-1];               
   	               j++;
                  }while(k+(linf[j-1]-1)-dl[j-1]<i);          
                        
                  if(dl[j-1]<linf[j-1]-1)
                 {   
                       A3++;
                       do
                       {           /*Now j is the lobaluy from where the (local)link will emanate*/
	                    done=true;
	                    k=rand_int((who[j-1][0]),(who[j-1][0]+linf[j-1]-1));
  
  	                    if (k==j)
  	 	            {	
  		                 done=false;
  		            }
		            else if(dl[j-1]>=1)
		            {	
		                 for(i=0; i<dl[j-1]; i++)
	 		         {  
	 		             if(cl[j-1][i]==k)
	 		             {	
	 		 	         done=false;
	 	                     }
	 	 	         }
	 	            }     
                       }while(done==false);
                 
                       award(lpropd,gpropd,gpropc,lpropc,dl,dg,cg,cl,n,m);
                           
 	               lpropd[j-1]++;  /*Register the extra link*/     
 	               lpropc[j-1][dl[j-1]] = k; /*Draw the extra link*/

                       if(check(lpropc,gpropc,lpropd,gpropd,n)==1)
	               {   
                            A4++;  
       	                    acc=exp(lambdaL*infper)-1;
		            /*acc = mgf((-1.0d0)*lambdaL)-1*/
	                    acc = acc*((double)(Nllinks-sum(dl,n)))/((double)(sum(dl,n)+1)); 
	                    u=randomInteger01();
	                    if((u)<acc)
                            {    A++;
	                         award(dl,dg,cg,cl,lpropd,gpropd,gpropc,lpropc,n,m);
		            }
		       }      
                  }      	     
             }	
             else if((sum(dg,n))<n*(n-1))  /*Choose to add a global edge*/
	     {	 
		  /*pick a random vertex to add an edge to*/
                  B2++;
		  i=randomIntegern(n*(n-1)-(sum(dg,n)));               //  2  //
		  j=1;
		  k=0;
	          do
		  {    
		       k=k+(n-1)-dg[j-1];
    		       j++;
		  }while (k+(n-1)-dg[j-1]<i);    /*pick a random vertex to link to*/
                  
                  do
                  {	
		       done=true; 
     		       k=randomIntegern(n);
		       if((k==j)||(dg[j-1]==n-1))
		       {   
			    done=false;
		       }
		       else if(dg[j-1]>=1)
		       {  
			    for(i=0; i<dg[j-1]; i++)
			    {
				 if(cg[j-1][i]==k)   
				 {	
				      done=false;
				 }
			    }    
		       }		
                  }while(done==false);
                   
                                     /*Set up the proposed new graph*/  
	          award(lpropd,gpropd,gpropc,lpropc,dl,dg,cg,cl,n,m);  

		  gpropd[j-1]++;   /*Register the extra link*/
                  gpropc[j-1][dg[j-1]] = k;  /*Draw the extra link*/
                  B3++;
		  if (check(lpropc,gpropc,lpropd,gpropd,n)==1)
                  {    B4++;
                      
		       acc=exp((lambdaG*infper)/((double)pop)) -1;
		       acc=acc*((double)((n*(n-1))-(sum(dg,n))))/((double)(sum(dg,n)+1));
		       u=randomInteger01();
		       if((u)<acc)
                       {    B++;
			    award(dl,dg,cg,cl,lpropd,gpropd,gpropc,lpropc,n,m);
                       }
                  }     
             }
        }    
        else if(p==2)   /*Delete an edge*/
	{    
	     up=randomInteger01();
	     if((up<cholprob) && sum(dl,n)>0)  /*Choose to delete local*/
	     {    C2++;
                  i=randomIntegern(sum(dl,n));     /*pick a random edge*/    //  3  //
		  j=orio[1]+1; 
		  k=0;     
		  do
		  {   
		       k=k+dl[j-1];
    	               j++;
		  }while (k+dl[j-1]<i);
		  k=i -k; 
		  award(lpropd,gpropd,gpropc,lpropc,dl,dg,cg,cl,n,m); 
	          if (k<dl[j-1] && k>0)
		  {    
		       for(i=(k-1); i<(dl[j-1]-1); i++) 
 		       {
			    lpropc[j-1][i]= lpropc[j-1][i+1];
		       }
		  }
                  C3++;
 		  lpropc[j-1][dl[j-1]-1]=0;
                  lpropd[j-1]--;
 		  loglikelihood(numer,lambdaL,lambdaG,lpropc,gpropc,lpropd,gpropd,grsize,n,pop);
		  if (numer[1]==1)
 	 	  {    C4++;
		       loglikelihood(denom,lambdaL,lambdaG,cl,cg,dl,dg,grsize,n,pop);
		       acc= numer[0] -denom[0];
		       u=randomInteger01();
  		       if (log(u)<acc)            /*Update Graph*/
  		       {    C++;
                            award(dl,dg,cg,cl,lpropd,gpropd,gpropc,lpropc,n,m);    
                       }
                  }
             }		    
             else if(sum(dg,n)>(sumdat-1))   /*Choose global*/
             {    D2++;
	 	  i=randomIntegern(sum(dg,n));  /*pick a random edge*/      //  4  //
		  j=1;
  		  k=0;
  		  do
		  {    
		       k=k+dg[j-1];
  		       j++;
		  }while(k+dg[j-1]<i);
		  k=i-k;
		  award(lpropd,gpropd,gpropc,lpropc,dl,dg,cg,cl,n,m); 
                  //look2(gpropd,gpropc,n,n,j);

		  if (k<dg[j-1] && k>0)
		  {    
		       for(i=(k-1); i<(dg[j-1]-1); i++)
		       {
	                    gpropc[j-1][i]= gpropc[j-1][i+1];
 		       }
		  }
                  D3++;
   		  gpropc[j-1][dg[j-1]-1]=0;
 		  gpropd[j-1]--;
	          loglikelihood(numer,lambdaL,lambdaG,lpropc,gpropc,lpropd,gpropd,grsize,n,pop);
 	          if(numer[1]==1)
		  {    D4++;
	               loglikelihood(denom,lambdaL,lambdaG,cl,cg,dl,dg,grsize,n,pop);
 		       acc=numer[0] -denom[0];
		       u=randomInteger01();
 		       if ((log(u))<acc) 
		       {    D++;
		            award(dl,dg,cg,cl,lpropd,gpropd,gpropc,lpropc,n,m);
		       }
  		   }
             }
        }      
        if ((counter>0) && (counter%samgap==0)) 
        {     
              run++;
              lmean=lmean+lambdaL;
              gmean=gmean+lambdaG;
              fprintf(fd," %d  %d  %d   %d  %f %f  %d    \n",sum(dl,n),sum(dg,n),max_array(dl,n),max_array(dg,n),lambdaL,lambdaG,counter); 
       

        }            
   } 

   // look(dl,cl,50,m);
   //look(dg,cg,80,30);
   printf("\nlocal mean=%f,global mean=%f run=%d\n",(lmean/(run)),(gmean/(run)),run);
   printf("add l: A2=%d->A=%d~%f\nadd g: B2=%d->B=%d~%f\ndel l: C2=%d->C=%d~%f\ndel g: D2=%d->D=%d~%f\n",A2,A,((double)(((double)(A))/((double)(A2)))),B2,B,((double)(((double)(B))/((double)(B2)))),C2,C,((double)(((double)(C))/((double)(C2)))),D2,D,((double)(((double)(D))/((double)(D2)))));
   printf("sunolo eginan %d protinomenes allages,dektike tis %d\n",A2+B2+C2+D2,A+B+C+D);
    printf("A2=%d A3=%d A4=%d A=%d\n",A2,A3,A4,A);
    printf("B2=%d B3=%d B4=%d B=%d\n",B2,B3,B4,B);
    printf("C2=%d C3=%d C4=%d C=%d\n",C2,C3,C4,C);
    printf("D2=%d D3=%d D4=%d D=%d\n",D2,D3,D4,D);
   printf(" LOCAL: add/delete= %d / %d = %f\n Global: add/delete= %d / %d = %f\n",A,C,((double)(((double)(A))/((double)(C)))),B,D,((double)(((double)(B))/((double)(D)))));

   printf("from the %d times,we accept the %d ,α=%f\n",(-burn+looplen),miko,((double)(((double)(miko))/((double)(-burn+looplen)))));
	mean=mikro1/mikro2;
	printf("%f\n",mean);
   free(linf);
   free(orio);
   free(dg);             
   free(dl);            
   free(gpropd);  
   free(lpropd);
   free(grsize);

clock_t end = clock();
float seconds = (float)(end - start) / CLOCKS_PER_SEC;
printf("sec=%f\n",seconds);

   return (0);
}
              
/* the end of program start of function */

/* Generate a random integer number sto(1,n) */
int randomIntegern(int n)
{
    int k;
    double number;
   
    number = (double)rand()/((double) RAND_MAX+1);
    k=(int)(number*n);
    return (k+1);
}

int rand_int(int a,int b)
{
    if (a > b)
        return((rand() % (a-b+1)) + b);
    else if (b > a)
        return((rand() % (b-a+1)) + a);
    else
        return a;
}

/* epistrefei apo mia omoiomorfh se [0,1] */

double randomInteger01(void)
{
    double number;
    
    number = (double)rand()/((double) RAND_MAX+1);
    return (number);
}

double randn (double mu, double sigma)
{
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;
    
    if (call == 1)
    {
        call = !call;
        return (mu + sigma* (double) X2);
    }
    
    do
    {
        U1 = -1 + ((double) rand () / RAND_MAX) * 2;
        U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
    }
    while (W >= 1 || W == 0);
    
    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;
    
    call = !call;
    
    return (mu + sigma * (double) X1);
}


void loglikelihood(double *result,double lamdaL,double lamdaG,int **lpropc,int **gpropc,int*ldegree,int*gdegree,int *locsize,int n,int pop)
{
     double llike1,sum6,sum7,sum8,sum9,I=4.1,total;
     int llike2,i;
     llike1=0;    
     llike2=0;
     total=0;
     if(check(lpropc,gpropc,ldegree,gdegree,n)==1)
     {  
        llike2=1;
	for(i=0; i<n; i++) /*gia kathe individual*/
	{
		sum6=(log(1-exp((-1)*lamdaL*(I))))*((double)(ldegree[i]));      
		sum7=(log(1-exp((-1)*lamdaG*(I)/pop)))*((double)(gdegree[i]));
		sum8=((-1)*lamdaL*(I))*((double)(locsize[i]-ldegree[i]-1));
  		sum9=((-1)*lamdaG*(I)/pop)*((double)(pop-gdegree[i]-1));
		
		total=total+sum6+sum7+sum8+sum9; 
	}
        llike1=total;
     }
     result[0]=llike1;  /*return the value of the likelihood if calculated*/
     result[1]=llike2;  /*return if graph is valid */
}


int check(int **lpropc,int **gpropc,int *ldegree ,int *gdegree,int n)
{
    int checkout,*hit,start;
    
    hit=(int*)calloc(n,sizeof(int));
    start=1;
    hit[start-1]=1;
    dohits(lpropc,gpropc,ldegree,gdegree,start,n,hit);
    if(sum(hit,n)==n)
    {
        checkout=1;
    }
    else
    {
        checkout=0;
    }

    return checkout;
}

void dohits(int **lpropc,int **gpropc,int *ldegree,int *gdegree,int start,int n,int *hit)
{
    int i;
    if(gdegree[start-1]>0)
    {
        for(i=0; i<gdegree[start-1]; i++)  
        {
            if(hit[gpropc[start-1][i]-1]==0)  /*Look for global hits */
            {
                hit[gpropc[start-1][i]-1]=1;
                dohits(lpropc,gpropc,ldegree,gdegree,gpropc[start-1][i],n,hit); 
            }
        }
    }
    if(ldegree[start-1]>0)
    {
        for(i=0; i<ldegree[start-1]; i++)  
        {                               /* Look for local hits */
            if(hit[lpropc[start-1][i]-1]==0)  
            {
                hit[lpropc[start-1][i]-1]=1;
                dohits(lpropc,gpropc,ldegree,gdegree,lpropc[start-1][i],n,hit); 
            }
        }
    }

}

long int choose(int n,int k)
{
    double result;
    int i;
    
    if(k>n) return(0);
    
    if(k>n/2) k=n-k;
    
    result=1;
    for(i=1;i<=k;i++)
        result*=(double)(n+1-i)/(i);
    
    return ((long int) result);
}

int sum(int *pa, int num_elements)
{
    int i, sum1=0;
    for (i=0; i<num_elements; i++)
    {
        sum1 = sum1 + pa[i];
    }
    return(sum1);
}
void award(int*A,int*B,int**C,int**D,int*A2,int*B2,int**C2,int**D2,int n,int m)
{
    int i,j;
    for(i=0; i<n; i++)           
	{
	  A[i]=A2[i];
	  B[i]=B2[i];  
      	      for(j=0; j<n; j++)
              {
            	   C[i][j]=C2[i][j];  /*(nxn)*/
              }
	      for(j=0; j<m; j++)
	      {
	      	   D[i][j]=D2[i][j]; 
	      }
         }
}


double gamratio(double x1,double x2,double y1,double y2)
{
  double ratio,lamb=0.0001,ni=0.0001;
  ratio=(lamb*(y1+y2-x1-x2))*((ni-1)*(log(x1*x2)-log(y1*y2)));
  return(ratio);
}

  
int ** Init_Matrix_2D(int nx,int ny) 
{ 
 int **X; 
 int  i; 
 X = (int **)calloc(nx,sizeof(int *)); 
 if (X == NULL) 
 	{ 
	     return NULL; 
	} 
        int  s = nx*ny;  
	int* p = (int *)calloc(s, sizeof(int)); 
	if (p == NULL)  
	{ 
	     free(X); 
	     return NULL; 
	} 	 
	for (i = 0; i < nx; i++) 
	{ 
	     X[i] = p; 
	     p += ny; 
	} 		 
 return X; 
}

int max_array(int *a, int num_elements)
{
   int i, max=0;
   for (i=0; i<num_elements; i++)
   {
	 if (a[i]>max)
	 {
	    max=a[i];
	 }
   }
   return(max);
}

int look(int*vector,int**matrix,int row,int col)
{
     int i,j;
     printf("\nthe matrix:\n");       
     for(i=0; i<row;i++)
          {
               for(j=0; j<col;j++)
               {  
                    printf(" %d ",matrix[i][j]);
               }
               printf("\n");   
	  } 
          printf("\nand the vector -->:\n");
          for(j=0; j<row;j++)
          {  
               printf(" %d ",vector[j]);
          }
return 0;
}

int look2(int*vector,int**matrix,int vn,int mn,int j)
{
     int i;
     printf("\nthe j's matrix row:\n");       
     for(i=0; i<mn; i++)
          {
              printf(" %d ",matrix[j-1][i]);
          }
          printf("\n");   
          printf("and the vector -->:\n");
          for(i=0; i<vn;i++)
          {  
               printf(" %d ",vector[i]);
          }
return 0;
}

