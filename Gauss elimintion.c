//Gauss elimination using partial pivoting
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
int main(){
int i,j,k,n,t; //declaring variables
printf("No. of rows="); // to write no. of rows of augmented matrix
scanf("%d",&n); // n = no. of rows
float a[n][n+1],x[n],b[n],sum,constant,temp,max;
FILE *fp;
fp=fopen("input_3.txt","r"); //to read augmented input matrix from a file
while(!feof(fp))
{
  for(i=0;i<n;i++)
  {
    for(j=0;j<n+1;j++)
     {
       fscanf(fp,"%f",&a[i][j]);
     }
  }
}
fclose(fp);
for (j=0;j<n;j++)
 {
   max=a[j][j];
   for(i=j+1;i<n;i++)
    {
       if(a[i][j]>max) // partial pivoting
        {
            max=a[i][j]; //to find maximum value in the column
            t=i; //row in which maximum value is present
        }
     }
    for(k=0;k<n+1;k++) //interchanging row with maximum value in the column and jth row(pivot row)
      {
        temp=a[j][k]; //storing a jth row element temporarily for interchanging
        a[j][k]=a[t][k];
        a[t][k]=temp;
      }
   for(i=0;i<n;i++)
     {
      if ( i>j) // make all the values below pivot to zero by row operations
       {
         constant=a[i][j]/a[j][j];
         for(k=0;k<n+1;k++)
          {
            a[i][k]=a[i][k]-constant*a[j][k]; //changes in values of component in a row due to the row operation
          }
        }

      }
   }
 FILE *fptr;
 fptr=fopen("output_3.txt","w"); // to write the out put in a file named output_3
 fprintf(fptr,"The final augmented matrix is:\n");
 for (int i=0;i<n;i++)
    {
        for (int j=0;j<n+1;j++)
         {
            fprintf(fptr,"%.2f\t",a[i][j]); // writing the final augmented matrix
         }
        fprintf(fptr,"\n");
    }
  if (a[n-1][n-1]!=0) // unique solution exist
  {
    x[n-1] = (a[n-1][n])/(a[n-1][n-1]); //the value of z
    for(i=n-1;i>=0;i--)
    {
      sum=0;
      for(j=i+1;j<n;j++)
      {
        sum=sum+a[i][j]*x[j];
      }
      x[i]=(a[i][n]-sum)/a[i][i]; //back substitution to get the values of y and x
    }
    fprintf(fptr,"The solution is:\n");
    for(i=0;i<n;i++)
    {
      fprintf(fptr,"%.2f\t",x[i]); // writing the solution
    }
 }
 fclose (fptr);
 printf("output is written on a file named output_3.txt"); // to know where the output is
 return 0;
}
