#include <stdlib.h>
#include <stdio.h>
#include <GL/freeglut.h>
#include <gtk/gtk.h>
#include <math.h>

#define ALPHA 0
#define CONTRA 1
#define MMSE 2
#define ARITHM 3
#define MTM 4

typedef struct {
   unsigned short int type;                 /* Magic identifier            */
   unsigned int size;                       /* File size in bytes          */
   unsigned short int reserved1, reserved2;
   unsigned int offset;                     /* Offset to image data, bytes */
} HEADER;

typedef struct {
   unsigned int size;               /* Header size in bytes      */
   int width,height;                /* Width and height of image */
   unsigned short int planes;       /* Number of colour planes   */
   unsigned short int bits;         /* Bits per pixel            */
   unsigned int compression;        /* Compression type          */
   unsigned int imagesize;          /* Image size in bytes       */
   int xresolution,yresolution;     /* Pixels per meter          */
   unsigned int ncolours;           /* Number of colours         */
   unsigned int importantcolours;   /* Important colours         */
} INFOHEADER;

typedef struct{
    unsigned char blue,green,red;
}COLOR;

FILE * fp_read,*fp_write;
HEADER* getHeader()
{
  HEADER *header;
  header=(HEADER*) malloc(sizeof(HEADER));

    fread(&header->type,sizeof(header->type),1,fp_read);
    if(header->type!=19778)
    {printf("Not a BMP Image\n");
     exit(0);}
    fread(&header->size,sizeof(header->size),1,fp_read);
    fread(&header->reserved1,sizeof(header->reserved1),1,fp_read);
     fread(&header->reserved2,sizeof(header->reserved2),1,fp_read);
     fread(&header->offset,sizeof(header->offset),1,fp_read);
     return header;
}

INFOHEADER* getInfoHeader()
{
    INFOHEADER *infoheader;
    infoheader=(INFOHEADER*) malloc(sizeof(INFOHEADER));
    fread(&infoheader->size,sizeof(infoheader->size),1,fp_read);
    fread(&infoheader->width,sizeof(infoheader->width),1,fp_read);
    fread(&infoheader->height,sizeof(infoheader->height),1,fp_read);
    fread(&infoheader->planes,sizeof(infoheader->planes),1,fp_read);
    fread(&infoheader->bits,sizeof(infoheader->bits),1,fp_read);
    if(infoheader->bits!=24)
    {printf("Not a 24 bit Image\n");
     exit(0);}
    fread(&infoheader->compression,sizeof(infoheader->compression),1,fp_read);
    fread(&infoheader->imagesize,sizeof(infoheader->imagesize),1,fp_read);
    fread(&infoheader->xresolution,sizeof(infoheader->xresolution),1,fp_read);
    fread(&infoheader->yresolution,sizeof(infoheader->yresolution),1,fp_read);
    fread(&infoheader->ncolours,sizeof(infoheader->ncolours),1,fp_read);
    fread(&infoheader->importantcolours,sizeof(infoheader->importantcolours),1,fp_read);



    return infoheader;
}

unsigned char *matrix;
unsigned char *new_mat;
unsigned char *noisy_mat;
unsigned char *gausi_mat;
HEADER *head;
INFOHEADER *infohead;
COLOR *color;

int window1,window2,window3,window4,window5,window6,window7;
int mat_size,i,j;
size_t num;
char* x;
int salt_flag=0,gauss_flag=0;

float val_msesp[5]={0,0,0,0,0};
float val_msegu[5]={0,0,0,0,0};
float val_psnrsp[5]={0,0,0,0,0};
float val_psnrgu[5]={0,0,0,0,0};
float val_snrsp[5]={0,0,0,0,0};
float val_snrgu[5]={0,0,0,0,0};

void mtm_N(void);
void arithm_mean_filter_N(void);
void mmse_N(void);
void alpha_N(void);
void contra_N(void);
void mtm_G(void);
void arithm_mean_filter_G(void);
void mmse_G(void);
void alpha_G(void);
void contra_G(void);
void read_matrix();

long getrand(int l)
{
    long num=0;
    int i;
    for( i=0;i<l;i++ )
    {
        num+=rand();

    }
    if(num>infohead->height*infohead->width)
    {
        num=num%(rand()+1);
    }
    return num;
}

void create_salt_pepper()
{
int i=0,j=0, k=0,z,zcount=0;
int limit;

limit=infohead->height*infohead->width*3;
int mat_size=infohead->height*infohead->width*3;
noisy_mat=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
int pixels=infohead->height*infohead->width;
int fivep=pixels*5/100;
long narr[fivep];
if(pixels>RAND_MAX)
{
 z=(pixels/RAND_MAX)+1;
}

for(i=0;i<fivep;i++)
{
    narr[i]=getrand(z);
    zcount++;
}

for(i=0;i<infohead->height;i++)
{
    for(j=0;j<infohead->width;j++)
    {
        noisy_mat[i*infohead->width*3+j*3]=matrix[i*infohead->width*3+j*3];
        noisy_mat[i*infohead->width*3+j*3+1]=matrix[i*infohead->width*3+j*3+1];
        noisy_mat[i*infohead->width*3+j*3+2]=matrix[i*infohead->width*3+j*3+2];
    }
}

for(i=0;i<zcount;i++)
{
    unsigned char p=0+(rand()%2)*255;
    noisy_mat[(int)(narr[i]/infohead->width)*infohead->width*3+(int)((narr[i]%infohead->width)-1)*3]=p;
    noisy_mat[(int)(narr[i]/infohead->width)*infohead->width*3+(int)((narr[i]%infohead->width)-1)*3+1]=p;
    noisy_mat[(int)(narr[i]/infohead->width)*infohead->width*3+(int)((narr[i]%infohead->width)-1)*3+2]=p;
}

}

void create_gauss_noise()
{
 int i=0,j=0, k=0,count=0,z,zcount=0;
 int limit;
 limit=infohead->height*infohead->width*3;
 int mat_size=infohead->height*infohead->width*3;
 unsigned char* noisy_mat;
 gausi_mat=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
 int mean=0;
 float var=0.10,vroot;

 vroot=pow(var,0.5);
 unsigned char temp;
 for(i=0;i<infohead->height;i++)
  {
     for(j=0;j<infohead->width;j++)
     {
         temp=vroot*(rand()%75)+mean;
         gausi_mat[i*infohead->width*3+j*3]=matrix[i*infohead->width*3+j*3]+temp;
         gausi_mat[i*infohead->width*3+j*3+1]=matrix[i*infohead->width*3+j*3+1]+temp;
         gausi_mat[i*infohead->width*3+j*3+2]=matrix[i*infohead->width*3+j*3+2]+temp;
     }
  }
}

long calc_mse(int n)
{
int count=0;
long long se=0;
double mse;
long num=(infohead->height-n)*(infohead->width-n);
count=0;

for(i=n/2;i<infohead->height-n/2;i++)
{
    for(j=n/2;j<infohead->width-n/2;j++)
    {
        se+=(matrix[i*(infohead->width*3)+j*3]-new_mat[count])*(matrix[i*(infohead->width*3)+j*3]-new_mat[count]);
        count+=3;
    }
}

mse=(double)se/num;
return mse;
}

float calc_psnr(double mse)
{
    //10log(base 10)Rsqr/mse
    int R=255,R2=pow(R,2);

    float v=R2/mse;
    float psnr=10*log10f(v);
    return psnr;
}

float calc_snr(int n)
{
int count;
long long ssum=0,nsum=0;
float snr,snrdb;
long num=(infohead->height-n)*(infohead->width-n);
count=0;

for(i=n/2;i<infohead->height-n/2;i++)
{
    for(j=n/2;j<infohead->width-n/2;j++)
    {
        ssum+=matrix[i*(infohead->width*3)+j*3];
        nsum+=abs(matrix[i*(infohead->width*3)+j*3]-new_mat[count]);
        count+=3;
    }
}
snr=(ssum/num)/(nsum/num);
snrdb=10*log10(snr);
return snrdb;
}

void display_noisy()
{
    int i=0,j=0, k=0,count=0;
int limit;
//int r=0,g=0,b=0;
//int p[2] = {0, 0};
glClear(GL_COLOR_BUFFER_BIT);
limit=infohead->height*infohead->width*3;

for(i=0;i<infohead->height;i++)     //-2 for mean filter,-4 for MTM
    {
        count=(i*((infohead->width*3)));    //-6 for mean filter,-12 for MTM
        for(j=0;j<infohead->width;j++)      //-2 for mean filter,-4 for MTM
        {
          k=3*j;
          glBegin(GL_POINTS);
          glColor3ub(noisy_mat[count+k+2],noisy_mat[count+k+1],noisy_mat[count+k]);
          glVertex2i(j,i);
          glEnd();
        }
    }

glFlush();
}

void display_gaussi()
{

    int i=0,j=0, k=0,count=0;
int limit;
//int r=0,g=0,b=0;
//int p[2] = {0, 0};
glClear(GL_COLOR_BUFFER_BIT);
limit=infohead->height*infohead->width*3;

for(i=0;i<infohead->height;i++)     //-2 for mean filter,-4 for MTM
    {
        count=(i*((infohead->width*3)));    //-6 for mean filter,-12 for MTM
        for(j=0;j<infohead->width;j++)      //-2 for mean filter,-4 for MTM
        {
          k=3*j;
          glBegin(GL_POINTS);
          glColor3ub(gausi_mat[count+k+2],gausi_mat[count+k+1],gausi_mat[count+k]);
          glVertex2i(j,i);
          glEnd();
        }
    }

glFlush();
}

void display(void) {
int i=0,j=0, k=0,count=0;
int limit;
glClear(GL_COLOR_BUFFER_BIT);
limit=infohead->height*infohead->width*3;

for(i=0;i<infohead->height;i++)     //-2 for mean filter,-4 for MTM
    {
        count=(i*((infohead->width*3)));    //-6 for mean filter,-12 for MTM
        for(j=0;j<infohead->width;j++)      //-2 for mean filter,-4 for MTM
        {
          k=3*j;
          glBegin(GL_POINTS);
          glColor3ub(matrix[count+k+2],matrix[count+k+1],matrix[count+k]);
          glVertex2i(j,i);
          glEnd();
        }
    }
glFlush();
}

 void display2(){
int i=0,j=0, k=0,count=0;
int limit;
char grey;
glClear(GL_COLOR_BUFFER_BIT);
limit=infohead->height*infohead->width*3;
arithm_mean_filter_N();

for(i=0;i<infohead->height-2;i++)     //-2 for mean filter,-4 for MTM and MMSE
    {
        count=(i*((infohead->width*3-6)));    //-6 for mean filter,-12 for MTM and MMSE
        for(j=0;j<infohead->width-2;j++)      //-2 for mean filter,-4 for MTM and MMSE
        {
          k=3*j;
          glBegin(GL_POINTS);
          glColor3ub(new_mat[count+k+2],new_mat[count+k+1],new_mat[count+k]);
          glVertex2i(j,i);
          glEnd();
        }
    }
glFlush();
 }

void display2_G(){
int i=0,j=0, k=0,count=0;
int limit;
char grey;

glClear(GL_COLOR_BUFFER_BIT);
limit=infohead->height*infohead->width*3;
arithm_mean_filter_G();

for(i=0;i<infohead->height-2;i++)     //-2 for mean filter,-4 for MTM and MMSE
    {
        count=(i*((infohead->width*3-6)));    //-6 for mean filter,-12 for MTM and MMSE
        for(j=0;j<infohead->width-2;j++)      //-2 for mean filter,-4 for MTM and MMSE
        {
          k=3*j;
          glBegin(GL_POINTS);
          glColor3ub(new_mat[count+k+2],new_mat[count+k+1],new_mat[count+k]);
          glVertex2i(j,i);
          glEnd();
        }
    }
glFlush();
 }

 void display3(void){
int i=0,j=0, k=0,count=0;
int limit;
char grey;
glClear(GL_COLOR_BUFFER_BIT);
limit=infohead->height*infohead->width*3;
mtm_N();
for(i=0;i<infohead->height-4;i++)     //-2 for mean filter,-4 for MTM and MMSE
    {
        count=(i*((infohead->width*3-12)));    //-6 for mean filter,-12 for MTM and MMSE
        for(j=0;j<infohead->width-4;j++)      //-2 for mean filter,-4 for MTM and MMSE
        {
          k=3*j;
          glBegin(GL_POINTS);
          glColor3ub(new_mat[count+k+2],new_mat[count+k+1],new_mat[count+k]);
          glVertex2i(j,i);
          glEnd();
        }
    }
glFlush();}

void display3_G(void){
int i=0,j=0, k=0,count=0;
int limit;
char grey;
//int r=0,g=0,b=0;
//int p[2] = {0, 0};
glClear(GL_COLOR_BUFFER_BIT);
limit=infohead->height*infohead->width*3;
mtm_G();
for(i=0;i<infohead->height-4;i++)     //-2 for mean filter,-4 for MTM and MMSE
    {
        count=(i*((infohead->width*3-12)));    //-6 for mean filter,-12 for MTM and MMSE
        for(j=0;j<infohead->width-4;j++)      //-2 for mean filter,-4 for MTM and MMSE
        {
          k=3*j;
          glBegin(GL_POINTS);
          glColor3ub(new_mat[count+k+2],new_mat[count+k+1],new_mat[count+k]);
      //    glColor3ub(matrix[count+k+2],matrix[count+k+1],matrix[count+k]);
          glVertex2i(j,i);
          glEnd();
        }
    }

glFlush();}

void display4(void){
int i=0,j=0, k=0,count=0;
glClear(GL_COLOR_BUFFER_BIT);
mmse_N();

for(i=0;i<infohead->height-4;i++)     //-2 for mean filter,-4 for MTM and MMSE
    {
        count=(i*((infohead->width*3-12)));    //-6 for mean filter,-12 for MTM and MMSE
        for(j=0;j<infohead->width-4;j++)      //-2 for mean filter,-4 for MTM and MMSE
        {
          k=3*j;
          glBegin(GL_POINTS);
          glColor3ub(new_mat[count+k+2],new_mat[count+k+1],new_mat[count+k]);
      //    glColor3ub(matrix[count+k+2],matrix[count+k+1],matrix[count+k]);
          glVertex2i(j,i);
          glEnd();
        }
    }

glFlush();}

void display4_G(void){
int i=0,j=0, k=0,count=0;
glClear(GL_COLOR_BUFFER_BIT);
mmse_G();

for(i=0;i<infohead->height-4;i++)     //-2 for mean filter,-4 for MTM and MMSE
    {
        count=(i*((infohead->width*3-12)));    //-6 for mean filter,-12 for MTM and MMSE
        for(j=0;j<infohead->width-4;j++)      //-2 for mean filter,-4 for MTM and MMSE
        {
          k=3*j;
          glBegin(GL_POINTS);
          glColor3ub(new_mat[count+k+2],new_mat[count+k+1],new_mat[count+k]);
      //    glColor3ub(matrix[count+k+2],matrix[count+k+1],matrix[count+k]);
          glVertex2i(j,i);
          glEnd();
        }
    }

glFlush();}

void display5(void){
int i=0,j=0, k=0,dip=0;

glClear(GL_COLOR_BUFFER_BIT);

alpha_N();
for(i=0;i<infohead->height-4;i++)     //-2 for mean filter,-4 for MTM and MMSE
    {
        dip=(i*((infohead->width*3-12)));    //-6 for mean filter,-12 for MTM and MMSE
        for(j=0;j<infohead->width-4;j++)      //-2 for mean filter,-4 for MTM and MMSE
        {
          k=3*j;
          glBegin(GL_POINTS);
          glColor3ub(new_mat[dip+k+2],new_mat[dip+k+1],new_mat[dip+k]);
      //    glColor3ub(matrix[count+k+2],matrix[count+k+1],matrix[count+k]);
          glVertex2i(j,i);
          glEnd();
        }
    }

glFlush();}

void display5_G(void){
int i=0,j=0, k=0,dip=0;

glClear(GL_COLOR_BUFFER_BIT);

alpha_G();
for(i=0;i<infohead->height-4;i++)     //-2 for mean filter,-4 for MTM and MMSE
    {
        dip=(i*((infohead->width*3-12)));    //-6 for mean filter,-12 for MTM and MMSE
        for(j=0;j<infohead->width-4;j++)      //-2 for mean filter,-4 for MTM and MMSE
        {
          k=3*j;
          glBegin(GL_POINTS);
          glColor3ub(new_mat[dip+k+2],new_mat[dip+k+1],new_mat[dip+k]);
      //    glColor3ub(matrix[count+k+2],matrix[count+k+1],matrix[count+k]);
          glVertex2i(j,i);
          glEnd();
        }
    }

glFlush();}

void display6(void){
int i=0,j=0, k=0,dip=0;

glClear(GL_COLOR_BUFFER_BIT);

contra_N();
for(i=0;i<infohead->height-4;i++)     //-2 for mean filter,-4 for MTM and MMSE
    {
        dip=(i*((infohead->width*3-12)));    //-6 for mean filter,-12 for MTM and MMSE
        for(j=0;j<infohead->width-4;j++)      //-2 for mean filter,-4 for MTM and MMSE
        {
          k=3*j;
          glBegin(GL_POINTS);
          glColor3ub(new_mat[dip+k+2],new_mat[dip+k+1],new_mat[dip+k]);
      //    glColor3ub(matrix[count+k+2],matrix[count+k+1],matrix[count+k]);
          glVertex2i(j,i);
          glEnd();
        }
    }

glFlush();}

void myinit(void) {
glClearColor(1.0, 1.0, 1.0, 1.0);
glColor3b(255, 0, 0);
glMatrixMode(GL_PROJECTION);
glLoadIdentity();
gluOrtho2D(0.0,infohead->width, 0.0,infohead->height);
glMatrixMode(GL_MODELVIEW);
}
int lb=-30, ub=30;
double cpoints=9999999;

void init (void)
{
    glClearColor(1,1,1,1);
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(0,1450,0,700);
    glMatrixMode(GL_MODELVIEW);
}

void display6_G(void){
int i=0,j=0, k=0,dip=0;

glClear(GL_COLOR_BUFFER_BIT);

contra_G();
for(i=0;i<infohead->height-4;i++)     //-2 for mean filter,-4 for MTM and MMSE
    {
        dip=(i*((infohead->width*3-12)));    //-6 for mean filter,-12 for MTM and MMSE
        for(j=0;j<infohead->width-4;j++)      //-2 for mean filter,-4 for MTM and MMSE
        {
          k=3*j;
          glBegin(GL_POINTS);
          glColor3ub(new_mat[dip+k+2],new_mat[dip+k+1],new_mat[dip+k]);
      //    glColor3ub(matrix[count+k+2],matrix[count+k+1],matrix[count+k]);
          glVertex2i(j,i);
          glEnd();
        }
    }

glFlush();}


void arithm_mean_filter_N(void){

int sum_col1=0,sum_col2=0,sum_col3=0,div=9,mat_size,count=0,i,j;

mat_size=(infohead->width-2)*(infohead->height-2)*3;

new_mat=(unsigned char*)calloc(mat_size,sizeof(unsigned char));

for(i=0;i<infohead->height-2;i++)
 {
        for(j=0;j<infohead->width-2;j++)
        {
            sum_col1=noisy_mat[i*infohead->width*3+j*3]+noisy_mat[(i+1)*infohead->width*3+j*3]+noisy_mat[(i+2)*infohead->width*3+j*3];
            sum_col2=noisy_mat[i*infohead->width*3+(j+1)*3]+noisy_mat[(i+1)*infohead->width*3+(j+1)*3]+noisy_mat[(i+2)*infohead->width*3+(j+1)*3];
            sum_col3=noisy_mat[i*infohead->width*3+(j+2)*3]+noisy_mat[(i+1)*infohead->width*3+(j+2)*3]+noisy_mat[(i+2)*infohead->width*3+(j+2)*3];

            new_mat[count]=(sum_col1+sum_col2+sum_col3)/div;
            count++;
            new_mat[count]=(sum_col1+sum_col2+sum_col3)/div;
            count++;
            new_mat[count]=(sum_col1+sum_col2+sum_col3)/div;
            count++;
        }
 }
double mse=calc_mse(2);
val_msesp[ARITHM]=mse;
float psnr=calc_psnr(mse);
val_psnrsp[ARITHM]=psnr;
float snr=calc_snr(2);
val_snrsp[ARITHM]=snr;

}

void arithm_mean_filter_G(void){

int sum_col1=0,sum_col2=0,sum_col3=0,div=9,mat_size,count=0,i,j;

mat_size=(infohead->width-2)*(infohead->height-2)*3;

new_mat=(unsigned char*)calloc(mat_size,sizeof(unsigned char));

for(i=0;i<infohead->height-2;i++)
 {
        for(j=0;j<infohead->width-2;j++)
        {
            sum_col1=gausi_mat[i*infohead->width*3+j*3]+gausi_mat[(i+1)*infohead->width*3+j*3]+gausi_mat[(i+2)*infohead->width*3+j*3];
            sum_col2=gausi_mat[i*infohead->width*3+(j+1)*3]+gausi_mat[(i+1)*infohead->width*3+(j+1)*3]+gausi_mat[(i+2)*infohead->width*3+(j+1)*3];
            sum_col3=gausi_mat[i*infohead->width*3+(j+2)*3]+gausi_mat[(i+1)*infohead->width*3+(j+2)*3]+gausi_mat[(i+2)*infohead->width*3+(j+2)*3];

            new_mat[count]=(sum_col1+sum_col2+sum_col3)/div;
            count++;
            new_mat[count]=(sum_col1+sum_col2+sum_col3)/div;
            count++;
            new_mat[count]=(sum_col1+sum_col2+sum_col3)/div;
            count++;
        }
 }
double mse=calc_mse(2);
val_msegu[ARITHM]=mse;
float psnr=calc_psnr(mse);
val_psnrgu[ARITHM]=psnr;
float snr=calc_snr(2);
val_snrgu[ARITHM]=snr;

}

void mtm_N(){

int mat_size,median_arr[9],m,n,t,temp,l,median,sum,grey,k,count=0;
float K=1.0,stddev=40.98;

mat_size=(infohead->width-4)*(infohead->height-4)*3;

new_mat=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
for(i=2;i<infohead->height-2;i++)
 {
        for(j=2;j<infohead->width-2;j++)
        {
            t=0;
            for(m=0;m<=2;m++)
            {
              for(n=0;n<=2;n++)
              {
                median_arr[t]=noisy_mat[(long)(i+m)*infohead->width*3+j*3+(n)*3];
                t++;
              }
            }

            for(k=1;k<9;k++)
            {
                temp=median_arr[k];
                l=k-1;
                while(l>=0 && median_arr[l]>temp)
                {
                   median_arr[l+1]=median_arr[l];
                   l=l-1;
                }
                median_arr[l+1]=temp;
            }

            median=median_arr[4];
            sum=0;
            t=0;
            for(m=-2;m<=2;m++)
            {
                for(n=-2;n<=2;n++)
                {
                    grey=noisy_mat[(i+m)*(infohead->width)*3+(j+n)*3];
                    if(grey>=median-(K*stddev))
                    {
                        if(grey<=median+(K*stddev))
                        {
                            sum+=grey;
                            t+=1;
                        }
                    }
                }
            }
             new_mat[count]=(unsigned char)((float)sum/(float)t);
             new_mat[count+1]=(unsigned char)((float)sum/(float)t);
             new_mat[count+2]=(unsigned char)((float)sum/(float)t);
             count+=3;

        }
 }
double mse=calc_mse(4);
float psnr=calc_psnr(mse);
float snr=calc_snr(4);

val_msesp[MTM]=mse;
val_psnrsp[MTM]=psnr;
val_snrsp[MTM]=snr;
}

void mtm_G(){

int mat_size,median_arr[9],m,n,t,temp,l,median,sum,grey,k,count=0;
float K=1.0,stddev=40.98;

mat_size=(infohead->width-2)*(infohead->height-2)*3;

new_mat=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
for(i=2;i<infohead->height-2;i++)
 {
        for(j=2;j<infohead->width-2;j++)
        {
            t=0;
            for(m=0;m<=2;m++)
            {
              for(n=0;n<=2;n++)
              {
                median_arr[t]=gausi_mat[(long)(i+m)*infohead->width*3+j*3+(n)*3];
                t++;
              }
            }

            for(k=1;k<9;k++)
            {
                temp=median_arr[k];
                l=k-1;
                while(l>=0 && median_arr[l]>temp)
                {
                   median_arr[l+1]=median_arr[l];
                   l=l-1;
                }
                median_arr[l+1]=temp;
            }

            median=median_arr[4];
            sum=0;
            t=0;
            for(m=-2;m<=2;m++)
            {
                for(n=-2;n<=2;n++)
                {
                    grey=gausi_mat[(i+m)*(infohead->width)*3+(j+n)*3];
                    if(grey>=median-(K*stddev))
                    {
                        if(grey<=median+(K*stddev))
                        {
                            sum+=grey;
                            t+=1;
                        }
                    }
                }
            }
             new_mat[count]=(unsigned char)((float)sum/(float)t);
             new_mat[count+1]=(unsigned char)((float)sum/(float)t);
             new_mat[count+2]=(unsigned char)((float)sum/(float)t);
             count+=3;

        }
 }
double mse=calc_mse(4);
float psnr=calc_psnr(mse);
float snr=calc_snr(4);

val_msegu[MTM]=mse;
val_psnrgu[MTM]=psnr;
val_snrgu[MTM]=snr;
}

 void mmse_N(){

int mat_size,median_arr[9],m,n,t,temp,l,sum,sumsq,g,count;
float meansq,mean,nv=20,lv;
n=5/2;
mat_size=(infohead->width-2)*(infohead->height-2)*3;

new_mat=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
unsigned char *p;
count=0;
p=noisy_mat;
for(i=n;i<infohead->height-n;i++)
{
    for(j=n;j<infohead->width-n;j++)
    {
        sum=0;
        sumsq=0;
        t=0;
         for(m=-n;m<=n;m++)
            {
              for(l=-n;l<=n;l++)
              {
                sum+=p[(long)(i+m)*infohead->width*3+j*3+(l)*3];
                sumsq+=(p[(long)(i+m)*infohead->width*3+j*3+(l)*3]*p[(long)(i+m)*infohead->width*3+j*3+(l)*3]);
                t++;
              }
            }
         meansq= (float)sumsq/t;
         mean=(float)sum/t;
         lv=meansq-mean*mean;
         if(lv==0.0)
           g=(int)(mean+0.5);
         else
           g=(int)((1-nv/lv)*p[i*infohead->width*3+j*3]+nv/lv*mean+0.5);
         if(g>255)
            g=255;
         else if(g<0)
            g=0;
        new_mat[count]=g;
        new_mat[count+1]=g;
        new_mat[count+2]=g;
        count+=3;
    }
}

double mse=calc_mse(4);
float psnr=calc_psnr(mse);
float snr=calc_snr(4);

val_msesp[MMSE]=mse;
val_psnrsp[MMSE]=psnr;
val_snrsp[MMSE]=snr;
}

void mmse_G(){

int mat_size,median_arr[9],m,n,t,temp,l,sum,sumsq,g,count;
float meansq,mean,nv=20,lv;
n=5/2;
mat_size=(infohead->width-2)*(infohead->height-2)*3;

new_mat=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
unsigned char *p;
count=0;
p=gausi_mat;
for(i=n;i<infohead->height-n;i++)
{
    for(j=n;j<infohead->width-n;j++)
    {
        sum=0;
        sumsq=0;
        t=0;
         for(m=-n;m<=n;m++)
            {
              for(l=-n;l<=n;l++)
              {
                sum+=p[(long)(i+m)*infohead->width*3+j*3+(l)*3];
                sumsq+=(p[(long)(i+m)*infohead->width*3+j*3+(l)*3]*p[(long)(i+m)*infohead->width*3+j*3+(l)*3]);
                t++;
              }
            }
         meansq= (float)sumsq/t;
         mean=(float)sum/t;
         lv=meansq-mean*mean;
         if(lv==0.0)
           g=(int)(mean+0.5);
         else
           g=(int)((1-nv/lv)*p[i*infohead->width*3+j*3]+nv/lv*mean+0.5);
         if(g>255)
            g=255;
         else if(g<0)
            g=0;
        new_mat[count]=g;
        new_mat[count+1]=g;
        new_mat[count+2]=g;
        count+=3;
    }
}

double mse=calc_mse(4);
float psnr=calc_psnr(mse);
float snr=calc_snr(4);

val_msegu[MMSE]=mse;
val_psnrgu[MMSE]=psnr;
val_snrgu[MMSE]=snr;
}

 void alpha_N(void){

int mat_size,m,n,z,x,y,i,j,arr[121],count,sum,temp,p=3;
float mean;

n=5;
mat_size=(infohead->width-2)*(infohead->height-2)*3;

new_mat=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
count=0;
for(y=n/2;y<infohead->height-n/2;y++){
    for(x=n/2;x<infohead->width-n/2;x++){
        z=0;
        for(j=-n/2;j<=n/2;j++){
            for(i=-n/2;i<=n/2;i++){
                arr[z]=noisy_mat[(y+j)*infohead->width*3+(x+i)*3];
                z++;
            }
        }

    for(j=1;j<n*n;j++)
            {
                temp=arr[j];
                i=j-1;
                while(i>=0 && arr[i]>temp)
                {
                   arr[i+1]=arr[i];
                   i=i-1;
                }
               arr[i+1]=temp;
            }

       sum=0;z=0;
       for(j=p;j<n*n-p;j++){
            sum+=arr[j];
            z++;
       }
mean=(sum/z)+0.5;
  new_mat[count]=(unsigned char)mean;
  new_mat[count+1]=(unsigned char)mean;
  new_mat[count+2]=(unsigned char)mean;
  count+=3;
  }
}

double mse=calc_mse(4);
float psnr=calc_psnr(mse);
float snr=calc_snr(4);

val_msesp[ALPHA]=mse;
val_psnrsp[ALPHA]=psnr;
val_snrsp[ALPHA]=snr;
}

void alpha_G(void){

int mat_size,m,n,z,x,y,i,j,arr[121],count,sum,temp,p=3;
float mean;

n=5;
mat_size=(infohead->width-2)*(infohead->height-2)*3;

new_mat=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
count=0;
for(y=n/2;y<infohead->height-n/2;y++){
    for(x=n/2;x<infohead->width-n/2;x++){
        z=0;
        for(j=-n/2;j<=n/2;j++){
            for(i=-n/2;i<=n/2;i++){
                arr[z]=gausi_mat[(y+j)*infohead->width*3+(x+i)*3];
                z++;
            }
        }

    for(j=1;j<n*n;j++)
            {
                temp=arr[j];
                i=j-1;
                while(i>=0 && arr[i]>temp)
                {
                   arr[i+1]=arr[i];
                   i=i-1;
                }
               arr[i+1]=temp;
            }

       sum=0;z=0;
       for(j=p;j<n*n-p;j++){
            sum+=arr[j];
            z++;
       }
mean=(sum/z)+0.5;
  new_mat[count]=(unsigned char)mean;
  new_mat[count+1]=(unsigned char)mean;
  new_mat[count+2]=(unsigned char)mean;
  count+=3;
  }
}

double mse=calc_mse(4);
float psnr=calc_psnr(mse);
float snr=calc_snr(4);

val_msegu[ALPHA]=mse;
val_psnrgu[ALPHA]=psnr;
val_snrgu[ALPHA]=snr;
}

void contra_N(void){

int mat_size,m,n,z,x,y,i,j,count,temp;
float sum,sum1,p=0.059;

n=5;
mat_size=(infohead->width-4)*(infohead->height-4)*3;

new_mat=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
count=0;
for(y=n/2;y<infohead->height-n/2;y++){
    for(x=n/2;x<infohead->width-n/2;x++){
        z=0;
        sum=0.0;
        sum1=0.0;
        for(j=-n/2;j<=n/2;j++){
            for(i=-n/2;i<=n/2;i++){
             temp=noisy_mat[(y+j)*infohead->width*3+(x+i)*3];
        if(temp==0&&p<0){
             z=1;}
            else{
                sum+=pow((double)temp,(double)p+1);
                sum1+=pow((double)temp,(double)p);
            }

            }
        }
   if(z==1){
    new_mat[count]=0;
    new_mat[count+1]=0;
    new_mat[count+2]=0;
    count+=3;
   }
   else{
    if(sum1==0.0){
       new_mat[count]=0;
       new_mat[count+1]=0;
       new_mat[count+2]=0;
       count+=3;
    }
    else{
        m=(int)sum/sum1;
        if(m>255){
            m=255;
        }
       new_mat[count]=m;
       new_mat[count+1]=m;
       new_mat[count+2]=m;
       count+=3;
    }
   }
  }
 }
 double mse=calc_mse(4);
 float psnr=calc_psnr(mse);
 float snr=calc_snr(4);

val_msesp[CONTRA]=mse;
val_psnrsp[CONTRA]=psnr;
val_snrsp[CONTRA]=snr;

}

void contra_G(void){

int mat_size,m,n,z,x,y,i,j,count,temp;
float sum,sum1,p=0.059;

n=5;
mat_size=(infohead->width-4)*(infohead->height-4)*3;

new_mat=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
count=0;
for(y=n/2;y<infohead->height-n/2;y++){
    for(x=n/2;x<infohead->width-n/2;x++){
        z=0;
        sum=0.0;
        sum1=0.0;
        for(j=-n/2;j<=n/2;j++){
            for(i=-n/2;i<=n/2;i++){
             temp=gausi_mat[(y+j)*infohead->width*3+(x+i)*3];
        if(temp==0&&p<0){
             z=1;}
            else{
                sum+=pow((double)temp,(double)p+1);
                sum1+=pow((double)temp,(double)p);
            }

            }
        }
   if(z==1){
    new_mat[count]=0;
    new_mat[count+1]=0;
    new_mat[count+2]=0;
    count+=3;
   }
   else{
    if(sum1==0.0){
       new_mat[count]=0;
       new_mat[count+1]=0;
       new_mat[count+2]=0;
       count+=3;
    }
    else{
        m=(int)sum/sum1;
        if(m>255){
            m=255;
        }
       new_mat[count]=m;
       new_mat[count+1]=m;
       new_mat[count+2]=m;
       count+=3;
    }
   }
  }
 }
 double mse=calc_mse(4);
 float psnr=calc_psnr(mse);
 float snr=calc_snr(4);

val_msegu[CONTRA]=mse;
val_psnrgu[CONTRA]=psnr;
val_snrgu[CONTRA]=snr;

}

static void mmsed(GtkWidget* widget, gpointer data)
{
	if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget)))
	{
		   fp_read = fopen (x,"rb+");
    if(fp_read==NULL)
    {printf("file not found \nexiting....");
     exit(0);}
     head=getHeader();
     infohead=getInfoHeader();
     color=(COLOR*) malloc(sizeof(COLOR));
     mat_size=(infohead->width*3)*infohead->height;
     matrix=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
     read_matrix();

    char** argv_m = 0;
    int argc_m = 0;
    glutInit(&argc_m,argv_m);
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize (infohead->width, infohead->height);
    glutInitWindowPosition (0, 0);
    window1=glutCreateWindow ("Original");
    myinit();
    glutDisplayFunc (display);

    if(salt_flag==0)
    {create_salt_pepper();
     salt_flag++;}

    glutInitWindowPosition (infohead->width+15, 0);
    window4=glutCreateWindow ("Noisy_MTM");
    myinit();
    glutDisplayFunc (display_noisy);

    if(gauss_flag==0)
    {create_gauss_noise();
     gauss_flag++;}

     glutInitWindowPosition (2*infohead->width+30, 0);
    window5=glutCreateWindow ("Gaussy_MMSE");
    myinit();
    glutDisplayFunc (display_gaussi);

    glutInitWindowPosition (3*infohead->width+35, 0);
    window5=glutCreateWindow ("NoiseR_MMSE");
    myinit();
    glutDisplayFunc (display4);

    glutInitWindowPosition (4*infohead->width+50, 0);
    window5=glutCreateWindow ("GaussR_MMSE");
    myinit();
    glutDisplayFunc (display4_G);

    glutMainLoop();
	}else
	{
		g_print(" Grey Code Deactivated\n");
	}
}

static void contrad(GtkWidget* widget, gpointer data)
{
	if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget)))
	{
		   fp_read = fopen (x,"rb+");
    if(fp_read==NULL)
    {printf("file not found \nexiting....");
     exit(0);}
     head=getHeader();
     infohead=getInfoHeader();
     color=(COLOR*) malloc(sizeof(COLOR));
     mat_size=(infohead->width*3)*infohead->height;
     matrix=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
     read_matrix();

char** argv_m = 0;
int argc_m = 0;
glutInit(&argc_m,argv_m);
glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS);
glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize (infohead->width, infohead->height);
    glutInitWindowPosition (0, 0);
    window1=glutCreateWindow ("Original");
    myinit();
    glutDisplayFunc (display);

    if(salt_flag==0)
    {create_salt_pepper();
     salt_flag++;}

    glutInitWindowPosition (infohead->width+15, 0);
    window4=glutCreateWindow ("Noisy_Contra");
    myinit();
    glutDisplayFunc (display_noisy);

    if(gauss_flag==0)
    {create_gauss_noise();
     gauss_flag++;}

    glutInitWindowPosition (2*infohead->width+30, 0);
    window5=glutCreateWindow ("Gaussy_Contra");
    myinit();
    glutDisplayFunc (display_gaussi);

    glutInitWindowPosition (3*infohead->width+35, 0);
    window5=glutCreateWindow ("NoiseR_Contra");
    myinit();
    glutDisplayFunc (display6);

    glutInitWindowPosition (4*infohead->width+50, 0);
    window5=glutCreateWindow ("GaussR_Contra");
    myinit();
    glutDisplayFunc (display6_G);

    glutMainLoop();
		}else
	{
		g_print(" Grey Code Deactivated\n");
	}
}

static void mtmd(GtkWidget* widget, gpointer data)
{
	if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget)))
	{
		   fp_read = fopen (x,"rb+");
    if(fp_read==NULL)
    {printf("file not found \nexiting....");
     exit(0);}
     head=getHeader();
     infohead=getInfoHeader();
     color=(COLOR*) malloc(sizeof(COLOR));
     mat_size=(infohead->width*3)*infohead->height;
     matrix=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
     read_matrix();


char** argv_l = 0;
int argc_l = 0;
glutInit(&argc_l,argv_l);
glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS);
glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize (infohead->width, infohead->height);
    glutInitWindowPosition (0, 0);
    window1=glutCreateWindow ("Original");
    myinit();
    glutDisplayFunc (display);

    if(salt_flag==0)
    {create_salt_pepper();
     salt_flag++;}

    glutInitWindowPosition (infohead->width+15, 0);
    window4=glutCreateWindow ("Noisy_MTM");
    myinit();
    glutDisplayFunc (display_noisy);

      if(gauss_flag==0)
    {create_gauss_noise();
     gauss_flag++;}

     glutInitWindowPosition (2*infohead->width+30, 0);
    window5=glutCreateWindow ("Gaussy_MTM");
    myinit();
    glutDisplayFunc (display_gaussi);

    glutInitWindowPosition (3*infohead->width+35, 0);
    window5=glutCreateWindow ("NoiseR_MTM");
    myinit();
    glutDisplayFunc (display3);

    glutInitWindowPosition (4*infohead->width+50, 0);
    window5=glutCreateWindow ("GaussR_MTM");
    myinit();
    glutDisplayFunc (display3_G);

    glutMainLoop();
		}else
	{
		g_print(" Grey Code Deactivated\n");
	}
}

static void arithmd(GtkWidget* widget, gpointer data)
{
	if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget)))
	{
		   fp_read = fopen (x,"rb+");
    if(fp_read==NULL)
    {printf("file not found \nexiting....");
     exit(0);}
     head=getHeader();
     infohead=getInfoHeader();
     color=(COLOR*) malloc(sizeof(COLOR));
     mat_size=(infohead->width*3)*infohead->height;
     matrix=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
     read_matrix();

char** argv_ = 0;
int argc_ = 0;
glutInit(&argc_,argv_);
glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize (infohead->width, infohead->height);
    glutInitWindowPosition (0, 0);
    window1=glutCreateWindow ("Original");
    myinit();
    glutDisplayFunc (display);

    if(salt_flag==0)
    {create_salt_pepper();
     salt_flag++;}

    glutInitWindowPosition (infohead->width+15, 0);
    window4=glutCreateWindow ("Noisy_Arithmetic");
    myinit();
    glutDisplayFunc (display_noisy);
    glutInitWindowPosition (2*infohead->width+30, 0);

      if(gauss_flag==0)
    {create_gauss_noise();
     gauss_flag++;}

     glutInitWindowPosition (2*infohead->width+30, 0);
    window5=glutCreateWindow ("Gaussy_Arithmetic");
    myinit();
    glutDisplayFunc (display_gaussi);

    glutInitWindowPosition (3*infohead->width+35, 0);
    window5=glutCreateWindow ("NoiseR_Arithmetic");
    myinit();
    glutDisplayFunc (display2);

    glutInitWindowPosition (4*infohead->width+50, 0);
    window5=glutCreateWindow ("GaussR_Arithmetic");
    myinit();
    glutDisplayFunc (display2_G);

    glutMainLoop();
    }else
	{
		g_print(" Grey Code Deactivated\n");
	}
}

static void alphad(GtkWidget* widget, gpointer data)
{
	if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget)))
	{
		   fp_read = fopen (x,"rb+");
   //fp_write= fopen("op.bmp","wb");
    if(fp_read==NULL)
    {printf("file not found \nexiting....");
     exit(0);}
     head=getHeader();
     infohead=getInfoHeader();
     color=(COLOR*) malloc(sizeof(COLOR));
     mat_size=(infohead->width*3)*infohead->height;
     matrix=(unsigned char*)calloc(mat_size,sizeof(unsigned char));
     read_matrix();

    char** argv_ = 0;
    int argc_ = 0;
    glutInit(&argc_,argv_);
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize (infohead->width, infohead->height);
    glutInitWindowPosition (0, 0);
    window1=glutCreateWindow ("Original");
    myinit();
    glutDisplayFunc (display);

    if(salt_flag==0)
    {create_salt_pepper();
     salt_flag++;}

    glutInitWindowPosition (infohead->width+15, 0);
    window4=glutCreateWindow ("Noisy_Alpha");
    myinit();
    glutDisplayFunc (display_noisy);
    glutInitWindowPosition (2*infohead->width+30, 0);

   if(gauss_flag==0)
    {create_gauss_noise();
     gauss_flag++;}

     glutInitWindowPosition (2*infohead->width+30, 0);
     window5=glutCreateWindow ("Gaussy_Alpha");
     myinit();
     glutDisplayFunc (display_gaussi);

     glutInitWindowPosition (3*infohead->width+35, 0);
     glutCreateWindow ("NoiseR_Alpha");
     myinit();
     glutDisplayFunc (display5);

     glutInitWindowPosition (4*infohead->width+50, 0);
     glutCreateWindow ("GaussR_Alpha");
     myinit();
     glutDisplayFunc (display5_G);

    glutMainLoop();
    }else
	{
		g_print("Toggle button Deactivated\n");
	}
}

void read_matrix()
{
    unsigned char temp,grey;
    int count_row=0,count_col=0;
    int count=0,i,j;
     while (count_row<infohead->height)
      {
        while(count_col<infohead->width)
        {
            fread(&color->blue,sizeof(color->blue),1,fp_read);
            matrix[count]=color->blue;
            count++;
            fread(&color->green,sizeof(color->green),1,fp_read);
            matrix[count]=color->green;
            count++;
            fread(&color->red,sizeof(color->red),1,fp_read);
            matrix[count]=color->red;

            grey=(color->red+color->blue+color->green)/3;
            matrix[count-2]=grey;
            matrix[count-1]=grey;
            matrix[count]=grey;

            count_col++;
            count++;
        }

            for(i=1;i<=((infohead->width)%4)%4;i++)
               {
                  fread(&temp,sizeof(temp),1,fp_read);
               }

        count_row++;
        count_col=0;
      }
}

float max_val(int x){
   int i;
   float temp;
   switch(x)
   {
   case 0:
       temp=val_msesp[0];
       for(i=1;i<5;i++)
       {
           if(val_msesp[i]>temp)
            temp=val_msesp[i];
       }
       return temp;
   break;
   case 1:
       temp=val_msegu[0];
       for(i=1;i<5;i++)
       {
           if(val_msegu[i]>temp)
            temp=val_msegu[i];
       }
       return temp;
    break;
   case 2:
       temp=val_psnrsp[0];
       for(i=1;i<5;i++)
       {
           if(val_psnrsp[i]>temp)
            temp=val_psnrsp[i];
       }
       return temp;
   break;
   case 3:
       temp=val_psnrgu[0];
       for(i=1;i<5;i++)
       {
           if(val_psnrgu[i]>temp)
            temp=val_psnrgu[i];
       }
       return temp;
   break;
   case 4:
       temp=val_snrsp[0];
       for(i=1;i<5;i++)
       {
           if(val_snrsp[i]>temp)
            temp=val_snrsp[i];
       }
       return temp;
   break;
   case 5:
       temp=val_snrgu[0];
       for(i=1;i<5;i++)
       {
           if(val_snrgu[i]>temp)
            temp=val_snrgu[i];
       }
       return temp;
   break;
   }
}

void graphfunct2D(void)
{
    double dx, x, y;
    glClear(GL_COLOR_BUFFER_BIT);
    int i;

        glBegin(GL_POINTS);
        y=5;
        glColor3ub(0,0,0);

        float msesp_max=max_val(0);
        float msegu_max=max_val(1);
        float psnrsp_max=max_val(2);
        float psnrgu_max=max_val(3);
        float snrsp_max=max_val(4);
        float snrgu_max=max_val(5);

        for(y=0;y<700;y++){
            for(x=0;x<1450;x++){

                 if(x==1200 || x==1201 || x==1202)
                 {
                   glColor3ub(0,0,0);
                   glVertex2i(x,y);
                 }
                if(y==350 && x<1200){
                   glColor3ub(0,0,0);
                   glVertex2i(x,y);
                }
                if(x==400 || x==800){
                   glColor3ub(0,0,0);
                   glVertex2i(x,y);
                }
                // End of Screen Partition, Start of pallet
                if((x>1210 && x<1230) && (y>60 && y<80))
                {
                    glColor3ub(217,244,66);
                    glVertex2i(x,y);
                }
                if((x>1210 && x<1230) && (y>200 && y<220))
                {
                   glColor3ub(10,5,55);
                    glVertex2i(x,y);
                }
                 if((x>1210 && x<1230) && (y>340 && y<360))
                {
                   glColor3ub(200,56,55);
                    glVertex2i(x,y);
                }
                 if((x>1210 && x<1230) && (y>480 && y<500))
                {
                   glColor3ub(0,156,55);
                    glVertex2i(x,y);
                }
                 if((x>1210 && x<1230) && (y>620 && y<640))
                {
                   glColor3ub(0,56,255);
                    glVertex2i(x,y);
                }
                //End of Screen Pallet
                if((x>=0 && x<400)&&(y>350 && y<=700)){

             if(!(x<0+25 || x>375) && y==375)
              { glColor3ub(0,0,0);
                  glVertex2i(x,y);   }

          if(!(x<25 || x>375) && y==376)
              {  glColor3ub(0,0,0);
                  glVertex2i(x,y);      }

           if(!(y<350+25 || y>675) && x==25)
              { glColor3ub(0,0,0);
                   glVertex2i(x,y);      }

           if(!(y<350+25 || y>675) && x==26)  //end of axis
              {  glColor3ub(0,0,0);
                  glVertex2i(x,y);      }

        if(val_msesp[ALPHA]!=0)
           if((x>35 && x<93) && (y<((300/log(msesp_max))*log(val_msesp[ALPHA]))+350 && y>376))
           {
               glColor3ub(0,56,255);
               glVertex2i(x,y);
           }

         if(val_msesp[CONTRA]!=0)
           if((x>103 && x<161) && (y<(300/log(msesp_max))*log(val_msesp[CONTRA])+350 && y>376))
           {
               glColor3ub(0,156,55);
               glVertex2i(x,y);
           }

         if(val_msesp[MMSE]!=0)
           if((x>171 && x<229) && (y<(300/log(msesp_max))*log(val_msesp[MMSE])+350 && y>376))
           {
               glColor3ub(200,56,55);
               glVertex2i(x,y);
           }

           if(val_msesp[ARITHM]!=0)
           if((x>239 && x<297) && (y<(300/log(msesp_max))*log(val_msesp[ARITHM])+350 && y>376))
           {
               glColor3ub(10,5,55);
               glVertex2i(x,y);
           }

          if(val_msesp[MTM]!=0)
           if((x>307 && x<365) && (y<(300/log(msesp_max))*log(val_msesp[MTM])+350 && y>376))
           {
               glColor3ub(217,244,66);
               glVertex2i(x,y);
           }
                }
                // 2 quadrant
                if((x>=400 && x<800)&&(y>350 && y<=700)){

             if(!(x<425 || x>775) && y==375)
              { glColor3ub(0,0,0);
                  glVertex2i(x,y);   }

          if(!(x<425 || x>775) && y==376)
              {  glColor3ub(0,0,0);
                  glVertex2i(x,y);      }

           if(!(y<350+25 || y>675) && x==425)
              { glColor3ub(0,0,0);
                   glVertex2i(x,y);      }

           if(!(y<350+25 || y>675) && x==426)  //end of axis
              {  glColor3ub(0,0,0);
                  glVertex2i(x,y);      }

        if(val_psnrsp[ALPHA]!=0)
           if((x>435 && x<493) && (y<(300/psnrsp_max)*val_psnrsp[ALPHA]+350 && y>376))
           {
               glColor3ub(0,56,255);
               glVertex2i(x,y);
           }

         if(val_psnrsp[CONTRA]!=0)
           if((x>503 && x<561) && (y<(300/psnrsp_max)*val_psnrsp[CONTRA]+350 && y>376))
           {
               glColor3ub(0,156,55);
               glVertex2i(x,y);
           }

         if(val_psnrsp[MMSE]!=0)
           if((x>571 && x<629) && (y<(300/psnrsp_max)*val_psnrsp[MMSE]+350 && y>376))
           {
               glColor3ub(200,56,55);
               glVertex2i(x,y);
           }

           if(val_psnrsp[ARITHM]!=0)
           if((x>639 && x<697) && (y<(300/psnrsp_max)*val_psnrsp[ARITHM]+350 && y>376))
           {
               glColor3ub(10,5,55);
               glVertex2i(x,y);
           }

          if(val_psnrsp[MTM]!=0)
           if((x>707 && x<765) && (y<(300/psnrsp_max)*val_psnrsp[MTM]+350 && y>376))
           {
               glColor3ub(217,244,66);
               glVertex2i(x,y);
           }
                }
                // 3 quadrant
                if((x>=800 && x<1200)&&(y>350 && y<=700)){

             if(!(x<825 || x>1175) && y==375)
              { glColor3ub(0,0,0);
                  glVertex2i(x,y);   }

          if(!(x<825 || x>1175) && y==376)
              {  glColor3ub(0,0,0);
                  glVertex2i(x,y);      }

           if(!(y<350+25 || y>675) && x==825)
              { glColor3ub(0,0,0);
                   glVertex2i(x,y);      }

           if(!(y<350+25 || y>675) && x==826)  //end of axis
              {  glColor3ub(0,0,0);
                  glVertex2i(x,y);      }

        if(val_snrsp[ALPHA]!=0)
           if((x>835 && x<893) && (y<(300/snrsp_max)*val_snrsp[ALPHA]+350 && y>376))
           {
               glColor3ub(0,56,255);
               glVertex2i(x,y);
           }

         if(val_snrsp[CONTRA]!=0)
           if((x>903 && x<961) && (y<(300/snrsp_max)*val_snrsp[CONTRA]+350 && y>376))
           {
               glColor3ub(0,156,55);
               glVertex2i(x,y);
           }

         if(val_snrsp[MMSE]!=0)
           if((x>971 && x<1029) && (y<(300/snrsp_max)*val_snrsp[MMSE]+350 && y>376))
           {
               glColor3ub(200,56,55);
               glVertex2i(x,y);
           }

           if(val_snrsp[ARITHM]!=0)
           if((x>1039 && x<1097) && (y<(300/snrsp_max)*val_snrsp[ARITHM]+350 && y>376))
           {
               glColor3ub(10,5,55);
               glVertex2i(x,y);
           }

          if(val_snrsp[MTM]!=0)
           if((x>1107 && x<1165) && (y<(300/snrsp_max)*val_snrsp[MTM]+350 && y>376))
           {
               glColor3ub(217,244,66);
               glVertex2i(x,y);
           }
                }

                // 4th quadrant
                if((x>=0 && x<400)&&(y>25 && y<=325)){

             if(!(x<25 || x>375) && y==26)
              { glColor3ub(0,0,0);
                  glVertex2i(x,y);   }

          if(!(x<25 || x>375) && y==27)
              {  glColor3ub(0,0,0);
                  glVertex2i(x,y);      }

           if(!(y<25 || y>375) && x==25)
              { glColor3ub(0,0,0);
                   glVertex2i(x,y);      }

           if(!(y<25 || y>375) && x==26)  //end of axis
              {  glColor3ub(0,0,0);
                  glVertex2i(x,y);      }

        if(val_msegu[ALPHA]!=0)
           if((x>35 && x<93) && (y<(300/log(msegu_max))*log(val_msegu[ALPHA]) && y>27))
           {
               glColor3ub(0,56,255);
               glVertex2i(x,y);
           }

         if(val_msegu[CONTRA]!=0)
           if((x>103 && x<161) && (y<(300/log(msegu_max))*log(val_msegu[CONTRA]) && y>27))
           {
               glColor3ub(0,156,55);
               glVertex2i(x,y);
           }

         if(val_msegu[MMSE]!=0)
           if((x>171 && x<229) && (y<(300/log(msegu_max))*log(val_msegu[MMSE]) && y>27))
           {
               glColor3ub(200,56,55);
               glVertex2i(x,y);
           }

           if(val_msegu[ARITHM]!=0)
           if((x>239 && x<297) && (y<(300/log(msegu_max))*log(val_msegu[ARITHM])&& y>27))
           {
               glColor3ub(10,5,55);
               glVertex2i(x,y);
           }

          if(val_msegu[MTM]!=0)
           if((x>307 && x<365) && (y<(300/log(msegu_max))*log(val_msegu[MTM]) && y>27))
           {
               glColor3ub(217,244,66);
               glVertex2i(x,y);
           }
                }
                // 5 quadrant
                if((x>=400 && x<800)&&(y>25 && y<=325)){

             if(!(x<425 || x>775) && y==26)
              { glColor3ub(0,0,0);
                  glVertex2i(x,y);   }

          if(!(x<425 || x>775) && y==27)
              {  glColor3ub(0,0,0);
                  glVertex2i(x,y);      }

           if(!(y<25 || y>375) && x==425)
              { glColor3ub(0,0,0);
                   glVertex2i(x,y);      }

           if(!(y<25 || y>375) && x==426)  //end of axis
              {  glColor3ub(0,0,0);
                  glVertex2i(x,y);      }

        if(val_psnrgu[ALPHA]!=0)
           if((x>435 && x<493) && (y<(300/psnrgu_max)*val_psnrgu[ALPHA] && y>27))
           {
               glColor3ub(0,56,255);
               glVertex2i(x,y);
           }

         if(val_psnrgu[CONTRA]!=0)
           if((x>503 && x<561) && (y<(300/psnrgu_max)*val_psnrgu[CONTRA] && y>27))
           {
               glColor3ub(0,156,55);
               glVertex2i(x,y);
           }

         if(val_psnrgu[MMSE]!=0)
           if((x>571 && x<629) && (y<(300/psnrgu_max)*val_psnrgu[MMSE] && y>27))
           {
               glColor3ub(200,56,55);
               glVertex2i(x,y);
           }

           if(val_psnrgu[ARITHM]!=0)
           if((x>639 && x<697) && (y<(300/psnrgu_max)*val_psnrgu[ARITHM] && y>27))
           {
               glColor3ub(10,5,55);
               glVertex2i(x,y);
           }

          if(val_psnrgu[MTM]!=0)
           if((x>707 && x<765) && (y<(300/psnrgu_max)*val_psnrgu[MTM] && y>27))
           {
               glColor3ub(217,244,66);
               glVertex2i(x,y);
           }
                }

                // 6 quadrant
                if((x>=800 && x<1200)&&(y>25 && y<=325)){

             if(!(x<825 || x>1175) && y==26)
              { glColor3ub(0,0,0);
                  glVertex2i(x,y);   }

          if(!(x<825 || x>1175) && y==27)
              {  glColor3ub(0,0,0);
                  glVertex2i(x,y);      }

           if(!(y<25 || y>375) && x==825)
              { glColor3ub(0,0,0);
                   glVertex2i(x,y);      }

           if(!(y<25 || y>375) && x==826)  //end of axis
              {  glColor3ub(0,0,0);
                  glVertex2i(x,y);      }

        if(val_snrgu[ALPHA]!=0)
           if((x>835 && x<893) && (y<(300/snrgu_max)*val_snrgu[ALPHA] && y>27))
           {
               glColor3ub(0,56,255);
               glVertex2i(x,y);
           }

         if(val_snrgu[CONTRA]!=0)
           if((x>903 && x<961) && (y<(300/snrgu_max)*val_snrgu[CONTRA] && y>27))
           {
               glColor3ub(0,156,55);
               glVertex2i(x,y);
           }

         if(val_snrgu[MMSE]!=0)
           if((x>971 && x<1029) && (y<(300/snrgu_max)*val_snrgu[MMSE] && y>27))
           {
               glColor3ub(200,56,55);
               glVertex2i(x,y);
           }

           if(val_snrgu[ARITHM]!=0)
           if((x>1039 && x<1097) && (y<(300/snrgu_max)*val_snrgu[ARITHM] && y>27))
           {
               glColor3ub(10,5,55);
               glVertex2i(x,y);
           }

          if(val_snrgu[MTM]!=0)
           if((x>1107 && x<1165) && (y<(300/snrgu_max)*val_snrgu[MTM] && y>27))
           {
               glColor3ub(217,244,66);
               glVertex2i(x,y);
           }
                }
            }
        }

        glEnd();
      //pallet text start
    glRasterPos2i(1240,65);
	glutBitmapString(GLUT_BITMAP_HELVETICA_18,"MTM");
	glRasterPos2i(1240,205);
	glutBitmapString(GLUT_BITMAP_HELVETICA_18,"Arithmetic Mean");
    glRasterPos2i(1240,345);
	glutBitmapString(GLUT_BITMAP_HELVETICA_18,"MMSE");
	glRasterPos2i(1240,485);
	glutBitmapString(GLUT_BITMAP_HELVETICA_18,"Contra Harmonic");
	glRasterPos2i(1240,625);
	glutBitmapString(GLUT_BITMAP_HELVETICA_18,"Alpha Trimmed Mean");
        //end of pallet text,start of bar values
    glRasterPos2i(25, 360);
	glutBitmapString(GLUT_BITMAP_HELVETICA_12,"MSE value for Salt and Pepper Noise (Image Quality ~ 1/MSE)");
	if(val_msesp[ALPHA]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_msesp[ALPHA]);
            glRasterPos2i(35,(300/log(msesp_max))*log(val_msesp[ALPHA])+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
     if(val_msesp[CONTRA]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_msesp[CONTRA]);
            glRasterPos2i(103,(300/log(msesp_max))*log(val_msesp[CONTRA])+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
     if(val_msesp[MMSE]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_msesp[MMSE]);
            glRasterPos2i(171,(300/log(msesp_max))*log(val_msesp[MMSE])+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
     if(val_msesp[ARITHM]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_msesp[ARITHM]);
            glRasterPos2i(239,(300/log(msesp_max))*log(val_msesp[ARITHM])+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
    if(val_msesp[MTM]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_msesp[MTM]);
            glRasterPos2i(307,(300/log(msesp_max))*log(val_msesp[MTM])+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
	glRasterPos2i(425, 360);
	glutBitmapString(GLUT_BITMAP_HELVETICA_12,"PSNR value for Salt and Pepper Noise (Image Quality ~ PSNR)");
	if(val_psnrsp[ALPHA]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_psnrsp[ALPHA]);
            glRasterPos2i(435,(300/psnrsp_max)*val_psnrsp[ALPHA]+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
     if(val_psnrsp[CONTRA]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_psnrsp[CONTRA]);
            glRasterPos2i(503,(300/psnrsp_max)*val_psnrsp[CONTRA]+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
      if(val_psnrsp[MMSE]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_psnrsp[MMSE]);
            glRasterPos2i(571,(300/psnrsp_max)*val_psnrsp[MMSE]+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
       if(val_psnrsp[ARITHM]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_psnrsp[ARITHM]);
            glRasterPos2i(639,(300/psnrsp_max)*val_psnrsp[ARITHM]+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
       if(val_psnrsp[MTM]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_psnrsp[MTM]);
            glRasterPos2i(707,(300/psnrsp_max)*val_psnrsp[MTM]+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
	glRasterPos2i(825, 360);
	glutBitmapString(GLUT_BITMAP_HELVETICA_12,"SNR value for Salt and Pepper Noise (Image Quality ~ SNR)");
		if(val_snrsp[ALPHA]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_snrsp[ALPHA]);
            glRasterPos2i(835,(300/snrsp_max)*val_snrsp[ALPHA]+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
     if(val_snrsp[CONTRA]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_snrsp[CONTRA]);
            glRasterPos2i(903,(300/snrsp_max)*val_snrsp[CONTRA]+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
      if(val_snrsp[MMSE]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_snrsp[MMSE]);
            glRasterPos2i(971,(300/snrsp_max)*val_snrsp[MMSE]+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
       if(val_snrsp[ARITHM]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_snrsp[ARITHM]);
            glRasterPos2i(1039,(300/snrsp_max)*val_snrsp[ARITHM]+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
       if(val_snrsp[MTM]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_snrsp[MTM]);
            glRasterPos2i(1107,(300/snrsp_max)*val_snrsp[MTM]+350+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
	glRasterPos2i(25, 10);
	glutBitmapString(GLUT_BITMAP_HELVETICA_12,"MSE value for Gaussian Noise (Image Quality ~ 1/MSE)");
		if(val_msegu[ALPHA]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_msegu[ALPHA]);
            glRasterPos2i(35,(300/log(msegu_max))*log(val_msegu[ALPHA])+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
     if(val_msegu[CONTRA]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_msegu[CONTRA]);
            glRasterPos2i(103,(300/log(msegu_max))*log(val_msegu[CONTRA])+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
      if(val_msegu[MMSE]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_msegu[MMSE]);
            glRasterPos2i(171,(300/log(msegu_max))*log(val_msegu[MMSE])+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
       if(val_msegu[ARITHM]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_msegu[ARITHM]);
            glRasterPos2i(239,(300/log(msegu_max))*log(val_msegu[ARITHM])+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
       if(val_msegu[MTM]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_msegu[MTM]);
            glRasterPos2i(307,(300/log(msegu_max))*log(val_msegu[MTM])+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
	glRasterPos2i(425, 10);
	glutBitmapString(GLUT_BITMAP_HELVETICA_12,"PSNR value for Gaussian Noise (Image Quality ~ PSNR)");
		if(val_psnrgu[ALPHA]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_psnrgu[ALPHA]);
            glRasterPos2i(435,(300/psnrgu_max)*val_psnrgu[ALPHA]+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
     if(val_psnrgu[CONTRA]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_psnrgu[CONTRA]);
            glRasterPos2i(503,(300/psnrgu_max)*val_psnrgu[CONTRA]+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
      if(val_psnrgu[MMSE]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_psnrgu[MMSE]);
            glRasterPos2i(571,(300/psnrgu_max)*val_psnrgu[MMSE]+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
       if(val_psnrgu[ARITHM]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_psnrgu[ARITHM]);
            glRasterPos2i(639,(300/psnrgu_max)*val_psnrgu[ARITHM]+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
       if(val_psnrgu[MTM]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_psnrgu[MTM]);
            glRasterPos2i(707,(300/psnrgu_max)*val_psnrgu[MTM]+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
	glRasterPos2i(825, 10);
	glutBitmapString(GLUT_BITMAP_HELVETICA_12,"SNR value for Gaussian Noise (Image Quality ~ SNR)");
    if(val_snrgu[ALPHA]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_snrgu[ALPHA]);
            glRasterPos2i(835,(300/snrgu_max)*val_snrgu[ALPHA]+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
     if(val_snrgu[CONTRA]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_snrgu[CONTRA]);
            glRasterPos2i(903,(300/snrgu_max)*val_snrgu[CONTRA]+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
      if(val_snrgu[MMSE]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_snrgu[MMSE]);
            glRasterPos2i(971,(300/snrgu_max)*val_snrgu[MMSE]+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
       if(val_snrgu[ARITHM]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_snrgu[ARITHM]);
            glRasterPos2i(1039,(300/snrgu_max)*val_snrgu[ARITHM]+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }
       if(val_snrgu[MTM]!=0)
       {
           char al[10];
           sprintf(al,"%0.3f",val_snrgu[MTM]);
            glRasterPos2i(1107,(300/snrgu_max)*val_snrgu[MTM]+5);
            glutBitmapString(GLUT_BITMAP_HELVETICA_12,al);
       }

    glFinish();
}
void statdisp(void){
char** argv_ = 0;
int argc_ = 0;
glutInit(&argc_,argv_);
glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS);
glutInitWindowSize(1450,700);
glutInitWindowPosition (0, 0);
    glutCreateWindow("Statistical Analysis");
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    init();
    glutDisplayFunc(graphfunct2D);
    glutMainLoop();}

void open_dialog(GtkWidget* button, gpointer window)
{       gtk_widget_hide((GtkWidget*)window);
        GtkWidget *dialog,*image,*loc,*win,*window_1, *vbox, *toggle; // setting a widget for show casing a file selecting space
        dialog = gtk_file_chooser_dialog_new("Choose a file for processing", GTK_WINDOW(window), GTK_FILE_CHOOSER_ACTION_OPEN, GTK_STOCK_OK, GTK_RESPONSE_OK, GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,NULL);
        // for above line we created a special function to accept any file format and show the case using Action_NEW and ok and cancel with creation and response
        gtk_widget_show_all(dialog); // show all the widgets on monitor

        gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(dialog), g_get_home_dir());// it would first display the home directory directly
        gtk_widget_show_all(dialog);// prints all the things under dialog to display on screen
        gint resp = gtk_dialog_run(GTK_DIALOG(dialog));// creating var resp to know the type of response given the user whether yes or no
        if(resp == GTK_RESPONSE_OK){
                x=gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
                gtk_main_quit();

    gtk_widget_destroy(dialog);

    window_1 = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    g_signal_connect(window, "delete-event", G_CALLBACK(gtk_main_quit), NULL);
    gtk_window_set_default_size(GTK_WINDOW(window_1), 250, 100);

	gtk_window_set_title((GtkWindow*)window_1,"Select Option");
	gtk_window_set_position((GtkWindow*)window_1, GTK_WIN_POS_CENTER);
	vbox = gtk_vbox_new(0,0);

	toggle = gtk_toggle_button_new_with_mnemonic("MMSE");
	gtk_box_pack_start(GTK_BOX(vbox), toggle, 0,0,0);
	g_signal_connect(toggle, "toggled", G_CALLBACK(mmsed), NULL);

	toggle = gtk_toggle_button_new_with_mnemonic("Arithmetic");
	gtk_box_pack_start(GTK_BOX(vbox), toggle, 0,0,0);
	g_signal_connect(toggle, "toggled", G_CALLBACK(arithmd), NULL);

	toggle = gtk_toggle_button_new_with_mnemonic("MTM");
	gtk_box_pack_start(GTK_BOX(vbox), toggle, 0,0,0);
	g_signal_connect(toggle, "toggled", G_CALLBACK(mtmd), NULL);

	toggle = gtk_toggle_button_new_with_mnemonic("Contra harmonic");
	gtk_box_pack_start(GTK_BOX(vbox), toggle, 0,0,0);
	g_signal_connect(toggle, "toggled", G_CALLBACK(contrad), NULL);

	toggle = gtk_toggle_button_new_with_mnemonic("Alpha Trim");
	gtk_box_pack_start(GTK_BOX(vbox), toggle, 0,0,0);
	g_signal_connect(toggle, "toggled", G_CALLBACK(alphad), NULL);

	toggle = gtk_toggle_button_new_with_mnemonic("Statistics");
	gtk_box_pack_start(GTK_BOX(vbox), toggle, 0,0,0);
    g_signal_connect_swapped (toggle, "clicked",G_CALLBACK (statdisp),NULL);

	toggle = gtk_toggle_button_new_with_mnemonic("Exit");
	gtk_box_pack_start(GTK_BOX(vbox), toggle, 0,0,0);
    g_signal_connect_swapped (toggle, "clicked",G_CALLBACK (gtk_widget_destroy),window_1);

	gtk_container_add(GTK_CONTAINER(window_1), vbox);
    gtk_widget_show_all(window_1);
    gtk_main(); //start the main loop

        }
         else
                g_print("You pressed Cancel\n");
        gtk_widget_destroy(dialog);
}

int main(int argc, char **argv) {

        gtk_init(&argc, &argv);  // intialize the gtk
        GtkWidget *window, *button;  //declare the needed var
        window = gtk_window_new(GTK_WINDOW_TOPLEVEL); //creating the window
        g_signal_connect(window, "delete-event", G_CALLBACK(gtk_main_quit), NULL); // if 'X' is clicked then exit the prog, command
            //our program starts from here
        gtk_window_set_title((GtkWindow*)window,"Noisy Image Processor");
        button = gtk_button_new_with_label("Select File"); //you will see it
        gtk_container_set_border_width(GTK_CONTAINER(window), 100);//dimensions of the button
        gtk_window_set_default_size(GTK_WINDOW(window), 100, 150);
        gtk_container_add(GTK_CONTAINER(window), button);// add a button icon to click
        g_signal_connect(button, "clicked", G_CALLBACK(open_dialog), window);//after click tell user it is clicked and show the file choose
        gtk_widget_show_all(window);  // show all widgets on to the window
        gtk_main(); //start the main loop

//toggle buttons
   free(infohead);
   free(head);
   free(color);
   fclose(fp_write);
   fclose(fp_read);
   return 0;
}
