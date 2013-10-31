#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct
{
	float	x,y,z;
}float3D;
typedef struct
{
	int	a,b,c;
}int3D;
typedef struct
{
    int a,b;
}int2D;

float norm3D(float3D a)
{
	return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}
float3D sub3D(float3D a, float3D b)
{
	return (float3D){a.x-b.x,a.y-b.y,a.z-b.z};
}
float3D add3D(float3D a, float3D b)
{
	return (float3D){a.x+b.x,a.y+b.y,a.z+b.z};
}
float3D sca3D(float3D a, float t)
{
	return (float3D){a.x*t,a.y*t,a.z*t};
}
#define kTOL	0.00001
#define kSVG	1
int eq3D(float3D a, float3D b)
{
	return (fabs(a.x-b.x)<kTOL && fabs(a.y-b.y)<kTOL && fabs(a.z-b.z)<kTOL);
}
int eq2D(int2D a, int2D b)
{
	return (a.a==b.a && a.b==b.b);
}

int load_mesh(char *path, float **p, int **t, int *np, int *nt)
{
    FILE	*f;
    int		i;
    char	str[1024];
    
    f=fopen(path,"r");
	
	if(f==NULL)
		return 1;
    
	fgets(str,1024,f);
	sscanf(str," %i %i ",np,nt);
    
	// read vertices
	*p=(float*)calloc(*np,sizeof(float)*3);
	for(i=0;i<*np;i++)
	{
		fgets(str,1024,f);
		sscanf(str," %f %f %f ",&((*p)[3*i+0]),&((*p)[3*i+1]),&((*p)[3*i+2]));
	}
	// read triangles
	*t=(int*)calloc(*nt,sizeof(int)*3);
	for(i=0;i<*nt;i++)
	{
		fgets(str,1024,f);
		sscanf(str," %i %i %i ",&((*t)[3*i+0]),&((*t)[3*i+1]),&((*t)[3*i+2]));
    }
	fclose(f);
	return 0;
}
void write_svg_header(FILE *f, float xoff, float yoff, float w, float h, int n, int m, float dw, float dh)
{
    fprintf(f,"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
    fprintf(f,"<svg\n");
    fprintf(f," version=\"1.1\"\n");
    fprintf(f," id=\"Brain contours\"\n");
    fprintf(f," xmlns=\"http://www.w3.org/2000/svg\"\n");
    fprintf(f," width=\"%.3fmm\"\n",w);
    fprintf(f," height=\"%.3fmm\"\n",h);
    fprintf(f," viewBox=\"%.3f %.3f %.3f %.3f\"\n",xoff,yoff,w,h);
    fprintf(f," xml:space=\"preserve\"\n");
    fprintf(f,">\n");

    // grid
    int i;
    for(i=0;i<=m;i++)
        fprintf(f,
                "<path fill=\"none\" stroke=\"#ff0000\" stroke-width=\"0.01mm\" stroke-miterlimit=\"10\" d=\"M%f,%f L%f,%f\"/>\n",
                xoff,yoff+i*dh,xoff+n*dw,yoff+i*dh);
    for(i=0;i<=n;i++)
        fprintf(f,
                "<path fill=\"none\" stroke=\"#ff0000\" stroke-width=\"0.01mm\" stroke-miterlimit=\"10\" d=\"M%f,%f L%f,%f\"/>\n",
                xoff+i*dw,yoff,xoff+i*dw,yoff+m*dh);
    /*
     fprintf(f,"<svg version=\"1.1\" id=\"Brain contours\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" x=\"0mm\" y=\"0mm\" width=\"%.3fmm\" height=\"%.3fmm\" viewBox=\"%.3f %.3f %.3f %.3f\" enable-background=\"new %.3f %.3f %.3f %.3f\" xml:space=\"preserve\">\n",xb-xa,zb-za,xa,za,xb-xa,zb-za,xa,za,xb-xa,zb-za);
     */
}

void write_svg_footer(FILE *f)
{
    fprintf(f,"</svg>\n");
    fclose(f);
}

int main(int argc, char *argv[])
{
	// argv[1]	input text format mesh
	// argv[2]	output svg file
	// argv[3]	input -n number_of_slices or -t thickness_of_slices_in_mm
	// argv[4]	(value for the previous switch)
	// argv[5]	input -w maximum_width_in_mm
	// argv[6]	(value for the previous switch)
	// argv[7]	input -h maximum_height_in_mm
	// argv[8]	(value for the previous switch)
	// argv[9]	input -s scale
	// argv[10]	(value for the previous switch)
	
	float3D	*p,*e,*c;
	int3D	*t;
    int2D   *e2,*c2;
	int		*cindex,ncindex;
	int		np,nt,ne,nc,nv;
	int		err;
	int		i,j,nslices,slicenum,found=0;
	float	y,ya,yb,dy,a,w,h;
	float	xa,xb,za,zb;
	float	xsum,zsum;
	float	scale;                      // scaling factor
	int		m;                          // number of slice rows per page
    int     n;                          // number of slice columns per page
    int     pages;                      // number of pages
    int     mm;                         // current row index
    int     nn;                         // current column index
    int     pp;                         // current page
	int		ntFlag;
    char    fileName[1024];
	FILE	*f;
	
	ntFlag=0;
	if(strcmp(argv[3],"-n")==0)
		ntFlag=1;
	else
        if(strcmp(argv[3],"-t")==0)
            ntFlag=2;
        else
        {
            printf("ERROR: Enter either number (-n) or thickness (-t) of slices\n");
            return 1;
        }
	
	// 1. load mesh
	err=load_mesh(argv[1],(float**)&p,(int**)&t,&np,&nt);
	if(err!=0)
	{
		printf("Error %i reading file\n",err);
		return 1;
	}
	printf("np:%i nt:%i\n",np,nt);
    
	// 2. find frontal- and occipital-most coordinate
	xa=xb=p[0].x;
	ya=yb=p[0].y;
	za=zb=p[0].z;
	for(i=1;i<np;i++)
	{
		if(p[i].x<xa)	xa=p[i].x;
		if(p[i].x>xb)	xb=p[i].x;
        
		if(p[i].y<ya)	ya=p[i].y;
		if(p[i].y>yb)	yb=p[i].y;
        
		if(p[i].z<za)	za=p[i].z;
		if(p[i].z>zb)	zb=p[i].z;
	}
	printf("xmin,xmax=%f,%f (%f)\n",xa,xb,xb-xa);
	printf("ymin,ymax=%f,%f (%f)\n",ya,yb,yb-ya);
	printf("zmin,zmax=%f,%f (%f)\n",za,zb,zb-za);
	
	// determine number and thickness of slices
	if(ntFlag==1) // input: number of slices
	{
		nslices=atoi(argv[4]);
		dy=(yb-ya)/(float)nslices;
	}
	if(ntFlag==2) //input: slice thickness
	{
		dy=atof(argv[4]);
		nslices=(yb-ya)/dy;
	}
    
	// determine number of slice columns
    w=atof(argv[6]);
	n=w/(xb-xa);

    // determine number of slice rows
    h=atof(argv[8]);
    m=h/(zb-za);
    
    // determine number of pages
    pages=nslices/(n*m)+0.5;

	printf("n:%i, m:%i, pages:%i\n",n,m,pages);
    
	// scale
	scale=1;
	if(strcmp(argv[7],"-s")==0)
		scale=atof(argv[8]);
	for(i=0;i<np;i++)
		p[i]=sca3D(p[i],scale);
    
	mm=0;
    nn=0;
    pp=1;
	xsum=0;
	zsum=0;
    
	// 3. cut the contours
	if(kSVG)
	{
        sprintf(fileName,"%s-%03i.svg",argv[2],pp);
        f=fopen(fileName,"w");
        write_svg_header(f, xa,za,w,h,n,m,(xb-xa),(zb-za));
	}

	slicenum=0;
	for(y=ya+dy*0.5;y<yb;y+=dy)
        //y=ya+dy;
	{
		printf("(%i %i) slice %i: ",slicenum/(n+1),nn,slicenum);
		// find number of segments
		ne=0;
		for(i=0;i<nt;i++)
		{
			if((p[t[i].a].y-y)*(p[t[i].b].y-y)<=0)
				ne++;
			if((p[t[i].b].y-y)*(p[t[i].c].y-y)<=0)
				ne++;
			if((p[t[i].c].y-y)*(p[t[i].a].y-y)<=0)
				ne++;
		}
		
		// compute the segments
		e=(float3D*)calloc(ne,sizeof(float3D));
		e2=(int2D*)calloc(ne,sizeof(int2D));
		ne=0;
		for(i=0;i<nt;i++)
		{
			if((p[t[i].a].y-y)*(p[t[i].b].y-y)<=0)
			{
				a=(y-p[t[i].a].y)/(p[t[i].b].y-p[t[i].a].y);
                e2[ne]=(int2D){(t[i].a<t[i].b)?t[i].a:t[i].b,(t[i].a<t[i].b)?t[i].b:t[i].a};
				e[ne++]=add3D(sca3D(p[t[i].a],1-a),sca3D(p[t[i].b],a));
			}
			if((p[t[i].b].y-y)*(p[t[i].c].y-y)<=0)
			{
				a=(y-p[t[i].b].y)/(p[t[i].c].y-p[t[i].b].y);
                e2[ne]=(int2D){(t[i].b<t[i].c)?t[i].b:t[i].c,(t[i].b<t[i].c)?t[i].c:t[i].b};
				e[ne++]=add3D(sca3D(p[t[i].b],1-a),sca3D(p[t[i].c],a));
			}
			if((p[t[i].c].y-y)*(p[t[i].a].y-y)<=0)
			{
				a=(y-p[t[i].c].y)/(p[t[i].a].y-p[t[i].c].y);
                e2[ne]=(int2D){(t[i].c<t[i].a)?t[i].c:t[i].a,(t[i].c<t[i].a)?t[i].a:t[i].c};
				e[ne++]=add3D(sca3D(p[t[i].c],1-a),sca3D(p[t[i].a],a));
			}
		}
		ne=ne/2;
		printf("%i vertices, ",ne);
		
		// join segments into a single contour
		cindex=(int*)calloc(1000,sizeof(int));
		cindex[0]=0;
		ncindex=0;
		
		nv=ne;
		c=(float3D*)calloc(nv+2000,sizeof(float3D));
		c2=(int2D*)calloc(nv+2000,sizeof(int2D));
		nc=0;
		c2[nc]=e2[0];
        c[nc++]=e[0];
        c2[nc]=e2[1];
		c[nc++]=e[1];
        
		
		ne--;
		e[0]=e[2*ne];
		e[1]=e[2*ne+1];
        e2[0]=e2[2*ne];
        e2[1]=e2[2*ne+1];
		
		for(i=1;i<nv;i++)
		{
			found=0;
			for(j=0;j<ne;j++)
			{
                if(eq2D(c2[nc-1],e2[2*j+0]))
				{
					c2[nc]=e2[2*j+1];
                    c[nc++]=e[2*j+1];
					
					ne--;
					e[2*j+0]=e[2*ne];
					e[2*j+1]=e[2*ne+1];
					e2[2*j+0]=e2[2*ne];
					e2[2*j+1]=e2[2*ne+1];
					found=1;
					break;
				}
				else
                    if(eq2D(c2[nc-1],e2[2*j+1]))
                    {
                        c2[nc]=e2[2*j+0];
                        c[nc++]=e[2*j+0];
                        
                        ne--;
                        e[2*j+0]=e[2*ne];
                        e[2*j+1]=e[2*ne+1];
                        e2[2*j+0]=e2[2*ne];
                        e2[2*j+1]=e2[2*ne+1];
                        found=1;
                        break;
                    }
			}
			if(found==0)  // end of closed loop reached
			{
				cindex[++ncindex]=nc;
				
				if(ne==0)
					break;
                
				c2[nc]=e2[0];
				c[nc++]=e[0];
				c2[nc]=e2[1];
				c[nc++]=e[1];
                
				ne--;
                
				e[0]=e[2*ne];
				e[1]=e[2*ne+1];
				e2[0]=e2[2*ne];
				e2[1]=e2[2*ne+1];
			}
			if(ne<0)
				break;
		}
		if(found)
			cindex[++ncindex]=nc;
		printf("%i contours\n",ncindex);
		
		// save the contour
		
		fprintf(f,"<g>\n");
		
		// contours
		for(j=0;j<ncindex;j++)
		{
			if(kSVG)
			{
				fprintf(f,"<path fill=\"none\" stroke=\"#ff0000\" stroke-width=\"0.01mm\" stroke-miterlimit=\"10\" id=\"%i_%i\" d=\"M",slicenum,j);
				for(i=cindex[j];i<cindex[j+1]-1;i++)
					fprintf(f,"%f,%f L",c[i].x+xsum,c[i].z+zsum);
                fprintf(f,"%f,%f ",c[i].x+xsum,c[i].z+zsum);
				fprintf(f,"\"/>\n");
			}
			else
			{
				for(i=cindex[j];i<cindex[j+1];i++)
					printf("%f,%f\n",c[i].x,c[i].z);
				printf("\n");
			}
		}
        
        // slice number
        if(kSVG)
        {
            fprintf(f,
                    "<text x=\"%f\" y=\"%f\" style=\"fill:none;stroke:#000000;stroke-width:0.01mm;stroke-miterlimit:10;font-size:3px\">%i</text>\n",
                    xa+xsum+1,za+zsum+4,slicenum+1);
        }
        
        
		// guides
		for(i=25;i<=75;i+=25)
            for(j=25;j<=75;j+=25)
                fprintf(f,"<path fill=\"none\" stroke=\"#ff0000\" stroke-width=\"0.01mm\" stroke-miterlimit=\"10\" d=\"M%f,%f L%f,%f L%f,%f L%f,%f L%f,%f\"/>\n",
						xsum+xa+i/100.0*(xb-xa)-0.5,zsum+za+j/100.0*(zb-za)-0.5,
						xsum+xa+i/100.0*(xb-xa)+0.5,zsum+za+j/100.0*(zb-za)-0.5,
						xsum+xa+i/100.0*(xb-xa)+0.5,zsum+za+j/100.0*(zb-za)+0.5,
						xsum+xa+i/100.0*(xb-xa)-0.5,zsum+za+j/100.0*(zb-za)+0.5,
						xsum+xa+i/100.0*(xb-xa)-0.5,zsum+za+j/100.0*(zb-za)-0.5);
        
		fprintf(f,"</g>\n");
		
		xsum+=(xb-xa);
		if(++nn==n)
		{
			xsum=0;
			nn=0;
			zsum+=(zb-za);
            if(++mm==m)
            {
                zsum=0;
                mm=0;
                
                // close current file
                write_svg_footer(f);
                
                // open new file
                pp++;
                sprintf(fileName,"%s-%03i.svg",argv[2],pp);
                f=fopen(fileName,"w");
                write_svg_header(f,xa,za,w,h,n,m,(xb-xa),(zb-za));
           }
		}
		slicenum++;
		
		free(e);
		free(e2);
        free(c);
        free(c2);
        
		free(cindex);
	}
    write_svg_footer(f);
	   
	free(p);
	free(t);
	
	return 0;
}
