#include <stdio.h>
#include<iostream>
#include<fstream>
#include<cmath>
#include"ppm.c"
#ifdef _WIN32
  #include <io.h>
  #include <fcntl.h>
#endif

using namespace std;
void bezier(double p1,double p2, double p3, double  p4, double bezier_array[],double derivative_array[]);
double ReadData();
void linearint(double J[],double p1, double  p2, double Lt[],double dt[]);
/* Global Variables*/
double X[100],Lxdt[100],dydx[1000],Lxdtcoord[1000],Lydt[100],Lydtcoord[1000],Lxcoord[1000],Lycoord[1000],Lx[100],Ly[100], Y[100],dx[1000], dy[1000], Xcoord[10000],Ycoord[10000],dxcoord[10000],dycoord[1000],points[100];
int number;
int length;
double Jojox[100], Jojoy[100],Jojo,Jojo2;
float arclength[10000];
int main()
{ 	 int result= _setmode(fileno(stdout),O_BINARY);
	double x[4],y[4],x0,x1,x2,x3,y0,y1,y2,y3;//initial x and y ps and last points on the curve
	int level=1;//helps find the distance along curve
	ReadData();//Reads Spline points from text file

	//initial coordinates	
	x[0]=points[0];
	y[0]=points[1];
	
	x[3]=points[2];
	y[3]=points[3];

	x[1]= x[0]+(x[3]-x[0])/4;
	y[1]=y[0]+(y[3]-y[0])/4;

	x[2]=x[0]+(x[3]-x[0])/2;
	y[2]=y[0]+(y[3]-y[0])/2;

	linearint(Jojox,x[0],x[3],Lx,Lxdt);
	linearint(Jojoy,-y[0],-y[3],Ly,Lydt);
//	Jojo=Jojox[99];
//	fprintf(stderr,"%lf\n",Jojo);
	bezier(x[0], x[1], x[2], x[3],X,dx);//finds x Bezier curve
	bezier(-y[0],-y[1],-y[2],-y[3],Y,dy);//finds y Bezier curve
	
	/*Arclength*/
int numberOfDivisions = 100;
int maxPoint = numberOfDivisions;

// _formula.GetPoint(t) is assumed to take a parameter t and return 
// a point using the Bézier formula.  The _formula object would 
// have the Bézier coefficients to calculate this.


	/*assigning values of first bezier to full bezier*/
	for(int curve_count=0;curve_count<100;curve_count++)
		{
		Xcoord[curve_count]=X[curve_count];
		Ycoord[curve_count]=Y[curve_count];
		dxcoord[curve_count]=dx[curve_count];
		dycoord[curve_count]=dy[curve_count];
		dydx[curve_count]=(-dy[curve_count])/dx[curve_count];
		Lxcoord[curve_count]=Lx[curve_count];
		Lycoord[curve_count]=Ly[curve_count];
		Lxdtcoord[curve_count]=Lxdt[curve_count];
		Lydtcoord[curve_count]=Lydt[curve_count];
		
		}
	//	fprintf(stderr, "%lf\n",dydx[5]);
	//	slope[curve_count]=-1/(Lydt[curve_count]/Lxdt[curve_count];)


/******ARCLENGTH**************/
double previousPointx = Xcoord[0];  // for Bézier == P0 control point
double previousPointy = Ycoord[0]; 

arclength[0]=0;
for (int i = 1; i < maxPoint; i++)
{
    double arcx = Xcoord[i];
		double arcy = Ycoord[i];
    arclength[i] += arclength[i-1]+sqrt(pow(arcx-previousPointx,2)+pow(arcy-previousPointy,2));
//	fprintf(stderr,"%d %lf\n",i,arclength[i]);
    previousPointx = arcx;
		previousPointy = arcy;
}
//fprintf(stderr,"%d %lf\n",99,arclength[99]);
		
	
	for(int Loop=0; Loop<number-2;Loop++)
		{
		/****Continuity****/
		/**X**/
		x1=x[1]; x2=x[2]; x3=x[3];//buffer
		x[0]=x3; /*C0 continuity pn=q0
			  	 A(1)=B(0)*/
		x[1]=2*x3-x2; /*C1 continuity Q1=2*Pn-Pn-1
			  		  A'(1)=B'(0)*/
		x[2]=4*x3-4*x2+x1;/*C2 cont Q2=2*Q1-2*Pn-1-Pn-2 */
		
		/**Y**/
		y0=y[0];y1=y[1];y2=y[2];y3=y[3];//buffer
		y[0]=y3;
		y[1]=2*y3-y2;
		y[2]=4*y3-4*y2+y1;
		

		x[3]=points[4+2*Loop];//input end of next curve
		y[3]=points[5+2*Loop];
		//printf("Bezier Curve #%d\n",c_number);
		bezier(x[0],x[1],x[2],x[3],X,dx);
		bezier(-y[0],-y[1],-y[2],-y[3],Y,dy);
		linearint(Jojox,x[0],x[3],Lx,Lxdt);
		linearint(Jojoy,-y[0],-y[3],Ly,Lydt);
	//	Jojo=Jojo+Jojox[99];

		/****assigning next value for bezeir****/
		for(int curve_count=0;curve_count<100;curve_count++)
			{
			Xcoord[100*level+curve_count]=X[curve_count];
			Ycoord[100*level+curve_count]=Y[curve_count];
			dxcoord[100*level+curve_count]=dx[curve_count];
			dycoord[100*level+curve_count]=dy[curve_count];
			dydx[100*level+curve_count]=-dy[curve_count]/dx[curve_count];
			Lxcoord[100*level+curve_count]=Lx[curve_count];
			Lycoord[100*level+curve_count]=Ly[curve_count];
			Lxdtcoord[100*level+curve_count]=Lxdt[curve_count];
			Lydtcoord[100*level+curve_count]=Lydt[curve_count];
	//		slope[100*level+curve_count]=-1/(Lydt[100*level+curve_count]/Lxdt[100*level+curve_count]);
			
			}	
	arclength[level*100]=arclength[(level*100)-1];
for (int i = 1; i < maxPoint; i++)
{
    double arcx = Xcoord[level*100+i];
		double arcy = Ycoord[level*100+i];
    arclength[level*100+i] += arclength[(level*100+i)-1]+sqrt(pow(arcx-previousPointx,2)+pow(arcy-previousPointy,2));
//	fprintf(stderr,"%d %lf\n",i,arclength[i]);
    previousPointx = arcx;
		previousPointy = arcy;
}
		level++;
		
		Xcoord[((100*level)-1)]=x[3];
		Ycoord[((100*level)-1)]=- y[3];
		Lxcoord[((100*level)-1)]=x[3];
		Lycoord[((100*level)-1)]=-y[3];
	}
//fprintf(stderr, "%lf\n", arclength[(100*level-1)]);
float max=arclength[(100*level)-1];
float u; // the parameterized value, assumed to be between 0 and 1 represent percentage along curve so .25=25% along the curve 
float finalt[100]; // we are attempting to find t for u
 int index[100];
// a pre-calculated list of arc-lengths along the Bézier 
for(u=0;u<100;u++)
{
// get the target arcLength of curve for parameter u
float targetArcLength; 
targetArcLength= (u/100)*max;

// the next function would be a binary search, for efficiency

for (int i=0; i<level*100;i++)
{
if(arclength[i]>targetArcLength)
{
	index[(int)u]=i-1;
	i=level*100;
//	fprintf(stderr,"1\n");
}
else if(arclength[i]==targetArcLength)
{
	index[(int)u]=i;
	i=level*100;
	//fprintf(stderr,"2\n");
}
else
{
index[(int)u]=i;
//fprintf(stderr,"3\n");
}
}
//fprintf(stderr,"%lf %d %lf\n", u, index[(int)u],targetArcLength);

//fprintf(stderr, "%d %lf\n", index[(int)u], targetArcLength);
//}
//fprintf(stderr,"%lf %lf\n",index,targetArcLength);.


// if exact match, return t based on exact index
/*if (index == targetArcLength)
{
    t = index / (float) (arcLengths.Length - 1);
}
else  // need to interpolate between two points
{
    float lengthBefore = _arcLengths[index];
    float lengthAfter = _arcLengths[index+1];
    float segmentLength = lengthAfter - lengthBefore;

    // determine where we are between the 'before' and 'after' points.
    float segmentFraction = (targetLength - lengthBefore) / segmentLength;
                          
    // add that fractional amount to t 
    t = (index + segmentFraction) / (float) (_arcLengths.Length -1);
}*/
}


// function to calculate straight-line distance between two points:
	//	fprintf(stderr,"%d %lf\n",i,arclength[i]);

	/*Image Part*/
	image original, distorted;
	pixel P,Px,Py,P_p,P_n;//positive pixel and negative pixel
	double Xco[1000],Yco[1000],Uco[1000],Vco[1000];
	int X_p,X_n,Y_p,Y_n;

	
	//fprintf(stderr,"%lf",Jojo2);
	//copy original image to distorted imagidwe
	//for t=0 to 100  (from 0% along the spline to 100% along the spline, head to tail)
	if(get_ppm( &original )&& copy_ppm(original,&distorted))

	{
	
	
	for(int t=0;t<100;t++) // 0 to 500 
										// 200
		{//array of x and why coordinates/ this is not how finding the percentage works, so its a placeholder
		//Xco[0]=Xcoord[0]; 
		//Yco[0]=Xcoord[0];
		//set <x,y> = spline at position t Question: point at pos t or spline formula at pos t or does that mean the same thing
		int one = index[t];

	//	Xco[t]=Xcoord[one];
	//	Yco[t]=-Ycoord[one];
		Xco[t]=Xcoord[one];
		Yco[t]=-Ycoord[one];
	//	Xco[99]=Xcoord[((100*level)-1)];
		// Yco[99]=-Ycoord[((100*level)-1)];
		//set <u,v> = vector perpendicular to spline/u,v= (-dy, dx) dy=y2-y1
//		Uco[t]=dxcoord[one];;
//		Vco[t]=-dycoord[one];; 
	if(dycoord[one]==0)
	{
		Uco[t]=0;
		Vco[t]=0; 
		
	}
	else
	{
		Uco[t]=-1/(dydx[one]);
		Vco[t]=-1/(dydx[one]); 
	}
	//	fprintf(stderr, "#%d X:%lf Y:%lf\n",t,Xco[t],Yco[t]);
	//	fprintf(stderr,"#%d dx:%lf dy:%lf\n",t, Uco[t],Vco[t]);
		//set pixel( t, 50 ) = original image pixel( x, y )
//		P=get_pixel(original,Xco[t],Yco[t]);
//		put_pixel( distorted, t, 50, P);

		for(int v=0; v<50;v++)
			{//for i=0..50:  
        	//set pixel( t, 50+i )= original image pixel(  x+iu, y+iv )  set pixel(Xco[t])
        	//set pixel( t, 50-i )=  original image pixel(  x-iu, y-iv ) 

				X_p=Xco[t]+v*Uco[t];//x+iu
				X_n=Xco[t]-v*Uco[t];//x-iu
				Y_p=v*Vco[t]+Yco[t];
				Y_n=-v*Vco[t]+Yco[t];
			if((X_p>0)&&(X_p<original.width)&&(Y_n>0)&&Y_n<original.height)
			{
			P_p=get_pixel(original, X_p,Y_n);
				put_pixel(distorted, t,50-v,P_p);
				fprintf(stderr, "%d X_p=%d X_n=%d Y_p=%d Y_n=%d\n",v,X_p,X_n,Y_p,Y_n);
			} 
			else
			{
				P_p= get_pixel(original, 0, 0);
				put_pixel(distorted, 0, 0,P_p);
			}
			if(((X_n>0)&&(X_n<original.width))&&((Y_p>0)&&((Y_p<original.height))))
			{
			P_n=get_pixel(original, X_n,Y_p);
			put_pixel(distorted, t,50+v,P_n);
		//		fprintf(stderr, "%d X_p=%d X_n=%d Y_p=%d Y_n=%d\n",v,X_p,X_n,Y_p,Y_n);
			} 
						else
			{
				P_p= get_pixel(original, 0,0);
				put_pixel(distorted, 0, 0,P_p);
			}
				fprintf(stderr, "%d X_p=%d X_n=%d Y_p=%d Y_n=%d\n",v,X_p,X_n,Y_p,Y_n);
			//	P_p=get_pixel(original, X_p,Y_n);
			//P_n=get_pixel(original, X_n,Y_p);
			//	put_pixel(distorted, t,50+v,P_p);
			//	put_pixel(distorted, t,50-v,P_n);
		}
		
		}
		
	
			put_ppm( distorted);
		dealloc_ppm( original);
		dealloc_ppm( distorted );
	}
	return 0;
}

void bezier(double p1,double p2, double p3, double  p4, double bezier_array[],double derivative_array[])
{
	double t=0;
	for (t = 0.0; t < 100; t ++)
	{//Bezier Caluculation
	//B(t) = (1 - t)3P0 + 3(1-t)2tP1 + 3(1-t)t2P2 + t3P3
		double Ct = pow(1 - (t/100), 3)*p1 + 3 * (t / 100)*pow(1 - (t / 100), 2)*p2 + 3 * pow((t / 100), 2)*(1 - (t / 100))*p3 + pow((t / 100), 3)*p4;
		double Dt = 3*pow(1-(t/100),2)*(p2-p1)+6*(1-(t/100))*t*(p3-p2)+3*pow((t/100),2)*(p4-p3);
		Dt=Dt;
		bezier_array[(int)t] = Ct;//used for outputting point
		derivative_array[(int)t]=Dt;		
	}
}
void linearint(double J[], double p1, double  p2, double Lt[],double dt[])
 { 
	 for(double t=0.0; t<100; t++)
	 {
	 Lt[(int)t]=(1-(t/100))*p1+t/100*p2;
	 dt[(int)t]=(p2-p1);

	 }
	
 }
double ReadData()
{
	std::fstream myfile("Values.txt", std::ios_base::in);
	
	
    myfile >> number;
	for(int input=0; input<number*2; input++)
	{
	myfile>>points[input];
	}
	return points[1000];
}
