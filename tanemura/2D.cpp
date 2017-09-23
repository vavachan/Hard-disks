#include<iostream>
#include<math.h>
#include<fstream>
using namespace std;
double box=0;
double twob;
double density;
struct atom;
struct delunay
{
	int A;
	int B;
	float circum_x=0.;
	float circum_y=0.;
	delunay *next=NULL;
};
class set_of_delunay
{
	public :
		delunay *initial=NULL;
};
struct atom 
{
	double x=0;
	double y=0;
	int neighlist[500];
	int contigous[500];
	int edge_index[500]={0};
	int neighbours=0;
	int conti=0;
	double radius=1.;
	set_of_delunay D;
 	//delunay D[12];
};

////////struct delunay
////////{
////////	atom A,B;
////////};
void update_neighbours(atom Atoms[],int nAtoms)
{
	double R_CUT;
	R_CUT=sqrt(200./(4*3.14*density));
	//cout<<R_CUT<<"\n";
	for(int i=0;i<nAtoms-1;i++)
	{
		for(int j=i+1;j<nAtoms;j++)
		{
			double drx,dry,dr;
			drx=Atoms[i].x-Atoms[j].x;
			dry=Atoms[i].y-Atoms[j].y;
			drx=(drx-(twob*lround(drx/twob)));
			dry=(dry-(twob*lround(dry/twob)));
			dr=drx*drx+dry*dry;
			if(dr<R_CUT*R_CUT)
			{
				Atoms[i].neighlist[Atoms[i].neighbours]=j;
				Atoms[j].neighlist[Atoms[j].neighbours]=i;
				Atoms[i].neighbours++;	
				Atoms[j].neighbours++;	
			}
		}
		//cout<<"nei="<<Atoms[i].neighbours<<"\n";
	}
}
void first_delunay(atom *ATOM,atom Atoms[])
{
	double drx,dry,dr;
	double min=2.*twob*twob;
	int nearest,flag;
	for(int i=0;i<ATOM->neighbours;i++)
	{
	////////drx=ATOM->x-Atoms[ATOM->neighlist[i]].x;
	////////dry=ATOM->y-Atoms[ATOM->neighlist[i]].y;
	////////drx=(drx-(twob*lround(drx/twob)));
	////////dry=(dry-(twob*lround(dry/twob)));
	////////dr=drx*drx+dry*dry;
	//////////cout<<dr<<"\t"<<i<<"\n";
	////////if(dr<min)
	////////{
	////////	min=dr;
	////////	nearest=ATOM->neighlist[i];
	////////}
		double X,Y,x,y;
		double rA,rS,DIS,l,dis_i,tan_sq;
		X=Atoms[ATOM->neighlist[i]].x-ATOM->x;
	        Y=Atoms[ATOM->neighlist[i]].y-ATOM->y;
	        X=(X-(twob*lround(X/twob)));
	        Y=(Y-(twob*lround(Y/twob)));
		DIS=sqrt(X*X+Y*Y);
		rA=Atoms[ATOM->neighlist[i]].radius;
		rS=ATOM->radius;
		l=0.5*(DIS+(rS*rS-rA*rA)/DIS);
		x=l/DIS*X;
		y=l/DIS*Y;
		dis_i=(x*x+y*y);
		tan_sq=dis_i-rS*rS;
		if(tan_sq<min)
		{
			min=tan_sq;
	        	nearest=ATOM->neighlist[i];
		}

	}
	ATOM->contigous[ATOM->conti]=nearest;
	ATOM->edge_index[ATOM->conti]++;
	ATOM->conti++;
////////cout<<"nea\t"<<nearest<<"\t"<<min<<"\n";
////////cout<<Atoms[nearest].x<<"\t"<<Atoms[nearest].y<<"\n";
	double DIS_MIN=box*box;
	int DIS_atom;
	double dis;
	double X,Y;
	double circx,circy;
	for(int i=0;i<ATOM->neighbours;i++)
	{
		if(ATOM->neighlist[i]!=nearest)
		{
			atom L=*ATOM;
			atom M=Atoms[ATOM->contigous[0]];
			atom R=Atoms[ATOM->neighlist[i]];
			double ax=M.x-L.x;
			double ay=M.y-L.y;
			double bx=R.x-L.x;
			double by=R.y-L.y;
			ax=(ax-(twob*lround(ax/twob)));
			ay=(ay-(twob*lround(ay/twob)));
			bx=(bx-(twob*lround(bx/twob)));
			by=(by-(twob*lround(by/twob)));
	////////	double A=ax*ax+ay*ay;				
	////////	double B=bx*bx+by*by;
	////////	double x=0.5*(ay*B-by*A)/(ay*bx-ax*by);
	////////	double denom=2.*(ay*bx-ax*by);
	////////	double y=-1.*ax/ay*x+A/(2.*ay);
	////////	dis=sqrt(pow(x-ax,2)+pow(y-ay,2));
			double rA,rS,rB;
			double XA,YA,XB,YB;
			double xA,yA,xB,yB;
			double l;
			double DISA;
			double DISB;
			double MA,MB,INMA,INMB;
			double CA,CB;
			double tan_sq;
			XA=ax;//Atoms[Atoms[SAM].contigous[D->A]].x-Atoms[SAM].x;
			YA=ay;//Atoms[Atoms[SAM].contigous[D->A]].y-Atoms[SAM].y;
			XB=bx;//Atoms[Atoms[SAM].contigous[D->B]].x-Atoms[SAM].x;
			YB=by;//Atoms[Atoms[SAM].contigous[D->B]].y-Atoms[SAM].y;
		////////XA=(XA-(twob*lround(XA/twob)));
		////////YA=(YA-(twob*lround(YA/twob)));
		////////XB=(XB-(twob*lround(XB/twob)));
		////////YB=(YB-(twob*lround(YB/twob)));
		////////cout<<"xa=";
		  //      cout<<XA<<"\t"<<YA<<"\n";
		////////cout<<"xb=";
		    //    cout<<XB<<"\t"<<YB<<"\n";
			DISA=sqrt(XA*XA+YA*YA);
			DISB=sqrt(XB*XB+YB*YB);
		////////cout<<"dis="<<DISA<<"\n";
		////////cout<<"dis="<<DISB<<"\n";
			MA=YA/XA;
			MB=YB/XB;
			INMA=-1./MA;
			INMB=-1./MB;
			rA=M.radius;
			rB=R.radius;
		////////cout<<"radiusA="<<rA<<"\n";
		////////cout<<"radiusb="<<rB<<"\n";
			rS=L.radius;
			//cout<<"radiusS="<<rS<<"\n";
			l=0.5*(DISA+(rS*rS-rA*rA)/DISA);
			xA=l/DISA*XA;
			yA=l/DISA*YA;
			l=0.5*(DISB+(rS*rS-rB*rB)/DISB);
			xB=l/DISB*XB;
			yB=l/DISB*YB;
			CA=yA-INMA*xA;
			CB=yB-INMB*xB;
			X=(CB-CA)/(INMA-INMB);
			Y=INMA*X+CA;
		//	cout<<"thisguy2=";
		//	cout<<X<<"\t"<<Y<<"\n";
			tan_sq=X*X+Y*Y-rS*rS;
			X=X+L.x;
			Y=Y+L.y;
////////  		cout<<Atoms[ATOM->neighlist[i]].x<<"\t"<<Atoms[ATOM->neighlist[i]].y<<"\n";
////////  		cout<<"DIS_MIN="<<DIS_MIN<<"tan="<<tan_sq<<"\n";
////////		cout<<X<<"\t"<<Y<<"\n";
			if(DIS_MIN > tan_sq)
			{
                		DIS_MIN=tan_sq;
                		DIS_atom=ATOM->neighlist[i];
	        		circx=X;
	        		circy=Y;
			}
		}

	}
////////cout<<Atoms[ATOM->neighlist[101]].x<<"\t"<<Atoms[ATOM->neighlist[101]].y<<"\n";
////////cout<<ATOM->neighlist[101]<<"\n";
	ATOM->contigous[ATOM->conti]=DIS_atom;
////////cout<<Atoms[ATOM->neighlist[101]].x<<"\t"<<Atoms[ATOM->neighlist[101]].y<<"\n";
////////cout<<ATOM->neighlist[101]<<"\n";
	ATOM->edge_index[ATOM->conti]++;
	ATOM->conti++;
	ATOM->D.initial= new delunay;
	ATOM->D.initial->A=0;
	ATOM->D.initial->B=1;
	ATOM->D.initial->circum_x=circx;
	ATOM->D.initial->circum_y=circy;
	//cout<<ATOM->D.initial->A<<"\n";
	//cout<<ATOM->D.initial<<"\n";
}
void print_delunay(atom *ATOM,delunay *D,atom Atoms[])
{
	float Sx,Sy;
	float Px,Py;
	Sx=ATOM->x-Atoms[ATOM->contigous[D->A]].x;
	Sy=ATOM->y-Atoms[ATOM->contigous[D->A]].y;
	Sx=(Sx-(twob*lround(Sx/twob)));
	Sy=(Sy-(twob*lround(Sy/twob)));
	cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\n";
	Sx=ATOM->x-Atoms[ATOM->contigous[D->B]].x;
	Sy=ATOM->y-Atoms[ATOM->contigous[D->B]].y;
	Sx=(Sx-(twob*lround(Sx/twob)));
	Sy=(Sy-(twob*lround(Sy/twob)));
	cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\n";
	Sx=ATOM->x-Atoms[ATOM->contigous[D->B]].x;
	Sy=ATOM->y-Atoms[ATOM->contigous[D->B]].y;
	Sx=(Sx-(twob*lround(Sx/twob)));
	Sy=(Sy-(twob*lround(Sy/twob)));
	Px=ATOM->x-Atoms[ATOM->contigous[D->A]].x;
	Py=ATOM->y-Atoms[ATOM->contigous[D->A]].y;
	Px=(Px-(twob*lround(Px/twob)));
	Py=(Py-(twob*lround(Py/twob)));
	cout<<ATOM->x-Px<<"\t"<<ATOM->y-Py<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\n";
////////Sx=Atoms[ATOM->contigous[D->A]].x-Atoms[ATOM->contigous[D->B]].x;
////////Sy=Atoms[ATOM->contigous[D->A]].y-Atoms[ATOM->contigous[D->B]].y;
////////Sx=(Sx-(twob*lround(Sx/twob)));
////////Sy=(Sy-(twob*lround(Sy/twob)));
}
void print_edge(atom *ATOM,delunay *D,atom Atoms[])
{


}

int calculate_line(double a,double b,double c,double d,double x,double y)
{
	double m=(d-b)/(c-a);
	double C=b-m*a;
	if((y-(m*x+C))<0.)
		return -1;
	else
		return 1;
}
void complete_del(atom *ATOM,atom Atoms[],int nAtoms)
{
	double Y_MIN=box*box;
	int DIS_atom;
	for(int i=0;i<ATOM->conti;i++)
	{
		Y_MIN=box*box;
		delunay *D;
		D=ATOM->D.initial;
		if(ATOM->edge_index[i]!=2)
		{
			while(1)
			{
				
				if(D->A==i)
				{
					break;
				}
				if(D->B==i)
				{
					break;
				}
				if(D->next)
					D=D->next;
				else 
				{
					cout<<"impossible\n";
					return;
				}
				
			}
			if(i==D->A)
			{
				//cout<<"hereA\n";
				double a=ATOM->x;
				double b=ATOM->y;
				double c=Atoms[ATOM->contigous[D->A]].x;	
				double d=Atoms[ATOM->contigous[D->A]].y;	
				double Sx=Atoms[ATOM->contigous[D->B]].x-a;
				double Sy=Atoms[ATOM->contigous[D->B]].y-b;
				Sx=(Sx-(twob*lround(Sx/twob)));
				Sy=(Sy-(twob*lround(Sy/twob)));
				double Y=((d-b)-(twob*lround((d-b)/twob)));
				double X=((c-a)-(twob*lround((c-a)/twob)));
				double m=Y/X;
				double C=0.;
				int sign;
				int flag=1;
				int k;
				double circx;
				double circy;
				if((Sy-(m*Sx+C))<0.)
					sign=-1;
				else
					sign=1;
				for(int j=0;j<ATOM->neighbours;j++)
				{
					//cout<<"this is a neighbour="<<ATOM->neighlist[j]<<"\n";
					
					int sign_N;
					double x=Atoms[ATOM->neighlist[j]].x-a;
					double y=Atoms[ATOM->neighlist[j]].y-b;
					x=(x-(twob*lround(x/twob)));
					y=(y-(twob*lround(y/twob)));
					if((y-(m*x+C))<0.)
						sign_N=-1;
					else
						sign_N=1;
				//		cout<<"sign_N="<<sign_N<<"\n";
					flag=1;
				        for(k=0;k<ATOM->conti;k++)
				        {
				        	if(ATOM->contigous[k]==ATOM->neighlist[j])
				        	{
				        		if(ATOM->edge_index[k]==2)
				        		{
				        			flag=0;
				        			  //cout<<"when\t"<<j<<"\n";
				        			  //cout<<Atoms[ATOM->neighlist[j]].x<<"\t"<<Atoms[ATOM->neighlist[j]].y<<"\n";
				        			break;
				        		}
				        	}
				        }
					//cout<<j<<"\t"<<flag<<"\n";
					if(sign!=sign_N && flag)
					{
						atom L=*ATOM;
						atom M=Atoms[ATOM->contigous[D->A]];
						atom R=Atoms[ATOM->neighlist[j]];
						double ax=M.x-L.x;
						double ay=M.y-L.y;
						double bx=R.x-L.x;
						double by=R.y-L.y;
						ax=(ax-(twob*lround(ax/twob)));
						ay=(ay-(twob*lround(ay/twob)));
						bx=(bx-(twob*lround(bx/twob)));
						by=(by-(twob*lround(by/twob)));
						double rA,rS,rB;
						double XA,YA,XB,YB;
						double xA,yA,xB,yB;
						double l;
						double DISA;
						double DISB;
						double MA,MB,INMA,INMB;
						double CA,CB;
						double tan_sq;
						XA=ax;//Atoms[Atoms[SAM].contigous[D->A]].x-Atoms[SAM].x;
						YA=ay;//Atoms[Atoms[SAM].contigous[D->A]].y-Atoms[SAM].y;
						XB=bx;//Atoms[Atoms[SAM].contigous[D->B]].x-Atoms[SAM].x;
						YB=by;//Atoms[Atoms[SAM].contigous[D->B]].y-Atoms[SAM].y;
					////////XA=(XA-(twob*lround(XA/twob)));
					////////YA=(YA-(twob*lround(YA/twob)));
					////////XB=(XB-(twob*lround(XB/twob)));
					////////YB=(YB-(twob*lround(YB/twob)));
					////////cout<<"xa=";
					  //      cout<<XA<<"\t"<<YA<<"\n";
					////////cout<<"xb=";
					    //    cout<<XB<<"\t"<<YB<<"\n";
						DISA=sqrt(XA*XA+YA*YA);
						DISB=sqrt(XB*XB+YB*YB);
					////////cout<<"dis="<<DISA<<"\n";
					////////cout<<"dis="<<DISB<<"\n";
						MA=YA/XA;
						MB=YB/XB;
						INMA=-1./MA;
						INMB=-1./MB;
						rA=M.radius;
						rB=R.radius;
					////////cout<<"radiusA="<<rA<<"\n";
					////////cout<<"radiusb="<<rB<<"\n";
						rS=L.radius;
						//cout<<"radiusS="<<rS<<"\n";
						l=0.5*(DISA+(rS*rS-rA*rA)/DISA);
						xA=l/DISA*XA;
						yA=l/DISA*YA;
						l=0.5*(DISB+(rS*rS-rB*rB)/DISB);
						xB=l/DISB*XB;
						yB=l/DISB*YB;
						CA=yA-INMA*xA;
						CB=yB-INMB*xB;
						X=(CB-CA)/(INMA-INMB);
						Y=INMA*X+CA;
					//	cout<<"thisguy2=";
					//	cout<<X<<"\t"<<Y<<"\n";
						tan_sq=X*X+Y*Y-rS*rS;
				////////	double A=ax*ax+ay*ay;				
				////////	double B=bx*bx+by*by;
				////////	double x=0.5*(ay*B-by*A)/(ay*bx-ax*by);
				////////	double denom=2.*(ay*bx-ax*by);
				////////	double y=-1.*ax/ay*x+A/(2.*ay);
				////////	double dis=sqrt(pow(x-ax,2)+pow(y-ay,2));
				////////	double Y=sqrt(pow(x-ax/2.,2)+pow(y-ay/2.,2));
						double Y_AXIS=sqrt(pow(X-xA,2)+pow(Y-yA,2));
						int sign_C;
						if((Y-(m*X+C))<0.)
							sign_C=-1;
						else
							sign_C=1;
						if(sign_C!=sign_N)
							Y_AXIS=-1.*Y_AXIS;
					//			cout<<"circum="<<dis<<"\n";
						X=X+L.x;
						Y=Y+L.y;
						if(Y_AXIS<Y_MIN)
						{
							;
							Y_MIN=Y_AXIS;
						        DIS_atom=ATOM->neighlist[j];
							circx=X;
							circy=Y;

						}
					}
					
				}//j loop 
				flag=1;
				for(k=0;k<ATOM->conti;k++)
				{
					if(ATOM->contigous[k]==DIS_atom)
					{
						flag=0;
						break;
					}
				}

				if(flag)
				{
					ATOM->contigous[ATOM->conti]=DIS_atom;
					ATOM->edge_index[D->A]++;
					ATOM->edge_index[ATOM->conti]++;
				}
				else 
				{
					ATOM->edge_index[D->A]++;
					ATOM->edge_index[k]++;
				}

				//cout<<DIS_atom<<"\n";
				//cout<<Atoms[DIS_atom].x<<"\t"<<Atoms[DIS_atom].y<<"\n";
				delunay *temp;
				temp=D;
				while(1)
				{
					if(temp->next)
						temp=temp->next;
					else 
						break;

				}
				temp->next= new delunay;
				if(flag)
					temp->next->A=ATOM->conti;
				else 
					temp->next->A=k;
				temp->next->B=D->A;
				temp->next->circum_x=circx;
				temp->next->circum_y=circy;
				if(flag)
					ATOM->conti++;		
			}// D->A loop
			else
			{
				//cout<<"hereb\n";
				double a=ATOM->x;
				double b=ATOM->y;
				double c=Atoms[ATOM->contigous[D->B]].x;	
				double d=Atoms[ATOM->contigous[D->B]].y;	
				double Sx=Atoms[ATOM->contigous[D->A]].x-a;
				double Sy=Atoms[ATOM->contigous[D->A]].y-b;
				Sx=(Sx-(twob*lround(Sx/twob)));
				Sy=(Sy-(twob*lround(Sy/twob)));
				double Y=((d-b)-(twob*lround((d-b)/twob)));
				double X=((c-a)-(twob*lround((c-a)/twob)));
				double m=Y/X;
				double C=0.;
				int sign;
				int flag=1;
				int k;
				double circx;
				double circy;
				if((Sy-(m*Sx+C))<0.)
					sign=-1;
				else
					sign=1;
				for(int j=0;j<ATOM->neighbours;j++)
				{
				//	cout<<"this is a neighbour="<<ATOM->neighlist[i]<<"\n";
					int sign_N;
					double x=Atoms[ATOM->neighlist[j]].x-a;
					double y=Atoms[ATOM->neighlist[j]].y-b;
					x=(x-(twob*lround(x/twob)));
					y=(y-(twob*lround(y/twob)));
					if((y-(m*x+C))<0.)
						sign_N=-1;
					else
						sign_N=1;
				//	cout<<"sign_N="<<sign_N<<"\n";
				
					flag=1;
				        for(k=0;k<ATOM->conti;k++)
				        {
				        	if(ATOM->contigous[k]==ATOM->neighlist[j])
				        	{
				        		if(ATOM->edge_index[k]==2)
				        		{
				        			flag=0;
				        			break;
				        		}
				        	}
				        }
					//cout<<j<<"\t"<<flag<<"\n";
					if(sign!=sign_N && flag)
					{
						atom L=*ATOM;
						atom M=Atoms[ATOM->contigous[D->B]];
						atom R=Atoms[ATOM->neighlist[j]];
						double ax=M.x-L.x;
						double ay=M.y-L.y;
						double bx=R.x-L.x;
						double by=R.y-L.y;
						ax=(ax-(twob*lround(ax/twob)));
						ay=(ay-(twob*lround(ay/twob)));
						bx=(bx-(twob*lround(bx/twob)));
						by=(by-(twob*lround(by/twob)));
						double rA,rS,rB;
						double XA,YA,XB,YB;
						double xA,yA,xB,yB;
						double l;
						double DISA;
						double DISB;
						double MA,MB,INMA,INMB;
						double CA,CB;
						double tan_sq;
						XA=ax;//Atoms[Atoms[SAM].contigous[D->A]].x-Atoms[SAM].x;
						YA=ay;//Atoms[Atoms[SAM].contigous[D->A]].y-Atoms[SAM].y;
						XB=bx;//Atoms[Atoms[SAM].contigous[D->B]].x-Atoms[SAM].x;
						YB=by;//Atoms[Atoms[SAM].contigous[D->B]].y-Atoms[SAM].y;
					////////XA=(XA-(twob*lround(XA/twob)));
					////////YA=(YA-(twob*lround(YA/twob)));
					////////XB=(XB-(twob*lround(XB/twob)));
					////////YB=(YB-(twob*lround(YB/twob)));
					////////cout<<"xa=";
					  //      cout<<XA<<"\t"<<YA<<"\n";
					////////cout<<"xb=";
					    //    cout<<XB<<"\t"<<YB<<"\n";
						DISA=sqrt(XA*XA+YA*YA);
						DISB=sqrt(XB*XB+YB*YB);
					////////cout<<"dis="<<DISA<<"\n";
					////////cout<<"dis="<<DISB<<"\n";
						MA=YA/XA;
						MB=YB/XB;
						INMA=-1./MA;
						INMB=-1./MB;
						rA=M.radius;
						rB=R.radius;
					////////cout<<"radiusA="<<rA<<"\n";
					////////cout<<"radiusb="<<rB<<"\n";
						rS=L.radius;
						//cout<<"radiusS="<<rS<<"\n";
						l=0.5*(DISA+(rS*rS-rA*rA)/DISA);
						xA=l/DISA*XA;
						yA=l/DISA*YA;
						l=0.5*(DISB+(rS*rS-rB*rB)/DISB);
						xB=l/DISB*XB;
						yB=l/DISB*YB;
						CA=yA-INMA*xA;
						CB=yB-INMB*xB;
						X=(CB-CA)/(INMA-INMB);
						Y=INMA*X+CA;
					//	cout<<"thisguy2=";
					//	cout<<X<<"\t"<<Y<<"\n";
						tan_sq=X*X+Y*Y-rS*rS;
				////////	double A=ax*ax+ay*ay;				
				////////	double B=bx*bx+by*by;
				////////	double x=0.5*(ay*B-by*A)/(ay*bx-ax*by);
				////////	double denom=2.*(ay*bx-ax*by);
				////////	double y=-1.*ax/ay*x+A/(2.*ay);
				////////	double dis=sqrt(pow(x-ax,2)+pow(y-ay,2));
				////////	double Y=sqrt(pow(x-ax/2.,2)+pow(y-ay/2.,2));
						double Y_AXIS=sqrt(pow(X-xA,2)+pow(Y-yA,2));
						int sign_C;
						if((Y-(m*X+C))<0.)
							sign_C=-1;
						else
							sign_C=1;
						if(sign_C!=sign_N)
							Y_AXIS=-1.*Y_AXIS;
					//			cout<<"circum="<<dis<<"\n";
						X=X+L.x;
						Y=Y+L.y;
						if(Y_AXIS<Y_MIN)
						{
							;
							Y_MIN=Y_AXIS;
						        DIS_atom=ATOM->neighlist[j];
							circx=X;
							circy=Y;

						}
					////////double A=ax*ax+ay*ay;				
					////////double B=bx*bx+by*by;
					////////double x=0.5*(ay*B-by*A)/(ay*bx-ax*by);
					////////double denom=2.*(ay*bx-ax*by);
					////////double y=-1.*ax/ay*x+A/(2.*ay);
					////////double dis=sqrt(pow(x-ax,2)+pow(y-ay,2));
					////////double Y=sqrt(pow(x-ax/2.,2)+pow(y-ay/2.,2));
					////////int sign_C;
					////////if((y-(m*x+C))<0.)
					////////	sign_C=-1;
					////////else
					////////	sign_C=1;
					////////if(sign_C!=sign_N)
					////////	Y=-1.*Y;

				        //////////cout<<j<<"\t"<<Atoms[ATOM->neighlist[j]].x<<"\t"<<Atoms[ATOM->neighlist[j]].y<<"\t";
					//////////cout<<Y<<"\n";

					////////if(Y<Y_MIN)
					////////{
					////////	;
					////////	Y_MIN=Y;
					////////        DIS_atom=ATOM->neighlist[j];
					////////	circx=x+L.x;
					////////	circy=y+L.y;

					////////}
					}
				}//j loop 
				flag=1;
				for(k=0;k<ATOM->conti;k++)
				{
					if(ATOM->contigous[k]==DIS_atom)
					{
						flag=0;
						break;
					}
				}

				if(flag)
				{
					ATOM->contigous[ATOM->conti]=DIS_atom;
					ATOM->edge_index[D->B]++;
					ATOM->edge_index[ATOM->conti]++;
				}
				else 
				{
					ATOM->edge_index[D->B]++;
					ATOM->edge_index[k]++;
				}

				//cout<<DIS_atom<<"\n";
				//cout<<Atoms[DIS_atom].x<<"\t"<<Atoms[DIS_atom].y<<"\n";
				delunay *temp;
				temp=D;
				while(1)
				{
					if(temp->next)
						temp=temp->next;
					else 
						break;

				}
				temp->next= new delunay;
				if(flag)
					temp->next->A=ATOM->conti;
				else 
					temp->next->A=k;
				temp->next->B=D->B;
				temp->next->circum_x=circx;
				temp->next->circum_y=circy;
				if(flag)
					ATOM->conti++;		
			}// D->B loop
		}//For atom i in conti
	}//loop  over conti
}//end 

int main()
{

	int nAtoms=0;
	atom *Atoms;
       	std::ifstream infile("dat");
        infile>>nAtoms;
        infile>>box;
	twob=2*box;
	density=nAtoms/(twob*twob);
	int SAM=0;
        Atoms = new (nothrow) atom[nAtoms];
        long double b,c,d;				
        nAtoms=0;
	ofstream vor;
	vor.open("vor");
        while(infile>>b>>c>>d) 
	{ 
            Atoms[nAtoms].x=b;
            Atoms[nAtoms].y=c;
	    Atoms[nAtoms].radius=d;
	    //cout<<nAtoms<<"\t"<<b<<"\t"<<c<<"\n";
            nAtoms++;
	}
	update_neighbours(Atoms,nAtoms);
	double area=0;
	for(SAM=0;SAM<nAtoms;SAM++)
	{
	    //  for(int i=0;i<Atoms[SAM].neighbours;i++)
	    //  	cout<<i<<"\t"<<Atoms[Atoms[SAM].neighlist[i]].x<<"\t"<<Atoms[Atoms[SAM].neighlist[i]].y<<"\n";
		first_delunay(&(Atoms[SAM]),Atoms);
		complete_del(&(Atoms[SAM]),Atoms,nAtoms);
		delunay *D;
		D=Atoms[SAM].D.initial;
		int count=0;
		while(1)
		{
			print_delunay(&(Atoms[SAM]),D,Atoms);
			count++;
			if(D->next)
			{
				D=D->next;
			}
			else
				break;
		}
		D=Atoms[SAM].D.initial;
	//cout<<"here\n";
		double area_s=0;
                for(int i=0;i<Atoms[SAM].conti;i++)
                {
                	D=Atoms[SAM].D.initial;
                	vor<<"\n";
                	//cout<<i<<"\n";
			double a,b,p,q,x,y;
			x=Atoms[SAM].x;
			y=Atoms[SAM].y;
			int flaga=1;
			int flagb=1;
			vor<<x<<"\t"<<y<<"\t";
                	while(1)
                	{
                		if(D->A==i)
                		{
                			vor<<D->circum_x<<"\t"<<D->circum_y<<"\t";
		////////		if(flaga)
		////////		{
		////////			a=D->circum_x;
		////////			b=D->circum_y;
		////////			flaga=0;
		////////		}
                		}
                		if(D->B==i)
                		{
                			vor<<D->circum_x<<"\t"<<D->circum_y<<"\t";
		////////		if(flagb)
		////////		{
		////////			p=D->circum_x;
		////////			q=D->circum_y;
		////////			flagb=0;
		///			}
                		}
                		if(D->next)
                			D=D->next;
                		else 
                		{
                			break;
                		}
                		
                	}
			//area_s=area_s+0.5*abs((x-p)*(b-y)-(x-a)*(q-y));
			//cout<<0.5*abs((x-p)*(b-y)-(x-a)*(q-y))<<"\n";
			//cout<<a<<"\t"<<b<<"\t"<<p<<"\t"<<q<<"\n";
                	
                }

		while(1)
		{
			delunay *temp;
			temp=D;
			if(D->next)
			{
				D=D->next;
				delete temp;
			}
			else
			{
				delete temp;
				break;
			}

		}
	}
	delete[] Atoms;
	return 0;
}
