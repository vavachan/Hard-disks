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
	int neighlist[200];
	int contigous[200];
	int edge_index[100]={0};
	int neighbours=0;
	int conti=0;
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
	//t adcout<<R_CUT<<"\n";
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
	}
}
void first_delunay(atom *ATOM,atom Atoms[])
{
	double drx,dry,dr;
	double min=2.*twob*twob;
	int nearest,flag;
	for(int i=0;i<ATOM->neighbours;i++)
	{
		drx=ATOM->x-Atoms[ATOM->neighlist[i]].x;
		dry=ATOM->y-Atoms[ATOM->neighlist[i]].y;
		drx=(drx-(twob*lround(drx/twob)));
		dry=(dry-(twob*lround(dry/twob)));
		dr=drx*drx+dry*dry;
		//cout<<dr<<"\t"<<i<<"\n";
		if(dr<min)
		{
			min=dr;
			nearest=ATOM->neighlist[i];
		}
	}
	ATOM->contigous[ATOM->conti]=nearest;
	ATOM->edge_index[ATOM->conti]++;
	ATOM->conti++;
	//cout<<"nea\t"<<nearest<<"\t"<<min<<"\n";
	//cout<<Atoms[nearest].x<<"\t"<<Atoms[nearest].y<<"\n";
	double DIS_MIN=box;
	int DIS_atom;
	double dis;
	double circx;
	double circy;
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
			double A=ax*ax+ay*ay;				
			double B=bx*bx+by*by;
			double x=0.5*(ay*B-by*A)/(ay*bx-ax*by);
			double denom=2.*(ay*bx-ax*by);
			double y=-1.*ax/ay*x+A/(2.*ay);
			dis=sqrt(pow(x-ax,2)+pow(y-ay,2));
			if(DIS_MIN > dis)
			{
				DIS_MIN=dis;
				DIS_atom=ATOM->neighlist[i];
				circx=x+L.x;
				circy=y+L.y;
	//			cout<<Atoms[ATOM->neighlist[i]].x<<"\t"<<Atoms[ATOM->neighlist[i]].y<<"\n";
	//			cout<<dis<<"\n";
			}
		}

	}
////////cout<<Atoms[ATOM->neighlist[101]].x<<"\t"<<Atoms[ATOM->neighlist[101]].y<<"\n";
////////cout<<ATOM->neighlist[101]<<"\n";
////////cout<<DIS_atom<<"\n";
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
						double A=ax*ax+ay*ay;				
						double B=bx*bx+by*by;
						double x=0.5*(ay*B-by*A)/(ay*bx-ax*by);
						double denom=2.*(ay*bx-ax*by);
						double y=-1.*ax/ay*x+A/(2.*ay);
						double dis=sqrt(pow(x-ax,2)+pow(y-ay,2));
						double Y=sqrt(pow(x-ax/2.,2)+pow(y-ay/2.,2));
						int sign_C;
						if((y-(m*x+C))<0.)
							sign_C=-1;
						else
							sign_C=1;
						if(sign_C!=sign_N)
							Y=-1.*Y;
					//			cout<<"circum="<<dis<<"\n";
						if(Y<Y_MIN)
						{
							;
							Y_MIN=Y;
						        DIS_atom=ATOM->neighlist[j];
							circx=x+L.x;
							circy=y+L.y;

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
						double A=ax*ax+ay*ay;				
						double B=bx*bx+by*by;
						double x=0.5*(ay*B-by*A)/(ay*bx-ax*by);
						double denom=2.*(ay*bx-ax*by);
						double y=-1.*ax/ay*x+A/(2.*ay);
						double dis=sqrt(pow(x-ax,2)+pow(y-ay,2));
						double Y=sqrt(pow(x-ax/2.,2)+pow(y-ay/2.,2));
						int sign_C;
						if((y-(m*x+C))<0.)
							sign_C=-1;
						else
							sign_C=1;
						if(sign_C!=sign_N)
							Y=-1.*Y;

				        	//cout<<j<<"\t"<<Atoms[ATOM->neighlist[j]].x<<"\t"<<Atoms[ATOM->neighlist[j]].y<<"\t";
						//cout<<Y<<"\n";

						if(Y<Y_MIN)
						{
							;
							Y_MIN=Y;
						        DIS_atom=ATOM->neighlist[j];
							circx=x+L.x;
							circy=y+L.y;

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
        long double b,c;				
        nAtoms=0;
	ofstream vor;
	vor.open("vor");
        while(infile>>b>>c) 
	{ 
            Atoms[nAtoms].x=b;
            Atoms[nAtoms].y=c;
            nAtoms++;
	}
	update_neighbours(Atoms,nAtoms);
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
////////	cout<<"here\n";
        	for(int i=0;i<Atoms[SAM].conti;i++)
        	{
        		D=Atoms[SAM].D.initial;
        		vor<<"\n";
        		while(1)
        		{
        			if(D->A==i)
        			{
        				vor<<D->circum_x<<"\t"<<D->circum_y<<"\t";
        			}
        			if(D->B==i)
        			{
        				vor<<D->circum_x<<"\t"<<D->circum_y<<"\t";
        			}
        			if(D->next)
        				D=D->next;
        			else 
        			{
        				break;
        			}
        			
        		}
        		
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
