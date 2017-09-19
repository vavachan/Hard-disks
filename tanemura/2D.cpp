#include<iostream>
#include<math.h>
#include<fstream>
using namespace std;
double box=0;
double twob;
struct half_edge 
{
	struct half_edge *twin=NULL;
	struct half_edge *next=NULL;
	struct vertice *origin=NULL;
	struct face *facing=NULL;
};
//definition of vertice
struct vertice
{
	double x=0;
	double y=0;
	struct half_edge *leaving=NULL;
};
//definition of face
struct face
{
	struct half_edge *edge=NULL;
};
struct atom;
struct delunay
{
	int A;
	int B;
	delunay *next=NULL;
};
class set_of_delunay
{
	public :
		delunay *initial=NULL;
		void update_list(delunay *,delunay *);
};
void set_of_delunay::update_list(delunay *start, delunay *D)
{
	if(start->next)
		update_list((start->next),D);
	else
		start->next=D;
}
struct atom 
{
	double x=0;
	double y=0;
	int neighlist[100];
	int contigous[100];
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
			if(dr<600.*600.)
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
	cout<<Atoms[nearest].x<<"\t"<<Atoms[nearest].y<<"\n";
	double DIS_MIN=box;
	int DIS_atom;
	double dis;
	for(int i=0;i<ATOM->neighbours;i++)
	{
		if(i!=nearest)
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
		}
		if(DIS_MIN > dis)
		{
			DIS_MIN=dis;
			DIS_atom=ATOM->neighlist[i];
		}

	}
	ATOM->contigous[ATOM->conti]=DIS_atom;
	ATOM->edge_index[ATOM->conti]++;
	ATOM->conti++;
	ATOM->D.initial= new delunay;
	ATOM->D.initial->A=0;
	ATOM->D.initial->B=1;
	//cout<<ATOM->D.initial->A<<"\n";
	//cout<<ATOM->D.initial<<"\n";
}
void print_delunay(atom *ATOM,delunay *D,atom Atoms[])
{
	cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<Atoms[ATOM->contigous[D->A]].x<<"\t"<<Atoms[ATOM->contigous[D->A]].y<<"\n";
	cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<Atoms[ATOM->contigous[D->B]].x<<"\t"<<Atoms[ATOM->contigous[D->B]].y<<"\n";
	cout<<Atoms[ATOM->contigous[D->A]].x<<"\t"<<Atoms[ATOM->contigous[D->A]].y<<"\t"<<Atoms[ATOM->contigous[D->B]].x<<"\t"<<Atoms[ATOM->contigous[D->B]].y<<"\n";;
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
				//		cout<<"sign_N="<<sign_N<<"\n";
					if(sign!=sign_N)
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


				//			cout<<"circum="<<dis<<"\n";
						if(Y<Y_MIN)
						{
							;
							Y_MIN=Y;
					////////	dis_MIN=dis;
						        DIS_atom=ATOM->neighlist[j];
			////////  			cout<<"circum="<<Y_MIN<<"\n";
			////////  			cout<<"circum_atom="<<DIS_atom<<"\n";

						}
					}
					
				}//j loop 
				int flag=1;
				int k;
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
				if(flag)
					ATOM->conti++;		
			}// D->A loop
			else
			{
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
					if(sign!=sign_N)
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


				//		cout<<"circum="<<dis<<"\n";
						if(Y<Y_MIN)
						{
							;
							Y_MIN=Y;
					////////	dis_MIN=dis;
						        DIS_atom=ATOM->neighlist[j];
				//			cout<<"circum="<<dis<<"\n";
				//			cout<<"circum_atom="<<DIS_atom;

						}
					}
				}//j loop 
				int flag=1;
				int k;
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
        Atoms = new (nothrow) atom[nAtoms];
        long double b,c;				
        nAtoms=0;
        while(infile>>b>>c) 
	{ 
            Atoms[nAtoms].x=b;
            Atoms[nAtoms].y=c;
          //  cout<<nAtoms<<"\t"<<Atoms[nAtoms].x<<"\t"<<Atoms[nAtoms].y<<"\n";
            nAtoms++;
	}
	update_neighbours(Atoms,nAtoms);
	//cout<<Atoms[12].x<<"\t"<<Atoms[12].y<<"\t";
////////for(int i=0;i<Atoms[83].neighbours;i++)
////////{
//////////	cout<<Atoms[83].neighlist[i]<<"\n";
////////	cout<<Atoms[Atoms[83].neighlist[i]].x<<"\t"<<Atoms[Atoms[83].neighlist[i]].y<<"\n";
////////}
	first_delunay(&(Atoms[84]),Atoms);
	//print_delunay(&(Atoms[83]),Atoms[83].D.initial,Atoms);
	complete_del(&(Atoms[84]),Atoms,nAtoms);
////////complete_delunay(&(Atoms[83]),Atoms[83].D.initial,Atoms,nAtoms);
////////complete_delunay(&(Atoms[83]),Atoms[83].D.initial->next->next,Atoms,nAtoms);
////////complete_delunay(&(Atoms[83]),Atoms[83].D.initial->next->next->next,Atoms,nAtoms);
	delunay *D;
	D=Atoms[84].D.initial;
	while(1)
	{
		print_delunay(&(Atoms[84]),D,Atoms);
		if(D->next)
		{
			D=D->next;
		}
		else
			break;
	}
////////print_delunay(&(Atoms[83]),Atoms[83].D.initial,Atoms);
////////print_delunay(&(Atoms[83]),Atoms[83].D.initial->next,Atoms);
////////print_delunay(&(Atoms[83]),Atoms[83].D.initial->next->next,Atoms);
////////print_delunay(&(Atoms[83]),Atoms[83].D.initial->next->next->next,Atoms);
	//print_delunay(&(Atoms[17]),Atoms[17].D.initial->next->next->next->next,Atoms);
	//cout<<"hehre\n";
//      cout<<Atoms[13].x<<"\t"<<Atoms[13].y<<"\t"<<Atoms[13].D.initial->A->x<<"\t"<<Atoms[13].D.initial->A->y<<"\n";
//      cout<<Atoms[13].x<<"\t"<<Atoms[13].y<<"\t"<<Atoms[13].D.initial->B->x<<"\t"<<Atoms[13].D.initial->B->y<<"\n";
//      cout<<Atoms[13].D.initial->A->x<<"\t"<<Atoms[13].D.initial->A->y<<"\t"<<Atoms[13].D.initial->B->x<<"\t"<<Atoms[13].D.initial->B->y<<"\n";
	
	//for(int i=0;i<nAtoms;i++)
		//cout<<i<<"\t"<<Atoms[i].neighbours<<"\n";
		//for(int j=0;j<Atoms[i].neighbours;j++)

	delete[] Atoms;
	return 0;
}
