#include<iostream>
#include<math.h>
#include<fstream>
using namespace std;
#define distance(x,y) x*x+y*y
int dim=3;
double box=0;
double twob;
double density;
double r_cut=0.0;
struct atom;
struct face;
struct vertice;
struct site
{
	double x=0;
	double y=0;
	struct face *F=NULL;
};
void display_SITE(struct site *p)
{
	cout<<p->x<<"\t"<<p->y<<"\n";
}
//definition of half-edge 
struct half_edge 
{
	struct half_edge *twin=NULL;
	struct half_edge *next=NULL;
	struct vertice *origin=NULL;
	struct face *facing=NULL;
};
//definition of vertice
//definition of face
struct face
{
	struct site *p=NULL;
	struct half_edge *edge=NULL;
};

class vert_list
{
	public:
		void delete_vertice(vertice *EV,vertice *v);
		vertice* insert_vertice(vertice *,vertice *);
		void display_vertex(vertice *);
		void display_conn(vertice *ve);
}*V;
struct delunay
{
	int A;
	int B;
	int C;
	double circum_x=0.;
	double circum_y=0.;
	double Ax=0.;
	double Ay=0.;
	double Bx=0.;
	double By=0.;
	delunay *next=NULL;
};
struct vertice
{
	struct site *p=NULL;
	struct half_edge *leaving=NULL;
	struct vertice *next=NULL;
	struct vertice *prev=NULL;
	int A,B,C;
	delunay *D;
	int is_void=0;
	int cluster_index=-1;
	//vert_list *V;
	vertice *neib_vert[10];
	int neib_ed[10];
	int v_neigh_count=0;
}*start;
void add_connected(vertice *focus,vertice *add,int binv)
{
	int flag=1;
	int i;
	for(i=0;i<focus->v_neigh_count;i++)
	{
		if(focus->neib_vert[i]==add)
		{
			flag=0;
			break;
		}
	}
	if(flag)
	{
		focus->neib_vert[focus->v_neigh_count]=add;
		focus->neib_ed[focus->v_neigh_count]=binv;
		focus->v_neigh_count=focus->v_neigh_count+1;
	}
	else 
		focus->neib_ed[i]=binv;

}
int compare(struct site *p1,struct site *p2)
{
        double DX,DY;
        DX=p1->x-p2->x;
        DY=p1->y-p2->y;
	//cout<<DX<<"\t"<<DY<<"\n";	
        DX=(DX-(twob*lround(DX/twob)));
        DY=(DY-(twob*lround(DY/twob)));
	//cout<<DX<<"\t"<<DY<<"\n";	
////////if( DY > 0.  )
////////	return 1;
////////else if ((  DY == 0. ) && ( DX > 0. ) )
////////	return 1;
////////else if ( ( DY==0.) && ( DX == 0. ) )
////////	return 0;
////////else 
////////	return -1;
	if( abs(DX) <0.0000001  && abs(DX) >0.0000001 )
	{
		cout<<p1->x<<"\t"<<p1->y<<"\n";
		cout<<p2->x<<"\t"<<p2->y<<"\n";
	}
	if( abs(DY) <0.0000001  && abs(DX) <0.0000001 )
		return 0;
	else if ( (p1->y == p2->y) && ( p1->x > p2->x ) )
		return 1;
	else if( p1->y > p2->y ) 
		return 1;
	else 
		return -1;
}
vertice* vert_list::insert_vertice(vertice *EV,vertice *v)
{
	int flag;
	flag=compare(EV->p,v->p);
///     cout<<"insert\n";
///     display_SITE(v->p);
///     cout<<"vetfoc\n";
///     display_SITE(EV->p);
///     cout<<"flag\t"<<flag<<"\n";
	if(flag==1)
	{
		if(EV->next)
		{
	//		cout<<"everytime?\n";
			insert_vertice(EV->next,v);
		}
		else 
		{
			EV->next=v;
			v->prev=EV;
	///		cout<<"here?\n";
			return v;
		}
	}
	else if (flag==-1)
	{
		v->next=EV;
		if(EV->prev)
		{
			v->prev=EV->prev;
			EV->prev->next=v;
		}
		else 
			start=v;
		EV->prev=v;
	//	cout<<"here2?\n";
		return v;
	}
	else if (flag==0)
	{
		//cout<<"vertice already exists\n";
		delete v->p;
		delete v;
		return EV;
	}
}

void vert_list::delete_vertice(vertice *EV,vertice *v)
{
    int flag;
    flag=compare(EV->p,v->p);
    if (flag==0)
    {
        if(EV->prev)
            EV->prev->next=EV->next;
        else
            start=EV->next;
        if(EV->next)
            EV->next->prev=EV->prev;
        else
            EV->prev->next=NULL;
        //cout<<"herein delete\n";
        //display_SITE(v->p);
      //for(int n=0; n<v->v_neigh_count; n++)
      //{
      //    if(v->neib_vert[n])
      //        for(int m=0; m<v->neib_vert[n]->v_neigh_count; m++)
      //        {
      //            if(v->neib_vert[n]->neib_vert[m])
      //                if(!(compare(v->neib_vert[n]->neib_vert[m]->p,v->p)))
      //                {
      //                    v->neib_vert[n]->neib_vert[m]=NULL;
      //                }
      //        }
      //}
        //delete v->p;
        //delete v;
    }
    else
    {
        if(EV->next)
            delete_vertice(EV->next,v);
        else
            return;
    }
}
void vert_list::display_vertex(vertice *start)
{
	display_SITE(start->p);
	if(start->next)
		display_vertex(start->next);
	else 
		return;
	;
}
void vert_list::display_conn(vertice *ve)
{
	if(ve->v_neigh_count==2)
	{
		cout<<"#"<<ve->v_neigh_count<<"\n";
		for(int i=0;i<ve->v_neigh_count;i++)
		{
			double Sx,Sy;
			Sx=ve->neib_vert[i]->p->x-ve->p->x;
			Sy=ve->neib_vert[i]->p->y-ve->p->y;
			Sx=(Sx-(twob*lround(Sx/twob)));
			Sy=(Sy-(twob*lround(Sy/twob)));
			cout<<ve->p->x<<"\t"<<ve->p->y<<"\t"<<Sx+ve->p->x<<"\t"<<Sy+ve->p->y<<"\n";
		}
	}
        if(ve->next)
        	display_conn(ve->next);
        else 
		return;

}
void display_a_edge(struct half_edge *first,struct half_edge *E)
{
	ofstream DEL;
	//cout<<"d therer ana\n";
	DEL.open("dcel",ios::app);
	//DEL<<"#"<<E->facing->p->x<<" "<<E->facing->p->y<<"\n";
	//cout<<E->origin->p<<"\n";
	DEL<<E->origin->p->x<<"\t"<<E->origin->p->y<<"\t"<<E->twin->origin->p->x<<"\t"<<E->twin->origin->p->y<<"\n";
	//cout<<"d trer ana\n";
	DEL.close();
	if(E->next)
	{
		if(first==E->next)
		{
			return;
		}
		display_a_edge(first,E->next);
	}
	else 
	{
		//cout<<"wehnt this way\n";
		return;
	}

}
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
	int bondinvoid[500];
	int edge_index[500]={0};
	int neighbours=0;
	int conti=0;
	struct face *F=NULL;
	double radius=1.;
	double ignore=0;
	set_of_delunay D;
};
struct atom3D 
{
	double x=0;
	double y=0;
	double z=0;
	int neighlist[500];
	int contigous[500];
	int bondinvoid[500];
	int edge_index[500]={0};
	int neighbours=0;
	int conti=0;
	struct face *F=NULL;
	double radius=1.;
	set_of_delunay D;
};
double perimeter(double Vx,double Vy,double E12x,double E12y,double Ax,double Ay,double Arad)
{
	E12x=E12x-Vx;
	E12y=E12y-Vy;
	Ax=Ax-Vx;
	Ay=Ay-Vy;
	double areac,areat;
	double perimeter;
	double dis,disint,dise,disv;
	dis=sqrt(distance((E12x-Ax),(E12y-Ay)));
	if(dis<Arad)
	{
		disint=sqrt(Arad*Arad-dis*dis);
		dise=sqrt(distance(E12x,E12y));
		disv=sqrt(distance(Ax,Ay));
		double Aintx,Ainty,X1,Y1,Aint1x,Aint1y,X2,Y2;
		Aintx=(dise-disint)/dise*E12x;
		Ainty=(dise-disint)/dise*E12y;
		Aint1x=(disv-Arad)/disv*Ax;
		Aint1y=(disv-Arad)/disv*Ay;
		X1=Aintx-Ax;
		Y1=Ainty-Ay;
		X2=Aint1x-Ax;
		Y2=Aint1y-Ay;
		double theta=acos((X2*X1+Y1*Y2)/(Arad*Arad));
		if((X2*X1+Y1*Y2)/(Arad*Arad)>1.)
			theta=0.;
		perimeter=2.*M_PI*Arad*(theta/(2.*M_PI));
		return perimeter;
	}
	else
	{
		disv=sqrt(distance(Ax,Ay));
		double Aintx,Ainty,X1,Y1,Aint1x,Aint1y,X2,Y2;
		Aint1x=(disv-Arad)/disv*Ax;
		Aint1y=(disv-Arad)/disv*Ay;
		X1=E12x-Ax;
		Y1=E12y-Ay;
		X2=Aint1x-Ax;
		Y2=Aint1y-Ay;
		double theta=acos((X2*X1+Y1*Y2)/(Arad*dis));
		if((X2*X1+Y1*Y2)/(Arad*dis)>1.)
			theta=0.;
		perimeter=2.*M_PI*Arad*(theta/(2.*M_PI));
		return perimeter;
	}

}
double area_trangle(double Vx,double Vy,double E12x,double E12y,double Ax,double Ay,double Arad,double originx=0.,double originy=0.)
{
	E12x=E12x-Vx;
	E12y=E12y-Vy;
	Ax=Ax-Vx;
	Ay=Ay-Vy;
	double areac,areat;
	double dis,disint,dise,disv;
	dis=sqrt(distance((E12x-Ax),(E12y-Ay)));
	if(dis<Arad)
	{
		//cout<<dis<<"\t"<<Arad<<"\n";
		disint=sqrt(Arad*Arad-dis*dis);
		dise=sqrt(distance(E12x,E12y));
		disv=sqrt(distance(Ax,Ay));
		//E12x=(E12x-(twob*lround(E12x/twob)));
		//E12y=(E12y-(twob*lround(E12y/twob)));
		double Aintx,Ainty,X1,Y1,Aint1x,Aint1y,X2,Y2;
		Aintx=(dise-disint)/dise*E12x;
		Ainty=(dise-disint)/dise*E12y;
		//cout<<Aintx+Vx+originx<<"\t"<<Ainty+Vy+originy<<"\n";
		Aint1x=(disv-Arad)/disv*Ax;
		Aint1y=(disv-Arad)/disv*Ay;
		X1=Aintx-Ax;
		Y1=Ainty-Ay;
		X2=Aint1x-Ax;
		Y2=Aint1y-Ay;
		double theta=acos((X2*X1+Y1*Y2)/(Arad*Arad));
		if((X2*X1+Y1*Y2)/(Arad*Arad)>1.)
			theta=0.;
		//double theta=acos(-1.0);

		//cout<<theta<<"=theta\t"<<(X2*X1+Y1*Y2)/(Arad*Arad)+Y1*Y2<<"\n";
		//cout<<Arad<<"\n";
		areac=M_PI*Arad*Arad*(theta/(2.*M_PI));
		double a,b,p,q;
		a=Ax;
		b=Ay;
		p=Aintx;
		q=Ainty;
	////////cout<<"#\n";
	////////cout<<Vx+originx<<"\t"<<Vy+originy<<"\t"<<a+Vx+originx<<"\t"<<b+Vy+originy<<"\n";
	////////cout<<Vx+originx<<"\t"<<Vy+originy<<"\t"<<p+Vx+originx<<"\t"<<q+Vy+originy<<"\n";
	////////cout<<a+Vx+originx<<"\t"<<b+Vy+originy<<"\t"<<p+Vx+originx<<"\t"<<q+Vy+originy<<"\n";
	////////cout<<"#\n";
		areat=0.5*abs(a*q-b*p);
		//cout<<areat<<"\t"<<areac<<"\n";
		return areat-areac;
	}
	else
	{
		disv=sqrt(distance(Ax,Ay));
		//E12x=(E12x-(twob*lround(E12x/twob)));
		//E12y=(E12y-(twob*lround(E12y/twob)));
		double Aintx,Ainty,X1,Y1,Aint1x,Aint1y,X2,Y2;
	////////Aintx=Arad/dis*E12x;
	////////Ainty=Arad/dis*E12y;
		Aint1x=(disv-Arad)/disv*Ax;
		Aint1y=(disv-Arad)/disv*Ay;
		//cout<<Aint1x+Vx+originx<<"\t"<<Aint1y+Vy+originy<<"\n";
		//cout<<Aint1x<<"\t"<<Aint1y<<"\n";
		X1=E12x-Ax;
		Y1=E12y-Ay;
		X2=Aint1x-Ax;
		Y2=Aint1y-Ay;
		//cout<<X1<<"\t"<<Y1<<"\n";
		//cout<<X2<<"\t"<<Y2<<"\n";
		//cout<<sqrt(distance(X1,Y1))<<"\t"<<dis<<"\n";;
		//cout<<(X1*X2+Y1*Y2)/(Arad*dis)<<"\n";
		double theta=acos((X2*X1+Y1*Y2)/(Arad*dis));
		//cout<<"theta="<<theta<<"\n";
	        if((X2*X1+Y1*Y2)/(Arad*dis)>1.)
		{
			//cout<<theta<<"\n";
	        	theta=0.;
		}
		areac=M_PI*Arad*Arad*(theta/(2.*M_PI));
		//cout<<"areacircle="<<areac<<"\n";
		double a,b,p,q;
		a=Ax;
		b=Ay;
		p=E12x;
		q=E12y;
		areat=0.5*abs(a*q-b*p);
		//cout<<"areat="<<areat<<"\n";
		return areat-areac;
	}

}
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
	int binv1=0;
	int binv2=0;
	double midax=0.;
	double miday=0.;
	double midbx=0.;
	double midby=0.;
	double lmin=0;
	for(int i=0;i<ATOM->neighbours;i++)
	{
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
			binv1=0;
			midax=x+ATOM->x;
			miday=y+ATOM->y;
			lmin=l;
			if(DIS>rS+rA+2*r_cut)
			{
				binv1=1;
			}
		}

	}
	ATOM->bondinvoid[ATOM->conti]=binv1;
	ATOM->contigous[ATOM->conti]=nearest;
	ATOM->edge_index[ATOM->conti]++;
	ATOM->conti++;
	if(lmin<0.)
	{
		ATOM->ignore=1;
		//return;
	}
////////cout<<lmin<<" =l\n";
////////cout<<Atoms[nearest].x<<"\t"<<Atoms[nearest].y<<"\n";
////////cout<<"nea\t"<<nearest<<"\t"<<min<<"\n";
////////cout<<Atoms[nearest].x<<"\t"<<Atoms[nearest].y<<"\n";
	double DIS_MIN=box*box;
	int DIS_atom;
	double dis;
	double X,Y;
	double circx,circy;
	int temp;
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
			double rA,rS,rB;
			double XA,YA,XB,YB;
			double xA,yA,xB,yB;
			double l;
			double DISA;
			double DISB;
			double MA,MB,INMA,INMB;
			double CA,CB;
			double tan_sq;
			XA=ax;
			YA=ay;
			XB=bx;
			YB=by;
			DISA=sqrt(XA*XA+YA*YA);
			DISB=sqrt(XB*XB+YB*YB);
			MA=YA/XA;
			MB=YB/XB;
			INMA=-1./MA;
			INMB=-1./MB;
			rA=M.radius;
			rB=R.radius;
			rS=L.radius;
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
			tan_sq=X*X+Y*Y-rS*rS;
			double Y_AXIS=sqrt(pow(X-xA,2)+pow(Y-yA,2));
			int sign_C;
			int sign_N;
			double m;
			m=YA/XA;
			if((by-(m*bx))<0.)
				sign_N=-1;
			else 
				sign_N=1;
			if((Y-(m*X))<0.)
				sign_C=-1;
			else
				sign_C=1;
			X=X+L.x;
			Y=Y+L.y;
			if(sign_C!=sign_N)
			{
			////////cout<<bx<<"\t"<<by<<"\t"<<by-(m*bx)<<"\n";
			////////cout<<X<<"\t"<<Y<<"\t"<<Y-(m*X)<<"\n";
			////////cout<<sign_N<<"\t"<<sign_C<<"\n";
			////////cout<<"this didn't happen\n";
				Y_AXIS=-1.*Y_AXIS;
			}
			if(DIS_MIN > tan_sq)
			{
                		DIS_MIN=tan_sq;
                		DIS_atom=ATOM->neighlist[i];
	        		circx=X;
	        		circy=Y;
				binv1=0;
				midbx=xB+L.x;
				midby=yB+L.y;
				lmin=l;
				//cout<<DIS_atom<<"\n";
				if(DISB>(rB+rS+2*r_cut))
				{
					//cout<<"here\n";
					binv1=1;
				}
			////////if(Y_AXIS < 0.)
			////////{
			////////	temp=1;
			////////	cout<<"here\n";
			////////	cout<<X<<"\t"<<Y<<"\n";
			////////	cout<<bx+L.x<<"\t"<<by+L.y<<"\n";
			////////	//ATOM->bondinvoid[ATOM->conti-1]=1;
			////////}
			////////else 
			////////{
			////////	//cout<<"here2\n";
			////////	temp=0;
			////////}
			}
		}

	}
////////if(temp)
////////	ATOM->bondinvoid[ATOM->conti-1]=1;
////////cout<<lmin<<" =l\n";
////////cout<<Atoms[DIS_atom].x<<"\t"<<Atoms[DIS_atom].y<<"\n";
	ATOM->bondinvoid[ATOM->conti]=binv1;
	ATOM->contigous[ATOM->conti]=DIS_atom;
	ATOM->edge_index[ATOM->conti]++;
	ATOM->conti++;
	ATOM->D.initial= new delunay;
	ATOM->D.initial->A=0;
	ATOM->D.initial->B=1;
	ATOM->D.initial->circum_x=circx;
	ATOM->D.initial->circum_y=circy;
	ATOM->D.initial->Ax=midax;
       	ATOM->D.initial->Ay=miday;
	ATOM->D.initial->Bx=midbx;
       	ATOM->D.initial->By=midby;
}
void print_delunay(atom *ATOM,delunay *D,atom Atoms[])
{
	//cout<<"what\n";
	double Sx,Sy;
	double Px,Py;
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
	int binv;
	int binv2;
	double lmin;
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
			        //cout<<i<<"\n";
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
				int temp1;
				double midbx=0.;
				double midby=0.;
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
						XA=ax;
						YA=ay;
						XB=bx;
						YB=by;
						DISA=sqrt(XA*XA+YA*YA);
						DISB=sqrt(XB*XB+YB*YB);
						MA=YA/XA;
						MB=YB/XB;
						INMA=-1./MA;
						INMB=-1./MB;
						rA=M.radius;
						rB=R.radius;
						rS=L.radius;
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
						tan_sq=X*X+Y*Y-rS*rS;
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
						int sign_circ;
						if(Y_AXIS<Y_MIN)
						{
							;
							Y_MIN=Y_AXIS;
						        DIS_atom=ATOM->neighlist[j];
							circx=X;
							circy=Y;
							binv=0;
							midbx=xB+L.x;
							midby=yB+L.y;
							lmin=l;
							if(DISB>rS+rB+2*r_cut )
							{
								binv=1;
							}
						////////cout<<"compa\n";
						////////cout<<X<<"\t"<<Y<<"\n";
						////////cout<<D->circum_x<<"\t"<<D->circum_y<<"\n";
	////////					if(((D->circum_y-L.y)-(m*(D->circum_x-L.x)+C))<0.)
	////////						sign_circ=-1;
	////////					else
	////////						sign_circ=1;
	////////					if(sign_C==sign_circ)	
	////////					{
	////////						temp1=1;
	////////					}
	////////					else
	////////					{
	////////						temp1=0;
	////////					}


						}
					}
					
				}//j loop 
				flag=1;
				//cout<<ATOM->contigous[D->A]<<"\t"<<temp1<<"\n";
				//cout<<DIS_atom<<"\t"<<binv<<"\n";
				if(lmin<0.)
				{
					ATOM->ignore=1;
					//return;
				}
	////////		if(temp1)
	////////			ATOM->bondinvoid[D->A]=1;
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
					ATOM->bondinvoid[ATOM->conti]=binv;
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
				{
					temp->next->A=ATOM->conti;
					temp->next->Ax=midbx;
					temp->next->Ay=midby;
				}
				else 
				{
					temp->next->A=k;
					temp->next->Ax=midbx;
					temp->next->Ay=midby;
				}
				temp->next->B=D->A;
				temp->next->Bx=D->Ax;
				temp->next->By=D->Ay;
				temp->next->circum_x=circx;
				temp->next->circum_y=circy;
				if(flag)
					ATOM->conti++;		
			}// D->A loop
			else
			{
			        //cout<<"hereb\n";
			        ///cout<<i<<"\n";
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
				int temp1;
				double midbx=0.;
				double midby=0.;
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
						XA=ax;
						YA=ay;
						XB=bx;
						YB=by;
						DISA=sqrt(XA*XA+YA*YA);
						DISB=sqrt(XB*XB+YB*YB);
						MA=YA/XA;
						MB=YB/XB;
						INMA=-1./MA;
						INMB=-1./MB;
						rA=M.radius;
						rB=R.radius;
						rS=L.radius;
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
						tan_sq=X*X+Y*Y-rS*rS;
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
						int sign_circ;
						if(Y_AXIS<Y_MIN)
						{
							;
							Y_MIN=Y_AXIS;
						        DIS_atom=ATOM->neighlist[j];
							circx=X;
							circy=Y;
							midbx=xB+L.x;
							midby=yB+L.y;
							binv=0;
							lmin=l;
							if(DISB>rS+rB+2*r_cut)
							{
								binv=1;
							}
	////////					if(((D->circum_y-L.y)-(m*(D->circum_x-L.x)+C))<0.)
	////////						sign_circ=-1;
	////////					else
	////////						sign_circ=1;
	////////					if(sign_C==sign_circ)	
	////////					{
	////////						temp1=1;
	////////					}
	////////					else
	////////					{
	////////						temp1=0;
	////////					}
						////////if(Y_AXIS < 0.)
						////////	temp1=1;
						////////else
						////////	temp1=0;
						////////if(Y_AXIS < 0. && !ATOM->bondinvoid[D->A])
						////////	ATOM->bondinvoid[D->A]=1;

						}
					}
				}//j loop 
				flag=1;
				//cout<<ATOM->contigous[D->B]<<"\n";
	////////		if(temp1)
	////////			ATOM->bondinvoid[D->B]=1;
				//cout<<lmin<<"\n";
				//cout<<Atoms[DIS_atom].x<<"\t"<<Atoms[DIS_atom].y<<"\n";
				if(lmin<0.)
				{
					ATOM->ignore=1;
					//return;
				}
				if(lmin<0.)
				{

				}
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
					ATOM->bondinvoid[ATOM->conti]=binv;
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
				{
					temp->next->A=ATOM->conti;
					temp->next->Ax=midbx;
					temp->next->Ay=midby;
				}
				else 
				{
					temp->next->A=k;
					temp->next->Ax=midbx;
					temp->next->Ay=midby;
				}
				temp->next->B=D->B;
				temp->next->Bx=D->Bx;
				temp->next->By=D->By;
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
        long double b,c,d,e,dummy;				
        nAtoms=0;
	ofstream vor;
	vor.open("vor");
	infile>>dummy;
	infile>>dummy;
        while(infile>>b>>c>>d>>e) 
	{ 
            Atoms[nAtoms].x=b;
            Atoms[nAtoms].y=c;
	    Atoms[nAtoms].radius=e;
	    //cout<<nAtoms<<"\t"<<b<<"\t"<<c<<"\n";
            nAtoms++;
	}
	update_neighbours(Atoms,nAtoms);
	double area=0;
	vertice *save;
		//cout<<"here\n";
	for(SAM=0;SAM<nAtoms;SAM++)
	//for(SAM=428;SAM<429;SAM++)
	{
		{	
	////////cout<<SAM<<"=SAM\n";
	        first_delunay(&(Atoms[SAM]),Atoms);
	////////cout<<"SAMhere\n";
	////////cout<<Atoms[SAM].ignore<<"\n";
////////	if(Atoms[SAM].ignore)
////////	{
////////		//cout<<">>\n";
////////		continue;
////////	}
		complete_del(&(Atoms[SAM]),Atoms,nAtoms);
		//cout<<"SAMhere2\n";
		delunay *D;
		D=Atoms[SAM].D.initial;
		int count=0;
	////////while(1)
	////////{
	////////	//print_delunay(&(Atoms[SAM]),D,Atoms);
        ////////        if(SAM==443 || SAM==441 || SAM==402 || SAM==418)
        ////////        {
	////////	print_delunay(&(Atoms[SAM]),D,Atoms);

	////////	}
	////////	count++;
	////////	if(D->next)
	////////	{
	////////		D=D->next;
	////////	}
	////////	else
	////////		break;
	////////}
		D=Atoms[SAM].D.initial;
		double area_s=0;
		//cout<<SAM<<"sam\n";
                for(int i=0;i<Atoms[SAM].conti;i++)
                {
			//cout<<"what\n";
                	D=Atoms[SAM].D.initial;
		        //cout<<Atoms[SAM].contigous[i]<<"\n";
		////////cout<<Atoms[Atoms[SAM].contigous[i]].x<<"\t"<<Atoms[Atoms[SAM].contigous[i]].y<<"\n";
			double a,b,p,q,x,y;
			x=Atoms[SAM].x;
			y=Atoms[SAM].y;
			int flaga=1;
			int flagb=1;
			delunay *D_ONE=NULL;
			delunay *D_TWO=NULL;
                	while(1)
                	{
                		if(D->A==i)
                		{
                			//vor<<D->circum_x<<"\t"<<D->circum_y<<"\t";
					if(!D_ONE)
					{
						D_ONE=D;
					}
					else if(!D_TWO)
						D_TWO=D;
                		}
                		if(D->B==i)
                		{
                			//vor<<D->circum_x<<"\t"<<D->circum_y<<"\t";
					if(!D_ONE)
					{
						D_ONE=D;
					}
					else if(!D_TWO)
						D_TWO=D;
                		}
                		if(D->next)
                			D=D->next;
                		else 
                		{
                			break;
                		}
                		
                	}
	        	//cout<<D_ONE->circum_x<<"\t"<<D_ONE->circum_y<<"\t";
	        	vertice *temp_vert_o;
	        	vertice *temp_vert_d;
	        	if(!start)
	        	{
	        		start=new vertice;
	        		start->p=new site;
	        		start->p->x=D_ONE->circum_x;
	        		start->p->y=D_ONE->circum_y;
				start->A=SAM;
				start->D=D_ONE;
	        		temp_vert_o=start;
	        	}
	        	else 
	        	{
	        		vertice *temp;
	        		temp=new vertice;
	        		temp->p=new site;
	        		temp->p->x=D_ONE->circum_x;
                                temp->p->y=D_ONE->circum_y;
				temp->A=SAM;
				temp->D=D_ONE;
	        		temp=V->insert_vertice(start,temp);
	        		//display_SITE(temp->p);
	        		temp_vert_o=temp;
	        	}
	        	//cout<<D_TWO->circum_x<<"\t"<<D_TWO->circum_y<<"\n";
	        	vertice *temp;
	        	temp=new vertice;
	        	temp->p=new site;
	        	temp->p->x=D_TWO->circum_x;
                        temp->p->y=D_TWO->circum_y;
			temp->A=SAM;
			temp->D=D_TWO;
	        	temp=V->insert_vertice(start,temp);
	                //display_SITE(tempert);
	                //display_SITE(temp_vert_d->p);
			double m;
			double X,Y,dis;
			X=Atoms[Atoms[SAM].contigous[i]].x-Atoms[SAM].x;
			Y=Atoms[Atoms[SAM].contigous[i]].y-Atoms[SAM].y;
			X=(X-(twob*lround(X/twob)));
			Y=(Y-(twob*lround(Y/twob)));
			dis=sqrt(distance(X,Y));
			m=Y/X;
			int sign_C;
			if((D_ONE->circum_y-Atoms[SAM].y)-(m*(D_ONE->circum_x-Atoms[SAM].x))<0.)
				sign_C=-1;
			else
				sign_C=1;
			int sign_N;
			if((D_TWO->circum_y-Atoms[SAM].y)-(m*(D_TWO->circum_x-Atoms[SAM].x))<0.)
				sign_N=-1;
			else
				sign_N=1;
			//if(Atoms[SAM].bondinvoid[i])
			if(dis>Atoms[Atoms[SAM].contigous[i]].radius+Atoms[SAM].radius)
				Atoms[SAM].bondinvoid[i]=1;
			else if(sign_N == sign_C)
				Atoms[SAM].bondinvoid[i]=1;
			else 
				Atoms[SAM].bondinvoid[i]=0;
			//if(Atoms[SAM].bondinvoid[i])
			{
				vor<<x<<"\t"<<y<<"\t";
                		vor<<D_ONE->circum_x<<"\t"<<D_ONE->circum_y<<"\t"<<D_TWO->circum_x<<"\t"<<D_TWO->circum_y<<"\n";
                		vor<<"\n";
			}
	        	//display_SITE(temp->p);
	        	temp_vert_d=temp;
			add_connected(temp_vert_o,temp_vert_d,Atoms[SAM].bondinvoid[i]);
			add_connected(temp_vert_d,temp_vert_o,Atoms[SAM].bondinvoid[i]);
////////                if(SAM==428 || SAM==441 || SAM==402 || SAM==418)
////////                {
////////			D=Atoms[SAM].D.initial;
////////			print_delunay(&(Atoms[SAM]),D,Atoms);
////////	      //	cout<<"####\n";
////////              //	cout<<Atoms[SAM].x<<"\t"<<Atoms[SAM].y<<"\n";
////////              //	cout<<Atoms[SAM].contigous[i]<<"\t"<<Atoms[SAM].bondinvoid[i]<<"\n";
////////              //	cout<<Atoms[Atoms[SAM].contigous[i]].x<<"\t"<<Atoms[Atoms[SAM].contigous[i]].y<<"\n";
////////              //	display_SITE(temp_vert_o->p);
////////              //	display_SITE(temp_vert_d->p);
////////               }

                }
		
		}
	}
	//return 0;

	vertice *temp_start;
	temp_start=start;
	double r_cut_sq=0.5*0.5;
	int void_vert_count=0;
	while(1)
	{
            if(temp_start->v_neigh_count == 2)
            {
                V->delete_vertice(start,temp_start);
	        vertice *temp;
	        temp=temp_start;
                temp_start=temp_start->next;
	        delete temp->p;
	        delete temp;
                continue;
            }
            int flag=1;
            double AX,AY,BX,BY,X,Y;
            BX=temp_start->p->x ;
            BY=temp_start->p->y;
            AX=Atoms[temp_start->A].x;
            AY=Atoms[temp_start->A].y;
            double dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
            if(dis<Atoms[temp_start->A].radius)
            {
                flag=0;
                //break;
            }
            AX=Atoms[Atoms[temp_start->A].contigous[temp_start->D->A]].x;
            AY=Atoms[Atoms[temp_start->A].contigous[temp_start->D->A]].y;
            dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
            if(dis<Atoms[Atoms[temp_start->A].contigous[temp_start->D->A]].radius)
            {
                flag=0;
                //break;
            }
            AX=Atoms[Atoms[temp_start->A].contigous[temp_start->D->B]].x;
            AY=Atoms[Atoms[temp_start->A].contigous[temp_start->D->B]].y;
            dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
            //cout<<dis<<"\t"<<AX<<"\t"<<AY<<"\t"<<temp_start->D->B<<"\n";;
            if(dis<Atoms[Atoms[temp_start->A].contigous[temp_start->D->B]].radius)
            {
                flag=0;
                //break;
            }
            if(flag )//&& temp_start->v_neigh_count == 3)
            {
                temp_start->is_void=1;
		temp_start->cluster_index=void_vert_count;
                void_vert_count=void_vert_count+1;
            }
            if(temp_start->next)
                temp_start=temp_start->next;
            else
                break;
	}
////////while(1)
////////{
////////	int flag=1;
////////	double AX,AY,BX,BY,X,Y;
////////	//cout<<"ver\t"<<temp_start->p->x<<"\t"<<temp_start->p->y<<"\n";
////////        for(int i=0;i<nAtoms;i++)
////////        {
////////        	//double AX,AY,BX,BY;
////////        	AX=Atoms[i].x;
////////        	AY=Atoms[i].y;
////////        	BX=temp_start->p->x ;
////////        	BY=temp_start->p->y;
////////		X=AX-BX;
////////		Y=AY-BY;
////////		X=(X-(twob*lround(X/twob)));
////////		Y=(Y-(twob*lround(Y/twob)));
////////        	double dis=sqrt((X)*(X)+(Y)*(Y));
////////        	if(dis<r_cut+Atoms[i].radius)
////////        	{
////////        		flag=0;
////////        		break;
////////        	}
////////        	
////////        }	
////////////////BX=temp_start->p->x ;
////////////////BY=temp_start->p->y;
////////////////AX=Atoms[temp_start->A].x;
////////////////AY=Atoms[temp_start->A].y;
////////////////double dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
////////////////if(dis<r_cut+Atoms[temp_start->A].radius)
////////////////{
////////////////	flag=0;
////////////////	//break;
////////////////}
////////////////AX=Atoms[temp_start->D->A].x;
////////////////AY=Atoms[temp_start->D->A].y;
////////////////dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
////////////////if(dis<r_cut+Atoms[temp_start->D->A].radius)
////////////////{
////////////////	flag=0;
////////////////	//break;
////////////////}
////////////////AX=Atoms[temp_start->D->B].x;
////////////////AY=Atoms[temp_start->D->B].y;
////////////////dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
////////////////if(dis<r_cut+Atoms[temp_start->D->B].radius)
////////////////{
////////////////	flag=0;
////////////////	//break;
////////////////}
////////	if(flag )//&& temp_start->v_neigh_count == 3)
////////	{
////////		temp_start->is_void=1;
////////		temp_start->cluster_index=void_vert_count;
////////		void_vert_count=void_vert_count+1;
////////		//cout<<"void\t"<<temp_start->p->x<<"\t"<<temp_start->p->y<<"\n";
////////	}
////////	//else 
////////		//cout<<"not void\t"<<temp_start->p->x<<"\t"<<temp_start->p->y<<"\n";
////////	if(temp_start->next)
////////		temp_start=temp_start->next;
////////	else
////////		break;
////////}
	vertice **cavity_list;
	cavity_list = new (nothrow) vertice*[void_vert_count];
	temp_start=start;
	int i=0;
	while(1)
	{
		if(temp_start->is_void)
		{
			cavity_list[i]=temp_start;
			i++;
		}
		if(temp_start->next)
			temp_start=temp_start->next;
		else
			break;

	}
      //for(int i=246;i<247;i++)
      //{
      //	//cout<<cavity_list[i]<<"\n";
      //	cout<<i<<"\t"<<cavity_list[i]->p->x<<"\t"<<cavity_list[i]->p->y<<"\n";
      //        for(int n=0;n<cavity_list[i]->v_neigh_count;n++)
      //        {
      //        	cout<<i<<"\t"<<cavity_list[i]->neib_vert[n]->p->x<<"\t"<<cavity_list[i]->neib_vert[n]->p->y<<"\t"<<cavity_list[i]->neib_ed[n]<<"\n";
      //        }
      //}
	//return 0;
	int change=1;
	int *old_label;
	old_label=new (nothrow) int[void_vert_count];
	int min;
//	cout<<void_vert_count<<"\n";
	while(change)
	{
			//cout<<cavity_list[0]<<"\n";
        	for(int i=0;i<void_vert_count;i++)
        	{
			//cout<<i<<"\t"<<cavity_list[0]<<"\n";
        		old_label[i]=cavity_list[i]->cluster_index;
        		//cout<<i<<"\t"<<old_label[i]<<"\n";
        	}
			//cout<<cavity_list[0]<<"\n";
		for(int i=0;i<void_vert_count;i++)
		{
			//cout<<i<<"\n";
			min=cavity_list[i]->cluster_index;
////////		cout<<i<<"\n";
////////		cout<<cavity_list[i]->v_neigh_count<<"\n";
			for(int n=0;n<cavity_list[i]->v_neigh_count;n++)
			{
				//cout<<"n="<<n<<"\n";
				//if(cavity_list[i]->neib_vert[n]->cluster_index==751)
				//cout<<i<<"\t"<<n<<"\t"<<cavity_list[i]->neib_ed[n]<<"\n";
				if(cavity_list[i]->neib_vert[n]->cluster_index != -1 && cavity_list[i]->neib_ed[n])
				{
					if(min>cavity_list[i]->neib_vert[n]->cluster_index)
						min=cavity_list[i]->neib_vert[n]->cluster_index;
				}
				else
				{
					//cout<<"these guys\n";
				}

			}
			for(int n=0;n<cavity_list[i]->v_neigh_count;n++)
			{
				if(cavity_list[i]->neib_vert[n]->cluster_index != -1 && cavity_list[i]->neib_ed[n])
					cavity_list[i]->neib_vert[n]->cluster_index=min;
			}

		}
		//cout<<"here1\n";
		change=0;
		int flag=0;
		for(int i=0; i<void_vert_count; i++) 
		{
		    flag=(old_label[i]!=cavity_list[i]->cluster_index);
		    if(flag) 
		    {
			change=1;
			break;
		    }
		}

	}
	ofstream cav;
	cav.open("cav");
        for(int i=0;i<void_vert_count;i++)
	{
		cav<<"#"<<i<<"\n\n";
		for(int j=0;j<void_vert_count;j++)
			if(cavity_list[j]->cluster_index==i)
			{
				//cout<<cavity_list[j]->A<<"\t"<<cavity_list[j]->D->A<<"\t"<<cavity_list[j]->D->B<<"\n";
		        	cav<<cavity_list[j]->p->x<<"\t"<<cavity_list[j]->p->y<<"\n";
			}

	}
	//V->display_conn(start);
	//CALCULATING THE VOID VOLUME IN A GIVEN CAVITY
	double *cav_area;
	cav_area= new (nothrow) double[void_vert_count]; 
	double *cav_lenght;
	cav_lenght= new (nothrow) double[void_vert_count]; 
	for(int i=0;i<void_vert_count;i++)
	//for(int i=11;i<12;i++)
	{
		cav_area[i]=0;
		cav_lenght[i]=0;
		for(int j=0;j<void_vert_count;j++)
		{
			if(cavity_list[j]->cluster_index==i)
			{
				//print_delunay(&Atoms[cavity_list[j]->A],cavity_list[j]->D,Atoms);
				double A1x,A1y,A2x,A2y,A3x,A3y,Vx,Vy;
				double E12x,E12y;
				double E13x,E13y;
				double E23x,E23y;
				double A1r,A2r,A3r;
				double m;
				int sign1,sign2;
				int S12,S23,S13;
				double C;
				A1x=Atoms[cavity_list[j]->A].x;
				A1y=Atoms[cavity_list[j]->A].y;
				A1r=Atoms[cavity_list[j]->A].radius;
				A2x=Atoms[Atoms[cavity_list[j]->A].contigous[cavity_list[j]->D->A]].x-A1x;
				A2y=Atoms[Atoms[cavity_list[j]->A].contigous[cavity_list[j]->D->A]].y-A1y;
				A2r=Atoms[Atoms[cavity_list[j]->A].contigous[cavity_list[j]->D->A]].radius;
				A3x=Atoms[Atoms[cavity_list[j]->A].contigous[cavity_list[j]->D->B]].x-A1x;
				A3y=Atoms[Atoms[cavity_list[j]->A].contigous[cavity_list[j]->D->B]].y-A1y;
				A3r=Atoms[Atoms[cavity_list[j]->A].contigous[cavity_list[j]->D->B]].radius;
				E12x=cavity_list[j]->D->Ax-A1x;
				E12y=cavity_list[j]->D->Ay-A1y;
				E13x=cavity_list[j]->D->Bx-A1x;
				E13y=cavity_list[j]->D->By-A1y;
				E12x=(E12x-(twob*lround(E12x/twob)));
				E12y=(E12y-(twob*lround(E12y/twob)));
				E13x=(E13x-(twob*lround(E13x/twob)));
				E13y=(E13y-(twob*lround(E13y/twob)));
				Vx=cavity_list[j]->p->x-A1x;
				Vy=cavity_list[j]->p->y-A1y;
				A2x=(A2x-(twob*lround(A2x/twob)));
				A2y=(A2y-(twob*lround(A2y/twob)));
				A3x=(A3x-(twob*lround(A3x/twob)));
				A3y=(A3y-(twob*lround(A3y/twob)));
				Vx=(Vx-(twob*lround(Vx/twob)));
				Vy=(Vy-(twob*lround(Vy/twob)));
				double X,Y,x,y;
				double rA,rS,DIS,l,dis_i,tan_sq;
				X=Atoms[Atoms[cavity_list[j]->A].contigous[cavity_list[j]->D->A]].x-Atoms[Atoms[cavity_list[j]->A].contigous[cavity_list[j]->D->B]].x;
				Y=Atoms[Atoms[cavity_list[j]->A].contigous[cavity_list[j]->D->A]].y-Atoms[Atoms[cavity_list[j]->A].contigous[cavity_list[j]->D->B]].y;
				X=(X-(twob*lround(X/twob)));
				Y=(Y-(twob*lround(Y/twob)));
				DIS=sqrt(X*X+Y*Y);
				rA=A2r;
				rS=A3r;
				l=0.5*(DIS+(rS*rS-rA*rA)/DIS);
				x=l/DIS*X;
				y=l/DIS*Y;
				x=x+Atoms[Atoms[cavity_list[j]->A].contigous[cavity_list[j]->D->B]].x;
				y=y+Atoms[Atoms[cavity_list[j]->A].contigous[cavity_list[j]->D->B]].y;
				E23x=x-A1x;
				E23y=y-A1y;
				E23x=(E23x-(twob*lround(E23x/twob)));
				E23y=(E23y-(twob*lround(E23y/twob)));
				//calculating the area of A1VA2
				m=A2y/A2x;	
				if((A3y-m*A3x) > 0. )
				{
					sign1=1;
				}
				else
					sign1=-1;
				if((Vy-m*Vx) > 0. )
				{
					sign2=1;
				}
				else
					sign2=-1;
				if(sign1!=sign2)
				{
					S12=-1;
				}
				else S12=1;
				m=A3y/A3x;	
				if((A2y-m*A2x) > 0. )
				{
					sign1=1;
				}
				else
					sign1=-1;
				if((Vy-m*Vx) > 0. )
				{
					sign2=1;
				}
				else
					sign2=-1;
				if(sign1!=sign2)
				{
					S13=-1;
				}
				else 
					S13=1;
				//m=Y/X;
				m=(A3y-A2y)/(A3x-A2x);
				C=A3y-m*A3x;
				if(-1.*C > 0. )
				{
					sign1=1;
				}
				else
					sign1=-1;
				//cout<<sign1<<"\n";
				if((Vy-m*Vx-C) > 0. )
				{
					sign2=1;
				}
				else
					sign2=-1;
				//cout<<sign2<<"\n";
				if(sign1!=sign2)
				{
					S23=-1;
				}
				else 
					S23=1;
				
			      //cout<<Vx+A1x<<"\t"<<Vy+A1y<<"\n";
			      //cout<<E12x+A1x<<"\t"<<E12y+A1y<<"\n";
			      //cout<<E13x+A1x<<"\t"<<E13y+A1y<<"\n";
			      //cout<<E23x+A1x<<"\t"<<E23y+A1y<<"\n";
			        //cout<<A1x<<"\t"<<A1y<<"\n";
			        //cout<<A2x+A1x<<"\t"<<A2y+A1y<<"\n";
			        //cout<<A3x+A1x<<"\t"<<A3y+A1y<<"\n";
				//cout<<"begin\n";
				cav_area[i]=cav_area[i]+S12*area_trangle(Vx,Vy,E12x,E12y,A2x,A2y,A2r,A1x,A1y);
			        cav_area[i]=cav_area[i]+S12*area_trangle(Vx,Vy,E12x,E12y,0.,0.,A1r,A1x,A1y);
			        cav_area[i]=cav_area[i]+S13*area_trangle(Vx,Vy,E13x,E13y,A3x,A3y,A3r,A1x,A1y);
			        cav_area[i]=cav_area[i]+S13*area_trangle(Vx,Vy,E13x,E13y,0.,0.,A1r,A1x,A1y);
			        cav_area[i]=cav_area[i]+S23*area_trangle(Vx,Vy,E23x,E23y,A3x,A3y,A3r,A1x,A1y);
			        cav_area[i]=cav_area[i]+S23*area_trangle(Vx,Vy,E23x,E23y,A2x,A2y,A2r,A1x,A1y);
				cav_lenght[i]=cav_lenght[i]+S12*perimeter(Vx,Vy,E12x,E12y,A2x,A2y,A2r);
			        cav_lenght[i]=cav_lenght[i]+S12*perimeter(Vx,Vy,E12x,E12y,0.,0.,A1r);
			        cav_lenght[i]=cav_lenght[i]+S13*perimeter(Vx,Vy,E13x,E13y,A3x,A3y,A3r);
			        cav_lenght[i]=cav_lenght[i]+S13*perimeter(Vx,Vy,E13x,E13y,0.,0.,A1r);
			        cav_lenght[i]=cav_lenght[i]+S23*perimeter(Vx,Vy,E23x,E23y,A3x,A3y,A3r);
			        cav_lenght[i]=cav_lenght[i]+S23*perimeter(Vx,Vy,E23x,E23y,A2x,A2y,A2r);

			//0.5*abs((x-p)*(b-y)-(x-a)*(q-y))
			      //  cout<<"S12="<<S12<<"\n";
			      //  cout<<"S13="<<S13<<"\n";
			      //  cout<<"S23="<<S23<<"\n";
			      //cout<<E12x+A1x<<"\t"<<E12y+A1y<<"\n";
			      //cout<<E13x+A1x<<"\t"<<E13y+A1y<<"\n";
			      //cout<<E23x+A1x<<"\t"<<E23y+A1y<<"\n";



			////////cout<<cavity_list[j]->D->Ax<<"\t"<<cavity_list[j]->D->Ay<<"\n";
			////////cout<<cavity_list[j]->D->Bx<<"\t"<<cavity_list[j]->D->By<<"\n";
			////////cout<<x<<"\t"<<y<<"\n";



				//cout<<cavity_list[j]->A<<"\t"<<cavity_list[j]->D->A<<"\t"<<cavity_list[j]->D->B<<"\n";
				//print_delunay(&Atoms[cavity_list[j]->A],cavity_list[j]->D,Atoms);
		        	//cout<<cavity_list[j]->p->x<<"\t"<<cavity_list[j]->p->y<<"\n";
			}

		}
	}
	double cav_tot=0.;
	double ca_per_tot=0.;
        for(int i=0;i<void_vert_count;i++)
        {
		cav_tot=cav_tot+cav_area[i];
		ca_per_tot=ca_per_tot+cav_lenght[i];
        	//cout<<i<<"\t"<<cav_area[i]<<"\t"<<cav_lenght[i]<<"\n";
        }
	cout<<cav_tot<<"\t"<<ca_per_tot<<"\n";
	temp_start=start;
	while(1)
	{
		vertice *temp;
		temp=temp_start;
		if(temp_start->next)
		{
			temp_start=temp_start->next;
			delete temp->p;
			delete temp;
		}
		else
		{
			delete temp->p;
			delete temp;
			break;
		}
		
	}
	delunay *D;
	for(SAM=0;SAM<nAtoms;SAM++)
	{
                D=Atoms[SAM].D.initial;
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
	delete[] cavity_list;
	delete[] old_label;
	return 0;
}
