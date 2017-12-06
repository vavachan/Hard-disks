#include<iostream>
#include<math.h>
#include<fstream>
#include<ostream>
#include<iomanip>
#include <set>
#include <limits>
using namespace std;
#define distance(x,y) x*x+y*y
int dim=3;
double box=0;
double twob;
double density;
double r_cut=0.0;
double DMIN=0.00000001;//std::numeric_limits<double>::min();
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
    void delete_vertice(vertice *,vertice *,int type);
    vertice* insert_vertice(vertice *,vertice *,int,int);
    void display_vertex(vertice *);
    void display_conn(vertice *ve);
}*V;
struct delunay
{
    int A;
    int B;
    int C;
    int a;
    int b;
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
    double r;
    //vert_list *V;
    vertice *neib_vert[10];
    int neib_ed[10];
    int v_neigh_count=0;
}**start,*s_temp,*sites;
struct container_vertice
{
    struct vertice *V;
    struct container_vertice *next=NULL;
    struct container_vertice *prev=NULL;
}**CSTART,*VOID_START;
int compare(struct site *p1,struct site *p2,int del=0)
{
    double DX,DY;
    DX=p1->x-p2->x;
    DY=p1->y-p2->y;
    //cout<<DX<<"\t"<<DY<<"\n";
    if(!del)
    {
	    DX=(DX-(twob*lround(DX/twob)));
	    DY=(DY-(twob*lround(DY/twob)));
    }
    //cout<<DX<<"\t"<<DY<<"\n";
////////if( DY > 0.  )
////////	return 1;
////////else if ((  DY == 0. ) && ( DX > 0. ) )
////////	return 1;
////////else if ( ( DY==0.) && ( DX == 0. ) )
////////	return 0;
////////else
////////	return -1;
    if( abs(DX) <DMIN && abs(DX) >DMIN )
    {
        cout<<p1->x<<"\t"<<p1->y<<"\n";
        cout<<p2->x<<"\t"<<p2->y<<"\n";
    }
    if( abs(DY) <DMIN  && abs(DX) <DMIN )
        return 0;
    else if ( (p1->y == p2->y) && ( p1->x > p2->x ) )
        return 1;
    else if( p1->y > p2->y )
        return 1;
    else
        return -1;
}
void add_connected(vertice *focus,vertice *add,int binv,int debug=0)
{
    int flag=1;
    int i;
    if(debug)
    {
	    display_SITE(focus->p);
	    display_SITE(add->p);
    }
    //if(!compare(focus->p,add->p),0)
	//cout<<"here22323\n";
    for(i=0; i<focus->v_neigh_count; i++)
    {
        //if(focus->neib_vert[i]==add)
	if(focus->neib_vert[i])	
        if(!compare(focus->neib_vert[i]->p,add->p))
        {
	        if(focus->neib_vert[i]->p->y-add->p->y)
	        {
	        	if(focus->neib_vert[i]->p->y<0.)
	        	{
	        		focus->neib_vert[i]=add;
	        	}
	        }

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
class set_of_delunay
{
public :
    delunay *initial=NULL;
};
struct atom3D
{
    double x=0;
    double y=0;
    double z=0;
    int neighlist[500];
    int contigous[500];
    int bondinvoid[500];
    int edge_index[500]= {0};
    int neighbours=0;
    int conti=0;
    struct face *F=NULL;
    double radius=1.;
    set_of_delunay D;
};
struct atom
{
    double x=0;
    double y=0;
    int neighlist[500];
    int *contigous[500];
    int *bondinvoid[500];
    int *edge_index[500];
    int neighbours=0;
    int *conti;
    struct face *F=NULL;
    double radius=1.;
    double ignore=0;
    container_vertice **Cstart;
    set_of_delunay *D;
    int save_neighlist[500];
    int *save_contigous[500];
    int *save_bondinvoid[500];
    int *save_edge_index[500];
    int save_neighbours=0;
    int *save_conti;
    int type;
    set_of_delunay *save_D;
};
void save_atom(atom *Atom,int TYPE)
{
    for(int i=0; i<Atom->neighbours; i++)
    {
        Atom->save_neighlist[i]=Atom->neighlist[i];
    }
    for(int i=0; i<Atom->conti[TYPE]; i++)
    {
        Atom->save_contigous[i][TYPE]=Atom->contigous[i][TYPE];
        Atom->save_bondinvoid[i][TYPE]=Atom->bondinvoid[i][TYPE];
        Atom->save_edge_index[i][TYPE]=Atom->edge_index[i][TYPE];
    }
    Atom->save_neighbours=Atom->neighbours;
    Atom->save_conti[TYPE]=Atom->conti[TYPE];

}
void reset_atom(atom *Atom,int TYPE)
{
    for(int i=0; i<Atom->save_neighbours; i++)
    {
        Atom->neighlist[i]=Atom->save_neighlist[i];
    }
    for(int i=0; i<Atom->save_conti[TYPE]; i++)
    {
        Atom->contigous[i][TYPE]=Atom->save_contigous[i][TYPE];
        Atom->bondinvoid[i][TYPE]=Atom->save_bondinvoid[i][TYPE];
        Atom->edge_index[i][TYPE]=Atom->save_edge_index[i][TYPE];
    }
    Atom->neighbours=Atom->save_neighbours;
    Atom->conti[TYPE]=Atom->save_conti[TYPE];
    Atom->D[TYPE]=Atom->save_D[TYPE];
}
int insert_cvertice(container_vertice *EV,container_vertice *v,container_vertice *&Cstart)
{
    int flag;
    flag=compare(EV->V->p,v->V->p);
    //cout<<"insert\n";
    //display_SITE(v->V->p);
    //cout<<"vetfoc\n";
    //display_SITE(EV->V->p);
    if(flag==1)
    {
        if(EV->next)
        {
            //		cout<<"everytime?\n";
            insert_cvertice(EV->next,v,Cstart);
        }
        else
        {
            EV->next=v;
            v->prev=EV;
            v->next=NULL;
            ///		cout<<"here?\n";
            return flag;
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
        {
            Cstart=v;
        }
        EV->prev=v;
        //	cout<<"here2?\n";
        return flag;
    }
    else if (flag==0)
    {
        //cout<<"vertice already exists\n";
        delete v;
        return flag;
    }
}
vertice* vert_list::insert_vertice(vertice *EV,vertice *v,int type,int debug=0)
{
    int flag;
    //cout<<"here\n";
    if(debug)
    {
        cout<<"insert\n";
        display_SITE(v->p);
        cout<<"vetfoc\n"<<flush;
        display_SITE(EV->p);
    }
    flag=compare(EV->p,v->p);
    if(debug)
    {
        cout<<"flag\t"<<flag<<"\n";
    }
    if(flag==1)
    {
        if(EV->next)
        {
            insert_vertice(EV->next,v,type,debug);
        }
        else
        {
            EV->next=v;
            v->prev=EV;
            v->next=NULL;
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
            start[type]=v;
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
void vert_list::delete_vertice(vertice *EV,vertice *v,int type)
{
    int flag;
    flag=compare(EV->p,v->p,1);
    ///cout<<flag<<"\n";
    //display_SITE(EV->p);
    //display_SITE(v->p);
    //cout<<"\n";
    if (flag==0)
    {
        if(EV->prev)
            EV->prev->next=EV->next;
        else
            start[type]=EV->next;
        if(EV->next)
            EV->next->prev=EV->prev;
        else
            EV->prev->next=NULL;
        //cout<<"herein delete\n";
//	cout<<v->v_neigh_count<<"\n"<<std::flush;
        for(int n=0; n<v->v_neigh_count; n++)
        {
	  //  cout<<n<<"\n";
            if(v->neib_vert[n])
	    {
        //display_SITE(v->p);
	//	    display_SITE(v->neib_vert[n]->p);
                for(int m=0; m<v->neib_vert[n]->v_neigh_count; m++)
                {
                    if(v->neib_vert[n]->neib_vert[m])
		    {
                        if(!(compare(v->neib_vert[n]->neib_vert[m]->p,v->p)))
                        {
                            v->neib_vert[n]->neib_vert[m]=NULL;
                        }
		    }
                }
	    }
        }
        //delete v->p;
        //delete v;
    }
    else
    {
        if(EV->next)
            delete_vertice(EV->next,v,type);
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
        for(int i=0; i<ve->v_neigh_count; i++)
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
void insert_site(vertice *EV,vertice *v,int debug=0)
{
    int flag;
    //cout<<"here\n";
    if(debug)
    {
        cout<<"insert\n";
        display_SITE(v->p);
        cout<<"vetfoc\n"<<flush;
        display_SITE(EV->p);
    }
    flag=compare(EV->p,v->p);
    if(debug)
    {
        cout<<"flag\t"<<flag<<"\n";
    }
    if(flag==1)
    {
        if(EV->next)
        {
            insert_site(EV->next,v,debug);
        }
        else
        {
            EV->next=v;
            v->prev=EV;
            v->next=NULL;
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
            sites=v;
        EV->prev=v;
        //	cout<<"here2?\n";
    }
    else if (flag==0)
    {
        //cout<<"vertice already exists\n";
        delete v->p;
        delete v;
    }
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
    for(int i=0; i<nAtoms-1; i++)
    {
        for(int j=i+1; j<nAtoms; j++)
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
void print_delunay(atom *ATOM,delunay *D,atom Atoms[],int TYPE)
{
    //cout<<"what\n";
    cout<<"#\n";
    double Sx,Sy;
    double Px,Py;
    Sx=ATOM->x-Atoms[D->a].x;
    Sy=ATOM->y-Atoms[D->a].y;
    Sx=(Sx-(twob*lround(Sx/twob)));
    Sy=(Sy-(twob*lround(Sy/twob)));
    cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\n";
    Sx=ATOM->x-Atoms[D->b].x;
    Sy=ATOM->y-Atoms[D->b].y;
    Sx=(Sx-(twob*lround(Sx/twob)));
    Sy=(Sy-(twob*lround(Sy/twob)));
    cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\n";
    Sx=ATOM->x-Atoms[D->b].x;
    Sy=ATOM->y-Atoms[D->b].y;
    Sx=(Sx-(twob*lround(Sx/twob)));
    Sy=(Sy-(twob*lround(Sy/twob)));
    Px=ATOM->x-Atoms[D->a].x;
    Py=ATOM->y-Atoms[D->a].y;
    Px=(Px-(twob*lround(Px/twob)));
    Py=(Py-(twob*lround(Py/twob)));
    cout<<ATOM->x-Px<<"\t"<<ATOM->y-Py<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\n";
}
void first_delunay(atom *ATOM,atom Atoms[],int TYPE)
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
    ATOM->conti[TYPE]=0;
    for(int i=0; i<ATOM->neighbours; i++)
    {
        double X,Y,x,y;
        double rA,rS,DIS,l,dis_i,tan_sq;
        X=Atoms[ATOM->neighlist[i]].x-ATOM->x;
        Y=Atoms[ATOM->neighlist[i]].y-ATOM->y;
        X=(X-(twob*lround(X/twob)));
        Y=(Y-(twob*lround(Y/twob)));
        DIS=sqrt(X*X+Y*Y);
        rA=Atoms[ATOM->neighlist[i]].radius+r_cut;
        rS=ATOM->radius+r_cut;
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
            if(DIS>rS+rA)
            {
                binv1=1;
            }
        }

    }
    ATOM->bondinvoid[ATOM->conti[TYPE]][TYPE]=binv1;
    ATOM->contigous[ATOM->conti[TYPE]][TYPE]=nearest;
    ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
    ATOM->conti[TYPE]++;
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
    for(int i=0; i<ATOM->neighbours; i++)
    {
        if(ATOM->neighlist[i]!=nearest)
        {
            atom L=*ATOM;
            atom M=Atoms[ATOM->contigous[0][TYPE]];
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
            rA=M.radius+r_cut;
            rB=R.radius+r_cut;
            rS=L.radius+r_cut;
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
    ATOM->bondinvoid[ATOM->conti[TYPE]][TYPE]=binv1;
    ATOM->contigous[ATOM->conti[TYPE]][TYPE]=DIS_atom;
    ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
    ATOM->conti[TYPE]++;
    ATOM->D[TYPE].initial= new delunay;
    ATOM->D[TYPE].initial->A=0;
    ATOM->D[TYPE].initial->B=1;
    ATOM->D[TYPE].initial->a=ATOM->contigous[0][TYPE];
    ATOM->D[TYPE].initial->b=ATOM->contigous[1][TYPE];
    ATOM->D[TYPE].initial->circum_x=circx;
    ATOM->D[TYPE].initial->circum_y=circy;
    ATOM->D[TYPE].initial->Ax=midax;
    ATOM->D[TYPE].initial->Ay=miday;
    ATOM->D[TYPE].initial->Bx=midbx;
    ATOM->D[TYPE].initial->By=midby;
    delunay *D;
    //delunay *D;
    //if(TYPE==1)
    //{
    //		cout<<"first\n";
    //	D=ATOM->D[0].initial;
    //	while(1)
    //	{
    //		print_delunay(ATOM,D,Atoms,0);
    //		if(D->next)
    //		{
    //			D=D->next;
    //		}
    //		else
    //			break;
    //	}
    //}
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
void complete_del(atom *ATOM,atom Atoms[],int nAtoms,int TYPE)
{
    double Y_MIN=box*box;
    int DIS_atom;
    int binv;
    int binv2;
    double lmin;
    //loop over all contiguos atom so far
    for(int i=0; i<ATOM->conti[TYPE]; i++)
    {
        Y_MIN=box*box;
        delunay *D;
        //delunay *D;
        //if(TYPE==1)
        //{
        //	cout<<i<<"\n";
        //	D=ATOM->D[0].initial;
        //	while(1)
        //	{
        //		print_delunay(ATOM,D,Atoms,0);
        //		if(D->next)
        //		{
        //			D=D->next;
        //		}
        //		else
        //			break;
        //	}
        //}
        D=ATOM->D[TYPE].initial;
        //make sure that atom is not taking parting two delunay triangles
        if(ATOM->edge_index[i][TYPE]!=2)
        {
            //finding the delunay triangle the atom contiguous atom takes part in
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
                //the next section is basically a copy with D->B
                //cout<<"hereA\n";
                //cout<<i<<"\n";
                double a=ATOM->x;
                double b=ATOM->y;
                double c=Atoms[ATOM->contigous[D->A][TYPE]].x;
                double d=Atoms[ATOM->contigous[D->A][TYPE]].y;
                double Sx=Atoms[ATOM->contigous[D->B][TYPE]].x-a;
                double Sy=Atoms[ATOM->contigous[D->B][TYPE]].y-b;
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
                for(int j=0; j<ATOM->neighbours; j++)
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
                    //to avoid atoms with have two delunay
                    for(k=0; k<ATOM->conti[TYPE]; k++)
                    {
                        if(ATOM->contigous[k][TYPE]==ATOM->neighlist[j])
                        {
                            if(ATOM->edge_index[k][TYPE]==2)
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
                        atom M=Atoms[ATOM->contigous[D->A][TYPE]];
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
                        rA=M.radius+r_cut;
                        rB=R.radius+r_cut;
                        rS=L.radius+r_cut;
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
                            if(DISB>rS+rB)
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
                for(k=0; k<ATOM->conti[TYPE]; k++)
                {
                    if(ATOM->contigous[k][TYPE]==DIS_atom)
                    {
                        flag=0;
                        break;
                    }
                }

                if(flag)
                {
                    ATOM->bondinvoid[ATOM->conti[TYPE]][TYPE]=binv;
                    ATOM->contigous[ATOM->conti[TYPE]][TYPE]=DIS_atom;
                    ATOM->edge_index[D->A][TYPE]++;
                    ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
                }
                else
                {
                    ATOM->edge_index[D->A][TYPE]++;
                    ATOM->edge_index[k][TYPE]++;
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
                    temp->next->A=ATOM->conti[TYPE];
                    temp->next->a=ATOM->contigous[ATOM->conti[TYPE]][TYPE];
                    temp->next->Ax=midbx;
                    temp->next->Ay=midby;
                }
                else
                {
                    temp->next->A=k;
                    temp->next->a=ATOM->contigous[k][TYPE];
                    temp->next->Ax=midbx;
                    temp->next->Ay=midby;
                }
                temp->next->B=D->A;
                temp->next->b=ATOM->contigous[D->A][TYPE];
                temp->next->Bx=D->Ax;
                temp->next->By=D->Ay;
                temp->next->circum_x=circx;
                temp->next->circum_y=circy;
                if(flag)
                    ATOM->conti[TYPE]++;
            }// D->A loop
            else
            {
                //cout<<"hereb\n";
                ///cout<<i<<"\n";
                double a=ATOM->x;
                double b=ATOM->y;
                double c=Atoms[ATOM->contigous[D->B][TYPE]].x;
                double d=Atoms[ATOM->contigous[D->B][TYPE]].y;
                double Sx=Atoms[ATOM->contigous[D->A][TYPE]].x-a;
                double Sy=Atoms[ATOM->contigous[D->A][TYPE]].y-b;
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
                for(int j=0; j<ATOM->neighbours; j++)
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
                    for(k=0; k<ATOM->conti[TYPE]; k++)
                    {
                        if(ATOM->contigous[k][TYPE]==ATOM->neighlist[j])
                        {
                            if(ATOM->edge_index[k][TYPE]==2)
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
                        atom M=Atoms[ATOM->contigous[D->B][TYPE]];
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
                        rA=M.radius+r_cut;
                        rB=R.radius+r_cut;
                        rS=L.radius+r_cut;
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
                            if(DISB>rS+rB)
                            {
                                binv=1;
                            }
                        }
                    }
                }//j loop
                flag=1;
                if(lmin<0.)
                {
                    ATOM->ignore=1;
                    //return;
                }
                if(lmin<0.)
                {

                }
                for(k=0; k<ATOM->conti[TYPE]; k++)
                {
                    if(ATOM->contigous[k][TYPE]==DIS_atom)
                    {
                        flag=0;
                        break;
                    }
                }

                if(flag)
                {
                    ATOM->bondinvoid[ATOM->conti[TYPE]][TYPE]=binv;
                    ATOM->contigous[ATOM->conti[TYPE]][TYPE]=DIS_atom;
                    ATOM->edge_index[D->B][TYPE]++;
                    ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
                }
                else
                {
                    ATOM->edge_index[D->B][TYPE]++;
                    ATOM->edge_index[k][TYPE]++;
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
                    temp->next->A=ATOM->conti[TYPE];
                    temp->next->a=ATOM->contigous[ATOM->conti[TYPE]][TYPE];
                    temp->next->Ax=midbx;
                    temp->next->Ay=midby;
                }
                else
                {
                    temp->next->A=k;
                    temp->next->a=ATOM->contigous[k][TYPE];
                    temp->next->Ax=midbx;
                    temp->next->Ay=midby;
                }
                temp->next->B=D->B;
                temp->next->b=ATOM->contigous[D->B][TYPE];
                temp->next->Bx=D->Bx;
                temp->next->By=D->By;
                temp->next->circum_x=circx;
                temp->next->circum_y=circy;
                if(flag)
                    ATOM->conti[TYPE]++;
            }// D->B loop
        }//For atom i in conti
    }//loop  over conti
}//end

int main()
{
    int nAtoms=0;
    atom *Atoms;
    int counter=0;
    int config_count=0;
    //cout << "Maximum value for dou: " << std::numeric_limits<double>::min() << '\n';
    std::ifstream infile("config_256_0.38_2_0.70.dat");//dat_save");//config_2000_0.38_2_0.70.dat");
    infile>>nAtoms;
    infile>>box;
    cout<<nAtoms<<"\t"<<box<<"\n";
    config_count=3;
    int ntypes=2;
    twob=2*box;
    density=nAtoms/(twob*twob);
    int SAM=0;
    double *radius=new (nothrow) double[ntypes];
    long double b,c,d,e;
    for(int t=0; t<ntypes; t++)
    {
        infile>>radius[t];
    }
    //nAtoms=0;
    char buffer[64];
    vertice *temp_site;
    temp_site=sites;
    double free_dist[1000]={0.};
    double max_free_area=0.;
    double width;
    snprintf(buffer,sizeof(char)*64,"free_dist");//_%d_%f.dat",int(nAtoms),Press);
    ofstream fdist;
    fdist.open(buffer);
    for(int nconfig=0;nconfig<config_count;nconfig++)
    {
	start = new (nothrow) vertice*[ntypes];
	CSTART = new (nothrow) container_vertice*[ntypes];
	Atoms = new (nothrow) atom[nAtoms];
	cout<<nconfig<<"\n"<<std::flush;
	counter=0;
	while(infile>>b>>c>>d>>e)
	{
	  //Atoms[nAtoms].x=b;
	  //Atoms[nAtoms].y=c;
	  //Atoms[nAtoms].radius=d;
	    //cout<<counter<<"\t"<<b<<"\t"<<c<<"\t"<<d<<"\n"<<std::flush;
	  //nAtoms++;
	    counter++;
	    //cout<<counter<<"\n";
	    if(sites==NULL)
	    {
		    cout<<"here\n";
	    	sites=new vertice;
	    	sites->p=new site;
	    	sites->p->x=b;
	    	sites->p->y=c;
	    	sites->r=d;
	    }
	    else 
	    {
	    	temp_site=new vertice;
	    	temp_site->p=new site;
	    	temp_site->p->x=b;
	    	temp_site->p->y=c;
	    	temp_site->r=d;
	    	insert_site(sites,temp_site);
	    }
	    if(counter==nAtoms)
	    {
	        break;
	    }
	}
	temp_site=sites;
	int cunt=0;
	while(1)
	{
	    //display_SITE(temp_site->p);
	    //cout<<temp_site->r<<"\n";
	    Atoms[cunt].x=temp_site->p->x;
	    Atoms[cunt].y=temp_site->p->y;
	    Atoms[cunt].radius=temp_site->r;
	    cunt++;
	    if(temp_site->next)
	    	temp_site=temp_site->next;
	    else
	    	break;
	}
     // for(int i=0;i<nAtoms;i++)
     // {
     //         cout<<Atoms[i].x<<"\t"<<Atoms[i].y<<"\t"<<Atoms[i].radius<<"\t"<<Atoms[i].radius<<"\n";
     // }
	//return 0;
	update_neighbours(Atoms,nAtoms);
	double area=0;
	vertice *save;
	//cout<<"here\n";
	for(int i=0; i<nAtoms; i++)
	{
	    for(int t=0; t<ntypes; t++)
	    {
	        if(Atoms[i].radius==radius[t])
	        {
	    	Atoms[i].type=t;
	    	break;
	        }
	    }
	    Atoms[i].Cstart= new(nothrow) container_vertice*[ntypes];
	    Atoms[i].D= new(nothrow) set_of_delunay[ntypes];
	    Atoms[i].conti = new (nothrow) int[ntypes];
	    Atoms[i].save_conti = new (nothrow) int[ntypes];
	    Atoms[i].save_D= new(nothrow) set_of_delunay[ntypes];
	    for(int t=0; t<500; t++)
	    {
	        Atoms[i].contigous[t]= new (nothrow) int[ntypes];
	        Atoms[i].edge_index[t]= new (nothrow) int[ntypes];
	        Atoms[i].bondinvoid[t]= new (nothrow) int[ntypes];
	        Atoms[i].save_contigous[t]= new (nothrow) int[ntypes];
	        Atoms[i].save_edge_index[t]= new (nothrow) int[ntypes];
	        Atoms[i].save_bondinvoid[t]= new (nothrow) int[ntypes];
	        for(int TYPE=0; TYPE<ntypes; TYPE++)
	        {
	    	Atoms[i].contigous[t][TYPE]=0;
	    	Atoms[i].edge_index[t][TYPE]=0;
	    	Atoms[i].bondinvoid[t][TYPE]=0;
	        }
	    }
	    for(int TYPE=0; TYPE<ntypes; TYPE++)
	    {
	        Atoms[i].Cstart[TYPE]=NULL;
	        start[TYPE]=NULL;
	        CSTART[TYPE]=NULL;
	    }
	}

	for(int TYPE=0; TYPE<ntypes; TYPE++)
	{
	    snprintf(buffer,sizeof(char)*64,"vor_%d",int(TYPE));//_%d_%f.dat",int(nAtoms),Press);
	    ofstream vor;
	    vor.open(buffer);
	    r_cut=radius[TYPE];
	    cout<<TYPE<<"\t"<<r_cut<<"\n";
	    for(SAM=0 ; SAM<nAtoms ; SAM++)
	    {
	        {
	    	first_delunay(&(Atoms[SAM]),Atoms,TYPE);
	    	complete_del(&(Atoms[SAM]),Atoms,nAtoms,TYPE);
	    	delunay *D;
	    	int count=0;
	    	double area_s=0;
	    	//cout<<SAM<<"\n";
	    	//if(SAM==492)
	    	//{
	    	//	cout<<TYPE<<"\n";
	    	//	D=Atoms[SAM].D[0].initial;
	    	//	cout<<"cuthere\n";
	    	//	while(1)
	    	//	{
	    	//		print_delunay(&(Atoms[SAM]),D,Atoms,0);
	    	//		count++;
	    	//		if(D->next)
	    	//		{
	    	//			D=D->next;
	    	//		}
	    	//		else
	    	//			break;
	    	//	}
	    	//}
	    	for(int i=0; i<Atoms[SAM].conti[TYPE]; i++)
	    	{
	    	    D=Atoms[SAM].D[TYPE].initial;
	    	    double a,b,p,q,x,y;
	    	    x=Atoms[SAM].x;
	    	    y=Atoms[SAM].y;
	    	    int flaga=1;
	    	    int flagb=1;
	    	    delunay *D_ONE=NULL;
	    	    delunay *D_TWO=NULL;
	    	    D=Atoms[SAM].D[TYPE].initial;

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
	    	    if(!start[TYPE])
	    	    {
	    		start[TYPE]=new vertice;
	    		start[TYPE]->p=new site;
	    		start[TYPE]->p->x=D_ONE->circum_x;
	    		start[TYPE]->p->y=D_ONE->circum_y;
	    		start[TYPE]->A=SAM;
	    		start[TYPE]->D=D_ONE;
	    		temp_vert_o=start[TYPE];
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
	    		temp=V->insert_vertice(start[TYPE],temp,TYPE);
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
	    	    temp=V->insert_vertice(start[TYPE],temp,TYPE);
	    	    //display_SITE(tempert);
	    	    //display_SITE(temp_vert_d->p);
	    	    double m;
	    	    double X,Y,dis;
	    	    X=Atoms[Atoms[SAM].contigous[i][TYPE]].x-Atoms[SAM].x;
	    	    Y=Atoms[Atoms[SAM].contigous[i][TYPE]].y-Atoms[SAM].y;
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
	    	    if(dis>Atoms[Atoms[SAM].contigous[i][TYPE]].radius+Atoms[SAM].radius+2.*r_cut)
	    		Atoms[SAM].bondinvoid[i][TYPE]=1;
	    	    else if(sign_N == sign_C)
	    		Atoms[SAM].bondinvoid[i][TYPE]=1;
	    	    else
	    		Atoms[SAM].bondinvoid[i][TYPE]=0;
	    	    //if(SAM!=585  && Atoms[SAM].contigous[i][TYPE] !=585)
	    	    {
	    		vor<<std::setprecision(15)<<x<<"\t"<<y<<"\t";
	    		vor<<std::setprecision(15)<<D_ONE->circum_x<<"\t"<<D_ONE->circum_y<<"\t"<<D_TWO->circum_x<<"\t"<<D_TWO->circum_y<<"\n";
	    		vor<<"\n";
	    	    }
	    	    //display_SITE(temp->p);
	    	    temp_vert_d=temp;
	    	    add_connected(temp_vert_o,temp_vert_d,Atoms[SAM].bondinvoid[i][TYPE]);
	    	    add_connected(temp_vert_d,temp_vert_o,Atoms[SAM].bondinvoid[i][TYPE]);
	    	    container_vertice *temp_cvert;
	    	    if(Atoms[SAM].Cstart[TYPE]==NULL)
	    	    {
	    		Atoms[SAM].Cstart[TYPE]=new container_vertice;
	    		Atoms[SAM].Cstart[TYPE]->V=temp_vert_o;
	    	    }
	    	    else
	    	    {
	    		temp_cvert=new container_vertice;
	    		temp_cvert->V=temp_vert_o;
	    		insert_cvertice(Atoms[SAM].Cstart[TYPE],temp_cvert,Atoms[SAM].Cstart[TYPE]);
	    	    }
	    	    temp_cvert=new container_vertice;
	    	    temp_cvert->V=temp_vert_d;
	    	    //cout<<Atoms[SAM].Cstart[TYPE]<<"\n";
	    	    //display_SITE(Atoms[SAM].Cstart[TYPE]->V->p);
	    	    insert_cvertice(Atoms[SAM].Cstart[TYPE],temp_cvert,Atoms[SAM].Cstart[TYPE]);

	    	}

	        }
	    }
	    vertice *temp_start;
	    temp_start=start[TYPE];
	    int void_vert_count=0;
	    //cout<<"after  first tessellation \t"<<TYPE<<"\n";;
	    while(1)
	    {
	        if(temp_start->v_neigh_count == 2)
	        {
		      //display_SITE(temp_start->p);
		      //cout<<std::flush;
	          V->delete_vertice(start[TYPE],temp_start,TYPE);
	          vertice *temp;
	          temp=temp_start;
	          temp_start=temp_start->next;
	          delete temp->p;
	          delete temp;
	          continue;
	        }
		
	        //if(temp_start->v_neigh_count == 3)
		  //    display_SITE(temp_start->p);
			
	      //if(TYPE==0)
	      //{
	          //display_SITE(temp_start->p);
	      //}
	        int flag=1;
	        double AX,AY,BX,BY,X,Y;
	        BX=temp_start->p->x ;
	        BY=temp_start->p->y;
	        AX=Atoms[temp_start->A].x;
	        AY=Atoms[temp_start->A].y;
	        double dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
	        if(dis<r_cut+Atoms[temp_start->A].radius)
	        {
	    	flag=0;
	    	//break;
	        }
	        AX=Atoms[Atoms[temp_start->A].contigous[temp_start->D->A][TYPE]].x;
	        AY=Atoms[Atoms[temp_start->A].contigous[temp_start->D->A][TYPE]].y;
	        dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
	        if(dis<r_cut+Atoms[Atoms[temp_start->A].contigous[temp_start->D->A][TYPE]].radius)
	        {
	    	flag=0;
	    	//break;
	        }
	        AX=Atoms[Atoms[temp_start->A].contigous[temp_start->D->B][TYPE]].x;
	        AY=Atoms[Atoms[temp_start->A].contigous[temp_start->D->B][TYPE]].y;
	        dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
	        //cout<<dis<<"\t"<<AX<<"\t"<<AY<<"\t"<<temp_start->D->B<<"\n";;
	        if(dis<r_cut+Atoms[Atoms[temp_start->A].contigous[temp_start->D->B][TYPE]].radius)
	        {
	    	flag=0;
	    	//break;
	        }
	        if(flag )//&& temp_start->v_neigh_count == 3)
	        {
	    	temp_start->is_void=1;
	    	void_vert_count=void_vert_count+1;
	    	//cout<<"void\t"<<temp_start->p->x<<"\t"<<temp_start->p->y<<"\n";
	        }
	        //else
	        //cout<<"not void\t"<<temp_start->p->x<<"\t"<<temp_start->p->y<<"\n";
	        if(temp_start->next)
	    		temp_start=temp_start->next;
	        else
	    		break;
	    }
	    //cout<<"first tessellation \t"<<TYPE<<"\n";;
	    vor.close();
	}
	//CODE BEGINS FOR CALCULATING FREE VOLUME
	//return 0;
	cout<<"are we here\n"<<std::flush;
	ofstream vor;
	for(int t=0; t<ntypes; t++)
	{
	    snprintf(buffer,sizeof(char)*64,"dat_%d",int(t));//_%d_%f.dat",int(nAtoms),Press);
	    ofstream config;
	    config.open(buffer);
	    r_cut=radius[t];
	    for(int n=0; n<nAtoms; n++)
	    {
	        config<<n<<"\t"<<Atoms[n].x<<"\t"<<Atoms[n].y<<"\t"<<Atoms[n].radius<<"\t"<<Atoms[n].radius+r_cut<<"\n";
	    }
	    config.close();
	}
	cout<<"are we here2\n"<<std::flush;
	//SAM=2;
	int TYPE=1;
	//cout<<Atoms[35].conti[TYPE]<<"\n";
	//for(int k=0; k<Atoms[35].conti[TYPE]; k++)
	//{
	//    cout<<k<<"\t"<<Atoms[35].contigous[k][TYPE]<<"\n";
	//}
	//LOOP OVER ALL ATOMS , CALCULATE THE FREE VOL FOR EACH ATOM
	double *freearea;
	double *freeperi;
	freearea=new (nothrow) double [nAtoms];
	freeperi=new (nothrow) double [nAtoms];
	for(int i=0; i<nAtoms; i++)
	{
		freearea[i]=0.;
		freeperi[i]=0.;
	}
	for(int i=0; i<nAtoms; i++)
	{

	    container_vertice *new_vert=NULL;
	    int Atom_in_foc=i;
	    //set_of_delunay D;
	    //cout<<"Atom co-ordiantes\n";
	    //cout<<Atoms[i].x<<"\t"<<Atoms[i].y<<"\n";
	    r_cut=Atoms[i].radius;
	    //cout<<i<<" Atoms NO \n";
	    TYPE=Atoms[i].type;
	    //FIND THE ATOM RADIUS
	    
	//  for(int t=0; t<ntypes; t++)
	//  {
	//      if(r_cut==radius[t])
	//      {
	//          TYPE=t;
	//          break;
	//      }
	//  }
	    //cout<<TYPE<<"\t"<<"TYPE"<<"\n";
	    container_vertice *ctemp;
	    ctemp=Atoms[i].Cstart[TYPE];
	    //REMOVE THE VERTICES WHICH BELONGED TO THE VORONOI CELL OF THE ATOM IN CONSIDERATIO (i)
	  //while(1)
	  //{
	  //    display_SITE(ctemp->V->p);
	  //    //V->delete_vertice(start[TYPE],ctemp->V,TYPE);
	  //    if(ctemp->next)
	  //	ctemp=ctemp->next;
	  //    else
	  //	break;
	  //}
	    ctemp=Atoms[i].Cstart[TYPE];
	    while(1)
	    {
	        //display_SITE(ctemp->V->p);
	        V->delete_vertice(start[TYPE],ctemp->V,TYPE);
	        if(ctemp->next)
	    	ctemp=ctemp->next;
	        else
	    	break;
	    }
	    snprintf(buffer,sizeof(char)*64,"vor_%d",int(TYPE));//_%d_%f.dat",int(nAtoms),Press);
	    vor.open(buffer,std::ios_base::app);
	    CSTART[TYPE]=NULL;
	    //cout<<"nearest\n";
	    //cout<<Atoms[Atoms[i].contigous[0][TYPE]].x<<"\t"<<Atoms[Atoms[i].contigous[0][TYPE]].y<<"\n";
	    //cout<<TYPE<<" TYPE\n";
	    //LOOP OVER ALL THE ATOMS THAT ARE CONTIGUOUS TO i
	    for(int j=0; j<Atoms[i].conti[TYPE]; j++)
	    {
	       //cout<<j<<"\t"<<Atoms[Atoms[i].contigous[j][TYPE]].x<<"\t"<<Atoms[Atoms[i].contigous[j][TYPE]].y<<"\n";
	        //cout<<"eye "<<j<<"\t"<<Atoms[35].conti[TYPE]<<"\n";
	        int flag=0;
	        //FIND THE ATOM INDEX OF THE CONTIGUOUS ATOM
	        int SAM=Atoms[i].contigous[j][TYPE];
	        container_vertice *ctemp;
	        //SAVE THE DETAILS OF THE ATOM 'SAM' BECAUSE THEY ARE GOING TO BE RETESSELLATED
	        ctemp=Atoms[SAM].Cstart[TYPE];
	        save_atom(&(Atoms[SAM]),TYPE);
	        delunay *D;
	        D=Atoms[SAM].D[TYPE].initial;
	        delunay *temp;
	        temp=Atoms[SAM].D[TYPE].initial;
	        Atoms[SAM].save_D[TYPE].initial=temp;
	        //cout<<"nearest\t"<<Atoms[SAM].contigous[0][TYPE]<<"\n";
	        //REMOVE THE ATOM i FROM THE NEIBHOUR LIST OF 'SAM'
	        for(int k=0; k<Atoms[Atoms[i].contigous[j][TYPE]].neighbours-1; k++)
	        {
	    	if(Atoms[Atoms[i].contigous[j][TYPE]].neighlist[k]==i)
	    	{
	    	    flag=1;
	    	}
	    	if(flag)
	    	{
	    	    Atoms[Atoms[i].contigous[j][TYPE]].neighlist[k]=Atoms[Atoms[i].contigous[j][TYPE]].neighlist[k+1];
	    	}
	        }
	        //cout<<"here\n";
	        Atoms[SAM].neighbours=Atoms[SAM].neighbours-1;
	        //FIND ALL THE DELUNAY TRIANGLES THIS 'SAM' TAKES PART IN AFTER REMOVING i
	        first_delunay(&(Atoms[SAM]),Atoms,TYPE);
	        complete_del(&(Atoms[SAM]),Atoms,nAtoms,TYPE);
	        int count=0;
	        double area_s=0;
	        //cout<<SAM<<"\t"<<Atoms[SAM].conti[TYPE]<<"=SAM\n";
	        for(int i=0; i<Atoms[SAM].conti[TYPE]; i++)
	        {
	    	//cout<<i<<"\n";
	    	D=Atoms[SAM].D[TYPE].initial;
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
	    		if(!D_ONE)
	    		{
	    		    D_ONE=D;
	    		}
	    		else if(!D_TWO)
	    		    D_TWO=D;
	    	    }
	    	    if(D->B==i)
	    	    {
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
	    	vertice *temp_vert_o;
	    	vertice *temp_vert_d;
	    	if(!start[TYPE])
	    	{
	    	    start[TYPE]=new vertice;
	    	    start[TYPE]->p=new site;
	    	    start[TYPE]->p->x=D_ONE->circum_x;
	    	    start[TYPE]->p->y=D_ONE->circum_y;
	    	    start[TYPE]->A=SAM;
	    	    start[TYPE]->D=D_ONE;
	    	    temp_vert_o=start[TYPE];
	    	}
	    	else
	    	{
	    	    vertice *temp;
	    	    vertice *temp_o;
	    	    temp=new vertice;
	    	    temp->p=new site;
	    	    temp->p->x=D_ONE->circum_x;
	    	    temp->p->y=D_ONE->circum_y;
	    	    temp->A=SAM;
	    	    temp->D=D_ONE;
	    	    temp_o=temp;
	    	    temp=V->insert_vertice(start[TYPE],temp,TYPE);
	    	    if(temp_o==temp)
	    	    {
	    		double Ax,Ay,Bx,By,Cx,Cy,px,py,m,Co;
	    		int sign_1,sign_2;
	    		int flag=1;
	    		Ax=Atoms[temp->A].x;
	    		Ay=Atoms[temp->A].y;
	    		Bx=Atoms[temp->D->a].x;
	    		By=Atoms[temp->D->a].y;
	    		Cx=Atoms[temp->D->b].x;
	    		Cy=Atoms[temp->D->b].y;
	    		Bx=Bx-Ax;
	    		By=By-Ay;
	    		Cx=Cx-Ax;
	    		Cy=Cy-Ay;
	    		px=Atoms[Atom_in_foc].x;
	    		py=Atoms[Atom_in_foc].y;
	    		px=px-Ax;
	    		py=py-Ay;
	    		Cx=(Cx-(twob*lround(Cx/twob)));
	    		Cy=(Cy-(twob*lround(Cy/twob)));
	    		Bx=(Bx-(twob*lround(Bx/twob)));
	    		By=(By-(twob*lround(By/twob)));
	    		px=(px-(twob*lround(px/twob)));
	    		py=(py-(twob*lround(py/twob)));
	    		m=By/Bx;
	    		if(Cy-m*Cx>0.)
	    		{
	    		    sign_1=1;
	    		}
	    		else
	    		{
	    		    sign_1=-1;
	    		}
	    		if(py-m*px>0.)
	    		{
	    		    sign_2=1;
	    		}
	    		else
	    		{
	    		    sign_2=-1;
	    		}
	    		if(sign_1==sign_2)
	    		{
	    		    m=Cy/Cx;
	    		    if(By-m*Bx>0.)
	    		    {
	    			sign_1=1;
	    		    }
	    		    else
	    		    {
	    			sign_1=-1;
	    		    }
	    		    if(py-m*px>0.)
	    		    {
	    			sign_2=1;
	    		    }
	    		    else
	    		    {
	    			sign_2=-1;
	    		    }
	    		    if(sign_1==sign_2)
	    		    {
	    			m=(Cy-By)/(Cx-Bx);
	    			Co=Cy-m*Cx;
	    			if(-1.*Co>0.)
	    			{
	    			    sign_1=1;
	    			}
	    			else
	    			{
	    			    sign_1=-1;
	    			}
	    			if(py-m*px-Co>0.)
	    			{
	    			    sign_2=1;
	    			}
	    			else
	    			{
	    			    sign_2=-1;
	    			}
	    			if(sign_1==sign_2)
	    			{
	    			    flag=1;
	    			}
	    			else
	    			    flag=0;

	    		    }
	    		    else
	    		    {
	    			flag=0;
	    		    }
	    		}
	    		else
	    		    flag=0;
	    		if(flag)
	    		{
	    		    container_vertice *ctemp;
	    		    ctemp = new container_vertice;
	    		    ctemp->V=temp;
	    		    if(!CSTART[TYPE])
	    			CSTART[TYPE]=ctemp;
	    		    else
	    			insert_cvertice(CSTART[TYPE],ctemp,CSTART[TYPE]);
	    		}
	    		container_vertice *ctemp;
	    		ctemp = new container_vertice;
	    		ctemp->V=temp;
	    		if(!new_vert)
	    		{
	    		    //cout<<"here>\n";
	    		    new_vert=ctemp;
	    		}
	    		else
	    		{
	    		    //cout<<"here2\n";
	    		    //display_SITE(ctemp->V->p);
	    		    insert_cvertice(new_vert,ctemp,new_vert);
	    		    //cout<<"here3\n";
	    		}
	    	    }
	    	    //display_SITE(temp->p);
	    	    temp_vert_o=temp;
	    	}
	    	vertice *temp;
	    	vertice *temp_o;
	    	temp=new vertice;
	    	temp->p=new site;
	    	temp->p->x=D_TWO->circum_x;
	    	temp->p->y=D_TWO->circum_y;
	    	temp->A=SAM;
	    	temp->D=D_TWO;
	    	temp_o=temp;
	    	//cout<<"here\n";
	    	//display_SITE(temp->p);
	    	temp=V->insert_vertice(start[TYPE],temp,TYPE);
	    	//cout<<"here2\n";
	    	if(temp_o==temp)
	    	{
	    	    double Ax,Ay,Bx,By,Cx,Cy,px,py,m,Co;
	    	    int sign_1,sign_2;
	    	    int flag=1;
	    	    Ax=Atoms[temp->A].x;
	    	    Ay=Atoms[temp->A].y;
	    	    Bx=Atoms[temp->D->a].x;
	    	    By=Atoms[temp->D->a].y;
	    	    Cx=Atoms[temp->D->b].x;
	    	    Cy=Atoms[temp->D->b].y;
	    	    px=Atoms[Atom_in_foc].x;
	    	    py=Atoms[Atom_in_foc].y;
	    	  //cout<<Ax<<"\t"<<Ay<<"\n";
	    	  //cout<<Bx<<"\t"<<By<<"\n";
	    	  //cout<<Cx<<"\t"<<Cy<<"\n";
	    	  //cout<<px<<"\t"<<py<<"\n";
	    	    Bx=Bx-Ax;
	    	    By=By-Ay;
	    	    Cx=Cx-Ax;
	    	    Cy=Cy-Ay;
	    	    px=px-Ax;
	    	    py=py-Ay;
	    	    Cx=(Cx-(twob*lround(Cx/twob)));
	    	    Cy=(Cy-(twob*lround(Cy/twob)));
	    	    Bx=(Bx-(twob*lround(Bx/twob)));
	    	    By=(By-(twob*lround(By/twob)));
	    	    px=(px-(twob*lround(px/twob)));
	    	    py=(py-(twob*lround(py/twob)));
	    	    m=By/Bx;
	    	    if(Cy-m*Cx>0.)
	    	    {
	    		sign_1=1;
	    	    }
	    	    else
	    	    {
	    		sign_1=-1;
	    	    }
	    	    if(py-m*px>0.)
	    	    {
	    		sign_2=1;
	    	    }
	    	    else
	    	    {
	    		sign_2=-1;
	    	    }
	    	  //  cout<<sign_1<<"\t"<<sign_2<<"\n";
	    	    if(sign_1==sign_2)
	    	    {
	    	//	    cout<<"here\n";
	    		m=Cy/Cx;
	    		if(By-m*Bx>0.)
	    		{
	    		    sign_1=1;
	    		}
	    		else
	    		{
	    		    sign_1=-1;
	    		}
	    		if(py-m*px>0.)
	    		{
	    		    sign_2=1;
	    		}
	    		else
	    		{
	    		    sign_2=-1;
	    		}
	    	    //cout<<sign_1<<"\t"<<sign_2<<"\n";
	    		if(sign_1==sign_2)
	    		{
	    		    m=(Cy-By)/(Cx-Bx);
	    		    Co=Cy-m*Cx;
	    		    if(-1.*Co>0.)
	    		    {
	    			sign_1=1;
	    		    }
	    		    else
	    		    {
	    			sign_1=-1;
	    		    }
	    		    if(py-m*px-Co>0.)
	    		    {
	    			sign_2=1;
	    		    }
	    		    else
	    		    {
	    			sign_2=-1;
	    		    }
	    	    //cout<<sign_1<<"\t"<<sign_2<<"\n";
	    		    if(sign_1==sign_2)
	    		    {
	    			flag=1;
	    		    }
	    		    else
	    			flag=0;

	    		}
	    		else
	    		{
	    		    flag=0;
	    		}
	    	    }
	    	    else
	    		flag=0;
	    	    //cout<<flag<<"\n";
	    	    if(flag)
	    	    {
	    		container_vertice *ctemp;
	    		ctemp = new container_vertice;
	    		ctemp->V=temp;
	    		if(!CSTART[TYPE])
	    		    CSTART[TYPE]=ctemp;
	    		else
	    		    insert_cvertice(CSTART[TYPE],ctemp,CSTART[TYPE]);
	    	    }
	    	    container_vertice *ctemp;
	    	    ctemp = new container_vertice;
	    	    ctemp->V=temp;
	    	    //display_SITE(ctemp->V->p);
	    	    if(!new_vert)
	    	    {
	    		//cout<<"here>\n";
	    		new_vert=ctemp;
	    	    }
	    	    else
	    	    {
	    		//cout<<"here2\n";
	    		//display_SITE(ctemp->V->p);
	    		insert_cvertice(new_vert,ctemp,new_vert);
	    		//cout<<"here3\n";
	    	    }
	    	}
	    	temp_vert_d=temp;
	    	double m;
	    	double X,Y,dis;
	    	X=Atoms[Atoms[SAM].contigous[i][TYPE]].x-Atoms[SAM].x;
	    	Y=Atoms[Atoms[SAM].contigous[i][TYPE]].y-Atoms[SAM].y;
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
	    	if(dis>Atoms[Atoms[SAM].contigous[i][TYPE]].radius+Atoms[SAM].radius+2.*r_cut)
	    	    Atoms[SAM].bondinvoid[i][TYPE]=1;
	    	else if(sign_N == sign_C)
	    	    Atoms[SAM].bondinvoid[i][TYPE]=1;
	    	else
	    	    Atoms[SAM].bondinvoid[i][TYPE]=0;
	    	//if(Atoms[SAM].bondinvoid[i])
	    	{
	    	    vor<<std::setprecision(15)<<x<<"\t"<<y<<"\t"<<std::flush;
	    	    vor<<std::setprecision(15)<<D_ONE->circum_x<<"\t"<<D_ONE->circum_y<<"\t"<<D_TWO->circum_x<<"\t"<<D_TWO->circum_y<<"\n"<<std::flush;
	    	    vor<<"\n"<<std::flush;
	    	}
	    	add_connected(temp_vert_o,temp_vert_d,Atoms[SAM].bondinvoid[i][TYPE]);
	    	add_connected(temp_vert_d,temp_vert_o,Atoms[SAM].bondinvoid[i][TYPE]);
	        }
	        //cout<<"here\n";
	    }
	    //display_SITE(new_vert->V->p);
	    //display_SITE(new_vert->next->V->p);
	    //cout<<"eye "<<Atoms[35].conti[TYPE]<<"\n";
	    //if(i==8)
	    //{
	    //      delunay *D;
	    //      D=Atoms[2].D[1].initial;
	    //      while(1)
	    //      {
	    //          print_delunay(&(Atoms[2]),D,Atoms,1);
	    //          if(D->next)
	    //          {
	    //              D=D->next;
	    //          }
	    //          else
	    //              break;
	    //      }
	    //  }
	    container_vertice *cstart;
	    cstart=CSTART[TYPE];
	    container_vertice *temp_new_vert;
	    temp_new_vert=new_vert;
	    while(1)
	    {
		   // /cout<<"####\n";
		   // display_SITE(temp_new_vert->V->p);
		  //  for(int n=0;n<temp_new_vert->V->v_neigh_count;n++)
		   // {
		//	    cout<<n<<"\t";
		//	display_SITE(temp_new_vert->V->neib_vert[n]->p);
		  //  }
		if(temp_new_vert->next)
		{
			temp_new_vert=temp_new_vert->next;
		}
		else
			break;
	    }
	    temp_new_vert=new_vert;
	    while(1)
	    {
		if(!(compare(temp_new_vert->V->p,cstart->V->p)))
		{
			if((temp_new_vert->V->p->y-cstart->V->p->y)||(temp_new_vert->V->p->x-cstart->V->p->x))
			{
				if(temp_new_vert->V->p->y>0.)
				{
					cstart->V=temp_new_vert->V;
				}
			}

		}		
		if(temp_new_vert->next)
		{
			temp_new_vert=temp_new_vert->next;
		}
		else
			break;
	    }
	    vertice *temp_start;
	    container_vertice *temp_new_vert1;
	    temp_new_vert=new_vert;
	    int rep;
	    while(1)
	    {
	    	    temp_new_vert1=new_vert;
		    while(1)
		    {
	//		    cout<<"####\n";
	//		                display_SITE(temp_new_vert->V->p);
	//		        	display_SITE(temp_new_vert1->V->p);
		        rep=compare(temp_new_vert->V->p,temp_new_vert1->V->p);
		        if(rep==0)
		        {
				if((temp_new_vert1->V->p->y-temp_new_vert->V->p->y)||(temp_new_vert1->V->p->x-temp_new_vert->V->p->x))
				{
					if(temp_new_vert1->V->p->y>0.)
					{
						for(int n=0;n<temp_new_vert->V->v_neigh_count;n++)
						{
							add_connected(temp_new_vert1->V,temp_new_vert->V->neib_vert[n],temp_new_vert->V->neib_ed[n]);
							add_connected(temp_new_vert->V->neib_vert[n],temp_new_vert1->V,temp_new_vert->V->neib_ed[n]);
						}
					        if(temp_new_vert->prev)
					        {
					        	if(temp_new_vert->next)
					        	{
					        		temp_new_vert->prev->next=temp_new_vert->next;
					        		temp_new_vert->next->prev=temp_new_vert->prev;
					        	}
					        	else
					        	{
					        		temp_new_vert->prev->next=NULL;
					        	}
					        }
					        else
					        {
					        	new_vert=temp_new_vert->next;
					        }

	        				V->delete_vertice(start[TYPE],temp_new_vert->V,TYPE);
					        delete temp_new_vert->V->p;
					        delete temp_new_vert->V;
					        delete temp_new_vert;
					}
					else
					{
						for(int n=0;n<temp_new_vert1->V->v_neigh_count;n++)
						{
	//						cout<<n<<"\n";
							if(temp_new_vert1->V->neib_vert[n])
							{
							add_connected(temp_new_vert->V,temp_new_vert1->V->neib_vert[n],temp_new_vert1->V->neib_ed[n]);
							add_connected(temp_new_vert1->V->neib_vert[n],temp_new_vert->V,temp_new_vert1->V->neib_ed[n],0);
							}
						}
					        if(temp_new_vert1->prev)
					        {
					        	if(temp_new_vert1->next)
					        	{
					        		temp_new_vert1->prev->next=temp_new_vert1->next;
					        		temp_new_vert1->next->prev=temp_new_vert1->prev;
					        	}
					        	else
					        	{
					        		temp_new_vert1->prev->next=NULL;
					        	}
					        }
					        else
					        {
					        	new_vert=temp_new_vert1->next;
					        }
						container_vertice *temp;
						temp=temp_new_vert1;
	        				V->delete_vertice(start[TYPE],temp_new_vert1->V,TYPE);
					        delete temp->V->p;
					        delete temp->V;
						
					}
				}
		        }
			if(temp_new_vert1->next)
			{
				temp_new_vert1=temp_new_vert1->next;
			}
			else 
				break;
		    }
		    if(temp_new_vert->next)
		    {
		        temp_new_vert=temp_new_vert->next;
		    }
		    else 
		    	break;
	    }
	    //cout<<r_cut<<"\n";
	    temp_new_vert=new_vert;
	    while(1)
	    {
	        temp_start=temp_new_vert->V;
		//display_SITE(temp_start->p);
	        int flag=1;
	        double AX,AY,BX,BY,X,Y;
	        BX=temp_start->p->x ;
	        BY=temp_start->p->y;
	        AX=Atoms[temp_start->A].x;
	        AY=Atoms[temp_start->A].y;
	        double dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
	        if(dis<r_cut+Atoms[temp_start->A].radius)
	        {
	    	flag=0;
	    	//break;
	        }
	        AX=Atoms[Atoms[temp_start->A].contigous[temp_start->D->A][TYPE]].x;
	        AY=Atoms[Atoms[temp_start->A].contigous[temp_start->D->A][TYPE]].y;
	        dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
	        if(dis<r_cut+Atoms[Atoms[temp_start->A].contigous[temp_start->D->A][TYPE]].radius)
	        {
	    	flag=0;
	    	//break;
	        }
	        AX=Atoms[Atoms[temp_start->A].contigous[temp_start->D->B][TYPE]].x;
	        AY=Atoms[Atoms[temp_start->A].contigous[temp_start->D->B][TYPE]].y;
	        dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
	        //cout<<dis<<"\t"<<AX<<"\t"<<AY<<"\t"<<temp_start->D->B<<"\n";;
	        if(dis<r_cut+Atoms[Atoms[temp_start->A].contigous[temp_start->D->B][TYPE]].radius)
	        {
	    	flag=0;
	    	//break;
	        }
	        AX=Atoms[i].x;
	        AY=Atoms[i].y;
	        BX=temp_start->p->x;
	        BY=temp_start->p->y;
	        X=AX-BX;
	        Y=AY-BY;
	        X=(X-(twob*lround(X/twob)));
	        Y=(Y-(twob*lround(Y/twob)));
	        dis=sqrt((X)*(X)+(Y)*(Y));
	        if(flag)//&& temp_start_start->v_neigh_count == 3)
	        {
		  //  cout<<"####\n";
	    //	display_SITE(temp_start->p);
		  //  cout<<temp_start->D->Ax<<"\t"<<temp_start->D->Ay<<"\n";
		  //  cout<<temp_start->D->Bx<<"\t"<<temp_start->D->By<<"\n";
	    	////////if(dis<Atoms[i].radius+r_cut)
	    	////////	cout<<i<<"\t"<<"this one has a vertice inside\n";
	    	temp_start->is_void=1;
	    	///temp_start->cluster_index=void_vert_count;
	        }
	        if(temp_new_vert->next)
	    	temp_new_vert=temp_new_vert->next;
	        else
	    	break;
	    }
	    cstart=CSTART[TYPE];
	////if(!new_vert)
	////	new_vert=cstart;
	    int void_vert_count=0;
	    temp_new_vert=new_vert;
	    ///cout<<"here4\n";
	  //while(1)
	  //{
	  //    //temp_new_vert->V=cstart->V;
	  //    //print_delunay(&Atoms[cstart->V->A],cstart->V->D,Atoms,TYPE);
	  //    display_SITE(cstart->V->p);
	  //    if(cstart->next)
	  //    {
	  //        cstart=cstart->next;
	  //        //  temp_new_vert->next=new container_vertice;
	  //        //  temp_new_vert=temp_new_vert->next;
	  //    }
	  //    else
	  //        break;
	  //}
	  //cout<<" ... \n";
	  //while(1)
	  //{
	  //    //temp_new_vert->V=cstart->V;
	  //    //print_delunay(&Atoms[cstart->V->A],cstart->V->D,Atoms,TYPE);
	  //    //display_SITE(cstart->V->p);
	  //    display_SITE(temp_new_vert->V->p);
	  //    if(temp_new_vert->next)
	  //    {
	  //        temp_new_vert=temp_new_vert->next;
	  //    }
	  //    else
	  //        break;
	  //}
	    cstart=CSTART[TYPE];

	    while(1)
	    {
	        if(cstart->V->is_void)
	        {
	    	//display_SITE(cstart->V->p);
	    	void_vert_count++;
	        }
	        if(cstart->next)
	    	cstart=cstart->next;
	        else
	    	break;
	    }
	    int change=1;
	    int void_vert_count_prev;
	    int tem_ind=0;
	    while(change)
	    {
	        cstart=CSTART[TYPE];
	        void_vert_count_prev=void_vert_count;
	        while(1)
	        {
	    	//cout<<cstart->V->v_neigh_count<<"\n";
	    	if(cstart->V->is_void)
	    	{
	    	    //display_SITE(cstart->V->p);
	    	    //cout<<"this\n";
	    	    for(int i=0; i<cstart->V->v_neigh_count; i++)
	    	    {
	    		//cout<<"here2\n";
	    		//cout<<"here\n";
	    		if(cstart->V->neib_vert[i])
	    		{
	    		    //display_SITE(cstart->V->neib_vert[i]->p);
	    		    //cout<<cstart->V->neib_ed[i]<<"\n";
	    		    if(cstart->V->neib_vert[i]->is_void && cstart->V->neib_ed[i])
	    		    {
	    			////////cout<<"here\n";
	    			////////display_SITE(cstart->V->neib_vert[i]->p);
	    			////////cout<<"here2\n";
	    			int flag;
	    			container_vertice *ctemp;
	    			ctemp=new container_vertice;
	    			ctemp->V=cstart->V->neib_vert[i];
	    			flag=insert_cvertice(CSTART[TYPE],ctemp,CSTART[TYPE]);
	    			if(flag)
	    			{
	    			    void_vert_count=void_vert_count+1;
	    			}
	    		    }
	    		}
	    		else
	    		{
	    		    //cout<<"who\n";
	    		}
	    	    }
	    	}
	    	if(cstart->next)
	    	{
	    	    cstart=cstart->next;
	    	}
	    	else
	    	    break;
	        }
	        //tem_ind++;
	        //if(tem_ind==1)
	        //	break;
	        if(void_vert_count_prev==void_vert_count)
	    	change=0;
	        else
	    	change=1;
	    }
	    cstart=CSTART[TYPE];
	    void_vert_count=0;
	    while(1)
	    {
	        if(cstart->V->is_void)
	        {
	    		void_vert_count++;
	        }
	        if(cstart->next)
	    		cstart=cstart->next;
	        else
	            break;
	    }
	    //cout<<void_vert_count<<"void_count\n"<<std::flush;
	    vertice **cavity_list;
	    cavity_list = new (nothrow) vertice*[void_vert_count];
	    cstart=CSTART[TYPE];
	    {
	        int i=0;
	        while(1)
	        {
	    	if(cstart->V->is_void)
	    	{
	    	    cavity_list[i]=cstart->V;
	    	    //cout<<i<<"\n";
	    	    //cout<<"here\n";
	    	    i++;
	    	    //cout<<i<<"\t"<<cavity_list[i]<<"\t"<<temp_start<<"\n";;
	    	}
	    	if(cstart->next)
	    	    cstart=cstart->next;
	    	else
	    	    break;

	        }

	    }
	    //return 0;
	    //CALCULATING THE VOID VOLUME IN A GIVEN CAVITY
	    double *void_area;
	    void_area= new (nothrow) double[void_vert_count];
	    double *void_length;
	    void_length= new (nothrow) double[void_vert_count];
	    //return 0;
	    for(int j=0 ; j<void_vert_count; j++)
	    {
	        void_area[j]=0;
	        void_length[j]=0;
	        //display_SITE(cavity_list[j]->p);
	        //cout<<cavity_list[j]->A<<"\t"<<cavity_list[j]->D->A<<"\t"<<cavity_list[j]->D->B<<"\n";
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
	        A1r=Atoms[cavity_list[j]->A].radius+r_cut;
	        A2x=Atoms[cavity_list[j]->D->a].x-A1x;
	        A2y=Atoms[cavity_list[j]->D->a].y-A1y;
	        A2r=Atoms[cavity_list[j]->D->a].radius+r_cut;
	        A3x=Atoms[cavity_list[j]->D->b].x-A1x;
	        A3y=Atoms[cavity_list[j]->D->b].y-A1y;
	        A3r=Atoms[cavity_list[j]->D->b].radius+r_cut;
		//cout<<cavity_list[j]->D->Bx<<"\t"<<cavity_list[j]->D->By<<"\n";;
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
	////////cout<<j<<"\n";
	////////cout<<A1x<<"\t"<<A1y<<"\t"<<A1r<<"\n";;
	////////cout<<A2x+A1x<<"\t"<<A2y+A1y<<"\t"<<A2r<<"\n";
	////////cout<<A3x+A1x<<"\t"<<A3y+A1y<<"\t"<<A3r<<"\n";
	        Vx=(Vx-(twob*lround(Vx/twob)));
	        Vy=(Vy-(twob*lround(Vy/twob)));
	        double X,Y,x,y;
	        double rA,rS,DIS,l,dis_i,tan_sq;
	        X=Atoms[cavity_list[j]->D->a].x-Atoms[cavity_list[j]->D->b].x;
	        Y=Atoms[cavity_list[j]->D->a].y-Atoms[cavity_list[j]->D->b].y;
	        X=(X-(twob*lround(X/twob)));
	        Y=(Y-(twob*lround(Y/twob)));
	        DIS=sqrt(X*X+Y*Y);
	        rA=A2r;
	        rS=A3r;
	        l=0.5*(DIS+(rS*rS-rA*rA)/DIS);
	        x=l/DIS*X;
	        y=l/DIS*Y;
	        x=x+Atoms[cavity_list[j]->D->b].x;
	        y=y+Atoms[cavity_list[j]->D->b].y;
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

		//cout<<j<<"\n";
	        //cout<<Vx+A1x<<"\t"<<Vy+A1y<<"\n";
	        //cout<<E12x+A1x<<"\t"<<E12y+A1y<<"\n";
	        //cout<<E13x+A1x<<"\t"<<E13y+A1y<<"\n";
	        //cout<<E23x+A1x<<"\t"<<E23y+A1y<<"\n";
	        //cout<<A1x<<"\t"<<A1y<<"\n";
	        //cout<<A2x+A1x<<"\t"<<A2y+A1y<<"\n";
	        //cout<<A3x+A1x<<"\t"<<A3y+A1y<<"\n";
	        //cout<<"begin\n";
	        //cout<<void_area[j]<<"\n";
	        void_area[j]=void_area[j]+S12*area_trangle(Vx,Vy,E12x,E12y,A2x,A2y,A2r,A1x,A1y);
	        //cout<<void_area[j]<<"\n";
	        void_area[j]=void_area[j]+S12*area_trangle(Vx,Vy,E12x,E12y,0.,0.,A1r,A1x,A1y);
	        //cout<<void_area[j]<<"\n";
	        void_area[j]=void_area[j]+S13*area_trangle(Vx,Vy,E13x,E13y,A3x,A3y,A3r,A1x,A1y);
	        //cout<<void_area[j]<<"\n";
	        void_area[j]=void_area[j]+S13*area_trangle(Vx,Vy,E13x,E13y,0.,0.,A1r,A1x,A1y);
	        //cout<<void_area[j]<<"\n";
	        void_area[j]=void_area[j]+S23*area_trangle(Vx,Vy,E23x,E23y,A3x,A3y,A3r,A1x,A1y);
	        //cout<<void_area[j]<<"\n";
	        void_area[j]=void_area[j]+S23*area_trangle(Vx,Vy,E23x,E23y,A2x,A2y,A2r,A1x,A1y);
	        //cout<<void_area[j]<<"here\n";
	        void_length[j]=void_length[j]+S12*perimeter(Vx,Vy,E12x,E12y,A2x,A2y,A2r);
	        void_length[j]=void_length[j]+S12*perimeter(Vx,Vy,E12x,E12y,0.,0.,A1r);
	        void_length[j]=void_length[j]+S13*perimeter(Vx,Vy,E13x,E13y,A3x,A3y,A3r);
	        void_length[j]=void_length[j]+S13*perimeter(Vx,Vy,E13x,E13y,0.,0.,A1r);
	        void_length[j]=void_length[j]+S23*perimeter(Vx,Vy,E23x,E23y,A3x,A3y,A3r);
	        void_length[j]=void_length[j]+S23*perimeter(Vx,Vy,E23x,E23y,A2x,A2y,A2r);
	        //cout<<j<<"\t"<<void_length[j]<<"\t"<<void_area[j]<<"\n";

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
	        //print_delunay(&Atoms[cavity_list[j]->A],cavity_list[j]->D,Atoms,TYPE);
	        //cout<<cavity_list[j]->p->x<<"\t"<<cavity_list[j]->p->y<<"\n";
	    }
	    //delete [] cavity_list;
	    double cav_tot=0.;
	    double ca_per_tot=0.;
	    for(int i=0; i<void_vert_count; i++)
	    {
	        cav_tot=cav_tot+void_area[i];
	        ca_per_tot=ca_per_tot+void_length[i];
	        //cout<<i<<"\t"<<void_area[i]<<"\t"<<void_length[i]<<"\n";
	    }
	    freearea[i]=cav_tot;
	    freeperi[i]=ca_per_tot;
	    //cout<<cav_tot<<"\t"<<ca_per_tot<<"\n";
	    temp_new_vert=new_vert;
	    while(1)
	    {
	        V->delete_vertice(start[TYPE],temp_new_vert->V,TYPE);
	        //display_SITE(temp_new_vert->V->p);
	        if(temp_new_vert->next)
	        {
	    	temp_new_vert=temp_new_vert->next;
	        }
	        else
	    	break;
	    }
	    ctemp=Atoms[i].Cstart[TYPE];
	    while(1)
	    {
	        V->insert_vertice(start[TYPE],ctemp->V,TYPE);
	      //for(int n=0;n<ctemp->V->v_neigh_count;n++)
	      //{
	      //	if(ctemp->V->neib_vert[n])
	      //		display_SITE(ctemp->V->neib_vert[n]->p);	
	      //	//if(
	      //	//add_connected(ctemp->V->neib_vert[n],ctemp->V,ctemp->V->neib_ed[n],0);
	      //}
	        //V->delete_vertice(start[TYPE],ctemp->V,TYPE);
	        if(ctemp->next)
	    	ctemp=ctemp->next;
	        else
	    	break;
	    }
	    //cout<<"here2\n"<<std::flush;
	    //container_vertice *ctemp;
	    //for(int j=0; j<Atoms[i].conti[TYPE]; j++)
	    //{
	    //    int SAM=Atoms[i].contigous[j][TYPE];
            delunay *D;
            for(int j=0; j<Atoms[i].conti[TYPE]; j++)
            {
        	  //  cout<<j<<"\n"<<std::flush;
                D=Atoms[i].D[TYPE].initial;
                delunay *D_ONE=NULL;
                delunay *D_TWO=NULL;

                while(1)
                {
            	if(D->A==j)
            	{
            	    //cout<<D->circum_x<<"\t"<<D->circum_y<<"\t";
            	    if(!D_ONE)
            	    {
            		D_ONE=D;
            	    }
            	    else if(!D_TWO)
            		D_TWO=D;
            	}
            	if(D->B==j)
            	{
            	    //cout<<D->circum_x<<"\t"<<D->circum_y<<"\t";
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
        	///cout<<"her1\n"<<std::flush;
        	vertice *temp1;
        	vertice *temp2;
        	temp1=new vertice;
        	temp2=new vertice;
        	temp1->p=new site;
        	temp2->p=new site;
        	//cout<<"here1\n"<<std::flush;
        	temp1->p->x=D_ONE->circum_x;
        	temp1->p->y=D_ONE->circum_y;
        	temp2->p->x=D_TWO->circum_x;
        	temp2->p->y=D_TWO->circum_y;
            	temp1=V->insert_vertice(start[TYPE],temp1,TYPE);
            	temp2=V->insert_vertice(start[TYPE],temp2,TYPE);
                add_connected(temp1,temp2,Atoms[i].bondinvoid[j][TYPE],0);
                add_connected(temp2,temp1,Atoms[i].bondinvoid[j][TYPE],0);
        //cout<<"\n";
            }
////////    cout<<"here\n"<<std::flush;
            ctemp=Atoms[i].Cstart[TYPE];
            while(1)
            {
                for(int n=0;n<ctemp->V->v_neigh_count;n++)
                {
                ////////if(ctemp->V->neib_vert[n])
                ////////	display_SITE(ctemp->V->neib_vert[n]->p);	
                        if(ctemp->V->neib_vert[n])
                		add_connected(ctemp->V->neib_vert[n],ctemp->V,ctemp->V->neib_ed[n],0);
                }
                //V->delete_vertice(start[TYPE],ctemp->V,TYPE);
                if(ctemp->next)
            	ctemp=ctemp->next;
                else
            	break;
            }
	    
	  //temp_start=start[TYPE];
	  //while(1)
	  //{
	  //        if(temp_start->v_neigh_count!=3)
	  //        {
	  //		//display_SITE(temp_start->p);
	  //    //cout<<temp_start->v_neigh_count<<"\n";
	  //    for(int n=0;n<temp_start->v_neigh_count;n++)
	  //    {
	  //    	if(temp_start->neib_vert[n])
	  //    	display_SITE(temp_start->neib_vert[n]->p);
	  //    }
	  //    cout<<"#####\n";
	  //        }
	  //	if(temp_start->next)
	  //	{
	  //		temp_start=temp_start->next;
	  //	}
	  //	else
	  //		break;
	  //}
	    //cout<<"after insertion\n";
	    //}

	    for(int j=0; j<Atoms[i].conti[TYPE]; j++)
	    {
	        int flag=0;
	        int SAM=Atoms[i].contigous[j][TYPE];
	        delunay *D;
	        D=Atoms[SAM].D[TYPE].initial;
	        ////////cout<<n<<"\n";
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
	        reset_atom(&(Atoms[SAM]),TYPE);
	    }
	    delete[] void_area;
	    delete[] void_length;
	    delete[] cavity_list;
	    if(new_vert)
	    {
	        cstart=new_vert;
	        while(1)
	        {
	    	container_vertice *ctemp;
	    	ctemp=cstart;
	    	if(cstart->next)
	    	{
	    	    cstart=cstart->next;
	    	    delete ctemp->V->p;
	    	    delete ctemp->V;
	    	    delete ctemp;
	    	}
	    	else
	    	{
	    	    delete ctemp->V->p;
	    	    delete ctemp->V;
	    	    delete ctemp;
	    	    break;
	    	}

	        }
	    }
	    if(CSTART[TYPE])
	    {
	        cstart=CSTART[TYPE];
	        while(1)
	        {
	    	container_vertice *ctemp;
	    	ctemp=cstart;
	    	if(cstart->next)
	    	{
	    	    cstart=cstart->next;
	    	    delete ctemp;
	    	}
	    	else
	    	{
	    	    delete ctemp;
	    	    break;
	    	}

	        }
	    }
	    vor.close();
	}
	//DISTRIBUTION CALCULATION
	cout<<"are we here3\n"<<std::flush;
	if(nconfig==0)
	{
		for(int f=0;f<nAtoms;f++)
		{
			if(max_free_area<freearea[f])
				max_free_area=freearea[f];
		}
		width=(max_free_area+1.0)/1000.;
	}
	for(int f=0;f<nAtoms;f++)
	{
		//cout<<f<<"\t"<<int(freearea[f]/width)<<std::flush<<"\n";
		if(int(freearea[f]/width)<1000)
			free_dist[int(freearea[f]/width)]++;
	}
	//cout<<"ahere\n";
	double *sum_freearea=new (nothrow) double [ntypes];
	double *sum_freeperi=new (nothrow) double [ntypes];
	double *sum_freeratio=new (nothrow) double [ntypes];
	int *count=new (nothrow) int [ntypes];
	for(int t=0;t<ntypes;t++)
	{
	        count[t]=0;
		sum_freearea[t]=0.;
		sum_freeperi[t]=0.;
		sum_freeratio[t]=0.;
	}
	for(int i=0;i<nAtoms;i++)
	{
	        sum_freearea[Atoms[i].type]=sum_freearea[Atoms[i].type]+freearea[i];
	        sum_freeperi[Atoms[i].type]=sum_freeperi[Atoms[i].type]+freeperi[i];
	        sum_freeratio[Atoms[i].type]=sum_freeratio[Atoms[i].type]+freeperi[i]/freearea[i];
	        count[Atoms[i].type]++;
	}
        for(int t=0;t<ntypes;t++)
        {
                cout<<t<<"\t"<<sum_freeratio[t]/count[t]<<"\t"<<std::flush;///<<sum_freeperi[t]<<"\t"<<count[t]<<"\t";
        }
	delete[] sum_freearea;
	delete[] sum_freeperi;
	delete[] sum_freeratio;
	delete[] count;
	cout<<"\n";

	//return 0;
	delete[] freearea;
	delete[] freeperi;
	for(int i=0;i<nAtoms;i++)
	{
	    ////////delete[] Atoms[i].D;
	        delete[] Atoms[i].conti;
	        delete[] Atoms[i].save_conti;
	        for(int t=0;t<500;t++)
	        {
	    	    delete [] Atoms[i].contigous[t];
	    	    delete [] Atoms[i].edge_index[t];
	    	    delete [] Atoms[i].bondinvoid[t];
	    	    delete [] Atoms[i].save_contigous[t];
	    	    delete [] Atoms[i].save_edge_index[t];
	    	    delete [] Atoms[i].save_bondinvoid[t];
	        }
	}
cout<<"are we here5\n"<<std::flush;
	for(int t=0; t<ntypes; t++)
	{
	    container_vertice *cstart;
	    vertice *temp_start;
	    temp_start=start[t];
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
	    for(SAM=0; SAM<nAtoms; SAM++)
	    {
	        cstart=Atoms[SAM].Cstart[t];
	        while(1)
	        {
	    	container_vertice *temp_c;
	    	temp_c=cstart;
	    	if(cstart->next)
	    	{
	    	    cstart=cstart->next;
	    	    delete temp_c;
	    	}
	    	else
	    	{
	    	    delete temp_c;
	    	    break;
	    	}
	        }
	    }
	}
	cout<<"are we here4\n"<<std::flush;
	for(int n=0; n<nAtoms; n++)
	{
	    delete[] Atoms[n].Cstart;
	    delunay *D;
	    for(int t=0; t<ntypes; t++)
	    {
	        D=Atoms[n].D[t].initial;
	        ////////cout<<n<<"\n";
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
	    delete [] Atoms[n].D;
	    delete [] Atoms[n].save_D;
	    ////while(1)
	    ////{

	    ////}
	}
	vertice	*temp_start=sites;
	while(1)
	{
	    vertice *temp;
	    temp=temp_start;
	    //display_SITE(temp->p);
	    if(temp_start->next)
	    {
		temp_start=temp_start->next;
		delete temp->p;
		delete temp;
		//temp=NULL;
	    }
	    else
	    {
		delete temp->p;
		delete temp;
		//temp=NULL;
		break;
	    }
	}
	sites=NULL;
        delete[] Atoms;
	for(int i=0;i<1000;i++)
	{
	        fdist<<i*width<<"\t"<<free_dist[i]/(config_count*nAtoms*width)<<"\n"<<std::flush;
	}
    }
    delete[] radius;

    return 0;
}
