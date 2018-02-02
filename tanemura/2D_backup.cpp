#include<iostream>
#include<math.h>
#include<fstream>
#include<ostream>
#include<iomanip>
#include<vector>
#include <set>
#include <limits>
#include <functional>
#include <numeric>
#include <algorithm>
using namespace std;
#define distance(x,y) x*x+y*y
int dim=3;
long double box=0;
long double twob;
long double density;
long double r_cut=0.0;
long double tilt=0.0;
long double DMIN=0.000000000001;//std::numeric_limits<long double>::min();
long double epsilon=0.00000000000;
struct atom;
struct face;
struct vertice;
struct site
{
    long double x=0;
    long double y=0;
    long double z=0;
    struct face *F=NULL;
};
void display_SITE(struct site *p)
{
    cout<<p->x<<"\t"<<p->y<<"\t"<<p->z<<"\n";
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
    void delete_vertice(vertice *,vertice *,int type,int);
    vertice* insert_vertice(vertice *,vertice *,int,int);
    void display_vertex(vertice *);
    void display_conn(vertice *ve);
}*V;
struct delunay
{
    int A;
    int B;
    int C;
    int D;
	int ABf=0;
	int BCf=0;
	int CAf=0;
	long double ABx=0;
	long double ABy=0;
	long double ABz=0;
	long double BCx=0;
	long double BCy=0;
	long double BCz=0;
	long double CAx=0;
	long double CAy=0;
	long double CAz=0;
    int a;
    int b;
    int c;
    long double circum_x=0.;
    long double circum_y=0.;
    long double circum_z=0.;
    long double Ax=0.;
    long double Ay=0.;
    long double Az=0.;
    long double Bx=0.;
    long double By=0.;
    long double Bz=0.;
    long double Cx=0.;
    long double Cy=0.;
    long double Cz=0.;
    delunay *next=NULL;
};
struct vertice
{
    struct site *p=NULL;
    struct half_edge *leaving=NULL;
    struct vertice *next=NULL;
    struct vertice *prev=NULL;
    int A,B,C,Da;
    delunay *D=nullptr;
    int is_void=0;
    int cluster_index=-1;
    long double r;
    //vert_list *V;
    vertice *neib_vert[10]= {NULL};
    int neib_ed[10];
    int v_neigh_count=0;
}**start,*s_temp,*sites;
struct container_vertice
{
    struct vertice *V=nullptr;
    struct container_vertice *next=NULL;
    struct container_vertice *prev=NULL;
}**CSTART,*VOID_START;
int compare(struct site *p1,struct site *p2,int del=0)
{
    long double DX,DY,DZ;
    DX=p1->x-p2->x;
    DY=p1->y-p2->y;
    DZ=p1->z-p2->z;
    if(!del)
    {
        DX=(DX-(tilt*lroundl(DY/twob)));
        DX=(DX-(twob*lroundl(DX/twob)));
        DY=(DY-(twob*lroundl(DY/twob)));
        DZ=(DZ-(twob*lroundl(DZ/twob)));
    }
    if( abs(DX) <DMIN && abs(DX) >DMIN )
    {
        cout<<p1->x<<"\t"<<p1->y<<"\n";
        cout<<p2->x<<"\t"<<p2->y<<"\n";
    }
    if( abs(DY) <DMIN  && abs(DX) <DMIN && abs(DZ) < DMIN )
    {
        return 0;
    }
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
        cout<<"begin\n";
        cout<<focus<<"\n";
        display_SITE(focus->p);
        display_SITE(add->p);
    }
    for(i=0; i<focus->v_neigh_count; i++)
    {
        if(focus->neib_vert[i])
            if(!compare(focus->neib_vert[i]->p,add->p))
            {
                std::vector<int> EV1 {focus->neib_vert[i]->A,focus->neib_vert[i]->D->a,focus->neib_vert[i]->D->b};
                std::vector<int> v1 {add->A,add->D->a,add->D->b};
                std::sort(EV1.begin(),EV1.end());
                std::sort(v1.begin(),v1.end());
                if(EV1[0]==v1[0] && EV1[1]==v1[1] && EV1[2]==v1[2])
                {
                    if(debug)
                    {
                        cout<<"here\n";
                    }
                    if(focus->neib_vert[i]->p->y-add->p->y)
                    {
                        if(focus->neib_vert[i]->p->y<0.)
                        {
                            if(debug)
                                cout<<"replaced\n";
                            focus->neib_vert[i]=add;
                        }
                    }
                    flag=0;
                    break;
                }
            }
    }
    if(flag)
    {
        if(debug)
            cout<<"added\n";
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
    long double x=0;
    long double y=0;
    long double z=0;
    int neighlist[500];
    int contigous[500];
    int bondinvoid[500];
    int edge_index[500]= {0};
    int neighbours=0;
    int conti=0;
    struct face *F=NULL;
    long double radius=1.;
    set_of_delunay D;
};
struct atom
{
    long double x=0;
    long double y=0;
    long double z=0;
    int neighlist[500];
    int *contigous[500];
	int *part_c[500][100];
    int *bondinvoid[500];
    int *edge_index[500];
    int neighbours=0;
    int *conti;
    struct face *F=NULL;
    long double radius=1.;
    int ignore=0;
    container_vertice **Cstart=nullptr;
    set_of_delunay *D=nullptr;
    int save_neighlist[500];
    int *save_contigous[500];
    int *save_bondinvoid[500];
    int *save_edge_index[500];
    int save_neighbours=0;
    int *save_conti=nullptr;
    int type;
    set_of_delunay *save_D=nullptr;
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
int insert_cvertice(container_vertice *EV,container_vertice *v,container_vertice *&Cstart,int debug=0)
{
    int flag;
    flag=compare(EV->V->p,v->V->p);
    if(debug)
    {
        cout<<"insert\n";
        display_SITE(v->V->p);
        cout<<"vetfoc\n";
        display_SITE(EV->V->p);
        cout<<flag<<"\n";
    }
    if(flag==0)
    {
        std::vector<int> EV1 {EV->V->A,EV->V->D->a,EV->V->D->b};
        std::vector<int> v1 {v->V->A,v->V->D->a,v->V->D->b};
        std::sort(EV1.begin(),EV1.end());
        std::sort(v1.begin(),v1.end());
        if(EV1[0]==v1[0] && EV1[1]==v1[1] && EV1[2]==v1[2])
        {
            flag=0;
            ////////delete v->p;
            ////////delete v;
            ////////return EV;
        }
        else
            flag=1;

    }
    if(flag==1)
    {
        if(EV->next)
        {
            //		cout<<"everytime?\n";
            insert_cvertice(EV->next,v,Cstart,debug);
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
        if(EV->V->p->y-v->V->p->y)
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
        }
        else
        {
            delete v;
        }
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
        cout<<flag<<"\n";
    }
    //else if (flag==0)
    //{
    //        flag=-1;
    //}

    if(flag==0)
    {
        if(debug)
        {
            cout<<"here\n";
        }
        std::vector<int> EV1 {EV->A,EV->D->a,EV->D->b};
        std::vector<int> v1 {v->A,v->D->a,v->D->b};
        std::sort(EV1.begin(),EV1.end());
        std::sort(v1.begin(),v1.end());
        if(debug)
        {
            cout<<"here\n";
            cout<<EV1[0]<<"\t"<<EV1[1]<<"\t"<<EV1[2]<<"\n";
            cout<<v1[0]<<"\t"<<v1[1]<<"\t"<<v1[2]<<"\n";
        }
        if(EV1[0]==v1[0] && EV1[1]==v1[1] && EV1[2]==v1[2])
        {
            delete v->p;
            delete v;
            return EV;
        }
        else
            flag=1;

    }

    //if(debug)
    //{
    //    cout<<"flag\t"<<flag<<"\n";
    //}
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
}
void vert_list::delete_vertice(vertice *EV,vertice *v,int type,int debug=0)
{
    int flag;
    flag=compare(EV->p,v->p,1);
    if(flag==0)
    {
        std::vector<int> EV1 {EV->A,EV->D->a,EV->D->b};
        std::vector<int> v1 {v->A,v->D->a,v->D->b};
        std::sort(EV1.begin(),EV1.end());
        std::sort(v1.begin(),v1.end());
        if(EV1[0]==v1[0] && EV1[1]==v1[1] && EV1[2]==v1[2])
        {
            flag=0;
            ////////delete v->p;
            ////////delete v;
            ////////return EV;
        }
        else
            flag=1;

    }
    if(debug)
    {
        cout<<flag<<"\n";
        display_SITE(EV->p);
        display_SITE(v->p);
        cout<<"\n"<<std::flush;
    }
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
        if(debug)
        {
            cout<<v->v_neigh_count<<"\n"<<std::flush;
        }
        for(int n=0; n<v->v_neigh_count; n++)
        {
            if(debug)
            {
                cout<<n<<"\n";
            }
            if(v->neib_vert[n])
            {
                //display_SITE(v->p);
                //	    display_SITE(v->neib_vert[n]->p);
                for(int m=0; m<v->neib_vert[n]->v_neigh_count; m++)
                {
                    if(v->neib_vert[n]->neib_vert[m])
                    {
                        if(!(compare(v->neib_vert[n]->neib_vert[m]->p,v->p,1)))
                        {

                            std::vector<int> EV1 {v->neib_vert[n]->neib_vert[m]->A,v->neib_vert[n]->neib_vert[m]->D->a,v->neib_vert[n]->neib_vert[m]->D->b};
                            std::vector<int> v1 {v->A,v->D->a,v->D->b};
                            std::sort(EV1.begin(),EV1.end());
                            std::sort(v1.begin(),v1.end());
                            if(EV1[0]==v1[0] && EV1[1]==v1[1] && EV1[2]==v1[2])
                            {
                                flag=0;
                                v->neib_vert[n]->neib_vert[m]=NULL;
                                if(m==0)
                                {
                                    v->neib_vert[n]->neib_vert[m]=v->neib_vert[n]->neib_vert[m+1];
                                    v->neib_vert[n]->neib_vert[m+1]=v->neib_vert[n]->neib_vert[m+2];
                                    v->neib_vert[n]->neib_vert[m+2]=v->neib_vert[n]->neib_vert[m+3];
                                }
                                if(m==1)
                                {
                                    v->neib_vert[n]->neib_vert[m]=v->neib_vert[n]->neib_vert[m+1];
                                    v->neib_vert[n]->neib_vert[m+1]=v->neib_vert[n]->neib_vert[m+2];
                                }
                                if(m==2)
                                {
                                    v->neib_vert[n]->neib_vert[m]=v->neib_vert[n]->neib_vert[m+1];
                                }
                                v->neib_vert[n]->v_neigh_count=v->neib_vert[n]->v_neigh_count-1;
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        if(EV->next)
            delete_vertice(EV->next,v,type,debug);
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
            long double Sx,Sy;
            Sx=ve->neib_vert[i]->p->x-ve->p->x;
            Sy=ve->neib_vert[i]->p->y-ve->p->y;
            Sx=(Sx-(tilt*lroundl(Sy/twob)));
            Sx=(Sx-(twob*lroundl(Sx/twob)));
            Sy=(Sy-(twob*lroundl(Sy/twob)));
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
long double perimeter(long double Vx,long double Vy,long double E12x,long double E12y,long double Ax,long double Ay,long double Arad)
{
    E12x=E12x-Vx;
    E12y=E12y-Vy;
    Ax=Ax-Vx;
    Ay=Ay-Vy;
    long double areac,areat;
    long double perimeter;
    long double dis,disint,dise,disv;
    dis=sqrtl(distance((E12x-Ax),(E12y-Ay)));
    if(dis<Arad)
    {
        disint=sqrtl(Arad*Arad-dis*dis);
        dise=sqrtl(distance(E12x,E12y));
        disv=sqrtl(distance(Ax,Ay));
        long double Aintx,Ainty,X1,Y1,Aint1x,Aint1y,X2,Y2;
        Aintx=(dise-disint)/dise*E12x;
        Ainty=(dise-disint)/dise*E12y;
        Aint1x=(disv-Arad)/disv*Ax;
        Aint1y=(disv-Arad)/disv*Ay;
        X1=Aintx-Ax;
        Y1=Ainty-Ay;
        X2=Aint1x-Ax;
        Y2=Aint1y-Ay;
        long double theta=acos((X2*X1+Y1*Y2)/(Arad*Arad));
        if((X2*X1+Y1*Y2)/(Arad*Arad)>1.)
            theta=0.;
        perimeter=2.*M_PI*Arad*(theta/(2.*M_PI));
        return perimeter;
    }
    else
    {
        disv=sqrtl(distance(Ax,Ay));
        long double Aintx,Ainty,X1,Y1,Aint1x,Aint1y,X2,Y2;
        Aint1x=(disv-Arad)/disv*Ax;
        Aint1y=(disv-Arad)/disv*Ay;
        X1=E12x-Ax;
        Y1=E12y-Ay;
        X2=Aint1x-Ax;
        Y2=Aint1y-Ay;
        long double theta=acos((X2*X1+Y1*Y2)/(Arad*dis));
        if((X2*X1+Y1*Y2)/(Arad*dis)>1.)
            theta=0.;
        perimeter=2.*M_PI*Arad*(theta/(2.*M_PI));
        return perimeter;
    }

}
long double area_trangle(long double Vx,long double Vy,long double E12x,long double E12y,long double Ax,long double Ay,long double Arad,long double originx=0.,long double originy=0.)
{
    E12x=E12x-Vx;
    E12y=E12y-Vy;
    Ax=Ax-Vx;
    Ay=Ay-Vy;
    long double areac,areat;
    long double dis,disint,dise,disv;
    dis=sqrtl(distance((E12x-Ax),(E12y-Ay)));
    if(dis<Arad)
    {
        disint=sqrtl(Arad*Arad-dis*dis);
        dise=sqrtl(distance(E12x,E12y));
        disv=sqrtl(distance(Ax,Ay));
        long double Aintx,Ainty,X1,Y1,Aint1x,Aint1y,X2,Y2;
        Aintx=(dise-disint)/dise*E12x;
        Ainty=(dise-disint)/dise*E12y;
        Aint1x=(disv-Arad)/disv*Ax;
        Aint1y=(disv-Arad)/disv*Ay;
        X1=Aintx-Ax;
        Y1=Ainty-Ay;
        X2=Aint1x-Ax;
        Y2=Aint1y-Ay;
        long double a,b,c,p,q,l,m;
        a=sqrtl(X1*X1+Y1*Y1);
        b=sqrtl(X2*X2+Y2*Y2);
        c=sqrtl((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2));
        long double theta=2*asinl(sqrtl((c*c-(a-b)*(a-b))/(4*a*b)));
        areac=M_PIl*Arad*Arad*(theta/(2.*M_PIl));
        a=Ax;
        b=Ay;
        p=Aintx;
        q=Ainty;
        l=Aint1x;
        m=Aint1y;
        areat=0.5*abs(a*q-b*p);
        return areat-areac;
    }
    else
    {
        disv=sqrtl(distance(Ax,Ay));
        long double Aintx,Ainty,X1,Y1,Aint1x,Aint1y,X2,Y2;
        Aint1x=(disv-Arad)/disv*Ax;
        Aint1y=(disv-Arad)/disv*Ay;
        X1=E12x-Ax;
        Y1=E12y-Ay;
        X2=Aint1x-Ax;
        Y2=Aint1y-Ay;
        long double theta=acos((X2*X1+Y1*Y2)/(Arad*dis));
        if((X2*X1+Y1*Y2)/(Arad*dis)>1.)
        {
            theta=0.;
        }
        areac=M_PI*Arad*Arad*(theta/(2.*M_PI));
        long double a,b,p,q;
        a=Ax;
        b=Ay;
        p=E12x;
        q=E12y;
        areat=0.5*abs(a*q-b*p);
        return areat-areac;
    }

}
void update_neighbours(atom Atoms[],int nAtoms)
{
    long double R_CUT;
    //R_CUT=sqrtl(200./(4*3.14*density));
    R_CUT=powl(20./((4./3.)*3.14*density),1./3.);
    for(int i=0; i<nAtoms-1; i++)
    {
        for(int j=i+1; j<nAtoms; j++)
        {
            long double drx,dry,drz,dr;
            drx=Atoms[i].x-Atoms[j].x;
            dry=Atoms[i].y-Atoms[j].y;
            drz=Atoms[i].z-Atoms[j].z;
            drx=(drx-(tilt*lroundl(dry/twob)));
            drx=(drx-(twob*lroundl(drx/twob)));
            dry=(dry-(twob*lroundl(dry/twob)));
            drz=(drz-(twob*lroundl(drz/twob)));
            dr=drx*drx+dry*dry+drz*drz;
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
void print_delunay(atom *ATOM,delunay *D,atom Atoms[],int TYPE)
{
    cout<<"#\t";
	cout<<D->a<<"\t"<<D->b<<"\t"<<D->c<<"\n";
	cout<<"#\t";
	cout<<D->A<<"\t"<<D->B<<"\t"<<D->C<<"\n";
    long double Sx,Sy,Sz;
    long double Px,Py,Pz;
    Sx=ATOM->x-Atoms[D->a].x;
    Sy=ATOM->y-Atoms[D->a].y;
    Sz=ATOM->z-Atoms[D->a].z;
    Sx=(Sx-(tilt*lroundl(Sy/twob)));
    Sx=(Sx-(twob*lroundl(Sx/twob)));
    Sy=(Sy-(twob*lroundl(Sy/twob)));
    Sz=(Sz-(twob*lroundl(Sz/twob)));
    cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<ATOM->z<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\t"<<ATOM->z-Sz<<"\n";
    Sx=ATOM->x-Atoms[D->b].x;
    Sy=ATOM->y-Atoms[D->b].y;
    Sz=ATOM->z-Atoms[D->b].z;
    Sx=(Sx-(tilt*lroundl(Sy/twob)));
    Sx=(Sx-(twob*lroundl(Sx/twob)));
    Sy=(Sy-(twob*lroundl(Sy/twob)));
    Sz=(Sz-(twob*lroundl(Sz/twob)));
    cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<ATOM->z<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\t"<<ATOM->z-Sz<<"\n";
    Sx=ATOM->x-Atoms[D->c].x;
    Sy=ATOM->y-Atoms[D->c].y;
    Sz=ATOM->z-Atoms[D->c].z;
    Sx=(Sx-(tilt*lroundl(Sy/twob)));
    Sx=(Sx-(twob*lroundl(Sx/twob)));
    Sy=(Sy-(twob*lroundl(Sy/twob)));
    Sz=(Sz-(twob*lroundl(Sz/twob)));
    cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<ATOM->z<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\t"<<ATOM->z-Sz<<"\n";
    Sx=ATOM->x-Atoms[D->b].x;
    Sy=ATOM->y-Atoms[D->b].y;
    Sz=ATOM->z-Atoms[D->b].z;
    Sx=(Sx-(tilt*lroundl(Sy/twob)));
    Sx=(Sx-(twob*lroundl(Sx/twob)));
    Sy=(Sy-(twob*lroundl(Sy/twob)));
    Sz=(Sz-(twob*lroundl(Sz/twob)));
    Px=ATOM->x-Atoms[D->a].x;
    Py=ATOM->y-Atoms[D->a].y;
    Pz=ATOM->z-Atoms[D->a].z;
    Px=(Px-(tilt*lroundl(Py/twob)));
    Px=(Px-(twob*lroundl(Px/twob)));
    Py=(Py-(twob*lroundl(Py/twob)));
    Pz=(Pz-(twob*lroundl(Pz/twob)));
    cout<<ATOM->x-Px<<"\t"<<ATOM->y-Py<<"\t"<<ATOM->z-Pz<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\t"<<ATOM->z-Sz<<"\n";
    Sx=ATOM->x-Atoms[D->c].x;
    Sy=ATOM->y-Atoms[D->c].y;
    Sz=ATOM->z-Atoms[D->c].z;
    Sx=(Sx-(tilt*lroundl(Sy/twob)));
    Sx=(Sx-(twob*lroundl(Sx/twob)));
    Sy=(Sy-(twob*lroundl(Sy/twob)));
    Sz=(Sz-(twob*lroundl(Sz/twob)));
    Px=ATOM->x-Atoms[D->a].x;
    Py=ATOM->y-Atoms[D->a].y;
    Pz=ATOM->z-Atoms[D->a].z;
    Px=(Px-(tilt*lroundl(Py/twob)));
    Px=(Px-(twob*lroundl(Px/twob)));
    Py=(Py-(twob*lroundl(Py/twob)));
    Pz=(Pz-(twob*lroundl(Pz/twob)));
    cout<<ATOM->x-Px<<"\t"<<ATOM->y-Py<<"\t"<<ATOM->z-Pz<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\t"<<ATOM->z-Sz<<"\n";
    Sx=ATOM->x-Atoms[D->b].x;
    Sy=ATOM->y-Atoms[D->b].y;
    Sz=ATOM->z-Atoms[D->b].z;
    Sx=(Sx-(tilt*lroundl(Sy/twob)));
    Sx=(Sx-(twob*lroundl(Sx/twob)));
    Sy=(Sy-(twob*lroundl(Sy/twob)));
    Sz=(Sz-(twob*lroundl(Sz/twob)));
    Px=ATOM->x-Atoms[D->c].x;
    Py=ATOM->y-Atoms[D->c].y;
    Pz=ATOM->z-Atoms[D->c].z;
    Px=(Px-(tilt*lroundl(Py/twob)));
    Px=(Px-(twob*lroundl(Px/twob)));
    Py=(Py-(twob*lroundl(Py/twob)));
    Pz=(Pz-(twob*lroundl(Pz/twob)));
    cout<<ATOM->x-Px<<"\t"<<ATOM->y-Py<<"\t"<<ATOM->z-Pz<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\t"<<ATOM->z-Sz<<"\n";
}
void first_delunay(atom *ATOM,atom Atoms[],int TYPE)
{
    long double drx,dry,drz,dr;
    long double min=2.*twob*twob;
    int nearest,flag;
    int binv1=0;
    int binv2=0;
    long double midax=0.;
    long double miday=0.;
    long double midaz=0.;
    long double midbx=0.;
    long double midby=0.;
    long double midbz=0.;
    long double lmin=0;
    ATOM->conti[TYPE]=0;
    for(int i=0; i<ATOM->neighbours; i++)
    {
        long double X,Y,Z,x,y,z;
        long double rA,rS,DIS,l,dis_i,tan_sq;
        X=Atoms[ATOM->neighlist[i]].x-ATOM->x;
        Y=Atoms[ATOM->neighlist[i]].y-ATOM->y;
        Z=Atoms[ATOM->neighlist[i]].z-ATOM->z;
        X=(X-(tilt*lroundl(Y/twob)));
        X=(X-(twob*lroundl(X/twob)));
        Y=(Y-(twob*lroundl(Y/twob)));
        Z=(Z-(twob*lroundl(Z/twob)));
        DIS=sqrtl(X*X+Y*Y+Z*Z);
        rA=Atoms[ATOM->neighlist[i]].radius+r_cut;
        rS=ATOM->radius+r_cut;
        l=0.5*(DIS+(rS*rS-rA*rA)/DIS);
        x=l/DIS*X;
        y=l/DIS*Y;
        z=l/DIS*Z;
        dis_i=(x*x+y*y+z*z);
        tan_sq=dis_i-rS*rS;
	////cout<<"first\t";
	////cout<<i<<"\t"<<tan_sq<<"\n";
        if(tan_sq<min)
        {
            min=tan_sq;
            nearest=ATOM->neighlist[i];
            binv1=0;
            midax=x+ATOM->x;
            miday=y+ATOM->y;
            midaz=z+ATOM->z;
            lmin=l;
            if(DIS>rS+rA)
            {
                binv1=1;
            }
        }

    }
//	cout<<"nea\t"<<nearest<<"\n";
    ATOM->bondinvoid[ATOM->conti[TYPE]][TYPE]=binv1;
    ATOM->contigous[ATOM->conti[TYPE]][TYPE]=nearest;
    //ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
    ATOM->conti[TYPE]++;
  //if(lmin<0.)
  //{
  //    ATOM->ignore=1;
  //}
    long double DIS_MIN=box*box;
    int DIS_atom;
    long double dis;
    long double X,Y,Z;
    long double circx,circy,circz;
	long double norm_x,norm_y,norm_z;
    int temp;
    for(int i=0; i<ATOM->neighbours; i++)
    {
        if(ATOM->neighlist[i]!=nearest)
        {
            atom L=*ATOM;
            atom M=Atoms[ATOM->contigous[0][TYPE]];
            atom R=Atoms[ATOM->neighlist[i]];
            long double ax=M.x-L.x;
            long double ay=M.y-L.y;
            long double az=M.z-L.z;
            long double bx=R.x-L.x;
            long double by=R.y-L.y;
            long double bz=R.z-L.z;
            ax=(ax-(tilt*lroundl(ay/twob)));
            ax=(ax-(twob*lroundl(ax/twob)));
            ay=(ay-(twob*lroundl(ay/twob)));
            az=(az-(twob*lroundl(az/twob)));
            bx=(bx-(tilt*lroundl(by/twob)));
            bx=(bx-(twob*lroundl(bx/twob)));
            by=(by-(twob*lroundl(by/twob)));
            bz=(bz-(twob*lroundl(bz/twob)));
            long double rA,rS,rB;
            long double XA,YA,ZA,XB,YB,ZB;
            long double xA,yA,zA,xB,yB,zB;
            long double l;
            long double DISA;
            long double DISB;
            long double MA,MB,INMA,INMB;
            long double CA,CB;
            long double tan_sq;
            XA=ax;
            YA=ay;
            ZA=az;
            XB=bx;
            YB=by;
            ZB=bz;
            DISA=sqrtl(XA*XA+YA*YA+ZA*ZA);
            DISB=sqrtl(XB*XB+YB*YB+ZB*ZB);
          //MA=YA/XA;
          //MB=YB/XB;
          //INMA=-1./MA;
          //INMB=-1./MB;
            rA=M.radius+r_cut;
            rB=R.radius+r_cut;
            rS=L.radius+r_cut;
            l=0.5*(DISA+(rS*rS-rA*rA)/DISA);
            xA=l/DISA*XA;
            yA=l/DISA*YA;
            zA=l/DISA*ZA;
//			cout<<i<<"\n";
			//cout<<XA+L.x<<"\t"<<YA+L.y<<"\t"<<ZA+L.z<<"\n";
//			cout<<l<<"\t"<<rA<<"\t"<<DISA<<"\n";
//			cout<<xA+L.x<<"\t"<<yA+L.y<<"\t"<<zA+L.z<<"\n";
            l=0.5*(DISB+(rS*rS-rB*rB)/DISB);
//			cout<<l<<"\t"<<rB<<"\t"<<DISB<<"\n";
            xB=l/DISB*XB;
            yB=l/DISB*YB;
            zB=l/DISB*ZB;
			//cout<<XB+L.x<<"\t"<<YB+L.y<<"\t"<<ZB+L.z<<"\n";
//			cout<<xB+L.x<<"\t"<<yB+L.y<<"\t"<<zB+L.z<<"\n";
			long double a1,b1,c1,a2,b2,c2,a,b,c;
			a=zB*yA-zA*yB;
			b=zA*xB-xA*zB;
			c=yB*xA-yA*xB;
			a1=c*yA-zA*b;
			b1=zA*a-xA*c;
			c1=b*xA-yA*a;
			a2=c*yB-zB*b;
			b2=zB*a-xB*c;
			c2=b*xB-yB*a;
			long double t1,t2; 	
			t2=(b1*xA-a1*yA-b1*xB+a1*yB)/(b1*a2-a1*b2);
			t1=(xB-xA+a2*t2)/a1;
			X=xB+a2*t2;
			Y=yB+b2*t2;
			Z=zB+c2*t2;
          //CA=yA-INMA*xA;
          //CB=yB-INMB*xB;
          //X=(CB-CA)/(INMA-INMB);
          //Y=INMA*X+CA;
            tan_sq=X*X+Y*Y+Z*Z-rS*rS;
            //long double Y_AXIS=sqrtl(powl(X-xA,2)+powl(Y-yA,2)+powl(Z-zA,2));
          //int sign_C;
          //int sign_N;
          //long double m;
          //m=YA/XA;
          //if((by-(m*bx))<0.)
          //    sign_N=-1;
          //else
          //    sign_N=1;
          //if((Y-(m*X))<0.)
          //    sign_C=-1;
          //else
          //    sign_C=1;
          //X=X+L.x;
          //Y=Y+L.y;
          //if(sign_C!=sign_N)
          //{
          //    Y_AXIS=-1.*Y_AXIS;
          //}
		////cout<<"sec\t";
		////cout<<i<<"\t"<<tan_sq<<"\n";
            if(DIS_MIN > tan_sq)
            {
                DIS_MIN=tan_sq;
                DIS_atom=ATOM->neighlist[i];
                circx=X;
                circy=Y;
                circz=Z;
				norm_x=a;
				norm_y=b;
				norm_z=c;
                binv1=0;
                midbx=xB+L.x;
                midby=yB+L.y;
                midbz=zB+L.z;
                lmin=l;
                if(DISB>(rB+rS+2*r_cut))
                {
                    binv1=1;
                }
            }
        }

    }
//	cout<<"dis\t"<<DIS_atom<<"\n";
    ATOM->bondinvoid[ATOM->conti[TYPE]][TYPE]=binv1;
    ATOM->contigous[ATOM->conti[TYPE]][TYPE]=DIS_atom;
    //ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
    //ATOM->edge_index[ATOM->conti[TYPE]-1][TYPE]++;
    ATOM->conti[TYPE]++;
    DIS_MIN=box*box*box;
	long double circx_s=circx;
	long double circy_s=circy;
	long double circz_s=circz;
//	cout<<DIS_atom<<"\n";
    atom L=*ATOM;
    atom M=Atoms[ATOM->contigous[0][TYPE]];
    atom N=Atoms[ATOM->contigous[1][TYPE]];
////cout<<M.x<<"\t"<<M.y<<"\t"<<M.z<<"\n";
////cout<<N.x<<"\t"<<N.y<<"\t"<<N.z<<"\n";
////cout<<L.x<<"\t"<<L.y<<"\t"<<L.z<<"\n";
    for(int i=0; i<ATOM->neighbours; i++)
    {
        if(ATOM->neighlist[i]!=nearest && ATOM->neighlist[i]!=ATOM->contigous[1][TYPE])
        {
            atom R=Atoms[ATOM->neighlist[i]];
            long double ax=M.x-L.x;
            long double ay=M.y-L.y;
            long double az=M.z-L.z;
            long double bx=R.x-L.x;
            long double by=R.y-L.y;
            long double bz=R.z-L.z;
            long double cx=N.x-L.x;
            long double cy=N.y-L.y;
            long double cz=N.z-L.z;
            ax=(ax-(tilt*lroundl(ay/twob)));
            ax=(ax-(twob*lroundl(ax/twob)));
            ay=(ay-(twob*lroundl(ay/twob)));
            az=(az-(twob*lroundl(az/twob)));
            bx=(bx-(tilt*lroundl(by/twob)));
            bx=(bx-(twob*lroundl(bx/twob)));
            by=(by-(twob*lroundl(by/twob)));
            bz=(bz-(twob*lroundl(bz/twob)));
            cx=(cx-(tilt*lroundl(cy/twob)));
            cx=(cx-(twob*lroundl(cx/twob)));
            cy=(cy-(twob*lroundl(cy/twob)));
            cz=(cz-(twob*lroundl(cz/twob)));
            long double rA,rS,rB;
            long double XA,YA,ZA,XB,YB,ZB;
            long double xA,yA,zA,xB,yB,zB;
            long double l;
            long double DISA;
            long double DISB;
            long double MA,MB,INMA,INMB;
            long double CA,CB;
            long double tan_sq;
            XA=ax;
            YA=ay;
            ZA=az;
            XB=bx;
            YB=by;
            ZB=bz;
            DISA=sqrtl(XA*XA+YA*YA+ZA*ZA);
            DISB=sqrtl(XB*XB+YB*YB+ZB*ZB);
          //MA=YA/XA;
          //MB=YB/XB;
          //INMA=-1./MA;
          //INMB=-1./MB;
            rA=M.radius+r_cut;
            rB=R.radius+r_cut;
            rS=L.radius+r_cut;
            l=0.5*(DISA+(rS*rS-rA*rA)/DISA);
            xA=l/DISA*XA;
            yA=l/DISA*YA;
            zA=l/DISA*ZA;
            l=0.5*(DISB+(rS*rS-rB*rB)/DISB);
            xB=l/DISB*XB;
            yB=l/DISB*YB;
            zB=l/DISB*ZB;
			long double a1,b1,c1,a2,b2,c2,a,b,c;
			a=zB*yA-zA*yB;
			b=zA*xB-xA*zB;
			c=yB*xA-yA*xB;
			a1=c*yA-zA*b;
			b1=zA*a-xA*c;
			c1=b*xA-yA*a;
			a2=c*yB-zB*b;
			b2=zB*a-xB*c;
			c2=b*xB-yB*a;
			long double t1,t2,X1,Y1,Z1; 	
			t2=(b1*xA-a1*yA-b1*xB+a1*yB)/(b1*a2-a1*b2);
			t1=(xB-xA+a2*t2)/a1;
			X=xB+a2*t2;
			Y=yB+b2*t2;
			Z=zB+c2*t2;
			t2=(norm_y*circx_s-norm_x*circy_s-norm_y*X+norm_x*Y)/(norm_y*a-norm_x*b);
			t1=(X-circx_s+a*t2)/norm_x;
			X1=X+a*t2;
			Y1=Y+b*t2;
			Z1=Z+c*t2;
			X1=circx_s+norm_x*t1;
			Y1=circy_s+norm_y*t1;
			Z1=circz_s+norm_z*t1;
			//cout<<"ceneter\t";
			//cout<<X1<<"\t"<<Y1<<"\t"<<Z1<<"\n";
		////cout<<ATOM->neighlist[i]<<"\n";
		////cout<<ATOM->contigous[0][TYPE]<<"\t";
        ////cout<<ATOM->contigous[1][TYPE]<<"\n";
          //CA=yA-INMA*xA;
          //CB=yB-INMB*xB;
          //X=(CB-CA)/(INMA-INMB);
          //Y=INMA*X+CA;
            tan_sq=X1*X1+Y1*Y1+Z1*Z1-rS*rS;
			//cout<<"chec  \t";
			//cout<<(XA)*(X1-circx_s)+(YA)*(Y1-circy_s)+(ZA)*(Z1-circz_s)<<"\n";
            //long double Y_AXIS=sqrtl(powl(X-xA,2)+powl(Y-yA,2)+powl(Z-zA,2));
          //int sign_C;
          //int sign_N;
          //long double m;
          //m=YA/XA;
          //if((by-(m*bx))<0.)
          //    sign_N=-1;
          //else
          //    sign_N=1;
          //if((Y-(m*X))<0.)
          //    sign_C=-1;
          //else
          //    sign_C=1;
          //X1=X1+L.x;
          //Y1=Y1+L.y;
          //Z1=Z1+L.z;
          //X1=X1;
          //Y1=Y1;
          //Z1=Z1;
          //if(sign_C!=sign_N)
          //{
          //    Y_AXIS=-1.*Y_AXIS;
          //}
		////cout<<"third\t";
		////cout<<i<<"\t"<<tan_sq<<"\n";
            if(DIS_MIN > tan_sq)
            {
                DIS_MIN=tan_sq;
                DIS_atom=ATOM->neighlist[i];
                circx=X1+L.x;
                circy=Y1+L.y;
                circz=Z1+L.z;
                binv1=0;
                midbx=xB+L.x;
                midby=yB+L.y;
                midbz=zB+L.z;
                lmin=l;
                if(DISB>(rB+rS+2*r_cut))
                {
                    binv1=1;
                }
            }
        }

    }
////cout<<"this\t";
////cout<<(circx-circx_s-L.x)*(M.x-L.x)+(circy-circy_s-L.y)*(M.y-L.y)+(circz-circz_s-L.z)*(M.z-L.z)<<"\n";
////cout<<circx_s+L.x<<"\t"<<circy_s+L.y<<"\t"<<circz_s+L.z<<"\n";
    ATOM->bondinvoid[ATOM->conti[TYPE]][TYPE]=binv1;
    ATOM->contigous[ATOM->conti[TYPE]][TYPE]=DIS_atom;
    //ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
    ATOM->edge_index[0][TYPE]++;
    ATOM->edge_index[0][TYPE]++;
    ATOM->edge_index[1][TYPE]++;
    ATOM->edge_index[1][TYPE]++;
    ATOM->edge_index[2][TYPE]++;
    ATOM->edge_index[2][TYPE]++;
    ATOM->conti[TYPE]++;
    ATOM->D[TYPE].initial= new delunay;
    ATOM->D[TYPE].initial->A=0;
    ATOM->D[TYPE].initial->B=1;
    ATOM->D[TYPE].initial->C=2;
	ATOM->D[TYPE].initial->ABf=1;
	ATOM->D[TYPE].initial->BCf=1;
	ATOM->D[TYPE].initial->CAf=1;
	ATOM->part_c[0][1][TYPE]=1;
	ATOM->part_c[0][2][TYPE]=1;
	ATOM->part_c[1][2][TYPE]=1;
	ATOM->part_c[1][0][TYPE]=1;
	ATOM->part_c[2][1][TYPE]=1;
	ATOM->part_c[2][0][TYPE]=1;
    ATOM->D[TYPE].initial->a=ATOM->contigous[0][TYPE];
    ATOM->D[TYPE].initial->b=ATOM->contigous[1][TYPE];
    ATOM->D[TYPE].initial->c=ATOM->contigous[2][TYPE];
	//cout<<ATOM->D[TYPE].initial->a<<"\t"<<ATOM->D[TYPE].initial->b<<"\n";
    ATOM->D[TYPE].initial->circum_x=circx;
    ATOM->D[TYPE].initial->circum_y=circy;
    ATOM->D[TYPE].initial->circum_z=circz;
    ATOM->D[TYPE].initial->Ax=midax;
    ATOM->D[TYPE].initial->Ay=miday;
    ATOM->D[TYPE].initial->Az=midaz;
    ATOM->D[TYPE].initial->Bx=midbx;
    ATOM->D[TYPE].initial->By=midby;
    ATOM->D[TYPE].initial->Bz=midbz;
	//delunay *D;
}
void constr_del(atom *ATOM,atom Atoms[],int TYPE,long double p,long double q,long double r,long double Sx,long double Sy,long double Sz,int sign,int A,int B,int O,long double rA,long double rB,delunay *D,int debug=0)
{
//		cout<<"hereasdasdasn\n";
    long double a=ATOM->x;
    long double b=ATOM->y;
    long double c=ATOM->z;
	long double rS=ATOM->radius+r_cut;
	long double DIS,X,Y,Z;
	long double Y_MIN=box*box*box;
	long double circx,circy,circz;
	int DIS_atom;
//	cout<<p+a<<"\t"<<q+b<<"\t"<<r+c<<"\n";
///	cout<<Sx+a<<"\t"<<Sy+b<<"\t"<<Sz+c<<"\n";
//	cout<<a<<"\t"<<b<<"\t"<<c<<"\n";
	for(int j=0; j<ATOM->neighbours; j++)
	{
		if(debug)
				cout<<j<<"\t"<<ATOM->neighlist[j]<<"\n";
		int sign_N;
		long double x=Atoms[ATOM->neighlist[j]].x-a;
		long double y=Atoms[ATOM->neighlist[j]].y-b;
		long double z=Atoms[ATOM->neighlist[j]].z-c;
		long double rN=Atoms[ATOM->neighlist[j]].radius+r_cut;
		x=(x-(tilt*lroundl(y/twob)));
		x=(x-(twob*lroundl(x/twob)));
		y=(y-(twob*lroundl(y/twob)));
		z=(z-(twob*lroundl(z/twob)));
		long double ax=(Sy*r-Sz*q);
		long double ay=(Sz*p-Sx*r);
		long double az=(Sx*q-Sy*p);
		long double overlap=(x*ax+y*ay+z*az);	
		if(overlap<0.)
				sign_N=1;
		else
				sign_N=-1;
////////if((y-(m*x+C))<0.)
////////	sign_N=-1;
////////else
////////	sign_N=1;
		int flag=1;
		//to avoid atoms with have two delunay
		for(int k=0; k<ATOM->conti[TYPE]; k++)
		{
			if(ATOM->contigous[k][TYPE]==ATOM->neighlist[j])
			{
					//cout<<k<<"\n";
			        if(k==O ||  k==A || k==B )
			        {
			        		flag=0;
			        		break;
			        }

				int s=0;
				for(int p=0;p<ATOM->conti[TYPE];p++)
				{
						s=s+ATOM->part_c[k][p][TYPE];
				}
				if(s==2*ATOM->edge_index[k][TYPE])
				{
					////cout<<k<<"\n";
						//cout<<"see here\n";
					    //cout<<s<<"\t"<<ATOM->edge_index[k][TYPE]<<"\n";
						//cout<<"nohspprn\n";
						flag=0;
						break;
				}
			////if(ATOM->edge_index[k][TYPE]==4)
			////{
			////	flag=0;
			////	break;
			////}
			}
		}
		
		if(debug)
				cout<<sign<<"\t"<<sign_N<<"\t"<<flag<<"\n";
		if(sign!=sign_N && flag)
		{
				//cout<<"ATOM=\t";
		   // cout<<x+a<<"\t"<<y+b<<"\t"<<z+c<<"\n";
		////atom L=*ATOM;
		////atom M=Atoms[ATOM->contigous[D->A][TYPE]];
		////atom R=Atoms[ATOM->neighlist[j]];
		////long double ax=M.x-L.x;
		////long double ay=M.y-L.y;
		////long double bx=R.x-L.x;
		////long double by=R.y-L.y;
		////ax=(ax-(tilt*lroundl(ay/twob)));
		////ax=(ax-(twob*lroundl(ax/twob)));
		////ay=(ay-(twob*lroundl(ay/twob)));
		////bx=(bx-(tilt*lroundl(by/twob)));
		////bx=(bx-(twob*lroundl(bx/twob)));
		////by=(by-(twob*lroundl(by/twob)));
			long double XA,YA,ZA,XB,YB,ZB;
			long double xA,yA,zA,xB,yB,zB,xC,yC,zC;
			long double l;
			long double DISA;
			long double DISB;
			long double MA,MB,INMA,INMB;
			long double CA,CB;
			long double tan_sq;
			DISA=sqrtl(p*p+q*q+r*r);
			DISB=sqrtl(Sx*Sx+Sy*Sy+Sz*Sz);
			DIS=sqrtl(x*x+y*y+z*z);
		////rA=M.radius+r_cut;
		////rB=R.radius+r_cut;
		////rS=L.radius+r_cut;
			l=0.5*(DISA+(rS*rS-rA*rA)/DISA);
			//cout<<l<<"\t"<<rA<<"\t"<<DISA<<"\n";
			xA=l/DISA*p;
			yA=l/DISA*q;
			zA=l/DISA*r;
			//cout<<"ths\n";
		////cout<<p+a<<"\t"<<q+b<<"\t"<<r+c<<"\n";
		////cout<<Sx+a<<"\t"<<Sy+b<<"\t"<<Sz+c<<"\n";
			//cout<<xA+a<<"\t"<<yA+b<<"\t"<<zA+c<<"\n";
			l=0.5*(DIS+(rS*rS-rN*rN)/DIS);
			xB=l/DIS*x;
			yB=l/DIS*y;
			zB=l/DIS*z;
			//cout<<xB+a<<"\t"<<yB+b<<"\t"<<zB+c<<"\n";
			l=0.5*(DISB+(rS*rS-rB*rB)/DISB);
			//cout<<l<<"\t"<<rB<<"\t"<<DISB<<"\n";
			xC=l/DISB*Sx;
			yC=l/DISB*Sy;
			zC=l/DISB*Sz;
			//cout<<xC+a<<"\t"<<yC+b<<"\t"<<zC+c<<"\n";
			////cout<<"ths\n";
			long double v1x,v1y,v1z;
			long double v2x,v2y,v2z;
			v1x=yA*zB-zA*yB;
			v1y=zA*xB-xA*zB;
			v1z=xA*yB-yA*xB;
			v2x=yC*zA-zC*yA;
			v2y=zC*xA-xC*zA;
			v2z=xC*yA-yC*xA;
			long double a1,b1,c1;
			long double a2,b2,c2;
			long double a3,b3,c3;
			long double a4,b4,c4;
			a1=v1y*zA-v1z*yA;
			b1=v1z*xA-v1x*zA;
			c1=v1x*yA-v1y*xA;
			a2=v1y*zB-v1z*yB;
			b2=v1z*xB-v1x*zB;
			c2=v1x*yB-v1y*xB;
			a3=v2y*zC-v2z*yC;
			b3=v2z*xC-v2x*zC;
			c3=v2x*yC-v2y*xC;
			a4=v2y*zA-v2z*yA;
			b4=v2z*xA-v2x*zA;
			c4=v2x*yA-v2y*xA;
			long double t1,t2,X1,Y1,Z1,X2,Y2,Z2; 	
			t2=(b1*xA-a1*yA-b1*xB+a1*yB)/(b1*a2-a1*b2);
			X1=xB+a2*t2;
			Y1=yB+b2*t2;
			Z1=zB+c2*t2;
			t2=(b3*xC-a3*yC-b3*xA+a3*yA)/(b3*a4-a3*b4);
			////cout<<t2<<" wdfd\n";
			X2=xA+a4*t2;
			Y2=yA+b4*t2;
			Z2=zA+c4*t2;
		    //cout<<X1+a<<"\t"<<Y1+b<<"\t"<<Z1+c<<"\n";
		    //cout<<xC*v2x+yC*v2y+zC*v2z<<" inn\n";
		    //cout<<a1*v1x+b1*v1y+c1*v1z<<" inn\n";
		    //cout<<a2*v1x+b2*v1y+c2*v1z<<" inn\n";
			t2=(v1y*X1-v1x*Y1-v1y*X2+v1x*Y2)/(v1y*v2x-v1x*v2y);
			////cout<<t2<<"\n";
			X=X2+v2x*t2;
			Y=Y2+v2y*t2;
			Z=Z2+v2z*t2;
			//cout<<p<<"\t"<<q<<"\t"<<r<<"\n";
			//cout<<v2x*p+v2y*q+v2z*r<<" check this\n";
			////cout<<rS<<"\n";
		////CA=yA-INMA*xA;
		////CB=yB-INMB*xB;
		////X=(CB-CA)/(INMA-INMB);
		////Y=INMA*X+CA;
			////cout<<"centers\n";
		////cout<<"vec =";
		////cout<<X2+v2x*0.0+a<<"\t"<<Y2+v2y*0.0+b<<"\t"<<Z2+v2z*0.0+c<<"\t";
		////cout<<X2+v2x*-1.0+a<<"\t"<<Y2+v2y*-1.0+b<<"\t"<<Z2+v2z*-1.0+c<<"\n";
			//cout<<"Xchekc ";
			//cout<<X2+a<<"\t"<<Y2+b<<"\t"<<Z2+c<<"\n";
			tan_sq=X*X+Y*Y+Z*Z-rS*rS;
			//cout<<tan_sq<<"\n";
			long double Y_AXIS=sqrtl(powl(X-X2,2)+powl(Y-Y2,2)+powl(Z-Z2,2));
			int sign_C;
			long double overlap=(X*ax+Y*ay+Z*az);	
			if(overlap<0.)
					sign_C=1;
			else
					sign_C=-1;
			if(sign_C!=sign_N)
				Y_AXIS=-1.*Y_AXIS;
			//cout<<(X-X2)*p+(Y-Y2)*q+(Z-Z2)*r<<"\n";
			X=X+a;
			Y=Y+b;
			Z=Z+c;
			//cout<<X<<"\t"<<Y<<"\t"<<Z<<"\n";

////	//cout<<"trythis 		\t";
////	//cout<<x+a<<"\t"<<y+b<<"\t"<<z+c<<"\n";
////	//cout<<x+a<<"\t"<<y+b<<"\t"<<z+c<<"\n";
////		//cout<<"trythis\t";
			////cout<<"tried this one\n";
		////cout<<"YAIS=";
    	////cout<<ATOM->neighlist[j]<<"\t";
		////cout<<Y_AXIS<<"\t"<<tan_sq<<"\n";
			if(debug)
			{
				cout<<ATOM->neighlist[j]<<"\t"<<Y_AXIS<<"\n";
			}
			if(Y_AXIS<Y_MIN)
			{
				;
				Y_MIN=Y_AXIS;
				DIS_atom=ATOM->neighlist[j];
				circx=X;
				circy=Y;
				circz=Z;
			}
		}
	}//j loop
////cout<<"here\n";
////cout<<circx-D->circum_x<<"\t"<<circy-D->circum_y<<"\t"<<circz-D->circum_z<<"\n";
////cout<<D->circum_x<<"\t"<<D->circum_y<<"\t"<<D->circum_z<<"\n";
////cout<<p<<"\t"<<q<<"\t"<<r<<"\n";
////cout<<(circx-D->circum_x)*p+(circy-D->circum_y)*q+(circz-D->circum_z)*r<<"\n";
////cout<<"Xchekc end";
////cout<<"trythis\t"<<DIS_atom<<"\t"<<Y_MIN;
	//cout<<DIS_atom<<"\n";
	int flag=1;
////if(lmin<0.)
////{
////	ATOM->ignore=1;
////}
    int k;
	//cout<<ATOM->edge_index[12][TYPE]<<" twelve\n";
	for(k=0; k<ATOM->conti[TYPE]; k++)
	{
		if(ATOM->contigous[k][TYPE]==DIS_atom)
		{
			flag=0;
			break;
		}
	}
	if(debug)
	{
			cout<<"debug\t";
			cout<<flag<<"\t"<<k<<"\n";;
			cout<<DIS_atom<<"\n";
	}

	if(flag)
	{
		//ATOM->bondinvoid[ATOM->conti[TYPE]][TYPE]=binv;
	//		cout<<"hereag\n";
		ATOM->contigous[ATOM->conti[TYPE]][TYPE]=DIS_atom;
		ATOM->edge_index[A][TYPE]++;
		ATOM->edge_index[B][TYPE]++;
		ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
		ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
	    ATOM->part_c[ATOM->conti[TYPE]][A][TYPE]++;
	    ATOM->part_c[A][B][TYPE]++;
	    ATOM->part_c[B][ATOM->conti[TYPE]][TYPE]++;
	    ATOM->part_c[A][ATOM->conti[TYPE]][TYPE]++;
	    ATOM->part_c[B][A][TYPE]++;
	    ATOM->part_c[ATOM->conti[TYPE]][B][TYPE]++;
	}
	else
	{
			if(ATOM->part_c[A][k][TYPE]==0)
			{
				ATOM->edge_index[A][TYPE]++;
				ATOM->edge_index[k][TYPE]++;
			}
			if(ATOM->part_c[B][k][TYPE]==0)
			{
				ATOM->edge_index[B][TYPE]++;
				ATOM->edge_index[k][TYPE]++;
			}
	    ATOM->part_c[k][A][TYPE]++;
	    ATOM->part_c[A][B][TYPE]++;
	    ATOM->part_c[B][k][TYPE]++;
	    ATOM->part_c[A][k][TYPE]++;
	    ATOM->part_c[B][A][TYPE]++;
	    ATOM->part_c[k][B][TYPE]++;
		//ATOM->edge_index[k][TYPE]++;
	}
	//cout<<"look\t";
	//cout<<ATOM->edge_index[8][TYPE]<<"\t"<<A<<"\n";;
	//cout<<ATOM->part_c[12][10][TYPE]<<" track this\n";
	delunay *temp=nullptr;
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
	////temp->next->Ax=midbx;
	////temp->next->Ay=midby;
	}
	else
	{
		temp->next->A=k;
		temp->next->a=ATOM->contigous[k][TYPE];
	////temp->next->Ax=midbx;
	////temp->next->Ay=midby;
	}
	temp->next->B=A;
	temp->next->C=B;
	temp->next->b=ATOM->contigous[A][TYPE];
	temp->next->c=ATOM->contigous[B][TYPE];
	temp->next->Bx=D->Ax;
	temp->next->By=D->Ay;
	temp->next->circum_x=circx;
	temp->next->circum_y=circy;
	temp->next->circum_z=circz;
	if(flag)
		ATOM->conti[TYPE]++;
}
void complete_del(atom *ATOM,atom Atoms[],int nAtoms,int TYPE)
{
    long double Y_MIN=box*box;
    int DIS_atom;
    int binv;
    int binv2;
    long double lmin;
	long double a;
	long double b;
	long double c;
	long double p;
	long double q;
	long double r;
	long double Sx;
	long double Sy;
	long double Sz;
	long double Px;
	long double Py;
	long double Pz;
	int flag=1;
	long double rA,rB,rC;
    for(int i=0; i<ATOM->conti[TYPE]; i++)
    {
	    flag=1;
        Y_MIN=box*box;
        delunay *D;
		int s=0;
		//cout<<ATOM->conti[TYPE]<<"\n";
        D=ATOM->D[TYPE].initial;
		//cout<<"before\t";
		//cout<<ATOM->edge_index[3][TYPE]<<"\n";;
        while(flag)
        {
			//cout<<D->A<<"\t"<<D->B<<"\t"<<D->C<<"\n";
			//cout<<ATOM->edge_index[12][TYPE]<<" twelve\n";
            if(D->A==i)
            {
			    a=ATOM->x;
			    b=ATOM->y;
			    c=ATOM->z;
			    p=Atoms[ATOM->contigous[D->A][TYPE]].x-a;
			    q=Atoms[ATOM->contigous[D->A][TYPE]].y-b;
			    r=Atoms[ATOM->contigous[D->A][TYPE]].z-c;
				rA=Atoms[ATOM->contigous[D->A][TYPE]].radius+r_cut;			
			    Sx=Atoms[ATOM->contigous[D->B][TYPE]].x-a;
			    Sy=Atoms[ATOM->contigous[D->B][TYPE]].y-b;
			    Sz=Atoms[ATOM->contigous[D->B][TYPE]].z-c;
				rB=Atoms[ATOM->contigous[D->B][TYPE]].radius+r_cut;			
			    Px=Atoms[ATOM->contigous[D->C][TYPE]].x-a;
			    Py=Atoms[ATOM->contigous[D->C][TYPE]].y-b;
			    Pz=Atoms[ATOM->contigous[D->C][TYPE]].z-c;
			//	cout<<"otherside\n";
			//	cout<<"s\n";
				rC=Atoms[ATOM->contigous[D->C][TYPE]].radius+r_cut;			
			    p=(p-(tilt*lroundl(q/twob)));
			    p=(p-(twob*lroundl(p/twob)));
			    q=(q-(twob*lroundl(q/twob)));
			    r=(r-(twob*lroundl(r/twob)));
			    Sx=(Sx-(tilt*lroundl(Sy/twob)));
			    Sx=(Sx-(twob*lroundl(Sx/twob)));
			    Sy=(Sy-(twob*lroundl(Sy/twob)));
			    Sz=(Sz-(twob*lroundl(Sz/twob)));
			    Px=(Px-(tilt*lroundl(Py/twob)));
			    Px=(Px-(twob*lroundl(Px/twob)));
			    Py=(Py-(twob*lroundl(Py/twob)));
			    Pz=(Pz-(twob*lroundl(Pz/twob)));
			    int sign;
				if(ATOM->part_c[D->A][D->B][TYPE]==1)
				//if(D->ABf==1)
				{
						long double ax=(Sy*r-Sz*q);
						long double ay=(Sz*p-Sx*r);
						long double az=(Sx*q-Sy*p);
						long double overlap=(Px*ax+Py*ay+Pz*az);	
						if(overlap<0.)
								sign=1;
						else
								sign=-1;
					//	cout<<"D->A\t"<<"1"<<"\n";
						constr_del(ATOM,Atoms,TYPE,p,q,r,Sx,Sy,Sz,sign,D->A,D->B,D->C,rA,rB,D);
				}
				else if(ATOM->part_c[D->A][D->C][TYPE]==1)
				//if(D->CAf==1)
				{
						long double ax=(Py*r-Pz*q);
						long double ay=(Pz*p-Px*r);
						long double az=(Px*q-Py*p);
						long double overlap=(Sx*ax+Sy*ay+Sz*az);	
						if(overlap<0.)
								sign=1;
						else
								sign=-1;
					//	cout<<"D->A\t"<<"2"<<"\n";
						constr_del(ATOM,Atoms,TYPE,p,q,r,Px,Py,Pz,sign,D->A,D->C,D->B,rA,rC,D);
				}
			}
            if(D->B==i)
            {
					//cout<<"thisguy\n";
			    a=ATOM->x;
			    b=ATOM->y;
			    c=ATOM->z;
			    p=Atoms[ATOM->contigous[D->B][TYPE]].x-a;
			    q=Atoms[ATOM->contigous[D->B][TYPE]].y-b;
			    r=Atoms[ATOM->contigous[D->B][TYPE]].z-c;
				rA=Atoms[ATOM->contigous[D->B][TYPE]].radius+r_cut;			
			    Sx=Atoms[ATOM->contigous[D->A][TYPE]].x-a;
			    Sy=Atoms[ATOM->contigous[D->A][TYPE]].y-b;
			    Sz=Atoms[ATOM->contigous[D->A][TYPE]].z-c;
				rB=Atoms[ATOM->contigous[D->A][TYPE]].radius+r_cut;			
			    Px=Atoms[ATOM->contigous[D->C][TYPE]].x-a;
			    Py=Atoms[ATOM->contigous[D->C][TYPE]].y-b;
			    Pz=Atoms[ATOM->contigous[D->C][TYPE]].z-c;
				rC=Atoms[ATOM->contigous[D->C][TYPE]].radius+r_cut;			
			    p=(p-(tilt*lroundl(q/twob)));
			    p=(p-(twob*lroundl(p/twob)));
			    q=(q-(twob*lroundl(q/twob)));
			    r=(r-(twob*lroundl(r/twob)));
			    Sx=(Sx-(tilt*lroundl(Sy/twob)));
			    Sx=(Sx-(twob*lroundl(Sx/twob)));
			    Sy=(Sy-(twob*lroundl(Sy/twob)));
			    Sz=(Sz-(twob*lroundl(Sz/twob)));
			    Px=(Px-(tilt*lroundl(Py/twob)));
			    Px=(Px-(twob*lroundl(Px/twob)));
			    Py=(Py-(twob*lroundl(Py/twob)));
			    Pz=(Pz-(twob*lroundl(Pz/twob)));
			    int sign;
				if(ATOM->part_c[D->B][D->A][TYPE]==1)
				//if(D->ABf==1)
				{
						long double ax=(Sy*r-Sz*q);
						long double ay=(Sz*p-Sx*r);
						long double az=(Sx*q-Sy*p);
						long double overlap=(Px*ax+Py*ay+Pz*az);	
						if(overlap<0.)
								sign=1;
						else
								sign=-1;
						//cout<<"D->B\t"<<"1"<<"\n";
						constr_del(ATOM,Atoms,TYPE,p,q,r,Sx,Sy,Sz,sign,D->B,D->A,D->C,rA,rB,D,0);
				}
				else if(ATOM->part_c[D->B][D->C][TYPE]==1)
				//if(D->CAf==1)
				{
						long double ax=(Py*r-Pz*q);
						long double ay=(Pz*p-Px*r);
						long double az=(Px*q-Py*p);
					////cout<<Sx+a<<"\t"<<Sy+b<<"\t"<<Sz+c<<"\n";
					//////cout<<ax<<"\t"<<ay<<"\t"<<az<<"\n";
					////cout<<a<<"\t"<<b<<"\t"<<c<<"\t";
					////cout<<a+ax*0.2<<"\t"<<b+ay*0.2<<"\t"<<c+az*0.2<<"\n";
					////cout<<a<<"\t"<<b<<"\t"<<c<<"\t";
					////cout<<Px+a<<"\t"<<Py+b<<"\t"<<Pz+c<<"\n";
					////cout<<a<<"\t"<<b<<"\t"<<c<<"\t";
					////cout<<p+a<<"\t"<<q+b<<"\t"<<r+c<<"\n";
						long double overlap=(Sx*ax+Sy*ay+Sz*az);	
						if(overlap<0.)
								sign=1;
						else
								sign=-1;
						//cout<<"D->B\t"<<"2"<<"\n";
						constr_del(ATOM,Atoms,TYPE,p,q,r,Px,Py,Pz,sign,D->B,D->C,D->A,rA,rC,D);
				}
            }
			if(D->C==i)
			{
					//cout<<"here\t";
	//cout<<ATOM->part_c[8][10][TYPE]<<" track this\n";
			    a=ATOM->x;
			    b=ATOM->y;
			    c=ATOM->z;
			    p=Atoms[ATOM->contigous[D->C][TYPE]].x-a;
			    q=Atoms[ATOM->contigous[D->C][TYPE]].y-b;
			    r=Atoms[ATOM->contigous[D->C][TYPE]].z-c;
				rA=Atoms[ATOM->contigous[D->C][TYPE]].radius+r_cut;			
			    Sx=Atoms[ATOM->contigous[D->A][TYPE]].x-a;
			    Sy=Atoms[ATOM->contigous[D->A][TYPE]].y-b;
			    Sz=Atoms[ATOM->contigous[D->A][TYPE]].z-c;
				rB=Atoms[ATOM->contigous[D->A][TYPE]].radius+r_cut;			
			    Px=Atoms[ATOM->contigous[D->B][TYPE]].x-a;
			    Py=Atoms[ATOM->contigous[D->B][TYPE]].y-b;
			    Pz=Atoms[ATOM->contigous[D->B][TYPE]].z-c;
				rC=Atoms[ATOM->contigous[D->B][TYPE]].radius+r_cut;			
			    p=(p-(tilt*lroundl(q/twob)));
			    p=(p-(twob*lroundl(p/twob)));
			    q=(q-(twob*lroundl(q/twob)));
			    r=(r-(twob*lroundl(r/twob)));
			    Sx=(Sx-(tilt*lroundl(Sy/twob)));
			    Sx=(Sx-(twob*lroundl(Sx/twob)));
			    Sy=(Sy-(twob*lroundl(Sy/twob)));
			    Sz=(Sz-(twob*lroundl(Sz/twob)));
			    Px=(Px-(tilt*lroundl(Py/twob)));
			    Px=(Px-(twob*lroundl(Px/twob)));
			    Py=(Py-(twob*lroundl(Py/twob)));
			    Pz=(Pz-(twob*lroundl(Pz/twob)));
			    int sign;
				if(ATOM->part_c[D->C][D->A][TYPE]==1)
				//if(D->ABf==1)
				{
						long double ax=(Sy*r-Sz*q);
						long double ay=(Sz*p-Sx*r);
						long double az=(Sx*q-Sy*p);
						long double overlap=(Px*ax+Py*ay+Pz*az);	
						if(overlap<0.)
								sign=1;
						else
								sign=-1;
						//cout<<"here\n";
						constr_del(ATOM,Atoms,TYPE,p,q,r,Sx,Sy,Sz,sign,D->C,D->A,D->B,rA,rB,D,0);
				}
				else if(ATOM->part_c[D->C][D->B][TYPE]==1)
				//if(D->CAf==1)
				{
						long double ax=(Py*r-Pz*q);
						long double ay=(Pz*p-Px*r);
						long double az=(Px*q-Py*p);
						long double overlap=(Sx*ax+Sy*ay+Sz*az);	
						if(overlap<0.)
								sign=1;
						else
								sign=-1;
						constr_del(ATOM,Atoms,TYPE,p,q,r,Px,Py,Pz,sign,D->C,D->B,D->A,rA,rC,D);
				}

			}
            if(D->next)
			{
                D=D->next;
				////if(D->next)
				////{
				////		cout<<"!\t";
				////		cout<<D->next->A<<"\t"<<D->next->B<<"\t"<<D->next->C<<"\n";
					//}
			}
            else
            {
					s=0;
					//cout<<i<<"\t";
					for(int p=0;p<ATOM->conti[TYPE];p++)
					{
							s=s+ATOM->part_c[i][p][TYPE];
					}
					//cout<<s<<"`\t"<<ATOM->edge_index[i][TYPE]<<"\n";
					if(s==2*ATOM->edge_index[i][TYPE])
					{
					//		cout<<"somehting BDSBSDBSHBD\n";
						//	cout<<"here\n";
							flag=0;
					}
					break;
            }
        }
	}
}//end

int main( int argc , char * argv[] )
{
    int nAtoms=0;
    atom *Atoms=NULL;
    int counter=0;
    int config_count=0;
    //The file with configurations
    std::ifstream infile(argv[1]);///dat_trial");//config_2000_0.38_2_0.70.dat");
    //No of Atoms
    nAtoms=64;
    cout<<std::setprecision(5);
    //No of configurations in the input file
    config_count=1;
    //No of types of particle
    int ntypes=2;
    int SAM=0;
    long double *radius=new (nothrow) long double[ntypes];
    long double b,c,d,e,f;
    //radiuses of the particle
    radius[0]=0.2;
    radius[1]=0.4;
    //nAtoms=0;
    char buffer[64];
    vertice *temp_site=nullptr;
    temp_site=sites;
    long long free_dist[10000]= {0};
    long double max_free_area=0.;
    long double width;
    long double dummy;
  //snprintf(buffer,sizeof(char)*64,"free_dist");//_%d_%f.dat",int(nAtoms),Press);
  //ofstream fdist;
  //fdist.open(buffer);
    //This loop is over all the configurations
	ofstream fdist;
    for(int nconfig=0; nconfig<config_count; nconfig++)
    {
        //start[TYPE]is a list that stores all the voronoi vertices calculate with disks of radius radius[TYPE]
        start = new (nothrow) vertice*[ntypes];
        //CSTART is another list that stores all new vertices when you remove an atom and retessellate
        CSTART = new (nothrow) container_vertice*[ntypes];
        //The array of atoms
        Atoms = new (nothrow) atom[nAtoms];
        cout<<nconfig<<"\n"<<std::flush;
        counter=0;
        infile>>dummy;
        infile>>box;
        twob=2.0*box;
        //infile>>density;
        //infile>>dummy;
        //infile>>dummy;
        //infile>>tilt;
		if(nconfig==0)
		{
				cout<<tilt<<"\n";
			snprintf(buffer,sizeof(char)*64,"free_dist_%f",float(tilt));//_%d_%f.dat",int(nAtoms),Press);
			fdist.open(buffer);
		}
        while(infile>>b>>c>>d>>e>>f)
        {
            /* The loops puts the first nAtoms lines in the input file to list in the descending order*/
            counter++;
            if(sites==NULL)
            {
                sites=new vertice;
                sites->p=new site;
                sites->p->x=b;
                sites->p->y=c;
                sites->p->z=d;
                sites->r=e-epsilon;
				//display_SITE(sites->p);
            }
            else
            {
                temp_site=new vertice;
                temp_site->p=new site;
                temp_site->p->x=b;
                temp_site->p->y=c;
                temp_site->p->z=d;
                temp_site->r=e-epsilon;
                insert_site(sites,temp_site);
            }
            if(counter==nAtoms)
            {
                break;
            }
        }
		cout<<"brea\n";
        temp_site=sites;
        int cunt=0;
        while(1)
        {
            /* this loop puts the list of atoms into an array*/
            Atoms[cunt].x=temp_site->p->x;
            Atoms[cunt].y=temp_site->p->y;
            Atoms[cunt].z=temp_site->p->z;
            Atoms[cunt].radius=temp_site->r;
			//cout<<Atoms[cunt].x<<"\t"<<Atoms[cunt].y<<"\t"<<Atoms[cunt].z<<"\t"<<Atoms[cunt].radius<<"\n";;
            cunt++;
            if(temp_site->next)
                temp_site=temp_site->next;
            else
                break;
        }
		cout<<"here\n";
        //Make neighbour list for all atoms
        update_neighbours(Atoms,nAtoms);
        //this loops look for overlaps
        for(int i=0; i<nAtoms; i++)
        {
            for(int j=0; j<Atoms[i].neighbours; j++)
            {
                long double drx,dry,drz,dr;
                drx=Atoms[i].x-Atoms[Atoms[i].neighlist[j]].x;
                dry=Atoms[i].y-Atoms[Atoms[i].neighlist[j]].y;
                drz=Atoms[i].z-Atoms[Atoms[i].neighlist[j]].z;
                drx=(drx-(tilt*lroundl(dry/twob)));
                drx=(drx-(twob*lroundl(drx/twob)));
                dry=(dry-(twob*lroundl(dry/twob)));
                drz=(drz-(twob*lroundl(drz/twob)));
                dr=sqrtl(drx*drx+dry*dry+drz*drz);
                if(dr<Atoms[i].radius+Atoms[Atoms[i].neighlist[j]].radius)
                {
                    //cout<<"error\t"<<dr<<"\n";
					Atoms[i].ignore=1;
					break;
                }
                //cout<<Atoms[Atoms[i].neighlist[j]].x<<"\t"<<Atoms[Atoms[i].neighlist[j]].y<<"\n";
            }
        }
        long double area=0;
        vertice *save=nullptr;
        /* here all the arrays are initialized*/
        for(int i=0; i<nAtoms; i++)
        {
            for(int t=0; t<ntypes; t++)
            {
                if(Atoms[i].radius==radius[t]-epsilon)
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
				for(int j=0;j<100;j++)
				{
						{
								Atoms[i].part_c[t][j]=new (nothrow) int[ntypes];
						}
				}
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

        //This a loop over the types ,  if you have 2 kinds of atoms you have to do voronoi tessellation twice
        for(int TYPE=0; TYPE<ntypes; TYPE++)
        {
            snprintf(buffer,sizeof(char)*64,"vor_%d",int(TYPE));//_%d_%f.dat",int(nAtoms),Press);
            ofstream vor;
            vor.open(buffer);
            r_cut=radius[TYPE];
            cout<<TYPE<<"\t"<<r_cut<<"\n";
            //This is the loop over all atoms: we construct the voronoi cell for each atom
            for(SAM=0 ; SAM<nAtoms; SAM++)
            //for(SAM=0 ; SAM<nAtoms; SAM++)
            {
                {
                    //the first function calculates the first delunay triangle for the atom
                    first_delunay(&(Atoms[SAM]),Atoms,TYPE);
                    //this completes the delunay triangles of the atoms: all the triangle the atom takes part in
                    complete_del(&(Atoms[SAM]),Atoms,nAtoms,TYPE);
                    delunay *D=nullptr;
				    D=Atoms[SAM].D[TYPE].initial;
				  //while(1)
				  //{
				  //	print_delunay(&(Atoms[SAM]),D,Atoms,TYPE);
				  //	if(D->next)
				  //			D=D->next;
				  //	else
				  //			break;
				  //}
				////print_delunay(&(Atoms[SAM]),Atoms[SAM].D[TYPE].initial,Atoms,TYPE);
				////print_delunay(&(Atoms[SAM]),Atoms[SAM].D[TYPE].initial->next,Atoms,TYPE);
				////print_delunay(&(Atoms[SAM]),Atoms[SAM].D[TYPE].initial->next->next,Atoms,TYPE);
				////print_delunay(&(Atoms[SAM]),Atoms[SAM].D[TYPE].initial->next->next->next,Atoms,TYPE);
				////print_delunay(&(Atoms[SAM]),Atoms[SAM].D[TYPE].initial->next->next->next->next,Atoms,TYPE);
				////print_delunay(&(Atoms[SAM]),Atoms[SAM].D[TYPE].initial->next->next->next->next->next,Atoms,TYPE);
				////print_delunay(&(Atoms[SAM]),Atoms[SAM].D[TYPE].initial->next->next->next->next->next->next,Atoms,TYPE);
				////print_delunay(&(Atoms[SAM]),Atoms[SAM].D[TYPE].initial->next->next->next->next->next->next->next,Atoms,TYPE);
				////print_delunay(&(Atoms[SAM]),Atoms[SAM].D[TYPE].initial->next->next->next->next->next->next->next->next,Atoms,TYPE);
				////print_delunay(&(Atoms[SAM]),Atoms[SAM].D[TYPE].initial->next->next->next->next->next->next->next->next->next,Atoms,TYPE);
				////print_delunay(&(Atoms[SAM]),Atoms[SAM].D[TYPE].initial->next->next->next->next->next->next->next->next->next->next,Atoms,TYPE);
				////print_delunay(&(Atoms[SAM]),Atoms[SAM].D[TYPE].initial->next->next->next->next->next->next->next->next->next->next->next,Atoms,TYPE);
				////print_delunay(&(Atoms[SAM]),Atoms[SAM].D[TYPE].initial->next->next->next->next->next->next->next->next->next->next->next->next,Atoms,TYPE);
					//break;
                    int count=0;
                    long double area_s=0;
					cout<<"####\t####\t"<<SAM<<"\n";
                    for(int i=0; i<Atoms[SAM].conti[TYPE]; i++)
                    {
                        //We loop over all the atoms the atom "SAM" make a voronoi edge with
                        long double a,b,p,q,x,y,z;
                        x=Atoms[SAM].x;
                        y=Atoms[SAM].y;
                        z=Atoms[SAM].z;
                        int flaga=1;
                        int flagb=1;
                        delunay *D_ONE=NULL;
                        delunay *D_TWO=NULL;
						cout<<"#\t"<<i<<"\n";;
                        //D is a the first delunay triangle : Associted with each atom is a list of delunay triangle the atom is part of
                        //D_ONE and D_TWO are the two delunay triangle the atom SAM and its i'th contiguous atom takes part in (they define one edge
						for(int k=0;k<Atoms[SAM].conti[TYPE];k++)
						{
							if(Atoms[SAM].part_c[i][k][TYPE])
							{
                        		D=Atoms[SAM].D[TYPE].initial;
								D_ONE=NULL;
								D_TWO=NULL;
								while(1)
								{
									if(D->A==i)
									{
											if(D->B==k)
											{
													if(!D_ONE)
													{
														D_ONE=D;
													}
													else
													{
														D_TWO=D;
														break;
													}
											}
											else if(D->C==k)
											{
													if(!D_ONE)
													{
														D_ONE=D;
													}
													else
													{
														D_TWO=D;
														break;
													}

											}
									}
									else if(D->B==i)
									{
											if(D->A==k)
											{
													if(!D_ONE)
													{
														D_ONE=D;
													}
													else
													{
														D_TWO=D;
														break;
													}
											}
											else if(D->C==k)
											{
													if(!D_ONE)
													{
														D_ONE=D;
													}
													else
													{
														D_TWO=D;
														break;
													}

											}
									}
									else if(D->C==i)
									{
											if(D->A==k)
											{
													if(!D_ONE)
													{
														D_ONE=D;
													}
													else
													{
														D_TWO=D;
														break;
													}
											}
											else if(D->B==k)
											{
													if(!D_ONE)
													{
														D_ONE=D;
													}
													else
													{
														D_TWO=D;
														break;
													}

											}
									}
								////if(D->A==i && D->B==k)
								////{


								////}
									if(D->next)
											D=D->next;
									else
											break;
								}
								if(D_ONE && D_TWO)
								{
										cout<<"#\t"<<k<<"\n";
										cout<<D_ONE->circum_x<<"\t"<<D_ONE->circum_y<<"\t"<<D_ONE->circum_z<<"\t"<<D_TWO->circum_x<<"\t"<<D_TWO->circum_y<<"\t"<<D_TWO->circum_z<<"\n";
								}
								vertice *temp_vert_o=nullptr;
								vertice *temp_vert_d=nullptr;
								//THe D_ONE->circum_x/y are the co-ordinates of the voronoi vertice that is defined by the delunay triangle D_ONE
								//And we add them to the list of vertices / again in descending order
								if(!start[TYPE])
								{
									start[TYPE]=new vertice;
									start[TYPE]->p=new site;
									start[TYPE]->p->x=D_ONE->circum_x;
									start[TYPE]->p->y=D_ONE->circum_y;
									start[TYPE]->p->z=D_ONE->circum_z;
									start[TYPE]->A=SAM;
									start[TYPE]->D=D_ONE;
									temp_vert_o=start[TYPE];
								}
								else
								{
									vertice *temp=nullptr;
									temp=new vertice;
									temp->p=new site;
									temp->p->x=D_ONE->circum_x;
									temp->p->y=D_ONE->circum_y;
									temp->p->z=D_ONE->circum_z;
									temp->A=SAM;
									temp->D=D_ONE;
									temp=V->insert_vertice(start[TYPE],temp,TYPE);
									temp_vert_o=temp;
								}
								vertice *temp=nullptr;
								temp=new vertice;
								temp->p=new site;
								temp->p->x=D_TWO->circum_x;
								temp->p->y=D_TWO->circum_y;
								temp->p->z=D_TWO->circum_z;
								temp->A=SAM;
								temp->D=D_TWO;
								temp=V->insert_vertice(start[TYPE],temp,TYPE);
							}
						}
                      ////Now that we have the two vertices that define an edge we need to see if the bond between them lies in a void.
                      //long double m;
                      //long double X,Y,dis;
                      //X=Atoms[Atoms[SAM].contigous[i][TYPE]].x-Atoms[SAM].x;
                      //Y=Atoms[Atoms[SAM].contigous[i][TYPE]].y-Atoms[SAM].y;
                      //X=(X-(tilt*lroundl(Y/twob)));
                      //X=(X-(twob*lroundl(X/twob)));
                      //Y=(Y-(twob*lroundl(Y/twob)));
                      //dis=sqrtl(distance(X,Y));
                      //m=Y/X;
                      //int sign_C;
                      ////this part check if the two voronoi vertices lie on the same side of the line connecting the two atoms that make the voronoi edge
                      //if((D_ONE->circum_y-Atoms[SAM].y)-(m*(D_ONE->circum_x-Atoms[SAM].x))<0.)
                      //    sign_C=-1;
                      //else
                      //    sign_C=1;
                      //int sign_N;
                      //if((D_TWO->circum_y-Atoms[SAM].y)-(m*(D_TWO->circum_x-Atoms[SAM].x))<0.)
                      //    sign_N=-1;
                      //else
                      //    sign_N=1;
                      ////If the distance between atoms is greater than the sum of radiuses then the bond is in void
                      ////else it could be on the same side then it is considered to be in void (this will do no effect if the vertices are not in void
                      //if(dis>Atoms[Atoms[SAM].contigous[i][TYPE]].radius+Atoms[SAM].radius+2.*r_cut)
                      //    Atoms[SAM].bondinvoid[i][TYPE]=1;
                      //else if(sign_N == sign_C)
                      //    Atoms[SAM].bondinvoid[i][TYPE]=1;
                      //else
                      //    Atoms[SAM].bondinvoid[i][TYPE]=0;
                      ////write the edges to a file
                      //{
                      //    vor<<std::setprecision(15)<<x<<"\t"<<y<<"\t";
                      //    vor<<std::setprecision(15)<<D_ONE->circum_x<<"\t"<<D_ONE->circum_y<<"\t"<<D_TWO->circum_x<<"\t"<<D_TWO->circum_y<<"\n";
                      //    vor<<"\n";
                      //}
                      ////adding vertice
                      //temp_vert_d=temp;
                      ////connecting the two vertices ( the vertices connceted to a vertice is stores as a neibhourng vertice
                      //add_connected(temp_vert_o,temp_vert_d,Atoms[SAM].bondinvoid[i][TYPE]);
                      //add_connected(temp_vert_d,temp_vert_o,Atoms[SAM].bondinvoid[i][TYPE]);
                      //container_vertice *temp_cvert=nullptr;
                      ////for each atom we maintain a list of vertices that the atom takes part in using a container
                      //if(Atoms[SAM].Cstart[TYPE]==NULL)
                      //{
                      //    Atoms[SAM].Cstart[TYPE]=new container_vertice;
                      //    Atoms[SAM].Cstart[TYPE]->V=temp_vert_o;
                      //}
                      //else
                      //{
                      //    temp_cvert=new container_vertice;
                      //    temp_cvert->V=temp_vert_o;
                      //    insert_cvertice(Atoms[SAM].Cstart[TYPE],temp_cvert,Atoms[SAM].Cstart[TYPE]);
                      //}
                      //temp_cvert=new container_vertice;
                      //temp_cvert->V=temp_vert_d;
                      //insert_cvertice(Atoms[SAM].Cstart[TYPE],temp_cvert,Atoms[SAM].Cstart[TYPE]);

                    }

                }
            }
			return 0;
            vertice *temp_start=nullptr;
            temp_start=start[TYPE];
            int void_vert_count=0;
            cout<<"after  first tessellation \t"<<TYPE<<"\n";;
            cout<<r_cut<<"\n";
            //this part checks if the vertices are in void
            while(1)
            {
                if(temp_start->v_neigh_count != 3)
                {
                    //this should only happen rarely (if a vertice has more than 4 neighbours
                    display_SITE(temp_start->p);
                    V->delete_vertice(start[TYPE],temp_start,TYPE);
                    vertice *temp=nullptr;
                    temp=temp_start;
                    temp_start=temp_start->next;
                    delete temp->p;
                    delete temp;
                    continue;
                }
                int flag=1;
                long double AX,AY,BX,BY,X,Y;
                BX=temp_start->p->x ;
                BY=temp_start->p->y;
                AX=Atoms[temp_start->A].x;
                AY=Atoms[temp_start->A].y;
                long double dis=sqrtl((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
                if(dis<r_cut+Atoms[temp_start->A].radius)
                {
                    flag=0;
                }
                AX=Atoms[Atoms[temp_start->A].contigous[temp_start->D->A][TYPE]].x;
                AY=Atoms[Atoms[temp_start->A].contigous[temp_start->D->A][TYPE]].y;
                dis=sqrtl((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
                if(dis<r_cut+Atoms[Atoms[temp_start->A].contigous[temp_start->D->A][TYPE]].radius)
                {
                    flag=0;
                }
                AX=Atoms[Atoms[temp_start->A].contigous[temp_start->D->B][TYPE]].x;
                AY=Atoms[Atoms[temp_start->A].contigous[temp_start->D->B][TYPE]].y;
                dis=sqrtl((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
                if(dis<r_cut+Atoms[Atoms[temp_start->A].contigous[temp_start->D->B][TYPE]].radius)
                {
                    flag=0;
                }
                //flag=1 if the vertice is not in the exclusion area of the atoms
                if(flag )
                {
                    temp_start->is_void=1;
                    //mark it as void
                    void_vert_count=void_vert_count+1;
                }
                if(temp_start->next)
                    temp_start=temp_start->next;
                else
                    break;
            }
            vor.close();
        }
        //CODE BEGINS FOR CALCULATING FREE VOLUME
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
        int TYPE=1;
        //LOOP OVER ALL ATOMS , CALCULATE THE FREE VOL FOR EACH ATOM
        long double *freearea=nullptr;
        long double *freeperi=nullptr;
        freearea=new (nothrow) long double [nAtoms];
        freeperi=new (nothrow) long double [nAtoms];
        for(int i=0; i<nAtoms; i++)
        {
            freearea[i]=0.;
            freeperi[i]=0.;
        }
        int percentage=1;
        //loop over all the atoms to calculate the free volume
        for(int i=0; i<nAtoms; i++)
        {
			if(Atoms[i].ignore==0)
			{
            vertice *temp_start=nullptr;
            container_vertice *new_vert=NULL;
            int Atom_in_foc=i;
            if(i>percentage*nAtoms/4.)
            {
                //to show the progress
                cout<<percentage/4.*100.<<" completed\n";
                percentage++;

            }
			//cout<<i<<" #\n";
            TYPE=Atoms[i].type;
			//cout<<TYPE<<"\n";
			//cout<<Atoms[i].x<<"\t"<<Atoms[i].y<<"\n";
            //FIND THE ATOM RADIUS
            r_cut=radius[TYPE];
            container_vertice *ctemp=nullptr;
            ctemp=Atoms[i].Cstart[TYPE];
            //REMOVE THE VERTICES WHICH BELONGED TO THE VORONOI CELL OF THE ATOM IN CONSIDERATIO (i)
            ctemp=Atoms[i].Cstart[TYPE];
            while(1)
            {
                V->delete_vertice(start[TYPE],ctemp->V,TYPE);
                if(ctemp->next)
                    ctemp=ctemp->next;
                else
                    break;
            }
            snprintf(buffer,sizeof(char)*64,"vor_%d",int(TYPE));//_%d_%f.dat",int(nAtoms),Press);
            vor.open(buffer,std::ios_base::app);
            CSTART[TYPE]=NULL;
            //LOOP OVER ALL THE ATOMS THAT ARE CONTIGUOUS TO i
            delunay *D=nullptr;
            for(int j=0; j<Atoms[i].conti[TYPE]; j++)
            {
                int flag=0;
                //FIND THE ATOM INDEX OF THE CONTIGUOUS ATOM
                int SAM=Atoms[i].contigous[j][TYPE];
                container_vertice *ctemp=nullptr;
                //SAVE THE DETAILS OF THE ATOM 'SAM' BECAUSE THEY ARE GOING TO BE RETESSELLATED
                ctemp=Atoms[SAM].Cstart[TYPE];
                save_atom(&(Atoms[SAM]),TYPE);
                delunay *D=nullptr;
                D=Atoms[SAM].D[TYPE].initial;
                delunay *temp=nullptr;
                temp=Atoms[SAM].D[TYPE].initial;
                Atoms[SAM].save_D[TYPE].initial=temp;
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
                Atoms[SAM].neighbours=Atoms[SAM].neighbours-1;
                //FIND ALL THE DELUNAY TRIANGLES THIS 'SAM' TAKES PART IN AFTER REMOVING i
                first_delunay(&(Atoms[SAM]),Atoms,TYPE);
                complete_del(&(Atoms[SAM]),Atoms,nAtoms,TYPE);
                int count=0;
                long double area_s=0;
                for(int i=0; i<Atoms[SAM].conti[TYPE]; i++)
                {
                    D=Atoms[SAM].D[TYPE].initial;
                    long double a,b,p,q,x,y;
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
                    vertice *temp_vert_o=nullptr;
                    vertice *temp_vert_d=nullptr;
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
                        vertice *temp=nullptr;
                        vertice *temp_o=nullptr;
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
                            long double Ax,Ay,Bx,By,Cx,Cy,px,py,m,Co;
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
                            Cx=(Cx-(tilt*lroundl(Cy/twob)));
                            Cx=(Cx-(twob*lroundl(Cx/twob)));
                            Cy=(Cy-(twob*lroundl(Cy/twob)));
                            Bx=(Bx-(tilt*lroundl(By/twob)));
                            Bx=(Bx-(twob*lroundl(Bx/twob)));
                            By=(By-(twob*lroundl(By/twob)));
                            px=(px-(tilt*lroundl(py/twob)));
                            px=(px-(twob*lroundl(px/twob)));
                            py=(py-(twob*lroundl(py/twob)));
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
                                container_vertice *ctemp=nullptr;
                                ctemp = new container_vertice;
                                ctemp->V=temp;
                                if(!CSTART[TYPE])
                                    CSTART[TYPE]=ctemp;
                            }
                            container_vertice *ctemp=nullptr;
                            ctemp = new container_vertice;
                            ctemp->V=temp;
                            if(!new_vert)
                            {
                                new_vert=ctemp;
                            }
                            else
                            {
                                insert_cvertice(new_vert,ctemp,new_vert,0);
                            }
                        }
                        temp_vert_o=temp;
                    }
                    vertice *temp=nullptr;
                    vertice *temp_o=nullptr;
                    temp=new vertice;
                    temp->p=new site;
                    temp->p->x=D_TWO->circum_x;
                    temp->p->y=D_TWO->circum_y;
                    temp->A=SAM;
                    temp->D=D_TWO;
                    temp_o=temp;
                    temp=V->insert_vertice(start[TYPE],temp,TYPE);
                    if(temp_o==temp)
                    {
                        long double Ax,Ay,Bx,By,Cx,Cy,px,py,m,Co;
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
                        Bx=Bx-Ax;
                        By=By-Ay;
                        Cx=Cx-Ax;
                        Cy=Cy-Ay;
                        px=px-Ax;
                        py=py-Ay;
                        Cx=(Cx-(tilt*lroundl(Cy/twob)));
                        Cx=(Cx-(twob*lroundl(Cx/twob)));
                        Cy=(Cy-(twob*lroundl(Cy/twob)));
                        Bx=(Bx-(tilt*lroundl(By/twob)));
                        Bx=(Bx-(twob*lroundl(Bx/twob)));
                        By=(By-(twob*lroundl(By/twob)));
                        px=(px-(tilt*lroundl(py/twob)));
                        px=(px-(twob*lroundl(px/twob)));
                        py=(py-(twob*lroundl(py/twob)));
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
                            container_vertice *ctemp=nullptr;
                            ctemp = new container_vertice;
                            ctemp->V=temp;
                            if(!CSTART[TYPE])
                                CSTART[TYPE]=ctemp;
                        }
                        container_vertice *ctemp=nullptr;
                        ctemp = new container_vertice;
                        ctemp->V=temp;

                        if(!new_vert)
                        {
                            new_vert=ctemp;
                        }
                        else
                        {
                            insert_cvertice(new_vert,ctemp,new_vert,0);
                        }
                    }
                    temp_vert_d=temp;
                    long double m;
                    long double X,Y,dis;
                    X=Atoms[Atoms[SAM].contigous[i][TYPE]].x-Atoms[SAM].x;
                    Y=Atoms[Atoms[SAM].contigous[i][TYPE]].y-Atoms[SAM].y;
                    X=(X-(tilt*lroundl(Y/twob)));
                    X=(X-(twob*lroundl(X/twob)));
                    Y=(Y-(twob*lroundl(Y/twob)));
                    dis=sqrtl(distance(X,Y));
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
                    if(dis>Atoms[Atoms[SAM].contigous[i][TYPE]].radius+Atoms[SAM].radius+2.*r_cut)
                        Atoms[SAM].bondinvoid[i][TYPE]=1;
                    else if(sign_N == sign_C)
                        Atoms[SAM].bondinvoid[i][TYPE]=1;
                    else
                        Atoms[SAM].bondinvoid[i][TYPE]=0;
                  //{
                        vor<<std::setprecision(15)<<x<<"\t"<<y<<"\t"<<std::flush;
                        vor<<std::setprecision(15)<<D_ONE->circum_x<<"\t"<<D_ONE->circum_y<<"\t"<<D_TWO->circum_x<<"\t"<<D_TWO->circum_y<<"\n"<<std::flush;
                        vor<<"\n"<<std::flush;
                  //}
                    add_connected(temp_vert_o,temp_vert_d,Atoms[SAM].bondinvoid[i][TYPE],0);
                    add_connected(temp_vert_d,temp_vert_o,Atoms[SAM].bondinvoid[i][TYPE],0);
                }
            }
            temp_start=start[TYPE];
          //while(1)
          //{
          //    {
          //        int count=0;
          //        for(int n=0; n<temp_start->v_neigh_count; n++)
          //        {
          //            if(temp_start->neib_vert[n])
          //            {
          //                count++;
          //            }
          //        }
          //        if(temp_start->v_neigh_count!=3)
          //        {
          //            cout<<"#\t"<<count<<"= i need to \n";
          //            cout<<temp_start<<"\n";
          //            display_SITE(temp_start->p);
          //            for(int n=0; n<temp_start->v_neigh_count; n++)
          //            {
          //                if(temp_start->neib_vert[n])
          //                {
          //                    display_SITE(temp_start->neib_vert[n]->p);
          //                }
          //            }
          //        }
          //    }
          //    if(temp_start->next)
          //    {
          //        temp_start=temp_start->next;
          //    }
          //    else
          //        break;
          //}
			//cout<<"end\n";
            container_vertice *cstart=nullptr;
            cstart=CSTART[TYPE];
            container_vertice *temp_new_vert=nullptr;
            temp_start=start[TYPE];
            temp_new_vert=new_vert;
            temp_new_vert=new_vert;
            while(1)
            {
                if(!(compare(temp_new_vert->V->p,cstart->V->p)))
                {
                    if((temp_new_vert->V->p->y-cstart->V->p->y)||(temp_new_vert->V->p->x-cstart->V->p->x))
                    {
                        if(temp_new_vert->V->p->y>0.)
                        {
                            std::vector<int> EV1 {temp_new_vert->V->A,temp_new_vert->V->D->a,temp_new_vert->V->D->b};
                            std::vector<int> v1 {cstart->V->A,cstart->V->D->a,cstart->V->D->b};
                            std::sort(EV1.begin(),EV1.end());
                            std::sort(v1.begin(),v1.end());
                            if(EV1[0]==v1[0] && EV1[1]==v1[1] && EV1[2]==v1[2])
                            {
                                cstart->V=temp_new_vert->V;
                            }
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
            container_vertice *temp_new_vert1=nullptr;
            temp_new_vert=new_vert;
            int flag1=1,flag2=1;
            int rep;
            while(1)
            {
                temp_new_vert1=new_vert;
                flag1=1;
                while(1)
                {
                    flag2=1;
                    rep=compare(temp_new_vert->V->p,temp_new_vert1->V->p);
                    if(rep==0)
                    {
                        std::vector<int> EV1 {temp_new_vert->V->A,temp_new_vert->V->D->a,temp_new_vert->V->D->b};
                        std::vector<int> v1 {temp_new_vert1->V->A,temp_new_vert1->V->D->a,temp_new_vert1->V->D->b};
                        std::sort(EV1.begin(),EV1.end());
                        std::sort(v1.begin(),v1.end());
                        if(EV1[0]==v1[0] && EV1[1]==v1[1] && EV1[2]==v1[2])
                        {
                            rep=0;
                        }
                        else
                            rep=1;

                    }

                    if(rep==0)
                    {
                        if((temp_new_vert1->V->p->y-temp_new_vert->V->p->y)||(temp_new_vert1->V->p->x-temp_new_vert->V->p->x))
                        {
                            if(temp_new_vert1->V->p->y>0.)
                            {
                                for(int n=0; n<temp_new_vert->V->v_neigh_count; n++)
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
                                flag1=0;
                            }
                            else
                            {
                                for(int n=0; n<temp_new_vert1->V->v_neigh_count; n++)
                                {
                                    if(temp_new_vert1->V->neib_vert[n])
                                    {
                                        add_connected(temp_new_vert->V,temp_new_vert1->V->neib_vert[n],temp_new_vert1->V->neib_ed[n]);
                                        add_connected(temp_new_vert1->V->neib_vert[n],temp_new_vert->V,temp_new_vert1->V->neib_ed[n]);
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
                                container_vertice *temp=nullptr;
                                temp=temp_new_vert1;
                                V->delete_vertice(start[TYPE],temp_new_vert1->V,TYPE);
                                delete temp->V->p;
                                delete temp->V;
                                flag2=0;

                            }
                        }
                    }
                    if(temp_new_vert1->next)
                    {
                        if(flag2==0)
                        {
                            container_vertice *temp=nullptr;
                            temp=temp_new_vert1;
                            temp_new_vert1=temp_new_vert1->next;
                            delete temp;
                        }
                        else
                            temp_new_vert1=temp_new_vert1->next;
                    }
                    else
                    {
                        if(flag2==0)
                        {
                            container_vertice *temp=nullptr;
                            temp=temp_new_vert1;
                            delete temp;
                        }
                        break;
                    }
                }
                if(temp_new_vert->next)
                {
                    if(flag1==0)
                    {
                        container_vertice *temp=nullptr;
                        temp=temp_new_vert;
                        temp_new_vert=temp_new_vert->next;
                        delete temp;
                    }
                    else
                        temp_new_vert=temp_new_vert->next;
                }
                else
                {
                    if(flag1==0)
                    {
                        container_vertice *temp=nullptr;
                        temp=temp_new_vert;
                        delete temp;
                    }
                    break;
                }
            }
            temp_start=start[TYPE];
            while(1)
            {
                {
                    int count=0;
                    for(int n=0; n<temp_start->v_neigh_count; n++)
                    {
                        if(temp_start->neib_vert[n])
                        {
                            count++;
                        }
                    }
                    if(temp_start->v_neigh_count!=3)
                    {
                        cout<<"#\t"<<count<<"= i need to \n";
                        cout<<temp_start<<"\n";
                        display_SITE(temp_start->p);
                        for(int n=0; n<temp_start->v_neigh_count; n++)
                        {
                            if(temp_start->neib_vert[n])
                            {
                                display_SITE(temp_start->neib_vert[n]->p);
                            }
                        }
                    }
                }
                if(temp_start->next)
                {
                    temp_start=temp_start->next;
                }
                else
                    break;
            }
            temp_new_vert=new_vert;
            while(1)
            {
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
                temp_start=temp_new_vert->V;
                int flag=1;
                long double AX,AY,BX,BY,X,Y;
                BX=temp_start->p->x ;
                BY=temp_start->p->y;
                AX=Atoms[temp_start->A].x-BX;
                AY=Atoms[temp_start->A].y-BY;
                AX=AX-tilt*lroundl(AY/twob);
                AX=AX-twob*lroundl(AX/twob);
                AY=AY-twob*lroundl(AY/twob);
                long double dis=sqrtl((AX)*(AX)+(AY)*(AY));
                if(dis<r_cut+Atoms[temp_start->A].radius)
                {
                    flag=0;
                }
                AX=Atoms[Atoms[temp_start->A].contigous[temp_start->D->A][TYPE]].x-BX;
                AY=Atoms[Atoms[temp_start->A].contigous[temp_start->D->A][TYPE]].y-BY;
                AX=AX-tilt*lroundl(AY/twob);
                AX=AX-twob*lroundl(AX/twob);
                AY=AY-twob*lroundl(AY/twob);
                dis=sqrtl((AX)*(AX)+(AY)*(AY));
                if(dis<r_cut+Atoms[Atoms[temp_start->A].contigous[temp_start->D->A][TYPE]].radius)
                {
                    flag=0;
                }
                AX=Atoms[Atoms[temp_start->A].contigous[temp_start->D->B][TYPE]].x-BX;
                AY=Atoms[Atoms[temp_start->A].contigous[temp_start->D->B][TYPE]].y-BY;
                AX=AX-tilt*lroundl(AY/twob);
                AX=AX-twob*lroundl(AX/twob);
                AY=AY-twob*lroundl(AY/twob);
                dis=sqrtl((AX)*(AX)+(AY)*(AY));
                if(dis<r_cut+Atoms[Atoms[temp_start->A].contigous[temp_start->D->B][TYPE]].radius)
                {
                    flag=0;
                }
                AX=Atoms[i].x;
                AY=Atoms[i].y;
                BX=temp_start->p->x;
                BY=temp_start->p->y;
                X=AX-BX;
                Y=AY-BY;
                X=(X-(tilt*lroundl(Y/twob)));
                X=(X-(twob*lroundl(X/twob)));
                Y=(Y-(twob*lroundl(Y/twob)));
                dis=sqrtl((X)*(X)+(Y)*(Y));
                if(flag)//&& temp_start_start->v_neigh_count == 3)
                {
                    temp_start->is_void=1;
                }
                if(temp_new_vert->next)
                    temp_new_vert=temp_new_vert->next;
                else
                    break;
            }
            cstart=CSTART[TYPE];
            int void_vert_count=0;
            temp_new_vert=new_vert;
            cstart=CSTART[TYPE];
			int flag=0;
            while(1)
            {
                //if(cstart->V->is_void)
                {
                    void_vert_count++;
                }
                if(!cstart->V->is_void)
			    {
			        cout<<"this guy \n";
			    }
			////else
			////{
			////	flag=1;
			////	break;
			////}
                if(cstart->next)
                    cstart=cstart->next;
                else
                    break;
            }
		////if(flag)
		////{
		////		continue;
		////}
            int change=1;
            int void_vert_count_prev;
            int tem_ind=0;
            while(change)
            {
                cstart=CSTART[TYPE];
                void_vert_count_prev=void_vert_count;
                while(1)
                {
                    //if(cstart->V->is_void)
                    {
                        for(int i=0; i<cstart->V->v_neigh_count; i++)
                        {
                            if(cstart->V->neib_vert[i])
                            {
                                if(cstart->V->neib_vert[i]->is_void && cstart->V->neib_ed[i])
                                {
                                    int flag;
                                    container_vertice *ctemp=nullptr;
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
                if(void_vert_count_prev==void_vert_count)
                    change=0;
                else
                    change=1;
            }
            cstart=CSTART[TYPE];
            void_vert_count=0;
            while(1)
            {
                //if(cstart->V->is_void)
                {
                    void_vert_count++;
                }
                if(cstart->next)
                    cstart=cstart->next;
                else
                    break;
            }
			//cout<<"void count ="<<void_vert_count<<"\n";
            vertice **cavity_list=nullptr;
            cavity_list = new (nothrow) vertice*[void_vert_count];
            cstart=CSTART[TYPE];
            {
                int i=0;
                while(1)
                {
                    //if(cstart->V->is_void)
                    {
                        cavity_list[i]=cstart->V;
                        i++;
                    }
                    if(cstart->next)
                        cstart=cstart->next;
                    else
                        break;

                }

            }
            //CALCULATING THE VOID VOLUME IN A GIVEN CAVITY
            long double *void_area=nullptr;
            void_area= new (nothrow) long double[void_vert_count];
            long double *void_length=nullptr;
            void_length= new (nothrow) long double[void_vert_count];
            //return 0;
            for(int j=0 ; j<void_vert_count; j++)
            {
                void_area[j]=0;
                void_length[j]=0;
                long double A1x,A1y,A2x,A2y,A3x,A3y,Vx,Vy;
                long double E12x,E12y;
                long double E13x,E13y;
                long double E23x,E23y;
                long double A1r,A2r,A3r;
                long double m;
                int sign1,sign2;
                int S12,S23,S13;
                long double C;
                A1x=Atoms[cavity_list[j]->A].x;
                A1y=Atoms[cavity_list[j]->A].y;
                A1r=Atoms[cavity_list[j]->A].radius+r_cut;
                A2x=Atoms[cavity_list[j]->D->a].x-A1x;
                A2y=Atoms[cavity_list[j]->D->a].y-A1y;
                A2r=Atoms[cavity_list[j]->D->a].radius+r_cut;
                A3x=Atoms[cavity_list[j]->D->b].x-A1x;
                A3y=Atoms[cavity_list[j]->D->b].y-A1y;
                A3r=Atoms[cavity_list[j]->D->b].radius+r_cut;
                E12x=cavity_list[j]->D->Ax-A1x;
                E12y=cavity_list[j]->D->Ay-A1y;
                E13x=cavity_list[j]->D->Bx-A1x;
                E13y=cavity_list[j]->D->By-A1y;
                E12x=(E12x-(tilt*lroundl(E12y/twob)));
                E12x=(E12x-(twob*lroundl(E12x/twob)));
                E12y=(E12y-(twob*lroundl(E12y/twob)));
                E13x=(E13x-(tilt*lroundl(E13y/twob)));
                E13x=(E13x-(twob*lroundl(E13x/twob)));
                E13y=(E13y-(twob*lroundl(E13y/twob)));
                Vx=cavity_list[j]->p->x-A1x;
                Vy=cavity_list[j]->p->y-A1y;
                A2x=(A2x-(tilt*lroundl(A2y/twob)));
                A2x=(A2x-(twob*lroundl(A2x/twob)));
                A2y=(A2y-(twob*lroundl(A2y/twob)));
                A3x=(A3x-(tilt*lroundl(A3y/twob)));
                A3x=(A3x-(twob*lroundl(A3x/twob)));
                A3y=(A3y-(twob*lroundl(A3y/twob)));
                Vx=(Vx-(tilt*lroundl(Vy/twob)));
                Vx=(Vx-(twob*lroundl(Vx/twob)));
                Vy=(Vy-(twob*lroundl(Vy/twob)));
                long double X,Y,x,y;
                long double rA,rS,DIS,l,dis_i,tan_sq;
                X=Atoms[cavity_list[j]->D->a].x-Atoms[cavity_list[j]->D->b].x;
                Y=Atoms[cavity_list[j]->D->a].y-Atoms[cavity_list[j]->D->b].y;
                X=(X-(tilt*lroundl(Y/twob)));
                X=(X-(twob*lroundl(X/twob)));
                Y=(Y-(twob*lroundl(Y/twob)));
                DIS=sqrtl(X*X+Y*Y);
                rA=A2r;
                rS=A3r;
                l=0.5*(DIS+(rS*rS-rA*rA)/DIS);
                x=l/DIS*X;
                y=l/DIS*Y;
                x=x+Atoms[cavity_list[j]->D->b].x;
                y=y+Atoms[cavity_list[j]->D->b].y;
                E23x=x-A1x;
                E23y=y-A1y;
                E23x=(E23x-(tilt*lroundl(E23y/twob)));
                E23x=(E23x-(twob*lroundl(E23x/twob)));
                E23y=(E23y-(twob*lroundl(E23y/twob)));
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
                m=(A3y-A2y)/(A3x-A2x);
                C=A3y-m*A3x;
                if(-1.*C > 0. )
                {
                    sign1=1;
                }
                else
                    sign1=-1;
                if((Vy-m*Vx-C) > 0. )
                {
                    sign2=1;
                }
                else
                    sign2=-1;
                if(sign1!=sign2)
                {
                    S23=-1;
                }
                else
                    S23=1;
                void_area[j]=void_area[j]+S12*area_trangle(Vx,Vy,E12x,E12y,A2x,A2y,A2r,A1x,A1y);
                void_area[j]=void_area[j]+S12*area_trangle(Vx,Vy,E12x,E12y,0.,0.,A1r,A1x,A1y);
                void_area[j]=void_area[j]+S13*area_trangle(Vx,Vy,E13x,E13y,A3x,A3y,A3r,A1x,A1y);
                void_area[j]=void_area[j]+S13*area_trangle(Vx,Vy,E13x,E13y,0.,0.,A1r,A1x,A1y);
                void_area[j]=void_area[j]+S23*area_trangle(Vx,Vy,E23x,E23y,A3x,A3y,A3r,A1x,A1y);
                void_area[j]=void_area[j]+S23*area_trangle(Vx,Vy,E23x,E23y,A2x,A2y,A2r,A1x,A1y);
                void_length[j]=void_length[j]+S12*perimeter(Vx,Vy,E12x,E12y,A2x,A2y,A2r);
                void_length[j]=void_length[j]+S12*perimeter(Vx,Vy,E12x,E12y,0.,0.,A1r);
                void_length[j]=void_length[j]+S13*perimeter(Vx,Vy,E13x,E13y,A3x,A3y,A3r);
                void_length[j]=void_length[j]+S13*perimeter(Vx,Vy,E13x,E13y,0.,0.,A1r);
                void_length[j]=void_length[j]+S23*perimeter(Vx,Vy,E23x,E23y,A3x,A3y,A3r);
                void_length[j]=void_length[j]+S23*perimeter(Vx,Vy,E23x,E23y,A2x,A2y,A2r);
            }
            long double cav_tot=0.;
            long double ca_per_tot=0.;
            for(int i=0; i<void_vert_count; i++)
            {
                cav_tot=cav_tot+void_area[i];
                ca_per_tot=ca_per_tot+void_length[i];
            }
			//cout<<cav_tot<<"\n";
            if(CSTART[TYPE]->V->is_void)
            {
                freearea[i]=cav_tot;
            }
            else
            {
                freearea[i]=0.;
            }
            if(CSTART[TYPE]->V->is_void)
            {
                freeperi[i]=ca_per_tot;
            }
            else
            {
                freeperi[i]=0.;
            }
            temp_new_vert=new_vert;
            while(1)
            {
                V->delete_vertice(start[TYPE],temp_new_vert->V,TYPE);
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
                if(ctemp->next)
                    ctemp=ctemp->next;
                else
                    break;
            }
            for(int j=0; j<Atoms[i].conti[TYPE]; j++)
            {
                D=Atoms[i].D[TYPE].initial;
                delunay *D_ONE=NULL;
                delunay *D_TWO=NULL;

                while(1)
                {
                    if(D->A==j)
                    {
                        if(!D_ONE)
                        {
                            D_ONE=D;
                        }
                        else if(!D_TWO)
                            D_TWO=D;
                    }
                    if(D->B==j)
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
                vertice *temp1=nullptr;
                vertice *temp2=nullptr;
                temp1=new vertice;
                temp2=new vertice;
                temp1->p=new site;
                temp2->p=new site;
                temp1->p->x=D_ONE->circum_x;
                temp1->p->y=D_ONE->circum_y;
                temp2->p->x=D_TWO->circum_x;
                temp2->p->y=D_TWO->circum_y;
                temp1->A=i;
                temp2->A=i;
                temp1->D=D_ONE;
                temp2->D=D_TWO;
                temp1=V->insert_vertice(start[TYPE],temp1,TYPE,0);
                temp2=V->insert_vertice(start[TYPE],temp2,TYPE,0);
                add_connected(temp1,temp2,Atoms[i].bondinvoid[j][TYPE],0);
                add_connected(temp2,temp1,Atoms[i].bondinvoid[j][TYPE],0);
            }
            ctemp=Atoms[i].Cstart[TYPE];
            while(1)
            {
                for(int n=0; n<ctemp->V->v_neigh_count; n++)
                {
                    if(ctemp->V->neib_vert[n])
                        add_connected(ctemp->V->neib_vert[n],ctemp->V,ctemp->V->neib_ed[n],0);
                }
                if(ctemp->next)
                    ctemp=ctemp->next;
                else
                    break;
            }

            temp_start=start[TYPE];
            while(1)
            {
                {
                    int count=0;
                    for(int n=0; n<temp_start->v_neigh_count; n++)
                    {
                        if(temp_start->neib_vert[n])
                        {
                            count++;
                        }
                    }
                    if(temp_start->v_neigh_count!=3)
                    {
                        cout<<count<<"= i need to \n";
                        cout<<temp_start->v_neigh_count<<"\n";
                        display_SITE(temp_start->p);
                    }
                }
                if(temp_start->next)
                {
                    temp_start=temp_start->next;
                }
                else
                    break;
            }

            for(int j=0; j<Atoms[i].conti[TYPE]; j++)
            {
                int flag=0;
                int SAM=Atoms[i].contigous[j][TYPE];
                delunay *D=nullptr;
                D=Atoms[SAM].D[TYPE].initial;
                while(1)
                {
                    delunay *temp=nullptr;
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
                    container_vertice *ctemp=nullptr;
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
                    container_vertice *ctemp=nullptr;
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
        }
        delete [] CSTART;
        for(int i=0; i<nAtoms; i++)
        {
            fdist<<i<<"\t"<<freearea[i]<<"\n"<<std::flush;
        }
        long double *sum_freearea=new (nothrow) long double [ntypes];
        long double *sum_freeperi=new (nothrow) long double [ntypes];
        long double *sum_freeratio=new (nothrow) long double [ntypes];
        int *count=new (nothrow) int [ntypes];
        for(int t=0; t<ntypes; t++)
        {
            count[t]=0;
            sum_freearea[t]=0.;
            sum_freeperi[t]=0.;
            sum_freeratio[t]=0.;
        }
        for(int i=0; i<nAtoms; i++)
        {
            sum_freearea[Atoms[i].type]=sum_freearea[Atoms[i].type]+freearea[i];
            sum_freeperi[Atoms[i].type]=sum_freeperi[Atoms[i].type]+freeperi[i];
            if(freearea[i]!=0.0)
            {
                sum_freeratio[Atoms[i].type]=sum_freeratio[Atoms[i].type]+freeperi[i]/freearea[i];
                count[Atoms[i].type]++;
            }
        }
        delete[] sum_freearea;
        delete[] sum_freeperi;
        delete[] sum_freeratio;
        delete[] count;
        cout<<"\n";

        delete[] freearea;
        delete[] freeperi;
        for(int i=0; i<nAtoms; i++)
        {
            delete[] Atoms[i].conti;
            delete[] Atoms[i].save_conti;
            for(int t=0; t<500; t++)
            {
                delete [] Atoms[i].contigous[t];
                delete [] Atoms[i].edge_index[t];
                delete [] Atoms[i].bondinvoid[t];
                delete [] Atoms[i].save_contigous[t];
                delete [] Atoms[i].save_edge_index[t];
                delete [] Atoms[i].save_bondinvoid[t];
            }
        }
        for(int t=0; t<ntypes; t++)
        {
            container_vertice *cstart=nullptr;
            vertice *temp_start=nullptr;
            temp_start=start[t];
            while(1)
            {
                vertice *temp=nullptr;
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
                    container_vertice *temp_c=nullptr;
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
        delete [] start;
        for(int n=0; n<nAtoms; n++)
        {
            delete[] Atoms[n].Cstart;
            delunay *D=nullptr;
            for(int t=0; t<ntypes; t++)
            {
                D=Atoms[n].D[t].initial;
                while(1)
                {
                    delunay *temp=nullptr;
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
            delete [] Atoms[n].D;
            delete [] Atoms[n].save_D;
        }
        vertice	*temp_start=sites;
        while(1)
        {
            vertice *temp=nullptr;
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
        sites=NULL;
        delete[] Atoms;
    }
    delete[] radius;

    return 0;
}
