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
long double LAYER_CUT;
long double tilt=0.0;
long double DMIN=0.000000000001;//std::numeric_limits<long double>::min();
long double epsilon=0.00000000000;
int void_vert_count=0;
int **contigous;
int PBC=0;
long double *radius;
long double convex_vol=0.;
int nAtoms=0;
struct atom;
struct face;
struct vertice;
template <class T>
struct container
{
    T *t=nullptr;
    container <T> *next=nullptr;
    container <T> *prev=nullptr;
};
template <class T>
class set_of_container
{
public:
    container <T> *initial=nullptr;
    container <T> *end=nullptr;
    void insert_element(T *t);
    void remove_element(container <T> *t);
};
template < class T >
void set_of_container < T > :: insert_element(T *t)
{
    container <T> *temp=nullptr;
    temp = new container <T>;
    temp->t=t;
    if(!initial)
    {
        initial=temp;
        end=temp;
    }
    else
    {
        end->next=temp;
        temp->prev=end;
        end=end->next;
    }
}
struct site
{
    long double x=0;
    long double y=0;
    long double z=0;
};
long double distancesq(site p,site q)
{
    long double dist=0.;
    dist=(p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y)+(p.z-q.z)*(p.z-q.z);
    return dist;
}
struct vect
{
    long double x=0;
    long double y=0;
    long double z=0;
};
void display_SITE(struct site *p)
{

    cout<<"draw sphere\t{";
    cout<<p->x<<"\t"<<p->y<<"\t"<<p->z<<"}\tradius\t0.1\n";
}
//definition of half-edge
struct edge
{
    vertice *origin=nullptr;
    vertice *destination=nullptr;
    int A1,A2,A3;
    int evoid=0;
    edge *next=nullptr;
    edge *prev=nullptr;
};
//definition of vertice
//definition of face
struct face
{
    site A1,A2,A3;
    site B;
    face *next=nullptr;
};

struct delunay
{
    int AT[4]= {0};
    site MID[4][4][4];
    site MIDP[4][4];
    int EDGE[4][4];
    int FACE[4][4][4];
    int DEL=0;
    vertice *v=nullptr;
    long double circum_x=0.;
    long double circum_y=0.;
    long double circum_z=0.;
    delunay *next=nullptr;
    delunay *prev=nullptr;
    int hull=0;
    //int solid=0;
};
struct container_delunay
{
    delunay *D=nullptr;
    struct container_delunay *next=nullptr;
    struct container_delunay *prev=nullptr;
};
struct vertice
{
    struct site *p=nullptr;
    struct half_edge *leaving=nullptr;
    struct vertice *next=nullptr;
    struct vertice *prev=nullptr;
    int dangling=0;
    int visited=0;
    delunay *D=nullptr;
    int is_void=0;
    int cluster_index=-1;
    long double r;
    long double ball_r=0.;
    //vert_list *V;
    vertice *neib_vert[10]= {nullptr};
    int neib_ed[10];
    int v_neigh_count=0;
}*start,*s_temp,*sites;
struct container_vertice
{
    struct vertice *V=nullptr;
    struct container_vertice *next=nullptr;
    struct container_vertice *prev=nullptr;
}**CSTART,*VOID_START;
int compare(struct site *p1,struct site *p2,int del=0)
{
    long double DX,DY,DZ;
    DX=p1->x-p2->x;
    DY=p1->y-p2->y;
    DZ=p1->z-p2->z;
    if(!del)
    {
        DX=(DX-(tilt*PBC*PBC*lroundl(DY/twob)));
        DX=(DX-(twob*PBC*PBC*lroundl(DX/twob)));
        DY=(DY-(twob*PBC*lroundl(DY/twob)));
        DZ=(DZ-(twob*PBC*lroundl(DZ/twob)));
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
        {
            if(!compare(focus->neib_vert[i]->p,add->p))
            {
                //std::vector<int> EV1 {focus->neib_vert[i]->A,focus->neib_vert[i]->D->a,focus->neib_vert[i]->D->b,focus->neib_vert[i]->D->c};
                //std::vector<int> v1 {add->A,add->D->a,add->D->b,add->D->c};
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
    delunay *initial=nullptr;
    delunay *end=nullptr;
    void insert_delunay(delunay *d1);
    void remove_delunay(delunay *d1);
} FULLSETD;
void set_of_delunay :: insert_delunay(delunay *d1)
{
    if(!initial)
    {
        initial=d1;
        end=d1;
    }
    else
    {
        end->next=d1;
        d1->prev=end;
        end=end->next;
    }
}
void set_of_delunay :: remove_delunay(delunay *d1)
{
    if(d1==initial)
    {
        initial=d1->next;
        initial->prev=nullptr;
    }
    else if (d1==end)
    {
        end=d1->prev;
        end->next=nullptr;
    }
    else
    {
        d1->next->prev=d1->prev;
        d1->prev->next=d1->next;
    }

}
class set_of_cvert
{
public :
    container_vertice *initial=nullptr;
    container_vertice *end=nullptr;
    void insert_cvert(container_vertice *d1);
    void delete_cvert(container_vertice *d1);
};

void set_of_cvert:: insert_cvert(container_vertice *d1)
{
    if(!initial)
    {
        initial=d1;
        end=d1;
    }
    else
    {
        end->next=d1;
        end=end->next;
    }
}
void set_of_cvert::delete_cvert(container_vertice *d1)
{
    container_vertice *temp;
    temp=initial;
    while(temp)
    {
        if(temp->V->D->AT[0]==d1->V->D->AT[0] && temp->V->D->AT[1]==d1->V->D->AT[1] && temp->V->D->AT[2]==d1->V->D->AT[2] && temp->V->D->AT[3]==d1->V->D->AT[3] )
        {
            if(temp==initial)
            {
                initial=temp->next;
                delete temp;
                return ;
            }
            else if(temp==end)
            {
                end=temp->prev;
                delete temp;
                return ;
            }
            else
            {
                temp->prev->next=temp->next;
                temp->next->prev=temp->prev;
                delete temp;
                return ;
            }
        }
        temp->next;
    }
    cout<<"not found\n";
}
class convex_hull
{
public :
    face *initial=nullptr;
    face *end=nullptr;
    void insert_face(face *f);
} CH,solid_wall;

void convex_hull :: insert_face(face *f1)
{
    if(!initial)
    {
        initial=f1;
        end=f1;
    }
    else
    {
        end->next=f1;
        end=end->next;
    }
}
class edge_list
{
public:
    edge *initial=nullptr;
    edge *end=nullptr;
    void insert_edge(edge *e);
};
void edge_list::insert_edge(edge *e)
{
    if(!initial)
    {
        initial=e;
        end=e;
    }
    else
    {
        end->next=e;
        e->prev=end;
        end=e;
    }
}
struct atom
{
    site p;
    int neighlist[50000];
    int contigous[50]= {-1};
    set_of_container < edge > F[50];
    int *conti_index;
    int part_c[50][50]= {0};
    int bond_c[50][50]= {0};
    int bondinvoid[50]= {0};
    int D3bondinvoid[50][50]= {0};
    int edge_index[50]= {0};
    int neighbours=0;
    int conti= {0};
    int bor=0;
    site MIDP[50][50];
    site RMID[50];
    long double radius=1.;
    int ignore=0;
    int index;
    long double vor_vol;
    container_vertice *Cstart;
    set_of_delunay *D;
    set_of_cvert *S;
    container_delunay *D_FIRST;
    int type;
};
class vert_list
{
public:
    vertice *initial=nullptr;
    vertice *end=nullptr;
    void delete_vertice(vertice *,int type,int);
    void insert_vertice(vertice *,int,int);
    void display_vertex(vertice *);
    void display_conn(vertice *ve);
}*V;
void vert_list::insert_vertice(vertice *v,int type,int debug=0)
{
    if(!initial)
    {
        initial=v;
        end=v;
    }
    else
    {
        end->next=v;
        v->prev=end;
        end=v;
    }
}
void vert_list::delete_vertice(vertice *v,int type,int debug=0)
{
    if(v->next)
    {
        if(v->prev)
        {
            v->prev->next=v->next;
            v->next->prev=v->prev;
        }
        else
        {
            initial=v->next;
            v->next->prev=nullptr;
        }
    }
    else
    {
        end=v->prev;
        v->prev->next=nullptr;
    }
    vertice *v_neib=nullptr;
    for(int i=0; i<v->v_neigh_count; i++)
    {
        v_neib=v->neib_vert[i];
        int temp_index=0;
        for(int j=0; j<v_neib->v_neigh_count; j++)
        {
            if(v_neib->neib_vert[j]==v)
            {
                v_neib->neib_vert[j]=nullptr;
            }
        }
        for(int j=0; j<v_neib->v_neigh_count; j++)
        {
            if(v_neib->neib_vert[j]==nullptr)
            {
                if(v_neib->neib_vert[j+1])
                {
                    v_neib->neib_vert[j]=v_neib->neib_vert[j+1];
                    v_neib->neib_ed[j]=v_neib->neib_ed[j+1];
                    v_neib->neib_vert[j+1]=nullptr;
                }
            }
        }
        v_neib->v_neigh_count=v_neib->v_neigh_count-1;
    }
    delete v->p;
    delete v;
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
            Sx=(Sx-(tilt*PBC*lroundl(Sy/twob)));
            Sx=(Sx-(twob*PBC*lroundl(Sx/twob)));
            Sy=(Sy-(twob*PBC*lroundl(Sy/twob)));
            cout<<ve->p->x<<"\t"<<ve->p->y<<"\t"<<Sx+ve->p->x<<"\t"<<Sy+ve->p->y<<"\n";
        }
    }
    if(ve->next)
        display_conn(ve->next);
    else
        return;

}
void update_neighbours(atom Atoms[],int nAtoms)
{
    long double R_CUT;
    //R_CUT=sqrtl(200./(4*3.14*density));
    //R_CUT=powl(20./((4./3.)*3.14*density),1./3.);
    R_CUT=100;
    for(int i=0; i<nAtoms-1; i++)
    {
        for(int j=i+1; j<nAtoms; j++)
        {
            long double drx,dry,drz,dr;
            drx=Atoms[i].p.x-Atoms[j].p.x;
            dry=Atoms[i].p.y-Atoms[j].p.y;
            drz=Atoms[i].p.z-Atoms[j].p.z;
            drx=(drx-(tilt*PBC*lroundl(dry/twob)));
            drx=(drx-(twob*PBC*lroundl(drx/twob)));
            dry=(dry-(twob*PBC*lroundl(dry/twob)));
            drz=(drz-(twob*PBC*lroundl(drz/twob)));
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
void print_face(face *f,int trans=0)
{
////if(trans)
////{
////	cout<<"mol new\n";
////	cout<<"draw material Transparent\n";
////}
////else
////{
////	cout<<"mol new\n";
////	cout<<"draw material Opaque\n";
////}
    cout<<"draw color blue\n";
    cout<<"draw triangle \t{";
    cout<<f->A1.x<<"\t"<<f->A1.y<<"\t"<<f->A1.z<<"}\t{";
    cout<<f->A2.x<<"\t"<<f->A2.y<<"\t"<<f->A2.z<<"}\t{";
    cout<<f->A3.x<<"\t"<<f->A3.y<<"\t"<<f->A3.z<<"}\n";
    cout<<"draw color black\n";
    cout<<"draw line\t{";
    cout<<f->A1.x<<"\t"<<f->A1.y<<"\t"<<f->A1.z<<"}\t{";
    cout<<f->A2.x<<"\t"<<f->A2.y<<"\t"<<f->A2.z<<"}\n";
    cout<<"draw line\t{";
    cout<<f->A1.x<<"\t"<<f->A1.y<<"\t"<<f->A1.z<<"}\t{";
    cout<<f->A3.x<<"\t"<<f->A3.y<<"\t"<<f->A3.z<<"}\n";
    cout<<"draw line\t{";
    cout<<f->A2.x<<"\t"<<f->A2.y<<"\t"<<f->A2.z<<"}\t{";
    cout<<f->A3.x<<"\t"<<f->A3.y<<"\t"<<f->A3.z<<"}\n";
}
void print_delunay_solid(delunay *D,atom Atoms[],int TYPE)
{
    cout<<"##\t"<<D->AT[0]<<"\t"<<D->AT[1]<<"\t"<<D->AT[2]<<"\t"<<D->AT[3]<<"\n";
    cout<<"mol new\n";
    cout<<"draw material Opaque\n";
    //if(D->FACE[0][1][2])
    {
        cout<<"draw color blue\n";
// }
// else
// 	cout<<"draw color green\n";
// {
        cout<<"draw triangle\t{";
        cout<<Atoms[D->AT[0]].p.x<<"\t"<<Atoms[D->AT[0]].p.y<<"\t"<<Atoms[D->AT[0]].p.z<<"}\t{";
        cout<<Atoms[D->AT[1]].p.x<<"\t"<<Atoms[D->AT[1]].p.y<<"\t"<<Atoms[D->AT[1]].p.z<<"}\t{";
        cout<<Atoms[D->AT[2]].p.x<<"\t"<<Atoms[D->AT[2]].p.y<<"\t"<<Atoms[D->AT[2]].p.z<<"}\n";
    }
    // if(D->FACE[0][1][3])
    {
        cout<<"draw color blue\n";
// }
// else
// 	cout<<"draw color green\n";
// {
        cout<<"draw triangle\t{";
        cout<<Atoms[D->AT[0]].p.x<<"\t"<<Atoms[D->AT[0]].p.y<<"\t"<<Atoms[D->AT[0]].p.z<<"}\t{";
        cout<<Atoms[D->AT[1]].p.x<<"\t"<<Atoms[D->AT[1]].p.y<<"\t"<<Atoms[D->AT[1]].p.z<<"}\t{";
        cout<<Atoms[D->AT[3]].p.x<<"\t"<<Atoms[D->AT[3]].p.y<<"\t"<<Atoms[D->AT[3]].p.z<<"}\n";
    }
    // if(D->FACE[1][2][3])
    {
        cout<<"draw color blue\n";
// }
// else
// 	cout<<"draw color green\n";
// {
        cout<<"draw triangle\t{";
        cout<<Atoms[D->AT[1]].p.x<<"\t"<<Atoms[D->AT[1]].p.y<<"\t"<<Atoms[D->AT[1]].p.z<<"}\t{";
        cout<<Atoms[D->AT[2]].p.x<<"\t"<<Atoms[D->AT[2]].p.y<<"\t"<<Atoms[D->AT[2]].p.z<<"}\t{";
        cout<<Atoms[D->AT[3]].p.x<<"\t"<<Atoms[D->AT[3]].p.y<<"\t"<<Atoms[D->AT[3]].p.z<<"}\n";
    }
    //  if(D->FACE[0][2][3])
    {
        cout<<"draw color blue\n";
// }
// else
// 	cout<<"draw color green\n";
// {
        cout<<"draw triangle\t{";
        cout<<Atoms[D->AT[0]].p.x<<"\t"<<Atoms[D->AT[0]].p.y<<"\t"<<Atoms[D->AT[0]].p.z<<"}\t{";
        cout<<Atoms[D->AT[2]].p.x<<"\t"<<Atoms[D->AT[2]].p.y<<"\t"<<Atoms[D->AT[2]].p.z<<"}\t{";
        cout<<Atoms[D->AT[3]].p.x<<"\t"<<Atoms[D->AT[3]].p.y<<"\t"<<Atoms[D->AT[3]].p.z<<"}\n";
    }
////if(D->EDGE[0][1])
////{
////	cout<<"draw line\t{";
////	cout<<Atoms[D->AT[0]].p.x<<"\t"<<Atoms[D->AT[0]].p.y<<"\t"<<Atoms[D->AT[0]].p.z<<"}\t{";
////	cout<<Atoms[D->AT[1]].p.x<<"\t"<<Atoms[D->AT[1]].p.y<<"\t"<<Atoms[D->AT[1]].p.z<<"}\n";
////}
////if(D->EDGE[0][2])
////{
////	cout<<"draw line\t{";
////	cout<<Atoms[D->AT[0]].p.x<<"\t"<<Atoms[D->AT[0]].p.y<<"\t"<<Atoms[D->AT[0]].p.z<<"}\t{";
////	cout<<Atoms[D->AT[2]].p.x<<"\t"<<Atoms[D->AT[2]].p.y<<"\t"<<Atoms[D->AT[2]].p.z<<"}\n";
////}
////if(D->EDGE[0][3])
////{
////	cout<<"draw line\t{";
////	cout<<Atoms[D->AT[0]].p.x<<"\t"<<Atoms[D->AT[0]].p.y<<"\t"<<Atoms[D->AT[0]].p.z<<"}\t{";
////	cout<<Atoms[D->AT[3]].p.x<<"\t"<<Atoms[D->AT[3]].p.y<<"\t"<<Atoms[D->AT[3]].p.z<<"}\n";
////}
////if(D->EDGE[1][2])
////{
////	cout<<"draw line\t{";
////	cout<<Atoms[D->AT[1]].p.x<<"\t"<<Atoms[D->AT[1]].p.y<<"\t"<<Atoms[D->AT[1]].p.z<<"}\t{";
////	cout<<Atoms[D->AT[2]].p.x<<"\t"<<Atoms[D->AT[2]].p.y<<"\t"<<Atoms[D->AT[2]].p.z<<"}\n";
////}
////if(D->EDGE[1][3])
////{
////	cout<<"draw line\t{";
////	cout<<Atoms[D->AT[1]].p.x<<"\t"<<Atoms[D->AT[1]].p.y<<"\t"<<Atoms[D->AT[1]].p.z<<"}\t{";
////	cout<<Atoms[D->AT[3]].p.x<<"\t"<<Atoms[D->AT[3]].p.y<<"\t"<<Atoms[D->AT[3]].p.z<<"}\n";
////}
////if(D->EDGE[2][3])
////{
////	cout<<"draw line\t{";
////	cout<<Atoms[D->AT[2]].p.x<<"\t"<<Atoms[D->AT[2]].p.y<<"\t"<<Atoms[D->AT[2]].p.z<<"}\t{";
////	cout<<Atoms[D->AT[3]].p.x<<"\t"<<Atoms[D->AT[3]].p.y<<"\t"<<Atoms[D->AT[3]].p.z<<"}\n";
////}
    cout<<"mol new\n";
    cout<<"draw material Transparent\n";
    cout<<"draw color red\n";
    cout<<"draw sphere\t{";
    cout<<Atoms[D->AT[0]].p.x<<"\t"<<Atoms[D->AT[0]].p.y<<"\t"<<Atoms[D->AT[0]].p.z<<"}\tradius\t"<<Atoms[D->AT[0]].radius+r_cut<<"\tresolution 10\n";
    cout<<"draw sphere\t{";
    cout<<Atoms[D->AT[1]].p.x<<"\t"<<Atoms[D->AT[1]].p.y<<"\t"<<Atoms[D->AT[1]].p.z<<"}\tradius\t"<<Atoms[D->AT[1]].radius+r_cut<<"\tresolution 10\n";
    cout<<"draw sphere\t{";
    cout<<Atoms[D->AT[2]].p.x<<"\t"<<Atoms[D->AT[2]].p.y<<"\t"<<Atoms[D->AT[2]].p.z<<"}\tradius\t"<<Atoms[D->AT[2]].radius+r_cut<<"\tresolution 10\n";
    cout<<"draw sphere\t{";
    cout<<Atoms[D->AT[3]].p.x<<"\t"<<Atoms[D->AT[3]].p.y<<"\t"<<Atoms[D->AT[3]].p.z<<"}\tradius\t"<<Atoms[D->AT[3]].radius+r_cut<<"\tresolution 10\n";
}
void display_atom(atom *ATOM)
{
    cout<<"draw sphere\t{";
    cout<<ATOM->p.x<<"\t"<<ATOM->p.y<<"\t"<<ATOM->p.z<<"}\tradius\t"<<ATOM->radius+r_cut<<"\tresolution 10\n";
}
void print_delunay(delunay *D,atom Atoms[],int TYPE)
{
    cout<<"##\t";
    for(int i=0; i<4; i++)
    {
        cout<<D->AT[i]<<"\t";
    }
    //for(int i=0;i<4;i++)
    //	for(int j=0;j<4;j++)
    //		for(int k=0;k<4;k++)
    //			for(uu
    cout<<"\n";
    //return ;
    atom *ATOM;
    ATOM=&(Atoms[D->AT[0]]);
    long double Sx,Sy,Sz;
    long double Px,Py,Pz;
    Sx=ATOM->p.x-Atoms[D->AT[1]].p.x;
    Sy=ATOM->p.y-Atoms[D->AT[1]].p.y;
    Sz=ATOM->p.z-Atoms[D->AT[1]].p.z;
    Sx=(Sx-(tilt*PBC*lroundl(Sy/twob)));
    Sx=(Sx-(twob*PBC*lroundl(Sx/twob)));
    Sy=(Sy-(twob*PBC*lroundl(Sy/twob)));
    Sz=(Sz-(twob*PBC*lroundl(Sz/twob)));
    cout<<"draw line\t{";
    cout<<ATOM->p.x<<"\t"<<ATOM->p.y<<"\t"<<ATOM->p.z<<"}\t{"<<ATOM->p.x-Sx<<"\t"<<ATOM->p.y-Sy<<"\t"<<ATOM->p.z-Sz<<"}\twidth	1\n";
    Sx=ATOM->p.x-Atoms[D->AT[2]].p.x;
    Sy=ATOM->p.y-Atoms[D->AT[2]].p.y;
    Sz=ATOM->p.z-Atoms[D->AT[2]].p.z;
    Sx=(Sx-(tilt*PBC*lroundl(Sy/twob)));
    Sx=(Sx-(twob*PBC*lroundl(Sx/twob)));
    Sy=(Sy-(twob*PBC*lroundl(Sy/twob)));
    Sz=(Sz-(twob*PBC*lroundl(Sz/twob)));
    cout<<"draw line\t{";
    cout<<ATOM->p.x<<"\t"<<ATOM->p.y<<"\t"<<ATOM->p.z<<"}\t{"<<ATOM->p.x-Sx<<"\t"<<ATOM->p.y-Sy<<"\t"<<ATOM->p.z-Sz<<"}\twidth	1\n";
    Sx=ATOM->p.x-Atoms[D->AT[3]].p.x;
    Sy=ATOM->p.y-Atoms[D->AT[3]].p.y;
    Sz=ATOM->p.z-Atoms[D->AT[3]].p.z;
    Sx=(Sx-(tilt*PBC*lroundl(Sy/twob)));
    Sx=(Sx-(twob*PBC*lroundl(Sx/twob)));
    Sy=(Sy-(twob*PBC*lroundl(Sy/twob)));
    Sz=(Sz-(twob*PBC*lroundl(Sz/twob)));
    cout<<"draw line\t{";
    cout<<ATOM->p.x<<"\t"<<ATOM->p.y<<"\t"<<ATOM->p.z<<"}\t{"<<ATOM->p.x-Sx<<"\t"<<ATOM->p.y-Sy<<"\t"<<ATOM->p.z-Sz<<"}\twidth	1\n";
    Sx=ATOM->p.x-Atoms[D->AT[2]].p.x;
    Sy=ATOM->p.y-Atoms[D->AT[2]].p.y;
    Sz=ATOM->p.z-Atoms[D->AT[2]].p.z;
    Sx=(Sx-(tilt*PBC*lroundl(Sy/twob)));
    Sx=(Sx-(twob*PBC*lroundl(Sx/twob)));
    Sy=(Sy-(twob*PBC*lroundl(Sy/twob)));
    Sz=(Sz-(twob*PBC*lroundl(Sz/twob)));
    Px=ATOM->p.x-Atoms[D->AT[1]].p.x;
    Py=ATOM->p.y-Atoms[D->AT[1]].p.y;
    Pz=ATOM->p.z-Atoms[D->AT[1]].p.z;
    Px=(Px-(tilt*PBC*lroundl(Py/twob)));
    Px=(Px-(twob*PBC*lroundl(Px/twob)));
    Py=(Py-(twob*PBC*lroundl(Py/twob)));
    Pz=(Pz-(twob*PBC*lroundl(Pz/twob)));
    cout<<"draw line\t{";
    cout<<ATOM->p.x-Px<<"\t"<<ATOM->p.y-Py<<"\t"<<ATOM->p.z-Pz<<"}\t{"<<ATOM->p.x-Sx<<"\t"<<ATOM->p.y-Sy<<"\t"<<ATOM->p.z-Sz<<"}\twidth	1\n";
    Sx=ATOM->p.x-Atoms[D->AT[3]].p.x;
    Sy=ATOM->p.y-Atoms[D->AT[3]].p.y;
    Sz=ATOM->p.z-Atoms[D->AT[3]].p.z;
    Sx=(Sx-(tilt*PBC*lroundl(Sy/twob)));
    Sx=(Sx-(twob*PBC*lroundl(Sx/twob)));
    Sy=(Sy-(twob*PBC*lroundl(Sy/twob)));
    Sz=(Sz-(twob*PBC*lroundl(Sz/twob)));
    Px=ATOM->p.x-Atoms[D->AT[1]].p.x;
    Py=ATOM->p.y-Atoms[D->AT[1]].p.y;
    Pz=ATOM->p.z-Atoms[D->AT[1]].p.z;
    Px=(Px-(tilt*PBC*lroundl(Py/twob)));
    Px=(Px-(twob*PBC*lroundl(Px/twob)));
    Py=(Py-(twob*PBC*lroundl(Py/twob)));
    Pz=(Pz-(twob*PBC*lroundl(Pz/twob)));
    cout<<"draw line\t{";
    cout<<ATOM->p.x-Px<<"\t"<<ATOM->p.y-Py<<"\t"<<ATOM->p.z-Pz<<"}\t{"<<ATOM->p.x-Sx<<"\t"<<ATOM->p.y-Sy<<"\t"<<ATOM->p.z-Sz<<"}\twidth	1\n";
    Sx=ATOM->p.x-Atoms[D->AT[2]].p.x;
    Sy=ATOM->p.y-Atoms[D->AT[2]].p.y;
    Sz=ATOM->p.z-Atoms[D->AT[2]].p.z;
    Sx=(Sx-(tilt*PBC*lroundl(Sy/twob)));
    Sx=(Sx-(twob*PBC*lroundl(Sx/twob)));
    Sy=(Sy-(twob*PBC*lroundl(Sy/twob)));
    Sz=(Sz-(twob*PBC*lroundl(Sz/twob)));
    Px=ATOM->p.x-Atoms[D->AT[3]].p.x;
    Py=ATOM->p.y-Atoms[D->AT[3]].p.y;
    Pz=ATOM->p.z-Atoms[D->AT[3]].p.z;
    Px=(Px-(tilt*PBC*lroundl(Py/twob)));
    Px=(Px-(twob*PBC*lroundl(Px/twob)));
    Py=(Py-(twob*PBC*lroundl(Py/twob)));
    Pz=(Pz-(twob*PBC*lroundl(Pz/twob)));
    cout<<"draw line\t{";
    cout<<ATOM->p.x-Px<<"\t"<<ATOM->p.y-Py<<"\t"<<ATOM->p.z-Pz<<"}\t{"<<ATOM->p.x-Sx<<"\t"<<ATOM->p.y-Sy<<"\t"<<ATOM->p.z-Sz<<"}\twidth	1\n";
////cout<<"draw sphere\t{";
////cout<<Atoms[D->AT[0]].p.x<<"\t"<<Atoms[D->AT[0]].p.y<<"\t"<<Atoms[D->AT[0]].p.z<<"}\tradius\t"<<Atoms[D->AT[0]].radius+r_cut<<"\n";
////cout<<"draw sphere\t{";
////cout<<Atoms[D->AT[1]].p.x<<"\t"<<Atoms[D->AT[1]].p.y<<"\t"<<Atoms[D->AT[1]].p.z<<"}\tradius\t"<<Atoms[D->AT[1]].radius+r_cut<<"\n";
////cout<<"draw sphere\t{";
////cout<<Atoms[D->AT[2]].p.x<<"\t"<<Atoms[D->AT[2]].p.y<<"\t"<<Atoms[D->AT[2]].p.z<<"}\tradius\t"<<Atoms[D->AT[2]].radius+r_cut<<"\n";
////cout<<"draw sphere\t{";
////cout<<Atoms[D->AT[3]].p.x<<"\t"<<Atoms[D->AT[3]].p.y<<"\t"<<Atoms[D->AT[3]].p.z<<"}\tradius\t"<<Atoms[D->AT[3]].radius+r_cut<<"\n";
}
site cross_product(site a1,site a2)
{
    site cp;
    cp.x=a2.z*a1.y-a1.z*a2.y;
    cp.y=a1.z*a2.x-a1.x*a2.z;
    cp.z=a2.y*a1.x-a1.y*a2.x;
    return cp;
}
long double determinant(long double a[3][3])
{
    return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
}
site cramer(long double a[3][3],long double b[3],int debug=0)
{
    long double x[3];
    long double det_a=determinant(a);
    long double det_a1=0.;
    long double a_1[3][3];
    if(debug)
    {
        cout<<"[";
        for(int i=0; i<3; i++)
        {
            cout<<"[";
            for(int j=0; j<3; j++)
            {
                cout<<a[i][j]<<",";
                //a_1[i][j]=a[i][j];
            }
            cout<<"]";
            //cout<<"\n";
        }
        cout<<"]\n";
        cout<<"[";
        for(int i=0; i<3; i++)
        {
            cout<<b[i]<<",";
        }
        cout<<"]\n";
    }
////cout<<det_a<<"\n";

    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++)
            {
                a_1[j][k]=a[j][k];
            }
        }
        for(int j=0; j<3; j++)
        {
            a_1[j][i]=b[j]/det_a;
        }
        //for(int i=0;i<3;i++)
        //{
        //    for(int j=0;j<3;j++)
        //    {
        //        cout<<a_1[i][j]<<"\t";
        //    }
        //    cout<<"\n";
        //}
        //cout<<"\n\n";
        det_a1=determinant(a_1);
        //cout<<det_a1<<"\n";
        x[i]=det_a1;
        ////for(int j=0;j<3;j++)
        ////{
        ////    a_1[j][i]=a[j][i];
        ////}
    }
    site p;
    p.x=x[0];
    p.y=x[1];
    p.z=x[2];
    return p;
}
site center_of_triangle(site A1,site A2,site A3,long double rS,long double rA,long double rB,int debug=0)
{
    long double ax,ay,az;
    long double bx,by,bz;
    long double X,Y,Z;
    site center;
    ax=A2.x-A1.x;
    ay=A2.y-A1.y;
    az=A2.z-A1.z;
    bx=A3.x-A1.x;
    by=A3.y-A1.y;
    bz=A3.z-A1.z;
    ax=(ax-(tilt*PBC*lroundl(ay/twob)));
    ax=(ax-(twob*PBC*lroundl(ax/twob)));
    ay=(ay-(twob*PBC*lroundl(ay/twob)));
    az=(az-(twob*PBC*lroundl(az/twob)));
    bx=(bx-(tilt*PBC*lroundl(by/twob)));
    bx=(bx-(twob*PBC*lroundl(bx/twob)));
    by=(by-(twob*PBC*lroundl(by/twob)));
    bz=(bz-(twob*PBC*lroundl(bz/twob)));
    //long double rA,rS,rB;
    long double XA,YA,ZA,XB,YB,ZB;
    long double DISA;
    long double DISB;
    XA=ax;
    YA=ay;
    ZA=az;
    XB=bx;
    YB=by;
    ZB=bz;
    DISA=sqrtl(XA*XA+YA*YA+ZA*ZA);
    DISB=sqrtl(XB*XB+YB*YB+ZB*ZB);
    long double a1,b1,c1,a2,b2,c2,a,b,c;
    a=ZB*YA-ZA*YB;
    b=ZA*XB-XA*ZB;
    c=YB*XA-YA*XB;
    long double B[3],A[3][3];
    B[0]=(DISA*DISA+rS*rS-rA*rA)/2.;
    B[1]=(DISB*DISB+rS*rS-rB*rB)/2.;
    B[2]=0.;
    A[0][0]=XA;
    A[0][1]=YA;
    A[0][2]=ZA;
    A[1][0]=XB;
    A[1][1]=YB;
    A[1][2]=ZB;
    A[2][0]=a;
    A[2][1]=b;
    A[2][2]=c;
    center=cramer(A,B,debug);
    center.x=center.x+A1.x;
    center.y=center.y+A1.y;
    center.z=center.z+A1.z;
    return center;
}
long double volume_delunay(delunay *D,atom Atoms[])
{
    long double a[3][3];
    a[0][0]=Atoms[D->AT[1]].p.x-Atoms[D->AT[0]].p.x;
    a[0][1]=Atoms[D->AT[1]].p.y-Atoms[D->AT[0]].p.y;
    a[0][2]=Atoms[D->AT[1]].p.z-Atoms[D->AT[0]].p.z;

    a[0][0]=(a[0][0]-(tilt*PBC*lroundl(a[0][0]/twob)));
    a[0][0]=(a[0][0]-(twob*PBC*lroundl(a[0][0]/twob)));
    a[0][1]=(a[0][1]-(twob*PBC*lroundl(a[0][1]/twob)));
    a[0][2]=(a[0][2]-(twob*PBC*lroundl(a[0][2]/twob)));

    a[1][0]=Atoms[D->AT[2]].p.x-Atoms[D->AT[0]].p.x;
    a[1][1]=Atoms[D->AT[2]].p.y-Atoms[D->AT[0]].p.y;
    a[1][2]=Atoms[D->AT[2]].p.z-Atoms[D->AT[0]].p.z;

    a[1][0]=(a[1][0]-(tilt*PBC*lroundl(a[1][0]/twob)));
    a[1][0]=(a[1][0]-(twob*PBC*lroundl(a[1][0]/twob)));
    a[1][1]=(a[1][1]-(twob*PBC*lroundl(a[1][1]/twob)));
    a[1][2]=(a[1][2]-(twob*PBC*lroundl(a[1][2]/twob)));

    a[2][0]=Atoms[D->AT[3]].p.x-Atoms[D->AT[0]].p.x;
    a[2][1]=Atoms[D->AT[3]].p.y-Atoms[D->AT[0]].p.y;
    a[2][2]=Atoms[D->AT[3]].p.z-Atoms[D->AT[0]].p.z;

    a[2][0]=(a[2][0]-(tilt*PBC*lroundl(a[2][0]/twob)));
    a[2][0]=(a[2][0]-(twob*PBC*lroundl(a[2][0]/twob)));
    a[2][1]=(a[2][1]-(twob*PBC*lroundl(a[2][1]/twob)));
    a[2][2]=(a[2][2]-(twob*PBC*lroundl(a[2][2]/twob)));

    long double det=determinant(a);
    return abs(det/6.);
}
void delete_delunay(atom Atoms[],delunay *D,int TYPE)
{
    FULLSETD.remove_delunay(D);
}
void create_delunay(atom Atoms[],int A1,int A2,int A3,int A4,delunay *D,int TYPE)
{
    //cout<<"ghere\n";
    convex_vol=convex_vol+volume_delunay(D,Atoms);
    vertice *temp_v=nullptr;
    temp_v=new vertice;
    temp_v->p=new site;
    temp_v->p->x=D->circum_x;
    temp_v->p->y=D->circum_y;
    temp_v->p->z=D->circum_z;
    temp_v->D=D;
    D->v=temp_v;
    V->insert_vertice(temp_v,TYPE);
    for(int a=0; a<4; a++)
    {
        for(int b=a+1; b<4; b++)
        {
            site a1,a2,a3,center;
            site midp;
            long double r1,r2,r3;
            a1=Atoms[D->AT[a]].p;
            r1=Atoms[D->AT[a]].radius+r_cut;
            a2=Atoms[D->AT[b]].p;
            r2=Atoms[D->AT[b]].radius+r_cut;
            //D->MIDP[a][b]=
            long double ax,ay,az;
            long double X,Y,Z,DISA,l;
            ax=a2.x-a1.x;
            ay=a2.y-a1.y;
            az=a2.z-a1.z;
            ax=(ax-(tilt*PBC*lroundl(ay/twob)));
            ax=(ax-(twob*PBC*lroundl(ax/twob)));
            ay=(ay-(twob*PBC*lroundl(ay/twob)));
            az=(az-(twob*PBC*lroundl(az/twob)));
            long double XA,YA,ZA;
            long double xA,yA,zA;
            XA=ax;
            YA=ay;
            ZA=az;
            DISA=sqrtl(XA*XA+YA*YA+ZA*ZA);
            l=0.5*(DISA+(r1*r1-r2*r2)/DISA);
            xA=l/DISA*XA;
            yA=l/DISA*YA;
            zA=l/DISA*ZA;
            xA=xA+a1.x;
            yA=yA+a1.y;
            zA=zA+a1.z;
            D->MIDP[a][b].x=xA;
            D->MIDP[a][b].y=yA;
            D->MIDP[a][b].z=zA;

            D->MIDP[b][a].x=xA;
            D->MIDP[b][a].y=yA;
            D->MIDP[b][a].z=zA;
            ////{
            ////	cout<<"B\t"<<a<<"\t"<<b<<"\n";;
            ////	display_SITE(&D->MIDP[a][b]);
            ////}
            //cout<<a<<"\t"<<b<<"\n";
            for(int c=b+1; c<4; c++)
            {
                a3=Atoms[D->AT[c]].p;

                r3=Atoms[D->AT[c]].radius+r_cut;

                //cout<<r1<<"\t"<<r2<<"\t"<<r3<<"\n";
                ////display_SITE(&a1);
                ////display_SITE(&a2);
                ////display_SITE(&a3);

                center=center_of_triangle(a1,a2,a3,r1,r2,r3);

                D->MID[a][b][c]=center;
                D->MID[a][c][b]=center;

                D->MID[b][a][c]=center;
                D->MID[b][c][a]=center;

                D->MID[c][a][b]=center;
                D->MID[c][b][a]=center;

            }
        }
    }
    for(int k=0; k<4; k++)
    {
        atom *ATOM;
        atom *Ai;
        atom *Aj;
        atom *Ak;
        ATOM=&(Atoms[D->AT[k]]);
        Ai=&(Atoms[D->AT[(k+1)%4]]);
        Aj=&(Atoms[D->AT[(k+2)%4]]);
        Ak=&(Atoms[D->AT[(k+3)%4]]);
        int a1=-1,a2=-1,a3=-1;
        {
            if(!contigous[ATOM->index][D->AT[(k+1)%4]])
            {
                ATOM->contigous[ATOM->conti]=D->AT[(k+1)%4];
                Ai->contigous[Ai->conti]=ATOM->index;
                ATOM->conti_index[D->AT[(k+1)%4]]=ATOM->conti;
                Ai->conti_index[ATOM->index]=Ai->conti;

                contigous[ATOM->index][D->AT[(k+1)%4]]=1;
                contigous[D->AT[(k+1)%4]][ATOM->index]=1;

                a1=ATOM->conti;
                ATOM->conti++;
                Ai->conti++;
            }
            else
            {
                a1=ATOM->conti_index[D->AT[(k+1)%4]];
            }

            if(!contigous[ATOM->index][D->AT[(k+2)%4]])
            {
                ATOM->contigous[ATOM->conti]=D->AT[(k+2)%4];
                Aj->contigous[Aj->conti]=ATOM->index;
                ATOM->conti_index[D->AT[(k+2)%4]]=ATOM->conti;
                Aj->conti_index[ATOM->index]=Aj->conti;

                contigous[ATOM->index][D->AT[(k+2)%4]]=1;
                contigous[D->AT[(k+2)%4]][ATOM->index]=1;

                a2=ATOM->conti;
                ATOM->conti++;
                Aj->conti++;
            }
            else
            {
                a2=ATOM->conti_index[D->AT[(k+2)%4]];
            }

            if(!contigous[ATOM->index][D->AT[(k+3)%4]])
            {
                ATOM->contigous[ATOM->conti]=D->AT[(k+3)%4];
                Ak->contigous[Ak->conti]=ATOM->index;
                ATOM->conti_index[D->AT[(k+3)%4]]=ATOM->conti;
                Ak->conti_index[ATOM->index]=Ak->conti;

                contigous[ATOM->index][D->AT[(k+3)%4]]=1;
                contigous[D->AT[(k+3)%4]][ATOM->index]=1;

                a3=ATOM->conti;
                ATOM->conti++;
                Ak->conti++;
            }
            else
            {
                a3=ATOM->conti_index[D->AT[(k+3)%4]];
            }

            if(ATOM->part_c[a1][a2]==0)
            {
                ATOM->edge_index[a1]++;
                ATOM->edge_index[a2]++;
            }
            if(ATOM->part_c[a1][a3]==0)
            {
                ATOM->edge_index[a1]++;
                ATOM->edge_index[a3]++;
            }
            if(ATOM->part_c[a2][a3]==0)
            {
                ATOM->edge_index[a2]++;
                ATOM->edge_index[a3]++;
            }

            ATOM->part_c[a1][a2]++;
            ATOM->part_c[a2][a1]++;

            ATOM->part_c[a1][a3]++;
            ATOM->part_c[a3][a1]++;

            ATOM->part_c[a2][a3]++;
            ATOM->part_c[a3][a2]++;

            ATOM->MIDP[a1][a2]=D->MID[k][(k+1)%4][(k+2)%4];
            ATOM->MIDP[a2][a1]=D->MID[k][(k+1)%4][(k+2)%4];

            ATOM->MIDP[a1][a3]=D->MID[k][(k+1)%4][(k+3)%4];
            ATOM->MIDP[a3][a1]=D->MID[k][(k+1)%4][(k+3)%4];

            ATOM->MIDP[a2][a3]=D->MID[k][(k+2)%4][(k+3)%4];
            ATOM->MIDP[a3][a2]=D->MID[k][(k+2)%4][(k+3)%4];

            ATOM->RMID[a1]=D->MIDP[k][(k+1)%4];
            ATOM->RMID[a2]=D->MIDP[k][(k+2)%4];
            ATOM->RMID[a3]=D->MIDP[k][(k+3)%4];

            if(!ATOM->D_FIRST)
            {
                ATOM->D_FIRST=new container_delunay;
                ATOM->D_FIRST->D=D;
            }
            else
            {
                container_delunay *temp;
                temp=ATOM->D_FIRST;
                while(temp->next)
                {
                    temp=temp->next;
                }
                temp->next=new container_delunay;
                temp->next->D=D;
            }
        }
    }
}
void print_edge(edge *e)
{
    cout<<"draw line\t{";
    cout<<e->origin->p->x<<"\t"<<e->origin->p->y<<"\t"<<e->origin->p->z<<"}\t{";
    cout<<e->destination->p->x<<"\t"<<e->destination->p->y<<"\t"<<e->destination->p->z<<"}\n";
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
    int A1,A2,A3,A4;
    for(int i=0; i<ATOM->neighbours; i++)
    {
        long double X,Y,Z,x,y,z;
        long double rA,rS,DIS,l,dis_i,tan_sq;
        X=Atoms[ATOM->neighlist[i]].p.x-ATOM->p.x;
        Y=Atoms[ATOM->neighlist[i]].p.y-ATOM->p.y;
        Z=Atoms[ATOM->neighlist[i]].p.z-ATOM->p.z;
        X=(X-(tilt*PBC*lroundl(Y/twob)));
        X=(X-(twob*PBC*lroundl(X/twob)));
        Y=(Y-(twob*PBC*lroundl(Y/twob)));
        Z=(Z-(twob*PBC*lroundl(Z/twob)));
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
            midax=x+ATOM->p.x;
            miday=y+ATOM->p.y;
            midaz=z+ATOM->p.z;
            lmin=l;
            if(DIS>rS+rA)
            {
                binv1=1;
            }
        }
    }
    A1=ATOM->index;
    A2=nearest;
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
            atom M=Atoms[nearest];
            atom R=Atoms[ATOM->neighlist[i]];
            site center;
            center=center_of_triangle(L.p,M.p,R.p,L.radius,M.radius,R.radius);
            long double tan_sq=distancesq(center,L.p)-L.radius*L.radius;
            if(DIS_MIN > tan_sq)
            {
                DIS_MIN=tan_sq;
                DIS_atom=ATOM->neighlist[i];
            }
        }
    }
//	cout<<"dis\t"<<DIS_atom<<"\n";
    A3=DIS_atom;
    DIS_MIN=box*box*box;
    long double circx_s=circx;
    long double circy_s=circy;
    long double circz_s=circz;
//	cout<<DIS_atom<<"\n";
    atom L=*ATOM;
    atom M=Atoms[nearest];
    atom N=Atoms[DIS_atom];
    site center;
    for(int i=0; i<ATOM->neighbours; i++)
    {
        if(ATOM->neighlist[i]!=nearest && ATOM->neighlist[i]!=A3)
        {
            atom R=Atoms[ATOM->neighlist[i]];
            long double ax=M.p.x-L.p.x;
            long double ay=M.p.y-L.p.y;
            long double az=M.p.z-L.p.z;

            long double bx=R.p.x-L.p.x;
            long double by=R.p.y-L.p.y;
            long double bz=R.p.z-L.p.z;

            long double cx=N.p.x-L.p.x;
            long double cy=N.p.y-L.p.y;
            long double cz=N.p.z-L.p.z;

            ax=(ax-(tilt*PBC*lroundl(ay/twob)));
            ax=(ax-(twob*PBC*lroundl(ax/twob)));
            ay=(ay-(twob*PBC*lroundl(ay/twob)));
            az=(az-(twob*PBC*lroundl(az/twob)));
            bx=(bx-(tilt*PBC*lroundl(by/twob)));
            bx=(bx-(twob*PBC*lroundl(bx/twob)));
            by=(by-(twob*PBC*lroundl(by/twob)));
            bz=(bz-(twob*PBC*lroundl(bz/twob)));
            cx=(cx-(tilt*PBC*lroundl(cy/twob)));
            cx=(cx-(twob*PBC*lroundl(cx/twob)));
            cy=(cy-(twob*PBC*lroundl(cy/twob)));
            cz=(cz-(twob*PBC*lroundl(cz/twob)));
            long double rA,rS,rB,rN;
            long double XA,YA,ZA,XB,YB,ZB,XC,YC,ZC;
            long double xA,yA,zA,xB,yB,zB;
            long double l;
            long double DISA;
            long double DISB;
            long double DISC;
            long double MA,MB,INMA,INMB;
            long double CA,CB;
            long double tan_sq;
            long double B[3],A[3][3];
            XA=ax;
            YA=ay;
            ZA=az;
            XB=bx;
            YB=by;
            ZB=bz;
            XC=cx;
            YC=cy;
            ZC=cz;
            rA=M.radius+r_cut;
            rB=R.radius+r_cut;
            rS=L.radius+r_cut;
            rN=N.radius+r_cut;
            DISA=sqrtl(XA*XA+YA*YA+ZA*ZA);
            DISB=sqrtl(XB*XB+YB*YB+ZB*ZB);
            DISC=sqrtl(XC*XC+YC*YC+ZC*ZC);
            B[0]=(DISA*DISA+rS*rS-rA*rA)/2.;
            B[1]=(DISB*DISB+rS*rS-rB*rB)/2.;
            B[2]=(DISC*DISC+rS*rS-rN*rN)/2.;
            A[0][0]=XA;
            A[0][1]=YA;
            A[0][2]=ZA;
            A[1][0]=XB;
            A[1][1]=YB;
            A[1][2]=ZB;
            A[2][0]=XC;
            A[2][1]=YC;
            A[2][2]=ZC;
            center=cramer(A,B);
            tan_sq=center.x*center.x+center.y*center.y+center.z*center.z-rS*rS;
            if(DIS_MIN > tan_sq)
            {
                DIS_MIN=tan_sq;
                DIS_atom=ATOM->neighlist[i];
                circx_s=center.x+L.p.x;
                circy_s=center.y+L.p.y;
                circz_s=center.z+L.p.z;

            }
        }
    }

    A4=DIS_atom;

    std::vector<int> EV1 {A1,A2,A3,A4};
    std::sort(EV1.begin(),EV1.end());
    delunay *D;

    D=new delunay;
    FULLSETD.insert_delunay(D);

    D->AT[0]=EV1[0];
    D->AT[1]=EV1[1];
    D->AT[2]=EV1[2];
    D->AT[3]=EV1[3];
    if(D->AT[0]!=ATOM->index)
    {
        //THis shouldn't happen because if this happens it means in first delunay you made a delunay which involves atoms which were analyzed p
        //previously. If this delunay was acceptable then this atoms delunay should have been found before and we wouldn't be here
        //I have to do this to get this thing working. A good progrmmer would find some other way
        // All that is left in front of me this ad hoc trick.
        //cout<<"error\n";
        return;
    }
    D->circum_x=circx_s;
    D->circum_y=circy_s;
    D->circum_z=circz_s;

    create_delunay(Atoms,A1,A2,A3,A4,D,TYPE);
}
delunay* constr_del(atom *ATOM,atom Atoms[],int TYPE,long double p,long double q,long double r,long double Sx,long double Sy,long double Sz,int sign,int A,int B,int O,long double rA,long double rB,delunay *D,int debug=0)
{
//		cout<<"hereasdasdasn\n";
    long double a=ATOM->p.x;
    long double b=ATOM->p.y;
    long double c=ATOM->p.z;
    long double rS=ATOM->radius+r_cut;
    long double DIS,X,Y,Z;
    long double Y_MIN=box*box*box;
    long double circx,circy,circz;
    long double X1_s,Y1_s,Z1_s;
    long double X2_s,Y2_s,Z2_s;
    long double X3_s,Y3_s,Z3_s;
    int A1,A2,A3,A4;
    int DIS_atom=-1;
    A1=ATOM->index;
    A2=ATOM->contigous[A];
    A3=ATOM->contigous[B];
    int atleastoneatom=0;
    for(int j=0; j<ATOM->neighbours; j++)
    {
        if(debug)
            cout<<j<<"\t"<<ATOM->neighlist[j]<<"\n";
        int sign_N;
        long double x=Atoms[ATOM->neighlist[j]].p.x-a;
        long double y=Atoms[ATOM->neighlist[j]].p.y-b;
        long double z=Atoms[ATOM->neighlist[j]].p.z-c;
        long double rN=Atoms[ATOM->neighlist[j]].radius+r_cut;
        x=(x-(tilt*PBC*lroundl(y/twob)));
        x=(x-(twob*PBC*lroundl(x/twob)));
        y=(y-(twob*PBC*lroundl(y/twob)));
        z=(z-(twob*PBC*lroundl(z/twob)));
        long double ax=(Sy*r-Sz*q);
        long double ay=(Sz*p-Sx*r);
        long double az=(Sx*q-Sy*p);
        long double overlap=(x*ax+y*ay+z*az);
        if(overlap<0.)
            sign_N=1;
        else
            sign_N=-1;
        int flag=1;
        //to avoid atoms with have two delunay
        for(int k=0; k<ATOM->conti; k++)
        {
            if(ATOM->contigous[k]==ATOM->neighlist[j])
            {
                if(k==O ||  k==A || k==B )
                {
                    flag=0;
                    break;
                }

                int s=0;
                for(int p=0; p<ATOM->conti; p++)
                {
                    s=s+ATOM->part_c[k][p];
                }
                if(s==2*ATOM->edge_index[k])
                {
                    flag=0;
                    break;
                }
            }
        }

        if(debug)
            cout<<sign<<"\t"<<sign_N<<"\t"<<flag<<"\n";
        if(sign!=sign_N && flag)
        {
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
            long double B[3],A[3][3];
            B[0]=(DISA*DISA+rS*rS-rA*rA)/2.;
            B[1]=(DISB*DISB+rS*rS-rB*rB)/2.;
            B[2]=(DIS*DIS+rS*rS-rN*rN)/2.;
            A[0][0]=p;
            A[0][1]=q;
            A[0][2]=r;
            A[1][0]=Sx;
            A[1][1]=Sy;
            A[1][2]=Sz;
            A[2][0]=x;
            A[2][1]=y;
            A[2][2]=z;
            site center;
            center=cramer(A,B);
            long double norm=sqrtl(ax*ax+ay*ay+az*az);
            X=center.x;
            Y=center.y;
            Z=center.z;
            long double Y_AXIS;//=sqrtl(powl(X-X2,2)+powl(Y-Y2,2)+powl(Z-Z2,2));
            ////cout<<Y_AXIS<<" 1\n";
            Y_AXIS=1./norm*(ax*X+ay*Y+az*Z);
            Y_AXIS=abs(Y_AXIS);
            //cout<<Y_AXIS<<" 2\n";
            int sign_C;
            long double overlap=(X*ax+Y*ay+Z*az);
            if(overlap<0.)
                sign_C=1;
            else
                sign_C=-1;
            if(sign_C!=sign_N)
            {
                Y_AXIS=-1.*Y_AXIS;
            }
            X=X+a;
            Y=Y+b;
            Z=Z+c;

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
    if(DIS_atom != -1)
    {
        A4=DIS_atom;
        //cout<<A1<<"\t"<<A2<<"\t"<<A3<<"\t"<<A4<<"\n";
        std::vector<int> EV1 {A1,A2,A3,A4};
        std::sort(EV1.begin(),EV1.end());
        D=new delunay;
        FULLSETD.insert_delunay(D);
        D->AT[0]=EV1[0];
        D->AT[1]=EV1[1];
        D->AT[2]=EV1[2];
        D->AT[3]=EV1[3];
        D->circum_x=circx;
        D->circum_y=circy;
        D->circum_z=circz;
        create_delunay(Atoms,A1,A2,A3,A4,D,TYPE);
        return D;
    }
    else
    {
        return D;
    }
}
void complete_del_2(atom *ATOM,atom Atoms[],int nAtoms,int TYPE)
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
    for(int i=0; i<ATOM->conti; i++)
    {
        int s=0;
        //cout<<Atoms[ATOM->contigous[i]].p.x<<"\t"<<Atoms[ATOM->contigous[i]].p.y<<"\t"<<Atoms[ATOM->contigous[i]].p.z<<"\t"<<Atoms[ATOM->contigous[i]].radius<<"\n";
        for(int p=0; p<ATOM->conti; p++)
        {
            s=s+ATOM->part_c[i][p];
        }

        for(int j=0; j<ATOM->conti; j++)
        {

            if(i!=j)
            {
                if(ATOM->part_c[i][j]==1)
                {
                    for(int k=0; k<ATOM->conti; k++)
                    {
                        if(k!=i && k!=j)
                        {
                            container_delunay *temp;
                            temp=ATOM->D_FIRST;
                            while(1)
                            {
                                std::vector<int> EV1 {ATOM->index,ATOM->contigous[i],ATOM->contigous[j],ATOM->contigous[k]};//,A3,A4;
                                std::sort(EV1.begin(),EV1.end());
                                if(temp->D->AT[0]==EV1[0] && temp->D->AT[1]==EV1[1] && temp->D->AT[2]==EV1[2] && temp->D->AT[3]==EV1[3])
                                {
                                    delunay *D_ONE=nullptr,*D_TWO=nullptr;
                                    D_ONE=temp->D;
                                    a=ATOM->p.x;
                                    b=ATOM->p.y;
                                    c=ATOM->p.z;
                                    p=Atoms[ATOM->contigous[i]].p.x-a;
                                    q=Atoms[ATOM->contigous[i]].p.y-b;
                                    r=Atoms[ATOM->contigous[i]].p.z-c;
                                    rA=Atoms[ATOM->contigous[i]].radius+r_cut;
                                    Sx=Atoms[ATOM->contigous[j]].p.x-a;
                                    Sy=Atoms[ATOM->contigous[j]].p.y-b;
                                    Sz=Atoms[ATOM->contigous[j]].p.z-c;
                                    rB=Atoms[ATOM->contigous[j]].radius+r_cut;
                                    Px=Atoms[ATOM->contigous[k]].p.x-a;
                                    Py=Atoms[ATOM->contigous[k]].p.y-b;
                                    Pz=Atoms[ATOM->contigous[k]].p.z-c;
                                    rC=Atoms[ATOM->contigous[k]].radius+r_cut;
                                    p=(p-(tilt*PBC*lroundl(q/twob)));
                                    p=(p-(twob*PBC*lroundl(p/twob)));
                                    q=(q-(twob*PBC*lroundl(q/twob)));
                                    r=(r-(twob*PBC*lroundl(r/twob)));
                                    Sx=(Sx-(tilt*PBC*lroundl(Sy/twob)));
                                    Sx=(Sx-(twob*PBC*lroundl(Sx/twob)));
                                    Sy=(Sy-(twob*PBC*lroundl(Sy/twob)));
                                    Sz=(Sz-(twob*PBC*lroundl(Sz/twob)));
                                    Px=(Px-(tilt*PBC*lroundl(Py/twob)));
                                    Px=(Px-(twob*PBC*lroundl(Px/twob)));
                                    Py=(Py-(twob*PBC*lroundl(Py/twob)));
                                    Pz=(Pz-(twob*PBC*lroundl(Pz/twob)));
                                    int sign;
                                    if(ATOM->part_c[i][j]==1)
                                    {
                                        long double ax=(Sy*r-Sz*q);
                                        long double ay=(Sz*p-Sx*r);
                                        long double az=(Sx*q-Sy*p);
                                        long double overlap=(Px*ax+Py*ay+Pz*az);
                                        if(overlap<0.)
                                            sign=1;
                                        else
                                            sign=-1;
                                        D_TWO=constr_del(ATOM,Atoms,TYPE,p,q,r,Sx,Sy,Sz,sign,i,j,k,rA,rB,D_TWO);
                                        if(D_TWO)
                                        {
                                            long double X1,Y1,Z1;
                                            long double X2,Y2,Z2;

                                            vertice *temp_vert_o=nullptr;
                                            vertice *temp_vert_d=nullptr;

                                            temp_vert_o=D_ONE->v;
                                            temp_vert_d=D_TWO->v;

                                            a=ATOM->p.x;
                                            b=ATOM->p.y;
                                            c=ATOM->p.z;

                                            X1=Atoms[ATOM->contigous[i]].p.x-ATOM->p.x;
                                            Y1=Atoms[ATOM->contigous[i]].p.y-ATOM->p.y;
                                            Z1=Atoms[ATOM->contigous[i]].p.z-ATOM->p.z;

                                            X2=Atoms[ATOM->contigous[j]].p.x-ATOM->p.x;
                                            Y2=Atoms[ATOM->contigous[j]].p.y-ATOM->p.y;
                                            Z2=Atoms[ATOM->contigous[j]].p.z-ATOM->p.z;

                                            long double V1x,V1y,V1z;
                                            long double V2x,V2y,V2z;

                                            V1x=D_ONE->circum_x-a;
                                            V1y=D_ONE->circum_y-b;
                                            V1z=D_ONE->circum_z-c;

                                            V2x=D_TWO->circum_x-a;
                                            V2y=D_TWO->circum_y-b;
                                            V2z=D_TWO->circum_z-c;

                                            V2x=(V2x-(tilt*PBC*lroundl(V2y/twob)));
                                            V2x=(V2x-(twob*PBC*lroundl(V2x/twob)));
                                            V2y=(V2y-(twob*PBC*lroundl(V2y/twob)));
                                            V2z=(V2z-(twob*PBC*lroundl(V2z/twob)));

                                            V1x=(V1x-(tilt*PBC*lroundl(V1y/twob)));
                                            V1x=(V1x-(twob*PBC*lroundl(V1x/twob)));
                                            V1y=(V1y-(twob*PBC*lroundl(V1y/twob)));
                                            V1z=(V1z-(twob*PBC*lroundl(V1z/twob)));

                                            X1=(X1-(tilt*PBC*lroundl(Y1/twob)));
                                            X1=(X1-(twob*PBC*lroundl(X1/twob)));
                                            Y1=(Y1-(twob*PBC*lroundl(Y1/twob)));
                                            Z1=(Z1-(twob*PBC*lroundl(Z1/twob)));

                                            X2=(X2-(tilt*PBC*lroundl(Y2/twob)));
                                            X2=(X2-(twob*PBC*lroundl(X2/twob)));
                                            Y2=(Y2-(twob*PBC*lroundl(Y2/twob)));
                                            Z2=(Z2-(twob*PBC*lroundl(Z2/twob)));

                                            //long double a,b,c;
                                            long double rS=ATOM->radius+r_cut;

                                            a=(Y1*Z2-Z1*Y2);
                                            b=(Z1*X2-X1*Z2);
                                            c=(X1*Y2-Y1*X2);

                                            long double overlap1,overlap2;
                                            int sign1,sign2;
                                            long double disV1,disV2,disA;

                                            disV1=sqrtl((V1x*V1x+V1y*V1y+V1z*V1z));
                                            if((disV1*disV1-rS*rS)>0.)
                                            {
                                                temp_vert_o->is_void=1;
                                                temp_vert_o->ball_r=sqrtl(disV1*disV1-rS*rS);
                                            }
                                            disV2=sqrtl((V2x*V2x+V2y*V2y+V2z*V2z));
                                            if((disV2*disV2-rS*rS)>0.)
                                            {
                                                temp_vert_d->is_void=1;
                                                temp_vert_d->ball_r=sqrtl(disV2*disV2-rS*rS);
                                            }

                                            disA=sqrtl((a*a+b*b+c*c));//V1x*V1x+V1y*V1y+V1z*V1z));
                                            overlap1=a*V1x+b*V1y+c*V1z;
                                            if(overlap1<0.)
                                                sign1=1;
                                            else
                                                sign1=-1;
                                            overlap2=a*V2x+b*V2y+c*V2z;
                                            if(overlap2<0.)
                                                sign2=1;
                                            else
                                                sign2=-1;
                                            if(sign1==sign2)
                                            {
                                                long double disv1,disv2;
                                                disv1=sqrtl(distancesq(ATOM->MIDP[i][j],*temp_vert_o->p));
                                                disv2=sqrtl(distancesq(ATOM->MIDP[i][j],*temp_vert_d->p));
                                                if(disv1<disv2)
                                                {
                                                    if(temp_vert_o->is_void==1)
                                                    {

                                                        ATOM->D3bondinvoid[i][j]=1;
                                                        ATOM->D3bondinvoid[j][i]=1;
                                                    }
                                                    ////cout<<"draw sphere\t{";
                                                    ////cout<<temp_vert_o->p->x<<"\t"<<temp_vert_o->p->y<<"\t"<<temp_vert_o->p->z<<"} radius 0.3\t resolution 10\n";
                                                }
                                                else
                                                {
                                                    if(temp_vert_d->is_void==1)
                                                    {
                                                        ATOM->D3bondinvoid[i][j]=1;
                                                        ATOM->D3bondinvoid[j][i]=1;
                                                    }
                                                    ////cout<<"draw sphere\t{";
                                                    ////cout<<temp_vert_d->p->x<<"\t"<<temp_vert_d->p->y<<"\t"<<temp_vert_d->p->z<<"} radius 0.3\t resolution 10\n";
                                                }

                                                ////cout<<"sameside\n";
                                                ////if(overlap1/(disV1*disA)<overlap2/(disV2*disA))
                                                ////{
                                                ////	if(temp_vert_o->is_void==1)
                                                ////	{
                                                ////
                                                ////		ATOM->D3bondinvoid[i][j]=1;
                                                ////		ATOM->D3bondinvoid[j][i]=1;
                                                ////	}
                                                ////}
                                                ////else
                                                ////{
                                                ////	if(temp_vert_d->is_void==1)
                                                ////	{
                                                ////		ATOM->D3bondinvoid[i][j]=1;
                                                ////		ATOM->D3bondinvoid[j][i]=1;
                                                ////	}
                                                ////}
                                            }
                                            else
                                            {
                                                long double midx,midy,midz;
                                                int p,q,r;
                                                //cout<<ATOM->MIDP[i][j].x<<"\t"<<ATOM->MIDP[i][j].y<<"\t"<<ATOM->MIDP[i][j].z<<"\n";
                                                midx=ATOM->MIDP[i][j].x-ATOM->p.x;
                                                midy=ATOM->MIDP[i][j].y-ATOM->p.y;
                                                midz=ATOM->MIDP[i][j].z-ATOM->p.z;
                                                midx=(midx-(tilt*PBC*lroundl(midy/twob)));
                                                midx=(midx-(twob*PBC*lroundl(midx/twob)));
                                                midy=(midy-(twob*PBC*lroundl(midy/twob)));
                                                midz=(midz-(twob*PBC*lroundl(midz/twob)));

                                                long double dismsq=midx*midx+midy*midy+midz*midz;
                                                if((dismsq-rS*rS)>0.)
                                                {
                                                    ATOM->D3bondinvoid[i][j]=1;
                                                    ATOM->D3bondinvoid[j][i]=1;
                                                }
                                            }
                                            atom* Ai;
                                            atom* Aj;
                                            atom* Ak;

                                            Ai=&(Atoms[ATOM->contigous[i]]);
                                            Aj=&(Atoms[ATOM->contigous[j]]);
                                            Ak=&(Atoms[ATOM->contigous[k]]);
                                            if(temp_vert_o->is_void && !temp_vert_d->is_void)
                                            {
                                                face *f=nullptr;
                                                f= new face;
                                                f->A1=ATOM->p;
                                                f->A2=Ai->p;
                                                f->A3=Aj->p;
                                                f->B=Ak->p;
                                                solid_wall.insert_face(f);

                                            }
                                            if(!temp_vert_o->is_void && temp_vert_d->is_void)
                                            {
                                                face *f=nullptr;
                                                f= new face;
                                                f->A1=ATOM->p;
                                                f->A2=Ai->p;
                                                f->A3=Aj->p;
                                                f->B=Ak->p;
                                                solid_wall.insert_face(f);
                                            }
                                            //if(temp_vert_o->is_void && temp_vert_d->is_void && temp_vert_d->D->hull)
                                            //{
                                            //	face *f=nullptr;
                                            //	f= new face;
                                            //	f->A1=ATOM->p;
                                            //	f->A2=Ai->p;
                                            //	f->A3=Aj->p;
                                            //	f->B=Ak->p;
                                            //	solid_wall.insert_face(f);
                                            //}
                                            ////if(temp_vert_o->is_void && temp_vert_d->is_void)
                                            ////{
                                            ////	cout<<"draw color red\n";
                                            ////	cout<<"draw line\t{";
                                            ////	cout<<temp_vert_o->p->x<<"\t"<<temp_vert_o->p->y<<"\t"<<temp_vert_o->p->z<<"}\t{";
                                            ////	cout<<temp_vert_d->p->x<<"\t"<<temp_vert_d->p->y<<"\t"<<temp_vert_d->p->z<<"}\n";
                                            ////}

                                            add_connected(temp_vert_o,temp_vert_d,ATOM->D3bondinvoid[i][j]);
                                            add_connected(temp_vert_d,temp_vert_o,ATOM->D3bondinvoid[i][j]);
                                            edge *temp_e=nullptr;
                                            temp_e=new edge;
                                            temp_e->origin=temp_vert_o;
                                            temp_e->destination=temp_vert_d;
                                            temp_e->evoid=ATOM->D3bondinvoid[i][j];
                                            Ai=&(Atoms[ATOM->contigous[i]]);
                                            Aj=&(Atoms[ATOM->contigous[j]]);
                                            temp_e->A1=ATOM->index;
                                            temp_e->A2=Ai->index;
                                            temp_e->A3=Aj->index;
                                            ATOM->F[i].insert_element(temp_e);
                                            ATOM->F[j].insert_element(temp_e);
                                            Ai->F[Ai->conti_index[ATOM->index]].insert_element(temp_e);
                                            Aj->F[Aj->conti_index[ATOM->index]].insert_element(temp_e);
                                            Ai->F[Ai->conti_index[Aj->index]].insert_element(temp_e);
                                            Aj->F[Aj->conti_index[Ai->index]].insert_element(temp_e);

                                            if(ATOM->D3bondinvoid[i][j])// && !temp_vert_o->is_void)
                                            {
                                                if(temp_vert_o->is_void!=temp_vert_d->is_void)
                                                    cout<<"####\t"<<temp_vert_o->is_void<<"\t"<<temp_vert_d->is_void<<"\n";
                                                ////cout<<"draw color blue\n";
                                                ////print_delunay(D_ONE,Atoms,nAtoms);
                                                ////print_delunay(D_TWO,Atoms,nAtoms);
                                                ////cout<<"draw color red\n";
                                                ////cout<<"draw line\t{";
                                                ////cout<<temp_vert_o->p->x<<"\t"<<temp_vert_o->p->y<<"\t"<<temp_vert_o->p->z<<"}\t{";
                                                ////cout<<temp_vert_d->p->x<<"\t"<<temp_vert_d->p->y<<"\t"<<temp_vert_d->p->z<<"}\n";
                                            }

                                            ////atom* Ai;
                                            ////atom* Aj;
                                            ATOM->bond_c[i][j]=1;
                                            Ai->bond_c[Ai->conti_index[ATOM->index]][Ai->conti_index[Aj->index]]=1;
                                            Aj->bond_c[Aj->conti_index[ATOM->index]][Aj->conti_index[Ai->index]]=1;

                                            ATOM->bond_c[j][i]=1;
                                            Ai->bond_c[Ai->conti_index[Aj->index]][Ai->conti_index[ATOM->index]]=1;
                                            Aj->bond_c[Aj->conti_index[Ai->index]][Aj->conti_index[ATOM->index]]=1;
                                        }
                                        else
                                        {
                                            atom* Ai;
                                            atom* Aj;
                                            atom* Ak;

                                            Ai=&(Atoms[ATOM->contigous[i]]);
                                            Aj=&(Atoms[ATOM->contigous[j]]);
                                            Ak=&(Atoms[ATOM->contigous[k]]);

                                            ATOM->part_c[i][j]++;
                                            Ai->part_c[Ai->conti_index[ATOM->index]][Ai->conti_index[Aj->index]]++;
                                            Aj->part_c[Aj->conti_index[ATOM->index]][Aj->conti_index[Ai->index]]++;

                                            ATOM->part_c[j][i]++;
                                            Ai->part_c[Ai->conti_index[Aj->index]][Ai->conti_index[ATOM->index]]++;
                                            Aj->part_c[Aj->conti_index[Ai->index]][Aj->conti_index[ATOM->index]]++;

                                            ATOM->bond_c[i][j]=1;
                                            Ai->bond_c[Ai->conti_index[ATOM->index]][Ai->conti_index[Aj->index]]=1;
                                            Aj->bond_c[Aj->conti_index[ATOM->index]][Aj->conti_index[Ai->index]]=1;

                                            ATOM->bond_c[j][i]=1;
                                            Ai->bond_c[Ai->conti_index[Aj->index]][Ai->conti_index[ATOM->index]]=1;
                                            Aj->bond_c[Aj->conti_index[Ai->index]][Aj->conti_index[ATOM->index]]=1;
                                            face *f=nullptr;
                                            face *f_sw=nullptr;
                                            f= new face;
                                            f->A1=ATOM->p;
                                            f->A2=Ai->p;
                                            f->A3=Aj->p;
                                            f_sw= new face;
                                            f_sw->A1=ATOM->p;
                                            f_sw->A2=Ai->p;
                                            f_sw->A3=Aj->p;
                                            f->B=Ak->p;
                                            if(!ATOM->bor)
                                            {
                                                ATOM->bor=1;
                                            }
                                            if(!Ai->bor)
                                            {
                                                Ai->bor=1;
                                            }
                                            if(!Aj->bor)
                                            {
                                                Aj->bor=1;
                                            }
                                            CH.insert_face(f);
                                            if(D_ONE->v->is_void==0)
                                            {
                                                solid_wall.insert_face(f_sw);
                                            }
                                            else
                                            {
                                                //print_face(f,1);
                                            }
                                            ////if(!convex_hull)
                                            ////{
                                            ////	convex_hull=new face;
                                            ////	convex_hull->A1=
                                            ////}
                                            D_ONE->hull=1;
                                        }
                                    }
                                    break;
                                }
                                else if (temp->next)
                                {
                                    temp=temp->next;
                                }
                                else
                                    break;
                            }
                        }
                    }
                }
                else if (ATOM->part_c[i][j]==2)
                {
                    //cout<<ATOM->index<<"\t"<<ATOM->contigous[i]<<"\t"<<ATOM->contigous[j]<<"\n";
                    if(!ATOM->bond_c[i][j])
                    {
                        //break;
                        ////	Atoms[SAM].MIDP[i][j]=center_of_triangle(ATOM->p,Atoms[ATOM->contigous[i]].p,Atoms[ATOM->contigous[j]].p,ATOM->radius+r_cut,Atoms[ATOM->contigous[i]].radius+r_cut,Atoms[ATOM->contigous[j]].radius+r_cut);
                        //cout<<ATOM->MIDP[i][j].x<<"\t"<<ATOM->MIDP[i][j].y<<"\t"<<ATOM->MIDP[i][j].z<<"\n";
                        delunay *D_ONE=nullptr,*D_TWO=nullptr;
                        for(int k=0; k<ATOM->conti; k++)
                        {
                            if(D_ONE && D_TWO)
                            {
                                break;
                            }
                            if(k!=i && k!=j)
                            {
                                container_delunay *temp;
                                temp=ATOM->D_FIRST;
                                //int flag=1;
                                //cout<<i<<"\t"<<j<<"\t"<<k<<"\n";
                                while(1)
                                {
                                    std::vector<int> EV1 {ATOM->index,ATOM->contigous[i],ATOM->contigous[j],ATOM->contigous[k]};//,A3,A4;
                                    std::sort(EV1.begin(),EV1.end());
                                    if(temp->D->AT[0]==EV1[0] && temp->D->AT[1]==EV1[1] && temp->D->AT[2]==EV1[2] && temp->D->AT[3]==EV1[3])
                                    {
                                        if(!D_ONE)
                                        {
                                            D_ONE=temp->D;
                                            break;
                                        }
                                        else
                                        {
                                            D_TWO=temp->D;
                                            break;
                                        }

                                    }
                                    if(temp->next)
                                        temp=temp->next;
                                    else
                                        break;
                                }
                            }
                        }
                        if(D_ONE && D_TWO)
                        {
                            long double V1x,V1y,V1z;
                            long double V2x,V2y,V2z;
                            long double X1,Y1,Z1;
                            long double X2,Y2,Z2;

                            a=ATOM->p.x;
                            b=ATOM->p.y;
                            c=ATOM->p.z;


                            V1x=D_ONE->circum_x-a;
                            V1y=D_ONE->circum_y-b;
                            V1z=D_ONE->circum_z-c;

                            V2x=D_TWO->circum_x-a;
                            V2y=D_TWO->circum_y-b;
                            V2z=D_TWO->circum_z-c;

                            X1=Atoms[ATOM->contigous[i]].p.x-ATOM->p.x;
                            Y1=Atoms[ATOM->contigous[i]].p.y-ATOM->p.y;
                            Z1=Atoms[ATOM->contigous[i]].p.z-ATOM->p.z;

                            X2=Atoms[ATOM->contigous[j]].p.x-ATOM->p.x;
                            Y2=Atoms[ATOM->contigous[j]].p.y-ATOM->p.y;
                            Z2=Atoms[ATOM->contigous[j]].p.z-ATOM->p.z;

                            V2x=(V2x-(tilt*PBC*lroundl(V2y/twob)));
                            V2x=(V2x-(twob*PBC*lroundl(V2x/twob)));
                            V2y=(V2y-(twob*PBC*lroundl(V2y/twob)));
                            V2z=(V2z-(twob*PBC*lroundl(V2z/twob)));

                            V1x=(V1x-(tilt*PBC*lroundl(V1y/twob)));
                            V1x=(V1x-(twob*PBC*lroundl(V1x/twob)));
                            V1y=(V1y-(twob*PBC*lroundl(V1y/twob)));
                            V1z=(V1z-(twob*PBC*lroundl(V1z/twob)));

                            X1=(X1-(tilt*PBC*lroundl(Y1/twob)));
                            X1=(X1-(twob*PBC*lroundl(X1/twob)));
                            Y1=(Y1-(twob*PBC*lroundl(Y1/twob)));
                            Z1=(Z1-(twob*PBC*lroundl(Z1/twob)));

                            X2=(X2-(tilt*PBC*lroundl(Y2/twob)));
                            X2=(X2-(twob*PBC*lroundl(X2/twob)));
                            Y2=(Y2-(twob*PBC*lroundl(Y2/twob)));
                            Z2=(Z2-(twob*PBC*lroundl(Z2/twob)));

                            //long double a,b,c;
                            long double rS=ATOM->radius+r_cut;

                            vertice *temp_vert_o;
                            vertice *temp_vert_d;
                            temp_vert_o=D_ONE->v;
                            temp_vert_d=D_TWO->v;

                            a=(Y1*Z2-Z1*Y2);
                            b=(Z1*X2-X1*Z2);
                            c=(X1*Y2-Y1*X2);

                            long double overlap1,overlap2;
                            int sign1,sign2;
                            long double disV1,disV2,disA;

                            disV1=sqrtl((V1x*V1x+V1y*V1y+V1z*V1z));
                            disV2=sqrtl((V2x*V2x+V2y*V2y+V2z*V2z));

                            disA=sqrtl((a*a+b*b+c*c));//V1x*V1x+V1y*V1y+V1z*V1z));
                            overlap1=a*V1x+b*V1y+c*V1z;
                            if(overlap1<0.)
                                sign1=1;
                            else
                                sign1=-1;
                            overlap2=a*V2x+b*V2y+c*V2z;
                            if(overlap2<0.)
                                sign2=1;
                            else
                                sign2=-1;
                            if(sign1==sign2)
                            {
                                long double disv1,disv2;
                                disv1=sqrtl(distancesq(ATOM->MIDP[i][j],*temp_vert_o->p));
                                disv2=sqrtl(distancesq(ATOM->MIDP[i][j],*temp_vert_d->p));
                                if(disv1<disv2)
                                {
                                    if(temp_vert_o->is_void==1)
                                    {

                                        ATOM->D3bondinvoid[i][j]=1;
                                        ATOM->D3bondinvoid[j][i]=1;
                                    }
                                    ////cout<<"draw sphere\t{";
                                    ////cout<<temp_vert_o->p->x<<"\t"<<temp_vert_o->p->y<<"\t"<<temp_vert_o->p->z<<"} radius 0.3\t resolution 10\n";
                                }
                                else
                                {
                                    if(temp_vert_d->is_void==1)
                                    {
                                        ATOM->D3bondinvoid[i][j]=1;
                                        ATOM->D3bondinvoid[j][i]=1;
                                    }
                                    ////cout<<"draw sphere\t{";
                                    ////cout<<temp_vert_d->p->x<<"\t"<<temp_vert_d->p->y<<"\t"<<temp_vert_d->p->z<<"} radius 0.3\t resolution 10\n";
                                }
                            }
                            else
                            {
                                long double midx,midy,midz;
                                int p,q,r;
                                midx=ATOM->MIDP[i][j].x-ATOM->p.x;
                                midy=ATOM->MIDP[i][j].y-ATOM->p.y;
                                midz=ATOM->MIDP[i][j].z-ATOM->p.z;
                                midx=(midx-(tilt*PBC*lroundl(midy/twob)));
                                midx=(midx-(twob*PBC*lroundl(midx/twob)));
                                midy=(midy-(twob*PBC*lroundl(midy/twob)));
                                midz=(midz-(twob*PBC*lroundl(midz/twob)));
                                long double dismsq=midx*midx+midy*midy+midz*midz;
                                ////cout<<ATOM->index<<"\n";
                                ////cout<<"here\t"<<sqrtl(dismsq)<<"\t"<<rS<<"\n";
                                if((dismsq-rS*rS)>0.)
                                {

                                    ATOM->D3bondinvoid[i][j]=1;
                                    ATOM->D3bondinvoid[j][i]=1;
                                }
                            }
                            atom* Ai;
                            atom* Aj;
                            atom* Ak;

                            Ai=&(Atoms[ATOM->contigous[i]]);
                            Aj=&(Atoms[ATOM->contigous[j]]);
                            //Ak=&(Atoms[ATOM->contigous[k]]);
                            if(temp_vert_o->is_void && !temp_vert_d->is_void)
                            {
                                face *f=nullptr;
                                f= new face;
                                f->A1=ATOM->p;
                                f->A2=Ai->p;
                                f->A3=Aj->p;
                                //f->B=Ak->p;
                                solid_wall.insert_face(f);

                            }
                            if(!temp_vert_o->is_void && temp_vert_d->is_void)
                            {
                                face *f=nullptr;
                                f= new face;
                                f->A1=ATOM->p;
                                f->A2=Ai->p;
                                f->A3=Aj->p;
                                //f->B=Ak->p;
                                solid_wall.insert_face(f);
                            }
                            ////if(temp_vert_o->is_void && temp_vert_d->is_void && ATOM->D3bondinvoid[i][j]==0)
                            ////{
                            ////	face *f=nullptr;
                            ////	f= new face;
                            ////	f->A1=ATOM->p;
                            ////	f->A2=Ai->p;
                            ////	f->A3=Aj->p;
                            ////	//f->B=Ak->p;
                            ////	solid_wall.insert_face(f);
                            ////}
                            //if(D_ONE->AT[0]==60 && D_ONE->AT[1]==128 && D_ONE->AT[2]==146 && D_ONE->AT[3]==150)
                            //{
                            //	if(D_TWO->AT[0]==62 && D_TWO->AT[1]==128 && D_TWO->AT[2]==146 && D_TWO->AT[3]==150)
                            //	{
                            //        cout<<ATOM->D3bondinvoid[i][j]<<"\n";
                            //		cout<<"draw color blue\n";
                            //		cout<<"draw sphere\t{";
                            //		cout<<ATOM->MIDP[i][j].x<<"\t"<<ATOM->MIDP[i][j].y<<"\t"<<ATOM->MIDP[i][j].z<<"} radius 0.1\n";

                            //	}
                            //}
                            ////if(temp_vert_o->is_void && !temp_vert_d->is_void)
                            ////{
                            ////	cout<<"draw color red\n";
                            ////	cout<<"draw line\t{";
                            ////	cout<<temp_vert_o->p->x<<"\t"<<temp_vert_o->p->y<<"\t"<<temp_vert_o->p->z<<"}\t{";
                            ////	cout<<temp_vert_d->p->x<<"\t"<<temp_vert_d->p->y<<"\t"<<temp_vert_d->p->z<<"}\n";
                            ////}
                            if(ATOM->D3bondinvoid[i][j])// && !temp_vert_o->is_void)
                            {
                                if(temp_vert_o->is_void!=temp_vert_d->is_void)
                                    cout<<"####\t"<<temp_vert_o->is_void<<"\t"<<temp_vert_d->is_void<<"\n";
                                ////cout<<"draw color blue\n";
                                ////print_delunay(D_ONE,Atoms,nAtoms);
                                ////print_delunay(D_TWO,Atoms,nAtoms);
                                ////cout<<"draw color red\n";
                                ////cout<<"draw line\t{";
                                ////cout<<temp_vert_o->p->x<<"\t"<<temp_vert_o->p->y<<"\t"<<temp_vert_o->p->z<<"}\t{";
                                ////cout<<temp_vert_d->p->x<<"\t"<<temp_vert_d->p->y<<"\t"<<temp_vert_d->p->z<<"}\n";
                            }
                            add_connected(temp_vert_o,temp_vert_d,ATOM->D3bondinvoid[i][j]);
                            add_connected(temp_vert_d,temp_vert_o,ATOM->D3bondinvoid[i][j]);
                            edge *temp_e=nullptr;
                            temp_e=new edge;
                            temp_e->origin=temp_vert_o;
                            temp_e->destination=temp_vert_d;
                            temp_e->evoid=ATOM->D3bondinvoid[i][j];
                            Ai=&(Atoms[ATOM->contigous[i]]);
                            Aj=&(Atoms[ATOM->contigous[j]]);
                            temp_e->A1=ATOM->index;
                            temp_e->A2=Ai->index;
                            temp_e->A3=Aj->index;
                            ATOM->F[i].insert_element(temp_e);
                            ATOM->F[j].insert_element(temp_e);
                            Ai->F[Ai->conti_index[ATOM->index]].insert_element(temp_e);
                            Aj->F[Aj->conti_index[ATOM->index]].insert_element(temp_e);
                            Ai->F[Ai->conti_index[Aj->index]].insert_element(temp_e);
                            Aj->F[Aj->conti_index[Ai->index]].insert_element(temp_e);

                            //atom* Ai;
                            //atom* Aj;

                            ATOM->bond_c[i][j]=1;
                            Ai->bond_c[Ai->conti_index[ATOM->index]][Ai->conti_index[Aj->index]]=1;
                            Aj->bond_c[Aj->conti_index[ATOM->index]][Aj->conti_index[Ai->index]]=1;

                            ATOM->bond_c[j][i]=1;
                            Ai->bond_c[Ai->conti_index[Aj->index]][Ai->conti_index[ATOM->index]]=1;
                            Aj->bond_c[Aj->conti_index[Ai->index]][Aj->conti_index[ATOM->index]]=1;
                        }
                    }
                }
            }
        }
        s=0;
        for(int p=0; p<ATOM->conti; p++)
        {
            s=s+ATOM->part_c[i][p];
        }
    }
}

long double volume_tetrahedron(long double Ax,long double Ay,long double Az,long double Ex,long double Ey,long double Ez,long double Bx,long double By,long double Bz,long double vx,long double vy,long double vz,long double r,int compliment=0,int debug=0)
{
    long double ABx,ABy,ABz,EBx,EBy,EBz,VEx,VEy,VEz,DISB,DISBE,DISVE;

    ABx=Ax-Bx;
    ABy=Ay-By;
    ABz=Az-Bz;
    EBx=Ex-Bx;
    EBy=Ey-By;
    EBz=Ez-Bz;
    VEx=Ex-vx;
    VEy=Ey-vy;
    VEz=Ez-vz;

    ABx=(ABx-(tilt*PBC*lroundl(ABy/twob)));
    ABx=(ABx-(twob*PBC*lroundl(ABx/twob)));
    ABy=(ABy-(twob*PBC*lroundl(ABy/twob)));
    ABz=(ABz-(twob*PBC*lroundl(ABz/twob)));
    EBx=(EBx-(tilt*PBC*lroundl(EBy/twob)));
    EBx=(EBx-(twob*PBC*lroundl(EBx/twob)));
    EBy=(EBy-(twob*PBC*lroundl(EBy/twob)));
    EBz=(EBz-(twob*PBC*lroundl(EBz/twob)));
    VEx=(VEx-(tilt*PBC*lroundl(VEy/twob)));
    VEx=(VEx-(twob*PBC*lroundl(VEx/twob)));
    VEy=(VEy-(twob*PBC*lroundl(VEy/twob)));
    VEz=(VEz-(twob*PBC*lroundl(VEz/twob)));

    DISB=ABx*ABx+ABy*ABy+ABz*ABz;
    DISBE=EBx*EBx+EBy*EBy+EBz*EBz;
    DISVE=VEx*VEx+VEy*VEy+VEz*VEz;
    Ax=0.;
    Ay=0.;
    Az=0.;
    Bx=sqrtl(DISB);
    By=0.;
    Bz=0.;
    Ex=sqrtl(DISB);
    Ey=sqrtl(DISBE);
    Ez=0.;
	site* ABs;
	site* AEs;
	site* AVs;
	site* BVs;
	site* BEs;
	site* VEs;
    if(debug)
    {
        cout<<"###\n";
        cout<<"mol new\n";
        cout<<"draw material Transparent\n";
        cout<<"draw sphere\t{";
        cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\tradius\t"<<r<<"\t"<<"resolution 10\n";
        ////cout<<"draw sphere\t{";
        ////cout<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\tradius 0.3\n";
        ////cout<<"draw sphere\t{";
        ////cout<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\tradius 0.3\n";
        cout<<"mol new\n";
        cout<<"draw material Opaque\n";
        cout<<"draw color red\n";
        cout<<"draw line\t{";
        cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\n";
        cout<<"draw line\t{";
        cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\n";
        cout<<"draw line\t{";
        cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<vx<<"\t"<<vy<<"\t"<<vz<<"}\n";
        cout<<"draw line\t{";
        cout<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\t{"<<vx<<"\t"<<vy<<"\t"<<vz<<"}\n";
        cout<<"draw line\t{";
        cout<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\t{"<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\n";
        cout<<"draw line\t{";
        cout<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\t{"<<vx<<"\t"<<vy<<"\t"<<vz<<"}\n";
        cout<<"draw sphere\t{";
        cout<<vx<<"\t"<<vy<<"\t"<<vz<<"}\tradius 0.3\n";
        cout<<"###\n";
    }
    vx=sqrtl(DISB);
    vy=sqrtl(DISBE);
    vz=sqrtl(DISVE);
    long double x0,y0,z0;
    long double x0sq,y0sq,z0sq;
    x0=Bx;
    y0=Ey;
    z0=vz;
    x0sq=DISB;
    y0sq=DISBE;
    z0sq=DISVE;
    long double rB,rV,rE;
    rV=sqrtl(vx*vx+vy*vy+vz*vz);
    rE=sqrtl(Ex*Ex+Ey*Ey+Ez*Ez);
    rB=sqrtl(DISB);
    long double theta,x2,y2;
    x2=r*x0/rE;
    y2=r*y0/rE;
    theta=atanl(z0/y0);
    long double Vc,Vt;
    Vt=(x0*y0*z0)/6.;
    Vc=0.;
    if(r<rB)
    {
        Vc=(r*r*r/6.)*(2*theta-M_PI/2.-asinl((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq))));
        //cout<<Vc<<"\n";
    }
    if(rB<r && r<rE)
    {
        Vc=theta/2.*(r*r*x0-x0*x0*x0/3.)-(r*r*r/6.)*(M_PI/2.+asinl((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq))));
        //cout<<Vc<<"\n";
    }
    if(rE<r && r<rV)
    {
        Vc=0.5*(theta-M_PI/2.+asin(y0/sqrtl(r*r-x0sq)))*(r*r*x0-x0*x0*x0/3.)+x0*y0/6.*sqrtl(r*r-rE*rE)+r*r*r/6.*asinl((x2*x2-y2*y2-x0sq)/(r*r-x0sq))-r*r*r/6.*asinl((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq)));
        //cout<<Vc<<"\n";
    }
    if(Vt<Vc)
    {
        //cout<<"help \t"<<Vt<<"\t"<<Vc<<"\t"<<Vt-Vc<<"\n";
    }
    if(!compliment)
        return Vt-Vc;
    else
    {
        //cout<<"###\t"<<Vt<<"\t"<<Vc<<"\t"<<r<<"\t"<<rB<<"\t"<<rE<<"\t"<<rV<<"\n";
        if(Vc==0.)
            return Vt;
        else
            return Vc;
    }
}
long double volume_segment(site *a1,site *a2,site *a3,long double r)
{
	long double theta=acosl((a1->x*a2->x+a1->y*a2->y+a1->z*a2->z)/(r*r));		
	long double phi=acosl((a1->x*a3->x+a1->y*a3->y+a1->z*a3->z)/(r*r));		
	long double rm=r-r_cut;
	//cout<<theta<<"\t"<<phi<<"\n";
	return r*r*r/6.*theta*phi-rm*rm*rm/6.*theta*phi;
}
long double volume_of_roling_probe(long double ra,long double theta,long double alpha)
{
	return (2.*r_cut*r_cut*r_cut*sinl(alpha)*theta/3.)/2.;	
}
long double removesurface(long double Ax,long double Ay,long double Az,long double Ex,long double Ey,long double Ez,long double Bx,long double By,long double Bz,long double vx,long double vy,long double vz,long double r,int compliment=0,int debug=0)
{
    long double ABx,ABy,ABz,EBx,EBy,EBz,VEx,VEy,VEz,DISB,DISBE,DISVE;

    ABx=Ax-Bx;
    ABy=Ay-By;
    ABz=Az-Bz;
    EBx=Ex-Bx;
    EBy=Ey-By;
    EBz=Ez-Bz;
    VEx=Ex-vx;
    VEy=Ey-vy;
    VEz=Ez-vz;

    ABx=(ABx-(tilt*PBC*lroundl(ABy/twob)));
    ABx=(ABx-(twob*PBC*lroundl(ABx/twob)));
    ABy=(ABy-(twob*PBC*lroundl(ABy/twob)));
    ABz=(ABz-(twob*PBC*lroundl(ABz/twob)));
    EBx=(EBx-(tilt*PBC*lroundl(EBy/twob)));
    EBx=(EBx-(twob*PBC*lroundl(EBx/twob)));
    EBy=(EBy-(twob*PBC*lroundl(EBy/twob)));
    EBz=(EBz-(twob*PBC*lroundl(EBz/twob)));
    VEx=(VEx-(tilt*PBC*lroundl(VEy/twob)));
    VEx=(VEx-(twob*PBC*lroundl(VEx/twob)));
    VEy=(VEy-(twob*PBC*lroundl(VEy/twob)));
    VEz=(VEz-(twob*PBC*lroundl(VEz/twob)));

    DISB=ABx*ABx+ABy*ABy+ABz*ABz;
    DISBE=EBx*EBx+EBy*EBy+EBz*EBz;
    DISVE=VEx*VEx+VEy*VEy+VEz*VEz;
    Ax=0.;
    Ay=0.;
    Az=0.;
    Bx=sqrtl(DISB);
    By=0.;
    Bz=0.;
    Ex=sqrtl(DISB);
    Ey=sqrtl(DISBE);
    Ez=0.;
    vx=sqrtl(DISB);
    vy=sqrtl(DISBE);
    vz=sqrtl(DISVE);
    long double x0,y0,z0;
    long double x0sq,y0sq,z0sq;
    x0=Bx;
    y0=Ey;
    z0=vz;
    x0sq=DISB;
    y0sq=DISBE;
    z0sq=DISVE;
    long double rB,rV,rE;
    rV=sqrtl(vx*vx+vy*vy+vz*vz);
    rE=sqrtl(Ex*Ex+Ey*Ey+Ez*Ez);
    rB=sqrtl(DISB);
    long double theta,x2,y2;
    x2=r*x0/rE;
    y2=r*y0/rE;
    theta=atanl(z0/y0);
    long double Vc,Vt;
    Vt=(x0*y0*z0)/6.;
    Vc=0.;
	site* ABs;
	site* AEs;
	site* AVs;
	site* BVs;
	site* BEs;
	site* EVs;
    if(r<rB)
    {
        Vc=Vc+(r*r*r/6.)*(2*theta-M_PI/2.-asinl((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq))));
		r=r-r_cut;
        Vc=Vc-(r*r*r/6.)*(2*theta-M_PI/2.-asinl((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq))));
        //cout<<Vc<<"\n";
    }
    if(rB<r && r<rE)
    {
		long double ra=sqrtl(r*r-rB*rB);
		long double l=r/rV;
		AVs=new site;	
		BVs=new site;	
		BEs=new site;	
		AEs=new site;	
		AVs->x=vx*l;
		AVs->y=vy*l;
		AVs->z=vz*l;
		BVs->x=x0;
		BVs->y=ra*cosl(theta);
		BVs->z=ra*sinl(theta);
		BEs->x=x0;
		BEs->y=ra;
		BEs->z=0.;
		l=r/rE;
		AEs->x=Ex*l;
		AEs->y=Ey*l;
		AEs->z=Ez*l;
		long double alpha=M_PI/2.-acosl((BEs->x*Bx+BEs->y*By+BEs->z*Bz)/(r*rB));		
		alpha=abs(alpha);
		long double V_rp=volume_of_roling_probe(ra,theta,alpha);	
		Vc=Vc+volume_segment(AVs,BVs,BEs,r);
		Vc=Vc+volume_segment(AVs,AEs,BEs,r);
		Vc=Vc+V_rp;
        if(debug)
        {
        	cout<<"###\n";
        	cout<<"mol new\n";
        	cout<<"draw material Transparent\n";
        	cout<<"draw sphere\t{";
        	cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\tradius\t"<<r<<"\t"<<"resolution 25\n";
        	////cout<<"draw sphere\t{";
        	////cout<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\tradius 0.3\n";
        	////cout<<"draw sphere\t{";
        	////cout<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\tradius 0.3\n";
        	cout<<"mol new\n";
        	cout<<"draw material Opaque\n";
        	cout<<"draw color red\n";
        	cout<<"draw line\t{";
        	cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\n";
        	cout<<"draw line\t{";
        	cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\n";
        	cout<<"draw line\t{";
        	cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<vx<<"\t"<<vy<<"\t"<<vz<<"}\n";
        	cout<<"draw line\t{";
        	cout<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\t{"<<vx<<"\t"<<vy<<"\t"<<vz<<"}\n";
        	cout<<"draw line\t{";
        	cout<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\t{"<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\n";
        	cout<<"draw line\t{";
        	cout<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\t{"<<vx<<"\t"<<vy<<"\t"<<vz<<"}\n";
        	cout<<"draw sphere\t{";
        	cout<<vx<<"\t"<<vy<<"\t"<<vz<<"}\tradius 0.3\n";
        	cout<<"###\n";
        	cout<<"draw color green\n";
        	display_SITE(AVs);
        	display_SITE(BVs);
        	display_SITE(BEs);
        	display_SITE(AEs);
        }
	
        //Vc=theta/2.*(r*r*x0-x0*x0*x0/3.)-(r*r*r/6.)*(M_PI/2.+asinl((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq))));
        //cout<<Vc<<"\n";
    }
    if(rE<r && r<rV)
    {
		long double ra=sqrtl(r*r-rB*rB);
		long double l=r/rV;
		AVs=new site;	
		BVs=new site;	
		EVs=new site;	
		AVs->x=vx*l;
		AVs->y=vy*l;
		AVs->z=vz*l;
		BVs->x=x0;
		BVs->y=ra*cosl(theta);
		BVs->z=ra*sinl(theta);
		long double theta_1=acosl(Ey/ra);
		theta_1=abs(theta_1);
		EVs->x=Ex;
		EVs->y=Ey;
		EVs->z=ra*sinl(theta_1);
        if(debug)
        {
        	cout<<"###\n";
        	cout<<"mol new\n";
        	cout<<"draw material Transparent\n";
        	cout<<"draw sphere\t{";
        	cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\tradius\t"<<r<<"\t"<<"resolution 25\n";
        	////cout<<"draw sphere\t{";
        	////cout<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\tradius 0.3\n";
        	////cout<<"draw sphere\t{";
        	////cout<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\tradius 0.3\n";
        	cout<<"mol new\n";
        	cout<<"draw material Opaque\n";
        	cout<<"draw color red\n";
        	cout<<"draw line\t{";
        	cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\n";
        	cout<<"draw line\t{";
        	cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\n";
        	cout<<"draw line\t{";
        	cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<vx<<"\t"<<vy<<"\t"<<vz<<"}\n";
        	cout<<"draw line\t{";
        	cout<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\t{"<<vx<<"\t"<<vy<<"\t"<<vz<<"}\n";
        	cout<<"draw line\t{";
        	cout<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\t{"<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\n";
        	cout<<"draw line\t{";
        	cout<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\t{"<<vx<<"\t"<<vy<<"\t"<<vz<<"}\n";
        	cout<<"draw sphere\t{";
        	cout<<vx<<"\t"<<vy<<"\t"<<vz<<"}\tradius 0.3\n";
        	cout<<"###\n";
        	cout<<"draw color green\n";
        	display_SITE(AVs);
        	display_SITE(BVs);
        	display_SITE(EVs);
        }
		Vc=Vc+volume_segment(AVs,BVs,EVs,r);
		long double alpha=M_PI/2.-acosl((BVs->x*Bx+BVs->y*By+BVs->z*Bz)/(r*rB));		
		alpha=abs(alpha);
		long double V_rp=volume_of_roling_probe(ra,theta-theta_1,alpha);	
		Vc=Vc+r_cut*r_cut*r_cut/6.*alpha*theta_1;
		Vc=Vc+V_rp;
        //Vc=0.5*(theta-M_PI/2.+asin(y0/sqrtl(r*r-x0sq)))*(r*r*x0-x0*x0*x0/3.)+x0*y0/6.*sqrtl(r*r-rE*rE)+r*r*r/6.*asinl((x2*x2-y2*y2-x0sq)/(r*r-x0sq))-r*r*r/6.*asinl((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq)));
        //cout<<Vc<<"\n";
    }
	if(Vc<0.)
		cout<<"error\n";
	return -1.*Vc;
 // if(Vt<Vc)
 // {
 //     //cout<<"help \t"<<Vt<<"\t"<<Vc<<"\t"<<Vt-Vc<<"\n";
 // }
 // if(!compliment)
 //     return Vt-Vc;
 // else
 // {
 //     //cout<<"###\t"<<Vt<<"\t"<<Vc<<"\t"<<r<<"\t"<<rB<<"\t"<<rE<<"\t"<<rV<<"\n";
 //     if(Vc==0.)
 //         return Vt;
 //     else
 //         return Vc;
 // }
}
int sign_aaav(long double A1x,long double A1y,long double A1z,long double A2x,long double A2y,long double A2z,long double A3x,long double A3y,long double A3z,long double A4x,long double A4y,long double A4z,long double Vx,long double Vy,long double Vz)
{
    long double a1x,a1y,a1z;
    long double a2x,a2y,a2z;
    long double a3x,a3y,a3z;
    long double vx,vy,vz;
    a1x=A2x-A1x;
    a1y=A2y-A1y;
    a1z=A2z-A1z;
    a2x=A3x-A1x;
    a2y=A3y-A1y;
    a2z=A3z-A1z;
    a3x=A4x-A1x;
    a3y=A4y-A1y;
    a3z=A4z-A1z;
    vx=Vx-A1x;
    vy=Vy-A1y;
    vz=Vz-A1z;
    a3x=(a3x-(tilt*PBC*lroundl(a3y/twob)));
    a3x=(a3x-(twob*PBC*lroundl(a3x/twob)));
    a3y=(a3y-(twob*PBC*lroundl(a3y/twob)));
    a3z=(a3z-(twob*PBC*lroundl(a3z/twob)));
    vx=(vx-(tilt*PBC*lroundl(vy/twob)));
    vx=(vx-(twob*PBC*lroundl(vx/twob)));
    vy=(vy-(twob*PBC*lroundl(vy/twob)));
    vz=(vz-(twob*PBC*lroundl(vz/twob)));
    a1x=(a1x-(tilt*PBC*lroundl(a1y/twob)));
    a1x=(a1x-(twob*PBC*lroundl(a1x/twob)));
    a1y=(a1y-(twob*PBC*lroundl(a1y/twob)));
    a1z=(a1z-(twob*PBC*lroundl(a1z/twob)));
    a2x=(a2x-(tilt*PBC*lroundl(a2y/twob)));
    a2x=(a2x-(twob*PBC*lroundl(a2x/twob)));
    a2y=(a2y-(twob*PBC*lroundl(a2y/twob)));
    a2z=(a2z-(twob*PBC*lroundl(a2z/twob)));
    long double a,b,c;
    long double overlap1,overlap2;
    int sign1,sign2;
    a=a2z*a1y-a1z*a2y;
    b=a1z*a2x-a1x*a2z;
    c=a2y*a1x-a1y*a2x;
    overlap1=a*a3x+b*a3y+c*a3z;
    if(overlap1<0.)
        sign1=1;
    else
        sign1=-1;
    overlap2=a*vx+b*vy+c*vz;
    if(overlap2<0.)
        sign2=1;
    else
        sign2=-1;
    if(sign1==sign2)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}
int sign_aaae(long double A1x,long double A1y,long double A1z,long double A2x,long double A2y,long double A2z,long double A3x,long double A3y,long double A3z,long double Ex,long double Ey,long double Ez)
{
    long double a1x,a1y,a1z;
    long double a2x,a2y,a2z;
    long double ex,ey,ez;
    a1x=A1x-A2x;
    a1y=A1y-A2y;
    a1z=A1z-A2z;
    a2x=A3x-A2x;
    a2y=A3y-A2y;
    a2z=A3z-A2z;
    ex=Ex-A2x;
    ey=Ey-A2y;
    ez=Ez-A2z;
    a1x=(a1x-(tilt*PBC*lroundl(a1y/twob)));
    a1x=(a1x-(twob*PBC*lroundl(a1x/twob)));
    a1y=(a1y-(twob*PBC*lroundl(a1y/twob)));
    a1z=(a1z-(twob*PBC*lroundl(a1z/twob)));
    a2x=(a2x-(tilt*PBC*lroundl(a2y/twob)));
    a2x=(a2x-(twob*PBC*lroundl(a2x/twob)));
    a2y=(a2y-(twob*PBC*lroundl(a2y/twob)));
    a2z=(a2z-(twob*PBC*lroundl(a2z/twob)));
    ex=(ex-(tilt*PBC*lroundl(ey/twob)));
    ex=(ex-(twob*PBC*lroundl(ex/twob)));
    ey=(ey-(twob*PBC*lroundl(ey/twob)));
    ez=(ez-(twob*PBC*lroundl(ez/twob)));
    long double a,b,c;
    long double a1,b1,c1;
    a=a2z*a1y-a1z*a2y;
    b=a1z*a2x-a1x*a2z;
    c=a2y*a1x-a1y*a2x;
    a1=c*a1y-a1z*b;
    b1=a1z*a-a1x*c;
    c1=b*a1x-a1y*a;
    long double overlap1,overlap2;
    int sign1,sign2;
    overlap1=a1*ex+b1*ey+c1*ez;
    if(overlap1<0.)
    {
        sign1=1;
    }
    else
    {
        sign1=-1;
    }
    overlap2=a1*a2x+b1*a2y+c1*a2z;
    if(overlap2<0.)
    {
        sign2=1;
    }
    else
    {
        sign2=-1;
    }
    if(sign1==sign2)
    {
        return 1;
    }
    else
        return -1;
}
void initialization(atom Atoms[],int nAtoms, int ntypes)
{
    radius=new (nothrow) long double[ntypes];
    contigous = new (nothrow) int* [nAtoms];
    V=new vert_list;
    for(int i=0; i<nAtoms; i++)
    {
        Atoms[i].conti_index = new (nothrow) int  [nAtoms];
        contigous[i] = new (nothrow) int [nAtoms];
    }
}
void delete_everything(atom Atoms[],int nAtoms,int ntypes)
{
    for(int i=0; i<nAtoms; i++)
    {
        delete[] Atoms[i].conti_index;
    }
}
bool inside_delunay(vertice *v,delunay *D,atom Atoms[],int nAtoms)
{
    site V1;
    site V2;
    site V;
    site A;
    int sign1,sign2;
    long double overlap1,overlap2;
    atom *a1=&(Atoms[D->AT[0]]);
    atom *a2=&(Atoms[D->AT[1]]);;
    atom *a3=&(Atoms[D->AT[2]]);;
    atom *a4=&(Atoms[D->AT[3]]);;
    V1.x=a2->p.x-a1->p.x;
    V1.y=a2->p.y-a1->p.y;
    V1.z=a2->p.z-a1->p.z;
    V2.x=a3->p.x-a1->p.x;
    V2.y=a3->p.y-a1->p.y;
    V2.z=a3->p.z-a1->p.z;
    V.x=v->p->x-a1->p.x;
    V.y=v->p->y-a1->p.y;
    V.z=v->p->z-a1->p.z;
    A.x=a4->p.x-a1->p.x;
    A.y=a4->p.y-a1->p.y;
    A.z=a4->p.z-a1->p.z;
    site cross_12;
    cross_12=cross_product(V1,V2);
    overlap1=cross_12.x*A.x+cross_12.y*A.y+cross_12.z*A.z;
    overlap2=cross_12.x*V.x+cross_12.y*V.y+cross_12.z*V.z;
    if(overlap1<0.)
        sign1=1;
    else
        sign1=-1;

    if(overlap2<0.)
        sign2=1;
    else
        sign2=-1;

    if(sign1==sign2)
    {
        V1.x=a2->p.x-a1->p.x;
        V1.y=a2->p.y-a1->p.y;
        V1.z=a2->p.z-a1->p.z;
        V2.x=a4->p.x-a1->p.x;
        V2.y=a4->p.y-a1->p.y;
        V2.z=a4->p.z-a1->p.z;
        V.x=v->p->x-a1->p.x;
        V.y=v->p->y-a1->p.y;
        V.z=v->p->z-a1->p.z;
        A.x=a3->p.x-a1->p.x;
        A.y=a3->p.y-a1->p.y;
        A.z=a3->p.z-a1->p.z;
        cross_12=cross_product(V1,V2);
        overlap1=cross_12.x*A.x+cross_12.y*A.y+cross_12.z*A.z;
        overlap2=cross_12.x*V.x+cross_12.y*V.y+cross_12.z*V.z;
        if(overlap1<0.)
            sign1=1;
        else
            sign1=-1;

        if(overlap2<0.)
            sign2=1;
        else
            sign2=-1;
        if(sign1==sign2)
        {
            V1.x=a3->p.x-a1->p.x;
            V1.y=a3->p.y-a1->p.y;
            V1.z=a3->p.z-a1->p.z;
            V2.x=a4->p.x-a1->p.x;
            V2.y=a4->p.y-a1->p.y;
            V2.z=a4->p.z-a1->p.z;
            V.x=v->p->x-a1->p.x;
            V.y=v->p->y-a1->p.y;
            V.z=v->p->z-a1->p.z;
            A.x=a2->p.x-a1->p.x;
            A.y=a2->p.y-a1->p.y;
            A.z=a2->p.z-a1->p.z;
            cross_12=cross_product(V1,V2);
            overlap1=cross_12.x*A.x+cross_12.y*A.y+cross_12.z*A.z;
            overlap2=cross_12.x*V.x+cross_12.y*V.y+cross_12.z*V.z;
            if(overlap1<0.)
                sign1=1;
            else
                sign1=-1;

            if(overlap2<0.)
                sign2=1;
            else
                sign2=-1;
            if(sign1==sign2)
            {
                V1.x=a3->p.x-a2->p.x;
                V1.y=a3->p.y-a2->p.y;
                V1.z=a3->p.z-a2->p.z;
                V2.x=a4->p.x-a2->p.x;
                V2.y=a4->p.y-a2->p.y;
                V2.z=a4->p.z-a2->p.z;
                V.x=v->p->x-a2->p.x;
                V.y=v->p->y-a2->p.y;
                V.z=v->p->z-a2->p.z;
                A.x=a1->p.x-a2->p.x;
                A.y=a1->p.y-a2->p.y;
                A.z=a1->p.z-a2->p.z;
                cross_12=cross_product(V1,V2);
                overlap1=cross_12.x*A.x+cross_12.y*A.y+cross_12.z*A.z;
                overlap2=cross_12.x*V.x+cross_12.y*V.y+cross_12.z*V.z;
                if(overlap1<0.)
                    sign1=1;
                else
                    sign1=-1;

                if(overlap2<0.)
                    sign2=1;
                else
                    sign2=-1;
                if(sign1==sign2)
                    return true;
                else
                    return false;
            }
            else
                return false;
        }
        else
            return false;
    }
    else
        return false;
}
long double void_vol(vertice *v,atom Atoms[])
{
    long double vol=0.;
    site A1,A2,A3,A4;
    long double E1x,E1y,E1z;
    long double r1,r2,r3,r4;
    int S123,S124,S234,S134;
    int A1A2E123,A1A3E123,A2A3E123;
    int A1A2E124,A1A4E124,A2A4E124;
    int A1A3E134,A1A4E134,A3A4E134;
    int A2A3E234,A2A4E234,A3A4E234;
    //long double Vx,Vy,Vz;
    long double Vx=v->p->x;
    long double Vy=v->p->y;
    long double Vz=v->p->z;


    A1=Atoms[v->D->AT[0]].p;
    A2=Atoms[v->D->AT[1]].p;
    A3=Atoms[v->D->AT[2]].p;
    A4=Atoms[v->D->AT[3]].p;

    r1=Atoms[v->D->AT[0]].radius+r_cut;
    r2=Atoms[v->D->AT[1]].radius+r_cut;
    r3=Atoms[v->D->AT[2]].radius+r_cut;
    r4=Atoms[v->D->AT[3]].radius+r_cut;

    site E123,E124,E134,E234;
    site B12,B13,B14,B23,B34,B24;
    long double MIDP234x,MIDP234y,MIDP234z;
    site a1,a2,a3,a4;
    long double X,Y,Z;

    S123=sign_aaav(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,Vx,Vy,Vz);
    S124=sign_aaav(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,A3.x,A3.y,A3.z,Vx,Vy,Vz);
    S234=sign_aaav(A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A1.x,A1.y,A1.z,Vx,Vy,Vz);
    S134=sign_aaav(A1.x,A1.y,A1.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A2.x,A2.y,A2.z,Vx,Vy,Vz);

    E123=v->D->MID[0][1][2];
    E124=v->D->MID[0][1][3];
    E134=v->D->MID[0][2][3];
    E234=v->D->MID[1][2][3];

    B12=v->D->MIDP[0][1];
    B13=v->D->MIDP[0][2];
    B14=v->D->MIDP[0][3];
    B23=v->D->MIDP[1][2];
    B24=v->D->MIDP[1][3];
    B34=v->D->MIDP[2][3];

    A1A2E123=sign_aaae(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,E123.x,E123.y,E123.z);
    A1A3E123=sign_aaae(A1.x,A1.y,A1.z,A3.x,A3.y,A3.z,A2.x,A2.y,A2.z,E123.x,E123.y,E123.z);
    A2A3E123=sign_aaae(A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A1.x,A1.y,A1.z,E123.x,E123.y,E123.z);

    A1A2E124=sign_aaae(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,E124.x,E124.y,E124.z);
    A1A4E124=sign_aaae(A1.x,A1.y,A1.z,A4.x,A4.y,A4.z,A2.x,A2.y,A2.z,E124.x,E124.y,E124.z);
    A2A4E124=sign_aaae(A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,A1.x,A1.y,A1.z,E124.x,E124.y,E124.z);

    A1A3E134=sign_aaae(A1.x,A1.y,A1.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,E134.x,E134.y,E134.z);
    A1A4E134=sign_aaae(A1.x,A1.y,A1.z,A4.x,A4.y,A4.z,A3.x,A3.y,A3.z,E134.x,E134.y,E134.z);
    A3A4E134=sign_aaae(A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A1.x,A1.y,A1.z,E134.x,E134.y,E134.z);

    A2A3E234=sign_aaae(A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,E234.x,E234.y,E234.z);
    A2A4E234=sign_aaae(A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,A3.x,A3.y,A3.z,E234.x,E234.y,E234.z);
    A3A4E234=sign_aaae(A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A2.x,A2.y,A2.z,E234.x,E234.y,E234.z);

    vol=vol+S123*A1A2E123*volume_tetrahedron(A1.x,A1.y,A1.z,E123.x,E123.y,E123.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r1);
    vol=vol+S123*A1A3E123*volume_tetrahedron(A1.x,A1.y,A1.z,E123.x,E123.y,E123.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r1);
    vol=vol+S123*A1A2E123*volume_tetrahedron(A2.x,A2.y,A2.z,E123.x,E123.y,E123.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r2);
    vol=vol+S123*A2A3E123*volume_tetrahedron(A2.x,A2.y,A2.z,E123.x,E123.y,E123.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r2);
    vol=vol+S123*A1A3E123*volume_tetrahedron(A3.x,A3.y,A3.z,E123.x,E123.y,E123.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r3);
    vol=vol+S123*A2A3E123*volume_tetrahedron(A3.x,A3.y,A3.z,E123.x,E123.y,E123.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r3);

    vol=vol+S124*A1A2E124*volume_tetrahedron(A1.x,A1.y,A1.z,E124.x,E124.y,E124.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r1);
    vol=vol+S124*A1A4E124*volume_tetrahedron(A1.x,A1.y,A1.z,E124.x,E124.y,E124.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r1);
    vol=vol+S124*A1A2E124*volume_tetrahedron(A2.x,A2.y,A2.z,E124.x,E124.y,E124.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r2);
    vol=vol+S124*A2A4E124*volume_tetrahedron(A2.x,A2.y,A2.z,E124.x,E124.y,E124.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r2);
    vol=vol+S124*A1A4E124*volume_tetrahedron(A4.x,A4.y,A4.z,E124.x,E124.y,E124.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r4);
    vol=vol+S124*A2A4E124*volume_tetrahedron(A4.x,A4.y,A4.z,E124.x,E124.y,E124.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r4);

    vol=vol+S134*A1A3E134*volume_tetrahedron(A1.x,A1.y,A1.z,E134.x,E134.y,E134.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r1);
    vol=vol+S134*A1A4E134*volume_tetrahedron(A1.x,A1.y,A1.z,E134.x,E134.y,E134.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r1);
    vol=vol+S134*A1A3E134*volume_tetrahedron(A3.x,A3.y,A3.z,E134.x,E134.y,E134.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r3);
    vol=vol+S134*A3A4E134*volume_tetrahedron(A3.x,A3.y,A3.z,E134.x,E134.y,E134.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r3);
    vol=vol+S134*A1A4E134*volume_tetrahedron(A4.x,A4.y,A4.z,E134.x,E134.y,E134.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r4);
    vol=vol+S134*A3A4E134*volume_tetrahedron(A4.x,A4.y,A4.z,E134.x,E134.y,E134.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r4);

    vol=vol+S234*A2A3E234*volume_tetrahedron(A2.x,A2.y,A2.z,E234.x,E234.y,E234.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r2);
    vol=vol+S234*A2A4E234*volume_tetrahedron(A2.x,A2.y,A2.z,E234.x,E234.y,E234.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r2);
    vol=vol+S234*A2A3E234*volume_tetrahedron(A3.x,A3.y,A3.z,E234.x,E234.y,E234.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r3);
    vol=vol+S234*A3A4E234*volume_tetrahedron(A3.x,A3.y,A3.z,E234.x,E234.y,E234.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r3);
    vol=vol+S234*A2A4E234*volume_tetrahedron(A4.x,A4.y,A4.z,E234.x,E234.y,E234.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r4);
    vol=vol+S234*A3A4E234*volume_tetrahedron(A4.x,A4.y,A4.z,E234.x,E234.y,E234.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r4);

    return vol;
}
void display_cell(int SAM,atom Atoms[])
{
    for(int i=0; i<Atoms[SAM].conti; i++)
    {
        cout<<SAM<<"\t"<<Atoms[SAM].contigous[i]<<"\n";
        atom *ATOM;
        site *p=nullptr;
        ATOM=&(Atoms[SAM]);
        container <edge> *temp_ec=nullptr;
        temp_ec=ATOM->F[i].initial;
        int flag=1;
        while(temp_ec)
        {
            //print_edge(temp_ec->t);
            if(!temp_ec->t->evoid)
            {
                flag=0;
            }
            temp_ec=temp_ec->next;
        }
        if(flag)
        {
            container <edge> *temp_ec=nullptr;
            temp_ec=ATOM->F[i].initial;
            ////display_atom(&(Atoms[SAM]));
            ////display_atom(&(Atoms[Atoms[SAM].contigous[i]]));
            long double dis=sqrtl(distancesq(Atoms[SAM].p,Atoms[Atoms[SAM].contigous[i]].p));
            ////if(dis<(Atoms[SAM].radius+Atoms[Atoms[SAM].contigous[i]].radius+2.*r_cut))
            ////while(temp_ec)
            ////{
            ////	print_edge(temp_ec->t);
            ////    cout<<temp_ec->t->A1<<"\t"<<temp_ec->t->A2<<"\t"<<temp_ec->t->A3<<"\n";
            ////    display_atom(&(Atoms[temp_ec->t->A1]));
            ////    display_atom(&(Atoms[temp_ec->t->A2]));
            ////    display_atom(&(Atoms[temp_ec->t->A3]));
            ////    temp_ec=temp_ec->next;
            ////}
        }
        ////if(!temp_ec)
        ////{
        ////
        ////}
        ////cout<<"##compare\n";
        ////for(int j=0;j<Atoms[SAM].conti && !Atoms[ATOM->contigous[i]].bor;j++)
        ////{
        ////	if(Atoms[SAM].part_c[i][j]==2)
        ////	{
        ////		delunay *D_ONE=nullptr,*D_TWO=nullptr;
        ////		for(int k=0; k<ATOM->conti; k++)
        ////		{
        ////			if(D_ONE && D_TWO)
        ////			{
        ////				break;
        ////			}
        ////			if(k!=i && k!=j)
        ////			{
        ////				container_delunay *temp;
        ////				temp=ATOM->D_FIRST;
        ////				//int flag=1;
        ////				//cout<<i<<"\t"<<j<<"\t"<<k<<"\n";
        ////				while(1)
        ////				{
        ////					std::vector<int> EV1 {ATOM->index,ATOM->contigous[i],ATOM->contigous[j],ATOM->contigous[k]};//,A3,A4;
        ////					std::sort(EV1.begin(),EV1.end());
        ////					if(temp->D->AT[0]==EV1[0] && temp->D->AT[1]==EV1[1] && temp->D->AT[2]==EV1[2] && temp->D->AT[3]==EV1[3])
        ////					{
        ////						if(!D_ONE)
        ////						{
        ////							D_ONE=temp->D;
        ////							break;
        ////						}
        ////						else
        ////						{
        ////							D_TWO=temp->D;
        ////							break;
        ////						}

        ////					}
        ////					if(temp->next)
        ////						temp=temp->next;
        ////					else
        ////						break;
        ////				}
        ////			}
        ////		}
        ////		if(!p)
        ////		{
        ////			p=new site;
        ////			p->x=D_ONE->circum_x;
        ////			p->y=D_ONE->circum_y;
        ////			p->z=D_ONE->circum_z;
        ////		}
        ////		if(D_ONE && D_TWO)
        ////		{
        ////			long double dis=sqrtl(distancesq(ATOM->RMID[i],Atoms[ATOM->contigous[j]].p));
        ////			if(dis<Atoms[ATOM->contigous[j]].radius+r_cut)
        ////			{
        ////			////for(int p=0;p<D_ONE->v->v_neigh_count;p++)
        ////			////{
        ////			////	if(D_ONE->v->neib_vert[p]==D_TWO->v)
        ////			////	{
        ////			////		//cout<<D_ONE->v->neib_ed[p]<<" prev\n";
        ////			////		D_ONE->v->neib_ed[p]=0;
        ////			////	}
        ////			////}
        ////			////for(int p=0;p<D_TWO->v->v_neigh_count;p++)
        ////			////{
        ////			////	if(D_TWO->v->neib_vert[p]==D_ONE->v)
        ////			////	{
        ////			////		//cout<<D_TWO->v->neib_ed[p]<<" prev\n";
        ////			////		D_TWO->v->neib_ed[p]=0;
        ////			////	}
        ////			////}
        ////			////cout<<"mol new\n";
        ////			////cout<<"draw material Opaque\n";
        ////			////cout<<"draw color white\n";
        ////			////cout<<"draw line\t{\t";
        ////			////cout<<ATOM->p.x<<"\t"<<ATOM->p.y<<"\t"<<ATOM->p.z<<"}\t{";
        ////			////cout<<Atoms[ATOM->contigous[i]].p.x<<"\t"<<Atoms[ATOM->contigous[i]].p.y<<"\t"<<Atoms[ATOM->contigous[i]].p.z<<"}\twidth 2\n";
        ////			}
        ////		    cout<<"draw line \t{";
        ////		    cout<<D_ONE->circum_x<<"\t"<<D_ONE->circum_y<<"\t"<<D_ONE->circum_z<<"}\t{";
        ////		    cout<<D_TWO->circum_x<<"\t"<<D_TWO->circum_y<<"\t"<<D_TWO->circum_z<<"}\n";
        ////		    //cout<<p->x<<"\t"<<p->y<<"\t"<<p->z<<"}\n";
        ////		}
        ////	}
        ////}
    }
}
int main( int argc, char * argv[] )
{
    atom *Atoms=nullptr;
    int counter=0;
    int config_count=0;
    //The file with configurations
    std::ifstream infile(argv[1]);///dat_trial");//config_2000_0.38_2_0.70.dat");
    //No of Atoms
    infile>>nAtoms;
    //No of configurations in the input file
    config_count=1;
    //No of types of particle
    int ntypes=5;
    int SAM=0;

    long double b,c,d,e,f;
    //radiuses of the particle
    //nAtoms=0;
    Atoms = new (nothrow) atom[nAtoms];
    initialization(Atoms,nAtoms,ntypes);

// radius[0]=0.5;
// radius[1]=0.7;

    radius[0]=1.1;
    radius[1]=1.52;
    radius[2]=1.55;
    radius[3]=1.7;
    radius[4]=1.8;

    if(3*2*radius[1]<box)
    {
        LAYER_CUT=3*2*radius[1]-box;
    }
    else
        LAYER_CUT=0.;

    char buffer[64];
    vertice *temp_site=nullptr;
    //long long free_dist[10000]= {0};
    long double dummy;
    //ofstream fdist;
    //fdist.open(buffer);
    //This loop is over all the configurations
    ofstream fdist;
    int resno;
    for(int nconfig=0; nconfig<config_count; nconfig++)
    {
        counter=0;
        infile>>box;
        infile>>resno;
        twob=2.0*box;
        while(infile>>b>>c>>d>>e>>f)
        {
            Atoms[counter].p.x=b;
            Atoms[counter].p.y=c;
            Atoms[counter].p.z=d;
            Atoms[counter].radius=e-epsilon;
            Atoms[counter].index=counter;
            Atoms[counter].neighbours=0;
            counter++;
            if(counter==nAtoms)
            {
                break;
            }
        }
        update_neighbours(Atoms,nAtoms);
        for(int i=0; i<nAtoms; i++)
        {
            int TYPE=0;
            for(int t=0; t<ntypes; t++)
            {
                if(abs(Atoms[i].radius-radius[t]-epsilon)<0.00001)
                {
                    Atoms[i].type=t;
                    break;
                }
            }

            Atoms[i].Cstart=nullptr;
            Atoms[i].D_FIRST=nullptr;
            Atoms[i].conti=0;

            for(int j=0; j<nAtoms; j++)
            {
                Atoms[i].conti_index[j]=-1;
                contigous[i][j]=0;
            }
        }
        long double area=0;
        vertice *save=nullptr;
        /* here all the arrays are initialized*/

        {
            int TYPE=0;
            snprintf(buffer,sizeof(char)*64,"vor_%d",int(TYPE));//_%d_%f.dat",int(nAtoms),Press);
            ofstream vor;
            vor.open(buffer);
            //r_cut=radius;
            //r_cut=0.0;
            r_cut=1.4;
            int max_conti=0;
            int max_conti_index=-1;
            for(SAM=0; SAM< nAtoms; SAM++)
            {
                {
                    //cout<<"mol new\n";
                    //cout<<"draw material Transparent\n";
                    // if(SAM==146 || SAM==150 || SAM==128)
                    // {
                    //   cout<<"draw color blue\n";
                    //cout<<"##\t"<<SAM<<"\n";
                    //   cout<<"draw sphere\t{";
                    //   cout<<Atoms[SAM].p.x<<"\t"<<Atoms[SAM].p.y<<"\t"<<Atoms[SAM].p.z<<"}\tradius\t"<<Atoms[SAM].radius+r_cut<<"\tresolution 15\n";
                    // }
                    ////cout<<Atoms[SAM].p.x<<"\t"<<Atoms[SAM].p.y<<"\t"<<Atoms[SAM].p.z<<"}\tradius\t"<<0.3<<"\tresolution 25\n";
                    ////if(!Atoms[SAM].D_FIRST)
                    ////{
                    ////   cout<<"whewq\n";
                    ////}

                    cout<<std::flush;
                    if(!Atoms[SAM].D_FIRST)
                    {
                        first_delunay(&(Atoms[SAM]),Atoms,TYPE);
                    }
                    if(Atoms[SAM].D_FIRST)
                        complete_del_2(&(Atoms[SAM]),Atoms,nAtoms,TYPE);
                    if(max_conti<Atoms[SAM].conti)
                    {
                        max_conti=Atoms[SAM].conti;
                        max_conti_index=SAM;
                    }
                    //cout<<Atoms[SAM].conti<<"\n";
                    //max_conti=max_conti+Atoms[SAM].conti;
                    container_delunay *temp_d=nullptr;
                    temp_d=Atoms[SAM].D_FIRST;
                    Atoms[SAM].S = new set_of_cvert;
                    while(temp_d)
                    {
                        //if(SAM==666)
                        //	print_delunay(temp_d->D,Atoms,nAtoms);
                        container_vertice *temp=nullptr;
                        temp = new container_vertice;
                        temp->V=temp_d->D->v;
                        Atoms[SAM].S->insert_cvert(temp);
                        temp_d=temp_d->next;
                    }

                }
            }
            //cout<<"draw color red\n";
            vertice *temp;
            temp=V->initial;
            void_vert_count=0;
            int vert_count=0;
            int color;
            while(1)
            {
                if(temp->is_void==1)
                {
                    temp->cluster_index=void_vert_count;
                    ////cout<<"draw sphere \t{";
                    ////cout<<temp->p->x<<"\t"<<temp->p->y<<"\t"<<temp->p->z<<"}\t radius 0.05\t resolution 100\n";
                    void_vert_count=void_vert_count+1;
                }
                vert_count++;

                if(temp->next)
                {
                    temp=temp->next;
                }
                else
                    break;
            }
            //cout<<void_vert_count<<"\n";
            vertice **cavity_list=nullptr;
            cavity_list = new (nothrow) vertice*[void_vert_count];
            vertice **vertex_list=nullptr;
            vertex_list = new (nothrow) vertice*[vert_count];
            int *pocket;
            pocket = new (nothrow) int[void_vert_count];
            ////cout<<"here\n";
            vertice *temp_start;
            temp_start=V->initial;
            {
                int i=0;
                int j=0;
                while(1)
                {
                    if(temp_start->is_void)
                    {
                        cavity_list[i]=temp_start;
                        pocket[i]=0;
                        i++;
                    }
                    vertex_list[j]=temp_start;
                    j++;
                    if(temp_start->next)
                        temp_start=temp_start->next;
                    else
                        break;
                }
            }
            //cout<<max_conti<<"\n";
            double *cav_vol;
            cav_vol= new (nothrow) double[void_vert_count];
            double *cav_area;
            cav_area= new (nothrow) double[void_vert_count];
            //cout<<"draw material transparent\n";
            //cout<<"mol new\n";
            //cout<<"draw material Opaque\n";
            //for(int i=0; i<void_vert_count; i++)
            {
                //cav_vol[i]=0;
                //cav_area[i]=0;
                //cout<<"##\t"<<i<<"\n";
                for(int j=0; j<vert_count; j++)
                {
                    //if(cavity_list[j]->cluster_index==i)
                    {
                        ////    if(i==111)
                        ////    {
                        ////    	cout<<"draw color "<<i%16<<"\n";;
                        ////    	print_delunay_solid(cavity_list[j]->D,Atoms,nAtoms);
                        ////    ////for(int p=0;p<cavity_list[j]->v_neigh_count;p++)
                        ////    ////{
                        ////    ////	cout<<cavity_list[j]->neib_ed[p]<<"\n";
                        ////    ////	print_delunay_solid(cavity_list[j]->neib_vert[p]->D,Atoms,nAtoms);
                        ////    ////}
                        ////    }
                        site A1,A2,A3,A4;
                        long double E1x,E1y,E1z;
                        long double r1,r2,r3,r4;
                        int S123,S124,S234,S134;
                        int A1A2E123,A1A3E123,A2A3E123;
                        int A1A2E124,A1A4E124,A2A4E124;
                        int A1A3E134,A1A4E134,A3A4E134;
                        int A2A3E234,A2A4E234,A3A4E234;
                        //long double Vx,Vy,Vz;
                        long double Vx=vertex_list[j]->p->x;
                        long double Vy=vertex_list[j]->p->y;
                        long double Vz=vertex_list[j]->p->z;

                        A1=Atoms[vertex_list[j]->D->AT[0]].p;
                        A2=Atoms[vertex_list[j]->D->AT[1]].p;
                        A3=Atoms[vertex_list[j]->D->AT[2]].p;
                        A4=Atoms[vertex_list[j]->D->AT[3]].p;

                        r1=Atoms[vertex_list[j]->D->AT[0]].radius+r_cut;
                        r2=Atoms[vertex_list[j]->D->AT[1]].radius+r_cut;
                        r3=Atoms[vertex_list[j]->D->AT[2]].radius+r_cut;
                        r4=Atoms[vertex_list[j]->D->AT[3]].radius+r_cut;

                        site E123,E124,E134,E234;
                        site B12,B13,B14,B23,B34,B24;
                        long double MIDP234x,MIDP234y,MIDP234z;
                        site a1,a2,a3,a4;
                        long double X,Y,Z;

                        S123=sign_aaav(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,Vx,Vy,Vz);
                        S124=sign_aaav(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,A3.x,A3.y,A3.z,Vx,Vy,Vz);
                        S234=sign_aaav(A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A1.x,A1.y,A1.z,Vx,Vy,Vz);
                        S134=sign_aaav(A1.x,A1.y,A1.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A2.x,A2.y,A2.z,Vx,Vy,Vz);

                        E123=vertex_list[j]->D->MID[0][1][2];
                        E124=vertex_list[j]->D->MID[0][1][3];
                        E134=vertex_list[j]->D->MID[0][2][3];
                        E234=vertex_list[j]->D->MID[1][2][3];

                        B12=vertex_list[j]->D->MIDP[0][1];
                        B13=vertex_list[j]->D->MIDP[0][2];
                        B14=vertex_list[j]->D->MIDP[0][3];
                        B23=vertex_list[j]->D->MIDP[1][2];
                        B24=vertex_list[j]->D->MIDP[1][3];
                        B34=vertex_list[j]->D->MIDP[2][3];


                        A1A2E123=sign_aaae(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,E123.x,E123.y,E123.z);
                        A1A3E123=sign_aaae(A1.x,A1.y,A1.z,A3.x,A3.y,A3.z,A2.x,A2.y,A2.z,E123.x,E123.y,E123.z);
                        A2A3E123=sign_aaae(A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A1.x,A1.y,A1.z,E123.x,E123.y,E123.z);

                        A1A2E124=sign_aaae(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,E124.x,E124.y,E124.z);
                        A1A4E124=sign_aaae(A1.x,A1.y,A1.z,A4.x,A4.y,A4.z,A2.x,A2.y,A2.z,E124.x,E124.y,E124.z);
                        A2A4E124=sign_aaae(A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,A1.x,A1.y,A1.z,E124.x,E124.y,E124.z);

                        A1A3E134=sign_aaae(A1.x,A1.y,A1.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,E134.x,E134.y,E134.z);
                        A1A4E134=sign_aaae(A1.x,A1.y,A1.z,A4.x,A4.y,A4.z,A3.x,A3.y,A3.z,E134.x,E134.y,E134.z);
                        A3A4E134=sign_aaae(A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A1.x,A1.y,A1.z,E134.x,E134.y,E134.z);

                        A2A3E234=sign_aaae(A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,E234.x,E234.y,E234.z);
                        A2A4E234=sign_aaae(A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,A3.x,A3.y,A3.z,E234.x,E234.y,E234.z);
                        A3A4E234=sign_aaae(A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A2.x,A2.y,A2.z,E234.x,E234.y,E234.z);

                        //   cav_vol[i]=cav_vol[i]+S123*A1A2E123*volume_tetrahedron(A1.x,A1.y,A1.z,E123.x,E123.y,E123.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r1);
                        //   cav_vol[i]=cav_vol[i]+S123*A1A3E123*volume_tetrahedron(A1.x,A1.y,A1.z,E123.x,E123.y,E123.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r1);
                        //   cav_vol[i]=cav_vol[i]+S123*A1A2E123*volume_tetrahedron(A2.x,A2.y,A2.z,E123.x,E123.y,E123.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r2);
                        //   cav_vol[i]=cav_vol[i]+S123*A2A3E123*volume_tetrahedron(A2.x,A2.y,A2.z,E123.x,E123.y,E123.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r2);
                        //   cav_vol[i]=cav_vol[i]+S123*A1A3E123*volume_tetrahedron(A3.x,A3.y,A3.z,E123.x,E123.y,E123.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r3);
                        //   cav_vol[i]=cav_vol[i]+S123*A2A3E123*volume_tetrahedron(A3.x,A3.y,A3.z,E123.x,E123.y,E123.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r3);

                        //   cav_vol[i]=cav_vol[i]+S124*A1A2E124*volume_tetrahedron(A1.x,A1.y,A1.z,E124.x,E124.y,E124.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r1);
                        //   cav_vol[i]=cav_vol[i]+S124*A1A4E124*volume_tetrahedron(A1.x,A1.y,A1.z,E124.x,E124.y,E124.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r1);
                        //   cav_vol[i]=cav_vol[i]+S124*A1A2E124*volume_tetrahedron(A2.x,A2.y,A2.z,E124.x,E124.y,E124.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r2);
                        //   cav_vol[i]=cav_vol[i]+S124*A2A4E124*volume_tetrahedron(A2.x,A2.y,A2.z,E124.x,E124.y,E124.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r2);
                        //   cav_vol[i]=cav_vol[i]+S124*A1A4E124*volume_tetrahedron(A4.x,A4.y,A4.z,E124.x,E124.y,E124.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r4);
                        //   cav_vol[i]=cav_vol[i]+S124*A2A4E124*volume_tetrahedron(A4.x,A4.y,A4.z,E124.x,E124.y,E124.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r4);

                        //   cav_vol[i]=cav_vol[i]+S134*A1A3E134*volume_tetrahedron(A1.x,A1.y,A1.z,E134.x,E134.y,E134.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r1);
                        //   cav_vol[i]=cav_vol[i]+S134*A1A4E134*volume_tetrahedron(A1.x,A1.y,A1.z,E134.x,E134.y,E134.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r1);
                        //   cav_vol[i]=cav_vol[i]+S134*A1A3E134*volume_tetrahedron(A3.x,A3.y,A3.z,E134.x,E134.y,E134.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r3);
                        //   cav_vol[i]=cav_vol[i]+S134*A3A4E134*volume_tetrahedron(A3.x,A3.y,A3.z,E134.x,E134.y,E134.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r3);
                        //   cav_vol[i]=cav_vol[i]+S134*A1A4E134*volume_tetrahedron(A4.x,A4.y,A4.z,E134.x,E134.y,E134.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r4);
                        //   cav_vol[i]=cav_vol[i]+S134*A3A4E134*volume_tetrahedron(A4.x,A4.y,A4.z,E134.x,E134.y,E134.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r4);

                        //   cav_vol[i]=cav_vol[i]+S234*A2A3E234*volume_tetrahedron(A2.x,A2.y,A2.z,E234.x,E234.y,E234.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r2);
                        //   cav_vol[i]=cav_vol[i]+S234*A2A4E234*volume_tetrahedron(A2.x,A2.y,A2.z,E234.x,E234.y,E234.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r2);
                        //   cav_vol[i]=cav_vol[i]+S234*A2A3E234*volume_tetrahedron(A3.x,A3.y,A3.z,E234.x,E234.y,E234.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r3);
                        //   cav_vol[i]=cav_vol[i]+S234*A3A4E234*volume_tetrahedron(A3.x,A3.y,A3.z,E234.x,E234.y,E234.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r3);
                        //   cav_vol[i]=cav_vol[i]+S234*A2A4E234*volume_tetrahedron(A4.x,A4.y,A4.z,E234.x,E234.y,E234.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r4);
                        //   cav_vol[i]=cav_vol[i]+S234*A3A4E234*volume_tetrahedron(A4.x,A4.y,A4.z,E234.x,E234.y,E234.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r4);

                        Atoms[vertex_list[j]->D->AT[0]].vor_vol=Atoms[vertex_list[j]->D->AT[0]].vor_vol+S123*A1A2E123*volume_tetrahedron(A1.x,A1.y,A1.z,E123.x,E123.y,E123.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r1,1);
                        Atoms[vertex_list[j]->D->AT[0]].vor_vol=Atoms[vertex_list[j]->D->AT[0]].vor_vol+S123*A1A3E123*volume_tetrahedron(A1.x,A1.y,A1.z,E123.x,E123.y,E123.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r1,1);
                        Atoms[vertex_list[j]->D->AT[0]].vor_vol=Atoms[vertex_list[j]->D->AT[0]].vor_vol+S124*A1A2E124*volume_tetrahedron(A1.x,A1.y,A1.z,E124.x,E124.y,E124.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r1,1);
                        Atoms[vertex_list[j]->D->AT[0]].vor_vol=Atoms[vertex_list[j]->D->AT[0]].vor_vol+S124*A1A4E124*volume_tetrahedron(A1.x,A1.y,A1.z,E124.x,E124.y,E124.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r1,1);
                        Atoms[vertex_list[j]->D->AT[0]].vor_vol=Atoms[vertex_list[j]->D->AT[0]].vor_vol+S134*A1A3E134*volume_tetrahedron(A1.x,A1.y,A1.z,E134.x,E134.y,E134.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r1,1);
                        Atoms[vertex_list[j]->D->AT[0]].vor_vol=Atoms[vertex_list[j]->D->AT[0]].vor_vol+S134*A1A4E134*volume_tetrahedron(A1.x,A1.y,A1.z,E134.x,E134.y,E134.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r1,1);

                        Atoms[vertex_list[j]->D->AT[1]].vor_vol=Atoms[vertex_list[j]->D->AT[1]].vor_vol+S123*A1A2E123*volume_tetrahedron(A2.x,A2.y,A2.z,E123.x,E123.y,E123.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r2,1);
                        Atoms[vertex_list[j]->D->AT[1]].vor_vol=Atoms[vertex_list[j]->D->AT[1]].vor_vol+S123*A2A3E123*volume_tetrahedron(A2.x,A2.y,A2.z,E123.x,E123.y,E123.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r2,1);
                        Atoms[vertex_list[j]->D->AT[1]].vor_vol=Atoms[vertex_list[j]->D->AT[1]].vor_vol+S124*A1A2E124*volume_tetrahedron(A2.x,A2.y,A2.z,E124.x,E124.y,E124.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r2,1);
                        Atoms[vertex_list[j]->D->AT[1]].vor_vol=Atoms[vertex_list[j]->D->AT[1]].vor_vol+S124*A2A4E124*volume_tetrahedron(A2.x,A2.y,A2.z,E124.x,E124.y,E124.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r2,1);
                        Atoms[vertex_list[j]->D->AT[1]].vor_vol=Atoms[vertex_list[j]->D->AT[1]].vor_vol+S234*A2A3E234*volume_tetrahedron(A2.x,A2.y,A2.z,E234.x,E234.y,E234.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r2,1);
                        Atoms[vertex_list[j]->D->AT[1]].vor_vol=Atoms[vertex_list[j]->D->AT[1]].vor_vol+S234*A2A4E234*volume_tetrahedron(A2.x,A2.y,A2.z,E234.x,E234.y,E234.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r2,1);

                        Atoms[vertex_list[j]->D->AT[2]].vor_vol=Atoms[vertex_list[j]->D->AT[2]].vor_vol+S123*A1A3E123*volume_tetrahedron(A3.x,A3.y,A3.z,E123.x,E123.y,E123.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r3,1);
                        Atoms[vertex_list[j]->D->AT[2]].vor_vol=Atoms[vertex_list[j]->D->AT[2]].vor_vol+S123*A2A3E123*volume_tetrahedron(A3.x,A3.y,A3.z,E123.x,E123.y,E123.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r3,1);
                        Atoms[vertex_list[j]->D->AT[2]].vor_vol=Atoms[vertex_list[j]->D->AT[2]].vor_vol+S134*A1A3E134*volume_tetrahedron(A3.x,A3.y,A3.z,E134.x,E134.y,E134.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r3,1);
                        Atoms[vertex_list[j]->D->AT[2]].vor_vol=Atoms[vertex_list[j]->D->AT[2]].vor_vol+S134*A3A4E134*volume_tetrahedron(A3.x,A3.y,A3.z,E134.x,E134.y,E134.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r3,1);
                        Atoms[vertex_list[j]->D->AT[2]].vor_vol=Atoms[vertex_list[j]->D->AT[2]].vor_vol+S234*A2A3E234*volume_tetrahedron(A3.x,A3.y,A3.z,E234.x,E234.y,E234.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r3,1);
                        Atoms[vertex_list[j]->D->AT[2]].vor_vol=Atoms[vertex_list[j]->D->AT[2]].vor_vol+S234*A3A4E234*volume_tetrahedron(A3.x,A3.y,A3.z,E234.x,E234.y,E234.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r3,1);

                        Atoms[vertex_list[j]->D->AT[3]].vor_vol=Atoms[vertex_list[j]->D->AT[3]].vor_vol+S124*A1A4E124*volume_tetrahedron(A4.x,A4.y,A4.z,E124.x,E124.y,E124.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r4,1);
                        Atoms[vertex_list[j]->D->AT[3]].vor_vol=Atoms[vertex_list[j]->D->AT[3]].vor_vol+S124*A2A4E124*volume_tetrahedron(A4.x,A4.y,A4.z,E124.x,E124.y,E124.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r4,1);
                        Atoms[vertex_list[j]->D->AT[3]].vor_vol=Atoms[vertex_list[j]->D->AT[3]].vor_vol+S134*A1A4E134*volume_tetrahedron(A4.x,A4.y,A4.z,E134.x,E134.y,E134.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r4,1);
                        Atoms[vertex_list[j]->D->AT[3]].vor_vol=Atoms[vertex_list[j]->D->AT[3]].vor_vol+S134*A3A4E134*volume_tetrahedron(A4.x,A4.y,A4.z,E134.x,E134.y,E134.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r4,1);
                        Atoms[vertex_list[j]->D->AT[3]].vor_vol=Atoms[vertex_list[j]->D->AT[3]].vor_vol+S234*A2A4E234*volume_tetrahedron(A4.x,A4.y,A4.z,E234.x,E234.y,E234.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r4,1);
                        Atoms[vertex_list[j]->D->AT[3]].vor_vol=Atoms[vertex_list[j]->D->AT[3]].vor_vol+S234*A3A4E234*volume_tetrahedron(A4.x,A4.y,A4.z,E234.x,E234.y,E234.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r4,1);



                    }
                }
            }
            long double uni_sphere=0.;
            long double sum_sphere=0.;
            int borat=0;
            for(int i=0; i<nAtoms; i++)
            {
                if(Atoms[i].radius!=1.0)
                {
                    uni_sphere=uni_sphere+Atoms[i].vor_vol;
                    sum_sphere=sum_sphere+4./3.*M_PI*powl(Atoms[i].radius+r_cut,3);
                    //cout<<Atoms[i].p.x<<"\t"<<Atoms[i].p.y<<"\t"<<Atoms[i].p.z<<"\t"<<Atoms[i].radius+r_cut<<"\n";
                    cout<<i<<"\t"<<4./3.*M_PI*powl(Atoms[i].radius+r_cut,3)<<"\t"<<Atoms[i].vor_vol<<"\n";
                }
                else
                {
                    borat=borat+1;
                }
                //cout<<4./3.*M_PI*powl(Atoms[i].radius,3)<<"\t"<<Atoms[i].vor_vol<<"\n";
            }
            cout<<std::setprecision(16);
            cout<<resno<<"\t"<<uni_sphere<<"\t"<<convex_vol<<"\t"<<sum_sphere<<"\t"<<uni_sphere/convex_vol<<"\t"<<borat<<"\t"<<borat*1./nAtoms<<"\n";
            for(int j=0; j<void_vert_count; j++)
            {
                site A1,A2,A3,A4;
                long double E1x,E1y,E1z;
                long double r1,r2,r3,r4;
                int S123,S124,S234,S134;
                int A1A2E123,A1A3E123,A2A3E123;
                int A1A2E124,A1A4E124,A2A4E124;
                int A1A3E134,A1A4E134,A3A4E134;
                int A2A3E234,A2A4E234,A3A4E234;
                //long double Vx,Vy,Vz;
                long double Vx=cavity_list[j]->p->x;
                long double Vy=cavity_list[j]->p->y;
                long double Vz=cavity_list[j]->p->z;

                A1=Atoms[cavity_list[j]->D->AT[0]].p;
                A2=Atoms[cavity_list[j]->D->AT[1]].p;
                A3=Atoms[cavity_list[j]->D->AT[2]].p;
                A4=Atoms[cavity_list[j]->D->AT[3]].p;

                //print_delunay(cavity_list[j]->D,Atoms,nAtoms);

                r1=Atoms[cavity_list[j]->D->AT[0]].radius+r_cut;
                r2=Atoms[cavity_list[j]->D->AT[1]].radius+r_cut;
                r3=Atoms[cavity_list[j]->D->AT[2]].radius+r_cut;
                r4=Atoms[cavity_list[j]->D->AT[3]].radius+r_cut;

                site E123,E124,E134,E234;
                site B12,B13,B14,B23,B34,B24;
                long double MIDP234x,MIDP234y,MIDP234z;
                site a1,a2,a3,a4;
                long double X,Y,Z;

                S123=sign_aaav(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,Vx,Vy,Vz);
                S124=sign_aaav(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,A3.x,A3.y,A3.z,Vx,Vy,Vz);
                S234=sign_aaav(A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A1.x,A1.y,A1.z,Vx,Vy,Vz);
                S134=sign_aaav(A1.x,A1.y,A1.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A2.x,A2.y,A2.z,Vx,Vy,Vz);

                E123=cavity_list[j]->D->MID[0][1][2];
                E124=cavity_list[j]->D->MID[0][1][3];
                E134=cavity_list[j]->D->MID[0][2][3];
                E234=cavity_list[j]->D->MID[1][2][3];

            //  display_SITE(&E123);
            //  display_SITE(&E124);
            //  display_SITE(&E134);
            //  display_SITE(&E234);
			////cout<<"draw color red\n";
            ////for(int p=0; p<cavity_list[j]->v_neigh_count; p++)
            ////{
            ////	cout<<"draw line\t{";
            ////	cout<<cavity_list[j]->p->x<<"\t"<<cavity_list[j]->p->y<<"\t"<<cavity_list[j]->p->z<<"}\t{";
            ////	cout<<cavity_list[j]->neib_vert[p]->p->x<<"\t"<<cavity_list[j]->neib_vert[p]->p->y<<"\t"<<cavity_list[j]->neib_vert[p]->p->z<<"}\n";
            ////}
                B12=cavity_list[j]->D->MIDP[0][1];
                B13=cavity_list[j]->D->MIDP[0][2];
                B14=cavity_list[j]->D->MIDP[0][3];
                B23=cavity_list[j]->D->MIDP[1][2];
                B24=cavity_list[j]->D->MIDP[1][3];
                B34=cavity_list[j]->D->MIDP[2][3];


                A1A2E123=sign_aaae(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,E123.x,E123.y,E123.z);
                A1A3E123=sign_aaae(A1.x,A1.y,A1.z,A3.x,A3.y,A3.z,A2.x,A2.y,A2.z,E123.x,E123.y,E123.z);
                A2A3E123=sign_aaae(A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A1.x,A1.y,A1.z,E123.x,E123.y,E123.z);

                A1A2E124=sign_aaae(A1.x,A1.y,A1.z,A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,E124.x,E124.y,E124.z);
                A1A4E124=sign_aaae(A1.x,A1.y,A1.z,A4.x,A4.y,A4.z,A2.x,A2.y,A2.z,E124.x,E124.y,E124.z);
                A2A4E124=sign_aaae(A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,A1.x,A1.y,A1.z,E124.x,E124.y,E124.z);

                A1A3E134=sign_aaae(A1.x,A1.y,A1.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,E134.x,E134.y,E134.z);
                A1A4E134=sign_aaae(A1.x,A1.y,A1.z,A4.x,A4.y,A4.z,A3.x,A3.y,A3.z,E134.x,E134.y,E134.z);
                A3A4E134=sign_aaae(A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A1.x,A1.y,A1.z,E134.x,E134.y,E134.z);

                A2A3E234=sign_aaae(A2.x,A2.y,A2.z,A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,E234.x,E234.y,E234.z);
                A2A4E234=sign_aaae(A2.x,A2.y,A2.z,A4.x,A4.y,A4.z,A3.x,A3.y,A3.z,E234.x,E234.y,E234.z);
                A3A4E234=sign_aaae(A3.x,A3.y,A3.z,A4.x,A4.y,A4.z,A2.x,A2.y,A2.z,E234.x,E234.y,E234.z);

                //cav_vol[i]=cav_vol[i]+S123*A1A2E123*volume_tetrahedron(A1.x,A1.y,A1.z,E123.x,E123.y,E123.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r1);
                //cav_vol[i]=cav_vol[i]+S123*A1A3E123*volume_tetrahedron(A1.x,A1.y,A1.z,E123.x,E123.y,E123.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r1);
                //cav_vol[i]=cav_vol[i]+S123*A1A2E123*volume_tetrahedron(A2.x,A2.y,A2.z,E123.x,E123.y,E123.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r2);
                //cav_vol[i]=cav_vol[i]+S123*A2A3E123*volume_tetrahedron(A2.x,A2.y,A2.z,E123.x,E123.y,E123.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r2);
                //cav_vol[i]=cav_vol[i]+S123*A1A3E123*volume_tetrahedron(A3.x,A3.y,A3.z,E123.x,E123.y,E123.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r3);
                //cav_vol[i]=cav_vol[i]+S123*A2A3E123*volume_tetrahedron(A3.x,A3.y,A3.z,E123.x,E123.y,E123.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r3);

                //cav_vol[i]=cav_vol[i]+S124*A1A2E124*volume_tetrahedron(A1.x,A1.y,A1.z,E124.x,E124.y,E124.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r1);
                //cav_vol[i]=cav_vol[i]+S124*A1A4E124*volume_tetrahedron(A1.x,A1.y,A1.z,E124.x,E124.y,E124.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r1);
                //cav_vol[i]=cav_vol[i]+S124*A1A2E124*volume_tetrahedron(A2.x,A2.y,A2.z,E124.x,E124.y,E124.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r2);
                //cav_vol[i]=cav_vol[i]+S124*A2A4E124*volume_tetrahedron(A2.x,A2.y,A2.z,E124.x,E124.y,E124.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r2);
                //cav_vol[i]=cav_vol[i]+S124*A1A4E124*volume_tetrahedron(A4.x,A4.y,A4.z,E124.x,E124.y,E124.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r4);
                //cav_vol[i]=cav_vol[i]+S124*A2A4E124*volume_tetrahedron(A4.x,A4.y,A4.z,E124.x,E124.y,E124.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r4);

                //cav_vol[i]=cav_vol[i]+S134*A1A3E134*volume_tetrahedron(A1.x,A1.y,A1.z,E134.x,E134.y,E134.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r1);
                //cav_vol[i]=cav_vol[i]+S134*A1A4E134*volume_tetrahedron(A1.x,A1.y,A1.z,E134.x,E134.y,E134.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r1);
                //cav_vol[i]=cav_vol[i]+S134*A1A3E134*volume_tetrahedron(A3.x,A3.y,A3.z,E134.x,E134.y,E134.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r3);
                //cav_vol[i]=cav_vol[i]+S134*A3A4E134*volume_tetrahedron(A3.x,A3.y,A3.z,E134.x,E134.y,E134.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r3);
                //cav_vol[i]=cav_vol[i]+S134*A1A4E134*volume_tetrahedron(A4.x,A4.y,A4.z,E134.x,E134.y,E134.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r4);
                //cav_vol[i]=cav_vol[i]+S134*A3A4E134*volume_tetrahedron(A4.x,A4.y,A4.z,E134.x,E134.y,E134.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r4);

                //cav_vol[i]=cav_vol[i]+S234*A2A3E234*volume_tetrahedron(A2.x,A2.y,A2.z,E234.x,E234.y,E234.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r2);
                //cav_vol[i]=cav_vol[i]+S234*A2A4E234*volume_tetrahedron(A2.x,A2.y,A2.z,E234.x,E234.y,E234.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r2);
                //cav_vol[i]=cav_vol[i]+S234*A2A3E234*volume_tetrahedron(A3.x,A3.y,A3.z,E234.x,E234.y,E234.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r3);
                //cav_vol[i]=cav_vol[i]+S234*A3A4E234*volume_tetrahedron(A3.x,A3.y,A3.z,E234.x,E234.y,E234.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r3);
                //cav_vol[i]=cav_vol[i]+S234*A2A4E234*volume_tetrahedron(A4.x,A4.y,A4.z,E234.x,E234.y,E234.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r4);
                //cav_vol[i]=cav_vol[i]+S234*A3A4E234*volume_tetrahedron(A4.x,A4.y,A4.z,E234.x,E234.y,E234.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r4);

                Atoms[cavity_list[j]->D->AT[0]].vor_vol=Atoms[cavity_list[j]->D->AT[0]].vor_vol+S123*A1A2E123*removesurface(A1.x,A1.y,A1.z,E123.x,E123.y,E123.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r1,1,0);
                Atoms[cavity_list[j]->D->AT[0]].vor_vol=Atoms[cavity_list[j]->D->AT[0]].vor_vol+S123*A1A3E123*removesurface(A1.x,A1.y,A1.z,E123.x,E123.y,E123.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r1,1,0);
                Atoms[cavity_list[j]->D->AT[0]].vor_vol=Atoms[cavity_list[j]->D->AT[0]].vor_vol+S124*A1A2E124*removesurface(A1.x,A1.y,A1.z,E124.x,E124.y,E124.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r1,1,0);
                Atoms[cavity_list[j]->D->AT[0]].vor_vol=Atoms[cavity_list[j]->D->AT[0]].vor_vol+S124*A1A4E124*removesurface(A1.x,A1.y,A1.z,E124.x,E124.y,E124.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r1,1,0);
                Atoms[cavity_list[j]->D->AT[0]].vor_vol=Atoms[cavity_list[j]->D->AT[0]].vor_vol+S134*A1A3E134*removesurface(A1.x,A1.y,A1.z,E134.x,E134.y,E134.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r1,1,0);
                Atoms[cavity_list[j]->D->AT[0]].vor_vol=Atoms[cavity_list[j]->D->AT[0]].vor_vol+S134*A1A4E134*removesurface(A1.x,A1.y,A1.z,E134.x,E134.y,E134.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r1,1,0);

                Atoms[cavity_list[j]->D->AT[1]].vor_vol=Atoms[cavity_list[j]->D->AT[1]].vor_vol+S123*A1A2E123*removesurface(A2.x,A2.y,A2.z,E123.x,E123.y,E123.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r2,1,0);
                Atoms[cavity_list[j]->D->AT[1]].vor_vol=Atoms[cavity_list[j]->D->AT[1]].vor_vol+S123*A2A3E123*removesurface(A2.x,A2.y,A2.z,E123.x,E123.y,E123.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r2,1,0);
                Atoms[cavity_list[j]->D->AT[1]].vor_vol=Atoms[cavity_list[j]->D->AT[1]].vor_vol+S124*A1A2E124*removesurface(A2.x,A2.y,A2.z,E124.x,E124.y,E124.z,B12.x,B12.y,B12.z,Vx,Vy,Vz,r2,1,0);
                Atoms[cavity_list[j]->D->AT[1]].vor_vol=Atoms[cavity_list[j]->D->AT[1]].vor_vol+S124*A2A4E124*removesurface(A2.x,A2.y,A2.z,E124.x,E124.y,E124.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r2,1,0);
                Atoms[cavity_list[j]->D->AT[1]].vor_vol=Atoms[cavity_list[j]->D->AT[1]].vor_vol+S234*A2A3E234*removesurface(A2.x,A2.y,A2.z,E234.x,E234.y,E234.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r2,1,0);
                Atoms[cavity_list[j]->D->AT[1]].vor_vol=Atoms[cavity_list[j]->D->AT[1]].vor_vol+S234*A2A4E234*removesurface(A2.x,A2.y,A2.z,E234.x,E234.y,E234.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r2,1,0);

                Atoms[cavity_list[j]->D->AT[2]].vor_vol=Atoms[cavity_list[j]->D->AT[2]].vor_vol+S123*A1A3E123*removesurface(A3.x,A3.y,A3.z,E123.x,E123.y,E123.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r3,1,0);
                Atoms[cavity_list[j]->D->AT[2]].vor_vol=Atoms[cavity_list[j]->D->AT[2]].vor_vol+S123*A2A3E123*removesurface(A3.x,A3.y,A3.z,E123.x,E123.y,E123.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r3,1,0);
                Atoms[cavity_list[j]->D->AT[2]].vor_vol=Atoms[cavity_list[j]->D->AT[2]].vor_vol+S134*A1A3E134*removesurface(A3.x,A3.y,A3.z,E134.x,E134.y,E134.z,B13.x,B13.y,B13.z,Vx,Vy,Vz,r3,1,0);
                Atoms[cavity_list[j]->D->AT[2]].vor_vol=Atoms[cavity_list[j]->D->AT[2]].vor_vol+S134*A3A4E134*removesurface(A3.x,A3.y,A3.z,E134.x,E134.y,E134.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r3,1,0);
                Atoms[cavity_list[j]->D->AT[2]].vor_vol=Atoms[cavity_list[j]->D->AT[2]].vor_vol+S234*A2A3E234*removesurface(A3.x,A3.y,A3.z,E234.x,E234.y,E234.z,B23.x,B23.y,B23.z,Vx,Vy,Vz,r3,1,0);
                Atoms[cavity_list[j]->D->AT[2]].vor_vol=Atoms[cavity_list[j]->D->AT[2]].vor_vol+S234*A3A4E234*removesurface(A3.x,A3.y,A3.z,E234.x,E234.y,E234.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r3,1,0);

                Atoms[cavity_list[j]->D->AT[3]].vor_vol=Atoms[cavity_list[j]->D->AT[3]].vor_vol+S124*A1A4E124*removesurface(A4.x,A4.y,A4.z,E124.x,E124.y,E124.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r4,1,0);
                Atoms[cavity_list[j]->D->AT[3]].vor_vol=Atoms[cavity_list[j]->D->AT[3]].vor_vol+S124*A2A4E124*removesurface(A4.x,A4.y,A4.z,E124.x,E124.y,E124.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r4,1,0);
                Atoms[cavity_list[j]->D->AT[3]].vor_vol=Atoms[cavity_list[j]->D->AT[3]].vor_vol+S134*A1A4E134*removesurface(A4.x,A4.y,A4.z,E134.x,E134.y,E134.z,B14.x,B14.y,B14.z,Vx,Vy,Vz,r4,1,0);
                Atoms[cavity_list[j]->D->AT[3]].vor_vol=Atoms[cavity_list[j]->D->AT[3]].vor_vol+S134*A3A4E134*removesurface(A4.x,A4.y,A4.z,E134.x,E134.y,E134.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r4,1,0);
                Atoms[cavity_list[j]->D->AT[3]].vor_vol=Atoms[cavity_list[j]->D->AT[3]].vor_vol+S234*A2A4E234*removesurface(A4.x,A4.y,A4.z,E234.x,E234.y,E234.z,B24.x,B24.y,B24.z,Vx,Vy,Vz,r4,1,0);
                Atoms[cavity_list[j]->D->AT[3]].vor_vol=Atoms[cavity_list[j]->D->AT[3]].vor_vol+S234*A3A4E234*removesurface(A4.x,A4.y,A4.z,E234.x,E234.y,E234.z,B34.x,B34.y,B34.z,Vx,Vy,Vz,r4,1,0);

            }
          uni_sphere=0.;
          sum_sphere=0.;
          //int borat=0;
            for(int i=0; i<nAtoms; i++)
            {
                if(Atoms[i].radius!=1.0)
                {
                    uni_sphere=uni_sphere+Atoms[i].vor_vol;
                    sum_sphere=sum_sphere+4./3.*M_PI*powl(Atoms[i].radius+r_cut,3);
                    //cout<<Atoms[i].p.x<<"\t"<<Atoms[i].p.y<<"\t"<<Atoms[i].p.z<<"\t"<<Atoms[i].radius+r_cut<<"\n";
                    cout<<i<<"\t"<<4./3.*M_PI*powl(Atoms[i].radius+r_cut,3)<<"\t"<<Atoms[i].vor_vol<<"\n";
                }
                else
                {
                    borat=borat+1;
                }
                //cout<<4./3.*M_PI*powl(Atoms[i].radius,3)<<"\t"<<Atoms[i].vor_vol<<"\n";
            }
            cout<<std::setprecision(16);
            cout<<resno<<"\t"<<uni_sphere<<"\t"<<convex_vol<<"\t"<<sum_sphere<<"\t"<<uni_sphere/convex_vol<<"\t"<<borat<<"\t"<<borat*1./nAtoms<<"\n";

            //double cav_tot=0.;
            //double ca_per_tot=0.;
            //for(int i=0; i<void_vert_count; i++)
            //{
            //    if(cav_vol[i])
            //        cout<<i<<"\t"<<cav_vol[i]<<"\t"<<resno<<"\t"<<pocket[i]<<"\n";
            //    cav_tot=cav_tot+cav_vol[i];
            //    ca_per_tot=ca_per_tot+cav_area[i];
            //    //cout<<i<<"\t"<<cav_area[i]<<"\t"<<cav_lenght[i]<<"\n";
            //}
            //cout<<r_cut<<"\t"<<cav_tot<<"\t"<<ca_per_tot<<"\n";
            delete [] pocket;
            delete [] cav_vol;
            delete [] cav_area;
            delete [] cavity_list;
            return 0;
        }

        //for(int t=0; t<ntypes; t++)
        {
            container_vertice *cstart=nullptr;
            vertice *temp_start=nullptr;
            temp_start=V->initial;
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
        }
        delete [] start;
        delunay *D;
        delunay *temp;
        //for(int t=0; t<ntypes; t++)
        {
            temp=FULLSETD.initial;
            while(1)
            {
                D=temp;
                if(temp->next)
                {
                    temp=temp->next;
                    delete D;
                }
                else
                {
                    delete D;
                    break;
                }
            }
            FULLSETD.initial=nullptr;
            FULLSETD.end=nullptr;
        }
        for(int n=0; n<nAtoms; n++)
        {
            //delete[] Atoms[n].Cstart;
            container_delunay *D=nullptr;
            //for(int t=0; t<ntypes; t++)
            {
                D=Atoms[n].D_FIRST;
                while(1)
                {
                    container_delunay *temp=nullptr;
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
        }
    }
    delete_everything(Atoms,nAtoms,ntypes);
    delete[] Atoms;
    delete[] radius;
}
