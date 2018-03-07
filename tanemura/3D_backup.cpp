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
    int AT[4];
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
}******DSET;
struct container_delunay
{
    delunay *D;
    struct container_delunay *next=nullptr;
    struct container_delunay *prev=nullptr;
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
        {
            if(!compare(focus->neib_vert[i]->p,add->p))
            {
                std::vector<int> EV1 {focus->neib_vert[i]->A,focus->neib_vert[i]->D->a,focus->neib_vert[i]->D->b,focus->neib_vert[i]->D->c};
                std::vector<int> v1 {add->A,add->D->a,add->D->b,add->D->c};
                std::sort(EV1.begin(),EV1.end());
                std::sort(v1.begin(),v1.end());
                if(EV1[0]==v1[0] && EV1[1]==v1[1] && EV1[2]==v1[2] && EV1[3]==v1[3])
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
    int *D3bondinvoid[50][50];
    int *edge_index[500];
    int neighbours=0;
    int *conti;
    site *MIDP[100][100];
    struct face *F=NULL;
    long double radius=1.;
    int ignore=0;
    int index;
    container_vertice **Cstart=nullptr;
    set_of_delunay *D=nullptr;
    container_delunay *D_FIRST;
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
    int flag_dup;
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

    if(v->p->y>LAYER_CUT)
    {
        if(EV->prev==NULL)
        {
            vertice *temp_vert=EV;
            while(1)
            {
                if(!compare(temp_vert->p,v->p))
                {
                    delete v->p;
                    delete v;
                    return temp_vert;
                }
                if(temp_vert->next)
                    temp_vert=temp_vert->next;
                else
                    break;
            }
        }
    }
    if(flag==0)
    {
        if(debug)
        {
            cout<<"here\n";
        }
        std::vector<int> EV1 {EV->A,EV->D->a,EV->D->b,EV->D->c};
        std::vector<int> v1 {v->A,v->D->a,v->D->b,v->D->c};
        std::sort(EV1.begin(),EV1.end());
        std::sort(v1.begin(),v1.end());
        if(debug)
        {
            cout<<"here\n";
            cout<<EV1[0]<<"\t"<<EV1[1]<<"\t"<<EV1[2]<<"\n";
            cout<<v1[0]<<"\t"<<v1[1]<<"\t"<<v1[2]<<"\n";
        }
        if(EV1[0]==v1[0] && EV1[1]==v1[1] && EV1[2]==v1[2] && EV1[3]==v1[3])
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
  //for(int i=0; i<4; i++)
  //{
  //    cout<<D->AT[i]<<"\t";
  //}
  //for(int i=0;i<4;i++)
  //	for(int j=0;j<4;j++)
  //		for(int k=0;k<4;k++)
  //			for(uu
  //cout<<"\n";
    //return ;
    cout<<"#\t";
    cout<<D->AT[1]<<"\t"<<D->AT[2]<<"\t"<<D->AT[3]<<"\n";
    cout<<"#\t";
    cout<<D->A<<"\t"<<D->B<<"\t"<<D->C<<"\n";
	ATOM=&(Atoms[D->AT[0]]);
    long double Sx,Sy,Sz;
    long double Px,Py,Pz;
    Sx=ATOM->x-Atoms[D->AT[1]].x;
    Sy=ATOM->y-Atoms[D->AT[1]].y;
    Sz=ATOM->z-Atoms[D->AT[1]].z;
    Sx=(Sx-(tilt*lroundl(Sy/twob)));
    Sx=(Sx-(twob*lroundl(Sx/twob)));
    Sy=(Sy-(twob*lroundl(Sy/twob)));
    Sz=(Sz-(twob*lroundl(Sz/twob)));
    cout<<"draw line\t{";
    cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<ATOM->z<<"}\t{"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\t"<<ATOM->z-Sz<<"}\twidth	1\n";
    Sx=ATOM->x-Atoms[D->AT[2]].x;
    Sy=ATOM->y-Atoms[D->AT[2]].y;
    Sz=ATOM->z-Atoms[D->AT[2]].z;
    Sx=(Sx-(tilt*lroundl(Sy/twob)));
    Sx=(Sx-(twob*lroundl(Sx/twob)));
    Sy=(Sy-(twob*lroundl(Sy/twob)));
    Sz=(Sz-(twob*lroundl(Sz/twob)));
    cout<<"draw line\t{";
    cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<ATOM->z<<"}\t{"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\t"<<ATOM->z-Sz<<"}\twidth	1\n";
    Sx=ATOM->x-Atoms[D->AT[3]].x;
    Sy=ATOM->y-Atoms[D->AT[3]].y;
    Sz=ATOM->z-Atoms[D->AT[3]].z;
    Sx=(Sx-(tilt*lroundl(Sy/twob)));
    Sx=(Sx-(twob*lroundl(Sx/twob)));
    Sy=(Sy-(twob*lroundl(Sy/twob)));
    Sz=(Sz-(twob*lroundl(Sz/twob)));
    cout<<"draw line\t{";
    cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<ATOM->z<<"}\t{"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\t"<<ATOM->z-Sz<<"}\twidth	1\n";
    Sx=ATOM->x-Atoms[D->AT[2]].x;
    Sy=ATOM->y-Atoms[D->AT[2]].y;
    Sz=ATOM->z-Atoms[D->AT[2]].z;
    Sx=(Sx-(tilt*lroundl(Sy/twob)));
    Sx=(Sx-(twob*lroundl(Sx/twob)));
    Sy=(Sy-(twob*lroundl(Sy/twob)));
    Sz=(Sz-(twob*lroundl(Sz/twob)));
    Px=ATOM->x-Atoms[D->AT[1]].x;
    Py=ATOM->y-Atoms[D->AT[1]].y;
    Pz=ATOM->z-Atoms[D->AT[1]].z;
    Px=(Px-(tilt*lroundl(Py/twob)));
    Px=(Px-(twob*lroundl(Px/twob)));
    Py=(Py-(twob*lroundl(Py/twob)));
    Pz=(Pz-(twob*lroundl(Pz/twob)));
    cout<<"draw line\t{";
    cout<<ATOM->x-Px<<"\t"<<ATOM->y-Py<<"\t"<<ATOM->z-Pz<<"}\t{"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\t"<<ATOM->z-Sz<<"}\twidth	1\n";
    Sx=ATOM->x-Atoms[D->AT[3]].x;
    Sy=ATOM->y-Atoms[D->AT[3]].y;
    Sz=ATOM->z-Atoms[D->AT[3]].z;
    Sx=(Sx-(tilt*lroundl(Sy/twob)));
    Sx=(Sx-(twob*lroundl(Sx/twob)));
    Sy=(Sy-(twob*lroundl(Sy/twob)));
    Sz=(Sz-(twob*lroundl(Sz/twob)));
    Px=ATOM->x-Atoms[D->AT[1]].x;
    Py=ATOM->y-Atoms[D->AT[1]].y;
    Pz=ATOM->z-Atoms[D->AT[1]].z;
    Px=(Px-(tilt*lroundl(Py/twob)));
    Px=(Px-(twob*lroundl(Px/twob)));
    Py=(Py-(twob*lroundl(Py/twob)));
    Pz=(Pz-(twob*lroundl(Pz/twob)));
    cout<<"draw line\t{";
    cout<<ATOM->x-Px<<"\t"<<ATOM->y-Py<<"\t"<<ATOM->z-Pz<<"}\t{"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\t"<<ATOM->z-Sz<<"}\twidth	1\n";
    Sx=ATOM->x-Atoms[D->AT[2]].x;
    Sy=ATOM->y-Atoms[D->AT[2]].y;
    Sz=ATOM->z-Atoms[D->AT[2]].z;
    Sx=(Sx-(tilt*lroundl(Sy/twob)));
    Sx=(Sx-(twob*lroundl(Sx/twob)));
    Sy=(Sy-(twob*lroundl(Sy/twob)));
    Sz=(Sz-(twob*lroundl(Sz/twob)));
    Px=ATOM->x-Atoms[D->AT[3]].x;
    Py=ATOM->y-Atoms[D->AT[3]].y;
    Pz=ATOM->z-Atoms[D->AT[3]].z;
    Px=(Px-(tilt*lroundl(Py/twob)));
    Px=(Px-(twob*lroundl(Px/twob)));
    Py=(Py-(twob*lroundl(Py/twob)));
    Pz=(Pz-(twob*lroundl(Pz/twob)));
    cout<<"draw line\t{";
    cout<<ATOM->x-Px<<"\t"<<ATOM->y-Py<<"\t"<<ATOM->z-Pz<<"}\t{"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\t"<<ATOM->z-Sz<<"}\twidth	1\n";
}
void create_delunay(atom Atoms[],int A1,int A2,int A3,int A4,delunay *D,int TYPE)
{
////D->AT[0]=A1;
////D->AT[1]=A2;
////D->AT[2]=A3;
////D->AT[3]=A4;
    for(int k=0; k<4; k++)
    {
        atom *ATOM;
        ATOM=&(Atoms[D->AT[k]]);
        int a1=-1,a2=-1,a3=-1;
        //cout<<D->AT[k]<<"\n";
        //for(int add=1;add<4;add++)
        {
            int i;
            for(i=0; i<ATOM->conti[TYPE]; i++)
            {
                if(ATOM->contigous[i][TYPE]==D->AT[(k+1)%4])
                {
                    a1=i;
                }
                if(ATOM->contigous[i][TYPE]==D->AT[(k+2)%4])
                {
                    a2=i;
                }
                if(ATOM->contigous[i][TYPE]==D->AT[(k+3)%4])
                {
                    a3=i;
                }
            }
            if(a1==-1)
            {
                ATOM->contigous[i][TYPE]=D->AT[(k+1)%4];
                a1=i;
                i++;
                ATOM->conti[TYPE]++;
            }
            if(a2==-1)
            {
                ATOM->contigous[i][TYPE]=D->AT[(k+2)%4];
                a2=i;
                i++;
                ATOM->conti[TYPE]++;
            }
            if(a3==-1)
            {
                ATOM->contigous[i][TYPE]=D->AT[(k+3)%4];
                a3=i;
                i++;
                ATOM->conti[TYPE]++;
            }
			//cout<<a1<<"\t"<<a2<<"\t"<<a3<<"\n";
			//cout<<ATOM->part_c[a1][a2]<<"\n";
            if(ATOM->part_c[a1][a2][TYPE]==0)
            {
					//cout<<"here\n";
                ATOM->edge_index[a1][TYPE]++;
                ATOM->edge_index[a2][TYPE]++;
                ////ATOM->part_c[a1][a2]++;
                ////ATOM->part_c[a2][a1]++;
            }
            if(ATOM->part_c[a1][a3][TYPE]==0)
            {
                ATOM->edge_index[a1][TYPE]++;
                ATOM->edge_index[a3][TYPE]++;
                ////ATOM->part_c[a1][a3]++;
                ////ATOM->part_c[a3][a1]++;
            }
            if(ATOM->part_c[a2][a3][TYPE]==0)
            {
                ATOM->edge_index[a2][TYPE]++;
                ATOM->edge_index[a3][TYPE]++;
                ////ATOM->part_c[a2][a3]++;
                ////ATOM->part_c[a3][a2]++;
            }

            ATOM->part_c[a1][a2][TYPE]++;
            ATOM->part_c[a2][a1][TYPE]++;

            ATOM->part_c[a1][a3][TYPE]++;
            ATOM->part_c[a3][a1][TYPE]++;

            ATOM->part_c[a2][a3][TYPE]++;
            ATOM->part_c[a3][a2][TYPE]++;

            if(!ATOM->D_FIRST)
            {
                ATOM->D_FIRST=new container_delunay;
                ATOM->D_FIRST->D=D;
            }
            else
            {
                container_delunay *temp;
                temp=ATOM->D_FIRST;
                //temp=ATOM->D_FIRST->next;
                while(temp->next)
                {
					//cout<<temp->D->AT[0]<<"\t"<<temp->D->AT[1]<<"\t"<<temp->D->AT[2]<<"\t"<<temp->D->AT[3]<<"\n";
                    temp=temp->next;
                }
			//	cout<<"here\n";
                temp->next=new container_delunay;
                temp->next->D=D;
            }
        }
    }
    //ATOM->contigous[ATOM->conti[TYPE]][TYPE]=DIS_atom;
    ////ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
    //ATOM->edge_index[0][TYPE]++;
    //ATOM->edge_index[0][TYPE]++;
    //ATOM->edge_index[1][TYPE]++;
    //ATOM->edge_index[1][TYPE]++;
    //ATOM->edge_index[2][TYPE]++;
    //ATOM->edge_index[2][TYPE]++;
    //ATOM->conti[TYPE]++;
    //ATOM->D[TYPE].initial= new delunay;
///////////////////////////////////////////////////////////////////////////////////////////////
    ////if(!ATOM->D_FIRST)
    //{
    //	ATOM->D_FIRST= new container_delunay;
    //	ATOM->D_FIRST->D=ATOM->D[TYPE].initial;
    //}
    //if(!Atoms[ATOM->contigous[0][TYPE]].D_FIRST)
    //{
    //	Atoms[ATOM->contigous[0][TYPE]].D_FIRST=new container_delunay;
    //	Atoms[ATOM->contigous[0][TYPE]].D_FIRST->D=ATOM->D[TYPE].initial;
    //}
    //else
    //{
    //	container_delunay *temp;
    //	temp=Atoms[ATOM->contigous[0][TYPE]].D_FIRST;
    //	while(temp->next)
    //	{
    //		temp=temp->next;
    //	}
    //	temp->next=new container_delunay;
    //	temp->D=ATOM->D[TYPE].initial;
    //}
    //if(!Atoms[ATOM->contigous[1][TYPE]].D_FIRST)
    //{
    //	Atoms[ATOM->contigous[1][TYPE]].D_FIRST=new container_delunay;
    //	Atoms[ATOM->contigous[1][TYPE]].D_FIRST->D=ATOM->D[TYPE].initial;
    //}
    //else
    //{
    //	container_delunay *temp;
    //	temp=Atoms[ATOM->contigous[1][TYPE]].D_FIRST;
    //	while(temp->next)
    //	{
    //		temp=temp->next;
    //	}
    //	temp->next=new container_delunay;
    //	temp->D=ATOM->D[TYPE].initial;
    //}
    //if(!Atoms[ATOM->contigous[2][TYPE]].D_FIRST)
    //{
    //	Atoms[ATOM->contigous[2][TYPE]].D_FIRST=new container_delunay;
    //	Atoms[ATOM->contigous[2][TYPE]].D_FIRST->D=ATOM->D[TYPE].initial;
    //}
    //else
    //{
    //	container_delunay *temp;
    //	temp=Atoms[ATOM->contigous[2][TYPE]].D_FIRST;
    //	while(temp->next)
    //	{
    //		temp=temp->next;
    //	}
    //	temp->next=new container_delunay;
    //	temp->D=ATOM->D[TYPE].initial;
    //}
    //ATOM->D[TYPE].initial->A=0;
    //ATOM->D[TYPE].initial->B=1;
    //ATOM->MIDP[0][1][TYPE].x=circx_s+L.x;
    //ATOM->MIDP[0][1][TYPE].y=circy_s+L.y;
    //ATOM->MIDP[0][1][TYPE].z=circz_s+L.z;
    //ATOM->MIDP[1][0][TYPE].x=circx_s+L.x;
    //ATOM->MIDP[1][0][TYPE].y=circy_s+L.y;
    //ATOM->MIDP[1][0][TYPE].z=circz_s+L.z;
////ATOM->D[TYPE].initial->ABx=circx_s+L.x;
////ATOM->D[TYPE].initial->ABy=circy_s+L.y;
////ATOM->D[TYPE].initial->ABz=circz_s+L.z;
    //ATOM->D[TYPE].initial->C=2;
    //ATOM->D[TYPE].initial->ABf=1;
    //ATOM->D[TYPE].initial->BCf=1;
    //ATOM->D[TYPE].initial->CAf=1;
    //ATOM->part_c[0][1][TYPE]=1;
    //ATOM->part_c[0][2][TYPE]=1;
    //ATOM->part_c[1][2][TYPE]=1;
    //ATOM->part_c[1][0][TYPE]=1;
    //ATOM->part_c[2][1][TYPE]=1;
    //ATOM->part_c[2][0][TYPE]=1;
    //ATOM->D[TYPE].initial->a=ATOM->contigous[0][TYPE];
    //ATOM->D[TYPE].initial->b=ATOM->contigous[1][TYPE];
    //ATOM->D[TYPE].initial->c=ATOM->contigous[2][TYPE];
    ////cout<<ATOM->D[TYPE].initial->a<<"\t"<<ATOM->D[TYPE].initial->b<<"\n";
    //ATOM->D[TYPE].initial->circum_x=circx;
    //ATOM->D[TYPE].initial->circum_y=circy;
    //ATOM->D[TYPE].initial->circum_z=circz;
    //ATOM->D[TYPE].initial->Ax=midax;
    //ATOM->D[TYPE].initial->Ay=miday;
    //ATOM->D[TYPE].initial->Az=midaz;
    //ATOM->D[TYPE].initial->Bx=midbx;
    //ATOM->D[TYPE].initial->By=midby;
    //ATOM->D[TYPE].initial->Bz=midbz;
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
    int A1,A2,A3,A4;
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
    A1=ATOM->index;
    A2=nearest;
    //ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
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
            atom M=Atoms[nearest];
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
    A3=DIS_atom;
    //ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
    //ATOM->edge_index[ATOM->conti[TYPE]-1][TYPE]++;
    DIS_MIN=box*box*box;
    long double circx_s=circx;
    long double circy_s=circy;
    long double circz_s=circz;
//	cout<<DIS_atom<<"\n";
    atom L=*ATOM;
    atom M=Atoms[nearest];
    atom N=Atoms[DIS_atom];
////cout<<M.x<<"\t"<<M.y<<"\t"<<M.z<<"\n";
////cout<<N.x<<"\t"<<N.y<<"\t"<<N.z<<"\n";
////cout<<L.x<<"\t"<<L.y<<"\t"<<L.z<<"\n";
    for(int i=0; i<ATOM->neighbours; i++)
    {
        if(ATOM->neighlist[i]!=nearest && ATOM->neighlist[i]!=A3)
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
    //ATOM->contigous[ATOM->conti[TYPE]][TYPE]=DIS_atom;
    A4=DIS_atom;
    //cout<<A1<<"\t"<<A2<<"\t"<<A3<<"\t"<<A4<<"\n";
    std::vector<int> EV1 {A1,A2,A3,A4};
    std::sort(EV1.begin(),EV1.end());
    delunay *D;
////DSET[EV1[0]][EV1[1]][EV1[2]][EV1[3]][TYPE]=*D;
////DSET[EV1[0]][EV1[1]][EV1[2]][EV1[3]][TYPE]=new delunay;
    D=new delunay;
    D->AT[0]=EV1[0];
    D->AT[1]=EV1[1];
    D->AT[2]=EV1[2];
    D->AT[3]=EV1[3];
    create_delunay(Atoms,A1,A2,A3,A4,D,TYPE);
    //ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
    //ATOM->edge_index[0][TYPE]++;
    //ATOM->edge_index[0][TYPE]++;
    //ATOM->edge_index[1][TYPE]++;
    //ATOM->edge_index[1][TYPE]++;
    //ATOM->edge_index[2][TYPE]++;
    //ATOM->edge_index[2][TYPE]++;
    //ATOM->conti[TYPE]++;
    //ATOM->D[TYPE].initial= new delunay;
///////////////////////////////////////////////////////////////////////////////////////////////
    ////if(!ATOM->D_FIRST)
    //{
    //	ATOM->D_FIRST= new container_delunay;
    //	ATOM->D_FIRST->D=ATOM->D[TYPE].initial;
    //}
    //if(!Atoms[ATOM->contigous[0][TYPE]].D_FIRST)
    //{
    //	Atoms[ATOM->contigous[0][TYPE]].D_FIRST=new container_delunay;
    //	Atoms[ATOM->contigous[0][TYPE]].D_FIRST->D=ATOM->D[TYPE].initial;
    //}
    //else
    //{
    //	container_delunay *temp;
    //	temp=Atoms[ATOM->contigous[0][TYPE]].D_FIRST;
    //	while(temp->next)
    //	{
    //		temp=temp->next;
    //	}
    //	temp->next=new container_delunay;
    //	temp->D=ATOM->D[TYPE].initial;
    //}
    //if(!Atoms[ATOM->contigous[1][TYPE]].D_FIRST)
    //{
    //	Atoms[ATOM->contigous[1][TYPE]].D_FIRST=new container_delunay;
    //	Atoms[ATOM->contigous[1][TYPE]].D_FIRST->D=ATOM->D[TYPE].initial;
    //}
    //else
    //{
    //	container_delunay *temp;
    //	temp=Atoms[ATOM->contigous[1][TYPE]].D_FIRST;
    //	while(temp->next)
    //	{
    //		temp=temp->next;
    //	}
    //	temp->next=new container_delunay;
    //	temp->D=ATOM->D[TYPE].initial;
    //}
    //if(!Atoms[ATOM->contigous[2][TYPE]].D_FIRST)
    //{
    //	Atoms[ATOM->contigous[2][TYPE]].D_FIRST=new container_delunay;
    //	Atoms[ATOM->contigous[2][TYPE]].D_FIRST->D=ATOM->D[TYPE].initial;
    //}
    //else
    //{
    //	container_delunay *temp;
    //	temp=Atoms[ATOM->contigous[2][TYPE]].D_FIRST;
    //	while(temp->next)
    //	{
    //		temp=temp->next;
    //	}
    //	temp->next=new container_delunay;
    //	temp->D=ATOM->D[TYPE].initial;
    //}
    //ATOM->D[TYPE].initial->A=0;
    //ATOM->D[TYPE].initial->B=1;
    //ATOM->MIDP[0][1][TYPE].x=circx_s+L.x;
    //ATOM->MIDP[0][1][TYPE].y=circy_s+L.y;
    //ATOM->MIDP[0][1][TYPE].z=circz_s+L.z;
    //ATOM->MIDP[1][0][TYPE].x=circx_s+L.x;
    //ATOM->MIDP[1][0][TYPE].y=circy_s+L.y;
    //ATOM->MIDP[1][0][TYPE].z=circz_s+L.z;
////ATOM->D[TYPE].initial->ABx=circx_s+L.x;
////ATOM->D[TYPE].initial->ABy=circy_s+L.y;
////ATOM->D[TYPE].initial->ABz=circz_s+L.z;
    //ATOM->D[TYPE].initial->C=2;
    //ATOM->D[TYPE].initial->ABf=1;
    //ATOM->D[TYPE].initial->BCf=1;
    //ATOM->D[TYPE].initial->CAf=1;
    //ATOM->part_c[0][1][TYPE]=1;
    //ATOM->part_c[0][2][TYPE]=1;
    //ATOM->part_c[1][2][TYPE]=1;
    //ATOM->part_c[1][0][TYPE]=1;
    //ATOM->part_c[2][1][TYPE]=1;
    //ATOM->part_c[2][0][TYPE]=1;
    //ATOM->D[TYPE].initial->a=ATOM->contigous[0][TYPE];
    //ATOM->D[TYPE].initial->b=ATOM->contigous[1][TYPE];
    //ATOM->D[TYPE].initial->c=ATOM->contigous[2][TYPE];
    ////cout<<ATOM->D[TYPE].initial->a<<"\t"<<ATOM->D[TYPE].initial->b<<"\n";
    //ATOM->D[TYPE].initial->circum_x=circx;
    //ATOM->D[TYPE].initial->circum_y=circy;
    //ATOM->D[TYPE].initial->circum_z=circz;
    //ATOM->D[TYPE].initial->Ax=midax;
    //ATOM->D[TYPE].initial->Ay=miday;
    //ATOM->D[TYPE].initial->Az=midaz;
    //ATOM->D[TYPE].initial->Bx=midbx;
    //ATOM->D[TYPE].initial->By=midby;
    //ATOM->D[TYPE].initial->Bz=midbz;
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
    long double X1_s,Y1_s,Z1_s;
    long double X2_s,Y2_s,Z2_s;
    long double X3_s,Y3_s,Z3_s;
	int A1,A2,A3,A4;
    int DIS_atom;
	A1=ATOM->index;
	A2=ATOM->contigous[A][TYPE];
	A3=ATOM->contigous[B][TYPE];
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
                for(int p=0; p<ATOM->conti[TYPE]; p++)
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
            long double v3x,v3y,v3z;
            v1x=yA*zB-zA*yB;
            v1y=zA*xB-xA*zB;
            v1z=xA*yB-yA*xB;
            v2x=yC*zA-zC*yA;
            v2y=zC*xA-xC*zA;
            v2z=xC*yA-yC*xA;
            v3x=yC*zB-zC*yB;
            v3y=zC*xB-xC*zB;
            v3z=xC*yB-yC*xB;
            long double a1,b1,c1;
            long double a2,b2,c2;
            long double a3,b3,c3;
            long double a4,b4,c4;
            long double a5,b5,c5;
            long double a6,b6,c6;
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
            a5=v3y*zC-v3z*yC;
            b5=v3z*xC-v3x*zC;
            c5=v3x*yC-v3y*xC;
            a6=v3y*zB-v3z*yB;
            b6=v3z*xB-v3x*zB;
            c6=v3x*yB-v3y*xB;
            long double t1,t2,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3;
            t2=(b1*xA-a1*yA-b1*xB+a1*yB)/(b1*a2-a1*b2);
            X1=xB+a2*t2;
            Y1=yB+b2*t2;
            Z1=zB+c2*t2;
            t2=(b3*xC-a3*yC-b3*xA+a3*yA)/(b3*a4-a3*b4);
            ////cout<<t2<<" wdfd\n";
            X2=xA+a4*t2;
            Y2=yA+b4*t2;
            Z2=zA+c4*t2;
            t2=(b5*xC-a5*yC-b5*xB+a5*yB)/(b5*a6-a5*b6);
            X3=xB+a6*t2;
            Y3=yB+b6*t2;
            Z3=zB+c6*t2;
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
                X1_s=X1+a;
                Y1_s=Y1+b;
                Z1_s=Z1+c;
                X2_s=X2+a;
                Y2_s=Y2+b;
                Z2_s=Z2+c;
                X3_s=X3+a;
                Y3_s=Y3+b;
                Z3_s=Z3+c;
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
  //int flag=1;
////if(lmin<0.)
////{
////	ATOM->ignore=1;
////}

    A4=DIS_atom;
    //cout<<A1<<"\t"<<A2<<"\t"<<A3<<"\t"<<A4<<"\n";
    std::vector<int> EV1 {A1,A2,A3,A4};
    std::sort(EV1.begin(),EV1.end());
    //delunay *D;
////DSET[EV1[0]][EV1[1]][EV1[2]][EV1[3]][TYPE]=*D;
////DSET[EV1[0]][EV1[1]][EV1[2]][EV1[3]][TYPE]=new delunay;
    D=new delunay;
    D->AT[0]=EV1[0];
    D->AT[1]=EV1[1];
    D->AT[2]=EV1[2];
    D->AT[3]=EV1[3];
    create_delunay(Atoms,A1,A2,A3,A4,D,TYPE);
  //int k;
  ////cout<<ATOM->edge_index[12][TYPE]<<" twelve\n";
  //for(k=0; k<ATOM->conti[TYPE]; k++)
  //{
  //    if(ATOM->contigous[k][TYPE]==DIS_atom)
  //    {
  //        flag=0;
  //        break;
  //    }
  //}
  //if(debug)
  //{
  //    cout<<"debug\t";
  //    cout<<flag<<"\t"<<k<<"\n";;
  //    cout<<DIS_atom<<"\n";
  //}
  ////cout<<"#"<<A<<"\t"<<B<<"\t"<<k<<"\n";
  ////cout<<"#"<<X1_s<<"\t"<<Y1_s<<"\t"<<Z1_s<<"\n";
  ////cout<<"#"<<X2_s<<"\t"<<Y2_s<<"\t"<<Z2_s<<"\n";

  //if(flag)
  //{
  //    //ATOM->bondinvoid[ATOM->conti[TYPE]][TYPE]=binv;
  //    //		cout<<"hereag\n";
  //    ATOM->contigous[ATOM->conti[TYPE]][TYPE]=DIS_atom;
  //    ATOM->edge_index[A][TYPE]++;
  //    ATOM->edge_index[B][TYPE]++;
  //    ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
  //    ATOM->edge_index[ATOM->conti[TYPE]][TYPE]++;
  //    ATOM->part_c[ATOM->conti[TYPE]][A][TYPE]++;
  //    ATOM->part_c[A][B][TYPE]++;
  //    ATOM->part_c[B][ATOM->conti[TYPE]][TYPE]++;
  //    ATOM->part_c[A][ATOM->conti[TYPE]][TYPE]++;
  //    ATOM->part_c[B][A][TYPE]++;
  //    ATOM->part_c[ATOM->conti[TYPE]][B][TYPE]++;
  //    ATOM->MIDP[A][B][TYPE].x=X2_s;
  //    ATOM->MIDP[A][B][TYPE].y=Y2_s;
  //    ATOM->MIDP[A][B][TYPE].z=Z2_s;
  //    ATOM->MIDP[B][A][TYPE].x=X2_s;
  //    ATOM->MIDP[B][A][TYPE].y=Y2_s;
  //    ATOM->MIDP[B][A][TYPE].z=Z2_s;
  //    ATOM->MIDP[A][ATOM->conti[TYPE]][TYPE].x=X1_s;
  //    ATOM->MIDP[A][ATOM->conti[TYPE]][TYPE].y=Y1_s;
  //    ATOM->MIDP[A][ATOM->conti[TYPE]][TYPE].z=Z1_s;
  //    ATOM->MIDP[ATOM->conti[TYPE]][A][TYPE].x=X1_s;
  //    ATOM->MIDP[ATOM->conti[TYPE]][A][TYPE].y=Y1_s;
  //    ATOM->MIDP[ATOM->conti[TYPE]][A][TYPE].z=Z1_s;
  //    ATOM->MIDP[B][ATOM->conti[TYPE]][TYPE].x=X3_s;
  //    ATOM->MIDP[B][ATOM->conti[TYPE]][TYPE].y=Y3_s;
  //    ATOM->MIDP[B][ATOM->conti[TYPE]][TYPE].z=Z3_s;
  //    ATOM->MIDP[ATOM->conti[TYPE]][B][TYPE].x=X3_s;
  //    ATOM->MIDP[ATOM->conti[TYPE]][B][TYPE].y=Y3_s;
  //    ATOM->MIDP[ATOM->conti[TYPE]][B][TYPE].z=Z3_s;
  //    //ATOM->MIDP[A][B][TYPE]->
  //}
  //else
  //{
  //    if(ATOM->part_c[A][k][TYPE]==0)
  //    {
  //        ATOM->edge_index[A][TYPE]++;
  //        ATOM->edge_index[k][TYPE]++;
  //    }
  //    if(ATOM->part_c[B][k][TYPE]==0)
  //    {
  //        ATOM->edge_index[B][TYPE]++;
  //        ATOM->edge_index[k][TYPE]++;
  //    }
  //    ATOM->part_c[k][A][TYPE]++;
  //    ATOM->part_c[A][B][TYPE]++;
  //    ATOM->part_c[B][k][TYPE]++;
  //    ATOM->part_c[A][k][TYPE]++;
  //    ATOM->part_c[B][A][TYPE]++;
  //    ATOM->part_c[k][B][TYPE]++;
  //    ATOM->MIDP[A][B][TYPE].x=X2_s;
  //    ATOM->MIDP[A][B][TYPE].y=Y2_s;
  //    ATOM->MIDP[A][B][TYPE].z=Z2_s;
  //    ATOM->MIDP[B][A][TYPE].x=X2_s;
  //    ATOM->MIDP[B][A][TYPE].y=Y2_s;
  //    ATOM->MIDP[B][A][TYPE].z=Z2_s;
  //    ATOM->MIDP[A][k][TYPE].x=X1_s;
  //    ATOM->MIDP[A][k][TYPE].y=Y1_s;
  //    ATOM->MIDP[A][k][TYPE].z=Z1_s;
  //    ATOM->MIDP[k][A][TYPE].x=X1_s;
  //    ATOM->MIDP[k][A][TYPE].y=Y1_s;
  //    ATOM->MIDP[k][A][TYPE].z=Z1_s;
  //    ATOM->MIDP[B][k][TYPE].x=X3_s;
  //    ATOM->MIDP[B][k][TYPE].y=Y3_s;
  //    ATOM->MIDP[B][k][TYPE].z=Z3_s;
  //    ATOM->MIDP[k][B][TYPE].x=X3_s;
  //    ATOM->MIDP[k][B][TYPE].y=Y3_s;
  //    ATOM->MIDP[k][B][TYPE].z=Z3_s;
  //    //ATOM->edge_index[k][TYPE]++;
  //}
  ////cout<<"look\t";
  ////cout<<ATOM->edge_index[8][TYPE]<<"\t"<<A<<"\n";;
  ////cout<<ATOM->part_c[12][10][TYPE]<<" track this\n";
  //delunay *temp=nullptr;
  //temp=D;
  //while(1)
  //{
  //    if(temp->next)
  //        temp=temp->next;
  //    else
  //        break;

  //}
  //temp->next= new delunay;
  //if(flag)
  //{
  //    temp->next->A=ATOM->conti[TYPE];
  //    temp->next->a=ATOM->contigous[ATOM->conti[TYPE]][TYPE];
  //    ////temp->next->Ax=midbx;
  //    ////temp->next->Ay=midby;
  //}
  //else
  //{
  //    temp->next->A=k;
  //    temp->next->a=ATOM->contigous[k][TYPE];
  //    ////temp->next->Ax=midbx;
  //    ////temp->next->Ay=midby;
  //}
  //temp->next->B=A;
  //temp->next->C=B;
  //temp->next->b=ATOM->contigous[A][TYPE];
  //temp->next->c=ATOM->contigous[B][TYPE];
  //temp->next->Bx=D->Ax;
  //temp->next->By=D->Ay;
  //temp->next->circum_x=circx;
  //temp->next->circum_y=circy;
  //temp->next->circum_z=circz;
  //if(flag)
  //    ATOM->conti[TYPE]++;
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
    for(int i=0; i<ATOM->conti[TYPE]; i++)
    //for(int i=0; i<1; i++)
    {
        int s=0;
        //cout<<i<<"\n";
        for(int p=0; p<ATOM->conti[TYPE]; p++)
        {
            s=s+ATOM->part_c[i][p][TYPE];
        }
        //cout<<s<<"\t"<<2*ATOM->edge_index[i][TYPE]<<"\n";
      //if(s==2*ATOM->edge_index[i][TYPE])
      //{
      //    //		cout<<"somehting BDSBSDBSHBD\n";
      //    //	cout<<"here\n";
      //    break;
      //}

        for(int j=0; j<ATOM->conti[TYPE]; j++)
        {
            if(i!=j)
            {
				if(ATOM->part_c[i][j][TYPE]==1)
				{
        			for(int k=0; k<ATOM->conti[TYPE]; k++)
					{
						if(k!=i && k!=j)
						{
							container_delunay *temp;
							temp=ATOM->D_FIRST;
							//int flag=1;
							//cout<<i<<"\t"<<j<<"\t"<<k<<"\n";
							while(1)
							{
								std::vector<int> EV1 {ATOM->index,ATOM->contigous[i][TYPE],ATOM->contigous[j][TYPE],ATOM->contigous[k][TYPE]};//,A3,A4;
								std::sort(EV1.begin(),EV1.end());
								if(temp->D->AT[0]==EV1[0] && temp->D->AT[1]==EV1[1] && temp->D->AT[2]==EV1[2] && temp->D->AT[3]==EV1[3])
								{
								  //cout<<"existing delunay \t";
								  //cout<<temp->D->AT[0]<<"\t"<<temp->D->AT[1]<<"\t"<<temp->D->AT[2]<<"\t"<<temp->D->AT[3]<<"\n";
								  //cout<<EV1[0]<<"\t"<<EV1[1]<<"\t"<<EV1[2]<<"\t"<<EV1[3]<<"\n";
								  //cout<<ATOM->index<<"\t"<<ATOM->contigous[i][TYPE]<<"\t"<<ATOM->contigous[j][TYPE]<<"\t";
								    a=ATOM->x;
								    b=ATOM->y;
								    c=ATOM->z;
								    p=Atoms[ATOM->contigous[i][TYPE]].x-a;
								    q=Atoms[ATOM->contigous[i][TYPE]].y-b;
								    r=Atoms[ATOM->contigous[i][TYPE]].z-c;
								    rA=Atoms[ATOM->contigous[i][TYPE]].radius+r_cut;
								    Sx=Atoms[ATOM->contigous[j][TYPE]].x-a;
								    Sy=Atoms[ATOM->contigous[j][TYPE]].y-b;
								    Sz=Atoms[ATOM->contigous[j][TYPE]].z-c;
								    rB=Atoms[ATOM->contigous[j][TYPE]].radius+r_cut;
								    Px=Atoms[ATOM->contigous[k][TYPE]].x-a;
								    Py=Atoms[ATOM->contigous[k][TYPE]].y-b;
								    Pz=Atoms[ATOM->contigous[k][TYPE]].z-c;
								    //cout<<ATOM->contigous[k][TYPE]<<"\n";
								    //	cout<<"otherside\n";
								    //	cout<<"s\n";
								    rC=Atoms[ATOM->contigous[k][TYPE]].radius+r_cut;
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
								    delunay *D=nullptr;
								    if(ATOM->part_c[i][j][TYPE]==1)
								    	//if(iBf==1)
								    {
								    	long double ax=(Sy*r-Sz*q);
								    	long double ay=(Sz*p-Sx*r);
								    	long double az=(Sx*q-Sy*p);
								    	long double overlap=(Px*ax+Py*ay+Pz*az);
								    	if(overlap<0.)
								    		sign=1;
								    	else
								    		sign=-1;
								    	//	cout<<"i\t"<<"1"<<"\n";
										//cout<<"new del\n";
								    	constr_del(ATOM,Atoms,TYPE,p,q,r,Sx,Sy,Sz,sign,i,j,k,rA,rB,D);
								    	//cout<<ATOM->part_c[i][][TYPE]<<"\n";
								    ////for(int p=0; p<ATOM->conti[TYPE]; p++)
								    ////{
								    ////	cout<<p<<"\t"<<ATOM->part_c[i][p][TYPE]<<",";
								    ////}
								    	//cout<<"\n";
								    }
									//cout<<"here\n";
								    break;

								}
								else if (temp->next)
								{
									temp=temp->next;
								}
								else
									break;	
							}
							//cout<<i<<"\t"<<j<<"\n";
							if(ATOM->part_c[i][j][TYPE]==2)
							{
								break;
							}
							//cout<<ATOM->index<<"\t"<<i<<"\t"<<j<<"\t"<<k<<"\n";
						}
					}
				}
            }
        }
		s=0;
        for(int p=0; p<ATOM->conti[TYPE]; p++)
        {
            s=s+ATOM->part_c[i][p][TYPE];
			//cout<<p<<"\t"<<ATOM->part_c[i][p][TYPE]<<"\n";
        }
        //cout<<s<<"\t"<<2*ATOM->edge_index[i][TYPE]<<"\n";
    }
}

void vert_void(vertice *v,atom Atoms[], int nAtoms,int TYPE,int void_vert_count)
{
    int flag=1;
    long double AX,AY,AZ,BX,BY,BZ;
    BX=v->p->x ;
    BY=v->p->y;
    BZ=v->p->z;
    AX=Atoms[v->A].x;
    AY=Atoms[v->A].y;
    AZ=Atoms[v->A].z;
    long double dis=sqrtl((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY)+(AZ-BZ)*(AZ-BZ));
    if(dis<r_cut+Atoms[v->A].radius)
    {
        flag=0;
    }
    AX=Atoms[Atoms[v->A].contigous[v->D->A][TYPE]].x;
    AY=Atoms[Atoms[v->A].contigous[v->D->A][TYPE]].y;
    AZ=Atoms[Atoms[v->A].contigous[v->D->A][TYPE]].z;
    dis=sqrtl((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY)+(AZ-BZ)*(AZ-BZ));
    if(dis<r_cut+Atoms[Atoms[v->A].contigous[v->D->A][TYPE]].radius)
    {
        flag=0;
    }
    AX=Atoms[Atoms[v->A].contigous[v->D->B][TYPE]].x;
    AY=Atoms[Atoms[v->A].contigous[v->D->B][TYPE]].y;
    AZ=Atoms[Atoms[v->A].contigous[v->D->B][TYPE]].z;
    dis=sqrtl((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY)+(AZ-BZ)*(AZ-BZ));
    if(dis<r_cut+Atoms[Atoms[v->A].contigous[v->D->B][TYPE]].radius)
    {
        flag=0;
    }
    AX=Atoms[Atoms[v->A].contigous[v->D->C][TYPE]].x;
    AY=Atoms[Atoms[v->A].contigous[v->D->C][TYPE]].y;
    AZ=Atoms[Atoms[v->A].contigous[v->D->C][TYPE]].z;
    dis=sqrtl((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY)+(AZ-BZ)*(AZ-BZ));
    if(dis<r_cut+Atoms[Atoms[v->A].contigous[v->D->C][TYPE]].radius)
    {
        flag=0;
    }
    //flag=1 if the vertice is not in the exclusion area of the atoms
    if(flag )
    {
        v->is_void=1;
        //mark it as void
        void_vert_count=void_vert_count+1;
    }
}
long double volume_tetrahedron(long double Ax,long double Ay,long double Az,long double Ex,long double Ey,long double Ez,long double Bx,long double By,long double Bz,long double vx,long double vy,long double vz,long double r)
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
    ABx=(ABx-(tilt*lroundl(ABy/twob)));
    ABx=(ABx-(twob*lroundl(ABx/twob)));
    ABy=(ABy-(twob*lroundl(ABy/twob)));
    ABz=(ABz-(twob*lroundl(ABz/twob)));
    EBx=(EBx-(tilt*lroundl(EBy/twob)));
    EBx=(EBx-(twob*lroundl(EBx/twob)));
    EBy=(EBy-(twob*lroundl(EBy/twob)));
    EBz=(EBz-(twob*lroundl(EBz/twob)));
    VEx=(VEx-(tilt*lroundl(VEy/twob)));
    VEx=(VEx-(twob*lroundl(VEx/twob)));
    VEy=(VEy-(twob*lroundl(VEy/twob)));
    VEz=(VEz-(twob*lroundl(VEz/twob)));
    //cout<<ABx+Bx<<"\t"<<ABy+By<<"\t"<<ABz+Bz<<"}\tradius\t"<<r<<"\tresolution\t"<<500<<"\n";
    //cout<<"mol new\n";
    //cout<<"draw material Opaque\n";
    //cout<<"draw color red\n";
    //cout<<"draw line\t{";
    //cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\n";
    //cout<<"draw line\t{";
    //cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<vx<<"\t"<<vy<<"\t"<<vz<<"}\n";
    //cout<<"draw line\t{";
    //cout<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\t{"<<vx<<"\t"<<vy<<"\t"<<vz<<"}\n";
    //cout<<"draw line\t{";
    //cout<<vx<<"\t"<<vy<<"\t"<<vz<<"}\t{"<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\n";
    //cout<<"draw line\t{";
    //cout<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\t{"<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\n";
////long double a,b,c;
////long double a1,b1,c1;
////a=(ABy*EBz-ABz*EBy);
////b=(ABz*EBx-ABx*EBz);
////c=(ABx*EBy-ABy*EBx);
////a1=(VEy*c-VEz*b);
////b1=(VEz*a-VEx*c);
////c1=(VEx*b-VEy*a);
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
    //cout<<"mol new\n";
    //cout<<"draw material Opaque\n";
    //cout<<"draw sphere\t{";
    //cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\tradius\t"<<r<<"\tresolution\t"<<500<<"\n";//{"<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\n";
    //cout<<"mol new\n";
    //cout<<"draw material Opaque\n";
    //cout<<"draw line\t{";
    //cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\twidth 5\n";
    //cout<<"draw line\t{";
    //cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\twidth 5\n";
    //cout<<"draw line\t{";
    //cout<<Ax<<"\t"<<Ay<<"\t"<<Az<<"}\t{"<<vx<<"\t"<<vy<<"\t"<<vz<<"}\twidth 5\n";
    //cout<<"draw line\t{";
    //cout<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\t{"<<vx<<"\t"<<vy<<"\t"<<vz<<"}\twidth 5\n";
    //cout<<"draw line\t{";
    //cout<<vx<<"\t"<<vy<<"\t"<<vz<<"}\t{"<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\twidth 5\n";
    //cout<<"draw line\t{";
    //cout<<Ex<<"\t"<<Ey<<"\t"<<Ez<<"}\t{"<<Bx<<"\t"<<By<<"\t"<<Bz<<"}\twidth 5\n";
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
////cout<<"x0="<<x0<<"\t"<<y0<<"\t"<<z0<<"\n";
////cout<<"ra="<<r<<"\t"<<rB<<"\t"<<rE<<"\t"<<rV<<"\n";
    long double theta,x2,y2;
    x2=r*x0/rE;
    y2=r*y0/rE;
    theta=atanl(z0/y0);
    //cout<<theta<<"\n";
    //cout<<"#\t";
    long double Vc,Vt;
    Vt=(x0*y0*z0)/6.;
    //cout<<(x0*y0*z0)/6.<<" vol\n";
    if(r<rB)
    {
        Vc=(r*r*r/6.)*(2*theta-M_PI/2.-asinl((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq))));
        //return (x0*y0*z0)/6.-(r*r*r/6.)*(2*theta-M_PI/2.-asin((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq))));
    }
    if(rB<r && r<rE)
    {
        Vc=theta/2.*(r*r*x0-x0*x0*x0/3.)-(r*r*r/6.)*(M_PI/2.+asinl((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq))));
        //return (x0*y0*z0)/6.-theta/2.*(r*r*x0-x0*x0*x0/3.)-(r*r*r/6.)*(M_PI/2.+asin((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq))));
    }
    if(rE<r && r<rV)
    {
        //cout<<"hea\n";
        Vc=0.5*(theta-M_PI/2.+asin(y0/sqrtl(r*r-x0sq)))*(r*r*x0-x0*x0*x0/3.)+x0*y0/6.*sqrtl(r*r-rE*rE)+r*r*r/6.*asinl((x2*x2-y2*y2-x0sq)/(r*r-x0sq))-r*r*r/6.*asinl((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq)));
        //return (x0*y0*z0)/6.-0.5*(theta-M_PI/2.+asin(y0/sqrtl(r*r-x0sq)))*(r*r*x0-x0*x0*x0/3.)+x0*y0/6.*sqrtl(r*r-rE*rE)+r*r*r/6.*((x2*x2-y2*y2-x0sq)/(r*r-x0sq))-r*r*r/6.*asin((z0sq*x0sq-y0sq*rV*rV)/(rE*rE*(y0sq+z0sq)));
    }
    //cout<<Vt<<"\t"<<Vc<<"\t"<<Vt-Vc<<"\n";
    if(Vt<Vc)
    {
        cout<<"help \t"<<Vt<<"\t"<<Vc<<"\t"<<Vt-Vc<<"\n";
    }
    return Vt-Vc;
    //cout<<a<<"\t"<<b<<"\t"<<c<<"\n";

    //cout<<ABx*EBx+ABy*EBy+ABz*EBz<<"\n";
////cout<<VEx*EBx+VEy*EBy+VEz*EBz<<"\n";
////cout<<VEx*ABx+VEy*ABy+VEz*ABz<<"\n";
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
    a1x=(a1x-(tilt*lroundl(a1y/twob)));
    a1x=(a1x-(twob*lroundl(a1x/twob)));
    a1y=(a1y-(twob*lroundl(a1y/twob)));
    a1z=(a1z-(twob*lroundl(a1z/twob)));
    a2x=(a2x-(tilt*lroundl(a2y/twob)));
    a2x=(a2x-(twob*lroundl(a2x/twob)));
    a2y=(a2y-(twob*lroundl(a2y/twob)));
    a2z=(a2z-(twob*lroundl(a2z/twob)));
    ex=(ex-(tilt*lroundl(ey/twob)));
    ex=(ex-(twob*lroundl(ex/twob)));
    ey=(ey-(twob*lroundl(ey/twob)));
    ez=(ez-(twob*lroundl(ez/twob)));
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

int main( int argc, char * argv[] )
{
    int nAtoms=0;
    atom *Atoms=NULL;
    int counter=0;
    int config_count=0;
    //The file with configurations
    std::ifstream infile(argv[1]);///dat_trial");//config_2000_0.38_2_0.70.dat");
    //No of Atoms
    nAtoms=216;
    cout<<std::setprecision(5);
    //No of configurations in the input file
    config_count=1;
    //No of types of particle
    int ntypes=2;
    int SAM=0;
	//delunay *****DSET;
////DSET=new (nothrow) delunay***** [nAtoms];
////for(int i=0;i<nAtoms;i++)
////	DSET[i]=new (nothrow) delunay**** [nAtoms];
////for(int i=0;i<nAtoms;i++)
////	for(int j=0;j<nAtoms;j++)
////		DSET[i][j]=new (nothrow) delunay*** [nAtoms];
////for(int i=0;i<nAtoms;i++)
////	for(int j=0;j<nAtoms;j++)
////		for(int k=0;k<nAtoms;k++)
////		{
////			cout<<i<<"\t"<<j<<"\t"<<k<<"\n";
////			DSET[i][j][k]=new (nothrow) delunay** [nAtoms];
////		}
//  for(int i=0;i<nAtoms;i++)
//  	for(int j=0;j<nAtoms;j++)
//  		for(int k=0;k<nAtoms;k++)
//  			for(int l=0;l<nAtoms;l++)
//  				DSET[i][j][k][l]=new (nothrow) delunay* [ntypes];
    long double *radius=new (nothrow) long double[ntypes];
    long double b,c,d,e,f;
    //radiuses of the particle
    radius[0]=0.5;
    radius[1]=0.7;
    if(3*2*radius[1]<box)
    {
        LAYER_CUT=3*2*radius[1]-box;
    }
    else
        LAYER_CUT=0.;
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
        //cout<<nconfig<<"\n"<<std::flush;
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
            //while(infile>>b>>c>>d>>e)
            {
                Atoms[counter].x=b;
                Atoms[counter].y=c;
                Atoms[counter].z=d;
                Atoms[counter].radius=e-epsilon;
                Atoms[counter].index=counter;
                //cout<<nAtoms<<"\t"<<b<<"\t"<<c<<"\t"<<d<<"\t"<<e<<"\n";;
                //nAtoms++;
            }
            /* The loops puts the first nAtoms lines in the input file to list in the descending order*/
            counter++;
            //if(sites==NULL)
            //{
            //    sites=new vertice;
            //    sites->p=new site;
            //    sites->p->x=b;
            //    sites->p->y=c;
            //    sites->p->z=d;
            //    sites->r=e-epsilon;
            //    //display_SITE(sites->p);
            //}
            //else
            //{
            //    temp_site=new vertice;
            //    temp_site->p=new site;
            //    temp_site->p->x=b;
            //    temp_site->p->y=c;
            //    temp_site->p->z=d;
            //    temp_site->r=e-epsilon;
            //    insert_site(sites,temp_site);
            //}
            if(counter==nAtoms)
            {
                break;
            }
        }
        //cout<<"brea\n";
        //temp_site=sites;
        //int cunt=0;
        //while(1)
        //{
        //    /* this loop puts the list of atoms into an array*/
        //    Atoms[cunt].x=temp_site->p->x;
        //    Atoms[cunt].y=temp_site->p->y;
        //    Atoms[cunt].z=temp_site->p->z;
        //    Atoms[cunt].radius=temp_site->r;
        //    cunt++;
        //    if(temp_site->next)
        //        temp_site=temp_site->next;
        //    else
        //        break;
        //}
        //cout<<"here\n";
        //Make neighbour list for all atoms
        update_neighbours(Atoms,nAtoms);
        //this loops look for overlaps
        for(int i=0; i<nAtoms; i++)
        {
            //cout<<Atoms[i].x<<"\t"<<Atoms[i].y<<"\t"<<Atoms[i].z<<"\t"<<Atoms[i].radius<<"\n";;

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
                //cout<<Atoms[i].radius<<"\t"<<radius[t]-epsilon<<"\n";
                if(abs(Atoms[i].radius-radius[t]-epsilon)<0.00001)
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
                for(int j=0; j<100; j++)
                {
                    {
                        Atoms[i].part_c[t][j]=new (nothrow) int[ntypes];
                        for(int TYPE=0; TYPE<ntypes; TYPE++)
                        {
                            Atoms[i].part_c[t][j][TYPE]=0;
                        }
                    }
                    if(t<100)
                    {
                        Atoms[i].MIDP[t][j]=new (nothrow) site [ntypes];
                    }
                }
                for(int j=0; j<50 && t<50 ; j++)
                {
                    Atoms[i].D3bondinvoid[t][j]= new (nothrow) int[ntypes];
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
        //for(int TYPE=0; TYPE<ntypes; TYPE++)
        for(int TYPE=0; TYPE<1; TYPE++)
        {
            snprintf(buffer,sizeof(char)*64,"vor_%d",int(TYPE));//_%d_%f.dat",int(nAtoms),Press);
            ofstream vor;
            vor.open(buffer);
            r_cut=radius[TYPE];
            //cout<<TYPE<<"\t"<<r_cut<<"\n";
            //This is the loop over all atoms: we construct the voronoi cell for each atom
            int void_vert_count=0;
            //for(SAM=0 ; SAM<1; SAM++)
            for(SAM=0 ; SAM<nAtoms; SAM++)
            {
					cout<<"##\t"<<SAM<<"\n";

                //   cout<<"mol new\n";
                //   cout<<"draw material Opaque\n";
                // ////////cout<<"here\n";
                // ////////cout<<Atoms[SAM].type<<"\n";
                // if(Atoms[SAM].type==1)
                // {
                //     ////cout<<"mol new\t";
                //     ////cout<<"draw material Transparent\t";
                //       cout<<"draw color blue\n";
                //       cout<<"draw sphere \t";
                //       cout<<"{\t"<<Atoms[SAM].x<<"\t"<<Atoms[SAM].y<<"\t"<<Atoms[SAM].z<<"}\t"<<"radius\t"<<Atoms[SAM].radius<<"\t"<<"resolution\t100\n";
                // }
                // if(Atoms[SAM].type==0)
                // {
                //       //cout<<"mol new\t";
                //       //cout<<"draw material Transparent\t";
                //       cout<<"draw color red\n";
                //       cout<<"draw sphere \t";
                //       cout<<"{\t"<<Atoms[SAM].x<<"\t"<<Atoms[SAM].y<<"\t"<<Atoms[SAM].z<<"}\t"<<"radius\t"<<Atoms[SAM].radius<<"\t"<<"resolution\t100\n";
                // }
                //cout<<"mol new\n";
                //cout<<"draw material Transparent\n";
                ////////cout<<"here\n";
                ////////cout<<Atoms[SAM].type<<"\n";
                //if(Atoms[SAM].type==1)
                //{
                //    ////cout<<"mol new\t";
                //    ////cout<<"draw material Transparent\t";
                //    cout<<"draw color blue\n";
                //    cout<<"draw sphere \t";
                //    cout<<"{\t"<<Atoms[SAM].x<<"\t"<<Atoms[SAM].y<<"\t"<<Atoms[SAM].z<<"}\t"<<"radius\t"<<Atoms[SAM].radius+r_cut<<"\t"<<"resolution\t100\n";
                //}
                //if(Atoms[SAM].type==0)
                //{
                //    ////cout<<"mol new\t";
                //    ////cout<<"draw material Transparent\t";
                //    cout<<"draw color red\n";
                //    cout<<"draw sphere \t";
                //    cout<<"{\t"<<Atoms[SAM].x<<"\t"<<Atoms[SAM].y<<"\t"<<Atoms[SAM].z<<"}\t"<<"radius\t"<<Atoms[SAM].radius+r_cut<<"\t"<<"resolution\t100\n";
                //}
                //cout<<"mol new\n";
                //cout<<"draw material Opaque\n";
                {
                    //the first function calculates the first delunay triangle for the atom
					//cout<<SAM<<"\n";
					if(!Atoms[SAM].D_FIRST)
					{
						//cout<<"here\n";
                    	first_delunay(&(Atoms[SAM]),Atoms,TYPE);
					}
                    //this completes the delunay triangles of the atoms: all the triangle the atom takes part in
                  	delunay *D=nullptr;
				////for(int i=0;i<Atoms[SAM].conti[TYPE];i++)
				////{
				////	cout<<i<<"\t"<<Atoms[SAM].contigous[i][TYPE]<<"\n";
				////}
                  //SAM=60;
                  //D=Atoms[SAM].D_FIRST->D;
                  //cout<<Atoms[SAM].conti[TYPE]<<"\n";
                  //print_delunay(&(Atoms[SAM]),D,Atoms,TYPE);
                  //SAM=134;
                  //cout<<Atoms[SAM].conti[TYPE]<<"\n";
                  //D=Atoms[SAM].D_FIRST->D;
                  //print_delunay(&(Atoms[SAM]),D,Atoms,TYPE);
                  //SAM=164;
                  //cout<<Atoms[SAM].conti[TYPE]<<"\n";
                  //D=Atoms[SAM].D_FIRST->D;
                  //print_delunay(&(Atoms[SAM]),D,Atoms,TYPE);
                    complete_del_2(&(Atoms[SAM]),Atoms,nAtoms,TYPE);
				////D=Atoms[SAM].D_FIRST->D;
				////print_delunay(&(Atoms[SAM]),D,Atoms,TYPE);
				    D=Atoms[SAM].D_FIRST->next->D;
					container_delunay *temp;
				    temp=Atoms[SAM].D_FIRST;
                    while(1)
                    {
				    	print_delunay(&(Atoms[SAM]),temp->D,Atoms,TYPE);
                    	if(temp->next)
                    			temp=temp->next;
                    	else
                    			break;
                    }
			   	}
            }
		}
	}
}
