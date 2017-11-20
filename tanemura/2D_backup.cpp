#include<iostream>
#include<math.h>
#include<fstream>
#include<iomanip>
#include <set>
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
    //vert_list *V;
    vertice *neib_vert[10];
    int neib_ed[10];
    int v_neigh_count=0;
}**start,*s_temp;
struct container_vertice
{
    struct vertice *V;
    struct container_vertice *next=NULL;
    struct container_vertice *prev=NULL;
}**CSTART,*VOID_START;
void add_connected(vertice *focus,vertice *add,int binv)
{
    int flag=1;
    int i;
    for(i=0; i<focus->v_neigh_count; i++)
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
    flag=compare(EV->p,v->p);
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
        //display_SITE(v->p);
        for(int n=0; n<v->v_neigh_count; n++)
        {
            if(v->neib_vert[n])
                for(int m=0; m<v->neib_vert[n]->v_neigh_count; m++)
                {
                    if(v->neib_vert[n]->neib_vert[m])
                        if(!(compare(v->neib_vert[n]->neib_vert[m]->p,v->p)))
                        {
                            v->neib_vert[n]->neib_vert[m]=NULL;
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
    Sx=ATOM->x-Atoms[ATOM->contigous[D->A][TYPE]].x;
    Sy=ATOM->y-Atoms[ATOM->contigous[D->A][TYPE]].y;
    Sx=(Sx-(twob*lround(Sx/twob)));
    Sy=(Sy-(twob*lround(Sy/twob)));
    cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\n";
    Sx=ATOM->x-Atoms[ATOM->contigous[D->B][TYPE]].x;
    Sy=ATOM->y-Atoms[ATOM->contigous[D->B][TYPE]].y;
    Sx=(Sx-(twob*lround(Sx/twob)));
    Sy=(Sy-(twob*lround(Sy/twob)));
    cout<<ATOM->x<<"\t"<<ATOM->y<<"\t"<<ATOM->x-Sx<<"\t"<<ATOM->y-Sy<<"\n";
    Sx=ATOM->x-Atoms[ATOM->contigous[D->B][TYPE]].x;
    Sy=ATOM->y-Atoms[ATOM->contigous[D->B][TYPE]].y;
    Sx=(Sx-(twob*lround(Sx/twob)));
    Sy=(Sy-(twob*lround(Sy/twob)));
    Px=ATOM->x-Atoms[ATOM->contigous[D->A][TYPE]].x;
    Py=ATOM->y-Atoms[ATOM->contigous[D->A][TYPE]].y;
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
    std::ifstream infile("dat");
    infile>>nAtoms;
    infile>>box;
    int ntypes=2;
    twob=2*box;
    density=nAtoms/(twob*twob);
    int SAM=0;
    double *radius=new (nothrow) double[ntypes];
    start = new (nothrow) vertice*[ntypes];
    CSTART = new (nothrow) container_vertice*[ntypes];
    Atoms = new (nothrow) atom[nAtoms];
    long double b,c,d,e;
    for(int t=0; t<ntypes; t++)
    {
        infile>>radius[t];
    }
    nAtoms=0;
    char buffer[64];
    while(infile>>b>>c>>d>>e)
    {
        Atoms[nAtoms].x=b;
        Atoms[nAtoms].y=c;
        Atoms[nAtoms].radius=d;
        //cout<<nAtoms<<"\t"<<b<<"\t"<<c<<"\n";
        nAtoms++;
    }
    update_neighbours(Atoms,nAtoms);
    double area=0;
    vertice *save;
    //cout<<"here\n";
    for(int i=0; i<nAtoms; i++)
    {
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
        for(SAM=0 ; SAM<nAtoms ; SAM++)
        {
            {
                first_delunay(&(Atoms[SAM]),Atoms,TYPE);
                complete_del(&(Atoms[SAM]),Atoms,nAtoms,TYPE);
                delunay *D;
                int count=0;
                double area_s=0;
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
                    if(SAM!=56  && Atoms[SAM].contigous[i][TYPE] !=56)
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
	cout<<"after  first tessellation \t"<<TYPE<<"\n";;
        while(1)
        {
	    if(temp_start->v_neigh_count == 2)
	    {
	    	V->delete_vertice(start[TYPE],temp_start,TYPE);
	    	temp_start=temp_start->next;
	    	continue;	
	    }
	    if(TYPE==0)
	    {
	    display_SITE(temp_start->p);
	    }
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
        vor.close();
    }
    //return 0;
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
    //SAM=2;
    int TYPE=1;
  //cout<<Atoms[35].conti[TYPE]<<"\n";
  //for(int k=0; k<Atoms[35].conti[TYPE]; k++)
  //{
  //    cout<<k<<"\t"<<Atoms[35].contigous[k][TYPE]<<"\n";
  //}
    for(int i=0; i<nAtoms; i++)
    {
        container_vertice *new_vert;
        //set_of_delunay D;
        //cout<<Atoms[i].x<<"\t"<<Atoms[i].y<<"\n";
        r_cut=Atoms[i].radius;
        cout<<i<<" Atoms NO \n";
        for(int t=0; t<ntypes; t++)
        {
            if(r_cut==radius[t])
            {
                TYPE=t;
                break;
            }
        }
        cout<<TYPE<<"\t"<<"TYPE"<<"\n";
        container_vertice *ctemp;
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
        //cout<<TYPE<<" TYPE\n";
        for(int j=0; j<Atoms[i].conti[TYPE]; j++)
        {
            //cout<<"eye "<<j<<"\t"<<Atoms[35].conti[TYPE]<<"\n";
            int flag=0;
            int SAM=Atoms[i].contigous[j][TYPE];
            container_vertice *ctemp;
            ctemp=Atoms[SAM].Cstart[TYPE];
            save_atom(&(Atoms[SAM]),TYPE);
            delunay *D;
            D=Atoms[SAM].D[TYPE].initial;
            delunay *temp;
            temp=Atoms[SAM].D[TYPE].initial;
            Atoms[SAM].save_D[TYPE].initial=temp;
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
                        container_vertice *ctemp;
                        ctemp = new container_vertice;
                        ctemp->V=temp;
                        if(!CSTART[TYPE])
                            CSTART[TYPE]=ctemp;
                        else
                            insert_cvertice(CSTART[TYPE],ctemp,CSTART[TYPE]);
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
                    container_vertice *ctemp;
                    ctemp = new container_vertice;
                    ctemp->V=temp;
                    if(!CSTART[TYPE])
                    {
                        //cout<<"here>\n";
                        CSTART[TYPE]=ctemp;
                    }
                    else
                    {
                        //cout<<"here2\n";
                        //display_SITE(ctemp->V->p);
                        insert_cvertice(CSTART[TYPE],ctemp,CSTART[TYPE]);
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
                    vor<<x<<"\t"<<y<<"\t";
                    vor<<D_ONE->circum_x<<"\t"<<D_ONE->circum_y<<"\t"<<D_TWO->circum_x<<"\t"<<D_TWO->circum_y<<"\n";
                    vor<<"\n";
                }
                add_connected(temp_vert_o,temp_vert_d,Atoms[SAM].bondinvoid[i][TYPE]);
                add_connected(temp_vert_d,temp_vert_o,Atoms[SAM].bondinvoid[i][TYPE]);
            }
	    //cout<<"here\n";
        }
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
        vertice *temp_start;
        //cout<<r_cut<<"\n";
        while(1)
        {
            temp_start=cstart->V;
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
            // BX=temp_start->p->x;
            // BY=temp_start->p->y;
            X=AX-BX;
            Y=AY-BY;
            X=(X-(twob*lround(X/twob)));
            Y=(Y-(twob*lround(Y/twob)));
            dis=sqrt((X)*(X)+(Y)*(Y));
            if(flag)//&& temp_start_start->v_neigh_count == 3)
            {
                //display_SITE(temp_start->p);
                temp_start->is_void=1;
                ///temp_start->cluster_index=void_vert_count;
            }
            if(cstart->next)
                cstart=cstart->next;
            else
                break;
        }
        cstart=CSTART[TYPE];
////////if(!new_vert)
////////	new_vert=cstart;
        int void_vert_count=0;
        new_vert=new container_vertice;
        container_vertice *temp_new_vert;
        temp_new_vert=new_vert;
        while(1)
        {
            temp_new_vert->V=cstart->V;
            //print_delunay(&Atoms[cstart->V->A],cstart->V->D,Atoms,TYPE);
            //display_SITE(cstart->V->p);
            if(cstart->next)
            {
                cstart=cstart->next;
                temp_new_vert->next=new container_vertice;
                temp_new_vert=temp_new_vert->next;
            }
            else
                break;
        }
        //cout<<" ... \n";
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
                //display_SITE(cstart->V->p);
                void_vert_count++;
            }
            if(cstart->next)
                cstart=cstart->next;
            else
                break;
        }
        cout<<void_vert_count<<"\n";
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
            //print_delunay(&Atoms[cavity_list[j]->A],cavity_list[j]->D,Atoms);
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
        cout<<cav_tot<<"\t"<<ca_per_tot<<"\n";
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
        //container_vertice *ctemp;
      //for(int j=0; j<Atoms[i].conti[TYPE]; j++)
      //{
      //    int SAM=Atoms[i].contigous[j][TYPE];
        ctemp=Atoms[i].Cstart[TYPE];
        while(1)
        {
            //display_SITE(ctemp->V->p);
            V->insert_vertice(start[TYPE],ctemp->V,TYPE);
            //V->delete_vertice(start[TYPE],ctemp->V,TYPE);
            if(ctemp->next)
                ctemp=ctemp->next;
            else
                break;
        }
	//cout<<"after insertion\n";
////////temp_start=start[TYPE];
////////while(1)
////////{
////////	display_SITE(temp_start->p);
////////	if(temp_start->next)
////////	{
////////		temp_start=temp_start->next;
////////	}
////////	else
////////		break;
////////}
      //}
        for(int j=0; j<Atoms[i].conti[TYPE]; j++)
        {
            int flag=0;
            int SAM=Atoms[i].contigous[j][TYPE];
            reset_atom(&(Atoms[SAM]),TYPE);
        }
        delete[] void_area;
        delete[] void_length;
        delete[] cavity_list;
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
        //for(int j=0; j<Atoms[8].conti[TYPE]; j++)
        //{
        //    delunay *D;
        //    SAM=Atoms[8].contigous[j][TYPE];
        //    D=Atoms[SAM].D[1].initial;
        //    cout<<"#"<<SAM<<"\n";
        //    while(1)
        //    {
        //        print_delunay(&(Atoms[SAM]),D,Atoms,1);
        //        if(D->next)
        //        {
        //            D=D->next;
        //        }
        //        else
        //            break;
        //    }
        //}
        //SAM=31;
        //delunay *D;
        //D=Atoms[SAM].D[0].initial;
        //cout<<"#"<<SAM<<"\n";
        //while(1)
        //{
        //        print_delunay(&(Atoms[SAM]),D,Atoms,0);
        //        if(D->next)
        //        {
        //    	    D=D->next;
        //        }
        //        else
        //    	    break;
        //}
////////    SAM=31;
////////    for(int j=0;j<Atoms[SAM].conti[0];j++)
////////	    {
////////		    cout<<Atoms[SAM].contigous[j][0]<<"\n";

////////	    }
        vor.close();
        //      if(i==7)
        //      {
        //          for(int j=0; j<Atoms[2].neighbours; j++)
        //          {
        //              cout<<Atoms[2].neighlist[j]<<"\n";

        //          }

        //      }
/////////SAM=2;
/////////for(int j=0;j<Atoms[SAM].conti[1];j++)
/////////   {
/////////           cout<<Atoms[SAM].contigous[j][1]<<"\n";

/////////   }
    }
    //cout<<"ahere\n";

    //return 0;
    for(int t=0; t<ntypes; t++)
    {
	cout<<t<<"\n";    
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
	cout<<"here\n";
      //for(SAM=0; SAM<nAtoms; SAM++)
      //{
      //	//cout<<SAM<<"\n";
      //    cstart=Atoms[SAM].Cstart[t];
      //    while(1)
      //    {
      //        container_vertice *temp_c;
      //        temp_c=cstart;
      //        if(cstart->next)
      //        {
      //            cstart=cstart->next;
      //            delete temp_c;
      //        }
      //        else
      //        {
      //            delete temp_c;
      //            break;
      //        }
      //    }
      //}
    }
    for(int n=0; n<nAtoms; n++)
    {
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
        ////while(1)
        ////{

        ////}
    }
    delete[] Atoms;

    return 0;
}
/*
vertice *temp_start;
for(int i=35; i<36; i++)
{
r_cut=Atoms[i].radius;
container_vertice *ctemp;
    temp_start=start;
    ctemp=Atoms[i].Cstart;
    while(1)
    {
        V->delete_vertice(start,ctemp->V);
        if(ctemp->next)
            ctemp=ctemp->next;
        else
            break;
    }
    //loop thru the atoms to remove it and ressellate to study free volume or other evil stuff
    cout<<Atoms[i].conti<<"\n";
    for(int j=0; j<Atoms[i].conti; j++)
    {
        int flag=0;
        for(int k=0; k<Atoms[Atoms[i].contigous[j]].neighbours; k++)
        {
            if(Atoms[Atoms[i].contigous[j]].neighlist[k]==i)
            {
                flag=1;
            }
            if(flag)
            {
                Atoms[Atoms[i].contigous[j]].neighlist[k]=Atoms[Atoms[i].contigous[j]].neighlist[k+1];
            }
        }
        int SAM=Atoms[i].contigous[j];
        first_delunay(&(Atoms[SAM]),Atoms);
        complete_del(&(Atoms[SAM]),Atoms,nAtoms);
        delunay *D;
        D=Atoms[SAM].D.initial;
        int count=0;
        D=Atoms[SAM].D.initial;
        double area_s=0;
        for(int i=0; i<Atoms[SAM].conti; i++)
        {
            D=Atoms[SAM].D.initial;
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
	    container_vertice *ctemp;
	    ctemp = new container_vertice;
	    ctemp->V=temp;
	    if(!CSTART)
	       	CSTART=ctemp;
	    else
  	        insert_cvertice(CSTART,ctemp,CSTART);
                //display_SITE(temp->p);
                temp_vert_o=temp;
            }
            vertice *temp;
            temp=new vertice;
            temp->p=new site;
            temp->p->x=D_TWO->circum_x;
            temp->p->y=D_TWO->circum_y;
            temp->A=SAM;
            temp->D=D_TWO;
            temp=V->insert_vertice(start,temp);
	container_vertice *ctemp;
	ctemp = new container_vertice;
	ctemp->V=temp;
	if(!CSTART)
	   	CSTART=ctemp;
	else
  	    insert_cvertice(CSTART,ctemp,CSTART);
            temp_vert_d=temp;
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
            if(dis>Atoms[Atoms[SAM].contigous[i]].radius+Atoms[SAM].radius+2.*r_cut)
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
            add_connected(temp_vert_o,temp_vert_d,Atoms[SAM].bondinvoid[i]);
            add_connected(temp_vert_d,temp_vert_o,Atoms[SAM].bondinvoid[i]);
        }
    }
    vertice *temp;
    temp_start=start;
int void_vert_count;
    void_vert_count=0;
container_vertice *cstart;
cstart=CSTART;
    while(1)
    {
	temp_start=cstart->V;
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
            AX=Atoms[temp_start->D->A].x;
            AY=Atoms[temp_start->D->A].y;
            dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
            if(dis<r_cut+Atoms[temp_start->D->A].radius)
            {
            	flag=0;
            	//break;
            }
            AX=Atoms[temp_start->D->B].x;
            AY=Atoms[temp_start->D->B].y;
            dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
            if(dis<r_cut+Atoms[temp_start->D->B].radius)
            {
            	flag=0;
            	//break;
            }
            AX=Atoms[i].x;
            AY=Atoms[i].y;
         // BX=temp_start->p->x;
         // BY=temp_start->p->y;
            X=AX-BX;
            Y=AY-BY;
            X=(X-(twob*lround(X/twob)));
            Y=(Y-(twob*lround(Y/twob)));
            dis=sqrt((X)*(X)+(Y)*(Y));
    	if(flag)//&& temp_start_start->v_neigh_count == 3)
    	{
    		temp_start->is_void=1;
    		temp_start->cluster_index=void_vert_count;
    	}
    	if(cstart->next)
    		cstart=cstart->next;
    	else
    		break;
    }
cstart=CSTART;
//////  temp_start=start;
cout<<"here\n";
    while(1)
    {
	if(cstart->V->is_void)
	display_SITE(cstart->V->p);
    	if(cstart->next)
    	{
    		cstart=cstart->next;
    	}
    	else
    		break;

    }
int change=1;
int void_vert_count_prev;
int tem_ind=0;
while(change)
{
	cstart=CSTART;
	void_vert_count_prev=void_vert_count;
	while(1)
	{
		//cout<<cstart->V->v_neigh_count<<"\n";
		//display_SITE(cstart->V->p);
		//cout<<"this\n";
		if(CSTART->V->is_void)
		{
			for(int i=0;i<cstart->V->v_neigh_count;i++)
			{
				//cout<<"here2\n";
				//if(cstart->V->neib_vert[i]->p->x)
				//cout<<"here\n";
				if(cstart->V->neib_vert[i])
				{
					//display_SITE(cstart->V->neib_vert[i]->p);
					if(cstart->V->neib_vert[i]->is_void )
					{
					////////cout<<"here\n";
					////////display_SITE(cstart->V->neib_vert[i]->p);
					////////cout<<"here2\n";
						int flag;
						container_vertice *ctemp;
						ctemp=new container_vertice;
						ctemp->V=cstart->V->neib_vert[i];
						flag=insert_cvertice(CSTART,ctemp,CSTART);
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
////////tem_ind++;
////////if(tem_ind==3)
////////	break;
        if(void_vert_count_prev==void_vert_count)
        	change=0;
        else
        	change=1;
}
cstart=CSTART;
//////  temp_start=start;
    while(1)
    {
    	//display_SITE(cstart->V->p);
    	if(cstart->next)
    	{
    		cstart=cstart->next;
    	}
    	else
    		break;

    }
}
return 0;
temp_start=start;
double r_cut_sq=0.5*0.5;
int void_vert_count=0;
while(1)
{
    int flag=1;
    double AX,AY,BX,BY,X,Y;
    //cout<<"ver\t"<<temp_start->p->x<<"\t"<<temp_start->p->y<<"\n";
    ////////for(int i=0;i<nAtoms;i++)
    ////////{
    ////////	//double AX,AY,BX,BY;
    ////////	AX=Atoms[i].x;
    ////////	AY=Atoms[i].y;
    ////////	BX=temp_start->p->x ;
    ////////	BY=temp_start->p->y;
    ////////	X=AX-BX;
    ////////	Y=AY-BY;
    ////////	X=(X-(twob*lround(X/twob)));
    ////////	Y=(Y-(twob*lround(Y/twob)));
    ////////	double dis=sqrt((X)*(X)+(Y)*(Y));
    ////////	if(dis<r_cut+Atoms[i].radius)
    ////////	{
    ////////		flag=0;
    ////////		break;
    ////////	}
    ////////
    ////////}
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
    AX=Atoms[temp_start->D->A].x;
    AY=Atoms[temp_start->D->A].y;
    dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
    if(dis<r_cut+Atoms[temp_start->D->A].radius)
    {
        flag=0;
        //break;
    }
    AX=Atoms[temp_start->D->B].x;
    AY=Atoms[temp_start->D->B].y;
    dis=sqrt((AX-BX)*(AX-BX)+(AY-BY)*(AY-BY));
    if(dis<r_cut+Atoms[temp_start->D->B].radius)
    {
        flag=0;
        //break;
    }
    if(flag )//&& temp_start->v_neigh_count == 3)
    {
        temp_start->is_void=1;
        temp_start->cluster_index=void_vert_count;
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
        //cout<<i<<"\t"<<cavity_list[i]<<"\t"<<temp_start<<"\n";;
    }
    if(temp_start->next)
        temp_start=temp_start->next;
    ealse
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
    for(int i=0; i<void_vert_count; i++)
    {
        //cout<<i<<"\t"<<cavity_list[0]<<"\n";
        old_label[i]=cavity_list[i]->cluster_index;
        //cout<<i<<"\t"<<old_label[i]<<"\n";
    }
    //cout<<cavity_list[0]<<"\n";
    for(int i=0; i<void_vert_count; i++)
    {
        //cout<<i<<"\n";
        min=cavity_list[i]->cluster_index;
////////		cout<<i<<"\n";
////////		cout<<cavity_list[i]->v_neigh_count<<"\n";
        for(int n=0; n<cavity_list[i]->v_neigh_count; n++)
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
        for(int n=0; n<cavity_list[i]->v_neigh_count; n++)
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
for(int i=0; i<void_vert_count; i++)
{
    cav<<"#"<<i<<"\n\n";
    for(int j=0; j<void_vert_count; j++)
        if(cavity_list[j]->cluster_index==i)
        {
            //cout<<cavity_list[j]->A<<"\t"<<cavity_list[j]->D->A<<"\t"<<cavity_list[j]->D->B<<"\n";
            cav<<cavity_list[j]->p->x<<"\t"<<cavity_list[j]->p->y<<"\n";
        }

}
//V->display_conn(start);
//CALCULATING THE VOID VOLUME IN A GIVEN CAVITY
double *void_area;
void_area= new (nothrow) double[void_vert_count];
double *void_length;
void_length= new (nothrow) double[void_vert_count];
for(int i=0; i<void_vert_count; i++)
    //for(int i=11;i<12;i++)
{
    void_area[i]=0;
    void_length[i]=0;
    for(int j=0; j<void_vert_count; j++)
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
            void_area[i]=void_area[i]+S12*area_trangle(Vx,Vy,E12x,E12y,A2x,A2y,A2r,A1x,A1y);
            void_area[i]=void_area[i]+S12*area_trangle(Vx,Vy,E12x,E12y,0.,0.,A1r,A1x,A1y);
            void_area[i]=void_area[i]+S13*area_trangle(Vx,Vy,E13x,E13y,A3x,A3y,A3r,A1x,A1y);
            void_area[i]=void_area[i]+S13*area_trangle(Vx,Vy,E13x,E13y,0.,0.,A1r,A1x,A1y);
            void_area[i]=void_area[i]+S23*area_trangle(Vx,Vy,E23x,E23y,A3x,A3y,A3r,A1x,A1y);
            void_area[i]=void_area[i]+S23*area_trangle(Vx,Vy,E23x,E23y,A2x,A2y,A2r,A1x,A1y);
            void_length[i]=void_length[i]+S12*perimeter(Vx,Vy,E12x,E12y,A2x,A2y,A2r);
            void_length[i]=void_length[i]+S12*perimeter(Vx,Vy,E12x,E12y,0.,0.,A1r);
            void_length[i]=void_length[i]+S13*perimeter(Vx,Vy,E13x,E13y,A3x,A3y,A3r);
            void_length[i]=void_length[i]+S13*perimeter(Vx,Vy,E13x,E13y,0.,0.,A1r);
            void_length[i]=void_length[i]+S23*perimeter(Vx,Vy,E23x,E23y,A3x,A3y,A3r);
            void_length[i]=void_length[i]+S23*perimeter(Vx,Vy,E23x,E23y,A2x,A2y,A2r);

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
 `a   }

}
delunay *D;
for(SAM=0; SAM<nAtoms; SAM++)
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
*/
