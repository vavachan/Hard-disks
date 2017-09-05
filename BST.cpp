#include<iostream>
#include<math.h>
#include<fstream>
#include <iomanip>
using namespace std;
//Site is a point on the 2D surface
struct face;
struct site
{
	float x=0;
	float y=0;
	struct face *F;
};
void display_SITE(struct site *p)
{
	cout<<"("<<p->x<<" "<<p->y<<")\n";
}
//definition of half-edge 
struct half_edge 
{
	struct half_edge *twin;
	struct half_edge *next;
	struct vertice *origin;
	struct face *facing;
};
//definition of vertice
struct vertice
{
	float x=0;
	float y=0;
	struct half_edge *leaving;
};
//definition of face
struct face
{
	struct site *p;
	struct half_edge *edge;
};
//CLASS TO IMPLEMENT THE DCEL
class DCEL
{
	public:
		struct vertice *origin=NULL;
		struct vertice create_vertex(float , float);
		void create_edge(struct vertice *,struct vertice *);
		
}*vor;


//Breakpoint definition
struct break_point
{
	struct site *p1;
	struct site *p2;
	struct vertice *V=NULL;
};



//node definition 
struct node 
{
        struct site *p;
        struct break_point *B;
        struct node *left; 
        struct node *right;
        struct node *adj_left; 
        struct node *adj_right;
	struct node *parent;
	struct node *twin;
	struct node *L_breakpoint;
	struct node *R_breakpoint;
	struct event *circle_event;
	int isBP=0;
}*root;
//definition of an event 
struct event
{
	struct site *p;
	struct site *center;
	int circle_event=0;
	struct node *circle_node;
	struct event *next=NULL;
	struct event *prev=NULL;
}*start;
void display_events(struct event *EV)
{
	if(EV->circle_event)
		cout<<"circle event= ";
	else
	{
		//display_SITE(EV->p);
	}
	display_SITE(EV->p);
	if(EV->next!=NULL)
	{
		display_events(EV->next);
	}
	else
		cout<<"events over\n";
}
/*
 * Defining a Class called Binary Search Tree
 */
class BST
{
	public:
		void insert(node *, node *);
		void delete_node(node *,site *);
		void display(node *,float,int);
		void find(node *,int);
		void disBeach(node *);
}bst;
//CLASS PRIORITY LIST
class priority_list
{
	public:
		void insert_event(event *,event *,node *);
		void read_events(event *);
		void delete_event(event *,event *);
}P;
/*
 * compare between two sites:
 * Two sites are compared using a convention adopted.
 */
int compare(struct site *p1,struct site *p2)
{
	if( p1->y > p2->y )
		return 1;
	else if ( (p1->y == p2->y) && ( p1->x > p2->x ) )
		return 1;
	else if ( (p1->y == p2->y) && ( p1->x == p2->x ) )
		return 0;
	else 
		return -1;
}

/* cal_break_point calculates
 * the break point defined by the variable break_point with respect to y
 */

struct site cal_break_point(struct break_point *B,float y)
{
	float h=B->p1->y-y;
	float h1=B->p2->y-y;
	float x1=B->p2->x-B->p1->x;
	int flag=compare(B->p1,B->p2);
	struct site BP;
      //cout<<B.p1.x<<"\t"<<B.p1.y<<"\n";
      //cout<<B.p2.x<<"\t"<<B.p2.y<<"\n";
      //cout<<x1<<"\t"<<h<<"\t"<<h1<<"\n";
////////if( abs(h1-h) < 0.01)
////////{
////////	cout<<abs(h-h1)<<"here?\n";
////////}
	if(abs(h1-h)==0.)
	{
		BP.x=x1/2.+B->p1->x;	
		if(B->p1->x<B->p2->x)
		{
			BP.y=(pow((BP.x-B->p1->x),2)+pow((B->p1->y),2)-pow(y,2))/(2.*(B->p1->y-y));
		}
		else
		{
			BP.x=-1.*pow(10,10);
			BP.y=pow(10,10);
		}
		return BP;

	}



	if(flag==1)
	{
		BP.x=(-2.*h*x1+sqrt(pow(2.*h*x1,2)-4.*(h1-h)*((h-h1)*h*h1-h*pow(x1,2))))/(2.*(h1-h))+B->p1->x;
		if(BP.x<B->p2->x)
		{
			BP.y=(pow((BP.x-B->p1->x),2)+pow((B->p1->y),2)-pow(y,2))/(2.*(B->p1->y-y));
			return BP;
		}
		else 
		{
			BP.x=(-2.*h*x1-sqrt(pow(2.*h*x1,2)-4.*(h1-h)*((h-h1)*h*h1-h*pow(x1,2))))/(2.*(h1-h))+B->p1->x;
			BP.y=(pow((BP.x-B->p1->x),2)+pow((B->p1->y),2)-pow(y,2))/(2.*(B->p1->y-y));
			return BP;
		}
	}
	else if (flag == -1)
	{
		BP.x=(-2.*h*x1+sqrt(pow(2.*h*x1,2)-4.*(h1-h)*((h-h1)*h*h1-h*pow(x1,2))))/(2.*(h1-h))+B->p1->x;
		if(BP.x>B->p1->x)
		{
			BP.y=(pow((BP.x-B->p1->x),2)+pow((B->p1->y),2)-pow(y,2))/(2.*(B->p1->y-y));
			return BP;
		}
		else 
		{
			BP.x=(-2.*h*x1-sqrt(pow(2.*h*x1,2)-4.*(h1-h)*((h-h1)*h*h1-h*pow(x1,2))))/(2.*(h1-h))+B->p1->x;
			BP.y=(pow((BP.x-B->p1->x),2)+pow((B->p1->y),2)-pow(y,2))/(2.*(B->p1->y-y));
			return BP;
		}
	}
}
void display_BP(struct break_point *B)
{
	cout<<"("<<B->p1->x<<" "<<B->p1->y<<")+";
	cout<<"("<<B->p2->x<<" "<<B->p2->y<<")\n";
}
void create_edge(struct node *N1,struct node *N2,float y)
{
	struct site p;
	p=cal_break_point(N1->B,y);
//	cout<<y<<"\n";
     // display_BP(N1->B);
     // display_BP(N2->B);
//	display_SITE(p);
	if(N1->B->V==NULL)
	{
		N1->B->V= new vertice;
	}
	N1->B->V->x=p.x;
	N1->B->V->y=p.y;
	p=cal_break_point(N2->B,y);
//	display_SITE(p);
	if(N2->B->V==NULL)
	{
		N2->B->V= new vertice;
	}
	N2->B->V->x=p.x;
	N2->B->V->y=p.y;
	if(vor==NULL)
	{
		vor = new DCEL;
		vor->origin=N1->B->V;
		//cout<<N1->B.V->x<<" htis "<<N1->B.V->y<<"\n";
	}
//	cout<<"origin mana origin "<<vor->origin->x<<" "<<vor->origin->y<<"\n";
	N1->B->V->leaving=new half_edge;
	N1->B->V->leaving->origin=N1->B->V;
	N2->B->V->leaving=new half_edge;
	N2->B->V->leaving->origin=N2->B->V;
	N1->B->V->leaving->twin=N2->B->V->leaving;
	N2->B->V->leaving->twin=N1->B->V->leaving;
	if(N1->B->p1->F)
	{
		N1->B->V->leaving->facing=N1->B->p1->F;

	}
	else 
		N1->B->V->leaving->facing=new face;
	if(N2->B->p1->F)
	{
		N2->B->V->leaving->facing=N2->B->p1->F;

	}
	else 
		N2->B->V->leaving->facing=new face;
	N1->B->V->leaving->facing->edge=N1->B->V->leaving;
	N2->B->V->leaving->facing->edge=N2->B->V->leaving;
	N1->B->V->leaving->facing->p=N1->B->p1;
	N2->B->V->leaving->facing->p=N2->B->p1;
	if(!(N1->B->p1->F))
		N1->B->p1->F=N1->B->V->leaving->facing;
	if(!(N2->B->p1->F))
		N2->B->p1->F=N2->B->V->leaving->facing;
}
/* See if there is a circle event given three arcs
 */
void circlevent(struct node *L,struct node *M,struct node *R)
{
	float ax=M->p->x-L->p->x;
	float ay=M->p->y-L->p->y;
	float bx=R->p->x-L->p->x;
	float by=R->p->y-L->p->y;
	float A=ax*ax+ay*ay;				

	float B=bx*bx+by*by;
	float x=0.5*(ay*B-by*A)/(ay*bx-ax*by);
	float denom=2.*(ay*bx-ax*by);
	float y=-1.*ax/ay*x+A/(2.*ay);
	float dis=sqrt(pow(x-ax,2)+pow(y-ay,2));
////////cout<<"begin\n";
////////display_SITE(L->p);
////////display_SITE(M->p);
////////display_SITE(R->p);
////////cout<<x+L->p->x<<" "<<y+L->p->y<<"\n";
////////cout<<ay<<"\n";
////////cout<<denom<<"\n";
        if(denom > 0. && y+L->p->y > 0.)
        {
		event *circle;
////////	cout<<"begin\n";
////////	display_SITE(L->p);
////////	display_SITE(M->p);
////////	display_SITE(R->p);
		circle = new event;
		circle->p = new site;
		circle->center = new site;
		circle->p->x=x+L->p->x;
		circle->p->y=y+L->p->y-dis;
		circle->center->x=x+L->p->x;
		circle->center->y=y+L->p->y;
////////////////cout<<"circle_event=";
////////////////display_SITE(circle->p);
////////////////display_SITE(circle->center);
////////////////cout<<dis<<"\n";
////////////////cout<<"end\n";
		circle->circle_event=1;
		circle->circle_node=M;
		//M->circle_event=circle;
		P.insert_event(start,circle,M);
//		display_events(start);
		
        }
}
/* function to inset a node into the tree
 */		
void BST::insert(node *tree, node *newnode)
{
	if(root==NULL)
	{
		root = new node;
		root->p=newnode->p;	
		root->left=NULL;
		root->right=NULL;
		root->adj_left=NULL;
		root->adj_right=NULL;
		root->parent=NULL;
		return;
	}
	if( tree->left==NULL && tree->right == NULL )
	{
		if(tree->circle_event!=NULL)
		{
			//cout<<"circle event is false alarm\n";	
			//cout<<tree->circle_event->center->x<<" "<<tree->circle_event->center->y<<"\n";
			P.delete_event(start,tree->circle_event);
		}
		//creating a new node at left
		tree->left=new node;
		//assigning the parent to be the tree
		tree->left->parent=tree;
		//Moving the node down
		tree->left->p=tree->p;
		//assigning the left adjacent nodes for the left child
		tree->left->adj_left=tree->adj_left;
		if(tree->adj_left!=NULL)
			tree->adj_left->adj_right=tree->left;
		//No further branches  for the left child
		tree->left->left=NULL;
		tree->left->right=NULL;
		//creating a breakpoint at the parent node
		tree->B = new break_point;
		tree->B->p1=tree->left->p;
		tree->B->p2=newnode->p;
		tree->isBP=1;
		//making the pointer for the left child to all the breakpoints it belongs to!
		tree->left->R_breakpoint=tree;
		tree->left->L_breakpoint=tree->L_breakpoint;
		//creating a child at the right
		tree->right=new node;
		tree->right->twin=tree->left;
		tree->left->twin=tree->right;
		//assignging its parent to be the tree
		tree->right->parent=tree;
		//the right node also a breakpoint
		//Here we need to put code which can make a new vertex whose edges should point at the vertex before. 
		tree->right->B =  new break_point;
		tree->right->B->p1=newnode->p;
		tree->right->B->p2=tree->p;
		create_edge(tree,tree->right,newnode->p->y);
		tree->right->isBP=1;
		//Creating children at left and right
		tree->right->left=new node;
		tree->right->right=new node;
		tree->right->left->twin=tree->right->right;
		tree->right->right->twin=tree->right->left;
		//assigning their parents to be the right child of the tree
		tree->right->left->parent=tree->right;
		tree->right->right->parent=tree->right;
		//assigning the new node to the left child
		tree->right->left->p=newnode->p;
		//creating pointers to the break point for left child of the right child of tree
		tree->right->left->L_breakpoint=tree;
		tree->right->left->R_breakpoint=tree->right;
		//assigning the left adj node to be the tree's left child
		tree->right->left->adj_left=tree->left;	
		//tree's left child's right node is the tree's right child's left child
  		tree->left->adj_right=tree->right->left;
		//assigning the right child of the right child of the tree to have the tree's node value
		tree->right->right->p=tree->p;
		//creating pointer to breakpoint for right child of the right child.
		tree->right->right->L_breakpoint=tree->right;
		tree->right->right->R_breakpoint=tree->R_breakpoint;
		//assiging the left child of the right child of the tree adj right to be right child's right child
		tree->right->left->adj_right=tree->right->right;	
		tree->right->right->adj_left=tree->right->left;	
		tree->right->right->adj_right=tree->adj_right;	
		if(tree->adj_right!=NULL)
			tree->adj_right->adj_left=tree->right->right;
		tree->right->left->left=NULL;
		tree->right->left->right=NULL;
		tree->right->right->left=NULL;
		tree->right->right->right=NULL;
		if(tree->right->left->adj_right->adj_right != NULL)
		{
			//cout<<"whodunnit1\n";
			circlevent(tree->right->left,tree->right->left->adj_right,tree->right->left->adj_right->adj_right);
		}
		if(tree->right->left->adj_left->adj_left != NULL)
		{
			//cout<<"whodunnit2\n";
			circlevent(tree->right->left->adj_left->adj_left,tree->right->left->adj_left,tree->right->left);
		}
		return;
	}
	else
	{
		struct site bp;
		bp=cal_break_point(tree->B,newnode->p->y);
	//	cout<<bp.x<<"\n";
		if(newnode->p->x<bp.x)
			insert(tree->left,newnode);
		else
			insert(tree->right,newnode);
	}
}
void BST::delete_node(node *leaf, site *center)
{
//	cout<<leaf->p.x<<" "<<leaf->p.y<<"\n";
//	display_BP(leaf->L_breakpoint->B);
//	display_SITE(leaf->adj_left->p);
	//cout<<"welcome to the final frontier\n";
	
	//display_BP(leaf->L_breakpoint->B);
	leaf->L_breakpoint->B->p2=leaf->adj_right->p;
	leaf->R_breakpoint->B->p1=leaf->adj_left->p;
	//DCEL STUFF
	struct vertice *temp_vert;
	struct half_edge *temp_he;
	struct half_edge *temp_he2;
	temp_he = new half_edge;
	temp_he2 = new half_edge;
	temp_vert = new vertice;
	temp_vert->x=center->x;
	temp_vert->y=center->y;
	temp_he->origin=temp_vert;
	temp_he2->origin=temp_vert;
	temp_vert->leaving= new half_edge;
	temp_vert->leaving->origin=temp_vert;
	temp_vert->leaving->twin=leaf->L_breakpoint->B->V->leaving->twin;
	temp_vert->leaving->facing=leaf->L_breakpoint->B->V->leaving->facing;
	temp_vert->leaving->next=leaf->L_breakpoint->B->V->leaving->next;
	//leaf->L_breakpoint->B->V->leaving->facing->edge=temp_vert->leaving;
////////if(vor)
////////	cout<<vor->origin->leaving->twin->origin->x<<" "<<vor->origin->leaving->twin->origin->y<<"\n";
	//cout<<"where are we \n";
	leaf->L_breakpoint->B->V->leaving->twin->twin=temp_vert->leaving;
	//cout<<"not here\n";
	temp_he->twin=leaf->R_breakpoint->B->V->leaving->twin;
////////cout<<"address of right\n";
////////cout<<leaf->R_breakpoint->B->V->leaving->facing<<"\n";
////////cout<<leaf->R_breakpoint->B->V->leaving->facing->p->x<<" "<<leaf->R_breakpoint->B->V->leaving->facing->p->y<<"\n";
	temp_he->facing=leaf->R_breakpoint->B->V->leaving->facing;
	temp_he->next = leaf->R_breakpoint->B->V->leaving->next;
	//leaf->R_breakpoint->B->V->leaving->facing->edge=temp_he;
	leaf->R_breakpoint->B->V->leaving->facing->edge=leaf->L_breakpoint->B->V->leaving->twin;
	leaf->R_breakpoint->B->V->leaving->twin->twin=temp_he;
	temp_he2->facing = leaf->R_breakpoint->B->V->leaving->twin->facing;
	////////leaf->L_breakpoint->B.V->leaving->origin=temp_vert;
////////leaf->L_breakpoint->B.V->leaving->origin->leaving=leaf->L_breakpoint->B.V->leaving;
////////leaf->R_breakpoint->B.V->leaving->origin=temp_vert;
//	display_BP(leaf->L_breakpoint->B);
     // if(vor)
     // {
     // 	cout<<vor->origin<<"\n";
     // 	cout<<"after "<<vor->origin->x<<" "<<vor->origin->y<<"\n";
     // }
	//cout<<leaf->L_breakpoint->B.V<<"\n";
	
	temp_vert->leaving->twin->next=temp_he;
	leaf->R_breakpoint->B->V->leaving->twin->next=temp_he2;
	//DCEL STUFF ENDS

	if(leaf->parent == leaf->adj_right->L_breakpoint)
	{	
		leaf->adj_right->L_breakpoint=leaf->L_breakpoint;
		//#################################//
		leaf->L_breakpoint->B->V->x=center->x;
		leaf->L_breakpoint->B->V->y=center->y;
	    //  if(vor)
	    //  {
	    //  	//cout<<vor->origin<<"\n";
	    //  	cout<<"after2 "<<vor->origin->x<<" "<<vor->origin->y<<"\n";
	    //  }
		leaf->L_breakpoint->B->V->leaving->origin=leaf->L_breakpoint->B->V;
		leaf->L_breakpoint->B->V->leaving->twin=temp_he2;
		temp_he2->twin = leaf->L_breakpoint->B->V->leaving;
		leaf->L_breakpoint->B->V->leaving->next=temp_vert->leaving;
		leaf->L_breakpoint->B->V->leaving->facing->edge=leaf->L_breakpoint->B->V->leaving;
	////////cout<<"insideleft\n";
	////////cout<<leaf->L_breakpoint->B->V->leaving->facing<<"\n";
	////////cout<<leaf->L_breakpoint->B->V->leaving->facing->p->x<<" "<<leaf->L_breakpoint->B->V->leaving->facing->p->y<<"\n";
	////////leaf->L_breakpoint->B->V->leaving->next=temp_vert->leaving;
	////////temp_vert->leaving->twin->next=temp_he;
	////////leaf->R_breakpoint->B->V->leaving->twin->next=temp_he2;
	}
	else 
	{
////////	cout<<"right\n";
////////cout<<"edge\n";
////////cout<<temp_he->facing<<"\n";
////////cout<<temp_he->facing->p->x<<" "<<temp_he->facing->p->y<<"\n";
////////cout<<temp_he->facing->edge->origin->x<<" "<<temp_he->facing->edge->origin->y<<" ";
////////cout<<temp_he->facing->edge->twin->origin->x<<" "<<temp_he->facing->edge->twin->origin->y<<"\n";
	/////////.//////////////////////////////////////////////////////////////////
		leaf->adj_left->R_breakpoint=leaf->R_breakpoint;
		//#################################//
		leaf->R_breakpoint->B->V->x=center->x;
		leaf->R_breakpoint->B->V->y=center->y;
		leaf->R_breakpoint->B->V->leaving->origin=leaf->R_breakpoint->B->V;
		leaf->R_breakpoint->B->V->leaving->twin=temp_he2;
		leaf->R_breakpoint->B->V->leaving->facing=leaf->L_breakpoint->B->V->leaving->facing;
		//leaf->R_breakpoint->B->V->leaving->facing->edge=temp_vert->leaving->twin;
		temp_he2->twin = leaf->R_breakpoint->B->V->leaving;
		leaf->R_breakpoint->B->V->leaving->next=temp_vert->leaving;
		leaf->L_breakpoint->B->V->leaving->facing->edge=leaf->R_breakpoint->B->V->leaving;
////////cout<<"after\n";
////////cout<<leaf->R_breakpoint->B->V->leaving->facing<<"\n";
////////cout<<leaf->R_breakpoint->B->V->leaving->facing->p->x<<" "<<leaf->R_breakpoint->B->V->leaving->facing->p->y<<"\n";
////////cout<<leaf->R_breakpoint->B->V->leaving->facing->edge->origin->x<<" "<<leaf->R_breakpoint->B->V->leaving->facing->edge->origin->y<<" ";
////////cout<<leaf->R_breakpoint->B->V->leaving->facing->edge->twin->origin->x<<" "<<leaf->R_breakpoint->B->V->leaving->facing->edge->twin->origin->y<<"\n";
	}
	if(leaf->parent->parent->left==leaf->parent)
	{
		//cout<<"this didn't happen\n";
//		display(root,0.,0);
		leaf->parent->parent->left=leaf->twin;
		leaf->twin->twin=leaf->parent->parent->right;
		leaf->parent->parent->right->twin=leaf->twin;
////////	if(leaf->twin)
////////	{
////////		cout<<"really\n";
////////		display_SITE(leaf->p);
////////		cout<<leaf->twin->isBP<<"\n";
////////		display_BP(leaf->twin->B);
////////		display_SITE(leaf->twin->left->p);
////////		display_SITE(leaf->twin->right->p);
////////		display_SITE(leaf->twin->p);
////////	}
////////	display(root,0.,0);
	}
	else
	{
		leaf->parent->parent->right=leaf->twin;
		leaf->twin->twin=leaf->parent->parent->left;
		leaf->parent->parent->left->twin=leaf->twin;
	}
	leaf->twin->parent=leaf->parent->parent;
	leaf->adj_left->adj_right=leaf->adj_right;
	leaf->adj_right->adj_left=leaf->adj_left;
//	cout<<"serio\n";
	if(leaf->adj_right->circle_event != NULL)
	{
//		display_SITE(start->p);
		P.delete_event(start,leaf->adj_right->circle_event);
		leaf->adj_right->circle_event=NULL;
	}
	if(leaf->adj_left->circle_event != NULL)
	{
	//cout<<"where1\n";
		P.delete_event(start,leaf->adj_left->circle_event);
		leaf->adj_left->circle_event=NULL;
	////cout<<"where2\n";
	}
        if(leaf->adj_left->adj_left != NULL)
        {
		//cout<<"orthisguy?\n";
        	circlevent(leaf->adj_left->adj_left,leaf->adj_left,leaf->adj_right);	
        }
        if(leaf->adj_right->adj_right!= NULL)
        {
		//cout<<"orthisguy2?\n";
        	circlevent(leaf->adj_left,leaf->adj_right,leaf->adj_right->adj_right);	
        }
}

void BST::display(node *tree,float y,int level=0)
{
	if(tree==NULL)
		return;
	if(tree->right!=NULL)
	{
		display(tree->right,y,level+1);
	}
	for(int i=0;i<level;i++)
	{
		cout<<"      		";
	}
	if( tree->left!=NULL && tree->right != NULL )
	{
		struct site p;
		p=cal_break_point(tree->B,y);
                cout<<"("<<tree->B->p1->x<<" "<<tree->B->p1->y<<")";
                cout<<"("<<tree->B->p2->x<<" "<<tree->B->p2->y<<")";
		cout<<"("<<p.x<<","<<p.y<<")\n";

        /////////out<<std::setprecision(3)<<"("<<tree->B->p1->x<<" "<<tree->B->p1->y<<")";
        /////////out<<std::setprecision(3)<<"("<<tree->B->p2->x<<" "<<tree->B->p2->y<<")";
	/////////out<<std::setprecision(3)<<"("<<p.x<<","<<p.y<<")\n";
		//cout<<level<<"(X,X)\n";
	}
	else
	{
		cout<<tree->circle_event<<" ";
		cout<<"N("<<tree->p->x<<","<<tree->p->y<<")\n";
////////	if(tree->twin)
////////	{
////////		if(!tree->twin->isBP)
////////			cout<<"N("<<tree->twin->p->x<<","<<tree->twin->p->y<<")\n";
////////		else
////////			cout<<"bp\n";
////////	}
	}
	if(tree->left!=NULL)
	{
		display(tree->left,y,level+1);
	}
}
void BST::disBeach(node *tree)
{
	struct site bp;
	if( tree->left==NULL && tree->right == NULL )
	{
//		cout<<"("<<tree->p.x<<","<<tree->p.y<<")\n";
		cout<<"#######\n";
		if(tree->L_breakpoint)
		display_BP(tree->L_breakpoint->B);
	        display_SITE(tree->p);

	//        display_BP(tree->parent->B);
////////        if(tree->adj_left!=NULL)
////////	{
////////		cout<<"reverse\n";
////////        	display_SITE(tree->adj_left->p);
////////		if( tree->adj_left->left != NULL && tree->adj_left->right != NULL )
////////			cout<<"nono\n";
	//	}
	////////display_SITE(tree->adj_right->p);
	////////display_SITE(tree->adj_right->adj_right->p);
	////////display_SITE(tree->adj_right->adj_right->adj_right->p);
	////////cout<<"####\n";
	////////bp=cal_break_point(tree->parent->B,1.);	
	////////cout<<bp.x<<"\n";
		if(tree->adj_right != NULL )
		{
////////		if( tree->adj_right->left != NULL && tree->adj_right->right != NULL )
////////		{
////////        		display_SITE(tree->adj_right->p);
////////        		display_SITE(tree->adj_right->adj_left->p);
////////        		display_BP(tree->adj_right->B);
////////			cout<<tree->adj_right->isBP<<"\n";
////////			cout<<"error\n";
////////			return;
////////		}
			disBeach(tree->adj_right);
		}
	}
	else 
	disBeach(tree->left);
}
void priority_list::insert_event(event *EV,event *newevent,node *M=NULL)
{
	int flag=0;
	struct event *Nevent;
////////cout<<newevent->circle_event<<"this is check\n";
      //if(start != NULL)
      //	display_SITE(EV->p);
	if(start==NULL)
	{
		start=new event;
		start->p=newevent->p;
		//display_SITE(newevent->p);
	}
        else
        {
        	flag=compare(EV->p,newevent->p);
		//display_SITE(newevent->p);
        	if(flag==-1)
        	{
        		Nevent=new event;
        		Nevent->p=newevent->p;
			if(newevent->p->y == EV->p->y)
			{
				//cout<<"is this really happenign?\n";
				if(EV->next)
        				insert_event(EV->next,newevent,M);
				else 
				{
					EV->next=new event;
					EV->next->p=newevent->p;
					EV->next->circle_event=newevent->circle_event;
					EV->next->circle_node=M;
					EV->next->center=newevent->center;
					EV->next->prev=EV;
					if(newevent->circle_event)
						M->circle_event=EV->next;
				}
				return;

			}
	        	if(EV->prev)
	        	{
                		EV->prev->next=Nevent;
				Nevent->prev=EV->prev;
	        	}
			else if (EV == start)
				start=Nevent;
	////////	else 
	////////	{

	////////		start=Nevent;
	////////	}
        		Nevent->next=EV;
        		EV->prev=Nevent;
			if(newevent->circle_event)
			{
				//cout<<"start(incide)=";
				//display_SITE(start->p);
        			//cout<<"this event?";
        			//display_SITE(Nevent->p);
				Nevent->circle_event=1;
				Nevent->circle_node=M;
				M->circle_event=Nevent;
				Nevent->center=newevent->center;
			}
        	}
        	else if (flag==1)
        	{
			if(EV->next != NULL)
        			insert_event(EV->next,newevent,M);
			else 
			{
				EV->next=new event;
				EV->next->p=newevent->p;
				EV->next->circle_event=newevent->circle_event;
				EV->next->circle_node=M;
				EV->next->center=newevent->center;
				EV->next->prev=EV;
				if(newevent->circle_event)
					M->circle_event=EV->next;
			}
        	}
        	else if (flag==0)
		{
        		return;
		}
        }

}
void priority_list::delete_event(event *EV,event *eventtd)
{
	int flag;
        flag=compare(EV->p,eventtd->p);
      //cout<<"event=";
      //display_SITE(EV->p);
      //cout<<"tobedeleted=";
      //display_SITE(eventtd->p);
	if(flag==0)
	{
		if(EV->next)
		{
			EV->next->prev=EV->prev;
			EV->prev->next=EV->next;
		}
		else 
			EV->prev->next=NULL;
//		cout<<"what\n";
	}
	else
		delete_event(EV->next,eventtd);
	return;
}
void update_dcel(node *tree,float y)
{
	if(tree->left!=NULL && tree->right!=NULL)
		update_dcel(tree->left,y);
	else 
	{
		struct site p;
		if(tree->R_breakpoint)
		{
			p=cal_break_point(tree->R_breakpoint->B,y);
			{
			tree->R_breakpoint->B->V->x=p.x;
			tree->R_breakpoint->B->V->y=p.y;
			//cout<<"check this out\n";
	//		display_SITE(&p);
			}
			if(tree->adj_right)
				update_dcel(tree->adj_right,y);
		}
	}
}
void display_a_edge(struct half_edge *first,struct half_edge *E)
{
	ofstream DEL;
	DEL.open("dcel",ios::app);
	DEL<<"#"<<E->facing->p->x<<" "<<E->facing->p->y<<"\n";
	DEL<<E->origin->x<<"\t"<<E->origin->y<<"\t"<<E->twin->origin->x<<"\t"<<E->twin->origin->y<<"\n";
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
		return;

}
void display_dcel(struct event *EV)
{
	if(!(EV->circle_event))
	{
////////	cout<<"#";
////////	//display_SITE(EV->p);
////////	cout<<EV->p->x<<" "<<EV->p->y<<"\n";
//		cout<<EV->p->F<<"\n";
		if(EV->p->F)
		{
			display_a_edge(EV->p->F->edge,EV->p->F->edge);
		}
	}
	if(EV->next)
		display_dcel(EV->next);
	else
		return;
		
}
void priority_list::read_events(event *EV)
{
	node *temp;
	temp = new node;
	temp->p=EV->p;
	//cout<<"start=";
	//display_SITE(start->p);
	if(root != NULL)
	{
	//	update_dcel(root,EV->p->y);
	        //bst.display(root,EV->p->y);
	//	cout<<EV->p->y<<"\n";
	//	display_dcel(start);
		//display_events(start);
		//cout<<"######\n";
	////////if(start->p->F)
	////////	display_a_edge(start->p->F->edge,start->p->F->edge);
		//display_SITE(EV->p);
	      //cout<<"###################\n";
	        //bst.disBeach(root);
	}
  //    if(vor)
  //    {
  //    	//cout<<vor->origin<<"\n";
  //    	//cout<<"before "<<vor->origin->x<<" "<<vor->origin->y<<"\n";
  //    }
	if(EV->circle_event)
	{
    	 // 	display_SITE(EV->p);
	////////cout<<"circle_event=";
	////////display_SITE(EV->center);
		//cout<<"here1\n";
		bst.delete_node(EV->circle_node,EV->center);
		//display_events(start);
		//cout<<"here2\n";
	}
	else 
	{
		//cout<<"he1\n";
		bst.insert(root,temp);
		//cout<<"he2\n";
	}
	if(EV->next)
	{
		read_events(EV->next);
	}
	else 
	{
	        if(vor)
		{
		////////cout<<vor->origin<<"\n";
	        ////////cout<<"origin "<<vor->origin->x<<" "<<vor->origin->y<<"\n";
		}
		
	//	cout<<"done bitches\n";
	  	update_dcel(root,EV->p->y);
	}
}


int main()
{
	node *temp;
	temp = new node;
	event *newevent;
	newevent=new event;	
	std::ifstream INFILE("dat");
	float x,y;
	while(INFILE>>x>>y)
	{
		newevent->p=new site;
		newevent->p->x=x;
		newevent->p->y=y;
		P.insert_event(start,newevent);
        ////////temp->p.x=x;
	////////temp->p.y=y;
//        	bst.insert(root,temp);
	}
////////cout<<"over\n";
////////cout<<start<<"\n";
        display_events(start);
	P.read_events(start);
  //      display_events(start);
        display_events(start);
//        bst.disBeach(root);	
	display_dcel(start);
//	bst.display(root,-4.34269);
	return 0;
}
