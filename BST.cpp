#include<iostream>
#include<math.h>
#include<fstream>
/*
 * Program to implement a binary search tree 
 */
using namespace std;
//Site is a point on the 2D surface
struct site
{
	float x=0;
	float y=0;
};
//Breakpoint definition
struct break_point
{
	struct site p1,p2;
};
//node definition 
struct node 
{
        struct site p;
        struct break_point B;
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

struct event
{
	struct site p;
	int circle_event=0;
	struct node *circle_node;
	struct event *next=NULL;
	struct event *prev=NULL;
}*start;
class BST
{
	public:
		void insert(node *, node *);
		void delete_node(node *);
		void display(node *,int);
		void find(node *,int);
		void disBeach(node *);
}bst;

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
int compare(struct site p1,struct site p2)
{
	if( p1.y > p2.y )
		return 1;
	else if ( (p1.y == p2.y) && ( p1.x > p2.x ) )
		return 1;
	else if ( (p1.y == p2.y) && ( p1.x = p2.x ) )
		return 0;
	else 
		return -1;
}


struct site cal_break_point(struct break_point B,float y)
{
	float h=B.p1.y-y;
	float h1=B.p2.y-y;
	float x1=B.p2.x-B.p1.x;
	int flag=compare(B.p1,B.p2);
	struct site BP;
////////cout<<B.p1.x<<"\t"<<B.p1.y<<"\n";
////////cout<<B.p2.x<<"\t"<<B.p2.y<<"\n";
////////cout<<x1<<"\t"<<h<<"\t"<<h1<<"\n";
	if(flag==1)
	{
		BP.x=(-2.*h*x1+sqrt(pow(2.*h*x1,2)-4.*(h1-h)*((h-h1)*h*h1-h*pow(x1,2))))/(2.*(h1-h))+B.p1.x;
		if(BP.x<B.p2.x)
		{
			BP.y=pow(BP.x,2)/(2.*h)+2.*h/2.;
			return BP;
		}
		else 
		{
			BP.x=(-2.*h*x1-sqrt(pow(2.*h*x1,2)-4.*(h1-h)*((h-h1)*h*h1-h*pow(x1,2))))/(2.*(h1-h))+B.p1.x;
			BP.y=pow(BP.x,2)/(2.*h)+2.*h/2.;
			return BP;
		}
	}
	else if (flag == -1)
	{
		BP.x=(-2.*h*x1+sqrt(pow(2.*h*x1,2)-4.*(h1-h)*((h-h1)*h*h1-h*pow(x1,2))))/(2.*(h1-h))+B.p1.x;
		if(BP.x>B.p1.x)
		{
			BP.y=pow(BP.x,2)/(2.*h)+2.*h/2.;
			return BP;
		}
		else 
		{
			BP.x=(-2.*h*x1-sqrt(pow(2.*h*x1,2)-4.*(h1-h)*((h-h1)*h*h1-h*pow(x1,2))))/(2.*(h1-h))+B.p1.x;
			BP.y=pow(BP.x,2)/(2.*h)+2.*h/2.;
			return BP;
		}
	}
}
/*
 * Defining a Class called Binary Search Tree
 */
void display_BP(struct break_point B)
{
	cout<<"##########################\n";
	cout<<B.p1.x<<" "<<B.p1.y<<"\n";
	cout<<B.p2.x<<" "<<B.p2.y<<"\n";
	cout<<"##########################\n";
}
void display_SITE(struct site p)
{
	cout<<"("<<p.x<<" "<<p.y<<")\n";
}

void circlevent(struct node *L,struct node *M,struct node *R)
{
	float ax=M->p.x-L->p.x;
	float ay=M->p.y-L->p.y;
	float bx=R->p.x-L->p.x;
	float by=R->p.y-L->p.y;
	float A=ax*ax+ay*ay;
	float B=bx*bx+by*by;
	float x=0.5*(ay*B-by*A)/(ay*bx-ax*by);
	float y=-1.*ax/ay*x+A/(2.*ay);
	float dis=sqrt(pow(x-ax,2)+pow(y-ay,2));
        if(y < ay)
        {
		event *circle;
		circle = new event;
		circle->p.x=x+L->p.x;
		circle->p.y=y+L->p.y-dis;
		display_SITE(circle->p);
		circle->circle_event=1;
		circle->circle_node=M;
		P.insert_event(start,circle,M);
        }
}
		
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
			cout<<"circle event is false alarm\n";	
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
		tree->B.p1=tree->left->p;
		tree->B.p2=newnode->p;
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
		tree->right->B.p1=newnode->p;
		tree->right->B.p2=tree->p;
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
			circlevent(tree->right->left,tree->right->left->adj_right,tree->right->left->adj_right->adj_right);
		if(tree->right->left->adj_left->adj_left != NULL)
			circlevent(tree->right->left->adj_left->adj_left,tree->right->left->adj_left,tree->right->left);
		return;
	}
	else
	{
		struct site bp;
		bp=cal_break_point(tree->B,newnode->p.y);
	//	cout<<bp.x<<"\n";
		if(newnode->p.x<bp.x)
			insert(tree->left,newnode);
		else
			insert(tree->right,newnode);
	}
}
void BST::delete_node(node *leaf)
{
//	cout<<leaf->p.x<<" "<<leaf->p.y<<"\n";
//	display_BP(leaf->L_breakpoint->B);
//	display_SITE(leaf->adj_left->p);
	leaf->L_breakpoint->B.p2=leaf->adj_right->p;
	leaf->R_breakpoint->B.p1=leaf->adj_left->p;
//	display_BP(leaf->L_breakpoint->B);
	if(leaf->parent == leaf->adj_right->L_breakpoint)
		leaf->adj_right->L_breakpoint=leaf->L_breakpoint;
	else 
		leaf->adj_left->R_breakpoint=leaf->R_breakpoint;
	if(leaf->parent->parent->left==leaf->parent)
		leaf->parent->parent->left=leaf->twin;
	else
		leaf->parent->parent->right=leaf->twin;
	leaf->twin->parent=leaf->parent->parent;
	leaf->adj_left->adj_right=leaf->adj_right;
	leaf->adj_right->adj_left=leaf->adj_left;
	if(leaf->adj_right->circle_event != NULL)
		P.delete_event(start,leaf->adj_right->circle_event);
	if(leaf->adj_left->circle_event != NULL)
		P.delete_event(start,leaf->adj_left->circle_event);
}


void BST::display(node *tree,int level=0)
{
	if(tree==NULL)
		return;
	if(tree->right!=NULL)
	{
		display(tree->right,level+1);
	}
	for(int i=0;i<level;i++)
	{
		cout<<"      ";
	}
	if( tree->left!=NULL && tree->right != NULL )
		cout<<level<<"(X,X)\n";
	else
		cout<<"("<<tree->p.x<<","<<tree->p.y<<")\n";
	if(tree->left!=NULL)
	{
		display(tree->left,level+1);
	}
}
void BST::disBeach(node *tree)
{
	struct site bp;
	if( tree->left==NULL && tree->right == NULL )
	{
//		cout<<"("<<tree->p.x<<","<<tree->p.y<<")\n";
		cout<<"#######\n";
	        display_SITE(tree->p);

	//        display_BP(tree->parent->B);
	        if(tree->adj_left!=NULL)
		{
			cout<<"reverse\n";
	        	display_SITE(tree->adj_left->p);
			if( tree->adj_left->left != NULL && tree->adj_left->right != NULL )
				cout<<"nono\n";
		}
	////////display_SITE(tree->adj_right->p);
	////////display_SITE(tree->adj_right->adj_right->p);
	////////display_SITE(tree->adj_right->adj_right->adj_right->p);
	////////cout<<"####\n";
	////////bp=cal_break_point(tree->parent->B,1.);	
	////////cout<<bp.x<<"\n";
		if(tree->adj_right != NULL )
		{
			if( tree->adj_right->left != NULL && tree->adj_right->right != NULL )
			{
	        		display_SITE(tree->adj_right->p);
	        		display_SITE(tree->adj_right->adj_left->p);
	        		display_BP(tree->adj_right->B);
				cout<<tree->adj_right->isBP<<"\n";
				cout<<"error\n";
				return;
			}
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
	if(start==NULL)
	{
		start=new event;
		start->p=newevent->p;
	}
        else
        {
        	flag=compare(EV->p,newevent->p);
        	if(flag==-1)
        	{
        		Nevent=new event;
        		Nevent->p=newevent->p;
        		EV->prev->next=Nevent;
        		Nevent->next=EV;
        		EV->prev=Nevent;
			if(newevent->circle_event)
			{
				Nevent->circle_event=1;
				Nevent->circle_node=M;
				M->circle_event=Nevent;
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
				EV->next->prev=EV;
			}
        	}
        	else if (flag==0)
        		return;
        }

}
void priority_list::delete_event(event *EV,event *eventtd)
{
	int flag;
        flag=compare(EV->p,eventtd->p);
	if(flag==0)
	{
		EV->next->prev=EV->prev;
		EV->prev->next=EV->next;
	}
	else
		delete_event(EV->next,eventtd);
	return;
}
void priority_list::read_events(event *EV)
{
	node *temp;
	temp = new node;
	temp->p=EV->p;
	if(EV->circle_event)
	{
		bst.delete_node(EV->circle_node);
	}
	else 
		bst.insert(root,temp);
	if(EV->next)
		read_events(EV->next);
	else 
		cout<<"done bitches\n";
}
void display_events(struct event *start)
{
	display_SITE(start->p);
	if(start->next!=NULL)
	{
		if(start->circle_event)
			cout<<"circle event= ";
		display_events(start->next);
	}
	else
		cout<<"events over\n";
}

int main()
{
	node *temp;
	temp = new node;
	event *newevent;
	newevent=new event;	
	std::ifstream INFILE("sites");
	float x,y;
	while(INFILE>>x>>y)
	{
		newevent->p.x=x;
		newevent->p.y=y;
		P.insert_event(start,newevent);
        ////////temp->p.x=x;
	////////temp->p.y=y;
//        	bst.insert(root,temp);
	}
	display_events(start);
	P.read_events(start);
//	display_events(start);
	return 0;
}
