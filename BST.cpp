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
	int isBP=0;
}*root;

struct event
{
	struct site p;
	struct event *next;
};
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
////////if(y < ay)
////////{
////////	include_circle_event();
////////}
}
		
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
class BST
{
	public:
		void insert(node *, node *);
		void delete_node(node *);
		void display(node *,int);
		void find(node *,int);
		void disBeach(node *);
};
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
int main()
{
        BST bst;
	node *temp;
	temp = new node;
	std::ifstream INFILE("sites");
	float x,y;
	while(INFILE>>x>>y)
	{
        	temp->p.x=x;
		temp->p.y=y;
        	bst.insert(root,temp);
	}

//	cout<<temp->data<<"\n";
      //temp->p.x=1.;
      //temp->p.y=5.;
////////cout<<"hello\n";
//    //temp->p.x=1.0;
//    //temp->p.y=1.0;
      //bst.insert(root,temp);
      //temp->p.x=0.;
      //temp->p.y=4.;
      //cout<<"hello\n";
      //bst.insert(root,temp);
      //temp->p.x=-1.;
      //temp->p.y=1.;
      //bst.insert(root,temp);
      //bst.display(root);
//	bst.disBeach(root);
	//cout<<root->left->adj_right->p.x<<"\t";
//	display_SITE(root->left->adj_right->adj_right->p);
//        bst.display(root);
        bst.display(root);
	cout<<"######################\n";
        bst.disBeach(root);
        bst.delete_node(root->left->adj_right);//->adj_right->adj_right);
        bst.display(root);
     // cout<<"######################\n";
        bst.disBeach(root);
	return 0;
}
