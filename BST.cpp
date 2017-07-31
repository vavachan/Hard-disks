#include<iostream>
#include<math.h>
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
}*root;
/*
 * Defining a Class called Binary Search Tree
 */
class BST
{
	public:
		void insert(node *, node *);
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
		tree->left=new node;
		tree->left->parent=tree;
		//Moving the node down
		tree->left->p=tree->p;
		//assigning the left and right adjacent nodes
		tree->left->adj_left=tree->adj_left;
//		tree->left->adj_right=newnode;
		//No further branches 
		tree->left->left=NULL;
		tree->left->right=NULL;
		//creating a breakpoint
		tree->B.p1=tree->p;
		tree->B.p2=newnode->p;
		tree->right=new node;
		tree->right->parent=tree;
		//the right node also a breakpoint
		tree->right->B.p1=newnode->p;
		tree->right->B.p2=tree->p;
		tree->right->left=new node;
		tree->right->right=new node;
		tree->right->left->parent=tree->right;
		tree->right->right->parent=tree->right;
		tree->right->left->p=newnode->p;
		tree->right->left->adj_left=tree->left;	
  		tree->left->adj_right=tree->right->left;
		tree->right->right->p=tree->p;
		tree->right->left->adj_right=tree->right->right;	
		tree->right->right->adj_left=tree->right->left;	
		tree->right->right->adj_right=tree->adj_right;	
		tree->right->left->left=NULL;
		tree->right->left->right=NULL;
		tree->right->right->left=NULL;
		tree->right->right->right=NULL;
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
		cout<<"("<<tree->p.x<<","<<tree->p.y<<")\n";
	////////bp=cal_break_point(tree->parent->B,1.);	
	////////cout<<bp.x<<"\n";
		disBeach(tree->adj_right);
	}
	else 
	disBeach(tree->left);
}
int main()
{
        BST bst;
	node *temp;
	temp = new node;
//	cout<<temp->data<<"\n";
        temp->p.x=1.;
	temp->p.y=5.;
////////cout<<"hello\n";
//	temp->p.x=1.0;
//	temp->p.y=1.0;
        bst.insert(root,temp);
        temp->p.x=0.;
	temp->p.y=4.;
        cout<<"hello\n";
        bst.insert(root,temp);
        temp->p.x=-1.;
	temp->p.y=1.;
        bst.insert(root,temp);
        bst.display(root);
	bst.disBeach(root);
	return 0;
}
