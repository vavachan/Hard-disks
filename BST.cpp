#include<iostream>
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
int compare(struct site *p1,struct site *p2)
{
	if( p1->y > p2->y )
		return 1;
	else if ( (p1->y == p2->y) && ( p1->x > p2->x ) )
		return 1;
	else if ( (p1->y == p2->y) && ( p1->x = p2->x ) )
		return 0;
	else 
		return -1;
}
//node definition 
struct node 
{
	int data=1;
        struct site p;
        struct break_point B;
        struct node *left; 
        struct node *right;
}*root;
/*
 * Defining a Class called Binary Search Tree
 */
class BST
{
	public:
		void insert(node *, node *);
		void display(node *,int);
};
void BST::insert(node *tree, node *newnode)
{
	if(root==NULL)
	{
		root = new node;
		root->data=newnode->data;	
		root->left=NULL;
		root->right=NULL;
		return;
	}
	if(tree->data>newnode->data)
	{
		if(tree->left==NULL)
		{
			tree->left = new node;
			tree->left->data=newnode->data;
			tree->left->left=NULL;
			tree->left->right=NULL;
			return;
		}
		else
		{
			insert(tree->left,newnode);
		}
	}
	else if(tree->data<newnode->data)
	{
		if(tree->right==NULL)
		{
			tree->right = new node;
			tree->right->data=newnode->data;
			tree->right->left=NULL;
			tree->right->right=NULL;
			return;
		}
		else
		{
			insert(tree->right,newnode);
		}
	}
	else
		cout<<"file already exists\n";
	return;
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
	cout<<tree->data<<"\n";
	if(tree->left!=NULL)
	{
		display(tree->left,level+1);
	}
}
int main()
{
        BST bst;
	node *temp;
	temp = new node;
//	cout<<temp->data<<"\n";
////////cout<<"hello\n";
        temp->data=2;
////////cout<<"hello\n";
//	temp->p.x=1.0;
//	temp->p.y=1.0;
        bst.insert(root,temp);
        temp->data=3;
        bst.insert(root,temp);
        temp->data=1;
        bst.insert(root,temp);
        temp->data=4;
        bst.insert(root,temp);
        temp->data=7;
        bst.insert(root,temp);
        temp->data=6;
        bst.insert(root,temp);
        bst.display(root);
	return 0;
}
