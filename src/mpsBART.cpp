#include "mpsBART.h"
#include <random>
#include <RcppArmadillo.h>
using namespace std;

// =====================================
// Statistics Function
// =====================================
// Calculating the log-density of a MVN(0, Sigma)
//[[Rcpp::export]]
double log_dmvn(arma::vec& x, arma::mat& Sigma){

        arma::mat L = arma::chol(Sigma,"lower");
        arma::vec D = L.diag();
        double p = Sigma.n_cols;

        arma::vec z(p);
        double out;
        double acc;

        for(int ip=0;ip<p;ip++){
                acc = 0.0;
                for(int ii = 0; ii < ip; ii++){
                        acc += z(ii)*L(ip,ii);
                }
                z(ip) = (x(ip)-acc)/D(ip);
        }
        out = (-0.5*sum(square(z))-( (p/2.0)*log(2.0*M_PI) +sum(log(D)) ));


        return out;

};

//[[Rcpp::export]]
arma::mat sum_exclude_col(arma::mat mat, int exclude_int){

        // Setting the sum matrix
        arma::mat m(mat.n_rows,1);

        if(exclude_int==0){
                m = sum(mat.cols(1,mat.n_cols-1),1);
        } else if(exclude_int == (mat.n_cols-1)){
                m = sum(mat.cols(0,mat.n_cols-2),1);
        } else {
                m = arma::sum(mat.cols(0,exclude_int-1),1) + arma::sum(mat.cols(exclude_int+1,mat.n_cols-1),1);
        }

        return m;
}

// =====================================
// SPLINES FUNCTIONS
// (those functions are with respect to
//every modification about splines)
// =====================================
// Creating the pos function
double pos(double x, double x_i){

        double dif = (x-x_i) ;

        // Getting the positive part only
        if( dif> 0 ){
                return dif;
        } else {
                return 0.0;
        }
}



// Creating the pos function
arma::vec pos_vec(arma::vec x, double x_i){

        arma::vec dif(x.size());

        for(int i = 0; i < x.size(); i++){

                // Getting the positive part only
                if( (x(i)-x_i)>0){
                        dif(i) = x(i)-x_i;
                } else {
                        dif(i) = 0.0;
                }
        }

        return dif;

}

double pos_val(double x, double x_i){

        double dif = (x-x_i) ;

        // Getting the positive part only
        if( dif> 0 ){
                return dif;
        } else {
                return 0.0;
        }
}

// Create a function to generate matrix D (for penalisation)
arma::mat D(modelParam data){

        // Creating the matrix elements
        arma::mat D_m((data.p-2),data.p,arma::fill::zeros);

        for(int i=0;i<(data.p-2);i++){
                D_m(i,i) = 1;
                D_m(i,i+1) = -2;
                D_m(i,i+2) = 1;
        }

        return D_m;
}

// Create a function to generate matrix D (for penalisation)
arma::mat D_first(modelParam data){

        // Creating the matrix elements
        arma::mat D_m((data.p-1),data.p,arma::fill::zeros);

        for(int i=0;i<(data.p-1);i++){
                D_m(i,i) = -1;
                D_m(i,i+1) = 1;
        }

        return D_m;
}

// Create a function to generate matrix D (identity version)
arma::mat D_diag(modelParam data){

        // Creating the matrix elements
        arma::mat D_m(data.p,data.p,arma::fill::zeros);

        for(int i=0;i<data.p;i++){
                D_m(i,i) = 1;
        }

        return D_m;
}


// Function to generate B
//[[Rcpp::export]]
arma::mat bspline(arma::vec x,
                  arma::vec x_obs){

        arma::mat B(x.size(), x_obs.n_rows+1, arma::fill::ones);

        // Storing a copy
        arma::vec x_obs_copy = x_obs;
        arma::uword max_ind = x_obs.index_max();
        x_obs_copy(max_ind) = -std::numeric_limits<double>::infinity();
        double x_n_1 = x_obs_copy.max();

        // Setting the values for all columns
        B.col(1) = x;
        double x_n = max(x_obs);

        for(int i = 0; i < x.size(); i++) {

                for(int j = 0 ; j < (x_obs.size()-1); j++){
                                B(i,j+2) = (pow(pos_val(x(i),x_obs(j)),3) - pow(pos_val(x(i),x_n),3))/(x_n-x_obs(j)) - (pow(pos_val(x(i),x_n_1),3)-pow(pos_val(x(i),x_n),3))/(x_n-x_n_1);


                }

        }
        return B;
}



// Initialising the model Param
modelParam::modelParam(arma::mat x_train_,
                       arma::vec y_,
                       arma::mat x_test_,
                       arma::cube Z_train_,
                       arma::cube Z_test_,
                       int n_tree_,
                       double alpha_,
                       double beta_,
                       double tau_mu_,
                       double tau_b_,
                       double tau_b_intercept_,
                       double tau_,
                       double a_tau_,
                       double d_tau_,
                       double nu_,
                       double delta_,
                       double a_delta_, double d_delta_,
                       double n_mcmc_,
                       double n_burn_,
                       arma::vec p_sample_,
                       arma::vec p_sample_levels_){

        // Assign the variables
        x_train = x_train_;
        y = y_;
        x_test = x_test_;
        Z_train = Z_train_;
        Z_test = Z_test_;
        n_tree = n_tree_;
        p = Z_train_.n_cols;
        d_var = Z_train_.n_slices;
        alpha = alpha_;
        beta = beta_;
        tau_mu = tau_mu_;
        tau_b = arma::vec(Z_train_.n_slices,arma::fill::ones)*tau_b_;
        tau_b_intercept = tau_b_intercept_;
        tau = tau_;
        a_tau = a_tau_;
        d_tau = d_tau_;
        nu = nu_;
        delta = arma::vec(Z_train_.n_slices,arma::fill::ones)*delta_;
        a_delta = a_delta_;
        d_delta = d_delta_;
        n_mcmc = n_mcmc_;
        n_burn = n_burn_;
        p_sample = p_sample_;
        p_sample_levels = p_sample_levels_;

}

// Initialising a node
Node::Node(modelParam &data){
        isLeaf = true;
        isRoot = true;
        left = NULL;
        right = NULL;
        parent = NULL;
        train_index = arma::zeros<arma::vec>(data.x_train.n_rows);
        test_index = arma::zeros<arma::vec>(data.x_test.n_rows) ;
        train_index.fill(-1);
        test_index.fill(-1);

        // Initialising those vectors
        for(int i = 0; i< data.x_train.n_rows; i ++ ) {
                train_index[i] = -1;
        }
        for(int i = 0; i< data.x_test.n_rows; i ++ ) {
                test_index[i] = -1;
        }


        var_split = 0;
        var_split_rule = 0.0;
        lower = 0.0;
        upper = 1.0;
        curr_weight = 0.0;
        mu = 0.0;
        n_leaf = 0.0;
        r_sq_sum = 0.0;
        r_sum = 0.0;
        log_likelihood = 0.0;
        depth = 0;


        // Initialising all the parameters
        betas = arma::mat(data.p,data.d_var,arma::fill::zeros);


}

Node::~Node() {
        if(!isLeaf) {
                delete left;
                delete right;
        }
}

// Initializing a stump
void Node::Stump(modelParam& data){

        // Changing the left parent and right nodes;
        left = this;
        right = this;
        parent = this;
        // n_leaf  = data.x_train.n_rows;

        // Updating the training index with the current observations
        for(int i=0; i<data.x_train.n_rows;i++){
                train_index[i] = i;
        }

        // Updating the same for the test observations
        for(int i=0; i<data.x_test.n_rows;i++){
                test_index[i] = i;
        }

}

void Node::addingLeaves(modelParam& data){

     // Create the two new nodes
     left = new Node(data); // Creating a new vector object to the
     right = new Node(data);
     isLeaf = false;

     // Modifying the left node
     left -> isRoot = false;
     left -> isLeaf = true;
     left -> left = left;
     left -> right = left;
     left -> parent = this;
     left -> var_split = 0;
     left -> var_split_rule = 0.0;
     left -> lower = 0.0;
     left -> upper = 1.0;
     left -> mu = 0.0;
     left -> r_sq_sum = 0.0;
     left -> r_sum = 0.0;
     left -> log_likelihood = 0.0;
     left -> n_leaf = 0.0;
     left -> depth = depth+1;
     left -> train_index = arma::zeros<arma::vec>(data.x_train.n_rows);
     left -> test_index = arma::zeros<arma::vec>(data.x_test.n_rows);
     left -> train_index.fill(-1);
     left -> test_index.fill(-1);

     right -> isRoot = false;
     right -> isLeaf = true;
     right -> left = right; // Recall that you are saving the address of the right node.
     right -> right = right;
     right -> parent = this;
     right -> var_split = 0;
     right -> var_split_rule = 0.0;
     right -> lower = 0.0;
     right -> upper = 1.0;
     right -> mu = 0.0;
     right -> r_sq_sum = 0.0;
     right -> r_sum = 0.0;
     right -> log_likelihood = 0.0;
     right -> n_leaf = 0.0;
     right -> depth = depth+1;
     right -> train_index = arma::zeros<arma::vec>(data.x_train.n_rows);
     right -> test_index = arma::zeros<arma::vec>(data.x_test.n_rows);
     right -> train_index.fill(-1);
     right -> test_index.fill(-1);


     return;

}

// Creating boolean to check if the vector is left or right
bool Node::isLeft(){
        return (this == this->parent->left);
}

bool Node::isRight(){
        return (this == this->parent->right);
}

// Sample var
void Node::sampleSplitVar(modelParam &data){

          // Sampling one index from 0:(p-1)
          int original_index = std::rand()%data.p_sample.size();
          int new_sample;
          if(data.p_sample_levels[original_index]==0){
               var_split = data.p_sample[original_index];
          } else {
               new_sample = std::rand()%(int)data.p_sample_levels[original_index];
               var_split = data.p_sample[original_index] + data.p_sample_levels[original_index];
          }

}
// This functions will get and update the current limits for this current variable
void Node::getLimits(){

        // Creating  a new pointer for the current node
        Node* x = this;
        // Already defined this -- no?
        lower = 0.0;
        upper = 1.0;
        // First we gonna check if the current node is a root or not
        bool tree_iter = x->isRoot ? false: true;
        while(tree_iter){
                bool is_left = x->isLeft(); // This gonna check if the current node is left or not
                x = x->parent; // Always getting the parent of the parent
                tree_iter = x->isRoot ? false : true; // To stop the while
                if(x->var_split == var_split){
                        tree_iter = false ; // This stop is necessary otherwise we would go up til the root, since we are always update there is no prob.
                        if(is_left){
                                upper = x->var_split_rule;
                                lower = x->lower;
                        } else {
                                upper = x->upper;
                                lower = x->var_split_rule;
                        }
                }
        }
}


void Node::displayCurrNode(){

                std::cout << "Node address: " << this << std::endl;
                std::cout << "Node parent: " << parent << std::endl;

                std::cout << "Cur Node is leaf: " << isLeaf << std::endl;
                std::cout << "Cur Node is root: " << isRoot << std::endl;
                std::cout << "Cur The split_var is: " << var_split << std::endl;
                std::cout << "Cur The split_var_rule is: " << var_split_rule << std::endl;

                return;
}


void Node::deletingLeaves(){

     // Should I create some warn to avoid memoery leak
     //something like it will only delete from a nog?
     // Deleting
     delete left; // This release the memory from the left point
     delete right; // This release the memory from the right point
     left = this;  // The new pointer for the left become the node itself
     right = this; // The new pointer for the right become the node itself
     isLeaf = true;

     return;

}
// Getting the leaves (this is the function that gonna do the recursion the
//                      function below is the one that gonna initialise it)
void get_leaves(Node* x,  std::vector<Node*> &leaves_vec) {

        if(x->isLeaf){
                leaves_vec.push_back(x);
        } else {
                get_leaves(x->left, leaves_vec);
                get_leaves(x->right,leaves_vec);
        }

        return;

}



// Initialising a vector of nodes in a standard way
std::vector<Node*> leaves(Node* x) {
        std::vector<Node*> leaves_init(0); // Initialising a vector of a vector of pointers of nodes of size zero
        get_leaves(x,leaves_init);
        return(leaves_init);
}

// Sweeping the trees looking for nogs
void get_nogs(std::vector<Node*>& nogs, Node* node){
        if(!node->isLeaf){
                bool bool_left_is_leaf = node->left->isLeaf;
                bool bool_right_is_leaf = node->right->isLeaf;

                // Checking if the current one is a NOGs
                if(bool_left_is_leaf && bool_right_is_leaf){
                        nogs.push_back(node);
                } else { // Keep looking for other NOGs
                        get_nogs(nogs, node->left);
                        get_nogs(nogs, node->right);
                }
        }
}

// Creating the vectors of nogs
std::vector<Node*> nogs(Node* tree){
        std::vector<Node*> nogs_init(0);
        get_nogs(nogs_init,tree);
        return nogs_init;
}



// Initializing the forest
Forest::Forest(modelParam& data){

        // Creatina vector of size of number of trees
        trees.resize(data.n_tree);
        for(int  i=0;i<data.n_tree;i++){
                // Creating the stump for each tree
                trees[i] = new Node(data);
                // Filling up each stump for each tree
                trees[i]->Stump(data);
        }
}

// Function to delete one tree
// Forest::~Forest(){
//         for(int  i=0;i<trees.size();i++){
//                 delete trees[i];
//         }
// }

// Selecting a random node
Node* sample_node(std::vector<Node*> leaves_){

        // Getting the number of leaves
        int n_leaves = leaves_.size();
        return(leaves_[std::rand()%n_leaves]);
}

// Grow a tree for a given rule
void grow(Node* tree, modelParam &data, arma::vec &curr_res){

        // Getting the number of terminal nodes
        std::vector<Node*> t_nodes = leaves(tree) ;
        std::vector<Node*> nog_nodes = nogs(tree);

        // Selecting one node to be sampled
        Node* g_node = sample_node(t_nodes);

        // Store all old quantities that will be used or not
        double old_lower = g_node->lower;
        double old_upper = g_node->upper;
        int old_var_split = g_node->var_split;
        double old_var_split_rule = g_node->var_split_rule;

        // Calculate current tree log likelihood
        double tree_log_like = 0;

        // Calculating the whole likelihood fo the tree
        for(int i = 0; i < t_nodes.size(); i++){
                // cout << "Error SplineNodeLogLike" << endl;
                t_nodes[i]->splineNodeLogLike(data, curr_res);
                tree_log_like = tree_log_like + t_nodes[i]->log_likelihood;
        }

        // cout << "LogLike Node ok Grow" << endl;

        // Adding the leaves
        g_node->addingLeaves(data);

        // Selecting the var
        g_node-> sampleSplitVar(data);
        // Updating the limits
        g_node->getLimits();


        // Selecting a rule
        g_node->var_split_rule = (g_node->upper-g_node->lower)*rand_unif()+g_node->lower;

        // Create an aux for the left and right index
        int train_left_counter = 0;
        int train_right_counter = 0;

        int test_left_counter = 0;
        int test_right_counter = 0;

        // Updating the left and the right nodes
        for(int i = 0;i<data.x_train.n_rows;i++){
                if(g_node -> train_index[i] == -1 ){
                        g_node->left->n_leaf = train_left_counter;
                        g_node->right->n_leaf = train_right_counter;
                        break;
                }
                if(data.x_train(g_node->train_index[i],g_node->var_split)<g_node->var_split_rule){
                        g_node->left->train_index[train_left_counter] = g_node->train_index[i];
                        train_left_counter++;
                } else {
                        g_node->right->train_index[train_right_counter] = g_node->train_index[i];
                        train_right_counter++;
                }

        }


        // Updating the left and right nodes for the
        for(int i = 0;i<data.x_test.n_rows; i++){
                if(g_node -> test_index[i] == -1){
                        g_node->left->n_leaf_test = test_left_counter;
                        g_node->right->n_leaf_test = test_right_counter;
                        break;
                }
                if(data.x_test(g_node->test_index[i],g_node->var_split)<g_node->var_split_rule){
                        g_node->left->test_index[test_left_counter] = g_node->test_index[i];
                        test_left_counter++;
                } else {
                        g_node->right->test_index[test_right_counter] = g_node->test_index[i];
                        test_right_counter++;
                }
        }

        // If is a root node
        if(g_node->isRoot){
                g_node->left->n_leaf = train_left_counter;
                g_node->right->n_leaf = train_right_counter;
                g_node->left->n_leaf_test = test_left_counter;
                g_node->right->n_leaf_test = test_right_counter;
        }

        // Avoiding nodes lower than the node_min
        if((g_node->left->n_leaf<5) || (g_node->right->n_leaf<5) ){

                // cout << " NODES" << endl;
                // Returning to the old values
                g_node->var_split = old_var_split;
                g_node->var_split_rule = old_var_split_rule;
                g_node->lower = old_lower;
                g_node->upper = old_upper;
                g_node->deletingLeaves();
                return;
        }


        // Updating the loglikelihood for those terminal nodes
        g_node->left->splineNodeLogLike(data, curr_res);
        g_node->right->splineNodeLogLike(data, curr_res);


        // Calculating the prior term for the grow
        double tree_prior = log(data.alpha*pow((1+g_node->depth),-data.beta)) +
                log(1-data.alpha*pow((1+g_node->depth+1),-data.beta)) + // Prior of left node being terminal
                log(1-data.alpha*pow((1+g_node->depth+1),-data.beta)) - // Prior of the right noide being terminal
                log(1-data.alpha*pow((1+g_node->depth),-data.beta)); // Old current node being terminal

        // Getting the transition probability
        double log_transition_prob = log((0.3)/(nog_nodes.size()+1)) - log(0.3/t_nodes.size()); // 0.3 and 0.3 are the prob of Prune and Grow, respectively

        // Calculating the loglikelihood for the new branches
        double new_tree_log_like = tree_log_like - g_node->log_likelihood + g_node->left->log_likelihood + g_node->right->log_likelihood ;

        // Calculating the acceptance ratio
        double acceptance = exp(new_tree_log_like - tree_log_like + log_transition_prob + tree_prior);

        // Keeping the new tree or not
        if(rand_unif()<-acceptance){
                // Do nothing just keep the new tree
        } else {
                // Returning to the old values
                g_node->var_split = old_var_split;
                g_node->var_split_rule = old_var_split_rule;
                g_node->lower = old_lower;
                g_node->upper = old_upper;
                g_node->deletingLeaves();
        }

        return;

}


// Pruning a tree
void prune(Node* tree, modelParam&data, arma::vec &curr_res){


        // Getting the number of terminal nodes
        std::vector<Node*> t_nodes = leaves(tree);

        // Can't prune a root
        if(t_nodes.size()==1){
                // cout << "Nodes size " << t_nodes.size() <<endl;
                t_nodes[0]->splineNodeLogLike(data, curr_res);
                return;
        }

        std::vector<Node*> nog_nodes = nogs(tree);

        // Selecting one node to be sampled
        Node* p_node = sample_node(nog_nodes);


        // Calculate current tree log likelihood
        double tree_log_like = 0;

        // Calculating the whole likelihood fo the tree
        for(int i = 0; i < t_nodes.size(); i++){
                t_nodes[i]->splineNodeLogLike(data, curr_res);
                tree_log_like = tree_log_like + t_nodes[i]->log_likelihood;
        }

        // cout << "Error C1" << endl;
        // Updating the loglikelihood of the selected pruned node
        p_node->splineNodeLogLike(data, curr_res);
        // cout << "Error C2" << endl;

        // Getting the loglikelihood of the new tree
        double new_tree_log_like = tree_log_like + p_node->log_likelihood - (p_node->left->log_likelihood + p_node->right->log_likelihood);

        // Calculating the transition loglikelihood
        double transition_loglike = log((0.3)/(t_nodes.size())) - log((0.3)/(nog_nodes.size()));

        // Calculating the prior term for the grow
        double tree_prior = log(1-data.alpha*pow((1+p_node->depth),-data.beta))-
                log(data.alpha*pow((1+p_node->depth),-data.beta)) -
                log(1-data.alpha*pow((1+p_node->depth+1),-data.beta)) - // Prior of left node being terminal
                log(1-data.alpha*pow((1+p_node->depth+1),-data.beta));  // Prior of the right noide being terminal
                 // Old current node being terminal


        // Calculating the acceptance
        double acceptance = exp(new_tree_log_like - tree_log_like + transition_loglike + tree_prior);

        if(rand_unif()<acceptance){
                p_node->deletingLeaves();
        } else {
                // p_node->left->splineNodeLogLike(data, curr_res);
                // p_node->right->splineNodeLogLike(data, curr_res);
        }

        return;
}


// // Creating the change verb
void change(Node* tree, modelParam &data, arma::vec &curr_res){


        // Getting the number of terminal nodes
        std::vector<Node*> t_nodes = leaves(tree) ;
        std::vector<Node*> nog_nodes = nogs(tree);

        // Selecting one node to be sampled
        Node* c_node = sample_node(nog_nodes);

        // Calculate current tree log likelihood
        double tree_log_like = 0;

        if(c_node->isRoot){
                // cout << " THAT NEVER HAPPENS" << endl;
               c_node-> n_leaf = data.x_train.n_rows;
               c_node-> n_leaf_test = data.x_test.n_rows;
        }

        // cout << " Change error on terminal nodes" << endl;
        // Calculating the whole likelihood fo the tree
        for(int i = 0; i < t_nodes.size(); i++){
                // cout << "Loglike error " << ed
                t_nodes[i]->splineNodeLogLike(data, curr_res);
                tree_log_like = tree_log_like + t_nodes[i]->log_likelihood;
        }
        // cout << " Other kind of error" << endl;
        // If the current node has size zero there is no point of change its rule
        if(c_node->n_leaf==0) {
                return;
        }

        // Storing all the old loglikelihood from left
        double old_left_log_like = c_node->left->log_likelihood;
        arma::cube old_left_z = c_node->left->Z;
        arma::cube old_left_z_t = c_node->left->Z_t;
        arma::cube old_left_z_test = c_node->left->Z_test;
        arma::mat old_left_z_t_ones = c_node->left->z_t_ones;
        arma::vec old_left_leaf_res = c_node->left->leaf_res;



        arma::vec old_left_train_index = c_node->left->train_index;
        c_node->left->train_index.fill(-1); // Returning to the original
        int old_left_n_leaf = c_node->left->n_leaf;


        // Storing all of the old loglikelihood from right;
        double old_right_log_like = c_node->right->log_likelihood;
        arma::cube old_right_z = c_node->right->Z;
        arma::cube old_right_z_t = c_node->right->Z_t;
        arma::cube old_right_z_test = c_node->right->Z_test;
        arma::mat old_right_z_t_ones = c_node->right->z_t_ones;
        arma::vec old_right_leaf_res = c_node->right->leaf_res;


        arma::vec old_right_train_index = c_node->right->train_index;
        c_node->right->train_index.fill(-1);
        int old_right_n_leaf = c_node->right->n_leaf;



        // Storing test observations
        arma::vec old_left_test_index = c_node->left->test_index;
        arma::vec old_right_test_index = c_node->right->test_index;
        c_node->left->test_index.fill(-1);
        c_node->right->test_index.fill(-1);

        int old_left_n_leaf_test = c_node->left->n_leaf_test;
        int old_right_n_leaf_test = c_node->right->n_leaf_test;


        // Storing the old ones
        int old_var_split = c_node->var_split;
        int old_var_split_rule = c_node->var_split_rule;
        int old_lower = c_node->lower;
        int old_upper = c_node->upper;

        // Selecting the var
        c_node-> sampleSplitVar(data);
        // Updating the limits
        c_node->getLimits();
        // Selecting a rule
        c_node -> var_split_rule = (c_node->upper-c_node->lower)*rand_unif()+c_node->lower;
        // c_node -> var_split_rule = 0.0;

        // Create an aux for the left and right index
        int train_left_counter = 0;
        int train_right_counter = 0;

        int test_left_counter = 0;
        int test_right_counter = 0;


        // Updating the left and the right nodes
        for(int i = 0;i<data.x_train.n_rows;i++){
                // cout << " Train indexeses " << c_node -> train_index[i] << endl ;
                if(c_node -> train_index[i] == -1){
                        c_node->left->n_leaf = train_left_counter;
                        c_node->right->n_leaf = train_right_counter;
                        break;
                }
                // cout << " Current train index " << c_node->train_index[i] << endl;

                if(data.x_train(c_node->train_index[i],c_node->var_split)<c_node->var_split_rule){
                        c_node->left->train_index[train_left_counter] = c_node->train_index[i];
                        train_left_counter++;
                } else {
                        c_node->right->train_index[train_right_counter] = c_node->train_index[i];
                        train_right_counter++;
                }
        }



        // Updating the left and the right nodes
        for(int i = 0;i<data.x_test.n_rows;i++){

                if(c_node -> test_index[i] == -1){
                        c_node->left->n_leaf_test = test_left_counter;
                        c_node->right->n_leaf_test = test_right_counter;
                        break;
                }

                if(data.x_test(c_node->test_index[i],c_node->var_split)<c_node->var_split_rule){
                        c_node->left->test_index[test_left_counter] = c_node->test_index[i];
                        test_left_counter++;
                } else {
                        c_node->right->test_index[test_right_counter] = c_node->test_index[i];
                        test_right_counter++;
                }
        }

        // If is a root node
        if(c_node->isRoot){
                c_node->left->n_leaf = train_left_counter;
                c_node->right->n_leaf = train_right_counter;
                c_node->left->n_leaf_test = test_left_counter;
                c_node->right->n_leaf_test = test_right_counter;
        }


        if((c_node->left->n_leaf<5) || (c_node->right->n_leaf)<5){

                // Returning to the previous values
                c_node->var_split = old_var_split;
                c_node->var_split_rule = old_var_split_rule;
                c_node->lower = old_lower;
                c_node->upper = old_upper;

                // Returning to the old ones
                c_node->left->Z = old_left_z;
                c_node->left->Z_t = old_left_z_t;
                c_node->left->Z_test = old_left_z_test;
                c_node->left->z_t_ones = old_left_z_t_ones;
                c_node->left->leaf_res = old_left_leaf_res;

                c_node->left->n_leaf = old_left_n_leaf;
                c_node->left->n_leaf_test = old_left_n_leaf_test;
                c_node->left->log_likelihood = old_left_log_like;
                c_node->left->train_index = old_left_train_index;
                c_node->left->test_index = old_left_test_index;

                // Returning to the old ones
                c_node->right->Z = old_right_z;
                c_node->right->Z_t = old_right_z_t;
                c_node->right->Z_test = old_right_z_test;
                c_node->right->z_t_ones = old_right_z_t_ones;
                c_node->right->leaf_res = old_right_leaf_res;

                c_node->right->n_leaf = old_right_n_leaf;
                c_node->right->n_leaf_test = old_right_n_leaf_test;
                c_node->right->log_likelihood = old_right_log_like;
                c_node->right->train_index = old_right_train_index;
                c_node->right->test_index = old_right_test_index;

                return;
        }

        // Updating the new left and right loglikelihoods
        c_node->left->splineNodeLogLike(data,curr_res);
        c_node->right->splineNodeLogLike(data,curr_res);

        // Calculating the acceptance
        double new_tree_log_like =  - old_left_log_like - old_right_log_like + c_node->left->log_likelihood + c_node->right->log_likelihood;

        double acceptance = exp(new_tree_log_like);

        if(rand_unif()<acceptance){
                // Keep all the treesi
        } else {

                // Returning to the previous values
                c_node->var_split = old_var_split;
                c_node->var_split_rule = old_var_split_rule;
                c_node->lower = old_lower;
                c_node->upper = old_upper;

                // Returning to the old ones
                c_node->left->Z = old_left_z;
                c_node->left->Z_t = old_left_z_t;
                c_node->left->Z_test = old_left_z_test;
                c_node->left->z_t_ones = old_left_z_t_ones;
                c_node->left->leaf_res = old_left_leaf_res;


                c_node->left->n_leaf = old_left_n_leaf;
                c_node->left->n_leaf_test = old_left_n_leaf_test;
                c_node->left->log_likelihood = old_left_log_like;
                c_node->left->train_index = old_left_train_index;
                c_node->left->test_index = old_left_test_index;

                // Returning to the old ones
                c_node->right->Z = old_right_z;
                c_node->right->Z_t = old_right_z_t;
                c_node->right->Z_test = old_right_z_test;
                c_node->right->z_t_ones = old_right_z_t_ones;
                c_node->right->leaf_res = old_right_leaf_res;

                c_node->right->n_leaf = old_right_n_leaf;
                c_node->right->n_leaf_test = old_right_n_leaf_test;
                c_node->right->log_likelihood = old_right_log_like;
                c_node->right->train_index = old_right_train_index;
                c_node->right->test_index = old_right_test_index;

        }

        return;
}

// // Creating the change verb
// void OLD_change(Node* tree, modelParam &data, arma::vec &curr_res){
//
//
//         // Setting the size of the tree
//         // if(tree->isRoot) {
//         //         return;
//         // }
//
//         // Getting the number of terminal nodes
//         std::vector<Node*> t_nodes = leaves(tree) ;
//         std::vector<Node*> nog_nodes = nogs(tree);
//
//         // Selecting one node to be sampled
//         Node* c_node = sample_node(nog_nodes);
//
//         // Calculate current tree log likelihood
//         double tree_log_like = 0;
//
//         // Calculating the whole likelihood fo the tree
//         for(int i = 0; i < t_nodes.size(); i++){
//                 // cout << "Loglike error " << ed
//                 t_nodes[i]->splineNodeLogLike(data, curr_res);
//                 tree_log_like = tree_log_like + t_nodes[i]->log_likelihood;
//         }
//
//         cout << "Error C1" << endl;
//
//         // If the current node has size zero there is no point of change its rule
//         if(c_node->n_leaf==0) {
//                 cout << "Error C2" << endl;
//                 return;
//         }
//
//
//
//         // Storing all the old loglikelihood from left
//         double old_left_log_like = c_node->left->log_likelihood;
//         double old_left_r_sum = c_node->left->r_sum;
//         double old_left_r_sq_sum = c_node->left->r_sq_sum;
//         double old_left_n_leaf = c_node->left->n_leaf;
//         int old_left_train_index[data.x_train.n_rows];
//
//         // cout <<" Old left train" << endl;
//
//         for(int i = 0; i < data.x_train.n_rows;i++){
//                 old_left_train_index[i] = c_node->left->train_index[i];
//                 c_node->left->train_index[i] = -1;
//                 // cout << c_node->left->train_index[i] << endl;
//         }
//
//         // Storing all of the old loglikelihood from right;
//         double old_right_log_like = c_node->right->log_likelihood;
//         double old_right_r_sum = c_node->right->r_sum;
//         double old_right_r_sq_sum = c_node->right->r_sq_sum;
//         double old_right_n_leaf = c_node->right->n_leaf;
//         int old_right_train_index[data.x_train.n_rows];
//
//         for(int i = 0; i< data.x_train.n_rows;i++){
//                 old_right_train_index[i] = c_node->right->train_index[i];
//                 c_node->right->train_index[i] = -1;
//         }
//
//
//         // Storing test observations
//         int old_left_test_index[data.x_test.n_rows];
//         int old_right_test_index[data.x_test.n_rows];
//
//         for(int i = 0; i< data.x_test.n_rows;i++){
//                 old_left_test_index[i] = c_node->left->test_index[i];
//                 c_node->left->test_index[i] = -1;
//         }
//
//         for(int i = 0; i< data.x_test.n_rows;i++){
//                 old_right_test_index[i] = c_node->right->test_index[i];
//                 c_node->right->test_index[i] = -1;
//         }
//
//         // Storing the old ones
//         int old_var_split = c_node->var_split;
//         int old_var_split_rule = c_node->var_split_rule;
//         int old_lower = c_node->lower;
//         int old_upper = c_node->upper;
//
//         // Selecting the var
//         c_node-> sampleSplitVar(data);
//         // Updating the limits
//         c_node->getLimits();
//         // Selecting a rule
//         c_node -> var_split_rule = (c_node->upper-c_node->lower)*rand_unif()+c_node->lower;
//         // c_node -> var_split_rule = 0.0;
//
//         // Create an aux for the left and right index
//         int train_left_counter = 0;
//         int train_right_counter = 0;
//
//         int test_left_counter = 0;
//         int test_right_counter = 0;
//
//
//         // Updating the left and the right nodes
//         for(int i = 0;i<data.x_train.n_rows;i++){
//                 // cout << " Train indexeses " << c_node -> train_index[i] << endl ;
//                 if(c_node -> train_index[i] == -1){
//                         c_node->left->n_leaf = train_left_counter;
//                         c_node->right->n_leaf = train_right_counter;
//                         break;
//                 }
//                 // cout << " Current train index " << c_node->train_index[i] << endl;
//
//                 if(data.x_train(c_node->train_index[i],c_node->var_split)<c_node->var_split_rule){
//                         c_node->left->train_index[train_left_counter] = c_node->train_index[i];
//                         train_left_counter++;
//                 } else {
//                         c_node->right->train_index[train_right_counter] = c_node->train_index[i];
//                         train_right_counter++;
//                 }
//         }
//
//
//
//         // Updating the left and the right nodes
//         for(int i = 0;i<data.x_train.n_rows;i++){
//                 if(c_node -> test_index[i] == -1){
//                         c_node->left->n_leaf_test = test_left_counter;
//                         c_node->right->n_leaf_test = test_right_counter;
//                         break;
//                 }
//                 if(data.x_test(c_node->test_index[i],c_node->var_split)<c_node->var_split_rule){
//                         c_node->left->test_index[test_left_counter] = c_node->test_index[i];
//                         test_left_counter++;
//                 } else {
//                         c_node->right->test_index[test_right_counter] = c_node->test_index[i];
//                         test_right_counter++;
//                 }
//         }
//
//         // If is a root node
//         if(c_node->isRoot){
//                 c_node->left->n_leaf = train_left_counter;
//                 c_node->right->n_leaf = train_right_counter;
//                 c_node->left->n_leaf_test = train_left_counter;
//                 c_node->right->n_leaf_test = test_left_counter;
//         }
//
//         cout << "Error C3" << endl;
//         // cout << c_node->left->n_leaf + c_node->right->n_leaf << endl;
//         // Updating the new left and right loglikelihoods
//         c_node->left->splineNodeLogLike(data,curr_res);
//         c_node->right->splineNodeLogLike(data,curr_res);
//
//         cout << "Error C4" << endl;
//
//         // Calculating the acceptance
//         double new_tree_log_like =  - old_left_log_like - old_right_log_like + c_node->left->log_likelihood + c_node->right->log_likelihood;
//
//         double acceptance = exp(new_tree_log_like);
//
//         if(rand_unif()<acceptance){
//                 // Keep all the treesi
//         } else {
//
//                 // Returning to the previous values
//                 c_node->var_split = old_var_split;
//                 c_node->var_split_rule = old_var_split_rule;
//                 c_node->lower = old_lower;
//                 c_node->upper = old_upper;
//
//                 // Returning to the old ones
//                 c_node->left->r_sum = old_left_r_sum;
//                 c_node->left->r_sq_sum = old_left_r_sq_sum;
//                 c_node->left->n_leaf = old_left_n_leaf;
//                 c_node->left->log_likelihood = old_left_log_like;
//                 c_node->left->train_index = old_left_train_index;
//                 c_node->left->test_index = old_left_test_index;
//
//                 c_node->right->r_sum = old_right_r_sum;
//                 c_node->right->r_sq_sum = old_right_r_sq_sum;
//                 c_node->right->n_leaf = old_right_n_leaf;
//                 c_node->right->log_likelihood = old_right_log_like;
//                 c_node->right->train_index = old_right_train_index;
//                 c_node->right->test_index = old_right_test_index;
//
//         }
//
//         return;
// }




// Calculating the Loglilelihood of a node
void Node::splineNodeLogLike(modelParam& data, arma::vec &curr_res){

        // Getting number of leaves in case of a root
        if(isRoot){
                // Updating the r_sum
                r_sum = 0;
                n_leaf = data.x_train.n_rows;
                n_leaf_test = data.x_test.n_rows;
        }


        // Case of an empty node
        if(train_index[0]==-1){
        // if(n_leaf < 100){
                r_sum = 0;
                r_sq_sum = 10000;
                n_leaf = 0;
                log_likelihood = -2000000; // Absurd value avoid this case
                return;
        }

        // Creating the B spline
        arma::cube Z_(n_leaf,data.p,data.d_var);
        arma::cube Z_t_(data.p,n_leaf,data.d_var);

        arma::cube Z_test_(n_leaf_test,data.p,data.d_var);
        arma::mat z_t_ones_(data.p,data.d_var);
        arma::vec leaf_res_ = arma::vec(n_leaf);

        // Some aux elements
        arma::mat ones_vec(n_leaf,1,arma::fill::ones);
        arma::mat diag_aux(n_leaf,n_leaf,arma::fill::eye);
        arma::mat res_cov(n_leaf,n_leaf,arma::fill::zeros);
        s_tau_beta_0 = (n_leaf + data.tau_b_intercept/data.tau);

        // cout << "Z dimensions are: " << Z.n_rows << " " << Z.n_cols << " "<< Z.n_slices << endl;
        // cout << "Z test dimensions are: " << Z_test.n_rows << " " << Z_test.n_cols << " "<< Z_test.n_slices << endl;


        // Need to iterate a d-level
        for(int k = 0; k < data.d_var;k++){
                // cout << "Error on Z_t" << endl;

                // Train elements
                for(int i = 0; i < n_leaf;i++){
                        for(int j = 0 ; j < data.Z_train.n_cols; j++){
                                Z_(i,j,k) = data.Z_train(train_index[i],j,k);
                        }
                        leaf_res_(i) = curr_res(train_index[i]);
                }

                // cout << "Error on Z" << endl;
                Z_t_.slice(k) = Z_.slice(k).t();
                // cout << "Error on Z_ones" << endl;

                z_t_ones_.col(k) = Z_t_.slice(k)*ones_vec; // Col-sums from B - gonna use this again to sample beta_0 (remember is a row vector)
                // cout << "Error on res_cov" << endl;

                res_cov = res_cov + (1/data.tau_b(k))*Z_.slice(k)*Z_t_.slice(k);

                // Test elements
                for(int i = 0 ; i < n_leaf_test;i++){
                        for(int j = 0 ; j < data.Z_test.n_cols; j++){
                                Z_test_(i,j,k) = data.Z_test(test_index[i],j,k);
                        }
                }

        }


        // Sotring the new quantities for the node
        Z = Z_;
        Z_t = Z_t_;
        z_t_ones = z_t_ones_;
        Z_test = Z_test_;
        leaf_res = leaf_res_;

        // Adding the remaining quantities
        res_cov  = ((1/data.tau)*diag_aux + (1/data.tau_b_intercept) + res_cov);


        // If is smaller then the node size still need to update theq quantities;
        // if(n_leaf < 15){
        //         log_likelihood = -2000000; // Absurd value avoid this case
        //         return;
        // }

        // Getting the log-likelihood;
        log_likelihood = log_dmvn(leaf_res,res_cov);

        return;

}


// Update betas
void updateBeta(Node* tree, modelParam &data){

        // Getting the terminal nodes
        std::vector<Node*> t_nodes = leaves(tree);
        arma::mat aux_diag(data.p,data.p,arma::fill::eye);

        // Iterating over the terminal nodes and updating the beta values
        for(int i = 0; i < t_nodes.size();i++){
                if(t_nodes[i]->n_leaf==0 ){
                        /// Skip nodes that doesn't have any observation within terminal node
                        cout << " SKIPPED" << endl;
                        continue;
                }

                // Calculating a cube with all {Z^(j)}\beta
                // cout << "error in Z_j_beta" << endl;
                arma::mat Z_j_beta(t_nodes[i]->n_leaf,data.d_var);
                // cout << "Dimensions of Z are: " << t_nodes[i]->Z.n_rows << " " <<  t_nodes[i]->Z.n_cols << " " <<  t_nodes[i]->Z.n_slices << endl;
                for(int k = 0; k<data.d_var;k++){
                        Z_j_beta.col(k) = t_nodes[i]->Z.slice(k)*t_nodes[i]->betas.col(k);
                }

                // cout << "error in Beta sampling" << endl;

                // Iterating ove each predictor
                for(int j=0;j<data.d_var;j++){
                        // Calculating elements exclusive to each predictor
                        arma::mat aux_precision_inv = arma::inv_sympd(t_nodes[i]->Z_t.slice(j)*t_nodes[i]->Z.slice(j)+(data.tau_b(j)/data.tau)*aux_diag);
                        arma::mat beta_mean = aux_precision_inv*(t_nodes[i]->Z_t.slice(j)*t_nodes[i]->leaf_res-t_nodes[i]->Z_t.slice(j)*(t_nodes[i]->beta_zero+sum_exclude_col(Z_j_beta,j)));
                        arma::mat beta_cov = (1/data.tau)*aux_precision_inv;

                        // cout << "Error sample BETA" << endl;
                        arma::mat sample = arma::randn<arma::mat>(beta_cov.n_cols);
                        // cout << "Error variance" << endl;
                        t_nodes[i]->betas.col(j) = arma::chol(beta_cov,"lower")*sample + beta_mean;
                }
        }
}


// Update betas
void updateGamma(Node* tree, modelParam &data){

        // Getting the terminal nodes
        std::vector<Node*> t_nodes = leaves(tree);

        // Iterating over the terminal nodes and updating the beta values
        for(int i = 0; i < t_nodes.size();i++){

                if(t_nodes[i]->n_leaf==0 ){
                        /// Skip nodes that doesn't have any observation within terminal node
                        cout << " SKIPPED" << endl;
                        continue;
                }
                // cout << "Error mean" << endl;
                double s_gamma = t_nodes[i]->n_leaf+(data.tau_b_intercept/data.tau);
                double sum_beta_z_one = 0;



                for(int j = 0; j<data.d_var;j++){
                        sum_beta_z_one = sum_beta_z_one + arma::as_scalar(t_nodes[i]->betas.col(j).t()*(t_nodes[i]->z_t_ones.col(j)));
                }

                // cout << "Beta_zero_mean " << endl;
                double gamma_mean = ((1/s_gamma)*(accu(t_nodes[i]->leaf_res)-sum_beta_z_one));
                // cout << "Error sample" << endl;
                double gamma_sd = sqrt(1/(data.tau*s_gamma));
                // cout << "Error variance" << endl;
                t_nodes[i]->beta_zero = arma::randn()*gamma_sd+gamma_mean;
        }
}

// Get the prediction
void getPredictions(Node* tree,
                    modelParam &data,
                    arma::vec& current_prediction_train,
                    arma::vec& current_prediction_test){

        // Getting the current prediction
        vector<Node*> t_nodes = leaves(tree);
        for(int i = 0; i<t_nodes.size();i++){

                // Skipping empty nodes
                if(t_nodes[i]->n_leaf==0){
                        cout << " THERE ARE EMPTY NODES" << endl;
                        continue;
                }

                // Calculating the sum over multiple predictors
                arma::vec betas_b_sum(t_nodes[i]->n_leaf,arma::fill::zeros);
                // cout << "Dimensions of Z are " << t_nodes[i] -> Z.n_rows << " " << t_nodes[i]->Z.n_cols << " " << t_nodes[i]->Z.n_slices << endl;

                for(int j = 0;j<data.d_var;j++){
                        betas_b_sum = betas_b_sum + t_nodes[i]->Z.slice(j)*t_nodes[i]->betas.col(j);
                }

                // Getting the vector of prediction for the betas and b for this node
                arma::vec leaf_y_hat = t_nodes[i]->beta_zero + betas_b_sum ;

                // Message to check if the dimensions are correct;
                if(leaf_y_hat.size()!=t_nodes[i]->n_leaf){
                        cout << " Pay attention something is wrong here" << endl;
                }

                // For the training samples
                for(int j = 0; j<data.x_train.n_rows; j++){


                        if((t_nodes[i]->train_index[j])==-1.0){
                                break;
                        }
                        current_prediction_train[t_nodes[i]->train_index[j]] = leaf_y_hat[j];
                }

                if(t_nodes[i]->n_leaf_test == 0 ){
                        continue;
                }


                // Calculating the sum over multiple predictors
                arma::vec betas_b_sum_test(t_nodes[i]->Z_test.n_rows,arma::fill::zeros);

                for(int j = 0;j<data.d_var;j++){
                        betas_b_sum_test = betas_b_sum + t_nodes[i]->Z_test.slice(j)*t_nodes[i]->betas.col(j);
                }


                // Creating the y_hattest
                arma::vec leaf_y_hat_test = t_nodes[i]->beta_zero + betas_b_sum_test;


                // Regarding the test samples
                for(int j = 0; j< data.x_test.n_rows;j++){

                        if(t_nodes[i]->test_index[j]==-1){
                                break;
                        }

                        current_prediction_test[t_nodes[i]->test_index[j]] = leaf_y_hat_test[j];

                }

        }
}

// Updating the tau parameter
void updateTau(arma::vec &y_hat,
               modelParam &data){

        // Getting the sum of residuals square
        double tau_res_sq_sum = dot((y_hat-data.y),(y_hat-data.y));

        data.tau = R::rgamma((0.5*data.y.size()+data.a_tau),1/(0.5*tau_res_sq_sum+data.d_tau));

        return;
}

// Updating tau b parameter
void updateTauB(Forest all_trees,
                modelParam &data){


        arma::vec beta_count_total(data.d_var,arma::fill::zeros);
        arma::vec beta_sq_sum_total(data.d_var,arma::fill::zeros);

        for(int t = 0; t< all_trees.trees.size();t++){

                Node* tree = all_trees.trees[t];

                // Getting tau_b
                vector<Node*> t_nodes = leaves(tree);


                // Iterating over terminal nodes
                for(int i = 0; i< t_nodes.size(); i++ ){

                        // Simple error test
                        if(t_nodes[i]->betas.size()<1) {
                                continue;
                        }

                        for(int j = 0; j < data.d_var;j++){
                                beta_sq_sum_total(j) = beta_sq_sum_total(j) + arma::accu(arma::square(t_nodes[i]->betas.col(j)));
                                beta_count_total(j) = beta_count_total(j) + data.p;
                        }
                }

        }

        for(int j = 0; j < data.d_var; j++){
                data.tau_b(j) = R::rgamma((0.5*beta_count_total(j) + 0.5*data.nu),1/(0.5*beta_sq_sum_total(j)+0.5*data.delta(j)*data.nu));
        }
        return;

}

// Updating the posterior for Delta;
void updateDelta(modelParam &data){
        // Iterating over the "p" predictors of delta
        for(int j = 0; j<data.d_var;j++){
                data.delta(j) = R::rgamma((0.5*data.nu + data.a_delta),1/(0.5*data.nu*data.tau_b(j)+data.d_delta));
        }
        return;
}

// Updating tau b parameter
void updateTauBintercept(Forest all_trees,
                modelParam &data,
                double a_tau_b,
                double d_tau_b){


        double beta_count_total = 0.0;
        double beta_sq_sum_total = 0.0;

        for(int t = 0; t< all_trees.trees.size();t++){

                Node* tree = all_trees.trees[t];

                // Getting tau_b
                vector<Node*> t_nodes = leaves(tree);


                // Iterating over terminal nodes
                for(int i = 0; i< t_nodes.size(); i++ ){

                        if(t_nodes[i]->betas.size()<1) {
                                continue;
                        }

                        // Getting only the intercept
                        beta_sq_sum_total = beta_sq_sum_total + t_nodes[i]->betas[0]*t_nodes[i]->betas[0];
                        beta_count_total ++;
                }

        }

        data.tau_b_intercept = R::rgamma((0.5*beta_count_total + a_tau_b),1/(0.5*beta_sq_sum_total+d_tau_b));


        return;

}

// Creating the BART function
// [[Rcpp::export]]
Rcpp::List sbart(arma::mat x_train,
          arma::vec y_train,
          arma::mat x_test,
          arma::cube Z_train,
          arma::cube Z_test,
          int n_tree,
          int n_mcmc,
          int n_burn,
          double tau, double mu,
          double tau_mu, double tau_b, double tau_b_intercept,
          double alpha, double beta,
          double a_tau, double d_tau,
          double nu, double delta,
          double a_delta, double d_delta,
          double a_tau_b, double d_tau_b,
          arma::vec p_sample, arma::vec p_sample_levels){

        // Posterior counter
        int curr = 0;

        // cout << "Error ModelParam init" << endl;

        // Creating the structu object
        modelParam data(x_train,
                        y_train,
                        x_test,
                        Z_train,
                        Z_test,
                        n_tree,
                        alpha,
                        beta,
                        tau_mu,
                        tau_b,
                        tau_b_intercept,
                        tau,
                        a_tau,
                        d_tau,
                        nu,
                        delta,
                        a_delta,
                        d_delta,
                        n_mcmc,
                        n_burn,
                        p_sample,
                        p_sample_levels);

        // cout << "Error ModelParam ok" << endl;

        // Getting the n_post
        int n_post = n_mcmc - n_burn;

        // Defining those elements
        arma::mat y_train_hat_post = arma::zeros<arma::mat>(data.x_train.n_rows,n_post);
        arma::mat y_test_hat_post = arma::zeros<arma::mat>(data.x_test.n_rows,n_post);
        arma::cube all_tree_post(y_train.size(),n_tree,n_post,arma::fill::zeros);
        arma::vec tau_post = arma::zeros<arma::vec>(n_post);
        arma::mat tau_b_post = arma::mat(data.d_var,n_post);
        arma::vec tau_b_post_intercept = arma::zeros<arma::vec>(n_post);


        // Defining other variables
        arma::vec partial_pred = arma::zeros<arma::vec>(data.x_train.n_rows);
        // arma::vec partial_pred = (data.y*(n_tree-1))/n_tree;
        arma::vec partial_residuals = arma::zeros<arma::vec>(data.x_train.n_rows);
        arma::mat tree_fits_store = arma::zeros<arma::mat>(data.x_train.n_rows,data.n_tree);
        arma::mat tree_fits_store_test = arma::zeros<arma::mat>(data.x_test.n_rows,data.n_tree);

        // Getting zeros
        arma::vec prediction_train_sum = arma::zeros<arma::vec>(data.x_train.n_rows);
        arma::vec prediction_test_sum = arma::zeros<arma::vec>(data.x_test.n_rows);
        arma::vec prediction_train = arma::zeros<arma::vec>(data.x_train.n_rows);
        arma::vec prediction_test = arma::zeros<arma::vec>(data.x_test.n_rows);

        double verb;

        // Defining progress bars parameters
        const int width = 70;
        double pb = 0;


        // cout << " Error one " << endl;

        // Selecting the train
        Forest all_forest(data);

        for(int i = 0;i<data.n_mcmc;i++){

                // Initialising PB
                std::cout << "[";
                int k = 0;
                // Evaluating progress bar
                for(;k<=pb*width/data.n_mcmc;k++){
                        std::cout << "=";
                }

                for(; k < width;k++){
                        std:: cout << " ";
                }

                std::cout << "] " << std::setprecision(5) << (pb/data.n_mcmc)*100 << "%\r";
                std::cout.flush();


                prediction_train_sum = arma::zeros<arma::vec>(data.x_train.n_rows);
                prediction_test_sum = arma::zeros<arma::vec>(data.x_test.n_rows);

                for(int t = 0; t<data.n_tree;t++){

                        // Updating the partial residuals
                        partial_residuals = data.y-partial_pred+tree_fits_store.col(t);

                        // Iterating over all trees
                        verb = rand_unif();
                        if(all_forest.trees[t]->isLeaf & all_forest.trees[t]->isRoot){
                                verb = 0.27;
                        }

                        verb = 0.27;
                        // Selecting the verb
                        if(verb < 0.3){
                                // cout << " Grow error" << endl;
                                grow(all_forest.trees[t],data,partial_residuals);
                        } else if(verb>=0.3 & verb <0.6) {
                                // cout << " Prune error" << endl;
                                prune(all_forest.trees[t], data, partial_residuals);
                        } else {
                                // cout << " Change error" << endl;
                                change(all_forest.trees[t], data, partial_residuals);
                                // std::cout << "Error after change" << endl;
                        }


                        // Updating the all the parameters
                        // cout << "Error on Beta" << endl;
                        updateBeta(all_forest.trees[t], data);
                        // cout << "Error on Gamma" << endl;
                        updateGamma(all_forest.trees[t],data);
                        // Getting predictions
                        // cout << " Error on Get Predictions" << endl;
                        getPredictions(all_forest.trees[t],data,prediction_train,prediction_test);



                        // Updating the partial pred
                        partial_pred = partial_pred - tree_fits_store.col(t) + prediction_train;

                        tree_fits_store.col(t) = prediction_train;

                        prediction_train_sum = prediction_train_sum + prediction_train;

                        prediction_test_sum = prediction_test_sum + prediction_test;


                }

                // Updating the Tau
                // std::cout << "Error TauB: " << data.tau_b << endl;
                // updateTauB(all_forest,data);
                // std::cout << "Error Delta: " << data.delta << endl;
                updateDelta(data);
                // std::cout << "Error Tau: " << data.tau<< endl;
                updateTau(partial_pred, data);
                // std::cout << " All good " << endl;
                if(i >= n_burn){
                        // Storing the predictions
                        y_train_hat_post.col(curr) = prediction_train_sum;
                        y_test_hat_post.col(curr) = prediction_test_sum;
                        all_tree_post.slice(curr) = tree_fits_store;
                        tau_post(curr) = data.tau;
                        tau_b_post.col(curr) = data.tau_b;
                        tau_b_post_intercept(curr) = data.tau_b_intercept;
                        curr++;
                }

                pb += 1;

        }
        // Initialising PB
        std::cout << "[";
        int k = 0;
        // Evaluating progress bar
        for(;k<=pb*width/data.n_mcmc;k++){
                std::cout << "=";
        }

        for(; k < width;k++){
                std:: cout << " ";
        }

        std::cout << "] " << std::setprecision(5) << 100 << "%\r";
        std::cout.flush();

        std::cout << std::endl;

        return Rcpp::List::create(y_train_hat_post,
                                  y_test_hat_post,
                                  tau_post,
                                  all_tree_post,
                                  tau_b_post,
                                  tau_b_post_intercept);
}


//[[Rcpp::export]]
arma::mat mat_init(int n){
        arma::mat A(n,1,arma::fill::ones);
        return A + 4.0;
}


//[[Rcpp::export]]
arma::vec vec_init(int n){
        arma::vec A(n);
        return A;
}


// Comparing matrix inversions in armadillo
//[[Rcpp::export]]
arma::mat std_inv(arma::mat A, arma::vec diag){

        arma::mat diag_aux = arma::diagmat(diag);
        return arma::inv(A.t()*A+diag_aux);
}

//[[Rcpp::export]]
arma::mat std_pinv(arma::mat A, arma::vec diag){

        arma::mat diag_aux = arma::diagmat(diag);
        return arma::inv_sympd(A.t()*A+diag_aux);
}

//[[Rcpp::export]]
arma::mat faster_simple_std_inv(arma::mat A, arma::vec diag){
        arma::mat diag_aux = arma::diagmat(diag);
        arma::mat L = chol(A.t()*A+diag_aux,"lower");
        return arma::inv(L.t()*L);
}

//[[Rcpp::export]]
arma::mat faster_std_inv(arma::mat A, arma::vec diag){
        arma::mat ADinvAt = A.t()*arma::diagmat(1.0/diag)*A;
        arma::mat L = arma::chol(ADinvAt + arma::eye(ADinvAt.n_cols,ADinvAt.n_cols),"lower");
        arma::mat invsqrtDA = arma::solve(A.t()/arma::diagmat(arma::sqrt(diag)),L.t());
        arma::mat Ainv = invsqrtDA *invsqrtDA.t()/(ADinvAt + arma::eye(ADinvAt.n_cols,ADinvAt.n_cols));
        return Ainv;
}


//[[Rcpp::export]]
arma::vec rMVN2(const arma::vec& b, const arma::mat& Q)
{
        arma::mat Q_inv = arma::inv(Q);
        arma::mat U = arma::chol(Q_inv, "lower");
        arma::vec z= arma::randn<arma::mat>(Q.n_cols);

        return arma::solve(U.t(), arma::solve(U, z, arma::solve_opts::no_approx), arma::solve_opts::no_approx) + b;
}

//[[Rcpp::export]]
arma::vec rMVNslow(const arma::vec& b, const arma::mat& Q){

        // cout << "Error sample BETA" << endl;
        arma::vec sample = arma::randn<arma::mat>(Q.n_cols);
        return arma::chol(Q,"lower")*sample + b;

}

//[[Rcpp::export]]
arma::mat matrix_mat(arma::cube array){
        return array.slice(1).t()*array.slice(2);
}


