#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

double getDist(NumericVector x, NumericVector y){
  double dist = 0;
  
  dist = sqrt( pow(x[0]-y[0],2) + pow(x[1]-y[1],2));
  return dist;
}

double edgeProb(double a, double d, double L, double dist){
  return a * exp(d*dist/L);
}

// [[Rcpp::export]]
List waxman(int n,double a, double d) {
  NumericMatrix x(n,2);
  x(_,0) = runif(n,0,1);
  x(_,1) = runif(n,0,1);
  
  IntegerVector node_states(n);
  NumericVector p=runif(1,0,n);
  int p2=p(0);
  node_states(p2-1)=1;
  
  double dist = 0;
  
  double L = sqrt(2);
  
  int Nedges = n * (n-1) / 2;  
  NumericMatrix edges(Nedges,6);
  edges(_,0) = runif(Nedges,0,1);
  NumericMatrix A(n,n);
  
  int row = 0;
  
  for(int i = 0; i < n; i++){
    for(int j = (i+1); j < n; j++){
      dist =  getDist(x(i,_),x(j,_)); 
      edges(row,1) = dist;
      edges(row,2) = edgeProb(a,d,L,dist);
      if( edges(row,0) < edges(row,2)){
        edges(row,3) = 1;
        A(i,j)=1;
        A(j,i)=1;
      } else {
        edges(row,3) = 0;
      }
      edges(row,4) = i+1;
      edges(row,5) = j+1;
      row++;
    }
  }
  return List::create(_["nodes"] = x, _["edges"] = edges, _["A"]=A,_["node_states"] = node_states);
}

// [[Rcpp::export]]
NumericMatrix apply_cluster(NumericMatrix nodes, NumericMatrix F, double f){
  int N = nodes.nrow();
  int m = F.nrow();
  
  NumericMatrix adjusted_nodes(N,2);
  
  for(int i = 0; i < N; i++){
    double min_dist = 10;
    for(int j = 0; j < m; j++){
      double dist = getDist(nodes(i,_),F(j,_));
      if(dist < min_dist){
        min_dist = dist;
        double xn = F(j,0) + f * (F(j,0) - nodes(i,0));
        double yn = F(j,1) + f * (F(j,1) - nodes(i,1));
        adjusted_nodes(i,0) = xn;
        adjusted_nodes(i,1) = yn;
      }
    }
  }
  return adjusted_nodes;
}


// [[Rcpp::export]]
List keeling(int n, int F, double f,double a, double d){
    NumericMatrix tmp(n,2);
  tmp(_,0) = runif(n,0,1);
  tmp(_,1) = runif(n,0,1);
  
  IntegerVector node_states(n);
  NumericVector p=runif(1,0,n);
  int p2=p(0);
  node_states(p2-1)=1;
  
   NumericMatrix FLoc(F,2);
  FLoc(_,0) = runif(F,0,1);
  FLoc(_,1) = runif(F,0,1);
  
  NumericMatrix x=apply_cluster(tmp,FLoc,f);
  
  double dist = 0;
  
  double L = sqrt(2);
  
  int Nedges = n * (n-1) / 2;  
  NumericMatrix edges(Nedges,6);
  edges(_,0) = runif(Nedges,0,1);
  NumericMatrix A(n,n);
  
  int row = 0;
  
  for(int i = 0; i < n; i++){
    for(int j = (i+1); j < n; j++){
      dist =  getDist(x(i,_),x(j,_)); 
      edges(row,1) = dist;
      edges(row,2) = edgeProb(a,d,L,dist);
      if( edges(row,0) < edges(row,2)){
        edges(row,3) = 1;
        A(i,j)=1;
        A(j,i)=1;
      } else {
        edges(row,3) = 0;
      }
      edges(row,4) = i+1;
      edges(row,5) = j+1;
      row++;
    }
  }
  return List::create(_["nodes"] = x, _["edges"] = edges, _["A"]=A,_["node_states"] = node_states);
}

void print_edges(IntegerMatrix x){
  int n = x.nrow();
  Rprintf("From\tTo\tState\tState meaning\n");
  for(int i = 0; i < n; i++){
    Rprintf("%i\t%i\t%i\t",x(i,0),x(i,1),x(i,2));
    if(x(i,2)==0){
      Rprintf("Susceptible\n");
    } else if(x(i,2) == 1){
      Rprintf("Infected\n");
    } else {
      Rprintf("Recovered\n");
    }
  }
}

double get_time_to_event(int n_infect, int n_infect_edges, double gamma, double beta, int &state){
  double infect_rate = beta * n_infect_edges;
  double recover_rate = n_infect * gamma;
  double rate = infect_rate + recover_rate;
  
  NumericVector x = rexp(1,rate);
  NumericVector u = runif(1);
  
  if(u(0) < (infect_rate / rate)){
    state = 0; //infection
  } else {
    state = 1; // recovery
  }
  return x(0);
  
}

int get_n_infect_edges(IntegerMatrix edges){
  int n = edges.nrow();
  int count = 0;
  for(int i = 0; i < n; i++){
    if(edges(i,2)==1){
      count++;
    }
  }
  return count;
}

int select_new_infection(IntegerMatrix edges){
  int n = edges.nrow();
  NumericVector u = runif(n);
  int new_infect = 0;
  double min = 1;
  for(int i = 0; i < n ; i++){
    if(edges(i,2) == 1 & u(i) < min){
      min = u(i);
      new_infect = edges(i,1);
    }
  }
  return new_infect;
}

int select_recovery(IntegerVector nodes_state){
  int n = nodes_state.size();
  NumericVector u = runif(n);
  int new_recovery;
  double min = 1.1;
  for(int i = 0; i < n ; i++){
    if(nodes_state(i) == 1 & u(i) < min){
      min = u(i);
      new_recovery = i;
    }
  }
  return new_recovery;
}

void update_infect_edges(IntegerMatrix &edges, IntegerVector nodes_states){
  int n = edges.nrow();
  for(int i = 0; i < n; i++){
    int from_state = nodes_states(edges(i,0));
    int to_state = nodes_states(edges(i,1));
    if(from_state == 1 & to_state == 0){
      edges(i,2) = 1;
    } else if(to_state == 1 | to_state == 2){
      edges(i,2) = 2;
    } else {
      edges(i,2) = 0;
    }
  }
}

void print_states(IntegerVector x){
  int n = x.size();
  Rprintf("Node states:\nNode\tState\tState meaning\n");
  for(int i = 0; i < n; i++){
    Rprintf("%i\t%i\t",i,x(i));
    if(x(i)==0){
      Rprintf("Susceptible\n");
    } else if(x(i) == 1){
      Rprintf("Infected\n");
    } else {
      Rprintf("Recovered\n");
    }
  }
  Rprintf("\n");
}

int get_n_infected(IntegerVector node_states){
  int count = 0;
  for(int i = 0; i < node_states.size(); i++){
    if(node_states(i) == 1){
      count++;
    }
  }
  return count;
}

void print_times(NumericVector x){
  Rprintf("Node\tTime\n");
  for(int i = 0; i < x.size();i++){
    Rprintf("%i\t%4.2f\n",i,x(i));
  }
}
// [[Rcpp::export]]
IntegerVector sim_epidemic(IntegerMatrix edges, 
IntegerVector node_states, 
double gamma, 
double Beta,
int trace) {

if (sum(node_states==0)!=99){
  Rprintf("Error: node_states must have 99 susceptibles. \n");
  return node_states;
}
if (sum(node_states==1)!=1){
  Rprintf("Error: node_states must have 1 infective. \n");
  return node_states;
}

  // initialise
  double beta=Beta/(node_states.length()-1);
  // update edges
  update_infect_edges(edges, node_states);
  // get number of infected
  int n_infected = get_n_infected(node_states);
  // get number of infected edges
  int n_infected_edges = get_n_infect_edges(edges);
  double time_to_event;
  int state;
  int step = 1;
  int n_nodes = node_states.size();
  NumericVector time_to_infect(n_nodes);
  NumericVector time_to_recover(n_nodes);
  double time = 0;
  
  if(trace == 1){
    // print edges
    print_edges(edges);
    //print node states
    print_states(node_states);
  }
  
  while(n_infected > 0){
    if(trace == 1){Rprintf("Step: %i\n",step);}
    step++;
    
    // Get time to next event and find if infection 0 or recovery 1
    time_to_event = get_time_to_event(n_infected, n_infected_edges, gamma, beta, state);
    time = time + time_to_event;
    
    
    if(state == 0){
      //infection
      // Choose new infection
      int new_infect = select_new_infection(edges);
      // Update state vector
      node_states(new_infect) = 1;
      // printed new infected
      if(trace == 1){
        Rprintf("New infection is %i\n",new_infect);
        Rprintf("Infection\n");
      }
      // Add time 
      time_to_infect(new_infect) = time;
    } else {
      //recovery
      // Choose new recovery
      int new_recovery = select_recovery(node_states);
      // update nodes state
      node_states(new_recovery) = 2;
      // print recovery
      if(trace == 1){
        Rprintf("Recovery\n");
        Rprintf("New recovery is %i\n",new_recovery);
      }
      // add time
      time_to_recover(new_recovery) = time;
    }
    // update edges
    update_infect_edges(edges, node_states);
    // get number of infected
    n_infected = get_n_infected(node_states);
    // get number of infected edges
    n_infected_edges = get_n_infect_edges(edges);
    
    if(trace == 1){
      //print node states
      print_states(node_states);
      //print times
      Rprintf("Infection times\n");
      print_times(time_to_infect);
      //print times
      Rprintf("Recovery times\n");
      print_times(time_to_recover);
      // print edges
      print_edges(edges);
    }
  }
  
  return node_states;
}

// [[Rcpp::export]]
NumericVector finalSize(double Beta, int n) {
  double beta=Beta/(n-1);
  double p1;
  NumericVector q(n+1);
  q(1)=1;
  for (int Z2=-1; Z2<(n-2); Z2++){
    for (int Z1=(Z2+1); Z1<(n-1); Z1++){
      p1 = 1/(1+1/(beta*(n-(Z1+1))));
      q(Z1+2) = q(Z1+2)+q(Z1+1)*p1;
      q(Z1+1) = q(Z1+1)*(1-p1);
    }
  }
  return q;
}


// [[Rcpp::export]]
double Prod(NumericVector a){
  double total=1;
  for (int i=0; i<a.length(); i++){
    total *= a[i];
  }
  return total;
}

// [[Rcpp::export]]
int myWhich(NumericVector a, double b){
  int out;
  for (int i=0; i<a.length(); i++){
    if (a(i)==b){
      out=i;
    }
  }
  return out;
}
 
 // [[Rcpp::export]]
double mle(NumericVector Ne,NumericVector Beta){
  NumericVector lik(Beta.length());
  NumericVector pNe;
  NumericVector tp(Ne.length());
  for (int i=0; i<Beta.length(); i++){
    pNe = finalSize(Beta(i),100);
    for (int j=0; j<Ne.length(); j++){
      tp(j)=pNe(Ne(j));
    }
    lik(i) = Prod(tp);
  }
  return Beta(myWhich(lik,max(lik)));
}
