#include <Rcpp.h>
#include <queue>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::depends(igraph)]]

struct real_key_int {
  double key; //priority of this vertex
  int id; //index of the vertex in the vertex array
  
  operator int() {return id;}
  
  operator double() {return key;}
  
  friend bool operator<(const real_key_int &lhs, const real_key_int &rhs);
  
  friend bool operator>(const real_key_int &lhs, const real_key_int &rhs);
  
};

bool operator<(const real_key_int &lhs, const real_key_int &rhs) {
    return (lhs.key > rhs.key);
};
bool operator>(const real_key_int &lhs, const real_key_int &rhs) {
  return (lhs < rhs);
};

struct node {
  bool mark; //marks if this node has been added to ball
  double dist; //marks distance from root to this node
  int id, parent; //index of element in vertex array
  std::vector<int> out_edges; //indices of elements in edge list arrays
};

std::vector<int> upper_ids(NumericVector v1) {
    int upper_id, v1_size;
  double upper;
  std::vector<int> outlist;
  
  v1_size = v1.size();
  upper_id = which_max(v1); //This is the index of the max element
  upper = floor(v1(upper_id));
  for (int ii=0;ii<v1_size;ii++) {
    if (v1(ii) >= upper) outlist.push_back(ii);
  }
  return outlist;
}

// [[Rcpp::export]]
NumericMatrix growBall(NumericVector vl,NumericVector elSources,NumericVector 
  elTargets,NumericVector ew, int rv, double r,bool verbose=false,
  bool distances=false,bool parents=false) {
  //vl - vertex id list
  //vw - vertex weight list - ignored for now 
  //elSources - list of source vertices of edges
  //elTargets - list of target vertices of directed edges
  //ew - array of weights of each edge
  //rv - root vertex id
  //r - radius of ball
  //note: since the input lists are assumed to come from R
  //      they will begin at 1 instead of 0. Need to take care.
  std::vector<int> ballVec,pVec; //stores ids of vertices in ball and parent nodes
  std::vector<double> dVec; //store distances from nodes in ball to the root
  std::priority_queue< real_key_int, std::vector<real_key_int>, 
                      std::greater<real_key_int> > vQueue;
  NumericMatrix outmat;
  real_key_int current_key, next_key, last_key;
  node current_node, next_node;
  std::vector<node> vVec(vl.size());
  int vl_size, el_size, node_el_size, eId,b_size, outrows=1, ri;
  double key_val, inf_dist, temp_dist;
  bool finished;
  
  inf_dist=99999.0;
  vl_size = vl.size();
  el_size = elSources.size();
  
  for (int ii=0;ii<vl_size;ii++) {
    vVec[ii].id=ii;
    vVec[ii].dist=inf_dist;
    vVec[ii].mark=false;
    vVec[ii].parent=ii;
  }
  
  //we will handle vertex weights by adding them to
  //the weights of all outgoing edges
  for (int ii=0;ii<el_size;ii++) {
    //Handle vertex weights here... not yet implemented
    vVec[elSources[ii]-1].out_edges.push_back(ii);
  }
  
  //vVec[rv].dist = 0.0;
  current_key.id = rv-1;
  current_key.key = 0.0;
  vQueue.push(current_key);
  //finished=(vQueue.empty() || current_key.key > r);
  while (!vQueue.empty() ) { //&& current_key.key <= r) {
    last_key = current_key;
    current_key = vQueue.top();
    if (verbose) Rprintf("popped key for vertex %d from stack\n",current_key.id);
    if (verbose) Rprintf("\t key distance was %4.2f, vertex at key id had distance %f\n",
                          current_key.key,vVec[current_key.id].dist);
    vQueue.pop();
    if (current_key.key < vVec[current_key.id].dist && \
        current_key.key < inf_dist) {
      if (verbose) Rprintf("\t\tcurrent key indicates lower distance than current\n");
      vVec[current_key.id].dist = current_key.key;
      current_node = vVec[current_key.id];
      if (current_node.mark == false) {
        if (verbose) Rprintf("\t\tnode %d not in ball yet, adding it to ball\n",current_key.id);
        ballVec.push_back(current_key.id);
        vVec[current_key.id].mark = true;
      }
      node_el_size = current_node.out_edges.size();
      key_val = current_key.key;
      if (verbose) Rprintf("\t\tEnqueing node %d's neighbors: ",current_key.id);
      for(int ii=0;ii<node_el_size;ii++) {
        eId = current_node.out_edges[ii];
        if (ew[eId] >= 0.0) {
          temp_dist = key_val + (ew[eId]);
          if (temp_dist <= r){
            next_node = vVec[elTargets[eId]-1];
            next_key.id = next_node.id;
            next_key.key = temp_dist;
            if (verbose) {
              Rprintf("\n\t\t\tdistance of %d set to %f",next_key.id,next_key.key);
            }
            if (parents) {
              if (temp_dist < vVec[next_key.id].dist) {
                vVec[next_key.id].parent = current_key.id;
                if (verbose) {
                  Rprintf("\nparent of %d set to %d",next_key.id,current_key.id);
                }
              }
            }
            vQueue.push(next_key);
            //if (!parents){ if (!distances){ if (verbose) Rprintf("%d ",next_key.id);}}
          }
        }
      }
      if (verbose) Rprintf("\n");
      //finished=(vQueue.empty() || current_key.key > r);
    }
    R_FlushConsole();
  }
  //We used an stl vector to store ids. It needs to be
  //converted to a NumericVector to match the return type.
  b_size = ballVec.size();
  ri=0;
  if (parents) outrows++;
  if (distances) outrows++;
  outmat = NumericMatrix(outrows,b_size);
  outmat(ri,_) = Rcpp::as<NumericVector>(wrap(ballVec));
  if (distances) {
    ri++;
    for(int ii=0;ii<b_size;ii++){
      dVec.push_back(vVec[ballVec[ii]].dist);
    }
    outmat(ri,_) = Rcpp::as<NumericVector>(wrap(dVec));
  }
  if (parents) {
    ri++;
    for(int ii=0;ii<b_size;ii++){
      pVec.push_back(vVec[ballVec[ii]].parent);
    }
    outmat(ri,_) = Rcpp::as<NumericVector>(wrap(pVec));
  }

  return outmat;
}

//expandBall -  returns matrix containing the vertices of the ball
//              by default, distances and parents of each vertex are in 2nd and 3rd
//              rows. These must be fed in to any subsequent calls as ballData.
//Unlike growBall, expandBall which assumes no starting ball
//expandBall assumes that a previous ball is being fed in and tries to continue it.
//[[Rcpp::export]]
NumericMatrix expandBall(NumericVector vl,NumericMatrix eData, int rv, double r, 
  NumericMatrix ballData, int pv=-1,
  bool verbose=false, bool distances=true,bool parents=true) {
  //vl - list of vertices to grow ball over, rv - root of ball, r - radius of ball
  //eData - 3 x |E| matrix
  //    Row 1 - edge sources, Row 2 - edge targets, Row 3 - edge weights
  //ballData 3 x |V(ball)| matrix
  //    Row 1 - vertex ids, Row 3 - vertex distances, Row 4 - vertex parents
  //pv -  special node used during petal growing methods. Represents ending node of
  //      the path. Needed because path distance may change but the changes may not
  //      would not be updated otherwise as it often is on the interior of the petal.
  NumericVector vBall = ballData(0,_);
  NumericVector dBall = ballData(1,_);
  NumericVector pBall = ballData(2,_);
  NumericMatrix outmat;
  std::vector<node> vVec(vl.length());
  std::vector<int> restart_list, ballVec, dVec,pVec;
  std::priority_queue< real_key_int, std::vector<real_key_int>, 
                      std::greater<real_key_int> > vQueue;
  real_key_int current_key, next_key, last_key;
  node current_node, next_node, temp_node;
  int vl_size, el_size, node_el_size, eId,b_size, outrows=1, ri;
  double key_val, temp_dist,inf_dist = 99999.0;
  
  vl_size = vl.size();
  b_size = ballData.ncol();
  el_size = eData.nrow();
  if (verbose) {Rprintf("initializing vertex node list\n");
  R_FlushConsole();}
  for (int ii=0;ii<vl_size;ii++) {
    vVec[ii].dist = inf_dist;
    vVec[ii].id = ii;
    vVec[ii].mark = false;
  }
  if (verbose) {Rprintf("loading previous ball Data\n\tball Vertex count = %d\n",b_size);
  R_FlushConsole();}
  for (int ii=0;ii<b_size;ii++) {
    if (verbose) Rprintf("\tloading vertex: id %d, dist %f, parent %d \n",
        int(ballData(0,ii)),ballData(1,ii),int(ballData(2,ii)));
    vVec[ballData(0,ii)].dist = ballData(1,ii);
    vVec[ballData(0,ii)].parent = ballData(2,ii);
    vVec[ballData(0,ii)].mark = true;
    ballVec.push_back(ballData(0,ii));
  }
  R_FlushConsole();
  
  if (verbose) {Rprintf("loading edges into vertex list\n");
    Rprintf("\tthere are a total of %d edges listed\n",el_size);
    R_FlushConsole();}
  for (int ii=0;ii<el_size;ii++) {
    //Handle vertex weights here... not yet implemented
    vVec[eData(ii,0)-1].out_edges.push_back(ii);
  }
/*  if (verbose) {
    for (int ii=0;ii<vl_size;ii++){
      node_el_size = vVec[ii].out_edges.size();
      Rprintf("\nlisting a total of %d edges for vertex %d: \n\t",node_el_size,ii);
      for (int jj=0;jj<node_el_size;jj++) {
        Rprintf("%1.0f--%f-->%1.0f ",eData(vVec[ii].out_edges[jj],0),eData(vVec[ii].out_edges[jj],2),
        eData(vVec[ii].out_edges[jj],1));
      }
    }
  } // This produces a LOT of printing!*/
  
  if (verbose) {Rprintf("\nloading ball vertices into queue\n");
  R_FlushConsole();
  Rprintf("vQueue: ");}
  restart_list = upper_ids(ballData(1,_));
  for (int ii=0;ii<restart_list.size();ii++) {
    temp_node = vVec[ballData(0,restart_list[ii])];
    node_el_size = temp_node.out_edges.size();
    for (int jj=0;jj<node_el_size;jj++) {
      current_key.id = eData(temp_node.out_edges[jj],1)-1;
      current_key.key = temp_node.dist + eData(temp_node.out_edges[jj],2);
      vQueue.push(current_key);
      if (verbose) Rprintf("%d ",current_key.id);
    }
  }
  //Load special vertex into list if needed. Only used for petal growing.
  if (pv >= 0) {
  current_key.id = pv-1;
  current_key.key = vVec[pv-1].dist;
  vVec[pv-1].dist = inf_dist;
  vQueue.push(current_key);
  temp_node = vVec[pv-1];
  node_el_size = temp_node.out_edges.size();
  for (int jj=0;jj<node_el_size;jj++) {
      current_key.id = eData(temp_node.out_edges[jj],1)-1;
      current_key.key = temp_node.dist + eData(temp_node.out_edges[jj],2);
      vQueue.push(current_key);
      if (verbose) Rprintf("%d ",current_key.id);
    }
  /*if (verbose) Rprintf("%d ",current_key.id);Rprintf("\n");*/ }
  
  while (!vQueue.empty()) {
    last_key = current_key;
    current_key = vQueue.top();
    if (verbose) Rprintf("popped key for vertex %d from stack\n",current_key.id);
    if (verbose) Rprintf("\t key distance was %4.2f, vertex at key id had distance %f\n",
                          current_key.key,vVec[current_key.id].dist);
    vQueue.pop();
    if (current_key.key < vVec[current_key.id].dist && \
        current_key.key < inf_dist) {
      if (verbose) Rprintf("\t\tcurrent key indicates lower distance than current\n");
      vVec[current_key.id].dist = current_key.key;
      current_node = vVec[current_key.id];
      if (current_node.mark == false) {
        if (verbose) Rprintf("\t\tnode %d not in ball yet, adding it to ball\n",current_key.id);
        ballVec.push_back(current_key.id);
        vVec[current_key.id].mark = true;
      }
      node_el_size = current_node.out_edges.size();
      key_val = current_key.key;
      if (verbose) Rprintf("\t\tnode %d has %d neighbors\n",current_key.id,node_el_size);
      if (verbose) Rprintf("\t\tEnqueing node %d's neighbors: ",current_key.id);
      for(int ii=0;ii<node_el_size;ii++) {
        eId = current_node.out_edges[ii];
        if (verbose) Rprintf("\n\t\t\tchecking vertex %d",int(eData(eId,1)-1));
        if (eData(eId,2) >= 0.0) {
          temp_dist = key_val + (eData(eId,2));
          if (temp_dist <= r){
            next_node = vVec[eData(eId,1)-1];
            next_key.id = next_node.id;
            next_key.key = temp_dist;
            if (verbose) {
              Rprintf("\n\t\t\tdistance of key %d set to %f",next_key.id,next_key.key);
            }
            if (parents) {
              if (temp_dist < vVec[next_key.id].dist) {
                vVec[next_key.id].parent = current_key.id;
                if (verbose) {
                  Rprintf("\n\t\t\tparent of %d set to %d",next_key.id,current_key.id);
                }
              }
            }
            vQueue.push(next_key);
            //if (!parents){ if (!distances){ if (verbose) Rprintf("%d ",next_key.id);}}
          }
        }
      }
      if (verbose) Rprintf("\n");
    }
    R_FlushConsole();
  }
  if (verbose) Rprintf("constructing output matrix\n");
  R_FlushConsole();
  //We used an stl vector to store ids. It needs to be
  //converted to a NumericVector to match the return type.
  b_size = ballVec.size();
  ri=0;
  if (parents) outrows++;
  if (distances) outrows++;
  outmat = NumericMatrix(outrows,b_size);
  outmat(ri,_) = Rcpp::as<NumericVector>(wrap(ballVec));
  if (distances) {
    ri++;
    for(int ii=0;ii<b_size;ii++){
      dVec.push_back(vVec[ballVec[ii]].dist);
    }
    outmat(ri,_) = Rcpp::as<NumericVector>(wrap(dVec));
  }
  if (parents) {
    ri++;
    for(int ii=0;ii<b_size;ii++){
      pVec.push_back(vVec[ballVec[ii]].parent);
    }
    outmat(ri,_) = Rcpp::as<NumericVector>(wrap(pVec));
  }
  if (verbose) Rprintf("returning output:\n");
  R_FlushConsole();
  return outmat;

}
//[[Rcpp::export]]
void numMatTester(NumericMatrix mat) {
  Rprintf("matrix has %d rows\n",mat.nrow());
  Rprintf("matrix has %d colums",mat.ncol());
  for (int ii=0;ii<mat.nrow();ii++){
   Rprintf("\n\trow %d: ",ii);
   for(int jj=0;jj<mat.ncol();jj++){
     Rprintf("%f ",mat(ii,jj));
   }
  }
}