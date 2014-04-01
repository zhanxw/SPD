#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <deque>
#include "IO.h"
using namespace std;

#undef DEBUG
//#define DEBUG

#define PLUS_ZERO 1e-5

struct Node {
  Node() {
    degree = 0;
    neib[0] = neib[1] = neib[2] = 0;
    dist[0] = dist[1] = dist[2] = 0;
    outerNode = false;
  };
  int addNeighbor(const int to, const double d) {
    // check pre-existing
    for (int i = 0; i < degree; i++) {
      if (neib[i] == to) {
        fprintf(stdout, "Duplicate connection: to %d\n", to);
      }
    }
    // add connection
    neib[degree] = to;
    dist[degree] = d;
    degree ++;

    // maximum degree is 3
    assert(degree <= 3);

    if (degree == 1) {
      outerNode = true;
    } else {
      outerNode = false;
    }
    return 0;
  };
  int degree;
  int neib[3];
  double dist[3];
  bool outerNode;
};

struct NodeSet{
  vector<Node> nodes;
  vector<string> labels;
  map<string, int> label2idx;

  // @return idx
  int addNode(const string& nodeName) {
    if (label2idx.find(nodeName) == label2idx.end()) {
      int idx = label2idx.size();
      Node n;
      nodes.push_back(n);
      labels.push_back(nodeName);
      label2idx[nodeName] = idx;
      return idx;
    }
    return label2idx[nodeName];
  };

  void addDistance(const int from, const int to, const double dist) {
    nodes[from].addNeighbor(to, dist);
    nodes[to].addNeighbor(from, dist);
  };
  int findIndex(const string& nodeName) {
    if (label2idx.find(nodeName) == label2idx.end()) {
      return -1;
    }
    return label2idx[nodeName];
  }
  /**
   * Find the minimum non-zero distance between two nodes, and count total pairs
   * e.g. suppose distances are 0, 0, 1, 1, 1, 3
   * the minimum distance is 1 - 0 = 1
   * the count is 2 + 3 = 5
   */
  double findMinimumDistanceDifference(int* tie) {
    int& s = *tie;
    s = 0;
    double minDist = 1e20;
    
    // get all non-zero distances
    map<double, int> dist;
    for (size_t i = 0; i < nodes.size(); ++i) {
      for (int j = 0; j < nodes[i].degree; ++j) {
        if (nodes[i].dist[j] == 0.0) continue;
        
        dist [ nodes[i].dist[j] ] ++;
      }
    }
    map<double, int>::const_iterator it = dist.begin();
    double last = it->first;
    double lastCount = it->second;
    double current;
    double diff;
    ++it;
    while (it != dist.end()) {
      // printf("current = %lf, count = %d\n", it->first, it->second);
      current = it->first;
      diff = current - last;
      if (minDist > diff) {
        minDist = diff;
        s = lastCount + it->second ;
      }
      last = it->first;
      lastCount = it->second;
      
      ++it;      
    }
    // printf("minDist = %lf, count = %d\n", minDist, s);
    return minDist;
  }
  void randomizeDistance() {
    int s = 0;
    double d = 1e6;
    d = findMinimumDistanceDifference(&s);
    d /= s;
    assert(d > 0);
    fprintf(stdout, "Randomize unit is %lf\n", d);

    for (size_t i = 0; i < nodes.size(); ++i) {
      for (int j = 0; j < nodes[i].degree; ++j) {
        double runif = 1.0 * rand() /RAND_MAX;
        nodes[i].dist[j] += d * runif;
      }
    }
  };
}; // end struct NodeSet

template <class T>
int whichMax(T* vec, int len) {
  int idx = -1;
  if (len <= 0) return idx;

  T max = vec[idx];

  for (int i = 1; i < len; i++) {
    if ( vec[i] > max) {
      idx = i;
      max = vec[idx];
    }
  };

  return idx;
};

double findMaxApartedNodes(const NodeSet& nodeSet, int* from, int* to) {
  const vector<Node>& nodes = nodeSet.nodes;
  assert(from && to);
  int l = nodes.size();
  if (l <= 0) {
    return -1.0;
  }

  int numNode = nodes.size();
  double* dist = new double[ numNode ];
  for (int i = 0; i < numNode; i++ ) dist[i] = -1.0;

  // find an outer node
  int outerNode = -1;
  for (int i = 0; i < numNode; i++) {
    if (nodes[i].outerNode) {
      outerNode = i;
      break;
    }
  }
#ifdef DEBUG
  fprintf(stdout, "Found outer node %d.\n", outerNode);
#endif

  deque<int> toExplore;
  toExplore.push_back(outerNode);
  dist[outerNode] = 0.0;
  while (toExplore.size() > 0) {
    int n = toExplore.front();
    double d = dist[n];
    for (int i = 0; i < nodes[n].degree; i ++ ) {
      int newNode = nodes[n].neib[i];
      if (dist[newNode] < -0.1) { // unvisit, orelse dist[] should be -1.
        dist[newNode] = d + nodes[n].dist[i];
        toExplore.push_back(newNode);
      }
    }
    toExplore.pop_front();
  };
  int maxFrom = whichMax(dist, numNode);
  while (!nodes[maxFrom].outerNode){
#ifdef DEBUG
    fprintf(stdout, "skip inner node %d [%s]\n", maxFrom, nodeSet.labels[maxFrom].c_str());
#endif
    dist[maxFrom] -= 1000.0;
    maxFrom = whichMax(dist, numNode);
  }

#ifdef DEBUG
  fprintf(stdout, "max dist is: %d -> %d [dist = %.2f] \n", outerNode, maxFrom, dist[maxFrom]);
#endif

  for (int i = 0; i < numNode; i++ ) dist[i] = -1.0;
  toExplore.push_back(maxFrom);
  dist[maxFrom] = 0.0;
  while (toExplore.size() > 0) {
    int n = toExplore.front();
    double d = dist[n];
    for (int i = 0; i < nodes[n].degree; i ++ ) {
      int newNode = nodes[n].neib[i];
      if (dist[newNode] < -0.1) { // unvisit, orelse dist[] should be -1.
        dist[newNode] = d + nodes[n].dist[i];
        toExplore.push_back(newNode);
      }
    }
    toExplore.pop_front();
  };
#ifdef DEBUG
  for (int i = 0; i < numNode; i ++ ) {
    fprintf(stdout, "dist[%d] %s = %.2f\n", i, nodeSet.labels[i].c_str(), dist[i]);
  }
#endif
  int maxTo = whichMax(dist, numNode);
  while (!nodes[maxTo].outerNode){
    dist[maxTo] -= 100.0;
    maxTo = whichMax(dist, numNode);
  }
#ifdef DEBUG
  fprintf(stdout, "max dist is: %d -> %d [dist = %.2f] \n", maxFrom, maxTo, dist[maxTo]);
#endif

  *from = maxFrom;
  *to = maxTo;
  double ret = dist[maxTo];
  delete [] dist;
  return ret;
};


// find distance from @param node to @param toReach, but not visa nodes in @param visited
double findDistToPath(const vector<Node>& nodes, const int node, const set<int>& dest, set<int>* visited, double* dist) {
  if (dest.find(node) != dest.end() ) {
    return *dist;
  }

  int degree = nodes[node].degree;
  const Node& n = nodes[node];
  for (int i = 0; i < degree; i ++) {
    int j = n.neib[i];
    if (visited->find(j) != visited->end()) {
      continue;
    }
    visited->insert(j);
    (*dist) += n.dist[i];

    if (dest.find(j) != dest.end()) { // dest is reached
      return *dist;
    }
    if (findDistToPath(nodes, j, dest, visited, dist) < 0) { // cannot reach dest via j
      visited->erase(j);
      (*dist) -= n.dist[i];
      continue;
    } else {
      return *dist;
    }
  }
  return -999;
};

void unionPath(set<int>* dest, const set<int>& src) {
  set<int>::const_iterator iter;
  for (iter = src.begin(); iter != src.end(); iter++){
    dest->insert( *iter);
  }
}


void usage(int argc, char** argv) {
  fprintf(stdout, "%s dist_file tree_file result_file: select most diverged populations from given tree file.\n", argv[0]);
};

void removeEmpty(vector<string>& vector) {
  int idx = 0;
  for (int i = 0; i < (int)vector.size(); i++){
    if ( vector[i].size() > 0) {
      vector[idx++] = vector[i];
    }
  }
  vector.resize(idx);
};


void findMaxDistNode(const char* fn, string* n1, string* n2, double* d) {
  LineReader lr(fn);
  vector< string> fd;
  int num = -1;
  vector<string> label;
  double maxDist = -1;
  int maxIndex = -1;
  while (lr.readLineBySep(&fd, "\t ")) {
    removeEmpty(fd);
    if (num < 0) {
      num = atoi(fd[0].c_str());
      fprintf(stderr, "There are %d nodes.\n", num);
      continue;
    }
    if ( (int)fd.size() != num + 1) {
      fprintf(stderr, "Cannot match lines, fd = %zu\n", fd.size());
      exit(1);
    }
    label.push_back(fd[0]);
    const int n = fd.size();
    for (int i = 1; i < n; ++i) {
      double d = atof(fd[i].c_str());
      if (d > maxDist) {
        (*n1) = fd[0];
        maxIndex = i - 1;
        maxDist = d;
      }
    }
  }
  (*n2) = label[maxIndex];
  (*d) = maxDist;
  fprintf(stdout, "Max dist = %lf between node %s and node %s\n", maxDist, n1->c_str(), n2->c_str());
  assert( maxDist > 0);
};


int main(int argc, char** argv) {
  if (argc != 4) {
    usage(argc, argv);
    exit(1);
  }
  string maxNode1;
  string maxNode2;
  double dist; 
  findMaxDistNode(argv[1], &maxNode1, &maxNode2, &dist);


  FileWriter fw(argv[3]);

  NodeSet nodeSet;

  LineReader lr(argv[2]);
  std::vector< std::string> fd;
  while (lr.readLineBySep(&fd, "\t ")) {
    removeEmpty(fd); // only keep non-empty columns
    if (fd.size() == 0) continue;
    assert(fd.size() == 3);

    double dist = atof(fd[2].c_str());
    int from = nodeSet.addNode(fd[0]);
    int to = nodeSet.addNode(fd[1]);

    nodeSet.addDistance(from, to, dist);
  }
  int numOuterNode = 0;
  int numTotalNode = nodeSet.nodes.size();
  for (int i = 0; i < numTotalNode ; i++ ){
    if (nodeSet.nodes[i].outerNode)
      numOuterNode ++;
  };

  fprintf(stdout, "Total %d nodes loaded, %d are outer nodes / species.\n", numTotalNode, numOuterNode);

  nodeSet.randomizeDistance();
  
#ifndef DEBUG
  int tempNode[] = {1, 17, 27};
  fprintf(stdout, "%d <- %s, %d <- %s, %d <- %s. \n",
          tempNode[0],
          nodeSet.labels[tempNode[0]].c_str(),
          tempNode[1],
          nodeSet.labels[tempNode[1]].c_str(),
          tempNode[2],
          nodeSet.labels[tempNode[2]].c_str());
#endif

  // find maximum aparted nodes from the tree we readed
  int maxFrom, maxTo;
  // double dist;
  // dist = findMaxApartedNodes(nodeSet, &maxFrom, &maxTo);
  // assert(dist > 0.0);
  maxFrom = nodeSet.findIndex(maxNode1);
  maxTo = nodeSet.findIndex(maxNode2);
  assert(maxFrom > 0 && maxTo > 0);
  
  fprintf(stdout, "Furthest aparted nodes: %s and %s with distance %.2f.\n", nodeSet.labels[maxFrom].c_str(), nodeSet.labels[maxTo].c_str(), dist);
  fw.writeLine(nodeSet.labels[maxFrom].c_str());
  fw.writeLine(nodeSet.labels[maxTo].c_str());

  // find path
  set<int> path;
  set<int> visited;
  path.insert(maxFrom);
  visited.insert(maxTo);
  assert(findDistToPath(nodeSet.nodes, maxFrom, visited, &path, &dist) > 0);
#ifdef DEBUG
  fprintf(stdout, "Max distance is %.2f\n", dist);
#endif
  unionPath(&visited, path);

  // gradually add furthest aparted nodes in.
  while ( (int)visited.size() < numTotalNode) {
    double maxDist = -1.0;
    int maxIdx = -1;
    set<int> maxPath;
    for (int i = 0 ; i < numTotalNode; i ++ ) {
      if (!nodeSet.nodes[i].outerNode) continue;
      if (visited.find(i) != visited.end()) { // previous visited
        continue;
      }
#ifdef DEBUG
      // fprintf(stdout, "Check node %s \n", nodeSet.labels[i].c_str());
#endif
      dist = 0.0;
      path.clear();
      path.insert(i);
      assert (findDistToPath(nodeSet.nodes, i, visited, &path, &dist) >= 0.0 );
      if (dist > maxDist) {
        maxIdx = i;
        maxDist = dist;
        maxPath = path;
      }
    };
    assert (maxIdx > 0);
    unionPath(&visited, maxPath);

    fprintf(stdout, "Add new node %s with distance %.2f.\n", nodeSet.labels[maxIdx].c_str(), maxDist);
    fw.writeLine(nodeSet.labels[maxIdx].c_str());
  }
  fw.close();
};
