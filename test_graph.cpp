// g++ -Wall -std=c++11 test_graph.cpp -o test_graph

#include <iostream>
#include <vector>
#include <map>
#include <cstdlib>

#include <queue>
#include <algorithm> // for string reverse

using namespace std;

//#define DEBUG
#define EXAMPLE_W2


inline double prob() {
  return (static_cast<double>(rand())/RAND_MAX);
}


class Graph {
friend ostream& operator<<(ostream& out, const Graph& graph);
public:
  Graph(int nodes = 50) : nodes(nodes) {
    reset();
  };

  void reset() {
    mat.resize(nodes, vector<double>(nodes, -1)); // matrice of size nodes * nodes initialized with 0 values
    closed_set = false; // whether closed set is known or not
    closed_set_src = -1; // src from which all shortest paths (stored in closed_set) are computed
    closed_set_dist.resize(nodes, -1);
    closed_set_prev.resize(nodes, -1);

    for (int i = 0; i < nodes; i++) {
      for (int j = 0; j < nodes; j++) {
        mat[i][j] = -1;
      }
    }
    for (int i = 0; i < nodes; i++) {
      add_edge(i, i, 0);
      closed_set_dist[i] = -1;
      closed_set_prev[i] = -1;
    }
  }

  double monte_carlo_simu(double density, double range) {
    //cout << "DBG:" << density << ", " << range << endl;
    reset();
    for (int i = 0; i < nodes; i++) {
      add_edge(i, i, 0);
      for (int j = i + 1; j < nodes; j++) { // from j = i + 1: undirected graph
        if (prob() < density) {
          add_edge(i, j, 1 + prob() * (range-1));
        }
      }
    }
    shortest_paths(0);
    //print_shortest_paths();
    return average_shortest_paths();
  }

  void add_edge(int x, int y, double cost) {
    if (x < nodes && y < nodes)
      mat[x][y] = mat[y][x] = cost;
  };

  string path(int src, int dst) {
    string path = to_string(dst);
    path += " >- "; // will become -> after reverse

    if (!closed_set)
      shortest_paths(src);

    int node = dst;
    while (closed_set_prev[node] >= 0 && closed_set_prev[node] != src) {
      node = closed_set_prev[node];
      path += to_string(node);
      path += " >- ";
    }
    path += to_string(src);
    reverse(path.begin(), path.end());
    return path;
  }

  double path_size(int src, int dst) {
    if (!closed_set)
      shortest_paths(src);

    return closed_set_dist[dst];
  }

  void print_shortest_paths() {
    if (closed_set) {
      cout << "Shortest paths from " << closed_set_src << endl;
      for (int i = 0; i < nodes; i++) {
        cout << "\tnode " << i << " => " << "dist=" << closed_set_dist[i] << " " << "prev=" << closed_set_prev[i] << endl;
      }
    }
  }

  double average_shortest_paths() {
    double sum = 0.0;
    int num = 0;

    if (closed_set) {
      for (int i = 0; i < nodes; i++) {
        if (closed_set_dist[i] > 0) {
          sum += closed_set_dist[i];
          num += 1;
        }
      }
      if (num > 0)
        return (sum/num);
    }
    return -1;
  }

private:
  void shortest_paths(int src) { // cf wikipedia Dijkstra: UniformCostSearch
    closed_set = true;
    closed_set_src = src;

    vector<double> open_set_dist(nodes, -1); // faster for 'belong to open_set' test
    typedef tuple<double, int, int> dii; // dist from src, previous node, node
    priority_queue<dii, vector<dii>, greater<dii>> open_set_pq; // faster for 'sort open_set'

    open_set_pq.push(make_tuple(0.0, src, src));

#ifdef DEBUG
    int step = 1;
#endif

    while (!open_set_pq.empty()) {
      dii node = open_set_pq.top();
      double node_dist = get<0>(node);
      int node_prev = get<1>(node);
      int node_id = get<2>(node);
      open_set_pq.pop();

      // if not already in closed set
      if (closed_set_dist[node_id] < 0) {
        closed_set_dist[node_id] = node_dist;
        closed_set_prev[node_id] = node_prev;

#ifdef DEBUG
        cout << "STEP " << step++ << endl;
        for (int i = 0; i < nodes; i++) {
          if (closed_set_dist[i] >= 0) {
            cout << "NODE " << i << ": " << "dist=" << closed_set_dist[i] << " prev=" << closed_set_prev[i] << endl;
          }
        }
#endif

        for (int i = 0; i < nodes; i++) {
          // if neighbor && not in closed_set
          if (mat[node_id][i] > 0 && closed_set_dist[i] < 0) {
            double new_dist = node_dist + mat[node_id][i];
            // if not in open_set or shorter path in open_set found
            if (open_set_dist[i] < 0 || new_dist < open_set_dist[i]) {
              open_set_dist[i] = new_dist;
              open_set_pq.push(make_tuple(new_dist, node_id, i));
            }
          }
        } // end for
      } // end if
    } // end while

  };

  int nodes; // number of nodes
  vector< vector<double> > mat; // connectivity matrix

  bool closed_set; // whether closed set is known or not
  int closed_set_src; // src from which all shortest paths (stored in closed_set) are computed
  vector<double> closed_set_dist;
  vector<int> closed_set_prev;
};

ostream& operator<<(ostream& out, const Graph& graph) {
  out << "Graph adjacency\n";
  for (int i = 0; i < graph.nodes; i++) {
    for (int j = i + 1; j < graph.nodes; j++) { // j = i + 1: undirected graph
      if (graph.mat[i][j] > 0)
        out << "\tbetween " << i << " and " << j << ": " << graph.mat[i][j] << endl;
    }
  }
  return out;
}

int main() {
  srand(time(0));

#ifdef EXAMPLE_W2
  // unitary test
  Graph g(9);
  g.add_edge(0, 1, 4.0);
  g.add_edge(0, 2, 3.0);
  g.add_edge(0, 4, 7.0);

  g.add_edge(1, 3, 1.0);

  g.add_edge(2, 4, 4.0);

  g.add_edge(3, 4, 3.0);
  g.add_edge(3, 5, 1.0);

  g.add_edge(4, 5, 1.0);
  g.add_edge(4, 6, 5.0);
  g.add_edge(4, 8, 3.0);

  g.add_edge(5, 7, 2.0);
  g.add_edge(5, 8, 4.0);

  g.add_edge(7, 8, 3.0);

  cout << g;
  cout << "shortest path from 0 to 8: " << g.path(0, 8) << endl;
  cout << "shortest path from 0 to 8: " << g.path_size(0, 8) << endl;
  g.print_shortest_paths();
  cout << "average_shortest_paths: " << g.average_shortest_paths() << endl << endl << endl;
#endif

  Graph g50(50);
  int nsimu = 1000;
  double path20 = 0;
  double path40 = 0;
  for (int i = 0; i < nsimu; i++) {
    path20 += g50.monte_carlo_simu(0.2, 10);
    path40 += g50.monte_carlo_simu(0.4, 10);
  }

  cout << "Average shortest path for a Graph with 50 nodes, density 20%, range 10: " << path20/nsimu << endl;
  cout << "Average shortest path for a Graph with 50 nodes, density 40%, range 10: " << path40/nsimu << endl;

#ifdef DEBUG
  Graph g(10);
  g.add_edge(0, 3, 2.81);
  g.add_edge(0, 2, 5.98);
  g.add_edge(2, 3, 9.39);
  g.add_edge(2, 7, 4.86);
  g.add_edge(7, 5, 7.12);
  g.add_edge(5, 1, 3.28);
  g.add_edge(3, 8, 4.94);
  g.add_edge(8, 1, 5.61);
  g.add_edge(1, 9, 3.64);
  g.add_edge(9, 4, 4.84);

  cout << g;
  cout << "shortest path from 0 to 8: " << g.path(0, 5) << endl;
  cout << "shortest path from 0 to 8: " << g.path_size(0, 5) << endl;
  g.print_shortest_paths();
  cout << "average_shortest_paths: " << g.average_shortest_paths() << endl;
#endif


  return 1;
}
