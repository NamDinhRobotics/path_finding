// g++ -g -Wall -std=c++11 test_graph.cpp -o test_graph

#include <iostream> // for cout
#include <vector> // for vector container (main container here)
#include <queue> // for priority_queue
#include <algorithm> // for string reverse operation

#include <fstream>
#include <iterator>
#include <cfloat> // for DBL_MAX

using namespace std;

// Debugging and unitary tests flags
#define TEST_MST


// computes a probability value between [0;1]
inline double prob() 
{
  return (static_cast<double>(rand())/RAND_MAX);
}

// Graph implementation: with methods for
//  constructing a graph
//  computing shortest paths with dijkstra algorithm
//  Monte Carlo simulation to compute average shortest path
//  During Monte Carlo simulation a Graph is characterized by:
//    its density: percentage of edges per node 
//    its range: path distance between [1; range] for connected nodes

// In terms of data storage and representation:
//    nodes are integer from 0 to N-1 (N being the number of nodes)
//    distances are double
//    a connectivity matrix is used to represent edges (alternatively an adjacency list could have been used)
//    a closed_set for shortest paths information storage is used
//    vector is the main container used: for direct O(1) access time
//    a priority_queue is used during Djikstra: for fast sort/dequeing operations 

// 2 small examples are provided for debugging / checking / unitary tests
// (cf DEBUG and EXAMPLE_W2)

// The Monte Carlo simulation results are for 1000 simulations average of average shortest path:
// Average shortest path for a Graph with 50 nodes, density 20%, range 10: 6.97412
// Average shortest path for a Graph with 50 nodes, density 40%, range 10: 4.70775


class Graph 
{
  friend ostream& operator<<(ostream& out, const Graph& graph);

  public:
    // constructor with 50 nodes by default
    Graph(int nodes = 50) : nodes(nodes) 
    {
      reset();
    };

    Graph(string fname)
    {
      ifstream data_file(fname);
      istream_iterator<int> start(data_file), end;
      vector<int> data(start, end);

      //reset();

      cout << data[0] << endl;
      //auto it = start;
      cout << *start << endl;
      cout << data.size() << endl;
    };
  
    // reset/init all Graph data execot number of nodes
    // typically before starting a new Monte Carlo simulation
    void reset() 
    {
      mat.resize(nodes, vector<double>(nodes)); // matrice of size nodes * nodes initialized with 0 values
      closed_set = false; // whether closed set is known or not
      closed_set_src = -1; // src from which all shortest paths (stored in closed_set) are computed
      closed_set_dist.resize(nodes);
      closed_set_prev.resize(nodes);
  
      for (int i = 0; i < nodes; i++) 
      {
        for (int j = 0; j < nodes; j++) 
        {
          mat[i][j] = -1;
        }
      }

      for (int i = 0; i < nodes; i++) 
      {
        add_edge(i, i, 0);
        closed_set_dist[i] = -1;
        closed_set_prev[i] = -1;
      }
    }

  
    // Monte Carlo simulation that computes average shortest path
    double monte_carlo_simu(double density, double range, int nsimu) 
    {
      double path = 0;

      for (int i = 0; i < nsimu; i++) 
      {
        reset();
        for (int i = 0; i < nodes; i++) 
        {
          add_edge(i, i, 0);
          for (int j = i + 1; j < nodes; j++) // from j = i + 1: undirected graph
          {
            if (prob() < density) 
              add_edge(i, j, 1 + prob() * (range-1));
          }
        }
        shortest_paths(0);
        //print_shortest_paths();
        path += average_shortest_path();
      }
      return path/nsimu;
    }
  
    // add an edge into an existing connectivity matrix
    void add_edge(int x, int y, double cost) 
    {
      if (x < nodes && y < nodes)
        mat[x][y] = mat[y][x] = cost;
      else
        cout << "WARNING: edge " << x << "," << y << " with non existing vertex" << endl;
    };
  
    // return the path from src to dst as a string
    string path(int src, int dst) 
    {
      string path = to_string(dst);
      path += " >- "; // will become -> after reverse
  
      if (!closed_set || closed_set_src != src)
        shortest_paths(src);
  
      int node = dst;
      while (closed_set_prev[node] >= 0 && closed_set_prev[node] != src) 
      {
        node = closed_set_prev[node];
        path += to_string(node);
        path += " >- ";
      }
      path += to_string(src);
      reverse(path.begin(), path.end());
      return path;
    }
  
    // return the shortest possible distance from src to dst
    double path_size(int src, int dst) 
    {
      if (!closed_set || closed_set_src != src)
        shortest_paths(src);
  
      return closed_set_dist[dst];
    }
  
    // print all shortest paths (mainly for debugging)
    void print_shortest_paths() 
    {
      if (closed_set) 
      {
        cout << "Shortest paths from " << closed_set_src << endl;
        for (int i = 0; i < nodes; i++) 
          cout << "\tnode " << i << " => " << "dist=" << closed_set_dist[i] << " " << "prev=" << closed_set_prev[i] << endl;
      }
    }
  
    // compute average shortest path (mainly for Monte Carlo simulation)
    double average_shortest_path() 
    {
      double sum = 0.0;
      int num = 0;
  
      if (closed_set) 
      {
        for (int i = 0; i < nodes; i++) 
        {
          if (closed_set_dist[i] > 0) 
          {
            sum += closed_set_dist[i];
            num += 1;
          }
        }
        if (num > 0)
          return (sum/num);
      }
      return -1;
    }

    // Compute Minimum Spanning Tree via Prim algorithm (cf https://en.wikipedia.org/wiki/Prim%27s_algorithm) 
    // using a priority queue for faster execution
    // The logic and code is very similar to Dijkstra implementation
    bool mst(void)
    {
      int mst_len = 0;
      mst_set = false;
      mst_cost = 0;
      mst_set_dist.resize(nodes, DBL_MAX);
      mst_set_prev.resize(nodes, -1);

      typedef tuple<double, int, int> dii; // dist from mst (must be first for sorting), edge from mst, node id
      priority_queue<dii, vector<dii>, greater<dii>> open_set_pq; // faster for 'sort/dequeue' on open_set

      open_set_pq.push(make_tuple(0.0, 0, 0));

      while (!open_set_pq.empty())
      {
        dii node = open_set_pq.top();
        double node_dist = get<0>(node);
        int node_prev = get<1>(node);
        int node_id = get<2>(node);
        open_set_pq.pop();

        // if not yet in MST
        if (mst_set_prev[node_id] < 0)
        {
          mst_set_dist[node_id] = node_dist;
          mst_set_prev[node_id] = node_prev;
          mst_len++;
          mst_cost += node_dist;

          for (int i = 0; i < nodes; i++)
          {
            // if connection to node i exists && node i not yet in MST => add node i to open_set
            if (mat[node_id][i] > 0 && mst_set_prev[i] < 0)
            {
              open_set_pq.push(make_tuple(mat[node_id][i], node_id, i));
            }
          }

          // for debug
          cout << "mst_len=" << mst_len << endl;
          print_mst();
        }
      }

      return mst_set = (mst_len == nodes);
    }

    // print Minimum Spanning Tree
    void print_mst() 
    {
      cout << "Minimum Spanning tree: " << mst_set << " mst_cost=" << mst_cost << endl;
      for (int i = 0; i < nodes; i++)
      {
        if (mst_set_prev[i] >= 0)
          cout << "\tnode " << i << " => " << "dist=" << mst_set_dist[i] << " " << "prev=" << mst_set_prev[i] << endl;
      }
    }
  
  private:
    // compute shortest paths as per wikipedia Dijkstra UniformCostSearch method
    void shortest_paths(int src) 
    {
      reset_closed_set();
      closed_set = true;
      closed_set_src = src;
  
      vector<double> open_set_dist(nodes, -1); // faster for 'belong to open_set' test
      typedef tuple<double, int, int> dii; // dist from src (must be first for sorting), previous node, node
      priority_queue<dii, vector<dii>, greater<dii>> open_set_pq; // faster for 'sort/dequeue' on open_set
  
      open_set_pq.push(make_tuple(0.0, src, src));
  
#ifdef DEBUG
      int step = 1;
#endif
  
      while (!open_set_pq.empty()) 
      {
        dii node = open_set_pq.top();
        double node_dist = get<0>(node);
        int node_prev = get<1>(node);
        int node_id = get<2>(node);
        open_set_pq.pop();
  
        // if not already in closed set
        if (closed_set_dist[node_id] < 0) 
        {
          closed_set_dist[node_id] = node_dist;
          closed_set_prev[node_id] = node_prev;
  
#ifdef DEBUG
          cout << "STEP " << step++ << endl;
          for (int i = 0; i < nodes; i++) 
          {
            if (closed_set_dist[i] >= 0) 
              cout << "NODE " << i << ": " << "dist=" << closed_set_dist[i] << " prev=" << closed_set_prev[i] << endl;
          }
#endif
  
          for (int i = 0; i < nodes; i++) 
          {
            // if neighbor && not in closed_set
            if (mat[node_id][i] > 0 && closed_set_dist[i] < 0) 
            {
              double new_dist = node_dist + mat[node_id][i];
              // if not in open_set or shorter path in open_set found
              if (open_set_dist[i] < 0 || new_dist < open_set_dist[i]) 
              {
                open_set_dist[i] = new_dist;
                open_set_pq.push(make_tuple(new_dist, node_id, i));
              }
            }
          } // end for
        } // end if
      } // end while
  
    };

    void reset_closed_set() 
    {
      closed_set = false; // whether closed set is known or not
      closed_set_src = -1; // src from which all shortest paths (stored in closed_set) are computed
      closed_set_dist.resize(nodes);
      closed_set_prev.resize(nodes);

      for (int i = 0; i < nodes; i++) 
      {
        closed_set_dist[i] = -1;
        closed_set_prev[i] = -1;
      }
    }
  
    int nodes; // number of nodes
    vector< vector<double> > mat; // connectivity matrix
  
    // closed_set stores information about all shortest paths
    bool closed_set; // whether closed set is known or not
    int closed_set_src; // src from which all shortest paths (stored in closed_set) are computed
    vector<double> closed_set_dist; // shortest distance from src to nodes
    vector<int> closed_set_prev; // previous node on the shortest path

    // Minnimum Spanning Tree
    bool mst_set; // whether MST exists or not
    double mst_cost;
    vector<double> mst_set_dist; // array of distances from the closest edge (cf below) to each vertex
    vector<int> mst_set_prev; // array indicating, for a given vertex, which vertex in the tree it is closest to
};

ostream& operator<<(ostream& out, const Graph& graph) 
{
  out << "Graph connectivity\n";
  for (int i = 0; i < graph.nodes; i++) 
  {
    for (int j = i + 1; j < graph.nodes; j++) // j = i + 1: undirected graph
    {
      if (graph.mat[i][j] > 0)
        out << "\tbetween " << i << " and " << j << ": " << graph.mat[i][j] << endl;
    }
  }
  return out;
}

// test1 from example W2
void test1_shortest_path()
{
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
  cout << "average_shortest_path: " << g.average_shortest_path() << endl << endl << endl;
}

// test2: more tricky
void test2_shortest_path()
{
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
  cout << "average_shortest_path: " << g.average_shortest_path() << endl;
}

void test1_mst()
{
  cout << "UNITARY TEST: " << __func__ << endl;
  Graph g(9);
  
  g.add_edge(0, 1, 4.0);
  g.add_edge(0, 7, 8.0);

  g.add_edge(1, 2, 8.0);
  g.add_edge(1, 7, 11.0);

  g.add_edge(2, 3, 7.0);
  g.add_edge(2, 5, 4.0);
  g.add_edge(2, 8, 2.0);

  g.add_edge(3, 4, 9.0);
  g.add_edge(3, 5, 14.0);

  g.add_edge(4, 5, 10.0);

  g.add_edge(5, 6, 2.0);

  g.add_edge(6, 7, 1.0);
  g.add_edge(6, 8, 6.0);

  g.add_edge(7, 8, 7.0);

  cout << g.mst() << endl << endl;
}

int main() {
  srand(time(0));

#ifdef TEST_SHORTEST_PATH
  test1_shortest_path();
  test2_shortest_path();
#endif

  Graph g50(50);

  cout << endl;
  cout << "Average shortest path for a Graph with 50 nodes, density 20%, range 10: " << g50.monte_carlo_simu(0.2, 10, 1000) << endl;
  cout << "Average shortest path for a Graph with 50 nodes, density 40%, range 10: " << g50.monte_carlo_simu(0.4, 10, 1000) << endl;
  cout << endl;

  //Graph gfile("data.txt");

#ifdef TEST_MST
  test1_mst();
#endif

  return 1;
}
