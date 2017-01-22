// g++ -Wall -std=c++11 dijkstra.cpp

#include <iostream>
#include <map> // header file needed for to use MAP STL
#include <queue>

using namespace std;

// Keypoint here is to understand the structure used to represent a graph:
// a map (nodes) of map (neighbors with associated distance)
// map<char, map<char,int>> g; // like a python dict (dict of dict)

class Graph {
public:
  void add_node(char node) 
  {
     g[node]; // empty initialization
  }

  void add_edge(char from, char to, int distance) {
    g[from][to] = distance;
    g[to];
  }

  int get_edge_value(char from, char to) 
  {
    return g[from][to];
  }

  void delete_edge(char from, char to) 
  {
    g[from].erase(to);
  }

  int shortest_path(char src, char dest) 
  {
    if (g.count(src) == 0 || g.count(dest) == 0) 
    {
      return -1;
    }

    // closed_set: nodes that have known shortest distance
    map<char, pair<int, char>> closed_set;
    closed_set[src] = make_pair(0, src); // with map: node -> dist, prev

    // open_set: what is reachable from closed_set nodes (nodes under exploration)
    map<char, pair<int, char>> open_set;

    // int (dist from src) char (previous node) char (node)
    typedef tuple<int, char, char> icc;
    // use a min priority queue for better efficiency of open_set handling
    priority_queue<icc, vector<icc>, greater<icc>> open_set_pq;

    char expand_node = src;
    int expand_dist = 0;
    int done = 0;
    int new_dist = -1;

    int step = 1;

    while (!done) 
    {
      cout << "step: " << step << endl;

      for (auto iter : g[expand_node]) 
      {
        auto neighbor = iter.first;
        auto dist = iter.second;
        if (closed_set.count(neighbor) == 0) // if neighbor not in closed_set
        {
          new_dist = expand_dist + dist;
          if (open_set.count(neighbor) == 0 || new_dist < open_set[neighbor].first)
          {
            open_set[neighbor] = make_pair(new_dist, expand_node);
            open_set_pq.push(make_tuple(new_dist, expand_node, neighbor));
          }
        }
      } 

      icc expand = open_set_pq.top();
      expand_node = get<2>(expand);
      expand_dist = get<0>(expand);
      // put it in closed set and delete it from open_set
      closed_set[expand_node] = make_pair(expand_dist, get<1>(expand));
      open_set_pq.pop();

      if (expand_node == dest)
      {
        done = 1;
        // // show content:
        for (auto node : closed_set)
        {
          cout << "shortest distance from " << src << " for node " << node.first << " = " << node.second.first;
          cout << " (previous node: " << node.second.second << ")\n";
        }
        return expand_dist;
      }

      new_dist = -1;
      step += 1;
    }

    return 0;
  }

  friend ostream& operator<<(ostream& out, const Graph& graph) 
  {
    cout << "Graph adjacency list\n";
    for (auto node : graph.g) // pythonic for
    {
      out << "\tnode " << node.first << " => ";
      for (auto neighbor : node.second) // pythonic for
      {
        out << node.first << ":" << neighbor.first << "=" << neighbor.second << " ";
      }
      out << endl;
    }
    return out;
  }

private:
  map<char, map<char,int>> g; // like a python dict (dict of dict)
};




int main(void)
{
  Graph g;

  g.add_node('Q');
  g.add_node('T');

  g.add_edge('P', 'Q', 666);
  g.delete_edge('P', 'Q');
  g.add_edge('S', 'A', 4);
  g.add_edge('S', 'D', 7);
  g.add_edge('S', 'B', 3);
  
  g.add_edge('A', 'C', 1);
  
  g.add_edge('B', 'S', 3);
  g.add_edge('B', 'D', 4);
  
  g.add_edge('C', 'D', 3);
  g.add_edge('C', 'E', 1);
  
  g.add_edge('D', 'E', 1);
  g.add_edge('D', 'T', 3);
  g.add_edge('D', 'F', 5);
  
  g.add_edge('E', 'G', 2);
  g.add_edge('E', 'T', 4);
  
  g.add_edge('G', 'E', 2);
  g.add_edge('G', 'T', 3);

  g.add_node('G');

  cout << g;
  //cout << "number_of_nodes: " << g.number_of_nodes() << endl;
  //cout << "number_of_edges: " << g.number_of_edges() << endl;
  //cout << g.adjacent('A', 'D') << endl;

  cout << g.shortest_path('S', 'T') << endl;

  return 0;
}
