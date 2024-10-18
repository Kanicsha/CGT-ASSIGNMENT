#include <iostream>
#include <vector>
#include<cstdlib>
#include<queue>
#include<utility>
#include <climits>
#include <queue>
#include <set>

using namespace std;

//removing first element
int rem(vector<int>& arr, int n, int pos) {
    arr[pos] = 0;
    return n - 1;
}

int isZero(vector<int>& arr, int n, int size) {
    for (int i = 0; i < n; i++) {
        if (arr[i] < 0) return -1;
        if (arr[i] != 0) return 0;
    }
    return 1;
}

int findStartVert(vector<vector<int> >& tempGraph,int NODE){
   for(int i = 0; i<NODE; i++){
      int deg = 0;
      for(int j = 0; j<NODE; j++){
         if(tempGraph[i][j])
         deg++; //increase degree, when connected edge found
      }
      if(deg % 2 != 0) //when degree of vertices are odd
      return i; //i is node with odd degree
   }
   return 0; //when all vertices have even degree, start from 0
}
bool isBridge(vector<vector<int> >& tempGraph,int u, int v, int NODE){
   int deg = 0;
   for(int i = 0; i<NODE; i++)
      if(tempGraph[v][i])
         deg++;
      if(deg>1){
         return false; //the edge is not forming bridge
      }
   return true; //edge forming a bridge
}
int edgeCount(vector<vector<int> >& tempGraph, int NODE){
   int count = 0;
   for(int i = 0; i<NODE; i++)
      for(int j = i; j<NODE; j++)
         if(tempGraph[i][j])
            count++;
   return count; //count nunber of edges in the graph
}
void fleuryAlgorithm(vector<vector<int> >& tempGraph,int start, int NODE){
   static int edge = edgeCount(tempGraph, NODE);
   for(int v = 0; v<NODE; v++){
      if(tempGraph[start][v]){ //when (u,v) edge is presnt and not forming bridge
         if(edge <= 1 || !isBridge(tempGraph,start, v,NODE)){
            cout << start << "--" << v << " ";
            tempGraph[start][v] = tempGraph[v][start] = 0; //remove edge from graph
            edge--; //reduce edge
            fleuryAlgorithm(tempGraph,v,NODE);
         }
      }
   }
}
void dfs(vector<vector<int> >& adj_matrix, int num_nodes, int start, bool *visited) {
    // Set the visit for starting node as true
    visited[start] = true;
    for (int i = 0; i < num_nodes; i++) {
        if (adj_matrix[start][i] == 1 && !visited[i]) {
            // Check if there's an edge between the two vertices and if it's visited or not
            dfs(adj_matrix, num_nodes, i, visited);
        }
    }
}

// Function to check if a graph is connected
bool is_connected_graph(vector<vector<int> >& adj_matrix, int num_nodes) {
    bool *visited = new bool[num_nodes];
    for (int i = 0; i < num_nodes; i++) {
        visited[i] = false;
    }

    dfs(adj_matrix, num_nodes, 0, visited);
    // After DFS if all the nodes aren't visited , then it's disconnected

    for (int i = 0; i < num_nodes; i++) {
        if (!visited[i]) {
            delete[] visited;
            return false; // Disconnected graph
        }
    }

    delete[] visited;
    return true; // Connected graph
}


vector<vector<int> > addWeight(vector<vector<int> >& adj_matrix, int n){
	int i,j,k=0,c=0;
	cout<<"BEFORE";
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			if(i<j){
			if(adj_matrix[i][j]==1) {
				cout<<"Enter weight in kgs:";
				cin>>adj_matrix[i][j];
				adj_matrix[j][i]=adj_matrix[i][j];
			}
		}
	}
	}
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			cout<<adj_matrix[i][j]<<"\t";
		}
		cout<<endl;
	}
	
	return adj_matrix;	
}
void sparseMatrix(int n,vector<vector<int> >& sparseMatrix){
	int i,j,k=0,c=0;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			if(sparseMatrix[i][j]==1) c++;
		}
	}
	vector<vector<int> > compactMatrix(3, vector<int>(c/2, 0));
	for ( i = 0; i < n; i++)
        for ( j = 0; j < n; j++)
            if (sparseMatrix[i][j] != 0 && i<=j)
            {
    			
                compactMatrix[0][k] = i;
                compactMatrix[1][k] = j;
                int randomNum = rand();

                compactMatrix[2][k] =22;
                k++;
            }
 
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<c/2; j++)
            printf("%d ", compactMatrix[i][j]);
 
        printf("\n");
    }
    
}
bool isEuler(vector<vector<int> >& adj_matrix, int num_nodes) {
    int i, degree;
    int odd_degree_count = 0;

    // Calculate the degree of each node
    for (i = 0; i < num_nodes; i++) {
        degree = 0;
        for (int j = 0; j < num_nodes; j++) {
            degree += adj_matrix[i][j];
        }
        if (degree % 2 != 0) {
            odd_degree_count++;
        }
    }
	
    // Check if the graph is Eulerian
    return (odd_degree_count == 0 || odd_degree_count == 2);
}
bool isEulerhh(vector<int> & hh, int n){
	int i;
	int even_deg=0;
	int odd_deg=0;
	for(i=0;i<n;i++){
		if(hh[i] & 1){
			cout<<"TOTO: "<<hh[i];
			odd_deg++;
		}
		else even_deg++;
	}
	cout<<"even::"<<even_deg<<endl<<"odd"<<odd_deg;
	
	if(odd_deg==2 || odd_deg== 0)
	{
			if(even_deg==n) 
	{
	cout<<"It's an Euler Circuit";

}
	if(odd_deg==2 && even_deg==n-2) {
		cout<<"It's an Euler Path";
}
	return true;
	
}
	return false;
	
}


// A utility function to find the vertex with minimum
// distance value, from the set of vertices not yet included
// in shortest path tree
int minDistance(vector<int> dist, vector<bool> sptSet,int V)
{

    // Initialize min value
    int min = INT_MAX, min_index;

    for (int v = 0; v < V; v++)
        if (sptSet[v] == false && dist[v] <= min)
            min = dist[v], min_index = v;

    return min_index;
}

// A utility function to print the constructed distance
// array
void printSolution(vector<int> dist, int V)
{
    cout << "Vertex \t Distance from Source" << endl;
    for (int i = 0; i < V; i++)
        cout << i << " \t\t\t\t" << dist[i] << endl;
}

// Function that implements Dijkstra's single source
// shortest path algorithm for a graph represented using
// adjacency matrix representation
void dijkstra(vector<vector<int> >& graph, int src, int V)
{
     vector<int> dist(V); // The output array.  dist[i] will hold the
                 // shortest
    // distance from src to i
 vector<bool> sptSet(V); // sptSet[i] will be true if vertex i is
                    // included in shortest
    // path tree or shortest distance from src to i is
    // finalized

    // Initialize all distances as INFINITE and stpSet[] as
    // false
    for (int i = 0; i < V; i++)
        dist[i] = INT_MAX, sptSet[i] = false;

    // Distance of source vertex from itself is always 0
    dist[src] = 0;

    // Find shortest path for all vertices
    for (int count = 0; count < V - 1; count++) {
        // Pick the minimum distance vertex from the set of
        // vertices not yet processed. u is always equal to
        // src in the first iteration.
        int u = minDistance(dist,sptSet,V);

        // Mark the picked vertex as processed
        sptSet[u] = true;

        // Update dist value of the adjacent vertices of the
        // picked vertex.
        for (int v = 0; v < V; v++)

            // Update dist[v] only if is not in sptSet,
            // there is an edge from u to v, and total
            // weight of path from src to  v through u is
            // smaller than current value of dist[v]
            if (!sptSet[v] && graph[u][v]
                && dist[u] != INT_MAX
                && dist[u] + graph[u][v] < dist[v])
                dist[v] = dist[u] + graph[u][v];
    }

    // print the constructed distance array
    printSolution(dist,V);
}
void selectionSort(vector<int>& arr, vector<int>& index, int n, int size) {
    for (int i = n - size; i < n - 1; i++) {
        int min_idx = i;
        for (int j = i + 1; j < n; j++) {
            if (arr[j] > arr[min_idx]) {
                min_idx = j;
            }
        }
        swap(arr[min_idx], arr[i]);
        swap(index[min_idx], index[i]);
    }
}

// Prim's algorithm

int minKey(int key[], bool mstSet[], int V)
{
    // Initialize min value
    int min = INT_MAX, min_index;

    for (int v = 0; v < V; v++)
        if (mstSet[v] == false && key[v] < min)
            min = key[v], min_index = v;

    return min_index;
}

// A utility function to print the
// constructed MST stored in parent[]
vector<pair<int, int> > printMST(int parent[],vector<vector<int> >& graph, int V)
{
	vector<vector<int> > spanning(V, vector<int>(V, 0));
		vector<vector<int> > spanningEdges(V-1, vector<int>(2, 0));
		vector<pair<int, int> > vec;
    cout << "Edge \tWeight\n";
    for (int i = 1; i < V; i++){
        cout << parent[i] << " - " << i << " \t"
             << graph[i][parent[i]] << " \n";
             spanning[parent[i]][i]=graph[i][parent[i]];
             spanning[i][parent[i]]=spanning[parent[i]][i];
              vec.push_back(make_pair(parent[i], i));
         }
    for(int i=0;i<V;i++){
    	for(int j=0;j<V;j++){
    		cout<< spanning[i][j]<<" ";
		}
		cout<<endl;
	}
	
for (size_t i = 0; i < vec.size(); ++i) {
    cout << "(" << vec[i].first << ", " << vec[i].second << ")" << endl;
}
return vec;
}
	/* Fundamental Cutsets*/


	/*End of Fundamental Cutsets*/

// Function to construct and print MST for
// a graph represented using adjacency
// matrix representation
vector<pair<int, int> > primMST(vector<vector<int> >& graph, int V)
{
    // Array to store constructed MST
    int parent[V];

    // Key values used to pick minimum weight edge in cut
    int key[V];

    // To represent set of vertices included in MST
    bool mstSet[V];

    // Initialize all keys as INFINITE
    for (int i = 0; i < V; i++)
        key[i] = INT_MAX, mstSet[i] = false;

    // Always include first 1st vertex in MST.
    // Make key 0 so that this vertex is picked as first
    // vertex.
    key[0] = 0;
  
    // First node is always root of MST
    parent[0] = -1;

    // The MST will have V vertices
    for (int count = 0; count < V - 1; count++) {
        
        // Pick the minimum key vertex from the
        // set of vertices not yet included in MST
        int u = minKey(key, mstSet,V);

        // Add the picked vertex to the MST Set
        mstSet[u] = true;

        // Update key value and parent index of
        // the adjacent vertices of the picked vertex.
        // Consider only those vertices which are not
        // yet included in MST
        for (int v = 0; v < V; v++)

            // graph[u][v] is non zero only for adjacent
            // vertices of m mstSet[v] is false for vertices
            // not yet included in MST Update the key only
            // if graph[u][v] is smaller than key[v]
            if (graph[u][v] && mstSet[v] == false
                && graph[u][v] < key[v])
                parent[v] = u, key[v] = graph[u][v];
    }

    // Print the constructed MST
    return printMST(parent, graph,V);
}
//    //
// Function to find fundamental circuits
vector<vector<int> > findFundamentalCircuits(vector<pair<int, int> >& mst_edges, vector<vector<int> >& graph, int V) {
    vector<vector<int> > circuits;

    // Find all edges not part of the MST (cotree edges)
    for (int i = 0; i < V; i++) {
        for (int j = i + 1; j < V; j++) {
            if (graph[i][j] != 0) {
                bool in_tree = false;

                // Check if the edge (i, j) is in the MST
                for (auto & edge:mst_edges) {
                    if ((edge.first == i && edge.second == j) || (edge.first == j && edge.second == i)) {
                        in_tree = true;
                        break;
                    }
                }

                // If edge (i, j) is not in the MST, it forms a fundamental circuit
                if (!in_tree) {
                    vector<int> circuit;
                    circuit.push_back(i);
                    circuit.push_back(j);

                    // Traverse the MST to find the path from i to j
                    // This path + edge (i, j) forms a cycle
                    vector<bool> visited(V, false);
                    vector<int> parent(V, -1);
                    queue<int> q;
                    q.push(i);
                    visited[i] = true;

                    while (!q.empty()) {
                        int u = q.front();
                        q.pop();

                        for (auto& edge : mst_edges) {
                            int v = -1;
                            if (edge.first == u)
                                v = edge.second;
                            else if (edge.second == u)
                                v = edge.first;

                            if (v != -1 && !visited[v]) {
                                visited[v] = true;
                                parent[v] = u;
                                q.push(v);

                                if (v == j)
                                    break;
                            }
                        }
                    }

                    // Reconstruct the path from i to j
                    int cur = j;
                    while (cur != i) {
                        circuit.push_back(cur);
                        cur = parent[cur];
                    }
                    circuit.push_back(i);

                    circuits.push_back(circuit);
                }
            }
        }
    }

    return circuits;
}

// Function to find fundamental cutsets
vector<vector<int> > findFundamentalCutsets(vector<pair<int, int> >& mst_edges, vector<vector<int> >& graph, int V) {
    vector<vector<int> > cutsets;

    // Fundamental cutsets are the edges between the two components
    for (int i = 0; i < V; i++) {
        for (int j = i + 1; j < V; j++) {
            if (graph[i][j] != 0) {
                bool in_tree = false;

                // Check if the edge (i, j) is in the MST
                for (auto & edge : mst_edges) {
                    if ((edge.first == i && edge.second == j) || (edge.first == j && edge.second == i)) {
                        in_tree = true;
                        break;
                    }
                }

                // If edge (i, j) is not in the MST, it defines a cutset
                if (!in_tree) {
                    vector<int> cutset;

                    // Find the connected components by removing the edge (i, j)
                    vector<bool> visited(V, false);
                    queue<int> q;
                    q.push(i);
                    visited[i] = true;

                    while (!q.empty()) {
                        int u = q.front();
                        q.pop();

                        for (auto& edge : mst_edges) {
                            int v = -1;
                            if (edge.first == u)
                                v = edge.second;
                            else if (edge.second == u)
                                v = edge.first;

                            if (v != -1 && !visited[v]) {
                                visited[v] = true;
                                q.push(v);
                            }
                        }
                    }

                    for (int k = 0; k < V; k++) {
                        if (visited[k])
                            cutset.push_back(k);
                    }

                    cutsets.push_back(cutset);
                }
            }
        }
    }

    return cutsets;
}
int main() {
    int n, ele, size;
    
	cout << "Enter number of vertices: ";
    cin >> n;
    size = n;

    vector<int> arr(n), index(n);
    for (int j = 0; j < n; j++) {
        index[j] = j;
    }

    vector<vector<int> > graph(n, vector<int>(n, 0));
    vector<vector<int> > graph2(n, vector<int>(n, 0));
    vector<vector<int> > graph3(n, vector<int>(n, 0));
    vector<pair<int, int> > vec;
   
    
    

    for (int i = 0; i < n; i++) {
        cout << "Enter degree: ";
        cin >> arr[i];
    }

    while (isZero(arr, n, size) == 0) {
        ele = arr[n - size];
        size = rem(arr, size, n - size);

        for (int j = 0; j < ele; j++) {
            arr[n - size + j] -= 1;
            graph[index[n - size - 1]][index[n - size + j]] = 1;
            graph[index[n - size + j]][index[n - size - 1]] = 1;
        }

        selectionSort(arr, index, n, size);
        cout << endl;
    }

    cout << "\nAdjacency Matrix:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << graph[i][j] << " ";
            graph2[i][j]=graph[i][j];
            graph3[i][j]=graph[i][j];
        }
        cout << endl;
    }
    
	if (isZero(arr, n, size) == 1){
        cout << "Sequence is graphical" << endl;
        cout<<"Using Adjacency matrix algo"<<endl;
       // if(isEuler(graph, n)) cout<<"IT's eulerian graph"<<endl;
        cout<<"Using Havel Hakimi algo"<<endl;
        isEulerhh(arr,n);
        if(1){
     	fleuryAlgorithm(graph,findStartVert(graph,n), n);
     	cout<<endl;
     	graph2=addWeight(graph2,n);
     	  for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << graph[i][j] << " ";
            graph3[i][j]=graph2[i][j];
        }
        cout << endl;
    }
     	int src;
     	cout<<"Enter source:";
     	cin>>src;
     	if(src<n)
     	dijkstra(graph2, src,n);
        cout << endl;
        vec=primMST(graph2,n);
         vector<vector<int> > circuits = findFundamentalCircuits(vec, graph3, n);
     vector<vector<int> > cutsets = findFundamentalCutsets(vec, graph3, n);
         cout << "Fundamental Circuits:\n";
    for (auto& circuit : circuits) {
        for (auto node : circuit) {
            cout << node << " ->";
        }
        cout << endl;
    }

    cout << "Fundamental Cutsets:\n";
    for (auto & cutset : cutsets) {
        for (auto node : cutset) {
            cout << node << " ->";
        }
        cout << endl;
    }
    }
}
        
    else
        {
		cout << "Cannot proceed further" << endl;
        sparseMatrix(n,graph);
}

    return 0;
	}
