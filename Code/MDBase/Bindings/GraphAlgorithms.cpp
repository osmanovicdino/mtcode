void DFUtil(int i, vector1<bool> &visited, matrix<int> &adj, vector1<int> &lens)
{
    visited[i] = true;

    cout << i << " ";

    for (int j = 0; j < lens[i]; j++)
    {
        int p = adj(i, j);
        if (!visited[p])
        {
            DFUtil(p, visited, adj, lens);
        }
    }
}

void DFUtilstore(int i, vector1<bool> &visited, matrix<int> &adj, vector1<int> &lens, vector1<int> &indexes, int &iter, int &len_dem)
{
    visited[i] = true;

    //cout << i << " ";
    indexes[iter] = i;
    iter++;
    len_dem++;

    for (int j = 0; j < lens[i]; j++)
    {
        int p = adj(i, j);
        if (!visited[p])
        {
            DFUtilstore(p, visited, adj, lens, indexes, iter, len_dem);
        }
    }
}

vector1<int> ConnectedComponents(matrix<int> &adj, vector1<int> &lens, vector1<int> &indexes)
{
    //find connected components of adjacency graph
    //return the vector of all the demarkers between one component and another
    //index returns

    int n = lens.getsize();

    if (indexes.getsize() != n || adj.getnrows() != n)
        error("error in sizes of input vectors in connected components");

    int iter = 0;
    int len_dem = 0;

    vector1<int> demarkus(n);

    int num_con_comp = 0;
    vector1<bool> already_accounted(n);

    for (int i = 0; i < n; i++)
    {
        if (already_accounted[i] == false)
        {
            DFUtilstore(i, already_accounted, adj, lens, indexes, iter, len_dem);

            demarkus[num_con_comp] = len_dem;

            num_con_comp++;
            len_dem = 0;

            // cout << "\n";
        }
    }

    vector1<int> nbins(num_con_comp + 1);
    nbins[0] = 0;
    int runningtotal = 0;
    for (int i = 0; i < num_con_comp; i++)
    {
        runningtotal += demarkus[i];
        nbins[i + 1] = runningtotal;
    }

    return nbins;
}