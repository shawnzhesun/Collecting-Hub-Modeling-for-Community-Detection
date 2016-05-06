#include <metis.h>

/* merge heterogeneous graphs into a single representation */
/* Yiye Ruan */
GraphType *mergeHeteroGraphs (GraphType * g1, GraphType * g2, int *heteroSource, int addSelfLoop)
{
    GraphType *graph = (GraphType *) malloc (sizeof (GraphType));
    InitGraph (graph);
    graph->nvtxs = g1->nvtxs;
    graph->vwgt = graph->adjwgt = NULL;
    idxtype *xadj = graph->xadj =
        idxsmalloc (graph->nvtxs + 1, 0, "MergeHeteroGraphs: xadj");
    /* assume that all adjacency lists have already been sorted
     * (except for self loop, if present)
     * merging multiple adjacency lists then becomes easier */
    idxtype *initAdjncy =
        idxmalloc (g1->nedges + g2->nedges, "MergeHeteroGraphs: init adjncy");

    idxtype realEdgeCount = 0;
    idxtype vtxId;
    for (vtxId = 0; vtxId < graph->nvtxs; vtxId++)
    {
        idxtype g1Index = g1->xadj[vtxId], g1End =
            g1->xadj[vtxId + 1], g2Index = g2->xadj[vtxId], g2End =
            g2->xadj[vtxId + 1];
        if (addSelfLoop == 1)
        {
            /* the first entry in adjncy[xadh[vtxId]] must be itself */
            assert (g1->adjncy[g1Index] == g2->adjncy[g2Index] && g1->adjncy[g1Index] == vtxId);
            initAdjncy[realEdgeCount++] = vtxId;
            g1Index++;
            g2Index++;
        }
        for (; g1Index < g1End && g2Index < g2End;)
        {
            idxtype diff = g1->adjncy[g1Index] - g2->adjncy[g2Index];
            if (diff == 0)
            {
                initAdjncy[realEdgeCount++] = g1->adjncy[g1Index];
                g1Index++;
                g2Index++;
            }
            else if (diff < 0)
                initAdjncy[realEdgeCount++] = g1->adjncy[g1Index++];
            else
                initAdjncy[realEdgeCount++] = g2->adjncy[g2Index++];
        }
        if (g1Index == g1End)
            for (; g2Index < g2End;)
                initAdjncy[realEdgeCount++] = g2->adjncy[g2Index++];
        else
            for (; g1Index < g1End;)
                initAdjncy[realEdgeCount++] = g1->adjncy[g1Index++];

        xadj[vtxId + 1] = realEdgeCount;
    }

    graph->nedges = realEdgeCount;
    idxtype *adjncy = graph->adjncy =
        idxmalloc (graph->nedges, "MergeHeteroGraphs: adjncy");
    memcpy (adjncy, initAdjncy, sizeof (*adjncy) * graph->nedges);
    free (initAdjncy);
    initAdjncy = NULL;

    return graph;
}

/* read node-associated lists from file */
/* Yiye Ruan */
void readList (GraphType * list, const char *fileName, idxtype offset)
{
    FILE *fpin;
    InitGraph (list);
    char *line = (char *) malloc (sizeof (char) * (MAXLINE + 1));
    if ((fpin = fopen (fileName, "r")) == NULL)
    {
        fprintf (stderr, "Failed to open file %s\n", fileName);
        exit (1);
    }

    /* the first line has number of nodes and total number of terms */
    fgets (line, MAXLINE, fpin);
    sscanf (line, "%d %d", &(list->nvtxs), &(list->nedges));
    idxtype *xadj = list->xadj =
        idxsmalloc (list->nvtxs + 1, 0, "ContentList: xadj");
    idxtype *adjncy = list->adjncy =
        idxmalloc (list->nedges, "ContentList: adjncy");
    idxtype vtxId, adjIndex;
    /* read consequent lines */
    for (vtxId = 0, adjIndex = 0; vtxId < list->nvtxs; vtxId++)
    {
        fgets (line, MAXLINE, fpin);
        if (strlen (line) == MAXLINE)
        {
            fprintf (stderr, "Buffer for fgets not big enough\n");
            exit (1);
        }
        char *oldStr = line, *newStr = NULL;
        idxtype termId;
        while ((termId = (idxtype) strtol (oldStr, &newStr, 10)) != 0)
        {
            /*adjncy[adjIndex++] = termId - 1 + offset;*/
            adjncy[adjIndex++] = termId - 1;
            oldStr = newStr;
        }
        xadj[vtxId + 1] = adjIndex;
    }

    assert (adjIndex == list->nedges);

    fclose (fpin);
    free (line);
    line = NULL;
}

/* estimate similarity between unions of multiple bags (e.g. adjacency
 * list, concept list) for nodes */
/* Yiye Ruan */

GraphType *getSingleBagSparsifiedGraph (GraphType * mergeGraph, GraphType * list1, GraphType * list2, int numHashes, float degreeExponent)
{
    int *randoms = (int *) malloc (sizeof (int) * 2 * numHashes);
    generateRandoms (2 * numHashes, randoms);

    Hashtable *ht = buildSingleBagHashtable (list1, list2, numHashes, randoms);
    free (randoms);
    randoms = NULL;

    printf ("Done minwise hashing.\n");
    fflush (stdout);

#ifdef DEBUG
    // check the result of minwise hashing
    for (int myIter = 0; myIter < 2; myIter++)
    {
        for (int myIndex = 0; myIndex < numHashes; myIndex++)
        {
            printf ("%d ", ht->hashes[myIter * numHashes + myIndex]);
        }
        printf ("\n");
    }
#endif

    GraphType *retGraph = (GraphType *) malloc (sizeof (GraphType));
    retGraph->nvtxs = mergeGraph->nvtxs;
    retGraph->nedges = mergeGraph->nedges;
    retGraph->xadj =
        (idxtype *) malloc (sizeof (idxtype) * (mergeGraph->nvtxs + 1));
    retGraph->adjncy =
        (idxtype *) malloc (sizeof (idxtype) * (mergeGraph->nedges));
    retGraph->adjwgt =
        (idxtype *) malloc (sizeof (idxtype) * (mergeGraph->nedges));

    int i, j, k;

    j = 0;
    for (i = 0; i < retGraph->nvtxs; i++)
    {
        retGraph->xadj[i] = j;
        for (k = mergeGraph->xadj[i]; k < mergeGraph->xadj[i + 1]; k++)
        {
            int numMatched = 0, l;
            int from = i, to = mergeGraph->adjncy[k];
            if (from == to)
            {
                // we won't consider self-loops
                continue;
            }
            idxtype *fromBase, *toBase;
            fromBase = ht->hashes + from * numHashes;
            toBase = ht->hashes + to * numHashes;
            for (l = 0; l < numHashes; l++)
            {
                if (fromBase[l] == toBase[l])
                    numMatched++;
            }
            retGraph->adjwgt[j] = numMatched;
            retGraph->adjncy[j] = to;
            j++;
        }
    }
    retGraph->xadj[retGraph->nvtxs] = j;

    printf ("Done estimating similarities for all edges.\n");
    //      printf("Total edges in sim. weighted graph:%d\n",
    //      retGraph->xadj[retGraph->nvtxs]);
    /*	if ( checkValidUndirectedGraph(retGraph) )
        printf("Sim. weighted graph is a valid undirected graph\n");
        else
        printf("Sim. weighted graph is a valid undirected graph\n");
        */
    fflush (stdout);

    retainTopNeighborsPerNode (retGraph, ht, degreeExponent);
    //      printf("Total edges in local sparsified graph:%d\n",
    //      retGraph->xadj[retGraph->nvtxs]);

    freeHashtable (ht);
    //      exactPruneGraph(retGraph, 0, threshold);
    return retGraph;
}

Hashtable *buildSingleBagHashtable (GraphType * list1, GraphType * list2, int numHashes, int *randoms)
{
    int i;
    Hashtable *ht = (Hashtable *) malloc (sizeof (Hashtable));
    idxtype *base;

    ht->numNodes = list1->nvtxs;
    ht->numHashes = numHashes;
    ht->sortedNodeIds = idxmalloc (list1->nvtxs,
            "buildHashtable:ht->sortedNodeIds");
    ht->hashes = idxmalloc (list1->nvtxs * numHashes,
            "buildHashtable:ht->hashes");

    int *mins = (int *) malloc (sizeof (int) * numHashes);

    base = ht->hashes;
    idxtype *adjList1 = list1->adjncy, *adjList2 = list2->adjncy;
    for (int j = 0; j < list1->nvtxs; j++)
    {
        int adjLen1 = list1->xadj[j + 1] - list1->xadj[j],
            adjLen2 = list2->xadj[j + 1] - list2->xadj[j];
        int adjlistLength = adjLen1 + adjLen2;
        idxtype *adjlist =
            (idxtype *) malloc (sizeof (*adjlist) * adjlistLength);
        memcpy (adjlist, adjList1, sizeof (*adjlist) * adjLen1);
        memcpy (adjlist + adjLen1, adjList2, sizeof (*adjlist) * adjLen2);
        if (adjlistLength == 1)
        {
            for (int i = 0; i < numHashes; i++)
                base[i] = adjlist[0];
        }
        else
        {
            getMinHashKeys (adjlist, adjlistLength, randoms,
                    numHashes, base, mins, numHashes);
        }

        base += numHashes;
        adjList1 += adjLen1;
        adjList2 += adjLen2;
        free (adjlist);
        adjlist = NULL;
    }

    free (mins);
    mins = NULL;
    /*	for ( i=0; i<graph->nvtxs; i++ )
        {
        ht->sortedNodeIds[i] = i;
        }
        */
    return ht;
}

GraphType *getSimWeightedSparsifiedGraph_cos (GraphType * mergeGraph, GraphType * list1, GraphType * list2, int numHash, float degreeExponent, double firstWeight, int normMode)
{
    /*
     * Construct LSH for cosine similarity
     */
    chunktype *ht1, *ht2;
    if (numHash != 0)
    {
        assert (numHash % (sizeof(chunktype) * 8) == 0);
        ht1 = generateHashTable(list1, numHash);
        ht2 = generateHashTable(list2, numHash);
    }

    GraphType *retGraph = (GraphType *) malloc (sizeof (GraphType));
    retGraph->nvtxs = mergeGraph->nvtxs;
    retGraph->nedges = mergeGraph->nedges;
    retGraph->xadj =
        (idxtype *) malloc (sizeof (idxtype) * (mergeGraph->nvtxs + 1));
    retGraph->adjncy =
        (idxtype *) malloc (sizeof (idxtype) * (mergeGraph->nedges));
    /*retGraph->adjwgt = (idxtype*)malloc(sizeof(idxtype)*(mergeGraph->nedges)); */

    /* where we temporarily store integrated similarity */
    IntDoublePair *tempPair =
        (IntDoublePair *) malloc (sizeof (*tempPair) * (mergeGraph->nedges));
    memset (tempPair, 0, sizeof (*tempPair) * (mergeGraph->nedges));

    int i, j, k;

    j = 0;
    for (i = 0; i < retGraph->nvtxs; i++)
    {
        retGraph->xadj[i] = j;

        double min1 = DBL_MAX, max1 = DBL_MIN, min2 = DBL_MAX, max2 = DBL_MIN;
        double *sim1 = (double *) malloc (sizeof (*sim1) * (mergeGraph->xadj[i+1] - mergeGraph->xadj[i]));
        double *sim2 = (double *) malloc (sizeof (*sim2) * (mergeGraph->xadj[i+1] - mergeGraph->xadj[i]));
        memset(sim1, 0, sizeof(*sim1) * (mergeGraph->xadj[i+1] - mergeGraph->xadj[i]));
        memset(sim2, 0, sizeof(*sim2) * (mergeGraph->xadj[i+1] - mergeGraph->xadj[i]));

        int jj = j;

        for (k = mergeGraph->xadj[i]; k < mergeGraph->xadj[i + 1]; k++)
        {
            int from = i, to = mergeGraph->adjncy[k];
            if (from == to)
                continue;

            double cosSim1, cosSim2;
            if (numHash != 0)
            {
                cosSim1 = cosSimByLSH(from, to, ht1, list1, numHash, numHash / (sizeof(chunktype) * 8));
                cosSim2 = cosSimByLSH(from, to, ht2, list2, numHash, numHash / (sizeof(chunktype) * 8));
            }
            else
            {
                cosSim1 = compCosSim(from, to, list1->adjncy, list1->xadj[from], list1->xadj[from+1], list1->xadj[to], list1->xadj[to+1]);
                cosSim2 = compCosSim(from, to, list2->adjncy, list2->xadj[from], list2->xadj[from+1], list2->xadj[to], list2->xadj[to+1]);
            }

            min1 = (cosSim1 < min1) ? cosSim1 : min1;
            max1 = (cosSim1 > max1) ? cosSim1 : max1;
            min2 = (cosSim2 < min2) ? cosSim2 : min2;
            max2 = (cosSim2 > max2) ? cosSim2 : max2;
            sim1[k - mergeGraph->xadj[i]] = cosSim1;
            sim2[k - mergeGraph->xadj[i]] = cosSim2;

            tempPair[j].id = to;
            j++;
        }

        double location1 = 0, location2 = 0, scale1 = 0, scale2 = 0;
        compLocationAndScale_double(i, sim1, mergeGraph->xadj[i], mergeGraph->xadj[i+1], mergeGraph->adjncy, normMode, min1, max1, &location1, &scale1);
        compLocationAndScale_double(i, sim2, mergeGraph->xadj[i], mergeGraph->xadj[i+1], mergeGraph->adjncy, normMode, min2, max2, &location2, &scale2);
#ifdef DEBUG
        printf
            ("min1 %f max1 %f min2 %f max2 %f location1 %f scale1 %f location2 %f scale2 %f\n",
             min1, max1, min2, max2, location1, scale1, location2, scale2);
#endif

        for (k = mergeGraph->xadj[i]; k < mergeGraph->xadj[i + 1]; k++)
        {
            if (mergeGraph->adjncy[k] == i)
                continue;
            tempPair[jj].value +=
                (sim1[k - mergeGraph->xadj[i]] - location1) / scale1 * firstWeight + 
                (sim2[k - mergeGraph->xadj[i]] - location2) / scale2 * (1 - firstWeight);
            jj++;
        }

        free (sim1);
        sim1 = NULL;
        free (sim2);
        sim2 = NULL;

#ifdef DEBUG
        printf("Done for node %d\n", i);
#endif

    }
    retGraph->xadj[retGraph->nvtxs] = j;
    if (numHash != 0)
    {
        free(ht1); 
        ht1 = NULL;
        free(ht2);
        ht2 = NULL;
    }

    printf ("Done estimating similarities for all edges.\n");
    //      printf("Total edges in sim. weighted graph:%d\n",
    //      retGraph->xadj[retGraph->nvtxs]);
    /*	if ( checkValidUndirectedGraph(retGraph) )
        printf("Sim. weighted graph is a valid undirected graph\n");
        else
        printf("Sim. weighted graph is a valid undirected graph\n");
        */
    fflush (stdout);

    retainSimWeightedTopNeighbors (retGraph, tempPair, degreeExponent);
    /*retainTopNeighborsPerNode(retGraph, ht, degreeExponent); */

    //      printf("Total edges in local sparsified graph:%d\n",
    //      retGraph->xadj[retGraph->nvtxs]);

    free (tempPair);
    tempPair = NULL;
    //      exactPruneGraph(retGraph, 0, threshold);
    return retGraph;
}

double compCosSim(const idxtype id1, const idxtype id2, const idxtype *adjncy, const idxtype startIdx1, const idxtype endIdx1, const idxtype startIdx2, const idxtype endIdx2)
{
    idxtype length1 = endIdx1 - startIdx1, length2 = endIdx2 - startIdx2;
    idxtype curIdx1 = startIdx1, curIdx2 = startIdx2;
    idxtype match = 0;
    for ( ; curIdx1 < endIdx1 && curIdx2 < endIdx2; )
    {
        if (adjncy[curIdx1] == id1 || adjncy[curIdx1] == id2)
        {
            curIdx1 ++;
            length1 --;
            continue;
        }
        if (adjncy[curIdx2] == id1 || adjncy[curIdx2] == id2)
        {
            curIdx2 ++;
            length2 --;
            continue;
        }
        idxtype diff = adjncy[curIdx1] - adjncy[curIdx2];
        if (diff == 0)
        {
            match ++;
            curIdx1 ++;
            curIdx2 ++;
        }
        else if (diff < 0)
            curIdx1 ++;
        else 
            curIdx2 ++;
    }

    if (match == 0 || length1 == 0 || length2 == 0)
        return 0.0;

    return match / sqrt(length1 * length2 * 1.0);
}

GraphType *getSimWeightedSparsifiedGraph_jac (GraphType * mergeGraph, GraphType * list1, GraphType * list2, int numHashes, float degreeExponent, double firstWeight, int normMode)
{
    int *randoms = (int *) malloc (sizeof (int) * 2 * numHashes);
    generateRandoms (2 * numHashes, randoms);

    Hashtable *ht1 = buildZeroFriendlyHashtable (list1, numHashes, randoms);
    Hashtable *ht2 = buildZeroFriendlyHashtable (list2, numHashes, randoms);
    free (randoms);
    randoms = NULL;

    printf ("Done minwise hashing.\n");
    fflush (stdout);

#ifdef DEBUG
    // check the result of minwise hashing
    for (int myIter = 0; myIter < 3; myIter++)
    {
        printf ("hashtable 1:\n");
        for (int myIndex = 0; myIndex < numHashes; myIndex++)
        {
            printf ("%d ", ht1->hashes[myIter * numHashes + myIndex]);
        }
        printf ("\n");
        printf ("hashtable 2:\n");
        for (int myIndex = 0; myIndex < numHashes; myIndex++)
        {
            printf ("%d ", ht2->hashes[myIter * numHashes + myIndex]);
        }
        printf ("\n");
    }
#endif

    GraphType *retGraph = (GraphType *) malloc (sizeof (GraphType));
    retGraph->nvtxs = mergeGraph->nvtxs;
    retGraph->nedges = mergeGraph->nedges;
    retGraph->xadj =
        (idxtype *) malloc (sizeof (idxtype) * (mergeGraph->nvtxs + 1));
    retGraph->adjncy =
        (idxtype *) malloc (sizeof (idxtype) * (mergeGraph->nedges));
    /*retGraph->adjwgt = (idxtype*)malloc(sizeof(idxtype)*(mergeGraph->nedges)); */

    /* where we temporarily store integrated similarity */
    IntDoublePair *tempPair =
        (IntDoublePair *) malloc (sizeof (*tempPair) * (mergeGraph->nedges));
    memset (tempPair, 0, sizeof (*tempPair) * (mergeGraph->nedges));

    int i, j, k;

    j = 0;
    for (i = 0; i < retGraph->nvtxs; i++)
    {
        retGraph->xadj[i] = j;

        int min1 = INT_MAX, max1 = INT_MIN, min2 = INT_MAX, max2 = INT_MIN;
        int *match1 = (int *) malloc (sizeof (*match1) * (mergeGraph->xadj[i + 1] - mergeGraph->xadj[i]));
        int *match2 = (int *) malloc (sizeof (*match2) * (mergeGraph->xadj[i + 1] - mergeGraph->xadj[i]));
        memset(match1, 0, sizeof(*match1) * (mergeGraph->xadj[i+1] - mergeGraph->xadj[i]));
        memset(match2, 0, sizeof(*match2) * (mergeGraph->xadj[i+1] - mergeGraph->xadj[i]));

        int jj = j;

        for (k = mergeGraph->xadj[i]; k < mergeGraph->xadj[i + 1]; k++)
        {
            int l;
            idxtype numMatched1 = 0, numMatched2 = 0;
            int from = i, to = mergeGraph->adjncy[k];
            if (from == to)
            {
                // we won't consider self-loops
                continue;
            }
            idxtype *fromBase, *toBase;
            fromBase = ht1->hashes + from * numHashes;
            toBase = ht1->hashes + to * numHashes;
            for (l = 0; l < numHashes; l++)
            {
                if (fromBase[l] == toBase[l])
                    numMatched1++;
            }
            fromBase = ht2->hashes + from * numHashes;
            toBase = ht2->hashes + to * numHashes;
            for (l = 0; l < numHashes; l++)
            {
                if (fromBase[l] == toBase[l])
                    numMatched2++;
            }

            min1 = (numMatched1 < min1) ? numMatched1 : min1;
            max1 = (numMatched1 > max1) ? numMatched1 : max1;
            min2 = (numMatched2 < min2) ? numMatched2 : min2;
            max2 = (numMatched2 > max2) ? numMatched2 : max2;
            match1[k - mergeGraph->xadj[i]] = numMatched1;
            match2[k - mergeGraph->xadj[i]] = numMatched2;

            /*!!!NEXT LINE DEPRECATED!!! */
            /*tempPair[j].value = firstWeight * numMatched1 + (1-firstWeight) * numMatched2; */
            tempPair[j].id = to;
            /*retGraph->adjwgt[j] = numMatched; */
            /*retGraph->adjncy[j] = to; */
            j++;
        }

        double location1 = 0, location2 = 0, scale1 = 0, scale2 = 0;
        compLocationAndScale_int(i, match1, mergeGraph->xadj[i], mergeGraph->xadj[i+1], mergeGraph->adjncy, normMode, min1, max1, &location1, &scale1);
        compLocationAndScale_int(i, match2, mergeGraph->xadj[i], mergeGraph->xadj[i+1], mergeGraph->adjncy, normMode, min2, max2, &location2, &scale2);
#ifdef DEBUG
        printf
            ("min1 %d max1 %d min2 %d max2 %d location1 %f scale1 %f location2 %f scale2 %f\n",
             min1, max1, min2, max2, location1, scale1, location2, scale2);
#endif

        for (k = mergeGraph->xadj[i]; k < mergeGraph->xadj[i + 1]; k++)
        {
            if (mergeGraph->adjncy[k] == i)
                continue;
            tempPair[jj].value +=
                (match1[k - mergeGraph->xadj[i]] - location1) / scale1 * firstWeight + 
                (match2[k - mergeGraph->xadj[i]] - location2) / scale2 * (1 - firstWeight);
            jj++;
        }

        free (match1);
        match1 = NULL;
        free (match2);
        match2 = NULL;

#ifdef DEBUG
        printf("Done for node %d\n", i);
#endif
    }
    retGraph->xadj[retGraph->nvtxs] = j;

    printf ("Done estimating similarities for all edges.\n");
    //      printf("Total edges in sim. weighted graph:%d\n",
    //      retGraph->xadj[retGraph->nvtxs]);
    /*	if ( checkValidUndirectedGraph(retGraph) )
        printf("Sim. weighted graph is a valid undirected graph\n");
        else
        printf("Sim. weighted graph is a valid undirected graph\n");
        */
    fflush (stdout);

    retainSimWeightedTopNeighbors (retGraph, tempPair, degreeExponent);
    /*retainTopNeighborsPerNode(retGraph, ht, degreeExponent); */

    //      printf("Total edges in local sparsified graph:%d\n",
    //      retGraph->xadj[retGraph->nvtxs]);

    freeHashtable (ht1);
    freeHashtable (ht2);
    free (tempPair);
    tempPair = NULL;
    //      exactPruneGraph(retGraph, 0, threshold);
    return retGraph;
}

/*Compute empirical location and scale of a series*/
void compLocationAndScale_double(const idxtype nodeId, const double *series, const idxtype xadjStart, const idxtype xadjEnd, const idxtype *adjncy, const idxtype normMode, const double min, const double max, double *location, double *scale){
#ifdef DEBUG
    printf("nodeId %d xadjStart %d xadjEnd %d normMode %d min %f max %f\n", nodeId, xadjStart, xadjEnd, normMode, min, max);
    int tempIdx;
    printf("Series: ");
    for (tempIdx = 0; tempIdx < xadjEnd - xadjStart; tempIdx++)
        printf("%f ", series[tempIdx]);
    printf("\n");
#endif
    if (min == max)
    {
        *location = min;
        *scale = 1;
    }
    else
    {
        assert (min < max);
        /*rescale to [0, 1]*/
        if (normMode == 1)
        {
            *location = min;
            *scale = max - min;
        }
        /*zero-mean-unit-variance*/
        else
        {
            double linSum = 0.0, sqSum = 0.0;
            double val;
            idxtype elementCount = xadjEnd - xadjStart;
            idxtype idx;
            for (idx = xadjStart; idx < xadjEnd; idx++)
            {
                if (adjncy[idx] == nodeId)
                {
                    elementCount --;
                    continue;
                }

                val = series[idx - xadjStart];
                linSum += val;
                sqSum += val * val;
            }
            *location = linSum / elementCount;
            *scale = sqrt((sqSum - linSum * linSum / elementCount) / (elementCount - 1));
        }
    }

#ifdef DEBUG
    printf("Location %f Scale %f\n", *location, *scale);
#endif

}

/*Compute empirical location and scale of a series*/
void compLocationAndScale_int(const idxtype nodeId, const idxtype *series, const idxtype xadjStart, const idxtype xadjEnd, const idxtype *adjncy, const idxtype normMode, const idxtype min, const idxtype max, double *location, double *scale){
#ifdef DEBUG
    printf("nodeId %d xadjStart %d xadjEnd %d normMode %d min %d max %d\n", nodeId, xadjStart, xadjEnd, normMode, min, max);
    int tempIdx;
    printf("Series: ");
    for (tempIdx = 0; tempIdx < xadjEnd - xadjStart; tempIdx++)
        printf("%d ", series[tempIdx]);
    printf("\n");
#endif
    if (min == max)
    {
        *location = min;
        *scale = 1;
    }
    else
    {
        assert (min < max);
        /*rescale to [0, 1]*/
        if (normMode == 1)
        {
            *location = min;
            *scale = max - min;
        }
        /*zero-mean-unit-variance*/
        else {
            idxtype linSum = 0, sqSum = 0, elementCount = xadjEnd - xadjStart;
            idxtype idx, matchCount;
            for (idx = xadjStart; idx < xadjEnd; idx++)
            {
                if (adjncy[idx] == nodeId)
                {
                    elementCount --;
                    continue;
                }

                matchCount = series[idx - xadjStart];
                linSum += matchCount;
                sqSum += matchCount * matchCount;
            }
            *location = linSum * 1.0 / elementCount;
            *scale = sqrt((sqSum - linSum * linSum * 1.0 / elementCount) / (elementCount - 1));
        }
    }
#ifdef DEBUG
    printf("Location %f Scale %f\n", *location, *scale);
#endif
}

/* compute the size of two sets' union */
idxtype getSetUnionSize (const idxtype * s1Start, const idxtype * s1End,
        const idxtype * s2Start, const idxtype * s2End)
{
    /*idxtype us1, us2; */
    /*us1 = getSetUnionSize(list1->adjncy + list1->xadj[from], list1->adjncy + list1->xadj[from+1], list1->adjncy + list1->xadj[to], list1->adjncy + list1->xadj[to+1]); */
    /*us2 = getSetUnionSize(list2->adjncy + list2->xadj[from], list2->adjncy + list2->xadj[from+1], list2->adjncy + list2->xadj[to], list2->adjncy + list2->xadj[to+1]); */

    /*if ( us1 != 0 ) */
    /*tempPair[j].value = firstWeight * numMatched1 + (1-firstWeight) * (us1 * 1.0 / us2) * numMatched2; */
    /*else if (us1 == 0) */
    /*{ */
    /*fprintf(stderr, "Puff, us1 is 0! Make sure addSelfLoop is on and try again"); */
    /*exit(1); */
    /*} */
    idxtype unionSize = 0;
    const idxtype *s1Index = s1Start, *s2Index = s2Start;
    for (; s1Index < s1End && s2Index < s2End;)
    {
        idxtype diff = *s1Index - *s2Index;
        unionSize++;
        if (diff == 0)
        {
            s1Index++;
            s2Index++;
        }
        else if (diff < 0)
            s1Index++;
        else
            s2Index++;
    }
    unionSize += (s1End - s1Index) + (s2End - s2Index);

    return unionSize;
}

/* modify Venu's buildHashtable function 
 * so that even if list for a vertex has length 0,
 * we can still deal with it */
Hashtable *buildZeroFriendlyHashtable (GraphType * graph, int numHashes, int *randoms)
{
    int i;
    Hashtable *ht = (Hashtable *) malloc (sizeof (Hashtable));
    idxtype *base;

    ht->numNodes = graph->nvtxs;
    ht->numHashes = numHashes;
    ht->sortedNodeIds = idxmalloc (graph->nvtxs,
            "buildHashtable:ht->sortedNodeIds");
    ht->hashes = idxmalloc (graph->nvtxs * numHashes,
            "buildHashtable:ht->hashes");

    int *mins = (int *) malloc (sizeof (int) * numHashes);

    base = ht->hashes;
    idxtype *adjlist = graph->adjncy;
    for (int j = 0; j < graph->nvtxs; j++)
    {
        int adjlistLength = graph->xadj[j + 1] - graph->xadj[j];
        if (adjlistLength == 0)
        {
            for (int i = 0; i < numHashes; i++)
                base[i] = -j;
        }
        else if (adjlistLength == 1)
        {
            for (int i = 0; i < numHashes; i++)
                base[i] = adjlist[0];
        }
        else
        {
            getMinHashKeys (adjlist, adjlistLength, randoms,
                    numHashes, base, mins, numHashes);
        }

        base += numHashes;
        adjlist += adjlistLength;
    }

    free (mins);
    /*	for ( i=0; i<graph->nvtxs; i++ )
        {
        ht->sortedNodeIds[i] = i;
        }
        */
    return ht;
}

/* Sorts the neighbor list of each node by edge weight, and
 * retains the top (log d + 1) neighbors for each node, where d
 * is the degree. */
/* Here we need to sort based on pair's value
 * as it is the weighted sum of similarity */
/* Yiye Ruan */
void retainSimWeightedTopNeighbors (GraphType * graph, IntDoublePair * pairArray, float exponent)
{
    int i, j, k;

    timer tmr;

    cleartimer (tmr);
    starttimer (tmr);

    sortAdjListsByAdjwgt (graph, pairArray);

    stoptimer (tmr);
    printf ("Time for sorting according to weighted similarities:%.2f\n", gettimer (tmr));
    printf ("Done sorting.\n");
    fflush (stdout);

    for (i = 0, j = 0; i < graph->nvtxs; i++)
    {
        k = graph->xadj[i];
        graph->xadj[i] = j;
        int d = graph->xadj[i + 1] - k;
        int numToRetain = 0;
        if (d > 0)
        {
            if (exponent == 0.5)
                numToRetain = (int) ceil (sqrt (d * 1.0));
            else
                numToRetain = (int) ceil (pow (d * 1.0, exponent));
            if (numToRetain <= 0)
                numToRetain = 1;
        }
        if (numToRetain > d)
            numToRetain = d;
        int skip = d - numToRetain;
        for (int l = 0; l < numToRetain; l++, j++)
        {
            // top neighbors are at the end of the list, since
            // parallelqsort sorts in ascending order.
            graph->adjncy[j] = graph->adjncy[k + skip + l];
            /*graph->adjwgt[j] = graph->adjwgt[k+skip+l]; */
        }
    }
    graph->xadj[graph->nvtxs] = j;
    printf ("Done retaining top neighbors.\n");
    fflush (stdout);
}

void sortAdjListsByAdjwgt (GraphType * graph, IntDoublePair * pairArray)
{
    for (idxtype vtxId = 0; vtxId < graph->nvtxs; vtxId++)
    {
        qsort (pairArray + graph->xadj[vtxId],
                graph->xadj[vtxId + 1] - graph->xadj[vtxId], sizeof (*pairArray),
                intDoublePairCmp);
        for (idxtype idx = graph->xadj[vtxId]; idx < graph->xadj[vtxId + 1]; idx++)
        {
            graph->adjncy[idx] = pairArray[idx].id;
        }

#ifdef DEBUG
        if (vtxId < 3)
        {
            printf ("vtxId %d\n", vtxId);
            for (idxtype idx = graph->xadj[vtxId]; idx < graph->xadj[vtxId + 1];
                    idx++)
            {
                printf ("id %d value %f\n", pairArray[idx].id,
                        pairArray[idx].value);
            }
            printf ("\n");
        }
#endif
    }
}

/*Dump the merged network */
void dumpNetwork (const char *fileName, GraphType * graph)
{
    FILE *dump = fopen (fileName, "w");
    int edgeCount = graph->nedges;
    for (int vtxId = 0; vtxId < graph->nvtxs; vtxId++)
    {
        if (graph->adjncy[graph->xadj[vtxId]] == vtxId)
            edgeCount--;
    }
    edgeCount /= 2;
    fprintf (dump, "%d %d\n", graph->nvtxs, edgeCount);
    for (int vtxId = 0; vtxId < graph->nvtxs; vtxId++)
    {
        for (int index = graph->xadj[vtxId]; index < graph->xadj[vtxId + 1];
                index++)
        {
            if (graph->adjncy[index] != vtxId)
                fprintf (dump, "%d ", graph->adjncy[index] + 1);
        }
        fprintf (dump, "\n");
    }
    fclose (dump);
}

int intDoublePairCmp (const void *pair1, const void *pair2)
{
    double diff =
        ((IntDoublePair *) pair1)->value - ((IntDoublePair *) pair2)->value;
    if (diff > 0)
        return 1;
    else if (diff < 0)
        return -1;
    else
        return 0;
}

