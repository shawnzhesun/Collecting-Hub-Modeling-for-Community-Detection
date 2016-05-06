#include <metis.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

chunktype * generateHashTable(const GraphType * const graph, const int numHash)
{
#ifndef DEBUG
    srand(time(0));
#endif

    // Each hash value is a bit, so we can pack multiple bits in one chunk
    size_t numChunkPerSignature = numHash / (sizeof(chunktype) * 8);
    // How many space do we need
    size_t numChunkPerTable = graph->nvtxs * numChunkPerSignature;
    chunktype *hashtable = (chunktype *) calloc(numChunkPerTable, sizeof(*hashtable));
    assert (hashtable != NULL);

    size_t universeSize = getUniverseSize(graph);
    size_t numProjPerTable = universeSize * numHash;
    /*
     * projtable will be indexed by numHash first
     * i.e. projtable[elementIdx * numHash + hashIdx] will be 
     * the hashIdx-th projection bit for the elementIdx-th element
     */
    projtype *projtable = (projtype *) calloc(numProjPerTable, sizeof(*projtable));
    assert (projtable != NULL);
    populateProjTable(projtable, numProjPerTable);
    populateHashTable(hashtable, numHash, numChunkPerSignature, projtable, graph);
    free(projtable);
    projtable = NULL;
    return hashtable;
}

void populateHashTable(chunktype * const hashtable, const int numHash, const int numChunkPerSignature, const projtype * const projtable, const GraphType * const graph)
{
    idxtype nodeIdx;
    for (nodeIdx = 0; nodeIdx < graph->nvtxs; nodeIdx++)
    {
        /*
         * Two cases for the length of the list:
         * 0: The list is empty, do nothing.
         * >0: Normal operation
         */
        int listLength = graph->xadj[nodeIdx+1] - graph->xadj[nodeIdx];
        if (listLength)
        {
            int *sumArray = (int *) calloc(numHash, sizeof(*sumArray));
            assert (sumArray != NULL);
            size_t elementIdx, hashIdx, elementBase;
            for (elementIdx = graph->xadj[nodeIdx]; elementIdx < graph->xadj[nodeIdx+1]; elementIdx++)
            {
                elementBase = graph->adjncy[elementIdx] * numHash;
                for (hashIdx = 0; hashIdx < numHash; hashIdx++)
                    sumArray[hashIdx] += projtable[elementBase++];
            }
            /*
             * Now we have counted the number of +1's in random projection
             * elements corresponding to the list elements.
             * We then need to compute the number -1's (list length -
             * +1's), and compare the two counts.
             */
            for (hashIdx = 0; hashIdx < numHash; hashIdx++)
                sumArray[hashIdx] = (2 * sumArray[hashIdx] >= listLength);
            packBits(sumArray, numHash, numChunkPerSignature, hashtable + nodeIdx * numChunkPerSignature);
            free(sumArray);
            sumArray = NULL;
        }
    }
}

void packBits(const int * const sumArray, const int numHash, const int numChunk, chunktype * const hashtable)
{
    size_t sumIdx = 0, inIdx, outIdx;
    for (outIdx = 0; outIdx < numChunk; outIdx++)
    {
        chunktype val = 0;
        for (inIdx = 0; inIdx < sizeof(chunktype) * 8; inIdx++)
            val = (val << 1) | sumArray[sumIdx++];
        hashtable[outIdx] = val;
    }
}

void populateProjTable(projtype * const projtable, const size_t tableSize)
{
    size_t projIdx = 0, outIdx, inIdx;
    int val;
    for (outIdx = 0; outIdx < tableSize / (sizeof(int) * 8); outIdx ++)
    {
        val = rand();
        for (inIdx = 0; inIdx < sizeof(int) * 8; inIdx++)
        {
            projtable[projIdx++] = val & 1;
            val >>= 1;
        }
    }
    val = rand();
    for (inIdx = 0; inIdx < tableSize % (sizeof(int) * 8); inIdx++)
    {
        projtable[projIdx++] = val & 1;
        val >>= 1;
    }
#ifdef DEBUG
    assert (projIdx == tableSize);
#endif
}

/*
 * Compute the universe size of the graph
 */
size_t getUniverseSize(const GraphType * const graph)
{
    size_t maxId = 0;
    /*
     * Assume each adjacency list in graph is already sorted
     * so we only need to compare the last element in each list
     */
    size_t nodeIdx;
    for (nodeIdx = 0; nodeIdx < graph->nvtxs; nodeIdx++)
    {
        size_t lastId = graph->adjncy[graph->xadj[nodeIdx+1]-1];
        maxId = (lastId > maxId) ? lastId : maxId;
    }
    // Index start from 0
    return maxId + 1;
}

/*
 * Count the number of 1-bits 
 * Reference: http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetKernighan
 */
size_t countBitSet(chunktype val)
{
    size_t ret;
    for (ret = 0; val; ret++)
        val &= val - 1; // clear the least significant bit set
    return ret;
}

double cosSimByLSH(const size_t id1, const size_t id2, const chunktype * const hashtable, const GraphType * const graph, const int numHash, const int numChunk)
{
    if (graph->xadj[id1 + 1] - graph->xadj[id1] == 0 || graph->xadj[id2 + 1] - graph->xadj[id2] == 0)
        return 0;
    int match = 0;
    size_t idx;
    size_t id1Base = id1 * numChunk, id2Base = id2 * numChunk;
    for (idx = 0; idx < numChunk; idx++)
        /*
         * We are using XOR here, which counts the number of bits that DO NOT MATCH
         * The final calculatin formula is accordingly modified 
         */
        match += countBitSet(hashtable[id1Base++] ^ hashtable[id2Base++]);
#ifdef DEBUG
    printf("from %lu to %lu\n", id1, id2);
    printf("list 1: ");
    for (idx = graph->xadj[id1]; idx < graph->xadj[id1+1]; idx++)
        printf("%d ", graph->adjncy[idx]);
    printf("\n");
    printf("list 2: ");
    for (idx = graph->xadj[id2]; idx < graph->xadj[id2+1]; idx++)
        printf("%d ", graph->adjncy[idx]);
    printf("\n");
    printf("signature 1: ");
    for (idx = 0; idx < numChunk; idx++)
        printf("%d ", hashtable[id1 * numChunk + idx]);
    printf("\n");
    printf("signature 2: ");
    for (idx = 0; idx < numChunk; idx++)
        printf("%d ", hashtable[id2 * numChunk + idx]);
    printf("\n");
    int sum1 = 0, sum2 = 0;
    for (idx = 0; idx < numChunk; idx++)
        sum1 += countBitSet(hashtable[id1 * numChunk + idx]);
    for (idx = 0; idx < numChunk; idx++)
        sum2 += countBitSet(hashtable[id2 * numChunk + idx]);
    printf("sum 1 %d sum 2 %d\n", sum1, sum2);
    printf("match is %d\n", numHash - match);
    printf("Return %f\n", cos(match * M_PI / numHash));
    /*exit(0);*/
#endif
    return cos(match * M_PI / numHash);
}
