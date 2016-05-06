/*
 *
 * $Id: sparsify.c,v 1.7 2010-11-30 23:17:20 venu Exp $
 *
 */

#include <metis.h>


/*************************************************************************
 * Let the game begin
 **************************************************************************/

void print_help (const char *program_name)
{
    /*fprintf(stderr, "Usage: %s <GraphFile> ",program_name); */
    fprintf (stderr,
            "Usage: %s <TopoGraphFile> <ContentGraphFile> <ContentListFile> ",
            program_name);
    fprintf (stderr, "[-f firstGraphWeight] [-n normMode] [-s simMetric] ");
    /*fprintf (stderr, "[-b sparsifyType] [-m numHashes] ");*/
    fprintf (stderr, "[-m numHashes] ");
    /*fprintf (stderr, "[-e sparsifyExponent] [-r globalFraction] [-o output file]");*/
    fprintf (stderr, "[-e sparsifyExponent] [-o output file]");
    fprintf (stderr, "\n");

    /*fprintf (stderr, "\nsparsifyType: 0 (default) - fast,");*/
    /*fprintf (stderr, "approximate, local sparsify\n ");*/
    /*fprintf (stderr, "\t\t1 - slow, exact local sparsify\n ");*/
    /*fprintf (stderr, "\t\t2 - fast, approximate global sparsify\n");*/
    /*fprintf (stderr, "\t\t3 - slow, exact global sim sparsify\n");*/
    fprintf (stderr, "numHashes: number of hashes for approx.");
    fprintf (stderr,
            " similarity (default: 30) (only for approx. similarity)\n");
    fprintf (stderr,
            "sparsifyExponent \"e\" (0 < e < 1, default 0.5): each node with degree d, d^e edges will be retained (only for local sparsification)\n");
    /*fprintf (stderr,*/
            /*"globalFraction (0 < r < 1): fraction of edges to retain in sparsified graph (default: 0.25) (only for global sparsification)\n");*/
    fprintf (stderr, \
            "firstGraphWeight \"alpha\" (0 < f < 1, default 0.5): the weight of first graph's edge similarity when combining scores\n");
    fprintf (stderr, "normMode: Normalization mode\n");
    fprintf (stderr, "\t\t 0 (default) - zero-mean-unit-variance\n");
    fprintf (stderr, "\t\t 1 - [0, 1] interval\n");
    fprintf (stderr, "simMetric: Similarity metric\n");
    fprintf (stderr, "\t\t 0 (default) - cosine\n");
    fprintf (stderr, "\t\t 1 - jaccard\n");
}

int main (int argc, char *argv[])
{

    GraphType topoGraph, contentGraph, contentList;
    char topoFileName[256], contentFileName[256], clistFileName[256],
         outputFile[256];
    int inputFileCount = 0;
    double firstGraphWeight = 0.5;
    int normMode = 0;
    int outputFileGiven = 0;
    int sparsifyType = 0, numHashes = 256, simMetric = 0;
    float globalSamplingRatio = 0.25, degreeExponent = 0.5;
    timer TOTALTmr, METISTmr, IOTmr;
    int topoWgtflag = 0, contentWgtflag = 0, addSelfLoop = 1, txtFormat = 0,
        minEdgeWeight = 0;

    if (argc < 2)
    {
        print_help (argv[0]);
        exit (0);
    }

    for (argv++; *argv != NULL; argv++)
    {
        if ((*argv)[0] == '-')
        {
            int temp;
            switch ((*argv)[1])
            {
                case 'o':
                case 'O':
                    outputFileGiven = 1;
                    strcpy (outputFile, *(++argv));
                    break;
                case 'b':
                case 'B':
                    sparsifyType = atoi (*(++argv));
                    break;
                case 'm':
                case 'M':
                    numHashes = atoi (*(++argv));
                    break;
                case 'e':
                case 'E':
                    degreeExponent = atof (*(++argv));
                    break;
                case 'w':
                case 'W':
                    minEdgeWeight = atoi (*(++argv));
                    break;
                case 'r':
                case 'R':
                    globalSamplingRatio = atof (*(++argv));
                    break;
                case 'f':
                case 'F':
                    firstGraphWeight = atof (*(++argv));
                    break;
                case 'n':
                case 'N':
                    normMode = atoi (*(++argv)); 
                    break;
                case 's':
                case 'S':
                    simMetric = atoi (*(++argv));
                    break;
            }
        }
        else
        {
            switch (inputFileCount)
            {
                case 0:
                    strcpy (topoFileName, *argv);
                    inputFileCount++;
                    break;
                case 1:
                    strcpy (contentFileName, *argv);
                    inputFileCount++;
                    break;
                case 2:
                    strcpy (clistFileName, *argv);
                    inputFileCount++;
                    break;
                default:
                    fprintf (stderr,
                            "Not supposed to take more than %d input files\n",
                            inputFileCount);
                    break;
            }
        }
    }

    if (inputFileCount != 3)
    {
        print_help (argv[0]);
        exit (0);
    }
    else
    {
        printf("Configurations -------------------------------------------------------\n");
        printf("\tTopo graph file: %s\n", topoFileName);
        printf("\tContent graph file: %s\n", contentFileName);
        printf("\tContent list file: %s\n", clistFileName);
        printf("\tSimilarity metric %d\n", simMetric);
        printf("\tNumber of hashes %d\n", numHashes);
        printf("\tNormalization mode: %d\n", normMode);
        printf("\tFirst graph weight: %f\n", firstGraphWeight);
    }

    if (!outputFileGiven)
    {
        fprintf (stderr, "Please specify output file using -o option\n");
        exit (0);
    }

    cleartimer (TOTALTmr);
    cleartimer (METISTmr);
    cleartimer (IOTmr);

    starttimer (TOTALTmr);
    starttimer (IOTmr);

    ReadGraph (&topoGraph, topoFileName, &topoWgtflag, addSelfLoop, txtFormat);

    ReadGraph (&contentGraph, contentFileName, &contentWgtflag, addSelfLoop, txtFormat);

    readList (&contentList, clistFileName, topoGraph.nvtxs);

    if (topoGraph.nvtxs <= 0 || contentGraph.nvtxs <= 0 || contentList.nvtxs <= 0)
    {
        fprintf (stderr, "Empty graph. Nothing to do.\n");
        exit (1);
    }
    stoptimer (IOTmr);

    printf ("Graph Information ----------------------------------------------------\n");
    printf ("  Name: %s, #Vertices: %d, #Edges: %d \n", topoFileName, topoGraph.nvtxs, topoGraph.nedges / 2);
    printf ("Graph Information ----------------------------------------------------\n");
    printf ("  Name: %s, #Vertices: %d, #Edges: %d \n", contentFileName, contentGraph.nvtxs, contentGraph.nedges / 2);
    printf ("Content List Information ---------------------------------------------\n");
    printf ("  Name: %s, #Vertices: %d, #Terms: %d \n", clistFileName, contentList.nvtxs, contentList.nedges);
    fflush (stdout);

    int myIndex, myIter;
#ifdef DEBUG
    printf ("First three nodes for contentList\n");
    for (myIter = 0; myIter < 3; myIter++)
    {
        for (myIndex = contentList.xadj[myIter];
                myIndex < contentList.xadj[myIter + 1]; myIndex++)
            printf ("%d ", contentList.adjncy[myIndex]);
        printf ("\n");
    }
    printf ("First three nodes for topoGraph\n");
    for (myIter = 0; myIter < 3; myIter++)
    {
        for (myIndex = topoGraph.xadj[myIter];
                myIndex < topoGraph.xadj[myIter + 1]; myIndex++)
            printf ("%d ", topoGraph.adjncy[myIndex]);
        printf ("\n");
    }
    printf ("First three nodes for contentGraph\n");
    for (myIter = 0; myIter < 3; myIter++)
    {
        for (myIndex = contentGraph.xadj[myIter];
                myIndex < contentGraph.xadj[myIter + 1]; myIndex++)
            printf ("%d ", contentGraph.adjncy[myIndex]);
        printf ("\n");
    }
#endif

    // Merge the heterogenous graphs into a single representation
    GraphType *inputg;
    inputg = mergeHeteroGraphs (&topoGraph, &contentGraph, NULL, addSelfLoop);
    // Free content graph as its no longer useful
    free (contentGraph.adjncy);
    contentGraph.adjncy = NULL;
    free (contentGraph.xadj);
    contentGraph.xadj = NULL;

    /*DumpNetwork(outputFile, inputg); */
    /*exit(0); */

#ifdef DEBUG
    printf
        ("Graph Information ---------------------------------------------------\n");
    printf ("  #Vertices: %d, #Edges: %d \n", inputg->nvtxs,
            inputg->nedges / 2);
    fflush (stdout);
#ifdef DEBUG
    printf ("First three nodes for inputg\n");
    for (myIter = 0; myIter < 3; myIter++)
    {
        for (myIndex = inputg->xadj[myIter]; myIndex < inputg->xadj[myIter + 1];
                myIndex++)
            printf ("%d ", inputg->adjncy[myIndex]);
        printf ("\n");
    }
    printf ("\n");
#endif
#endif

    starttimer (METISTmr);

    idxtype *new_xadj, *new_adjncy, *new_adjwgt;

    // For now let's just worry with the most basic one - Yiye
    GraphType *rg;

    // DO NOT USE
    // -----------------------------------
    /*Method 1 (sparification of union) */
    /*rg = getSingleBagSparsifiedGraph(inputg, &topoGraph, &contentList, numHashes, degreeExponent); */
    // -----------------------------------

    /*Method 2 (normalize and mix)*/
    if (simMetric == 0)
        rg = getSimWeightedSparsifiedGraph_cos (inputg, &topoGraph, &contentList, numHashes, degreeExponent, firstGraphWeight, normMode);
    else 
        rg = getSimWeightedSparsifiedGraph_jac (inputg, &topoGraph, &contentList, numHashes, degreeExponent, firstGraphWeight, normMode);
    /*Uncomment the following line if you want to do hop1 */
    /*rg = getSimWeightedSparsifiedGraph(&topoGraph, &topoGraph, &contentList, numHashes, degreeExponent, firstGraphWeight, normMode); */

    /*Venu's originial */
    /*rg = getHashSparsifiedGraph(inputTopoG, numHashes, degreeExponent); */
    new_xadj = rg->xadj;
    new_adjncy = rg->adjncy;
    usual_symmetrize_directed (topoGraph.nvtxs, rg->xadj[topoGraph.nvtxs],
            rg->xadj, rg->adjncy, NULL, &new_xadj,
            &new_adjncy, &new_adjwgt, 0);

    /*if ( sparsifyType == 0 ) // Estimate similarity using minhash */
    /*//and then retain similar edges. */
    /*{ */
    /*GraphType *rg; */
    /*rg = getHashSparsifiedGraph(inputg, numHashes, */
    /*degreeExponent); */
    /*//    new_xadj = rg->xadj; */
    /*//    new_adjncy = rg->adjncy; */

    /*new_xadj = rg->xadj; */
    /*new_adjncy = rg->adjncy; */
    /*usual_symmetrize_directed(topoGraph.nvtxs,rg->xadj[topoGraph.nvtxs], */
    /*rg->xadj, rg->adjncy, NULL, &new_xadj, &new_adjncy, */
    /*&new_adjwgt, 0); */
    /*} */
    /*else if (sparsifyType == 2 ) */
    /*{ */
    /*GraphType *rg; */
    /*rg = */
    /*getGlobalHashSparsifiedGraph(inputg,globalSamplingRatio, */
    /*numHashes); */
    /*//    new_xadj = rg->xadj; */
    /*//    new_adjncy = rg->adjncy; */

    /*usual_symmetrize_directed(inputg->nvtxs,rg->xadj[inputg->nvtxs], */
    /*rg->xadj, rg->adjncy, NULL, &new_xadj, &new_adjncy, */
    /*&new_adjwgt, 0); */

    /*} */
    /*else if ( sparsifyType == 1 ) */
    /*{ */
    /*GraphType *rg; */
    /*rg = getExactSimSparsifiedGraph(inputg, degreeExponent); */

    /*usual_symmetrize_directed(inputg->nvtxs,rg->xadj[inputg->nvtxs], */
    /*rg->xadj, rg->adjncy, NULL, &new_xadj, &new_adjncy, */
    /*&new_adjwgt, 0); */
    /*} */
    /*else if ( sparsifyType == 3 ) */
    /*{ */
    /*GraphType *rg; */
    /*rg = getExactGlobalSimSparsifiedGraph(inputg, */
    /*globalSamplingRatio); */

    /*usual_symmetrize_directed(inputg->nvtxs,rg->xadj[inputg->nvtxs], */
    /*rg->xadj, rg->adjncy, NULL, &new_xadj, &new_adjncy, */
    /*&new_adjwgt, 0); */
    /*} */

    stoptimer (METISTmr);

    printf ("No. of edges in sparsified graph:%d\n", (new_xadj[rg->nvtxs] / 2));

    starttimer (IOTmr);
    if (0 /*sparsifyType == 7 */ )
        WriteGraphWithWts (outputFile, rg->nvtxs, new_xadj,
                new_adjncy, new_adjwgt);
    else
        WriteGraph (outputFile, rg->nvtxs, new_xadj, new_adjncy);
    stoptimer (IOTmr);
    stoptimer (TOTALTmr);

    printf
        ("\nTiming Information ---------------------------------------------------\n");
    printf (" I/O:\t\t\t\t\t\t\t%7.3f\n", gettimer (IOTmr));
    printf (" Time for sparsifying:\t\t\t\t\t%7.3f\n", gettimer (METISTmr));
    printf (" Total:\t\t\t\t\t\t\t%7.3f\n", gettimer (TOTALTmr));
    printf
        ("**********************************************************************\n");




}
