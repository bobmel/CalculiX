/*  patchAndGoMPI.c  */

#include "../spoolesMPI.h"
#include "../../timings.h"

/*--------------------------------------------------------------------*/
int
main ( int argc, char *argv[] ) {
/*
   ------------------------------------------------------------
   all-in-one MPI program for each process

   order, factor and solve A X = Y
   We use a patch-and-go strategy 
   for the factorization without pivoting

   ( 1) read in matrix entries and form InpMtx object for A
   ( 2) order the system using minimum degree
   ( 3) permute the front tree
   ( 4) create the owners map IV object
   ( 5) permute the matrix A and redistribute
   ( 6) compute the symbolic factorization 
   ( 7) compute the numeric factorization
   ( 8) split the factors into submatrices
   ( 9) create the submatrix map and redistribute
   (10) read in right hand side entries 
        and form dense matrix DenseMtx object for Y
   (11) permute and redistribute Y
   (12) solve the linear system
   (13) gather X on processor 0

   created -- 98jun13, cca
   ------------------------------------------------------------
*/
/*--------------------------------------------------------------------*/
char            buffer[20] ;
Chv             *rootchv ;
ChvManager      *chvmanager ;
DenseMtx        *mtxX, *mtxY ;
double          fudge, toosmall ;
SubMtxManager   *mtxmanager, *solvemanager ;
FrontMtx        *frontmtx ;
InpMtx          *keepmtx, *mtxA ;
double          cutoff, imag, minops, real, tau = 100., value ;
double          cpus[20] ;
double          *opcounts ;
DV              *cumopsDV ;
ETree           *frontETree ;
FILE            *inputFile, *msgFile ;
Graph           *graph ;
int             error, firsttag, ient, irow, jcol, msglvl, myid, nedges,
                nent, neqns, nmycol, nproc, nrhs, nrow, patchAndGoFlag,
                root, seed, storeids, storevalues, symmetryflag, type ;
int             stats[20] ;
int             *rowind ;
IV              *oldToNewIV, *ownedColumnsIV, *ownersIV, 
                *newToOldIV, *vtxmapIV ;
IVL             *adjIVL, *symbfacIVL ;
SolveMap        *solvemap ;
/*--------------------------------------------------------------------*/
/*
   ---------------------------------------------------------------
   find out the identity of this process and the number of process
   ---------------------------------------------------------------
*/
MPI_Init(&argc, &argv) ;
MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
fprintf(stdout, "\n proc %d, argc = %d", myid, argc) ;
fflush(stdout) ;
/*--------------------------------------------------------------------*/
/*
   --------------------
   get input parameters
   --------------------
*/
if ( argc != 11 ) {
   fprintf(stdout, 
      "\n usage: %s msglvl msgFile type symmetryflag patchAndGoFlag"
      "\n        fudge toosmall storeids storevalues seed"
      "\n    msglvl -- message level"
      "\n    msgFile -- message file"
      "\n    type    -- type of entries"
      "\n      1 (SPOOLES_REAL)    -- real entries"
      "\n      2 (SPOOLES_COMPLEX) -- complex entries"
      "\n    symmetryflag -- type of matrix"
      "\n      0 (SPOOLES_SYMMETRIC)    -- symmetric entries"
      "\n      1 (SPOOLES_HERMITIAN)    -- Hermitian entries"
      "\n      2 (SPOOLES_NONSYMMETRIC) -- nonsymmetric entries"
      "\n    patchAndGoFlag -- flag for the patch-and-go strategy"
      "\n      0 -- none, stop factorization"
      "\n      1 -- optimization strategy"
      "\n      2 -- structural analysis strategy"
      "\n    fudge       -- perturbation parameter"
      "\n    toosmall    -- upper bound on a small pivot"
      "\n    storeids    -- flag to store ids of small pivots"
      "\n    storevalues -- flag to store perturbations"
      "\n    seed -- random number seed"
      "\n    "
      "\n   note: matrix entries are read in from patchMatrix.k.input"
      "\n         where k is the process number"
      "\n   note: rhs entries are read in from patchRhs.k.input"
      "\n         where k is the process number"
      "\n", argv[0]) ;
   return(0) ;
}
msglvl = atoi(argv[1]) ;
if ( strcmp(argv[2], "stdout") == 0 ) {
   msgFile = stdout ;
} else {
   sprintf(buffer, "res.%d", myid) ;
   if ( (msgFile = fopen(buffer, "w")) == NULL ) {
      fprintf(stderr, "\n fatal error in %s"
              "\n unable to open file %s\n",
              argv[0], buffer) ;
      return(-1) ;
   }
}
type           = atoi(argv[3]) ;
symmetryflag   = atoi(argv[4]) ;
patchAndGoFlag = atoi(argv[5]) ;
fudge          = atof(argv[6]) ;
toosmall       = atof(argv[7]) ;
storeids       = atoi(argv[8]) ;
storevalues    = atoi(argv[9]) ;
seed           = atoi(argv[10]) ;
IVzero(20, stats) ;
DVzero(20, cpus) ;
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------
   STEP 1: read the entries from the input file 
           and create the InpMtx object
   --------------------------------------------
*/
sprintf(buffer, "patchMatrix.%d.input", myid) ;
inputFile = fopen(buffer, "r") ;
fscanf(inputFile, "%d %d %d", &neqns, &neqns, &nent) ;
mtxA = InpMtx_new() ;
InpMtx_init(mtxA, INPMTX_BY_ROWS, type, nent, 0) ;
if ( type == SPOOLES_REAL ) {
   for ( ient = 0 ; ient < nent ; ient++ ) {
      fscanf(inputFile, "%d %d %le", &irow, &jcol, &value) ;
      InpMtx_inputRealEntry(mtxA, irow, jcol, value) ;
   }
} else if ( type == SPOOLES_COMPLEX ) {
   for ( ient = 0 ; ient < nent ; ient++ ) {
      fscanf(inputFile, "%d %d %le %le", &irow, &jcol, &real, &imag) ;
      InpMtx_inputComplexEntry(mtxA, irow, jcol, real, imag) ;
   }
}
fclose(inputFile) ;
InpMtx_sortAndCompress(mtxA) ;
InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n input matrix") ;
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   STEP 2 : find a low-fill ordering
   (1) create the Graph object
   (2) order the graph using multiple minimum degree
   (3) find out who has the best ordering w.r.t. op count,
       and broadcast that front tree object
   -------------------------------------------------------
*/
graph = Graph_new() ;
adjIVL = InpMtx_MPI_fullAdjacency(mtxA, stats, 
                                  msglvl, msgFile, MPI_COMM_WORLD) ;
nedges = IVL_tsize(adjIVL) ;
Graph_init2(graph, 0, neqns, 0, nedges, neqns, nedges, adjIVL,
            NULL, NULL) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n graph of the input matrix") ;
   Graph_writeForHumanEye(graph, msgFile) ;
   fflush(msgFile) ;
}
frontETree = orderViaMMD(graph, seed + myid, msglvl, msgFile) ;
Graph_free(graph) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n front tree from ordering") ;
   ETree_writeForHumanEye(frontETree, msgFile) ;
   fflush(msgFile) ;
}
opcounts = DVinit(nproc, 0.0) ;
opcounts[myid] = ETree_nFactorOps(frontETree, type, symmetryflag) ;
MPI_Allgather((void *) &opcounts[myid], 1, MPI_DOUBLE,
              (void *) opcounts, 1, MPI_DOUBLE, MPI_COMM_WORLD) ;
minops = DVmin(nproc, opcounts, &root) ;
DVfree(opcounts) ;
frontETree = ETree_MPI_Bcast(frontETree, root, 
                             msglvl, msgFile, MPI_COMM_WORLD) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n best front tree") ;
   ETree_writeForHumanEye(frontETree, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------------
   STEP 3: get the permutations and permute the front tree
   -------------------------------------------------------
*/
oldToNewIV = ETree_oldToNewVtxPerm(frontETree) ;
newToOldIV = ETree_newToOldVtxPerm(frontETree) ;
ETree_permuteVertices(frontETree, oldToNewIV) ;
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   STEP 4: generate the owners map IV object
           and the map from vertices to owners
   -------------------------------------------
*/
cutoff   = 1./(2*nproc) ;
cumopsDV = DV_new() ;
DV_init(cumopsDV, nproc, NULL) ;
ownersIV = ETree_ddMap(frontETree, 
                       type, symmetryflag, cumopsDV, cutoff) ;
DV_free(cumopsDV) ;
vtxmapIV = IV_new() ;
IV_init(vtxmapIV, neqns, NULL) ;
IVgather(neqns, IV_entries(vtxmapIV), 
         IV_entries(ownersIV), ETree_vtxToFront(frontETree)) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n map from fronts to owning processes") ;
   IV_writeForHumanEye(ownersIV, msgFile) ;
   fprintf(msgFile, "\n\n map from vertices to owning processes") ;
   IV_writeForHumanEye(vtxmapIV, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------
   STEP 5: permute the matrix and redistribute
   -------------------------------------------
*/
InpMtx_permute(mtxA, IV_entries(oldToNewIV), IV_entries(oldToNewIV)) ;
if (  symmetryflag == SPOOLES_SYMMETRIC 
   || symmetryflag == SPOOLES_HERMITIAN ) { 
   InpMtx_mapToUpperTriangle(mtxA) ;
}
InpMtx_changeCoordType(mtxA, INPMTX_BY_CHEVRONS) ;
firsttag = 0 ;
keepmtx = InpMtx_MPI_split(mtxA, vtxmapIV, stats, 
                           msglvl, msgFile, firsttag, MPI_COMM_WORLD) ;
InpMtx_free(mtxA) ;
mtxA = keepmtx ;
InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n permuted and split InpMtx") ;
   InpMtx_writeForHumanEye(mtxA, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   STEP 6: compute the symbolic factorization
   ------------------------------------------
*/
symbfacIVL = SymbFac_MPI_initFromInpMtx(frontETree, ownersIV, mtxA,
                     stats, msglvl, msgFile, firsttag, MPI_COMM_WORLD) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n local symbolic factorization") ;
   IVL_writeForHumanEye(symbfacIVL, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------
   STEP 7: initialize the front matrix object and the
      PatchAndGoInfo object to handle small pivots
   --------------------------------------------------
*/
mtxmanager = SubMtxManager_new() ;
SubMtxManager_init(mtxmanager, NO_LOCK, 0) ;
frontmtx = FrontMtx_new() ;
FrontMtx_init(frontmtx, frontETree, symbfacIVL, type, symmetryflag,
              FRONTMTX_DENSE_FRONTS, SPOOLES_NO_PIVOTING, NO_LOCK, myid,
              ownersIV, mtxmanager, msglvl, msgFile) ;
if ( patchAndGoFlag == 1 ) {
   frontmtx->patchinfo = PatchAndGoInfo_new() ;
   PatchAndGoInfo_init(frontmtx->patchinfo, 1, toosmall, fudge,
                       storeids, storevalues) ;
} else if ( patchAndGoFlag == 2 ) {
   frontmtx->patchinfo = PatchAndGoInfo_new() ;
   PatchAndGoInfo_init(frontmtx->patchinfo, 2, toosmall, fudge,
                       storeids, storevalues) ;
}
/*--------------------------------------------------------------------*/
/*
   ---------------------------------
   STEP 8: compute the factorization
   ---------------------------------
*/
chvmanager = ChvManager_new() ;
ChvManager_init(chvmanager, NO_LOCK, 0) ;
rootchv = FrontMtx_MPI_factorInpMtx(frontmtx, mtxA, tau, 0.0 ,
                     chvmanager, ownersIV, 0, &error, cpus, 
                     stats, msglvl, msgFile, firsttag, MPI_COMM_WORLD) ;
ChvManager_free(chvmanager) ;
if ( patchAndGoFlag == 1 ) {
   if ( frontmtx->patchinfo->fudgeIV != NULL ) {
      fprintf(msgFile, "\n small pivots found at these locations") ;
      IV_writeForHumanEye(frontmtx->patchinfo->fudgeIV, msgFile) ;
   }
   PatchAndGoInfo_free(frontmtx->patchinfo) ;
} else if ( patchAndGoFlag == 2 ) {
   if ( frontmtx->patchinfo->fudgeIV != NULL ) {
      fprintf(msgFile, "\n small pivots found at these locations") ;
      IV_writeForHumanEye(frontmtx->patchinfo->fudgeIV, msgFile) ;
   }
   if ( frontmtx->patchinfo->fudgeDV != NULL ) {
      fprintf(msgFile, "\n perturbations") ;
      DV_writeForHumanEye(frontmtx->patchinfo->fudgeDV, msgFile) ;
   }
   PatchAndGoInfo_free(frontmtx->patchinfo) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n numeric factorization") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
if ( error >= 0 ) {
   fprintf(stderr, 
          "\n proc %d : factorization error at front %d", myid, error) ;
   fprintf(msgFile, 
          "\n proc %d : factorization error at front %d", myid, error) ;
   MPI_Barrier(MPI_COMM_WORLD) ;
   fprintf(stderr, 
          "\n proc %d : calling MPI_Finalize()", myid) ;
   fprintf(msgFile, 
          "\n proc %d : calling MPI_Finalize()", myid) ;
   MPI_Finalize() ;
   return(-1) ;
}
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   STEP 9: post-process the factorization and split 
   the factor matrices into submatrices 
   ------------------------------------------------
*/
FrontMtx_MPI_postProcess(frontmtx, ownersIV, stats, msglvl,
                         msgFile, firsttag, MPI_COMM_WORLD) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n numeric factorization after post-processing");
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   -----------------------------------
   STEP 10: create the solve map object
   -----------------------------------
*/
solvemap = SolveMap_new() ;
SolveMap_ddMap(solvemap, frontmtx->symmetryflag, 
               FrontMtx_upperBlockIVL(frontmtx),
               FrontMtx_lowerBlockIVL(frontmtx),
               nproc, ownersIV, FrontMtx_frontTree(frontmtx), 
               seed, msglvl, msgFile);
if ( msglvl > 3 ) {
   SolveMap_writeForHumanEye(solvemap, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   ----------------------------------------------------
   STEP 11: redistribute the submatrices of the factors
   ----------------------------------------------------
*/
FrontMtx_MPI_split(frontmtx, solvemap, 
                   stats, msglvl, msgFile, firsttag, MPI_COMM_WORLD) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n numeric factorization after split") ;
   FrontMtx_writeForHumanEye(frontmtx, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   -------------------------------------------------
   STEP 12: read the entries from the rhs input file 
   and create the DenseMtx object for Y
   -------------------------------------------------
*/
sprintf(buffer, "patchRhs.%d.input", myid) ;
inputFile = fopen(buffer, "r") ;
fscanf(inputFile, "%d %d", &nrow, &nrhs) ;
mtxY = DenseMtx_new() ;
DenseMtx_init(mtxY, type, myid, -1, nrow, nrhs, 1, nrow) ;
DenseMtx_rowIndices(mtxY, &nrow, &rowind) ;
if ( type == SPOOLES_REAL ) {
   for ( irow = 0 ; irow < nrow ; irow++ ) {
      fscanf(inputFile, "%d", rowind + irow) ;
      for ( jcol = 0 ; jcol < nrhs ; jcol++ ) {
         fscanf(inputFile, "%le", &value) ;
         DenseMtx_setRealEntry(mtxY, irow, jcol, value) ;
      }
   }
} if ( type == SPOOLES_COMPLEX ) {
   for ( irow = 0 ; irow < nrow ; irow++ ) {
      fscanf(inputFile, "%d", rowind + irow) ;
      for ( jcol = 0 ; jcol < nrhs ; jcol++ ) {
         fscanf(inputFile, "%le %le", &real, &imag) ;
         DenseMtx_setComplexEntry(mtxY, irow, jcol, real, imag) ;
      }
   }
}
fclose(inputFile) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n rhs matrix in original ordering") ;
   DenseMtx_writeForHumanEye(mtxY, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------------
   STEP 13: permute and redistribute Y if necessary
   ------------------------------------------------
*/
DenseMtx_permuteRows(mtxY, oldToNewIV) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n rhs matrix in new ordering") ;
   DenseMtx_writeForHumanEye(mtxY, msgFile) ;
   fflush(msgFile) ;
}
if ( FRONTMTX_IS_PIVOTING(frontmtx) ) {
   IV   *rowmapIV ;
/*
   ----------------------------------------------------------
   pivoting has taken place, redistribute the right hand side
   to match the final rows and columns in the fronts
   ----------------------------------------------------------
*/
   rowmapIV = FrontMtx_MPI_rowmapIV(frontmtx, ownersIV, msglvl,
                                    msgFile, MPI_COMM_WORLD) ;
   mtxY = DenseMtx_MPI_splitByRows(mtxY, rowmapIV, stats, msglvl, 
                                   msgFile, firsttag, MPI_COMM_WORLD) ;
   IV_free(rowmapIV) ;
} else {
/*
   --------------------------------------------------------------
   pivoting has not taken place, redistribute the right hand side
   using the vertex map used to redistribute the matrix A
   --------------------------------------------------------------
*/
   mtxY = DenseMtx_MPI_splitByRows(mtxY, vtxmapIV, stats, msglvl, 
                                   msgFile, firsttag, MPI_COMM_WORLD) ;
}
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n rhs matrix after split") ;
   DenseMtx_writeForHumanEye(mtxY, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   ------------------------------------------
   STEP 14: create a solution DenseMtx object
            and solve the linear system
   ------------------------------------------
*/
ownedColumnsIV = FrontMtx_ownedColumnsIV(frontmtx, myid, ownersIV,
                                         msglvl, msgFile) ;
nmycol = IV_size(ownedColumnsIV) ;
mtxX = DenseMtx_new() ;
if ( nmycol > 0 ) {
   DenseMtx_init(mtxX, type, myid, -1, nmycol, nrhs, 1, nmycol) ;
   DenseMtx_rowIndices(mtxX, &nrow, &rowind) ;
   IVcopy(nmycol, rowind, IV_entries(ownedColumnsIV)) ;
}
solvemanager = SubMtxManager_new() ;
SubMtxManager_init(solvemanager, NO_LOCK, 0) ;
FrontMtx_MPI_solve(frontmtx, mtxX, mtxY, solvemanager, solvemap, cpus, 
                   stats, msglvl, msgFile, firsttag, MPI_COMM_WORLD) ;
SubMtxManager_free(solvemanager) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n solution in new ordering") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
}
/*--------------------------------------------------------------------*/
/*
   --------------------------------------------------------
   STEP 15: permute the solution into the original ordering
            and assemble the solution onto processor zero
   --------------------------------------------------------
*/
DenseMtx_permuteRows(mtxX, newToOldIV) ;
if ( msglvl > 2 ) {
   fprintf(msgFile, "\n\n solution in old ordering") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
   fflush(msgFile) ;
}
IV_fill(vtxmapIV, 0) ;
firsttag++ ;
mtxX = DenseMtx_MPI_splitByRows(mtxX, vtxmapIV, stats, msglvl, msgFile,
                                firsttag, MPI_COMM_WORLD) ;
if ( myid == 0 && msglvl > 0 ) {
   fprintf(msgFile, "\n\n complete solution in old ordering") ;
   DenseMtx_writeForHumanEye(mtxX, msgFile) ;
   fflush(msgFile) ;
}
/*--------------------------------------------------------------------*/
MPI_Finalize() ;

return(1) ; }
/*--------------------------------------------------------------------*/
