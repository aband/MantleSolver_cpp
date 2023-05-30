/* 
    Ouput with cgns.
    DMDA global vector.
*/

#include "cgns_io.h"
PetscErrorCode CGNSMeshWrite(DM dm, Vec * fullmesh){
    PetscFunctionBeginUser;
    PetscErrorCode    ierr;

    PetscInt         xs, ys, xm, ym, M, N;
    PetscInt         stencilWidth;
    int              size, rank;
    Vec              fmesh, lmesh;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    cgp_mpi_comm(PETSC_COMM_WORLD);

    fmesh = *fullmesh;

    ierr = DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL);                                          
           CHKERRQ(ierr);
    ierr = DMDAGetInfo(dm, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
           CHKERRQ(ierr);

    const cgsize_t zoneSize[3][2] = 
  // vertices   cells   boundary vertex
    {{M+1,N+1}, {M,N}, {0,0}};

    int index_file, index_base;
    int index_zone, index_grid, index_coordx, index_coordy;

    char basename[33], zonename[33];

    // Open CGNS file for write
    if (cgp_open("meshOut.cgns", CG_MODE_WRITE, &index_file)) cg_error_exit();

    // Create base
    strcpy(basename, "Base");
    int icelldim = 2;
    int iphysdim = 2;

    cg_base_write(index_file, basename, icelldim, iphysdim, &index_base);

    strcpy(zonename, "Mesh");
    // create zone
    if (cg_zone_write(index_file, index_base, zonename, (cgsize_t*)zoneSize,
                      CGNS_ENUMV(Structured), &index_zone)) cg_error_exit();

   if (cg_grid_write(index_file, index_base, index_zone, "GridCoordinates",
                     &index_grid)) cg_error_exit();

    // construct the grid coordinates nodes (user must use SIDS-standard names here)
   if (cgp_coord_write(index_file, index_base, index_zone,
                       CGNS_ENUMV(RealSingle), "CoordinateX",
                       &index_coordx)) cgp_error_exit();
   if (cgp_coord_write(index_file, index_base, index_zone,
                       CGNS_ENUMV(RealSingle), "CoordinateY",
                       &index_coordy)) cgp_error_exit();

    // Collective writing of file data
    int numLocalx, numLocaly;
    if (xs+xm == M){
        numLocalx = xm+1;
    }else {
        numLocalx = xm;
    }

    if (ys+ym == N){
        numLocaly = ym+1;
    }else {
        numLocaly = ym;
    }

    cgsize_t   s_rmin[2], s_rmax[2], m_dimvals[2], m_rmin[2], m_rmax[2];

    double *x = NULL;
    double *y = NULL;

    // Create gridpoints
    const int num_vertex = numLocalx*numLocaly;
    x = (double*)malloc(num_vertex*sizeof(double));
    y = (double*)malloc(num_vertex*sizeof(double));

    // Assign uniform grids to all coordinates
    Vector2D *data = (Vector2D *)malloc(numLocalx*numLocaly*sizeof(Vector2D));
    Vector2D **localcoords = (Vector2D **)malloc(numLocaly * sizeof(Vector2D **));

    for (int i=0; i<numLocaly; i++){
        localcoords[i] = &(data[numLocalx*i]);
    }

    // Get Coordinates
    Vector2D **coords;

    ierr = DMGetLocalVector(dm, &lmesh);
    ierr = DMGlobalToLocalBegin(dm, fmesh, INSERT_VALUES, lmesh);
    ierr = DMGlobalToLocalEnd(dm, fmesh, INSERT_VALUES, lmesh);

    ierr = DMDAVecGetArray(dm, lmesh, &coords);CHKERRQ(ierr);
    for (int j=0; j<numLocaly; j++){
    for (int i=0; i<numLocalx; i++){
        localcoords[j][i].x = coords[j+ys][i+xs].x;
        localcoords[j][i].y = coords[j+ys][i+xs].y;
    }}

    ierr = DMDAVecRestoreArray(dm, lmesh, &coords);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm, &lmesh);

    // assign 1d data thereafter
    for (int j=0; j<numLocaly; j++){
    for (int i=0; i<numLocalx; i++){
        int id = j*numLocalx+i;
        x[id] = localcoords[j][i].x;
        y[id] = localcoords[j][i].y;
    }}

    // Shape in file space
    s_rmin[0] = xs+1;
    s_rmax[0] = xs+numLocalx;
    s_rmin[1] = ys+1;
    s_rmax[1] = ys+numLocaly; 

    // Shape in memory
    m_dimvals[0] = numLocalx;
    m_rmin[0]    = 1;
    m_rmax[0]    = numLocalx;

    m_dimvals[1] = numLocaly;
    m_rmin[1]    = 1;
    m_rmax[1]    = numLocaly;

     if (cgp_coord_general_write_data(index_file, index_base, index_zone, 1, s_rmin, s_rmax, CGNS_ENUMV(RealDouble),2,m_dimvals, m_rmin, m_rmax, x)){
      cgp_error_exit();
    }

    if (cgp_coord_general_write_data(index_file, index_base, index_zone, 2, s_rmin, s_rmax, CGNS_ENUMV(RealDouble),2,m_dimvals, m_rmin, m_rmax, y)){
      cgp_error_exit();
    }

    free(x);
    free(y);

    PetscFunctionReturn(0);
}

/*
PetscErrorCode  DMDACgnsOut2D(DM dmMesh, Vec *fullmesh, DM dmCell, Vec *Sol, char *filename)
{
    // Cgns output code for 2D mesh  
    PetscErrorCode   ierr;
    PetscInt         xs, ys, xm, ym, M, N;
    int              size, rank;
    Vec              fmesh;
    PetscFunctionBeginUser;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    cgp_mpi_comm(PETSC_COMM_WORLD);

    fmesh = *fullmesh;

    ierr = DMDAGetCorners(dmMesh, &xs, &ys, NULL, &xm, &ym, NULL);                                       CHKERRQ(ierr);
    ierr = DMDAGetInfo(dmMesh, NULL, &M, &N, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

    const cgsize_t zoneSize[3][2] = {{M+1, N+1}, {M,N}, {0,0}};

    int index_file, index_base;
    int index_zone, index_grid, index_coordx, index_coordy;
    char basename[33], zonename[33];

    // Open cgns file for write
    if (cgp_open(filename, CG_MODE_WRITE, &index_file)){
        cg_error_exit();
    }

    // Create base
    strcpy(basename, "Base");
    int icelldim = 2;
    int iphysdim = 2;

    cg_base_write(index_file, basename, icelldim, iphysdim, &index_base);

    // Define zone name
    strcpy(zonename, "Grid");

    // create zone
    if (cg_zone_write(index_file, index_base, zonename, (cgsize_t*)zoneSize,CGNS_ENUMV(Structured), &index_zone)){
        cg_error_exit();
    }

    if (cg_grid_write(index_file, index_base, index_zone, "GridCoordinates", &index_grid)){
        cg_error_exit();
    }

    // Construct the grid coordinates nodes
    if (cgp_coord_write(index_file, index_base, index_zone,
                        CGNS_ENUMV(RealDouble), "CoordinateX",
                        &index_coordx)){
        cgp_error_exit();
    }

    if (cgp_coord_write(index_file, index_base, index_zone,
                        CGNS_ENUMV(RealDouble), "CoordinateY",
                        &index_coordy)){
        cgp_error_exit();
    }

    // Collective writing of coordinate data
    int numLocalx, numLocaly;
    if (xs+xm == M){
        numLocalx = xm+1;
    }else {
        numLocalx = xm;
    }

    if (ys+ym == N){
        numLocaly = ym+1;
    }else {
        numLocaly = ym;
    }

    cgsize_t   s_rmin[2], s_rmax[2], m_dimvals[2], m_rmin[2], m_rmax[2];

    double *x = NULL;
    double *y = NULL;

    // Create gridpoints
    const int num_vertex = numLocalx*numLocaly;
    x = (double*)malloc(num_vertex*sizeof(double));
    y = (double*)malloc(num_vertex*sizeof(double));

    // Assign uniform grids to all coordinates
    Vector2D *data = (Vector2D *)malloc(numLocalx*numLocaly*sizeof(Vector2D));
    Vector2D **localcoords = (Vector2D **)malloc(numLocaly * sizeof(Vector2D **));

    for (int i=0; i<numLocaly; i++){
        localcoords[i] = &(data[numLocalx*i]);
    }

    double hx = 1.0/(double)M;
    double hy = 1.0/(double)N;

    for (int j=0; j<numLocaly; j++){
    for (int i=0; i<numLocalx; i++){
        localcoords[j][i].x = (i+xs)*hx;
        localcoords[j][i].y = (j+ys)*hy;
    }}

    // Get Coordinates
    Vector2D **coords;

    ierr = DMDAVecGetArray(dmMesh, fmesh, &coords);CHKERRQ(ierr);
    for (int j=0; j<ym; j++){
    for (int i=0; i<xm; i++){
        localcoords[j][i].x = coords[j+ys][i+xs].x;
        localcoords[j][i].y = coords[j+ys][i+xs].y;
    }}

    ierr = DMDAVecRestoreArray(dmMesh, fmesh, &coords);CHKERRQ(ierr);

    // assign 1d data thereafter
    for (int j=0; j<numLocaly; j++){
    for (int i=0; i<numLocalx; i++){
        int id = j*numLocalx+i;
        x[id] = localcoords[j][i].x;
        y[id] = localcoords[j][i].y;
    }}


    // Shape in file space
    s_rmin[0] = xs+1;
    s_rmax[0] = xs+numLocalx;
    s_rmin[1] = ys+1;
    s_rmax[1] = ys+numLocaly; 

    // Shape in memory
    m_dimvals[0] = numLocalx;
    m_rmin[0]    = 1;
    m_rmax[0]    = numLocalx;

    m_dimvals[1] = numLocaly;
    m_rmin[1]    = 1;
    m_rmax[1]    = numLocaly;

    // Collectively write coordinates to file
    if (cgp_coord_general_write_data(index_file, index_base, index_zone, 1, s_rmin, s_rmax, CGNS_ENUMV(RealDouble),2,m_dimvals, m_rmin, m_rmax, x)){
      cgp_error_exit();
    }

    if (cgp_coord_general_write_data(index_file, index_base, index_zone, 2, s_rmin, s_rmax, CGNS_ENUMV(RealDouble),2,m_dimvals, m_rmin, m_rmax, y)){
      cgp_error_exit();
    }

    free(x);
    free(y);

    // Solution output
    Vec               gSol;
    int               index_flow, index_field;
    PetscScalar       **sol;
    PetscInt          cxs, cys, cxm, cym;
    char              solname[33];

    ierr = DMDAGetCorners(dmCell, &cxs, &cys, NULL, &cxm, &cym, NULL);CHKERRQ(ierr);

    gSol = *Sol;

    strcpy(solname, "FlowSolution");

    // Create flow solution node 
    cg_sol_write(index_file, index_base, index_zone, solname,
                 CGNS_ENUMV(CellCenter), &index_flow);

    // Go to position within tree at FlowSolution_t node 
    cg_goto(index_file, index_base, "Zone_t", index_zone, 
            "FlowSolution_t", index_flow, "end");

    if (cgp_field_write(index_file, index_base, index_zone,
                        index_flow, CGNS_ENUMV(RealDouble),
                        "Solution", &index_field))
        cgp_error_exit();

    ierr = DMDAVecGetArrayRead(dmCell, gSol, &sol);CHKERRQ(ierr);

    const int num_cell = ym*xm;
    double *soll;
    soll = (double*)malloc(num_cell*sizeof(double));

    for (int j=ys; j<ys+ym; j++){
    for (int i=xs; i<xs+xm; i++){
        int id = (j-ys)*xm+i-xs;
        soll[id] = (double)sol[j][i];
    }}

    ierr = DMDAVecRestoreArrayRead(dmCell, gSol, &sol);CHKERRQ(ierr);

    int c_s_rmin[2], c_s_rmax[2], c_m_dimvals[2], c_m_rmin[2], c_m_rmax[2];

    c_s_rmin[0] = cxs+1;
    c_s_rmax[0] = cxs+cxm;

    c_m_dimvals[0] = cxm;
    c_m_rmin[0]    = 1;
    c_m_rmax[0]    = cxm;
    
    c_s_rmin[1] = cys+1;
    c_s_rmax[1] = cys+cym;

    c_m_dimvals[1] = cym;
    c_m_rmin[1]    = 1;
    c_m_rmax[1]    = cym;

    if (cgp_field_general_write_data(index_file, index_base, index_zone,index_flow,1,c_s_rmin, c_s_rmax, CGNS_ENUMV(RealDouble),2, c_m_dimvals, c_m_rmin, c_m_rmax, soll)) cgp_error_exit();

    free(soll);

    cgp_close(index_file);

    PetscFunctionReturn(0);
}*/
