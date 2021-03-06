void triangulatePlane(int I0, int I1, int J0, int J1, int K,int coordmap[], int order, mxArray **vertexMatrix);
void triangulatePlaneSkip(int I0, int I1, int J0, int J1, int K,int coordmap[], int order, mxArray **vertexMatrix, int dI, int dJ);
void triangulateCuboid(int I0, int I1, int J0, int J1, int K0, int K1, mxArray **vertexMatrix);
void triangulateCuboidSkip(int I0, int I1, int J0, int J1, int K0, int K1, mxArray **vertexMatrix, int dI, int dJ, int dK);
void conciseTriangulateCuboid(int I0, int I1, int J0, int J1, int K0, int K1, 
			      mxArray **vertices, mxArray ** facets);
void conciseTriangulateCuboidSkip(int I0, int I1, int J0, int J1, int K0, int K1, 
				  int dI, int dJ, int dK,
			     mxArray **vertices, mxArray ** facets);
void pointsToVector(int p1[], int p2[], int *vector);
void crossProduct(int v1[], int v2[], int *v1crossv2);

void conciseCreateBoundary(int I0, int I1,int K0, int K1, mxArray **vertices, mxArray ** facets);
