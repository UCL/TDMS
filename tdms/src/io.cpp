#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "queuelinked.h"
#include "matio.h"

/*Returns:

-1 - could not open file
-2 - could not find particular matrix
-3 - too few matrices in mat file
-4 - either I_tot, J_tot or K_tot was not specified in the mat file
-5 - NULL pointer returned when getting matric pointer
*/
int openandorder(char *matfilename, char **matrixnames, const mxArray **matrixptrs, int nmatrices){
  char **mat_matrixnames;
  MATFile *tmatfile;
  int i,j, num_mats;
  mwSize ndims;
  //int dims[3];
  mwSize *dims;
  dims = (mwSize *)malloc(4*sizeof(mwSize));
  int I_tot, J_tot, K_tot;
  const char fdtdgrid_elements[][15] = {"Exy","Exz","Eyx","Eyz","Ezx","Ezy","Hxy","Hxz","Hyx","Hyz","Hzx","Hzy"};
  mxArray *element;


  tmatfile = matOpen(matfilename, "r");
  if(tmatfile==NULL){
    fprintf(stderr, "Error opening %s\n",matfilename);
    return(-1);
  }

  //get matrix names
  mat_matrixnames = matGetDir(tmatfile, &num_mats);

  if( num_mats < nmatrices){
    fprintf(stderr,"Not enough matrices in mat file (%d,%d)\n",num_mats,nmatrices);
    return(-3);
  }

  //now iterate through the matrix names and assign the pointer to each matrix
  //into the appropriate entry of matrixptrs
  for(i=0;i<nmatrices;i++){
    //now find the matrix name
    for(j=0;j<num_mats;j++){
      if(!strcmp(mat_matrixnames[j],matrixnames[i])){
        //matrix pointer found
        //fprintf(stderr,"Got %s/%s (%d)\n",mat_matrixnames[j],matrixnames[i],j);
        //	//	matrixptrs[i] = matGetArray(tmatfile, mat_matrixnames[j]);
        matrixptrs[i] = matGetVariable(tmatfile, mat_matrixnames[j]);
        if( matrixptrs[i] == NULL){
          fprintf(stderr, "Could not get pointer to %s\n",mat_matrixnames[j]);
          return -5;
        }
        break;
      }
      else if(j==(num_mats-1)){
        //matrix pointer NOT found
        fprintf(stderr, "Couldn't find matrix %s\n",matrixnames[i]);
        return(-2);
      }
    }
  }
  //fprintf(stderr, "Got all %d matrices\n",nmatrices);
  //fprintf(stderr, "%d\n",mxGetFieldNumber( (mxArray *)matrixptrs[0], "I_tot"));
  if(!strcmp(matrixnames[0],"fdtdgrid")){
    //get the dimensions of the fdtd grid
    if( mxGetFieldNumber( (mxArray *)matrixptrs[0], "I_tot") == -1 ){
      fprintf(stderr, "%s missing field fdtdgrid.I_tot\n",matfilename);
      return -4;
    }
    element = mxGetField( (mxArray *)matrixptrs[0], 0, "I_tot");
    I_tot = (int) (*mxGetPr(element));
    mxRemoveField((mxArray *)matrixptrs[0], mxGetFieldNumber(matrixptrs[0], "I_tot"));
    if( mxGetFieldNumber(matrixptrs[0], "J_tot") == -1 ){
      fprintf(stderr, "%s missing field fdtdgrid.J_tot\n",matfilename);
      return -4;
    }
    element = mxGetField( (mxArray *)matrixptrs[0], 0, "J_tot");
    J_tot = (int) *mxGetPr(element);
    mxRemoveField((mxArray *)matrixptrs[0], mxGetFieldNumber( (mxArray *)matrixptrs[0], "J_tot"));
    if( mxGetFieldNumber(matrixptrs[0], "K_tot") == -1 ){
      fprintf(stderr, "%s missing field fdtdgrid.K_tot\n",matfilename);
      return -4;
    }
    element = mxGetField( (mxArray *)matrixptrs[0], 0, "K_tot");
    K_tot = (int) *mxGetPr(element);
    mxRemoveField((mxArray *)matrixptrs[0], mxGetFieldNumber( (mxArray *)matrixptrs[0], "K_tot"));

    //fprintf(stderr,"%d %d %d\n",I_tot,J_tot,K_tot);

    //now need to add fields to fdtdgrid
    for(i=0;i<12;i++){
      mxAddField((mxArray *)(matrixptrs[0]),fdtdgrid_elements[i] );
      //now need to populate the fields with zero matrices
      element = mxGetField( (mxArray *)matrixptrs[0], 0, fdtdgrid_elements[i]);
      ndims = 3;

      dims[0] = I_tot + 1;
      dims[1] = J_tot + 1;
      dims[2] = K_tot + 1;

      element = mxCreateNumericArray( (const mwSize)ndims, (const mwSize *)dims, mxDOUBLE_CLASS, mxREAL); //this function initialises data to 0
      mxSetField(( mxArray *)matrixptrs[0], 0, fdtdgrid_elements[i],element);
    }
  }

  free(dims);
  matClose(tmatfile);
  mxFree(mat_matrixnames);
  return 0;
}


/*Save the resultant matrices
  Returns:
  -1 if can not open mat file for writing
  -2 error adding matrix to file


*/

int saveoutput(mxArray **plhs, int *matricestosave, char *matrixnames[], int nmatrices, char *outputfilename){

  MATFile *outfile;
  int i;
  int mpv_out;
  FILE *fp;

  //first open file
  outfile = matOpen(outputfilename, "w7.3");
  if(outfile==NULL){
    fprintf(stderr,"Couldn't open %s for writing\n",outputfilename);
    return -1;
  }

  //now iterate through the matrices, set names and add to matfile
  for(i=0;i<nmatrices;i++){
    mpv_out = matPutVariable(outfile, matrixnames[i], (mxArray *)plhs[matricestosave[i]]);

    if(mpv_out){
      fp = matGetFp(outfile);
      fprintf(stderr,"Could not write array %s to %s (%d,%d,%d)\n",matrixnames[i],outputfilename,mpv_out,feof(fp),ferror(fp));

    }
    //return -2;
  }
  //int matPutArray(MATFile *mfp, const mxArray *mp);
  //void mxSetName(mxArray *array_ptr, const char *name);
  matClose(outfile);
  return 0;

}

void freememory( int nmatrices, mxArray **matrixptrs){
  int i;
  for(i=0;i<nmatrices;i++)
    mxDestroyArray(matrixptrs[i]);
}
