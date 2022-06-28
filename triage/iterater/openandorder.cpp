/*****************************************************************
 *
 *  Project.....:  isotropic FDTD code
 *  Application.:  launch and file IO
 *  Module......:  openandorder.cpp
 *  Description.:  Code for processing command line arguments, 
 *                 opening input files,  passing matrices to 
 *                 the mexFunction and writing the output to the 
 *                 specified output file.
 *  Compiler....:  g++
 *  Written by..:  Peter Munro, Imperial College London, 2002-2008
 *  Environment.:  Linux
 *  Modified....:  Numerous times
 *
 ******************************************************************/

/*---------------------------------------------------------------*/
//                        INCLUDE section
/*---------------------------------------------------------------*/
#include <omp.h>
#include <stdio.h>
#include <string.h> 
#include <stdlib.h> 
#include "queuelinked.h"
#include "matio.h"

/*---------------------------------------------------------------*/
//                        DEFINE section
/*---------------------------------------------------------------*/

#define NMATRICES    36 //number of input matrices
#define NOUTMATRICES_WRITE 20  //number of output matrices to be written to output file
#define NOUTMATRICES_WRITE_ALL 22  //number of output matrices to be written to output file
#define NOUTMATRICES_PASSED 26  //number of output matrices passed by mexFunction
#define STRBUFLEN 256

int openandorder(char *matfilename, char **matrixnames, const mxArray **matrixptrs, int nmatrices);
int saveoutput(mxArray **plhs, int *matricestosave, char *matrixnames[], int nmatrics, char *outputfilename);
void freememory( int nmatrices, mxArray **matrixptrs);

int main(int nargs,char *argv[], char **envp){
  omp_set_num_threads(12);
  FILE *infile;

  /*
    There are two cases to consider, when the fdtdgrid matrix is specified in a separate mat file 
    and when it is in the same file as the other matrices.
  */

  const char *matrixnames[] = {"fdtdgrid","Cmaterial","Dmaterial","C","D","freespace","disp_params","delta","interface","Isource","Jsource","Ksource","grid_labels","omega_an","to_l","hwhm","Dxl","Dxu","Dyl","Dyu","Dzl","Dzu","Nt","dt","tind","sourcemode","runmode","exphasorsvolume","exphasorssurface","phasorsurface","phasorinc","dimension","conductive_aux","dispersive_aux","structure","f_ex_vec"};
  
  const char *matrixnames_infile[] = {"Cmaterial","Dmaterial","C","D","freespace","disp_params","delta","interface","Isource","Jsource","Ksource","grid_labels","omega_an","to_l","hwhm","Dxl","Dxu","Dyl","Dyu","Dzl","Dzu","Nt","dt","tind","sourcemode","runmode","exphasorsvolume","exphasorssurface","phasorsurface","phasorinc","dimension","conductive_aux","dispersive_aux","structure","f_ex_vec"};
  
  const char *matrixnames_gridfile[] = {"fdtdgrid"};

  const char *outputmatrices_all[] = {"Ex_out","Ey_out","Ez_out","Hx_out","Hy_out","Hz_out","x_out","y_out","z_out","Ex_i","Ey_i","Ez_i","Hx_i","Hy_i","Hz_i","x_i","y_i","z_i","vertices","camplitudes","facets","maxresfield"};
  const char *outputmatrices[] = {"Ex_out","Ey_out","Ez_out","Hx_out","Hy_out","Hz_out","x_out","y_out","z_out","Ex_i","Ey_i","Ez_i","Hx_i","Hy_i","Hz_i","x_i","y_i","z_i","camplitudes","maxresfield"};
  int  matricestosave_all[]   = {0,1,2,3,4,5,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
  int  matricestosave[]   = {0,1,2,3,4,5,10,11,12,13,14,15,16,17,18,19,20,21,23,25};
  int nmatrices = NMATRICES;
  mxArray *plhs[NOUTMATRICES_PASSED];
  const mxArray *matrixptrs[NMATRICES];
  /*********************************/
  int argi = 1;
  int save_geometry=1;
  int n_args_read = 0;
  char infile_str[STRBUFLEN], outfile_str[STRBUFLEN], gridfile_str[STRBUFLEN];
  /*********************************/
  Queue q = initQueue();
  ql_init(q);
  
  while( argi<nargs ){
    if( argv[argi][0]=='-'){
      //fprintf(stdout,"switch: %s\n",argv[i]);
      switch(argv[argi][1]){
      case 'h':
	fprintf(stdout,"Usage: openandorder [options] infile outfile\n       openandorder [options] infile gridfile outfile\nOptions:\n-h:\tDisplay this help message\n-m:\tMinimise output file size by not saving vertex and facet information\n\n");
	ql_deleteQueue(q);
	exit(0);
      case 'm':
	save_geometry=0;
      }
      argi++;
    }
    else{
      if(strlen(argv[argi])>(STRBUFLEN-1)){
	fprintf(stderr,"argument %d file name is too long (max length %d chars)\n",n_args_read+1,STRBUFLEN-1);
	ql_deleteQueue(q);
	exit(-1);
      }
      ql_push(q,argv[argi]);
      argi++;
      n_args_read++;
    }
  }
  
  if( (n_args_read!=2) &&  (n_args_read!=3) ){
    fprintf(stderr,"incorrect number of arguments\n");
    ql_deleteQueue(q);
    exit(-1);
  }
  else{
    
    ql_pop( q, infile_str, STRBUFLEN);
    if( n_args_read== 3)
      ql_pop( q, gridfile_str, STRBUFLEN);
    ql_pop( q, outfile_str, STRBUFLEN);
    if( n_args_read== 2)
      fprintf(stdout,"infile:[%s], outfile:[%s], m=%d\n",infile_str,outfile_str,save_geometry);
    else
      fprintf(stdout,"infile:[%s], gridfile:[%s], outfile:[%s], m=%d\n",infile_str,gridfile_str,outfile_str,save_geometry);
  }
  ql_deleteQueue(q);
  /*********************************/
  
  //for(temp=0;temp<100;temp++){
  infile = fopen(infile_str,"r");
  if(infile==NULL){
    fprintf(stderr,"Unable to open file %s\n",infile_str);
    exit(-1);
  }
  fclose(infile);
  
  infile = fopen(outfile_str,"a+");
  if(infile==NULL){
    fprintf(stderr,"Unable to open file %s\n",outfile_str);
    exit(0) ;
  }
  fclose(infile);

  if(n_args_read==3){
    infile = fopen(gridfile_str,"r");
    if(infile==NULL){
      fprintf(stderr,"Unable to open file %s\n",gridfile_str);
      exit(0) ;
    }
    fclose(infile);
  }

    
    //now it is safe to use matlab routines to open the file and order the matrices 
  if(n_args_read==2)
    openandorder(infile_str, (char **)matrixnames, matrixptrs, nmatrices);
  else{
    openandorder(infile_str, (char **)matrixnames_infile, matrixptrs+1, nmatrices-1);
    openandorder(gridfile_str, (char **)matrixnames_gridfile, matrixptrs, 1);
  }

    //now run the FDTD code
     mexFunction(NOUTMATRICES_PASSED, (mxArray **)plhs, NMATRICES, (const mxArray **)matrixptrs);
     /*
     MATFile *toutfile;
     toutfile  = matOpen("fdtdgrid_out.mat", "w");
     matPutVariable(toutfile, "fdtdgrid", (mxArray *) matrixptrs[0]);
     matClose(toutfile);
     */

    if( save_geometry ){//prints vertices and facets
      saveoutput(plhs, matricestosave_all, (char **)outputmatrices_all, NOUTMATRICES_WRITE_ALL,outfile_str );
    }
    else{//does not print vertices and facets
      saveoutput(plhs, matricestosave, (char **)outputmatrices, NOUTMATRICES_WRITE,outfile_str );
    freememory(NMATRICES,(mxArray **)matrixptrs);
    freememory(NOUTMATRICES_PASSED,(mxArray **)plhs);
    return 0;
    
  }
  //}
  
}
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
  dims = (mwSize *)malloc(3*sizeof(mwSize));
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
	//	matrixptrs[i] = matGetArray(tmatfile, mat_matrixnames[j]);
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
  
  //first open file
  outfile = matOpen(outputfilename, "w");
  if(outfile==NULL){
    fprintf(stderr,"Couldn't open %s for writing\n",outputfilename);
    return -1;
  }

  //now iterate through the matrices, set names and add to matfile
  for(i=0;i<nmatrices;i++){
    if(matPutVariable(outfile, matrixnames[i], (mxArray *)plhs[matricestosave[i]])){
      fprintf(stderr,"Could not write array %s to %s\n",matrixnames[i],outputfilename);
      return -2;
    }
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
