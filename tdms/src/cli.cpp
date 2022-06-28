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
#include "io.h"
#include "queuelinked.h"
#include "matio.h"

/*---------------------------------------------------------------*/
//                        DEFINE section
/*---------------------------------------------------------------*/

#define NMATRICES    49 //number of input matrices
#define NOUTMATRICES_WRITE 23  //number of output matrices to be written to output file
#define NOUTMATRICES_WRITE_ALL 25  //number of output matrices to be written to output file
#define NOUTMATRICES_PASSED 31  //number of output matrices passed by mexFunction
#define STRBUFLEN 256


int main(int nargs,char *argv[], char **envp){
  omp_set_num_threads(12);
  FILE *infile;

  /*
    There are two cases to consider, when the fdtdgrid matrix is specified in a separate mat file 
    and when it is in the same file as the other matrices.
  */

  const char *matrixnames[] = {"fdtdgrid","Cmaterial","Dmaterial","C","D","freespace","disp_params","delta","interface","Isource","Jsource","Ksource","grid_labels","omega_an","to_l","hwhm","Dxl","Dxu","Dyl","Dyu","Dzl","Dzu","Nt","dt","tind","sourcemode","runmode","exphasorsvolume","exphasorssurface","intphasorssurface","phasorsurface","phasorinc","dimension","conductive_aux","dispersive_aux","structure","f_ex_vec","exdetintegral","f_vec","Pupil","D_tilde","k_det_obs_global","air_interface","intmatprops","intmethod","tdfield","tdfdir","fieldsample","campssample"};
  
  const char *matrixnames_infile[] = {"Cmaterial","Dmaterial","C","D","freespace","disp_params","delta","interface","Isource","Jsource","Ksource","grid_labels","omega_an","to_l","hwhm","Dxl","Dxu","Dyl","Dyu","Dzl","Dzu","Nt","dt","tind","sourcemode","runmode","exphasorsvolume","exphasorssurface","intphasorssurface","phasorsurface","phasorinc","dimension","conductive_aux","dispersive_aux","structure","f_ex_vec","exdetintegral","f_vec","Pupil","D_tilde","k_det_obs_global","air_interface","intmatprops","intmethod","tdfield","tdfdir","fieldsample","campssample"};
  
  const char *matrixnames_gridfile[] = {"fdtdgrid"};

  const char *outputmatrices_all[] = {"Ex_out","Ey_out","Ez_out","Hx_out","Hy_out","Hz_out","x_out","y_out","z_out","Ex_i","Ey_i","Ez_i","Hx_i","Hy_i","Hz_i","x_i","y_i","z_i","vertices","camplitudes","facets","maxresfield","Id","fieldsample","campssample"};
  const char *outputmatrices[] = {"Ex_out","Ey_out","Ez_out","Hx_out","Hy_out","Hz_out","x_out","y_out","z_out","Ex_i","Ey_i","Ez_i","Hx_i","Hy_i","Hz_i","x_i","y_i","z_i","camplitudes","maxresfield","Id","fieldsample","campssample"};
  int  matricestosave_all[]   = {0,1,2,3,4,5,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28};
  int  matricestosave[]   = {0,1,2,3,4,5,10,11,12,13,14,15,16,17,18,19,20,21,23,25,26,27,28};
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
