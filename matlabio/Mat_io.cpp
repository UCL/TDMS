//---------------------------------------------------------------------------
#pragma hdrstop
#include<string>
#include<fstream>
#include<iostream>
using namespace std;
#include "Mat_io.h"
#include<time.h>
//#include<time>
//---------------------------------------------------------------------------

/*
"MATLAB 5.0 MAT-file, Platform: PCWIN, Created on: Tue Jun 22 22:08:51 2004"
fill spaces
00 01
49 4d (IM)

signed 8 bit 		 1	miINT8
unsigned 8 bit 		 2	miUINT8
signed 16-bit		 3	miINT16
unsigned 16-bit		 4	miUINT16
signed 32-bit		 5	miINT32
unsigned 32-bit		 6	miUINT32
single IEEE 754		 7	miSINGLE
Reserved			 8	--
double IEEE 754		 9	miDOUBLE
Reserved			10	--
Reserved			11	--
signed 64-bit		12	miINT64
unsigned 64-bit		13	miUINT64
MATLAB array		14	miMATRIX
Compressed Data		15	miCOMPRESSED
Unicode UTF-8 Data	16	miUTF8
Unicode UTF-16 Data	17	miUTF16
Unicode UTF-32 Data	18	miUTF32

Unknown						 0	mxUNKNOWN_CLASS
Cell array					 1	mxCELL_CLASS
Structure 					 2	mxSTRUCT_CLASS
Logical array				 3	mxLOGICAL_CLASS
Character array				 4	mxCHAR_CLASS
Void array					 5	mxVOID_CLASS
Double precision array		 6	mxDOUBLE_CLASS
Single precision array		 7	mxSINGLE_CLASS
8-bit, signed integer		 8	mxINT8_CLASS
8-bit, unsigned integer		 9	mxUINT8_CLASS
16-bit, signed integer		10	mxINT16_CLASS
16-bit, unsigned integer	11	mxUINT16_CLASS
32-bit, signed integer		12	mxINT32_CLASS
32-bit unsigned, integer	13	mxUINT32_CLASS
64-bit, signed integer		14	mxINT64_CLASS
64-bit unsigned, integer	15	mxUINT64_CLASS
Function 					16	mxFUNCTION_CLASS
Opaque class				17	mxOPAQUE_CLASS
Object	 					18	mxOBJECT_CLASS

*/

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Mat file format level 5, with or withour compression, not all functions are implemented and not checked for level 4 compatibility
MATFile* matOpen(const char* filename, const char* mode){
	MATFile* mf1;
	bool check1=true;
	switch(mode[0]){
		case 'r' : mf1 = new MATFile(filename,0); break;
		case 'u' : mf1 = new MATFile(filename,1); break;
		case 'w' : mf1 = new MATFile(filename,2); break;//w4 and w7.3 not implemented
	}
	if(!mf1->IsOpen())mf1=NULL;
	else
	{
		switch(mf1->GetMode()){
			case 0 : check1 = mf1->CheckHeader();
			         if(!check1){
						delete mf1;
						mf1 = NULL;
			         }
					 break;
			case 1 : 
			case 2 : mf1->WriteHeader(); break;
		}
	}
	return mf1;
}

//---------------------------------------------------------------------------

int matClose(MATFile* mf1){
	int ans=0;
	if(mf1!=NULL)delete mf1;
	else ans = 1;
	return ans;
}

//---------------------------------------------------------------------------

mxArray* matGetVariable(MATFile* mf1, const char* ArrayName){
  int AL = 0;
	while(ArrayName[AL]!=(char)0) AL++;
	miMatrix DataType, DataType_c;
	mxClassID ArrayType;
	int Bytes, Bytes_c;
	bool small, small_c;
	mxArray* mxA;
	mxA = new mxArray;
	bool found = false;
	bool Strm_ok = true;
	char* AN0, *AN1;
	int nc;
	int Npos;
	mf1->SetRWmode(true);
	mf1->SetPos(128);
	
	int count=1;
	while((!found)&&(!mf1->EoF())){
	  mf1->ReadTag(DataType, Bytes, small);
		Npos = mf1->GetPos() + Bytes;
		if(Bytes>0){
			if(DataType==miCOMPRESSED)Strm_ok = mf1->InitStream(Npos);
			if(Strm_ok)
			{
				if(DataType==miCOMPRESSED)
					mf1->ReadTag(DataType_c, Bytes_c, small_c);
				mf1->ReadMatrixHeader(mxA);
				nc=mxA->GetNc();
				if(nc>0){
					found=true;
					AN0 = mxA->GetName();
				}
				if(AL!=nc)found=false;
				else for(int i=0; i<nc; i++)if(ArrayName[i]!=AN0[i])found=false;
				if(!found){
					if(DataType==miCOMPRESSED)mf1->EndStream();
					mf1->SetPos(Npos);
				}
			}
		}
	}
	if(found){
	  mf1->ReadMatrixDataElement(mxA);
	  if(DataType==miCOMPRESSED)mf1->EndStream();
	  mf1->SetNindex(Npos);
	}
	else mxA = NULL;
	return mxA;
}


//---------------------------------------------------------------------------

mxArray* matGetNextVariable(MATFile* mf1, const char ** NamePtr){
	miMatrix DataType, DataType_c;
	mxClassID ArrayType;
	int Bytes, Bytes_c;
	int nc;
	bool small, small_c;
	bool Strm_ok = false;
	mxArray* mxA;
	mxA = new mxArray;
	int Npos;
	mf1->SetRWmode(true);
	mf1->SetPos(mf1->GetNindex());
	mf1->ReadTag(DataType, Bytes, small);
	if(!mf1->EoF()){
		Npos = mf1->GetPos() + Bytes;
		if(Bytes>0){
			if(DataType==miCOMPRESSED)Strm_ok = mf1->InitStream(Npos);
			if(Strm_ok)
			{
				if(DataType==miCOMPRESSED)
					mf1->ReadTag(DataType_c, Bytes_c, small_c);
				mf1->ReadMatrixHeader(mxA);
				mf1->ReadMatrixDataElement(mxA);
			}
			if(DataType==miCOMPRESSED)mf1->EndStream();
			nc=mxA->GetNc();
			if(nc>0)NamePtr[0] = mxA->GetName();
		}
		else mxA = NULL;
		mf1->SetNindex(Npos);
	}
	else{
		mxA = NULL;
		mf1->Reset();
	}
	return mxA;
}

//---------------------------------------------------------------------------

int matPutVariable(MATFile* mf1, char* Name, mxArray* mxA){
	mxArray* mxA1 = mxA;
	int nc = 0, Npos_c;
	while(Name[nc]!=(char)0) nc++;
	mxA1->SetName(Name,nc);
	mf1->SetPos(-1);
	bool small=false, small_c=false;
	bool Strm_ok = true;
	mf1->SetRWmode(false);
	miMatrix DataType_c = miMATRIX;
#if _ALLOW_COMPRESSION==1
	miMatrix DataType = miCOMPRESSED;
#else
	miMatrix DataType = miMATRIX;
#endif
	mf1->WriteTag(DataType,0,small);
	int Npos=mf1->GetPos();
//	cout << "Npos = " << Npos << endl;
	
	if((mxA1!=NULL)&&(!mxA1->GetEmpty())){
		if(DataType==miCOMPRESSED)Strm_ok = mf1->InitStream(Npos);
//		cout << "Strm_ok = " << Strm_ok << endl;
		if(Strm_ok)
		{
			if(DataType==miCOMPRESSED){
				mf1->WriteTag(DataType_c,0,small_c);
				Npos_c=mf1->GetPos();
//				cout << "Npos_c(1) = " << Npos_c << endl;
			}
			mf1->WriteMatrixHeader(mxA1);
			mf1->WriteMatrixDataElement(mxA1);
			if(DataType==miCOMPRESSED){
				small_c=true;
//				cout << "Npos_c(2) = " << Npos_c << endl;
				mf1->WriteTag(DataType_c,Npos_c,small_c);
			}
		}
		if(DataType==miCOMPRESSED)mf1->EndStream();
	}
	small=true;
	mf1->WriteTag(DataType,Npos,small);
	return 0;
}
//---------------------------------------------------------------------------
/*
char **matGetDir(MATFile *mf1, int *num){
  char **retval;
  char tempbuf[200][200];
  *num = 0;
  mxArray *mxA;
  mxA = matGetNextVariable(mf1, (const char **)&tempbuf[(num[0])++]);
  while( mxA != NULL ){
    fprintf(stderr,"[%s]\n",tempbuf[num[0]-1]);
    mxA = matGetNextVariable(mf1, (const char **)&tempbuf[(num[0])++]);
  }
  (*num)--;
  retval = (char **)malloc(sizeof(char *)*num[0]);
  for(int i=0;i<num[0];i++){
    retval[i] = (char *)malloc(sizeof(char)*(strlen(tempbuf[i])+1));
  }
    
}
*/


char **matGetDir(MATFile *mf1, int *num){
  char **retval;
  char **tempbuf;
  char **temp2;
  temp2 = (char **)malloc(sizeof(char *));
  temp2[0] = (char *)malloc(100*sizeof(char));
  tempbuf = (char **)malloc(sizeof(char *)*200);
  *num = 0;
  mxArray *mxA;
  //  mxA = matGetNextVariable(mf1, (const char **)&tempbuf[(num[0])++]);
  mxA = matGetNextVariable(mf1, (const char **)temp2);
  num[0]++;
  while( mxA != NULL ){
    tempbuf[num[0]-1] = (char *)malloc(sizeof(char)*100);
    strcpy(tempbuf[num[0]-1],temp2[0]);
    mxA = matGetNextVariable(mf1, (const char **)temp2);
    num[0]++;
    

  }
  
  (*num)--;
  retval = (char **)malloc(sizeof(char *)*num[0]);
  for(int i=0;i<num[0];i++){
    retval[i] = (char *)malloc(sizeof(char)*(strlen(tempbuf[i])+1));
    strcpy(retval[i],tempbuf[i]);
  }
  
  return retval;
    
}


/*
char **matGetDir(MATFile *mf1, int *num){

  char **retval;
  char tempbuf[200][200];
  miMatrix DataType, DataType_c;
  mxClassID ArrayType;
  int Bytes, Bytes_c;
  bool small, small_c;
  mxArray* mxA;
  mxA = new mxArray;
  bool Strm_ok = true;
  char* AN0, *AN1;
  int nc;
  int Npos;
  mf1->SetRWmode(true);
  mf1->SetPos(128);
  *num = 0;
  while(!mf1->EoF()){
    mf1->ReadTag(DataType, Bytes, small);
    Npos = mf1->GetPos() + Bytes;
    //if( (Bytes>0) && (DataType<=miUTF32) ){
      if( (Bytes>0) ){
      
      if(DataType==miCOMPRESSED)Strm_ok = mf1->InitStream(Npos);
      if(Strm_ok)
	{
	  if(DataType==miCOMPRESSED)
	    mf1->ReadTag(DataType_c, Bytes_c, small_c);
	  mf1->ReadMatrixHeader(mxA);
	  nc=mxA->GetNc();
	  	  if(nc>0){
	    AN0 = mxA->GetName();
	    strcpy(tempbuf[*num],AN0);
	    (*num)++;
	  }
	}
      if(DataType==miCOMPRESSED)
	mf1->EndStream();
      mf1->SetNindex(Npos);
    }
    
  }
  (*num)--;
  retval = (char **)malloc( (*num)*sizeof(char *) );
  for(int i=0;i<(*num);i++){
    retval[i] = (char *)malloc( sizeof(char *)*(strlen(tempbuf[i])+1) );
    strcpy(retval[i],tempbuf[i]);
  }
  mf1->Reset();
  return retval;
}
*/

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void mxDestroyArray(mxArray *mxA){
	delete mxA;
}

//---------------------------------------------------------------------------

void* mxGetData(mxArray *mxA){
	return mxA->GetDataP(0);
}

//---------------------------------------------------------------------------

double* mxGetPr(mxArray *mxA){
	return (double *)mxA->GetDataP(0);
}

//---------------------------------------------------------------------------

double* mxGetPi(mxArray *mxA){
	return (double *)mxA->GetDataP(1);
}

//---------------------------------------------------------------------------

mxLogical* mxGetLogicals(mxArray *mxA){
	return (mxLogical *)mxA->GetDataP(0);
}

//---------------------------------------------------------------------------

void mxFree(void *ptr){
  free(ptr);
}


//---------------------------------------------------------------------------
/*
char* mxGetName(mxArray *mxA){
	char *Name0;
	if(mxA->GetNc()>0)Name0 = mxA->GetName();
	else Name0 = NULL;
	return Name0;
}
*/
//---------------------------------------------------------------------------
/*
void mxSetName(mxArray *mxA, char* Name){
	int nc = 0;
	while(Name[nc]!=(char)0) nc++;
	mxA->SetName(Name,nc);
}
*/
//---------------------------------------------------------------------------
/*
void mxClearLogical(mxArray *mxA){
	mxA->SetLogical(false);
}
*/
//---------------------------------------------------------------------------
/*
void mxSetLogical(mxArray *mxA){
	mxA->SetLogical(true);
}
*/
//---------------------------------------------------------------------------

const mwSize* mxGetDimensions(mxArray *mxA){
	int nd = mxA->GetNd();
	mwSize *dims;
	dims = new mwSize[nd];
	for(int j=0; j<nd; j++){
		dims[j] = (mwSize)mxA->GetDim(j);
	}
	return (const mwSize *)dims;
}

//---------------------------------------------------------------------------

mwSize *mxGetNumberOfElements(const mxArray *pm){
  mwSize *numel;
  numel = (mwSize *)malloc(sizeof(mwSize));
  *numel = 1;
  int nd = ((mxArray *)pm)->GetNd();
  const mwSize *dims;
  dims = mxGetDimensions((mxArray *)pm);
  for(int i=0;i<nd;i++)
    *numel = (*numel)*dims[i];
  return numel;
}

//---------------------------------------------------------------------------

bool mxIsLogical(mxArray *mxA){
	return mxA->GetLogical();
}

//---------------------------------------------------------------------------

bool mxIsComplex(mxArray *mxA){
	return mxA->GetComplex();
}

//---------------------------------------------------------------------------

mxArray* mxGetField(mxArray *mxA, mwIndex n1, const char* Name1){
	int N1 = (int)n1;
	mxArray* mxA0;
	int NL1 = 0, NL2;
	while(Name1[NL1]!=(char)0)NL1++;
	char* Name2;
	int nf = mxA->GetNf();
	int n0=-1;
	for(int i=0; i<nf; i++){
		Name2 = mxA->GetFieldName(i);
		NL2 = 0;
		while(Name2[NL2]!=(char)0)NL2++;
		if(NL1==NL2){
			n0=i;
			for(int j=0; j<NL1; j++)
				if(Name1[j]!=Name2[j])n0=-1;
		}
		else n0=-1;
		if(n0>-1)i=nf;
	}
	if(n0>-1){
		n0=n0+N1*nf;
		mxA0 = mxA->GetMatrixElm(n0);
	}
	else mxA0 = NULL;
	return mxA0;
}

//---------------------------------------------------------------------------

void mxSetField(mxArray* mxA, mwIndex n1, const char* Name1, mxArray* mxA0){
	int N1 = int(n1);
	int NL1 = 0, NL2;
	while(Name1[NL1]!=(char)0)NL1++;
	char* Name2;
	int nf = mxA->GetNf();
	int n0=-1;
	for(int i=0; i<nf; i++){
		Name2 = mxA->GetFieldName(i);
		NL2 = 0;
		while(Name2[NL2]!=(char)0)NL2++;
		if(NL1==NL2){
			n0=i;
			for(int j=0; j<NL1; j++)
				if(Name1[j]!=Name2[j])n0=-1;
		}
		else n0=-1;
		if(n0>-1)i=nf;
	}
	if(n0>-1){
		n0=n0+N1*nf;
		mxA->SetMatrixElm(mxA0,n0);
	}
}

//---------------------------------------------------------------------------

extern int mxAddField(mxArray *pm, const char *fieldname){
  return (pm)->AddField(fieldname, 32);
  
  //mxArray* mxCreateStructMatrix(mwSize ndim, const mwSize *dims, int nf, const char** FieldNames)

}

//---------------------------------------------------------------------------
//mxArray *mxCreateStructArray(mwSize ndim, const mwSize *dims,int nfields, const char **fieldnames);
mxArray* mxCreateStructMatrix(mwSize ndim, const mwSize *dims, int nf, const char** FieldNames){
  //int dim[2];
  //dim[0] = (int)n;
  //dim[1] = (int)m;
	mxArray* mxA;
	mxA = new mxArray;
	mxA->SetEmpty();
	mxA->SetDim((int *)dims,ndim);
	mxA->SetArrayType(mxSTRUCT_CLASS);
	mxA->SetFieldNames(nf, 32, FieldNames);
	int Elements;
	mxA->MakeDataP(miINT8,0,Elements,0);//DataType, Bytes and mode irrelevant for struct
	return mxA;
}

//---------------------------------------------------------------------------

mxArray* mxCreateNumericMatrix(mwSize n, mwSize m, mxClassID ArrayType, mxComplexity RealComplex){
	int dim[2];
	dim[0] = (int)n;
	dim[1] = (int)m;
	mxArray* mxA;
	mxA = new mxArray;
	mxA->SetEmpty();
	mxA->SetDim(dim,2);
	mxA->SetArrayType(ArrayType);
	int Elements;
	mxA->MakeDataP(miINT8,dim[0]*dim[1],Elements,0);//DataType must match with Bytes like size and elements
	if(RealComplex==mxCOMPLEX){
		mxA->MakeDataP(miINT8,dim[0]*dim[1],Elements,1);
		mxA->SetComplex(true);
	}
	return mxA;
}

//---------------------------------------------------------------------------

mxArray *mxCreateDoubleMatrix(mwSize m, mwSize n, mxComplexity ComplexFlag){
  return mxCreateNumericMatrix( m, n, mxDOUBLE_CLASS, ComplexFlag);
}

//---------------------------------------------------------------------------

mxArray* mxCreateLogicalMatrix(mwSize n, mwSize m){
	int dim[2];
	dim[0] = (int)n;
	dim[1] = (int)m;
	mxArray* mxA;
	mxA = new mxArray;
	mxA->SetEmpty();
	mxA->SetDim(dim,2);
	mxA->SetArrayType(mxUINT8_CLASS);
	int Elements;
	mxA->MakeDataP(miUINT8,dim[0]*dim[1],Elements,0);//DataType must match with Bytes like size and elements
	mxA->SetLogical(true);
	return mxA;
}

//---------------------------------------------------------------------------

mxArray* mxCreateNumericArray(mwSize n, const mwSize* dim, mxClassID ArrayType, mxComplexity RealComplex){
	int N = (int)n;
	int* dim0;
	int Elm=1;
	dim0 = new int[N];
	for(int i=0; i<N; i++){
		dim0[i]=(int)dim[i];
		Elm=Elm*((int)dim[i]);
	}
	int Elements;
	mxArray* mxA;
	mxA = new mxArray;
	mxA->SetEmpty();
	mxA->SetDim(dim0,N);
	mxA->SetArrayType(ArrayType);
	mxA->MakeDataP(miINT8,Elm,Elements,0);//DataType must match with Bytes like size and elements
	if(RealComplex==mxCOMPLEX){
		mxA->MakeDataP(miINT8,Elm,Elements,1);
		mxA->SetComplex(true);
	}
	return mxA;
}

//---------------------------------------------------------------------------

int mxGetString(mxArray* mxA, char* Str, mwSize maxSize){
	int MAXSize = (int)maxSize;
	int fail=0;
	if(mxA->GetArrayType()!=mxCHAR_CLASS)fail=1;
	char* Data0 = (char *)mxA->GetDataP(0);
	int min=mxA->GetNdata();
	if(MAXSize<min){
		min=MAXSize;
		fail=1;
	}
	for(int i=0; i<min; i++)
		Str[i] = Data0[i];
	Str[min]=(char)0;
	return fail;
}

//---------------------------------------------------------------------------

mxArray* mxCreateString(const char* Str){
	int SL=0;
	while(Str[SL]!=(char)0)SL++;
	int dim[2];
	dim[0] = 1;
	dim[1] = SL;
	mxArray* mxA;
	mxA = new mxArray;
	mxA->SetEmpty();
	mxA->SetDim(dim,2);
	mxA->SetArrayType(mxCHAR_CLASS);
	int Elements;
	char *Data0;
	Data0 = (char *)mxA->MakeDataP(miINT8,SL,Elements,0);//DataType must match with Bytes like size and elements
	for(int i=0; i<SL; i++)
		Data0[i] = Str[i];
	return mxA;
}

//---------------------------------------------------------------------------

void* mxCalloc(size_t n,size_t size){
	void* Data0;
	Data0 = calloc(n,size);
	return Data0;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

int mexPrintf(const char* Buffer){
	cout << Buffer;
	return 0;
}

void mexErrMsgTxt(const char *errormsg){
  cerr << errormsg;
}

//---------------------------------------------------------------------------

int mexEvalString(const char*){
	return 0;
}//"drawnow;"


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

MATFile::MATFile(const char* filename, int mode1){
	mode = mode1;
	z_stream strm;
	strm_in = new unsigned char[CHUNK];
	strm_out = new unsigned char[CHUNK];
	ReadStrm=0;
	strm_error=Z_OK;
	Nindex = 128;
	Npos = 128;
	switch(mode){
		case 0 : File1.open(filename,ios_base::in|ios_base::binary);break;
		case 1 : File1.open(filename,ios_base::in|ios_base::out|ios_base::binary);break;
		case 2 : File1.open(filename,ios_base::in|ios_base::out|ios_base::trunc|ios_base::binary);break;
	}
	store_pos = new int * [MAX_STORED];
	for(int j=0; j<MAX_STORED; j++)
		store_pos[j] = new int [2];
}

//---------------------------------------------------------------------------

MATFile::~MATFile(){
	File1.close();
	if(ReadStrm>0)EndStream();
	delete[] strm_in;
	delete[] strm_out;
	for(int j=0; j<MAX_STORED; j++)
		delete[] store_pos[j];
	delete[] store_pos;
}

//---------------------------------------------------------------------------

bool MATFile::IsOpen(){
	return File1.is_open();
}

//---------------------------------------------------------------------------

bool MATFile::EoF(){
	return File1.eof();
}

//---------------------------------------------------------------------------

bool MATFile::Fail(){
	bool ret0;
	if(ReadStrm>0)ret0 = !((strm_error==Z_OK)||(strm_error==Z_STREAM_END));
	else ret0 = File1.fail();
	return ret0;
}

//---------------------------------------------------------------------------

int MATFile::GetMode(){
	return mode;
}


//---------------------------------------------------------------------------

void MATFile::Reset(){
	Nindex = 128;
	if(ReadStrm>0){
		cout << "should not reset in stream, error?" << endl;
		EndStream();
	}
	File1.clear();
	if(RWmode)File1.seekg(Nindex,ios_base::beg);
	else File1.seekp(Nindex,ios_base::beg);
}

//---------------------------------------------------------------------------

bool MATFile::InitStream(int end_pos){
	strm_begin=GetPos();
//	cout << "strm_begin = " << strm_begin << endl;
	Npos = end_pos;
	strm_idx=0;
	strm_pos=0;
	write_pos=0;
	if(RWmode)ReadStrm=1;
	else ReadStrm=2;
	if(ReadStrm==1){
		strm.zalloc=Z_NULL;
		strm.zfree=Z_NULL;
		strm.opaque=Z_NULL;
		strm.next_in=Z_NULL;
		strm.avail_in=0;
		strm_flush=Z_NO_FLUSH;
		strm_error=inflateInit(&strm);
		strm.next_in=strm_in;
		strm.avail_in=0;
		strm.next_out=strm_out;
		strm.avail_out=CHUNK;
	}
	else{
		strm.zalloc=Z_NULL;
		strm.zfree=Z_NULL;
		strm.opaque=Z_NULL;
		strm_flush=Z_NO_FLUSH;
		strm_error=deflateInit(&strm,Z_DEFAULT_COMPRESSION);
		strm.next_in=strm_in;
		strm.avail_in=0;
		strm.next_out=strm_out;
		strm.avail_out=CHUNK;
	}
	return ((strm_error == Z_OK)||(strm_error == Z_STREAM_END));
}

//---------------------------------------------------------------------------

void MATFile::EndStream(){
	if(ReadStrm==1)strm_error=inflateEnd(&strm);
	else{
		cout << "(end) strm_error = " << strm_error << endl;
		cout << "strm_idx = " << strm_idx << endl;
		cout << "strm.avail_in = " << strm.avail_in << endl;
		while(strm_error==Z_OK){
			cout << "strm.avail_out = " << strm.avail_out << endl;
					cout << "strm.next_in = " << (int)strm.next_in << endl;
					cout << "strm_in      = " << (int)strm_in << endl;
					cout << "strm.next_out = " << (int)strm.next_out << endl;
					cout << "strm_out      = " << (int)strm_out << endl;
			strm.next_in = strm_in; //why necessary??
			strm_error = deflate(&strm, Z_FINISH);
					cout << "strm.next_in = " << (int)strm.next_in << endl;
					cout << "strm_in      = " << (int)strm_in << endl;
					cout << "strm.next_out = " << (int)strm.next_out << endl;
					cout << "strm_out      = " << (int)strm_out << endl;
			cout << "strm.avail_in = " << strm.avail_in << endl;
			cout << "strm.avail_out = " << strm.avail_out << endl;
			strm_idx = strm.avail_in;
			strm.next_out = strm_out; //why necessary??
			cout << File1.is_open() << " ; " << File1.fail() << " ; " << File1.bad() << " ; "  << File1.eof() << endl;
			cout << "pos = " << File1.tellp() << endl;
			cout << "add = " << CHUNK-strm.avail_out << endl;
			File1.write(reinterpret_cast<char*>(strm_out),CHUNK-strm.avail_out);
			cout << "now = " << File1.tellp() << endl;
					cout << "strm.next_in = " << (int)strm.next_in << endl;
					cout << "strm_in      = " << (int)strm_in << endl;
					cout << "strm.next_out = " << (int)strm.next_out << endl;
					cout << "strm_out      = " << (int)strm_out << endl;
			cout << "strm.avail_out = " << strm.avail_out << endl;
			strm.avail_out = CHUNK;
		}
		if(write_pos>0){
			//set file to begin of the stream
			//read stream
			//make changes at stored positions
			//write stream
			//strm_error=deflateReset(&strm);
			for(int j=0; j<write_pos; j++)
				cout << "pos[" << j << "] = (" << store_pos[j][0] << "," << store_pos[j][1] << ")" << endl;
		}
		strm_error=deflateEnd(&strm);
	}
	ReadStrm=0;
	strm_error=Z_OK;
}

//---------------------------------------------------------------------------

int MATFile::GetPos(){
	int ans;
	if(ReadStrm>0)ans=strm_pos+strm_begin;
	else{
		if(RWmode)ans=File1.tellg();
		else ans=File1.tellp();
	}
//	cout << "Pos (" << ReadStrm << ") = " << ans << endl;
	return ans;
}
//---------------------------------------------------------------------------

void MATFile::SetPos(int pos){
	int pos_tmp;
	if(ReadStrm>0){
		if((strm_pos-strm_idx)>(pos-strm_begin)){
			if(RWmode){
				strm_error=inflateReset(&strm);
				strm_idx=0;
				strm_pos=0;
				strm.next_in=strm_in;
				strm.avail_in=0;
				strm.next_out=strm_out;
				strm.avail_out=CHUNK;
			}
			else{
				cout << "SetPos in write stream is not possible!" << endl;
			}
		}
		pos_tmp = pos-strm_begin-strm_pos;
		Read(NULL,pos_tmp);
	}
	else{
		if(RWmode){
			if(pos>=0)File1.seekg(pos,ios_base::beg);
			else File1.seekg(-pos-1,ios_base::end);
		}
		else{
			if(pos>=0)File1.seekp(pos,ios_base::beg);
			else File1.seekp(-pos-1,ios_base::end);
		}
	}
}


//---------------------------------------------------------------------------

void MATFile::SetNindex(int n1){
	Nindex = n1;
}

//---------------------------------------------------------------------------

void MATFile::SetRWmode(bool RW){
	RWmode = RW;
}

//---------------------------------------------------------------------------

int MATFile::GetNindex(){
	return Nindex;
}

//---------------------------------------------------------------------------

bool MATFile::CheckHeader(){
	bool check1=true;
	char buffer [117];
	char buffer2[] = "MATLAB 5.0 MAT-file,";
	File1 >> noskipws;
	File1.get(buffer,117,EOF);
	for (int i=0; i<20; i++)
		if(buffer[i]!=buffer2[i])check1=false;
	File1.get(buffer,9,EOF);//8 more bytes subsys data offset
	char *a1;
	a1 = new char [2*sizeof(unsigned short int)];
	unsigned short int *c1;
	File1.read( a1, 2*sizeof(unsigned short int) );
	c1 = reinterpret_cast<unsigned short int*>(a1);
	if((c1[0]==256)||(c1[0]==1)){
		if(c1[0]==1){
			if(c1[1]!=18765)check1=false;
			else miRead = true;
		}
		if(c1[0]==256){
			if(c1[1]!=19785)check1=false;
			else miRead = false;
		}
	}
	else check1=false;
	delete[] a1;
	return check1;
}

//---------------------------------------------------------------------------

void MATFile::WriteHeader(){
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	string tmp = "MATLAB 5.0 MAT-file, Platform: PCWIN, Created on: ";
	File1 << tmp;
	tmp = asctime(timeinfo);
	File1 << tmp;
	tmp = " Self-defined!";
	File1 << tmp;
	int p1 = File1.tellp();
	tmp = " ";
	for (int i=p1; i<124; i++)
		File1 << tmp;
	unsigned short int c1[2];
	c1[0] = (unsigned short int)(256);
	c1[1] = (unsigned short int)(19785);
	File1.write( reinterpret_cast<char*>(c1), 2*sizeof(unsigned short int) );
	if(mode==1)SetPos(-1);
}

//---------------------------------------------------------------------------

void MATFile::Read(char *data, int& Ndata){
	if(ReadStrm==1){
		//compressed  stream
		int data_end, data_total=0;
		while(data_total<Ndata){
			if((CHUNK-strm.avail_out-strm_idx)>=Ndata)
				data_end = strm_idx + Ndata;
			else
				data_end = CHUNK-strm.avail_out;
			if(data!=NULL)
				for(int j=strm_idx; j<data_end; j++)
					data[j-strm_idx] = (char)strm_out[j];
			strm_pos += (data_end-strm_idx);
			data_total += (data_end-strm_idx);
			strm_idx = data_end;
			if(data_total<Ndata){
				if(strm.avail_in==0){
					strm.avail_in = (Npos-File1.tellg());
					if(strm.avail_in>CHUNK)strm.avail_in=CHUNK;
					File1.read(reinterpret_cast<char*>(strm_in),strm.avail_in);
					strm.next_in = strm_in; //probably not necessary
				}
				strm.avail_out = CHUNK;
				strm.next_out = strm_out; //probably not necessary
				strm_error = inflate(&strm, Z_NO_FLUSH);
				strm_idx = 0;
				if((strm_error!=Z_OK)&&(strm_error!=Z_STREAM_END))
					Ndata = data_total;
			}
			else
				if(strm_error==Z_STREAM_END)
					Ndata = data_total;
		}
	}
	else{
		File1.read(data, Ndata);
	}
}

//---------------------------------------------------------------------------

void MATFile::Write(char *data, int Ndata){
	if(ReadStrm==2){
		//compressed  stream
		int data_end, data_total=0;
//		cout << "Ndata = " << Ndata << endl;
		while(data_total<Ndata){
			if((CHUNK-strm_idx)>=Ndata)
				data_end = strm_idx + Ndata;
			else
				data_end = CHUNK;
//			cout << "data_end = " << data_end << endl;
//			cout << "data_total = " << data_total << endl;
			for(int j=strm_idx; j<data_end; j++)
				strm_in[j] = (unsigned char)data[j-strm_idx];
			strm_pos += (data_end-strm_idx);
			data_total += (data_end-strm_idx);
			strm_idx = data_end;
			strm.avail_in = strm_idx;
			if(data_total<Ndata){
//				cout << "strm.next_in = " << (int)strm.next_in << endl;
//				cout << "strm_in      = " << (int)strm_in << endl;
				strm.next_in = strm_in; //probably not necessary
//				cout << "strm.avail_in = " << strm.avail_in << endl;
				while(strm.avail_in!=0){
//					cout << "strm.next_out = " << (int)strm.next_out << endl;
//					cout << "strm_out      = " << (int)strm_out << endl;
					strm_error = deflate(&strm, Z_NO_FLUSH);
//					cout << "strm.next_in = " << (int)strm.next_in << endl;
//					cout << "strm_in      = " << (int)strm_in << endl;
//					cout << "strm.next_out = " << (int)strm.next_out << endl;
//					cout << "strm_out      = " << (int)strm_out << endl;
					strm.avail_in = strm_idx;
					strm.next_out = strm_out; //probably not necessary??
//					cout << "strm.avail_out = " << strm.avail_out << endl;
					File1.write(reinterpret_cast<char*>(strm_out),CHUNK-strm.avail_out);
//					cout << "strm.next_out = " << (int)strm.next_out << endl;
//					cout << "strm_out      = " << (int)strm_out << endl;
//					cout << "strm.avail_out = " << strm.avail_out << endl;
					strm.avail_out = CHUNK;
					if((strm_error!=Z_OK)&&(strm_error!=Z_STREAM_END)){
						strm.avail_in = 0;
						Ndata = data_total;
					}
				}
				strm_idx = 0;
				strm.avail_in = strm_idx;
			}
			else
				if(strm_error==Z_STREAM_END)
					Ndata=data_total;
		}
	}
	else{
		File1.write(data, Ndata);
	}
}


//---------------------------------------------------------------------------

bool MATFile::ReadTag(miMatrix& DataType, int& Bytes, bool& small){
	Bytes = 0;
	DataType = miMATRIX;
	small = false;
	int tmp;
	unsigned int* c1;
	char a1[4], a2[4];
	bool failed=false;
	if(!Fail()){
		tmp=4;
		Read(a1,tmp);
		if(miRead)for(int i=0; i<4; i++)a2[i] = a1[3-i];
		else for(int i=0; i<4; i++)a2[i] = a1[i];
		c1 = reinterpret_cast<unsigned int*>(a2);
		
		if(c1[0]>65535){//small data element
			small=true;
			DataType = (miMatrix)(c1[0]%65536);
			Bytes = c1[0]/65536;
		}
		else{
			small=false;
			DataType = (miMatrix)c1[0];
			
			tmp=4;
			Read(a1,tmp);
			if(miRead)for(int i=0; i<4; i++)a2[i] = a1[3-i];
			else for(int i=0; i<4; i++)a2[i] = a1[i];
			c1 = reinterpret_cast<unsigned int*>(a2);
			Bytes = c1[0];
		}
	}
	else failed = true;
	return failed;
}

//---------------------------------------------------------------------------

bool MATFile::WriteTag(miMatrix DataType, int Bytes, bool& small){
	unsigned int* c1;
	bool failed=false;
	if(!Fail()){
		if((DataType==miCOMPRESSED)||(DataType==miMATRIX)){
			if(!small){
				c1 = new unsigned int [2];
				c1[0] = DataType;
				c1[1] = 0;
				Write(reinterpret_cast<char*>(c1),8);
			}
			else{
				int Npos = GetPos();
				if(ReadStrm>0){
					if((strm_pos-strm_idx)<(Bytes-4)){
						c1 = new unsigned int [1];
						c1[0] = Npos-Bytes;
						char *a1;
						a1 = reinterpret_cast<char*>(c1);
						int j_start = (Bytes-4)-(strm_pos-strm_idx)-strm_begin;
						for(int j=j_start; j<(j_start+4); j++){
							strm_in[j] = (unsigned char)a1[j-j_start];
						}
					}
					else{//Store and write later
//						cout << "Bytes = " << Bytes << endl;
						store_pos[write_pos][0]=Bytes-4-strm_begin;//deflated stream pos
						store_pos[write_pos][1]=Npos-Bytes;
						write_pos++;
					}
				}
				else{
					SetPos(Bytes-4);
					c1 = new unsigned int [1];
					c1[0] = Npos-Bytes;
					Write(reinterpret_cast<char*>(c1),4);
					SetPos(Npos);
				}
			}
		}
		else{
			if((Bytes>4)||(Bytes==0)){
				small=false;
				c1 = new unsigned int [2];
				c1[0] = DataType;
				c1[1] = Bytes;
				Write(reinterpret_cast<char*>(c1),8);
			}
			else{
				small=true;
				c1 = new unsigned int [1];
				c1[0] = DataType+Bytes*65536;
				Write(reinterpret_cast<char*>(c1),4);
			}
		}
		delete[] c1;
	}
	else failed = true;
	return failed;
}

//---------------------------------------------------------------------------

bool MATFile::ReadDataElements(mxClassID ArrayType, miMatrix DataType, void* Data, int Bytes, int size, bool small){
	mxLogical* d00;
	char* d01;
	double* d02;
	float* d03;
	int8s* d04;
	int8u* d05;
	int16s* d06;
	int16u* d07;
	int32s* d08;
	int32u* d09;
	int64s* d10;
	int64u* d11;
	switch(ArrayType){
		case mxUNKNOWN_CLASS : break; //not used in this function
		case mxCELL_CLASS    : break; //not used in this function
		case mxSTRUCT_CLASS  : break; //not used in this function
		case mxLOGICAL_CLASS : d00 = reinterpret_cast<mxLogical*>(Data); break;
		case mxCHAR_CLASS    : d01 = reinterpret_cast<char*>(Data); break;
		case mxVOID_CLASS    : d02 = reinterpret_cast<double*>(Data); break;
		case mxDOUBLE_CLASS  : d02 = reinterpret_cast<double*>(Data); break;
		case mxSINGLE_CLASS  : d03 = reinterpret_cast<float*>(Data); break;
		case mxINT8_CLASS    : d04 = reinterpret_cast<int8s*>(Data); break;
		case mxUINT8_CLASS   : d05 = reinterpret_cast<int8u*>(Data); break;
		case mxINT16_CLASS   : d06 = reinterpret_cast<int16s*>(Data); break;
		case mxUINT16_CLASS  : d07 = reinterpret_cast<int16u*>(Data); break;
		case mxINT32_CLASS   : d08 = reinterpret_cast<int32s*>(Data); break;
		case mxUINT32_CLASS  : d09 = reinterpret_cast<int32u*>(Data); break;
		case mxINT64_CLASS   : d10 = reinterpret_cast<int64s*>(Data); break;
		case mxUINT64_CLASS  : d11 = reinterpret_cast<int64u*>(Data); break;
		case mxFUNCTION_CLASS: break; //not used in this function
		case mxOPAQUE_CLASS  : break; //not used in this function
		case mxOBJECT_CLASS  : break; //not used in this function
	}

	int8s* c01;
	int8u* c02;
	int16s* c03;
	int16u* c04;
	int32s* c05;
	int32u* c06;
	float* c07;
	double* c09;
	int64s* c12;
	int64u* c13;
	int8u* c16;
	int16u* c17;
	int32u* c18;
	char *a1, *a2;
	a1 = new char[size];
	a2 = new char[size];
	int Elements = Bytes/size;
	for(int j=0; j<Elements; j++){
		Read(a1,size);
		if(miRead)for(int i=0; i<size; i++)a2[i] = a1[size-i-1];
		else for(int i=0; i<size; i++)a2[i] = a1[i];
		switch(DataType){
			case miINT8   : c01 = reinterpret_cast<int8s*>(a2);
			                switch(ArrayType){
			                	case mxLOGICAL_CLASS: d00[j] = (mxLogical)c01[0]; break;
			                	case mxCHAR_CLASS   : d01[j] = (char)c01[0]; break;
			                	case mxVOID_CLASS   :
			                	case mxDOUBLE_CLASS : d02[j] = (double)c01[0]; break;
			                	case mxSINGLE_CLASS : d03[j] = (float)c01[0]; break;
			                	case mxINT8_CLASS   : d04[j] = c01[0]; break;
			                	case mxUINT8_CLASS  : d05[j] = (int8u)c01[0]; break;
			                	case mxINT16_CLASS  : d06[j] = (int16s)c01[0]; break;
			                	case mxUINT16_CLASS : d07[j] = (int16u)c01[0]; break;
			                	case mxINT32_CLASS  : d08[j] = (int32s)c01[0]; break;
			                	case mxUINT32_CLASS : d09[j] = (int32u)c01[0]; break;
			                	case mxINT64_CLASS  : d10[j] = (int64s)c01[0]; break;
			                	case mxUINT64_CLASS : d11[j] = (int64u)c01[0]; break;
			                }; break;
			case miUINT8  : c02 = reinterpret_cast<int8u*>(a2);
			                switch(ArrayType){
			                	case mxLOGICAL_CLASS: d00[j] = (mxLogical)c02[0]; break;
			                	case mxCHAR_CLASS   : d01[j] = (char)c02[0]; break;
			                	case mxVOID_CLASS   :
			                	case mxDOUBLE_CLASS : d02[j] = (double)c02[0]; break;
			                	case mxSINGLE_CLASS : d03[j] = (float)c02[0]; break;
			                	case mxINT8_CLASS   : d04[j] = (int8s)c02[0]; break;
			                	case mxUINT8_CLASS  : d05[j] = c02[0]; break;
			                	case mxINT16_CLASS  : d06[j] = (int16s)c02[0]; break;
			                	case mxUINT16_CLASS : d07[j] = (int16u)c02[0]; break;
			                	case mxINT32_CLASS  : d08[j] = (int32s)c02[0]; break;
			                	case mxUINT32_CLASS : d09[j] = (int32u)c02[0]; break;
			                	case mxINT64_CLASS  : d10[j] = (int64s)c02[0]; break;
			                	case mxUINT64_CLASS : d11[j] = (int64u)c02[0]; break;
			                }; break;
			case miINT16  : c03 = reinterpret_cast<int16s*>(a2);
			                switch(ArrayType){
			                	case mxLOGICAL_CLASS: d00[j] = (mxLogical)c03[0]; break;
			                	case mxCHAR_CLASS   : d01[j] = (char)c03[0]; break;
			                	case mxVOID_CLASS   :
			                	case mxDOUBLE_CLASS : d02[j] = (double)c03[0]; break;
			                	case mxSINGLE_CLASS : d03[j] = (float)c03[0]; break;
			                	case mxINT8_CLASS   : d04[j] = (int8s)c03[0]; break;
			                	case mxUINT8_CLASS  : d05[j] = (int8u)c03[0]; break;
			                	case mxINT16_CLASS  : d06[j] = c03[0]; break;
			                	case mxUINT16_CLASS : d07[j] = (int16u)c03[0]; break;
			                	case mxINT32_CLASS  : d08[j] = (int32s)c03[0]; break;
			                	case mxUINT32_CLASS : d09[j] = (int32u)c03[0]; break;
			                	case mxINT64_CLASS  : d10[j] = (int64s)c03[0]; break;
			                	case mxUINT64_CLASS : d11[j] = (int64u)c03[0]; break;
			                }; break;
			case miUINT16 : c04 = reinterpret_cast<int16u*>(a2);
			                switch(ArrayType){
			                	case mxLOGICAL_CLASS: d00[j] = (mxLogical)c04[0]; break;
			                	case mxCHAR_CLASS   : d01[j] = (char)c04[0]; break;
			                	case mxVOID_CLASS   :
			                	case mxDOUBLE_CLASS : d02[j] = (double)c04[0]; break;
			                	case mxSINGLE_CLASS : d03[j] = (float)c04[0]; break;
			                	case mxINT8_CLASS   : d04[j] = (int8s)c04[0]; break;
			                	case mxUINT8_CLASS  : d05[j] = (int8u)c04[0]; break;
			                	case mxINT16_CLASS  : d06[j] = (int16s)c04[0]; break;
			                	case mxUINT16_CLASS : d07[j] = c04[0]; break;
			                	case mxINT32_CLASS  : d08[j] = (int32s)c04[0]; break;
			                	case mxUINT32_CLASS : d09[j] = (int32u)c04[0]; break;
			                	case mxINT64_CLASS  : d10[j] = (int64s)c04[0]; break;
			                	case mxUINT64_CLASS : d11[j] = (int64u)c04[0]; break;
			                }; break;
			case miINT32  : c05 = reinterpret_cast<int32s*>(a2);
			                switch(ArrayType){
			                	case mxLOGICAL_CLASS: d00[j] = (mxLogical)c05[0]; break;
			                	case mxCHAR_CLASS   : d01[j] = (char)c05[0]; break;
			                	case mxVOID_CLASS   :
			                	case mxDOUBLE_CLASS : d02[j] = (double)c05[0]; break;
			                	case mxSINGLE_CLASS : d03[j] = (float)c05[0]; break;
			                	case mxINT8_CLASS   : d04[j] = (int8s)c05[0]; break;
			                	case mxUINT8_CLASS  : d05[j] = (int8u)c05[0]; break;
			                	case mxINT16_CLASS  : d06[j] = (int16s)c05[0]; break;
			                	case mxUINT16_CLASS : d07[j] = (int16u)c05[0]; break;
			                	case mxINT32_CLASS  : d08[j] = c05[0]; break;
			                	case mxUINT32_CLASS : d09[j] = (int32u)c05[0]; break;
			                	case mxINT64_CLASS  : d10[j] = (int64s)c05[0]; break;
			                	case mxUINT64_CLASS : d11[j] = (int64u)c05[0]; break;
			                }; break;
			case miUINT32 : c06 = reinterpret_cast<int32u*>(a2);
			                switch(ArrayType){
			                	case mxLOGICAL_CLASS: d00[j] = (mxLogical)c06[0]; break;
			                	case mxCHAR_CLASS   : d01[j] = (char)c06[0]; break;
			                	case mxVOID_CLASS   :
			                	case mxDOUBLE_CLASS : d02[j] = (double)c06[0]; break;
			                	case mxSINGLE_CLASS : d03[j] = (float)c06[0]; break;
			                	case mxINT8_CLASS   : d04[j] = (int8s)c06[0]; break;
			                	case mxUINT8_CLASS  : d05[j] = (int8u)c06[0]; break;
			                	case mxINT16_CLASS  : d06[j] = (int16s)c06[0]; break;
			                	case mxUINT16_CLASS : d07[j] = (int16u)c06[0]; break;
			                	case mxINT32_CLASS  : d08[j] = (int32s)c06[0]; break;
			                	case mxUINT32_CLASS : d09[j] = c06[0]; break;
			                	case mxINT64_CLASS  : d10[j] = (int64s)c06[0]; break;
			                	case mxUINT64_CLASS : d11[j] = (int64u)c06[0]; break;
			                }; break;
			case miSINGLE : c07 = reinterpret_cast<float*>(a2);
			                switch(ArrayType){
			                	case mxLOGICAL_CLASS: d00[j] = (mxLogical)c07[0]; break;
			                	case mxCHAR_CLASS   : d01[j] = (char)c07[0]; break;
			                	case mxVOID_CLASS   :
			                	case mxDOUBLE_CLASS : d02[j] = (double)c07[0]; break;
			                	case mxSINGLE_CLASS : d03[j] = c07[0]; break;
			                	case mxINT8_CLASS   : d04[j] = (int8s)c07[0]; break;
			                	case mxUINT8_CLASS  : d05[j] = (int8u)c07[0]; break;
			                	case mxINT16_CLASS  : d06[j] = (int16s)c07[0]; break;
			                	case mxUINT16_CLASS : d07[j] = (int16u)c07[0]; break;
			                	case mxINT32_CLASS  : d08[j] = (int32s)c07[0]; break;
			                	case mxUINT32_CLASS : d09[j] = (int32u)c07[0]; break;
			                	case mxINT64_CLASS  : d10[j] = (int64s)c07[0]; break;
			                	case mxUINT64_CLASS : d11[j] = (int64u)c07[0]; break;
			                }; break;
			case miDOUBLE : c09 = reinterpret_cast<double*>(a2);
			                switch(ArrayType){
			                	case mxLOGICAL_CLASS: d00[j] = (mxLogical)c09[0]; break;
			                	case mxCHAR_CLASS   : d01[j] = (char)c09[0]; break;
			                	case mxVOID_CLASS   :
			                	case mxDOUBLE_CLASS : d02[j] = c09[0]; break;
			                	case mxSINGLE_CLASS : d03[j] = (float)c09[0]; break;
			                	case mxINT8_CLASS   : d04[j] = (int8s)c09[0]; break;
			                	case mxUINT8_CLASS  : d05[j] = (int8u)c09[0]; break;
			                	case mxINT16_CLASS  : d06[j] = (int16s)c09[0]; break;
			                	case mxUINT16_CLASS : d07[j] = (int16u)c09[0]; break;
			                	case mxINT32_CLASS  : d08[j] = (int32s)c09[0]; break;
			                	case mxUINT32_CLASS : d09[j] = (int32u)c09[0]; break;
			                	case mxINT64_CLASS  : d10[j] = (int64s)c09[0]; break;
			                	case mxUINT64_CLASS : d11[j] = (int64u)c09[0]; break;
			                }; break;
			case miINT64  : c12 = reinterpret_cast<int64s*>(a2);
			                switch(ArrayType){
			                	case mxLOGICAL_CLASS: d00[j] = (mxLogical)c12[0]; break;
			                	case mxCHAR_CLASS   : d01[j] = (char)c12[0]; break;
			                	case mxVOID_CLASS   :
			                	case mxDOUBLE_CLASS : d02[j] = (double)c12[0]; break;
			                	case mxSINGLE_CLASS : d03[j] = (float)c12[0]; break;
			                	case mxINT8_CLASS   : d04[j] = (int8s)c12[0]; break;
			                	case mxUINT8_CLASS  : d05[j] = (int8u)c12[0]; break;
			                	case mxINT16_CLASS  : d06[j] = (int16s)c12[0]; break;
			                	case mxUINT16_CLASS : d07[j] = (int16u)c12[0]; break;
			                	case mxINT32_CLASS  : d08[j] = (int32s)c12[0]; break;
			                	case mxUINT32_CLASS : d09[j] = (int32u)c12[0]; break;
			                	case mxINT64_CLASS  : d10[j] = c12[0]; break;
			                	case mxUINT64_CLASS : d11[j] = (int64u)c12[0]; break;
			                }; break;
			case miUINT64 : c13 = reinterpret_cast<int64u*>(a2);
			                switch(ArrayType){
			                	case mxLOGICAL_CLASS: d00[j] = (mxLogical)c13[0]; break;
			                	case mxCHAR_CLASS   : d01[j] = (char)c13[0]; break;
			                	case mxVOID_CLASS   :
			                	case mxDOUBLE_CLASS : d02[j] = (double)c13[0]; break;
			                	case mxSINGLE_CLASS : d03[j] = (float)c13[0]; break;
			                	case mxINT8_CLASS   : d04[j] = (int8s)c13[0]; break;
			                	case mxUINT8_CLASS  : d05[j] = (int8u)c13[0]; break;
			                	case mxINT16_CLASS  : d06[j] = (int16s)c13[0]; break;
			                	case mxUINT16_CLASS : d07[j] = (int16u)c13[0]; break;
			                	case mxINT32_CLASS  : d08[j] = (int32s)c13[0]; break;
			                	case mxUINT32_CLASS : d09[j] = (int32u)c13[0]; break;
			                	case mxINT64_CLASS  : d10[j] = (int64s)c13[0]; break;
			                	case mxUINT64_CLASS : d11[j] = c13[0]; break;
			                }; break;
			case miUTF8   : c16 = reinterpret_cast<int8u*>(a2);
			                switch(ArrayType){
			                	case mxLOGICAL_CLASS: d00[j] = (mxLogical)c16[0]; break;
			                	case mxCHAR_CLASS   : d01[j] = (char)c16[0]; break;
			                	case mxVOID_CLASS   :
			                	case mxDOUBLE_CLASS : d02[j] = (double)c16[0]; break;
			                	case mxSINGLE_CLASS : d03[j] = (float)c16[0]; break;
			                	case mxINT8_CLASS   : d04[j] = (int8s)c16[0]; break;
			                	case mxUINT8_CLASS  : d05[j] = c16[0]; break;
			                	case mxINT16_CLASS  : d06[j] = (int16s)c16[0]; break;
			                	case mxUINT16_CLASS : d07[j] = (int16u)c16[0]; break;
			                	case mxINT32_CLASS  : d08[j] = (int32s)c16[0]; break;
			                	case mxUINT32_CLASS : d09[j] = (int32u)c16[0]; break;
			                	case mxINT64_CLASS  : d10[j] = (int64s)c16[0]; break;
			                	case mxUINT64_CLASS : d11[j] = (int64u)c16[0]; break;
			                }; break;
			case miUTF16  : c17 = reinterpret_cast<int16u*>(a2);
			                switch(ArrayType){
			                	case mxLOGICAL_CLASS: d00[j] = (mxLogical)c17[0]; break;
			                	case mxCHAR_CLASS   : d01[j] = (char)c17[0]; break;
			                	case mxVOID_CLASS   :
			                	case mxDOUBLE_CLASS : d02[j] = (double)c17[0]; break;
			                	case mxSINGLE_CLASS : d03[j] = (float)c17[0]; break;
			                	case mxINT8_CLASS   : d04[j] = (int8s)c17[0]; break;
			                	case mxUINT8_CLASS  : d05[j] = (int8u)c17[0]; break;
			                	case mxINT16_CLASS  : d06[j] = (int16s)c17[0]; break;
			                	case mxUINT16_CLASS : d07[j] = c17[0]; break;
			                	case mxINT32_CLASS  : d08[j] = (int32s)c17[0]; break;
			                	case mxUINT32_CLASS : d09[j] = (int32u)c17[0]; break;
			                	case mxINT64_CLASS  : d10[j] = (int64s)c17[0]; break;
			                	case mxUINT64_CLASS : d11[j] = (int64u)c17[0]; break;
			                }; break;
			case miUTF32  : c18 = reinterpret_cast<int32u*>(a2);
			                switch(ArrayType){
			                	case mxLOGICAL_CLASS: d00[j] = (mxLogical)c18[0]; break;
			                	case mxCHAR_CLASS   : d01[j] = (char)c18[0]; break;
			                	case mxVOID_CLASS   :
			                	case mxDOUBLE_CLASS : d02[j] = (double)c18[0]; break;
			                	case mxSINGLE_CLASS : d03[j] = (float)c18[0]; break;
			                	case mxINT8_CLASS   : d04[j] = (int8s)c18[0]; break;
			                	case mxUINT8_CLASS  : d05[j] = (int8u)c18[0]; break;
			                	case mxINT16_CLASS  : d06[j] = (int16s)c18[0]; break;
			                	case mxUINT16_CLASS : d07[j] = (int16u)c18[0]; break;
			                	case mxINT32_CLASS  : d08[j] = (int32s)c18[0]; break;
			                	case mxUINT32_CLASS : d09[j] = c18[0]; break;
			                	case mxINT64_CLASS  : d10[j] = (int64s)c18[0]; break;
			                	case mxUINT64_CLASS : d11[j] = (int64u)c18[0]; break;
			                }; break;
		}
	}
	//padding
	int padding;
	if(!small){
		padding=(Elements*size)/8;
		if(((Elements*size)%8)!=0)padding++;
		padding=padding*8;
		padding=padding-(Elements*size);
	}
	else{
		padding=(Elements*size)/4;
		if(((Elements*size)%4)!=0)padding++;
		padding=padding*4;
		padding=padding-(Elements*size);
	}
	for(int i=0; i<(padding/size); i++)
		Read(a1,size);
	delete[] a1;
	delete[] a2;
	return !Fail();
}

//---------------------------------------------------------------------------
//for writing the compressed data format, this should be checked first and passed on to this function
bool MATFile::WriteDataElements(mxClassID ArrayType, miMatrix DataType, void* Data, int Bytes, int size, bool small){
	mxLogical* d00;
	char* d01;
	double* d02;
	float* d03;
	int8s* d04;
	int8u* d05;
	int16s* d06;
	int16u* d07;
	int32s* d08;
	int32u* d09;
	int64s* d10;
	int64u* d11;
	switch(ArrayType){
		case mxUNKNOWN_CLASS : break; //not used in this function
		case mxCELL_CLASS    : break; //not used in this function
		case mxSTRUCT_CLASS  : break; //not used in this function
		case mxLOGICAL_CLASS : d00 = reinterpret_cast<mxLogical*>(Data); break;
		case mxCHAR_CLASS    : d01 = reinterpret_cast<char*>(Data); break;
		case mxVOID_CLASS    : d02 = reinterpret_cast<double*>(Data); break;
		case mxDOUBLE_CLASS  : d02 = reinterpret_cast<double*>(Data); break;
		case mxSINGLE_CLASS  : d03 = reinterpret_cast<float*>(Data); break;
		case mxINT8_CLASS    : d04 = reinterpret_cast<int8s*>(Data); break;
		case mxUINT8_CLASS   : d05 = reinterpret_cast<int8u*>(Data); break;
		case mxINT16_CLASS   : d06 = reinterpret_cast<int16s*>(Data); break;
		case mxUINT16_CLASS  : d07 = reinterpret_cast<int16u*>(Data); break;
		case mxINT32_CLASS   : d08 = reinterpret_cast<int32s*>(Data); break;
		case mxUINT32_CLASS  : d09 = reinterpret_cast<int32u*>(Data); break;
		case mxINT64_CLASS   : d10 = reinterpret_cast<int64s*>(Data); break;
		case mxUINT64_CLASS  : d11 = reinterpret_cast<int64u*>(Data); break;
		case mxFUNCTION_CLASS: break; //not used in this function
		case mxOPAQUE_CLASS  : break; //not used in this function
		case mxOBJECT_CLASS  : break; //not used in this function
	}
	int8s c01[1];
	int8u c02[1];
	int16s c03[1];
	int16u c04[1];
	int32s c05[1];
	int32u c06[1];
	float c07[1];
	double c09[1];
	int64s c12[1];
	int64u c13[1];
	int8u c16[1];
	int16u c17[1];
	int32u c18[1];
	for(int i=0; i<(Bytes/size); i++){
		switch(ArrayType){
			case mxLOGICAL_CLASS  : switch(DataType){
			                      	case miINT8   : c01[0] = (int8s)d00[i]; Write(reinterpret_cast<char*>(c01),size); break;
			                      	case miUINT8  : c02[0] = (int8u)d00[i]; Write(reinterpret_cast<char*>(c02),size); break;
			                      	case miINT16  : c03[0] = (int16s)d00[i]; Write(reinterpret_cast<char*>(c03),size); break;
			                      	case miUINT16 : c04[0] = (int16u)d00[i]; Write(reinterpret_cast<char*>(c04),size); break;
			                      	case miINT32  : c05[0] = (int32s)d00[i]; Write(reinterpret_cast<char*>(c05),size); break;
			                      	case miUINT32 : c06[0] = (int32u)d00[i]; Write(reinterpret_cast<char*>(c06),size); break;
			                      	case miSINGLE : c07[0] = (float)d00[i]; Write(reinterpret_cast<char*>(c07),size); break;
			                      	case miDOUBLE : c09[0] = (double)d00[i]; Write(reinterpret_cast<char*>(c09),size); break;
			                      	case miINT64  : c12[0] = (int64s)d00[i]; Write(reinterpret_cast<char*>(c12),size); break;
			                      	case miUINT64 : c13[0] = (int64u)d00[i]; Write(reinterpret_cast<char*>(c13),size); break;
			                      	case miUTF8   : c16[0] = (int8u)d00[i]; Write(reinterpret_cast<char*>(c16),size); break;
			                      	case miUTF16  : c17[0] = (int16u)d00[i]; Write(reinterpret_cast<char*>(c17),size); break;
			                      	case miUTF32  : c18[0] = (int32u)d00[i]; Write(reinterpret_cast<char*>(c18),size); break;
			                      }; break;
			case mxCHAR_CLASS   : switch(DataType){
			                      	case miINT8   : c01[0] = (int8s)d01[i]; Write(reinterpret_cast<char*>(c01),size); break;
			                      	case miUINT8  : c02[0] = (int8u)d01[i]; Write(reinterpret_cast<char*>(c02),size); break;
			                      	case miINT16  : c03[0] = (int16s)d01[i]; Write(reinterpret_cast<char*>(c03),size); break;
			                      	case miUINT16 : c04[0] = (int16u)d01[i]; Write(reinterpret_cast<char*>(c04),size); break;
			                      	case miINT32  : c05[0] = (int32s)d01[i]; Write(reinterpret_cast<char*>(c05),size); break;
			                      	case miUINT32 : c06[0] = (int32u)d01[i]; Write(reinterpret_cast<char*>(c06),size); break;
			                      	case miSINGLE : c07[0] = (float)d01[i]; Write(reinterpret_cast<char*>(c07),size); break;
			                      	case miDOUBLE : c09[0] = (double)d01[i]; Write(reinterpret_cast<char*>(c09),size); break;
			                      	case miINT64  : c12[0] = (int64s)d01[i]; Write(reinterpret_cast<char*>(c12),size); break;
			                      	case miUINT64 : c13[0] = (int64u)d01[i]; Write(reinterpret_cast<char*>(c13),size); break;
			                      	case miUTF8   : c16[0] = (int8u)d01[i]; Write(reinterpret_cast<char*>(c16),size); break;
			                      	case miUTF16  : c17[0] = (int16u)d01[i]; Write(reinterpret_cast<char*>(c17),size); break;
			                      	case miUTF32  : c18[0] = (int32u)d01[i]; Write(reinterpret_cast<char*>(c18),size); break;
			                      }; break;
			case mxVOID_CLASS   :
			case mxDOUBLE_CLASS : switch(DataType){
			                      	case miINT8   : c01[0] = (int8s)d02[i]; Write(reinterpret_cast<char*>(c01),size); break;
			                      	case miUINT8  : c02[0] = (int8u)d02[i]; Write(reinterpret_cast<char*>(c02),size); break;
			                      	case miINT16  : c03[0] = (int16s)d02[i]; Write(reinterpret_cast<char*>(c03),size); break;
			                      	case miUINT16 : c04[0] = (int16u)d02[i]; Write(reinterpret_cast<char*>(c04),size); break;
			                      	case miINT32  : c05[0] = (int32s)d02[i]; Write(reinterpret_cast<char*>(c05),size); break;
			                      	case miUINT32 : c06[0] = (int32u)d02[i]; Write(reinterpret_cast<char*>(c06),size); break;
			                      	case miSINGLE : c07[0] = (float)d02[i]; Write(reinterpret_cast<char*>(c07),size); break;
			                      	case miDOUBLE : c09[0] = d02[i]; Write(reinterpret_cast<char*>(c09),size); break;
			                      	case miINT64  : c12[0] = (int64s)d02[i]; Write(reinterpret_cast<char*>(c12),size); break;
			                      	case miUINT64 : c13[0] = (int64u)d02[i]; Write(reinterpret_cast<char*>(c13),size); break;
			                      	case miUTF8   : c16[0] = (int8u)d02[i]; Write(reinterpret_cast<char*>(c16),size); break;
			                      	case miUTF16  : c17[0] = (int16u)d02[i]; Write(reinterpret_cast<char*>(c17),size); break;
			                      	case miUTF32  : c18[0] = (int32u)d02[i]; Write(reinterpret_cast<char*>(c18),size); break;
			                      }; break;
			case mxSINGLE_CLASS : switch(DataType){
			                      	case miINT8   : c01[0] = (int8s)d03[i]; Write(reinterpret_cast<char*>(c01),size); break;
			                      	case miUINT8  : c02[0] = (int8u)d03[i]; Write(reinterpret_cast<char*>(c02),size); break;
			                      	case miINT16  : c03[0] = (int16s)d03[i]; Write(reinterpret_cast<char*>(c03),size); break;
			                      	case miUINT16 : c04[0] = (int16u)d03[i]; Write(reinterpret_cast<char*>(c04),size); break;
			                      	case miINT32  : c05[0] = (int32s)d03[i]; Write(reinterpret_cast<char*>(c05),size); break;
			                      	case miUINT32 : c06[0] = (int32u)d03[i]; Write(reinterpret_cast<char*>(c06),size); break;
			                      	case miSINGLE : c07[0] = d03[i]; Write(reinterpret_cast<char*>(c07),size); break;
			                      	case miDOUBLE : c09[0] = (double)d03[i]; Write(reinterpret_cast<char*>(c09),size); break;
			                      	case miINT64  : c12[0] = (int64s)d03[i]; Write(reinterpret_cast<char*>(c12),size); break;
			                      	case miUINT64 : c13[0] = (int64u)d03[i]; Write(reinterpret_cast<char*>(c13),size); break;
			                      	case miUTF8   : c16[0] = (int8u)d03[i]; Write(reinterpret_cast<char*>(c16),size); break;
			                      	case miUTF16  : c17[0] = (int16u)d03[i]; Write(reinterpret_cast<char*>(c17),size); break;
			                      	case miUTF32  : c18[0] = (int32u)d03[i]; Write(reinterpret_cast<char*>(c18),size); break;
			                      }; break;
			case mxINT8_CLASS   : switch(DataType){
			                      	case miINT8   : c01[0] = d04[i]; Write(reinterpret_cast<char*>(c01),size); break;
			                      	case miUINT8  : c02[0] = (int8u)d04[i]; Write(reinterpret_cast<char*>(c02),size); break;
			                      	case miINT16  : c03[0] = (int16s)d04[i]; Write(reinterpret_cast<char*>(c03),size); break;
			                      	case miUINT16 : c04[0] = (int16u)d04[i]; Write(reinterpret_cast<char*>(c04),size); break;
			                      	case miINT32  : c05[0] = (int32s)d04[i]; Write(reinterpret_cast<char*>(c05),size); break;
			                      	case miUINT32 : c06[0] = (int32u)d04[i]; Write(reinterpret_cast<char*>(c06),size); break;
			                      	case miSINGLE : c07[0] = (float)d04[i]; Write(reinterpret_cast<char*>(c07),size); break;
			                      	case miDOUBLE : c09[0] = (double)d04[i]; Write(reinterpret_cast<char*>(c09),size); break;
			                      	case miINT64  : c12[0] = (int64s)d04[i]; Write(reinterpret_cast<char*>(c12),size); break;
			                      	case miUINT64 : c13[0] = (int64u)d04[i]; Write(reinterpret_cast<char*>(c13),size); break;
			                      	case miUTF8   : c16[0] = (int8u)d04[i]; Write(reinterpret_cast<char*>(c16),size); break;
			                      	case miUTF16  : c17[0] = (int16u)d04[i]; Write(reinterpret_cast<char*>(c17),size); break;
			                      	case miUTF32  : c18[0] = (int32u)d04[i]; Write(reinterpret_cast<char*>(c18),size); break;
			                      }; break;
			case mxUINT8_CLASS  : switch(DataType){
			                      	case miINT8   : c01[0] = (int8s)d05[i]; Write(reinterpret_cast<char*>(c01),size); break;
			                      	case miUINT8  : c02[0] = d05[i]; Write(reinterpret_cast<char*>(c02),size); break;
			                      	case miINT16  : c03[0] = (int16s)d05[i]; Write(reinterpret_cast<char*>(c03),size); break;
			                      	case miUINT16 : c04[0] = (int16u)d05[i]; Write(reinterpret_cast<char*>(c04),size); break;
			                      	case miINT32  : c05[0] = (int32s)d05[i]; Write(reinterpret_cast<char*>(c05),size); break;
			                      	case miUINT32 : c06[0] = (int32u)d05[i]; Write(reinterpret_cast<char*>(c06),size); break;
			                      	case miSINGLE : c07[0] = (float)d05[i]; Write(reinterpret_cast<char*>(c07),size); break;
			                      	case miDOUBLE : c09[0] = (double)d05[i]; Write(reinterpret_cast<char*>(c09),size); break;
			                      	case miINT64  : c12[0] = (int64s)d05[i]; Write(reinterpret_cast<char*>(c12),size); break;
			                      	case miUINT64 : c13[0] = (int64u)d05[i]; Write(reinterpret_cast<char*>(c13),size); break;
			                      	case miUTF8   : c16[0] = d05[i]; Write(reinterpret_cast<char*>(c16),size); break;
			                      	case miUTF16  : c17[0] = (int16u)d05[i]; Write(reinterpret_cast<char*>(c17),size); break;
			                      	case miUTF32  : c18[0] = (int32u)d05[i]; Write(reinterpret_cast<char*>(c18),size); break;
			                      }; break;
			case mxINT16_CLASS  : switch(DataType){
			                      	case miINT8   : c01[0] = (int8s)d06[i]; Write(reinterpret_cast<char*>(c01),size); break;
			                      	case miUINT8  : c02[0] = (int8u)d06[i]; Write(reinterpret_cast<char*>(c02),size); break;
			                      	case miINT16  : c03[0] = d06[i]; Write(reinterpret_cast<char*>(c03),size); break;
			                      	case miUINT16 : c04[0] = (int16u)d06[i]; Write(reinterpret_cast<char*>(c04),size); break;
			                      	case miINT32  : c05[0] = (int32s)d06[i]; Write(reinterpret_cast<char*>(c05),size); break;
			                      	case miUINT32 : c06[0] = (int32u)d06[i]; Write(reinterpret_cast<char*>(c06),size); break;
			                      	case miSINGLE : c07[0] = (float)d06[i]; Write(reinterpret_cast<char*>(c07),size); break;
			                      	case miDOUBLE : c09[0] = (double)d06[i]; Write(reinterpret_cast<char*>(c09),size); break;
			                      	case miINT64  : c12[0] = (int64s)d06[i]; Write(reinterpret_cast<char*>(c12),size); break;
			                      	case miUINT64 : c13[0] = (int64u)d06[i]; Write(reinterpret_cast<char*>(c13),size); break;
			                      	case miUTF8   : c16[0] = (int8u)d06[i]; Write(reinterpret_cast<char*>(c16),size); break;
			                      	case miUTF16  : c17[0] = (int16u)d06[i]; Write(reinterpret_cast<char*>(c17),size); break;
			                      	case miUTF32  : c18[0] = (int32u)d06[i]; Write(reinterpret_cast<char*>(c18),size); break;
			                      }; break;
			case mxUINT16_CLASS : switch(DataType){
			                      	case miINT8   : c01[0] = (int8s)d07[i]; Write(reinterpret_cast<char*>(c01),size); break;
			                      	case miUINT8  : c02[0] = (int8u)d07[i]; Write(reinterpret_cast<char*>(c02),size); break;
			                      	case miINT16  : c03[0] = (int16s)d07[i]; Write(reinterpret_cast<char*>(c03),size); break;
			                      	case miUINT16 : c04[0] = d07[i]; Write(reinterpret_cast<char*>(c04),size); break;
			                      	case miINT32  : c05[0] = (int32s)d07[i]; Write(reinterpret_cast<char*>(c05),size); break;
			                      	case miUINT32 : c06[0] = (int32u)d07[i]; Write(reinterpret_cast<char*>(c06),size); break;
			                      	case miSINGLE : c07[0] = (float)d07[i]; Write(reinterpret_cast<char*>(c07),size); break;
			                      	case miDOUBLE : c09[0] = (double)d07[i]; Write(reinterpret_cast<char*>(c09),size); break;
			                      	case miINT64  : c12[0] = (int64s)d07[i]; Write(reinterpret_cast<char*>(c12),size); break;
			                      	case miUINT64 : c13[0] = (int64u)d07[i]; Write(reinterpret_cast<char*>(c13),size); break;
			                      	case miUTF8   : c16[0] = (int8u)d07[i]; Write(reinterpret_cast<char*>(c16),size); break;
			                      	case miUTF16  : c17[0] = d07[i]; Write(reinterpret_cast<char*>(c17),size); break;
			                      	case miUTF32  : c18[0] = (int32u)d07[i]; Write(reinterpret_cast<char*>(c18),size); break;
			                      }; break;
			case mxINT32_CLASS  : switch(DataType){
			                      	case miINT8   : c01[0] = (int8s)d08[i]; Write(reinterpret_cast<char*>(c01),size); break;
			                      	case miUINT8  : c02[0] = (int8u)d08[i]; Write(reinterpret_cast<char*>(c02),size); break;
			                      	case miINT16  : c03[0] = (int16s)d08[i]; Write(reinterpret_cast<char*>(c03),size); break;
			                      	case miUINT16 : c04[0] = (int16u)d08[i]; Write(reinterpret_cast<char*>(c04),size); break;
			                      	case miINT32  : c05[0] = d08[i]; Write(reinterpret_cast<char*>(c05),size); break;
			                      	case miUINT32 : c06[0] = (int32u)d08[i]; Write(reinterpret_cast<char*>(c06),size); break;
			                      	case miSINGLE : c07[0] = (float)d08[i]; Write(reinterpret_cast<char*>(c07),size); break;
			                      	case miDOUBLE : c09[0] = (double)d08[i]; Write(reinterpret_cast<char*>(c09),size); break;
			                      	case miINT64  : c12[0] = (int64s)d08[i]; Write(reinterpret_cast<char*>(c12),size); break;
			                      	case miUINT64 : c13[0] = (int64u)d08[i]; Write(reinterpret_cast<char*>(c13),size); break;
			                      	case miUTF8   : c16[0] = (int8u)d08[i]; Write(reinterpret_cast<char*>(c16),size); break;
			                      	case miUTF16  : c17[0] = (int16u)d08[i]; Write(reinterpret_cast<char*>(c17),size); break;
			                      	case miUTF32  : c18[0] = (int32u)d08[i]; Write(reinterpret_cast<char*>(c18),size); break;
			                      }; break;
			case mxUINT32_CLASS : switch(DataType){
			                      	case miINT8   : c01[0] = (int8s)d09[i]; Write(reinterpret_cast<char*>(c01),size); break;
			                      	case miUINT8  : c02[0] = (int8u)d09[i]; Write(reinterpret_cast<char*>(c02),size); break;
			                      	case miINT16  : c03[0] = (int16s)d09[i]; Write(reinterpret_cast<char*>(c03),size); break;
			                      	case miUINT16 : c04[0] = (int16u)d09[i]; Write(reinterpret_cast<char*>(c04),size); break;
			                      	case miINT32  : c05[0] = (int32s)d09[i]; Write(reinterpret_cast<char*>(c05),size); break;
			                      	case miUINT32 : c06[0] = d09[i]; Write(reinterpret_cast<char*>(c06),size); break;
			                      	case miSINGLE : c07[0] = (float)d09[i]; Write(reinterpret_cast<char*>(c07),size); break;
			                      	case miDOUBLE : c09[0] = (double)d09[i]; Write(reinterpret_cast<char*>(c09),size); break;
			                      	case miINT64  : c12[0] = (int64s)d09[i]; Write(reinterpret_cast<char*>(c12),size); break;
			                      	case miUINT64 : c13[0] = (int64u)d09[i]; Write(reinterpret_cast<char*>(c13),size); break;
			                      	case miUTF8   : c16[0] = (int8u)d09[i]; Write(reinterpret_cast<char*>(c16),size); break;
			                      	case miUTF16  : c17[0] = (int16u)d09[i]; Write(reinterpret_cast<char*>(c17),size); break;
			                      	case miUTF32  : c18[0] = d09[i]; Write(reinterpret_cast<char*>(c18),size); break;
			                      }; break;
			case mxINT64_CLASS  : switch(DataType){
			                      	case miINT8   : c01[0] = (int8s)d10[i]; Write(reinterpret_cast<char*>(c01),size); break;
			                      	case miUINT8  : c02[0] = (int8u)d10[i]; Write(reinterpret_cast<char*>(c02),size); break;
			                      	case miINT16  : c03[0] = (int16s)d10[i]; Write(reinterpret_cast<char*>(c03),size); break;
			                      	case miUINT16 : c04[0] = (int16u)d10[i]; Write(reinterpret_cast<char*>(c04),size); break;
			                      	case miINT32  : c05[0] = (int32s)d10[i]; Write(reinterpret_cast<char*>(c05),size); break;
			                      	case miUINT32 : c06[0] = (int32u)d10[i]; Write(reinterpret_cast<char*>(c06),size); break;
			                      	case miSINGLE : c07[0] = (float)d10[i]; Write(reinterpret_cast<char*>(c07),size); break;
			                      	case miDOUBLE : c09[0] = (double)d10[i]; Write(reinterpret_cast<char*>(c09),size); break;
			                      	case miINT64  : c12[0] = d10[i]; Write(reinterpret_cast<char*>(c12),size); break;
			                      	case miUINT64 : c13[0] = (int64u)d10[i]; Write(reinterpret_cast<char*>(c13),size); break;
			                      	case miUTF8   : c16[0] = (int8u)d10[i]; Write(reinterpret_cast<char*>(c16),size); break;
			                      	case miUTF16  : c17[0] = (int16u)d10[i]; Write(reinterpret_cast<char*>(c17),size); break;
			                      	case miUTF32  : c18[0] = (int32u)d10[i]; Write(reinterpret_cast<char*>(c18),size); break;
			                      }; break;
			case mxUINT64_CLASS : switch(DataType){
			                      	case miINT8   : c01[0] = (int8s)d11[i]; Write(reinterpret_cast<char*>(c01),size); break;
			                      	case miUINT8  : c02[0] = (int8u)d11[i]; Write(reinterpret_cast<char*>(c02),size); break;
			                      	case miINT16  : c03[0] = (int16s)d11[i]; Write(reinterpret_cast<char*>(c03),size); break;
			                      	case miUINT16 : c04[0] = (int16u)d11[i]; Write(reinterpret_cast<char*>(c04),size); break;
			                      	case miINT32  : c05[0] = (int32s)d11[i]; Write(reinterpret_cast<char*>(c05),size); break;
			                      	case miUINT32 : c06[0] = (int32u)d11[i]; Write(reinterpret_cast<char*>(c06),size); break;
			                      	case miSINGLE : c07[0] = (float)d11[i]; Write(reinterpret_cast<char*>(c07),size); break;
			                      	case miDOUBLE : c09[0] = (double)d11[i]; Write(reinterpret_cast<char*>(c09),size); break;
			                      	case miINT64  : c12[0] = (int64s)d11[i]; Write(reinterpret_cast<char*>(c12),size); break;
			                      	case miUINT64 : c13[0] = d11[i]; Write(reinterpret_cast<char*>(c13),size); break;
			                      	case miUTF8   : c16[0] = (int8u)d11[i]; Write(reinterpret_cast<char*>(c16),size); break;
			                      	case miUTF16  : c17[0] = (int16u)d11[i]; Write(reinterpret_cast<char*>(c17),size); break;
			                      	case miUTF32  : c18[0] = (int32u)d11[i]; Write(reinterpret_cast<char*>(c18),size); break;
			                      }; break;
		}
	}
	//padding
	int padding;
	if(!small){
		padding=Bytes/8;
		if((Bytes%8)!=0)padding++;
		padding=padding*8;
		padding=padding-Bytes;
	}
	else{
		padding=Bytes/4;
		if((Bytes%4)!=0)padding++;
		padding=padding*4;
		padding=padding-Bytes;
	}
	char a1[1];
	a1[0] = 0;
	for(int i=0; i<padding; i++)
		Write(a1,1);
	return !Fail();
}

//---------------------------------------------------------------------------

bool MATFile::ReadMatrixHeader(mxArray* mxA){
	mxA->SetEmpty();
	bool failed=false;
	//Array Flags
	miMatrix DataType;
	int Bytes, MaxLength;
	bool small;
	//flag + class
	unsigned int *Data1;
	ReadTag(DataType, Bytes, small);
	Data1 = new unsigned int [Bytes/4];
	if(DataType==miUINT32)ReadDataElements(mxUINT32_CLASS,DataType, Data1, Bytes, 4, small);
	else failed=true;
	mxClassID ArrayType = (mxClassID)(Data1[0]%256);
	mxA->SetArrayType(ArrayType);
	if(Data1[0]==mxVOID_CLASS)mxA->SetNzmax(Data1[1]);
	Data1[0]=Data1[0]/512;
	if((Data1[0]%2)==1)mxA->SetLogical(true);
	else mxA->SetLogical(false);
	Data1[0]=Data1[0]/2;
	if((Data1[0]%2)==1)mxA->SetGlobal(true);
	else mxA->SetGlobal(false);
	Data1[0]=Data1[0]/2;
	if((Data1[0]%2)==1)mxA->SetComplex(true);
	else mxA->SetComplex(false);
	delete[] Data1;
	//Dimensions Array
	int *Data2;
	ReadTag(DataType, Bytes, small);
	Data2 = new int [Bytes/4];
	if(DataType==miINT32)ReadDataElements(mxINT32_CLASS,DataType, Data2, Bytes, 4, small);
	else failed=true;
	mxA->SetDim(Data2,Bytes/4);
	delete[] Data2;
	//Array Name
	char *Data3;
	ReadTag(DataType, Bytes, small);
	if(Bytes>0){
		Data3 = new char [Bytes];
		if(DataType==miINT8)ReadDataElements(mxCHAR_CLASS,DataType, Data3, Bytes, 1, small);
		else failed=true;
		mxA->SetName(Data3,Bytes);
		delete[] Data3;
	}
	
	if(ArrayType==mxSTRUCT_CLASS){
	  
		ReadTag(DataType, Bytes, small);
		if(DataType==miINT8){
		  Data3 = new char [Bytes];
		  ReadDataElements(mxCHAR_CLASS,DataType, Data3, Bytes, 1, small);
		  mxA->SetCname(Data3,Bytes);
		  delete[] Data3;
		  ReadTag(DataType, Bytes, small);
		}
		if(DataType==miINT32){
			Data2 = new int [Bytes/4];
			ReadDataElements(mxINT32_CLASS,DataType, Data2, Bytes, 4, small);
			MaxLength=Data2[0];
			delete[] Data2;
			ReadTag(DataType, Bytes, small);
			mxA->SetNf(Bytes,MaxLength);
			int Npos=0;
			for(int i=0; i<mxA->GetNf(); i++){
				Data3 = mxA->GetFieldName(i);
				Npos = GetPos();
				Read(Data3,MaxLength);
				SetPos(Npos+MaxLength);
			}
			int padding =0;
			char *temp = new char [8];
			if(!small){
			  padding= (MaxLength*(mxA->GetNf()))/8;
			  if( (MaxLength*(mxA->GetNf()))%8!=0 )
			    padding++;
			  padding=padding*8;
			  padding = padding-MaxLength*(mxA->GetNf());
			  			  
			}
			else{
			  padding=(MaxLength*(mxA->GetNf()))/4;
			  if( (MaxLength*(mxA->GetNf()))%4!=0 )padding++;
			  padding=padding*4;
			  padding=padding-(MaxLength*(mxA->GetNf()));
			}

			int size=1;
			for(int i=0; i<padding; i++)
			    Read(temp,size);
			delete[] temp;
		}
		else failed=true;
	}
	
	mxA->SetDataOffset(GetPos());
	return failed;
}

//---------------------------------------------------------------------------

bool MATFile::WriteMatrixHeader(mxArray* mxA){
	bool failed=false;
	bool small;
	mxClassID ArrayType = mxA->GetArrayType();
	//flag + class
	WriteTag(miUINT32, 8, small);
	unsigned int *Data1;
	Data1 = new unsigned int [2];
	Data1[0] = 0;
	if(mxA->GetLogical())Data1[0] = Data1[0]+2;
	if(mxA->GetGlobal())Data1[0] = Data1[0]+4;
	if(mxA->GetComplex())Data1[0] = Data1[0]+8;
	Data1[0] = Data1[0]*256+ArrayType;
	if(ArrayType==mxVOID_CLASS)Data1[1] = mxA->GetNzmax();
	else Data1[1]=0;
	WriteDataElements(mxUINT32_CLASS,miUINT32, Data1, 8, 4, small);
	delete[] Data1;
	//Dimensions Array
	WriteTag(miINT32, 4*mxA->GetNd(), small);
	int *Data2;
	Data2 = new int [mxA->GetNd()];
	for(int i=0; i<mxA->GetNd(); i++)
		Data2[i] = mxA->GetDim(i);
	WriteDataElements(mxINT32_CLASS,miINT32, Data2, 4*mxA->GetNd(), 4, small);
	delete[] Data2;
	//Array Name
	WriteTag(miINT8, mxA->GetNc(), small);
	char *Data3;
	if(mxA->GetNc()>0){
		Data3 = mxA->GetName();
		WriteDataElements(mxCHAR_CLASS,miINT8, Data3, mxA->GetNc(), 1, small);
	}
	if(ArrayType==mxSTRUCT_CLASS){
	  if(mxA->GetNe()>0){
	    //Class name
	    WriteTag(miINT8, mxA->GetNe(), small);
	    char *Data3;
	    if(mxA->GetNe()>0){
	      Data3 = mxA->GetCname();
	      WriteDataElements(mxCHAR_CLASS,miINT8, Data3, mxA->GetNe(), 1, small);
	    }
	  }
	  //Field Name Length
	  WriteTag(miINT32, 4, small);
	  Data2 = new int [1];
	  Data2[0] = 32;
	  WriteDataElements(mxINT32_CLASS,miINT32, Data2, 4, 4, small);
	  delete[] Data2;
	  //Field Names
	  WriteTag(miINT8, mxA->GetNf()*32, small);
	  for(int i=0; i<mxA->GetNf(); i++){
	    Data3 = mxA->GetFieldName(i);
	    WriteDataElements(mxCHAR_CLASS,miINT8, Data3, 32, 1, small);
	  }
	}
	return failed;
}

//---------------------------------------------------------------------------

bool MATFile::ReadMatrixDataElement(mxArray* mxA){
	int Bytes, Elements;
	bool small, failed=false;
	miMatrix DataType;
	mxClassID ArrayType = mxA->GetArrayType();
	SetPos(mxA->GetDataOffset());
	//numeric/char array
	void * Data4;
	if((ArrayType>mxSTRUCT_CLASS)&&(ArrayType<mxFUNCTION_CLASS)){
		//sparse array
		if(ArrayType==mxVOID_CLASS){
			//Ir
			ReadTag(DataType, Bytes, small);
			Data4 = mxA->MakeDataP(DataType,Bytes,Elements,2);
			if(Bytes>0)ReadDataElements(mxINT32_CLASS,DataType, Data4, Bytes, Bytes/Elements, small);
			//Jr
			ReadTag(DataType, Bytes, small);
			Data4 = mxA->MakeDataP(DataType,Bytes,Elements,3);
			if(Bytes>0)ReadDataElements(mxINT32_CLASS, DataType, Data4, Bytes, Bytes/Elements, small);
		}
		//Pr or Data
		ReadTag(DataType, Bytes, small);
		Data4 = mxA->MakeDataP(DataType,Bytes,Elements,0);
		if(Bytes>0)ReadDataElements(ArrayType,DataType, Data4, Bytes, Bytes/Elements, small);
		//Pi
		if(mxA->GetComplex()){
			ReadTag(DataType, Bytes, small);
			Data4 = mxA->MakeDataP(DataType,Bytes,Elements,1);
			if(Bytes>0)ReadDataElements(ArrayType,DataType, Data4, Bytes, Bytes/Elements, small);
		}
	}
	//cell array / struct array / object array
	if((ArrayType==mxCELL_CLASS)||(ArrayType==mxSTRUCT_CLASS)){
		//mxA0
		mxArray* mxA0;
		Data4 = mxA->MakeDataP(DataType,Bytes,Elements,0);
		for(int i=0; i<Elements; i++){
			mxA0 = mxA->MakeMatrixElm(i);
			bool ret = ReadTag(DataType, Bytes, small);
			if(Bytes>0){
			  ReadMatrixHeader(mxA0);
			  ReadMatrixDataElement(mxA0);
			}
			if(EoF()){failed=true;i=Elements;}
		}
	}
	return failed;
}

//---------------------------------------------------------------------------

bool MATFile::WriteMatrixDataElement(mxArray* mxA){
	int Bytes, size, Npos;
	bool small, failed=false;
	miMatrix DataType;
	mxClassID ArrayType = mxA->GetArrayType();
	void * Data4;
	//numeric array
	if((ArrayType>mxSTRUCT_CLASS)&&(ArrayType<mxFUNCTION_CLASS)){
		switch(ArrayType){
			case mxLOGICAL_CLASS: DataType=miUINT8; size=1; break;
			case mxCHAR_CLASS   : DataType=miUINT16; size=2; break;
			case mxVOID_CLASS   :
			case mxDOUBLE_CLASS : DataType=miDOUBLE; size=8; break;
			case mxSINGLE_CLASS : DataType=miSINGLE; size=4; break;
			case mxINT8_CLASS   : DataType=miINT8; size=1; break;
			case mxUINT8_CLASS  : DataType=miUINT8; size=1; break;
			case mxINT16_CLASS  : DataType=miINT16; size=2; break;
			case mxUINT16_CLASS : DataType=miUINT16; size=2; break;
			case mxINT32_CLASS  : DataType=miINT32; size=4; break;
			case mxUINT32_CLASS : DataType=miUINT32; size=4; break;
			case mxINT64_CLASS  : DataType=miINT64; size=8; break;
			case mxUINT64_CLASS : DataType=miUINT64; size=8; break;
		}
		if(ArrayType==mxVOID_CLASS){
			//Ir
			WriteTag(miINT32, 4*mxA->GetNdata(), small);
			Data4 = mxA->GetDataP(2);
			if(mxA->GetNdata()>0)WriteDataElements(ArrayType,miINT32, Data4, 4*mxA->GetNdata(), 4, small);
			//Jr
			WriteTag(miINT32, 4*mxA->GetNdata(), small);
			Data4 = mxA->GetDataP(3);
			if(mxA->GetNdata()>0)WriteDataElements(ArrayType,miINT32, Data4, 4*mxA->GetNdata(), 4, small);
		}
		//Pr or Data
		WriteTag(DataType, size*mxA->GetNdata(), small);
		Data4 = mxA->GetDataP(0);
		if(mxA->GetNdata()>0)WriteDataElements(ArrayType,DataType, Data4, size*mxA->GetNdata(), size, small);
		if(mxA->GetComplex()){
			//Pi
			WriteTag(DataType, size*mxA->GetNdata(), small);
			Data4 = mxA->GetDataP(1);
			if(mxA->GetNdata()>0)WriteDataElements(ArrayType,DataType, Data4, size*mxA->GetNdata(), size, small);
		}
	}
	//cell array / struct array / object array
	if((ArrayType==mxCELL_CLASS)||(ArrayType==mxSTRUCT_CLASS)){
		//Data
		mxArray* mxA0;
		for(int i=0; i<mxA->GetNdata(); i++){
			small=false;
			WriteTag(miMATRIX, 0, small);
			Npos = GetPos();
//			cout << "Npos(1) = " << Npos << endl;
			mxA0 = mxA->GetMatrixElm(i);
			if((mxA0!=NULL)&&(!mxA0->GetEmpty())){
				WriteMatrixHeader(mxA0);
				WriteMatrixDataElement(mxA0);
			}
			small=true;
//			cout << "Npos(2) = " << Npos << endl;
			WriteTag(miMATRIX, Npos, small);
		}
	}
	return failed;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

mxArray::mxArray(){
	ArrayType = mxUNKNOWN_CLASS;
	Ndata=0;
	nzmax=0;
	nd=0;
	nc=0;
	nf=0;
	ne=0;
	DataOffset=0;
	isComplex=false;
	isLogical=false;
	isGlobal=false;
	Empty=true;
}

//---------------------------------------------------------------------------

mxArray::~mxArray(){
	if(Ndata>0){
		if((ArrayType>mxSTRUCT_CLASS)&&(ArrayType<mxFUNCTION_CLASS)){
			if((ArrayType==mxDOUBLE_CLASS)||(ArrayType==mxVOID_CLASS)){
				delete[] _Pr;
				if(isComplex)delete[] _Pi;
			}
			else delete[] Data;
		}
		else{
			for(int i=0; i<Ndata; i++)
				delete mxA0[i];
			delete[] mxA0;
		}
	}
	if(nzmax>0){
		delete[] _Ir;
		delete[] _Jc;
	}
	if(nd>0)delete[] Dim;
	if(nc>0)delete[] Name;
	if(nf>0){
		for(int i=0; i<nf; i++)
			delete[] FieldNames[i];
		delete[] FieldNames;
	}
	if(ne>0)delete[] Cname;
}

//---------------------------------------------------------------------------

void mxArray::SetComplex(bool isC){
	isComplex=isC;
}

//---------------------------------------------------------------------------

void mxArray::SetLogical(bool isL){
	isLogical=isL;
}

//---------------------------------------------------------------------------

void mxArray::SetGlobal(bool isG){
	isGlobal=isG;
}

//---------------------------------------------------------------------------

void mxArray::SetArrayType(mxClassID ArrayType0){
	ArrayType = ArrayType0;
}

//---------------------------------------------------------------------------

void mxArray::SetName(char* Name1, int nc1){
	if(nc1>0){
		nc = nc1;
		Name = new char[nc];
		for(int i=0; i<nc; i++)
			Name[i] = Name1[i];
		Name[nc]=(char)0;
	}
}

//---------------------------------------------------------------------------

void mxArray::SetCname(char* Name1, int ne1){
	if(ne1>0){
		ne = ne1;
		Cname = new char[ne];
		for(int i=0; i<ne; i++)
			Cname[i] = Name1[i];
		Cname[ne] = (char)0;
	}
}

//---------------------------------------------------------------------------

void mxArray::SetNzmax(int nzmax0){
	nzmax = nzmax0;
}

//---------------------------------------------------------------------------

void mxArray::SetDim(int* dim1, int nd1){
	if(nd1>0){
		nd = nd1;
		Dim = new int[nd];
		for(int i=0; i<nd; i++)
			Dim[i] = dim1[i];
	}
}

//---------------------------------------------------------------------------

void* mxArray::MakeDataP(miMatrix DataType, int Bytes, int& Elements, int mode){
	int size;
	void *Data0;
	switch(DataType){
		case miINT8   : size=1; break;
		case miUINT8  : size=1; break;
		case miINT16  : size=2; break;
		case miUINT16 : size=2; break;
		case miINT32  : size=4; break;
		case miUINT32 : size=4; break;
		case miSINGLE : size=4; break;
		case miDOUBLE : size=8; break;
		case miINT64  : size=8; break;
		case miUINT64 : size=8; break;
		case miUTF8   : size=1; break;
		case miUTF16  : size=2; break;
		case miUTF32  : size=4; break;
	}
	if(ArrayType<mxLOGICAL_CLASS){
		if(ArrayType==mxCELL_CLASS)Elements=1;
		else Elements=nf;
		for(int i=0; i<nd; i++)
			Elements=Dim[i]*Elements;
	}
	else Elements=Bytes/size;
	if(mode<2)Ndata=Elements;
	if(Elements>0){		switch(ArrayType){
			case mxCELL_CLASS   :
			case mxSTRUCT_CLASS : mxA0 = new mxArray* [Elements]; for(int i=0; i<Elements; i++)mxA0[i] = NULL; break;
			case mxLOGICAL_CLASS: Data = new mxLogical[Elements]; break;
			case mxCHAR_CLASS   : Data = new char[Elements]; break;
			case mxVOID_CLASS   : switch(mode){
			                      	case 0 : _Pr = new double[Elements]; break;
			                      	case 1 : _Pi = new double[Elements]; break;
			                      	case 2 : _Ir = new int[Elements]; break;
			                      	case 3 : _Jc = new int[Elements]; break;
			                      };break;
			case mxDOUBLE_CLASS : if(mode==0)_Pr = new double[Elements];
			                      	if(mode==1)_Pi = new double[Elements];
			                      break;
			case mxSINGLE_CLASS : Data = new float[Elements]; break;
			case mxINT8_CLASS   : Data = new int8s[Elements]; break;
			case mxUINT8_CLASS  : Data = new int8u[Elements]; break;
			case mxINT16_CLASS  : Data = new int16s[Elements]; break;
			case mxUINT16_CLASS : Data = new int16u[Elements]; break;
			case mxINT32_CLASS  : Data = new int32s[Elements]; break;
			case mxUINT32_CLASS : Data = new int32u[Elements]; break;
			case mxINT64_CLASS  : Data = new int64s[Elements]; break;
			case mxUINT64_CLASS : Data = new int64u[Elements]; break;
		}
	}
	if((ArrayType==mxDOUBLE_CLASS)||(ArrayType==mxVOID_CLASS))
		switch(mode){
			case 0: Data0 = _Pr; break;
			case 1: Data0 = _Pi; break;
			case 2: Data0 = _Ir; break;
			case 3: Data0 = _Jc; break;
		}
	else Data0 = Data;
	return Data0;
}

//---------------------------------------------------------------------------

void* mxArray::GetDataP(int mode){
	void* Data0;
	if((ArrayType==mxDOUBLE_CLASS)||(ArrayType==mxVOID_CLASS))
		switch(mode){
			case 0: Data0 = _Pr; break;
			case 1: Data0 = _Pi; break;
			case 2: Data0 = _Ir; break;
			case 3: Data0 = _Jc; break;
		}
	else Data0 = Data;
	return Data0;
}

//---------------------------------------------------------------------------

void mxArray::SetDataP(int mode, void* Data0){
	if((ArrayType==mxDOUBLE_CLASS)||(ArrayType==mxVOID_CLASS))
		switch(mode){
			case 0: _Pr = (double *)Data0; break;
			case 1: _Pi = (double *)Data0; break;
			case 2: _Ir = (int *)Data0; break;
			case 3: _Jc = (int *)Data0; break;
		}
	else Data = Data0;
}

//---------------------------------------------------------------------------

bool mxArray::GetComplex(){
	return isComplex;
}

//---------------------------------------------------------------------------

bool mxArray::GetLogical(){
	return isLogical;
}

//---------------------------------------------------------------------------

bool mxArray::GetGlobal(){
	return isGlobal;
}

//---------------------------------------------------------------------------

mxClassID mxArray::GetArrayType(){
	return ArrayType;
}

//---------------------------------------------------------------------------

int mxArray::GetNc(){
	return nc;
}

//---------------------------------------------------------------------------

int mxArray::GetNd(){
	return nd;
}

//---------------------------------------------------------------------------

int mxArray::GetNf(){
	return nf;
}

//---------------------------------------------------------------------------

int mxArray::GetNe(){
	return ne;
}

//---------------------------------------------------------------------------

int mxArray::GetNzmax(){
	return nzmax;
}

//---------------------------------------------------------------------------

int mxArray::GetNdata(){
	return Ndata;
}

//---------------------------------------------------------------------------

char* mxArray::GetName(){
	return Name;
}

//---------------------------------------------------------------------------

char* mxArray::GetCname(){
	return Cname;
}

//---------------------------------------------------------------------------

int mxArray::GetDim(int index){
	return Dim[index];
}

//---------------------------------------------------------------------------

void mxArray::SetNf(int Bytes,int MaxLength){
	if(Bytes>0){
		nf = Bytes/MaxLength;
		FieldNames = new char * [nf];
		for(int i=0; i<nf; i++)
			FieldNames[i] = new char [MaxLength];
	}
}

//---------------------------------------------------------------------------

void mxArray::SetFieldNames(int nf1, int MaxLength, const char** FN1){
	nf = nf1;
	bool NL;
	FieldNames = new char * [nf];
	for(int i=0; i<nf; i++){
		FieldNames[i] = new char [MaxLength];
		NL=false;
		for(int j=0; j<MaxLength; j++){
			if(FN1[i][j]==(char)0)NL=true;
			if(!NL)FieldNames[i][j] = FN1[i][j];
			else FieldNames[i][j] = 0;
		}
	}
}
//---------------------------------------------------------------------------

int mxArray::AddField(const char *fieldname, int MaxLength){
  bool NL;
  nf+=1;
  //char **OldFieldNames = FieldNames;
  char **OldFieldNames = new char * [nf];
  for(int i=0; i<(nf-1); i++)
    OldFieldNames[i] = FieldNames[i];
  delete[] FieldNames;
  FieldNames = new char * [nf];
  	
  //add name
  for(int i=0; i<nf; i++){
    FieldNames[i] = new char [MaxLength];
    NL=false;
    if(i!=(nf-1)){
      for(int j=0; j<MaxLength; j++){
	if(OldFieldNames[i][j]==(char)0)NL=true;
	if(!NL){
	  FieldNames[i][j] = OldFieldNames[i][j];
	}
	else FieldNames[i][j] = 0;
      }
      delete[] OldFieldNames[i];
    }
    else
      for(int j=0; j<(strlen(fieldname)+1); j++){
	if(fieldname[j]==(char)0)NL=true;
	if(!NL)FieldNames[i][j] = fieldname[j];
	else FieldNames[i][j] = 0;
      }
    
  }
  delete[] OldFieldNames;
  int OldElements=Ndata;
  int NewElements=nf;
  

  for(int i=0; i<nd; i++)
    NewElements=Dim[i]*NewElements;
  Ndata=NewElements;
  mxArray** mxA0copy = new mxArray* [OldElements];
  memcpy(mxA0copy,mxA0,sizeof( mxArray *)*(OldElements));
  delete[] mxA0;
  mxA0 = new mxArray* [NewElements]; 
  for(int i=0; i<OldElements; i++)mxA0[i] = mxA0copy[i];
  for(int i=OldElements; i<NewElements; i++)mxA0[i] = NULL;
  

  delete[] mxA0copy;
  return nf;

}

//---------------------------------------------------------------------------

int mxGetFieldNumber(const mxArray *pm, const char *fieldname){
  return ((mxArray *)pm)->GetFieldNumber((char *)fieldname);
}

//---------------------------------------------------------------------------

int mxArray::GetFieldNumber(char *fieldname){
  int retval = -1;
  for(int i=0;i<nf;i++){
    if(strlen(fieldname)==strlen(FieldNames[i])){
      if(!strcmp(fieldname,FieldNames[i])){
	retval = i;
      }
    }
  }
  return retval;
}

//---------------------------------------------------------------------------

void mxRemoveField(mxArray *pm, int fieldnumber){
  pm->RemoveField(fieldnumber,32);
}

//---------------------------------------------------------------------------

void mxArray::RemoveField( int fieldnumber, int MaxLength){
  bool NL;
  nf-=1;
  //char **OldFieldNames = FieldNames;
  char **OldFieldNames = new char * [nf+1];
  for(int i=0; i<(nf+1); i++)
    OldFieldNames[i] = FieldNames[i];
  delete[] FieldNames;
  FieldNames = new char * [nf];
  	
  //add name
  int count=0;
  for(int i=0; i<(nf+1); i++){
    if( i!=fieldnumber ){
      FieldNames[count] = new char [MaxLength];
      NL=false;
      for(int j=0; j<MaxLength; j++){
	if(OldFieldNames[i][j]==(char)0)NL=true;
	if(!NL){
	  FieldNames[count][j] = OldFieldNames[i][j];
	}
	else FieldNames[count][j] = 0;
      }
      delete[] OldFieldNames[i];
      count++;
    }
  }
  delete[] OldFieldNames;
  int OldElements=Ndata;
  int NewElements=nf;
  

  for(int i=0; i<nd; i++)
    NewElements=Dim[i]*NewElements;
  Ndata=NewElements;
  mxArray** mxA0copy = new mxArray* [NewElements+1];
  memcpy(mxA0copy,mxA0,sizeof( mxArray *)*(NewElements+1));
  delete[] mxA0;
  mxA0 = new mxArray* [NewElements];
  count = 0;
  for(int i=0; i<(NewElements+1); i++){
    if( i!=fieldnumber ){
      mxA0[count] = mxA0copy[i];
      count++;
    }
  }

  delete[] mxA0copy;

}

//---------------------------------------------------------------------------

char* mxArray::GetFieldName(int index){
	return FieldNames[index];
}

//---------------------------------------------------------------------------

void mxArray::SetDataOffset(int DataOffset0){
	DataOffset = DataOffset0;
}

//---------------------------------------------------------------------------

int mxArray::GetDataOffset(){
	return DataOffset;
}
//---------------------------------------------------------------------------

int mxIsStruct(const mxArray *pm){
  return ((mxArray *)pm)->GetArrayType()==mxSTRUCT_CLASS;
}
//---------------------------------------------------------------------------

int mxIsDouble(const mxArray *pm){
  return ((mxArray *)pm)->GetArrayType()==mxDOUBLE_CLASS;
}

//---------------------------------------------------------------------------

bool mxIsEmpty(const mxArray *pm){
  return ((mxArray *)pm)->GetEmpty();
}

//---------------------------------------------------------------------------

int mxIsChar(const mxArray *pm){
  return ((mxArray *)pm)->GetArrayType()==mxCHAR_CLASS;
}

//---------------------------------------------------------------------------

int mxIsUint8(const mxArray *pm){
  return ((mxArray *)pm)->GetArrayType()==mxUINT8_CLASS;
}

//---------------------------------------------------------------------------

int mxGetNumberOfFields(const mxArray *pm){
  return ((mxArray *)pm)->GetNf();
}

//---------------------------------------------------------------------------

int mxGetNumberOfDimensions(const mxArray *pm){
  return ((mxArray *)pm)->GetNd();
}


//---------------------------------------------------------------------------

mxArray* mxArray::GetMatrixElm(int index){
	return mxA0[index];
}

//---------------------------------------------------------------------------

void mxArray::SetMatrixElm(mxArray* mxA1, int index){
	mxA0[index] = mxA1;
}

//---------------------------------------------------------------------------

mxArray* mxArray::MakeMatrixElm(int index){
	mxA0[index] = new mxArray;
	return mxA0[index];
}

//---------------------------------------------------------------------------

void mxArray::PrintArrayTypes(){
	for(int i=0; i<Ndata; i++)
		if(!mxA0[i]->GetEmpty())cout << "ArrayType(" << i << ") = " << mxA0[i]->GetArrayType() << endl;
		else cout << "ArrayType(" << i << ") = empty" << endl;
}

//---------------------------------------------------------------------------

void mxArray::SetEmpty(){
	Empty=false;
}

//---------------------------------------------------------------------------

bool mxArray::GetEmpty(){
	return Empty;
}

//---------------------------------------------------------------------------
