//---------------------------------------------------------------------------
#ifndef Mat_ioH
#define Mat_ioH
//---------------------------------------------------------------------------
#include <string>
#include <fstream>
#include "Mat_io_globals.h"
#include "zlib.h"
//---------------------------------------------------------------------------
#define _ALLOW_COMPRESSION 0

typedef size_t mwSize;         /* unsigned pointer-width integer */
typedef size_t mwIndex;        /* unsigned pointer-width integer */
//typedef ptrdiff_t mwSignedIndex;  /* a signed pointer-width integer */

//#define CHUNK 4194304
#define CHUNK 10194304
#define MAX_STORED 4096

typedef char int8s;
typedef unsigned char int8u;
typedef short int int16s;
typedef unsigned short int int16u;
typedef int int32s;
typedef unsigned int int32u;
typedef long long int64s;
typedef unsigned long long int64u;

typedef bool mxLogical;

enum mxClassID {
	mxUNKNOWN_CLASS = 0,
	mxCELL_CLASS,
	mxSTRUCT_CLASS,
	mxLOGICAL_CLASS,
	mxCHAR_CLASS,
	mxVOID_CLASS,
	mxDOUBLE_CLASS,
	mxSINGLE_CLASS,
	mxINT8_CLASS,
	mxUINT8_CLASS,
	mxINT16_CLASS,
	mxUINT16_CLASS,
	mxINT32_CLASS,
	mxUINT32_CLASS,
	mxINT64_CLASS,
	mxUINT64_CLASS,
	mxFUNCTION_CLASS,
	mxOPAQUE_CLASS,
	mxOBJECT_CLASS
};

//#define mxINDEX_CLASS ((sizeof(size_t)==4) ? mxUINT32_CLASS : mxUINT64_CLASS)

enum mxComplexity {mxREAL=0, mxCOMPLEX};

enum miMatrix {
	miINT8 = 1,
	miUINT8,
	miINT16,
	miUINT16,
	miINT32,
	miUINT32,
	miSINGLE,
	miDOUBLE = 9,
	miINT64 = 12,
	miUINT64,
	miMATRIX,
	miCOMPRESSED,
	miUTF8,
	miUTF16,
	miUTF32
};


class mxArray {
private:
	mxClassID ArrayType;
	void* Data;
	mxArray** mxA0;
	double *_Pr, *_Pi;
	int *_Ir, *_Jc;
	int Ndata, nzmax;
	bool isLogical, isComplex, isGlobal;
	int* Dim;
	char* Name;
	char* Cname;
	char ** FieldNames;
	/*
	  nd - number of dimensions
	  nc - number of chars in mat filename
	  nf - number of fields
	  ne - number of ??
	*/
	int nd, nc, nf, ne;
	int DataOffset;
	bool Empty;
public:
	mxArray();
	~mxArray();
	bool GetComplex();
	bool GetEmpty();
	bool GetGlobal();
	bool GetLogical();
	mxClassID GetArrayType();
	int GetDim(int);
	char* GetName();
	char* GetCname();
	char* GetFieldName(int);
	void* GetDataP(int);
	int GetDataOffset();
	mxArray* GetMatrixElm(int);
	int GetNf();
	int GetNzmax();
	int GetNc();
	int GetNd();
	int GetNe();
	int GetNdata();
	int GetFieldNumber(char *);
	void SetEmpty();
	void SetLogical(bool);
	void SetGlobal(bool);
	void SetComplex(bool);
	void SetArrayType(mxClassID);
	void SetName(char*,int);
	void SetCname(char*,int);
	void SetNf(int, int);
	void SetDim(int*,int);
	void SetDataP(int,void*);
	void SetNzmax(int);
	void SetDataOffset(int);
	void SetFieldNames(int, int, const char**);
	int AddField(const char *, int);
	void RemoveField( int , int);
	void SetMatrixElm(mxArray*,int);
	void* MakeDataP(miMatrix, int, int&, int);
	mxArray* MakeMatrixElm(int);
	void PrintArrayTypes();
};


class MATFile {
private:
	fstream File1;
	z_stream strm;
	int strm_flush, strm_error, strm_begin, strm_idx, strm_pos;
	unsigned char *strm_in, *strm_out;
	bool miRead, RWmode;
	int ReadStrm;
	int mode, Nindex, Npos;
	int write_pos, **store_pos;
public:
	MATFile(const char*, int);
	~MATFile();
	bool EoF();
	bool Fail();
	bool IsOpen();
	void Reset();
	int GetMode();
	int GetPos();
	int GetNindex();
	void SetRWmode(bool);
	void SetPos(int);
	void SetNindex(int);
	bool CheckHeader();
	void WriteHeader();
	bool ReadTag(miMatrix&, int&, bool&);
	bool ReadDataElements(mxClassID, miMatrix, void*, int, int, bool);
	bool ReadMatrixHeader(mxArray*);
	bool ReadMatrixDataElement(mxArray*);
	bool WriteTag(miMatrix, int, bool&);
	bool WriteDataElements(mxClassID, miMatrix, void*, int, int, bool);
	bool WriteMatrixHeader(mxArray*);
	bool WriteMatrixDataElement(mxArray*);
	bool InitStream(int);
	void EndStream();
	void Read(char*, int&);
	void Write(char*, int);
};


MATFile* matOpen(const char*, const char*);
int matClose(MATFile*);
mxArray* matGetVariable(MATFile*, const char*);
mxArray* matGetNextVariable(MATFile*, const char**);
int matPutVariable(MATFile*, char*, mxArray*);
char **matGetDir(MATFile *, int *);

void* mxCalloc(size_t,size_t);
void mxClearLogical(mxArray*);
mxArray* mxCreateLogicalMatrix(mwSize, mwSize);
mxArray* mxCreateNumericMatrix(mwSize, mwSize, mxClassID, mxComplexity);
mxArray* mxCreateNumericArray(mwSize, const mwSize*, mxClassID, mxComplexity);
mxArray* mxCreateString(const char*);
mxArray* mxCreateStructMatrix(mwSize, const mwSize *, int, const char**);
mxArray *mxCreateDoubleMatrix(mwSize m, mwSize n, mxComplexity ComplexFlag);
//mxArray* mxCreateStructMatrix(mwSize, mwSize, int, const char**);
void mxDestroyArray(mxArray*);
void* mxGetData(mxArray*);
mxArray* mxGetField(mxArray*, mwIndex, const char*);
mxLogical* mxGetLogicals(mxArray*);
const mwSize* mxGetDimensions(mxArray*);
mwSize *mxGetNumberOfElements(const mxArray *pm);
double* mxGetPr(mxArray*);
double* mxGetPi(mxArray*);
int mxGetString(mxArray*,char*,mwSize);
bool mxIsComplex(mxArray*);
bool mxIsLogical(mxArray*);
void mxSetField(mxArray*, mwIndex, const char*, mxArray*);
extern int mxAddField(mxArray *, const char *);
void mxRemoveField(mxArray *,int );
//char* mxGetName(mxArray*);
//void mxSetLogical(mxArray*);
//void mxSetName(mxArray*, char*);
int mxGetFieldNumber(const mxArray *, const char *);
int mxIsStruct(const mxArray *);
int mxIsDouble(const mxArray *);
int mxIsUint8(const mxArray *);
int mxGetNumberOfDimensions(const mxArray *);
int mxGetNumberOfFields(const mxArray *);
int mxIsChar(const mxArray *);
int mexPrintf(const char*);//evaluated directly
void mexErrMsgTxt(const char *);
int mexEvalString(const char*);//dummy
void mxFree(void *);
bool mxIsEmpty(const mxArray *array_ptr);
//---------------------------------------------------------------------------

#endif
