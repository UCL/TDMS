// own include
#include "hdf5_io.h"

// std
#include <string>

// external libraries
#include <H5Cpp.h>

// own headers


#define MAX_NAME_LENGTH 32
const std::string FileName("SimpleCompound.h5");
const std::string DatasetName("PersonalInformation");
const std::string member_age("Age");
const std::string member_sex("Sex");
const std::string member_name("Name");
const std::string member_height("Height");

typedef struct {
    int age;
    char sex;
    char name[MAX_NAME_LENGTH];
    float height;
} PersonalInformation;



void test_hdf5() {

    // Data to write
    PersonalInformation person_list[] = {
        { 18, 'M', "Mary",  152.0   },
        { 32, 'F', "Tom",   178.6   },
        { 29, 'M', "Tarou", 166.6   }
    };
    // the length of the data
    //int length = sizeof(person_list) / sizeof(PersonalInformation);
    // the array of each length of multidimentional data.
    hsize_t dim[1];
    dim[0] = sizeof(person_list) / sizeof(PersonalInformation);

    // the length of dim
    int rank = sizeof(dim) / sizeof(hsize_t);

    // defining the datatype to pass HDF55
    H5::CompType mtype(sizeof(PersonalInformation));
    mtype.insertMember(member_age, HOFFSET(PersonalInformation, age), H5::PredType::NATIVE_INT);
    mtype.insertMember(member_sex, HOFFSET(PersonalInformation, sex), H5::PredType::C_S1);
    mtype.insertMember(member_name, HOFFSET(PersonalInformation, name), H5::StrType(H5::PredType::C_S1, MAX_NAME_LENGTH));
    mtype.insertMember(member_height, HOFFSET(PersonalInformation, height), H5::PredType::NATIVE_FLOAT);

    // preparation of a dataset and a file.
    H5::DataSpace space(rank, dim);
    H5::H5File *file = new H5::H5File(FileName, H5F_ACC_TRUNC);
    H5::DataSet *dataset = new H5::DataSet(file->createDataSet(DatasetName, mtype, space));
    // Write
    dataset->write(person_list, mtype);

    delete dataset;
    delete file;
    return;

}
