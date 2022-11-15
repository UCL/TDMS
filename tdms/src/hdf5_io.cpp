// own include
#include "hdf5_io.h"

// std
#include <string>
#include <iostream>
#include <stdexcept>

// external
#include <H5Cpp.h>
#include <H5Fpublic.h>
#include <spdlog/spdlog.h>


/**
 * @brief Convert a HDF5FileMode to the #defined HDF5 file creation property.
 * @note Internal function, only used in this file by HDF5File::_open.
 */
unsigned int convert_to_h5f_global(HDF5FileMode mode) {
  switch (mode) {
    case READONLY:  return H5F_ACC_RDONLY;
    case READWRITE: return H5F_ACC_RDWR;
    case OVERWRITE: return H5F_ACC_TRUNC;
    default: return std::numeric_limits<int>::max();
  }
}

void HDF5File::_open() {
    spdlog::trace("Opening file: {}, in mode: {}", filename_, static_cast<int>(mode_));
    file_ = new H5::H5File(filename_.c_str(), convert_to_h5f_global(mode_));
    return;
}

HDF5File::~HDF5File() {
    file_->close();
    for (auto ds: datasets_)
        delete ds;
    delete file_;
}

void HDF5File::write() {
    for (int i=0; i<3; i++) {
        std::cout << convert_to_h5f_global(static_cast<HDF5FileMode>(i)) << std::endl;
    }
    return;
}
void HDF5File::read() {
    return;
}

bool HDF5File::isOK(bool print_debug) {

    if (print_debug) {
      // debug information
      spdlog::debug("File is named: {}", filename_);
      spdlog::debug("File mode: {}", static_cast<int>(mode_));
      spdlog::debug("Internal H5::H5File address: {:p}", (void*)file_);
    }

    // tests here
    if (file_ == nullptr) return false;

    if (file_->getFileSize() <= 0) return false;

    // passed all tests: it's ok
    return true;
}


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



void example_hdf5() {

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
