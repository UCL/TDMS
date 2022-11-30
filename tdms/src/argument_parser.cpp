#include "argument_parser.h"

#include <stdexcept>

#include <spdlog/spdlog.h>

using namespace std;


ArgumentNamespace ArgumentParser::parse_args(int n_args, char *arg_ptrs[]) {
  spdlog::debug("Parsing {} command line arguments", n_args);

  auto args = ArgumentNamespace(n_args, arg_ptrs);

  if (args.have_flag("-h")){
    print_help_message();
    exit(0);
  }

  if (args.have_flag("-q")){  // quiet operation
    spdlog::set_level(spdlog::level::off);
  }

  if (!args.have_correct_number_of_filenames()){
    fprintf(stderr,"Incorrect number of arguments. See below for help\n\n");
    print_help_message();
    exit(-1);
  }

  if( !args.has_grid_filename()){
    fprintf(stdout,"infile:[%s], outfile:[%s], m=%d\n",
            args.input_filename(),
            args.output_filename(),
            args.have_flag("-m"));
  }
  else{
    fprintf(stdout,"infile:[%s], gridfile:[%s], outfile:[%s], m=%d\n",
            args.input_filename(),
            args.grid_filename(),
            args.output_filename(),
            args.have_flag("-m"));
  }

  spdlog::debug("Finished parsing arguments");
  return args;
}

void ArgumentParser::print_help_message(){
  fprintf(stdout,"Usage:\n"
                 "openandorder [options] infile outfile\n"
                 "openandorder [options] infile gridfile outfile\n"
                 "Options:\n"
                 "-h:\tDisplay this help message\n"
                 "--finite-difference:\tUse the finite-difference solver, instead of the pseudo-spectral method.\n"
                 "-q:\tQuiet operation. Silence all logging\n"
                 "-m:\tMinimise output file size by not saving vertex and facet information\n\n");
}

ArgumentNamespace::ArgumentNamespace(int n_args, char *arg_ptrs[]) {

  arguments = std::vector<std::string>(arg_ptrs + 1, arg_ptrs + n_args);

  for (const auto& arg: arguments){

    if (!is_a_flag_argument(arg)){
      non_flag_arguments.push_back(arg);
      num_non_flag++;
    }
  }

}

bool ArgumentNamespace::is_a_flag_argument(std::string arg){
  return arg[0] == '-';
}

bool ArgumentNamespace::have_flag(std::string const &flag) const {

  for (const auto& arg : arguments){
    if(arg == flag) return true;
  }

  return false;
}

bool ArgumentNamespace::finite_difference() const {
  return this->have_flag("-fd") || this->have_flag("--finite-difference");
}

bool ArgumentNamespace::no_band_limited_schemes() const {
  return this->have_flag("-nbli") || this->have_flag("--no-band-limited");
}

const char* ArgumentNamespace::output_filename() {

  if (has_grid_filename()){
    return non_flag_arguments[2].c_str();
  } else if (num_non_flag == 2){
    return non_flag_arguments[1].c_str();
  }

  throw std::runtime_error("Failed to determine the output file from arguments");
}

const char* ArgumentNamespace::grid_filename() {

  if (num_non_flag == 3){
    return non_flag_arguments[1].c_str();
  }

  throw std::runtime_error("Failed to determine the grid file from arguments");
}

const char* ArgumentNamespace::input_filename() {

  if (num_non_flag > 0){
    return non_flag_arguments[0].c_str();
  }

  throw std::runtime_error("Failed to determine the input file from arguments");
}

bool ArgumentNamespace::has_grid_filename() const {
  return num_non_flag == 3;
}

vector<string> ArgumentNamespace::input_filenames() {

  vector<string> filenames;

  filenames.emplace_back(input_filename());
  if (has_grid_filename()){
    filenames.emplace_back(grid_filename());
  }

  return filenames;
}

bool ArgumentNamespace::have_correct_number_of_filenames() const {
  return num_non_flag == 2 || num_non_flag == 3;
}
