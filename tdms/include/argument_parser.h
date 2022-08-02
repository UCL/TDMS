#include "string"
#include "vector"


class ArgumentNamespace{

private:
    std::vector<std::string> arguments;
    std::vector<std::string> non_flag_arguments;

public:
    int num_non_flag = 0;

    explicit ArgumentNamespace(int nargs, char *argv[]);
    bool have_flag(std::string const &flag);
    bool has_grid_filename() const;
    bool have_correct_number_of_filenames() const;

    static bool is_a_flag_argument(std::string arg);

    const char* input_filename();
    const char* output_filename();
    const char* grid_filename();

    std::vector<std::string> input_filenames();
};


class ArgumentParser {

private:
    static void print_help_message();

public:
    explicit ArgumentParser();
    static ArgumentNamespace parse_args(int nargs, char *argv[]);
};
