
class Timer{

    double start_time; // in seconds
    double end_time;

public:
    void start();
    void end();
    void click();
    double delta_seconds() const;
};
