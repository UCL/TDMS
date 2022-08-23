
class Timer{

    double start_time; // in seconds
    double end_time;

public:
    void start();

    void end();

    /**
     * Log the difference in time and reset the timer
     */
    void click();

    double delta_seconds() const;
};
