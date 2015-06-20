#include <atomic>
#include <condition_variable>
#include <exception>
#include <mutex>

#ifndef MAXIS_SYNCHRONIZATION_H
#define MAXIS_SYNCHRONIZATION_H

namespace maxis {

class WorkerSynchronizer {
public:
    // The number of workers that will be executing is passed to the
    // constructor. The number of threads attempting to use the
    // synchronizer SHOULD BE EXACTLY THIS NUMBER.
    WorkerSynchronizer(unsigned int);

    // wait_on_cycle blocks until all workers are complete. Once it
    // returns, no new handles can be created until next_cycle is
    // called.
    void wait_on_cycle();

    // next_cycle is called to start a new work cycle.
    void next_cycle();

    // terminate signals to the workers that this cycle will be the
    // last, and that no actual work should be done, only clean up.
    // It blocks until all workers terminate.
    void terminate();

    // This class allows us to use RAII to synchronize the workers. A
    // worker constructs a handle for the current work cycle and
    // destroys it when it is done with that work cycle. The worker then
    // increments the number of cycles it has completed, and constructs
    // a handle for the next work cycle. That constructor will block
    // until the remaining threads complete the current work cycle and
    // next_cycle is called.
    class Handle {
    public:
        Handle(WorkerSynchronizer&, unsigned int);
        ~Handle();
        Handle(const Handle&) =delete;
        Handle& operator=(const Handle&) =delete;
        // Move constructors will be disabled automatically.

        // Returns true if this is the cleanup cycle.
        bool operator!();

    private:
        WorkerSynchronizer &sync;
    };

    struct SynchronizationError : public std::runtime_error {
        SynchronizationError(const char* s) : std::runtime_error(s) {};
    };

private:
    // available_handles is the total number of handles available for a
    // given work cycle
    const unsigned int available_handles;
    std::atomic<unsigned int> utilized_handles;  // number of successfully constructed handles
    std::atomic<unsigned int> discarded_handles; // number of successfully destructed handles
    std::atomic<unsigned int> current_cycle;
    std::atomic<bool> is_terminated;

    std::condition_variable worker_cv;
    std::mutex worker_lock;
    std::condition_variable manager_cv;
    std::mutex manager_lock;
};

} // namespace maxis

#endif // MAXIS_SYNCHRONIZATION_H
