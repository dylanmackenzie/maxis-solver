#include <atomic>
#include <condition_variable>
#include <exception>
#include <mutex>

#ifndef MAXIS_SYNCHRONIZATION_H
#define MAXIS_SYNCHRONIZATION_H

namespace maxis {

// This class is designed to facilitate a single manager thread sharing
// memory with many worker threads, each of which are responsible for
// contiguous, non-overlapping slices of that memory. It was implemented
// before I was aware of boost::shared_mutex, and instead uses atomics
// in combination with condition variables to implement what is
// basically a pair of semaphores (although they act are used slightly
// abnormally, see below).
//
// One semaphore keeps track of the number of initiated workers, and the
// other keeps track of the completed workers. Additionally, the
// synchronizer keeps an incrementing count of how many times a full set
// of workers has run to completion, and each worker keeps a count of how
// many cycles it has completed.
//
// To begin, a worker constructs a handle by passing a reference to the
// synchronizer by which it is managed and the number of work cycles it
// has completed. The handle constructor waits until the number of work
// cycles completed by the synchronizer is equal to the number of work
// cycles completed by the current worker. Once this is true, the
// constructor increments the initiated workers semaphore and allows the
// worker to execute. If incrementing the initiated workers semaphore
// results in it being greater than the max number of workers, the
// handle constructor throws.
//
// Once a worker finishes executing the current cycle, it lets its
// handle go out of scope. The destructor of the handle increments the
// completed workers semaphore. Notice that we do not decrement the
// initiated workers semaphore, otherwise accidentally creating too many
// handles for the same work cycle might not throw an error. This differs
// from the typical use of a semaphore.
//
// Once all workers finish running, the completed workers semaphore is
// equal to the max number of workers. When the manager sees this, it
// resumes with the knowledge that all of the worker threads are waiting
// on a handle for the next work cycle. It runs whatever task is
// required, then sets both semaphores back to 0. Finally, it increments
// the work cycle count of the synchronizer, allowing all of the workers
// to begin the next cycle.
//
// Because of this somewhat convoluted approach to concurrency, helgrind
// cannot meaningfully analyze it. Additionally, since we now lock the
// mutex for the condition variable whenever a worker or manager thread
// finishes a work cycle, it's probably not much faster than the
// shared_mutex implementation.
//
// TODO: use boost::shared_mutex instead of atomic "semaphores"
class WorkerSynchronizer {
public:
    // The number of workers that will be executing is passed to the
    // constructor. The number of threads attempting to use the
    // synchronizer SHOULD BE EXACTLY THIS NUMBER. The first work cycle
    // occurs as soon as a worker thread is spawned.
    WorkerSynchronizer(unsigned int);

    // wait_on_cycle blocks until all workers are complete. Once it
    // returns, no workers will execute until next_cycle is called.
    void wait_on_cycle();

    // next_cycle is called to start a new work cycle.
    void next_cycle();

    // terminate signals to the workers that this cycle will be the
    // last, and that no actual work should be done, only clean up.
    // It blocks until all workers terminate.
    void terminate();

    // This class allows us to use RAII to synchronize the workers. The
    // worker should check if !handle == true before it begins working
    // to see if is in the clean up cycle.
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
