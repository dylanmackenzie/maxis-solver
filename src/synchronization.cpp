#include <atomic>
#include <condition_variable>
#include <exception>
#include <mutex>

#include "maxis/synchronization.hpp"

namespace maxis {

WorkerSynchronizer::WorkerSynchronizer(
        unsigned int n
    ) : available_handles{n},
        utilized_handles{0},
        discarded_handles{0},
        current_cycle{0},
        is_terminated{false} {}

WorkerSynchronizer::Handle::Handle(WorkerSynchronizer &sync, unsigned int desired_cycle)
        : sync{sync} {

    std::unique_lock<std::mutex> l(sync.worker_lock);
    while(desired_cycle != sync.current_cycle) {
        sync.worker_cv.wait(l);
    }

    if (++sync.utilized_handles > sync.available_handles) {
        throw SynchronizationError("Too many workers attempting to reserve handles");
    }
}

WorkerSynchronizer::Handle::~Handle() {
    ++sync.discarded_handles;

    std::lock_guard<std::mutex> l(sync.manager_lock);
    sync.manager_cv.notify_all();
}

bool
WorkerSynchronizer::Handle::operator!() {
    return sync.is_terminated;
}

void
WorkerSynchronizer::wait_on_cycle() {
    std::unique_lock<std::mutex> l(manager_lock);
    while(discarded_handles != available_handles) {
        manager_cv.wait(l);
    }
}

void
WorkerSynchronizer::next_cycle() {
    utilized_handles = 0;
    discarded_handles = 0;
    ++current_cycle;

    std::lock_guard<std::mutex> l(worker_lock);
    worker_cv.notify_all();
}

void
WorkerSynchronizer::terminate() {
    is_terminated = true;
    next_cycle();
    wait_on_cycle();
}

}
