/* Copyright (c) 2021. The SimGrid Team.
 * All rights reserved.                                                     */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef SRC_KERNEL_TIMER_TIMER_HPP_
#define SRC_KERNEL_TIMER_TIMER_HPP_

#include <simgrid/forward.h>
#include <xbt/functional.hpp>
#include <xbt/utility.hpp>

#include <boost/heap/fibonacci_heap.hpp>

namespace simgrid {
namespace kernel {
namespace timer {

inline auto& kernel_timers() // avoid static initialization order fiasco
{
  using TimerQelt = std::pair<double, Timer*>;
  static boost::heap::fibonacci_heap<TimerQelt, boost::heap::compare<xbt::HeapComparator<TimerQelt>>> value;
  return value;
}

/** @brief Timer datatype */
class Timer {
public:
  const double date;
  std::remove_reference_t<decltype(kernel_timers())>::handle_type handle_;

  Timer(double date, xbt::Task<void()>&& callback) : date(date), callback(std::move(callback)) {}

  xbt::Task<void()> callback;
  void remove();

  template <class F> static inline Timer* set(double date, F callback)
  {
    return set(date, xbt::Task<void()>(std::move(callback)));
  }

  static Timer* set(double date, xbt::Task<void()>&& callback);
  static double next() { return kernel_timers().empty() ? -1.0 : kernel_timers().top().first; }
};

} // namespace timer
} // namespace kernel
} // namespace simgrid

#endif /* SRC_KERNEL_TIMER_TIMER_HPP_ */