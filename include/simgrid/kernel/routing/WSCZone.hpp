/* Copyright (c) 2013-2022. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef SURF_ROUTING_WSC_HPP_
#define SURF_ROUTING_WSC_HPP_

#include <simgrid/kernel/routing/RoutedZone.hpp>

#define OPENMP_THREADS 16

namespace simgrid {
namespace kernel {
namespace routing {

/** @ingroup ROUTING_API
 *  @brief NetZone with an explicit routing computed at initialization with Wafer-Scale Chip format
 *
 *  The path between components is computed at creation time from every one-hop links
 *
 *  This result in rather small platform file, slow initialization time,  and intermediate memory requirements
 *  (somewhere between the one of @{FloydZone} and the one of @{FullZone}).
 */
class XBT_PRIVATE WSCZone : public RoutedZone {
  /* vars to compute the WSC algorithm. */
  std::vector<std::vector<long>> predecessor_table_;
  std::vector<std::vector<std::unique_ptr<Route>>> link_table_;
  unsigned contiguous_x = 1; // First PE has ID of 1 due to host claiming ID 0
  unsigned wsc_x = 0;
  unsigned wsc_y = 0;
  bool dimensions_inferred = false;

  void init_tables(unsigned int table_size);
  void do_seal() override;

public:
  using RoutedZone::RoutedZone;
  WSCZone(const WSCZone&) = delete;
  WSCZone& operator=(const WSCZone&) = delete;

  void get_local_route(const NetPoint* src, const NetPoint* dst, Route* into, double* latency) override;
  void add_route(NetPoint* src, NetPoint* dst, NetPoint* gw_src, NetPoint* gw_dst,
                 const std::vector<s4u::LinkInRoute>& link_list, bool symmetrical) override;
};
} // namespace routing
} // namespace kernel
} // namespace simgrid

#endif /* SURF_ROUTING_WSC_HPP_ */
