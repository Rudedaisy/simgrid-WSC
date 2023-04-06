/* Copyright (c) 2013-2022. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef SURF_ROUTING_WSC_HPP_
#define SURF_ROUTING_WSC_HPP_

#include <cstdint>
#include <rocksdb/db.h>
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

  using COST_T = uint64_t; // Cost Table Type
  using PRED_T = int64_t;  // Predecessor Table Type

  /* vars to compute the WSC algorithm. */
  std::vector<COST_T> cost_table_;
  std::vector<PRED_T> predecessor_table_;
  std::vector<std::vector<std::unique_ptr<Route>>> link_table_;

  void init_tables(unsigned int table_size);
  void do_seal() override;

public:
  using RoutedZone::RoutedZone;
  WSCZone(const WSCZone&)            = delete;
  WSCZone& operator=(const WSCZone&) = delete;

  void get_local_route(const NetPoint* src, const NetPoint* dst, Route* into, double* latency) override;
  void add_route(NetPoint* src, NetPoint* dst, NetPoint* gw_src, NetPoint* gw_dst,
                 const std::vector<s4u::LinkInRoute>& link_list, bool symmetrical) override;

private:
  rocksdb::DB* ckpt_db{nullptr}; // Database to store checkpoint
  uint64_t table_size{0};
  uint64_t blocked_table_size{0};

  bool init_ckpt_db(const char* db_path);
  void floyd_warshall_blocked(const int ori_n, const int n, const int b, const int io_router_id, const int host_id);
  uint64_t floyd_warshall_blocked_init(const int n, const int block_size);
};
} // namespace routing
} // namespace kernel
} // namespace simgrid

#endif /* SURF_ROUTING_WSC_HPP_ */
