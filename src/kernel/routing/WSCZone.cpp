/* Copyright (c) 2009-2022. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <ios>
#include <omp.h>

#include <climits>
#include <iostream>
#include <rocksdb/compression_type.h>
#include <rocksdb/slice.h>
#include <rocksdb/status.h>
#include <simgrid/kernel/routing/NetPoint.hpp>
#include <simgrid/kernel/routing/WSCZone.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <xbt/string.hpp>

#include "rocksdb/db.h"
#include "rocksdb/options.h"

#include "src/kernel/resource/NetworkModel.hpp"

XBT_LOG_NEW_DEFAULT_SUBCATEGORY(ker_routing_wsc, ker_routing, "Kernel WSC Routing");

namespace simgrid {
namespace kernel::routing {

void WSCZone::init_tables(unsigned int table_size)
{
  if (link_table_.size() != table_size) {
    blocked_table_size = floyd_warshall_blocked_init(table_size, 256);
    /* Resize and initialize Cost, Predecessor and Link tables */
    link_table_.resize(table_size);
    for (auto& link : link_table_)
      link.resize(table_size); /* actual link between src and dst */
  }
}

void WSCZone::get_local_route(const NetPoint* src, const NetPoint* dst, Route* route, double* lat)
{
  get_route_check_params(src, dst);

  /* create a result route */
  std::vector<Route*> route_stack;
  unsigned long cur = dst->id();
  do {
    long pred = predecessor_table_[src->id() * blocked_table_size + cur];
    if (pred == -1)
      throw std::invalid_argument(xbt::string_printf("No route from '%s' to '%s'", src->get_cname(), dst->get_cname()));
    route_stack.push_back(link_table_[pred][cur].get());
    cur = pred;
  } while (cur != src->id());

  if (get_hierarchy() == RoutingMode::recursive) {
    route->gw_src_ = route_stack.back()->gw_src_;
    route->gw_dst_ = route_stack.front()->gw_dst_;
  }

  const NetPoint* prev_dst_gw = nullptr;
  while (not route_stack.empty()) {
    const Route* e_route = route_stack.back();
    route_stack.pop_back();
    if (get_hierarchy() == RoutingMode::recursive && prev_dst_gw != nullptr &&
        prev_dst_gw->get_cname() != e_route->gw_src_->get_cname()) {
      get_global_route(prev_dst_gw, e_route->gw_src_, route->link_list_, lat);
    }

    add_link_latency(route->link_list_, e_route->link_list_, lat);

    prev_dst_gw = e_route->gw_dst_;
  }
}

void WSCZone::add_route(NetPoint* src, NetPoint* dst, NetPoint* gw_src, NetPoint* gw_dst,
                        const std::vector<s4u::LinkInRoute>& link_list, bool symmetrical)
{
  /* set the size of table routing */
  auto table_size = get_table_size();
  init_tables(table_size);

  add_route_check_params(src, dst, gw_src, gw_dst, link_list, symmetrical);

  /* Check that the route does not already exist */
  if (gw_dst && gw_src) // netzone route (to adapt the error message, if any)
    xbt_assert(nullptr == link_table_[src->id()][dst->id()],
               "The route between %s@%s and %s@%s already exists (Rq: "
               "routes are symmetrical by default).",
               src->get_cname(), gw_src->get_cname(), dst->get_cname(), gw_dst->get_cname());
  else
    xbt_assert(nullptr == link_table_[src->id()][dst->id()],
               "The route between %s and %s already exists (Rq: routes are "
               "symmetrical by default).",
               src->get_cname(), dst->get_cname());

  link_table_[src->id()][dst->id()] = std::unique_ptr<Route>(
      new_extended_route(get_hierarchy(), gw_src, gw_dst, get_link_list_impl(link_list, false), true));
  predecessor_table_[src->id() * blocked_table_size + dst->id()] = src->id();
  cost_table_[src->id() * blocked_table_size + dst->id()]        = link_table_[src->id()][dst->id()]->link_list_.size();

  if (symmetrical) {
    if (gw_dst && gw_src) // netzone route (to adapt the error message, if any)
      xbt_assert(nullptr == link_table_[dst->id()][src->id()],
                 "The route between %s@%s and %s@%s already exists. You "
                 "should not declare the reverse path as symmetrical.",
                 dst->get_cname(), gw_dst->get_cname(), src->get_cname(), gw_src->get_cname());
    else
      xbt_assert(nullptr == link_table_[dst->id()][src->id()],
                 "The route between %s and %s already exists. You should "
                 "not declare the reverse path as symmetrical.",
                 dst->get_cname(), src->get_cname());

    if (gw_dst && gw_src) {
      NetPoint* gw_tmp = gw_src;
      gw_src           = gw_dst;
      gw_dst           = gw_tmp;
    }

    if (not gw_src || not gw_dst)
      XBT_DEBUG("Load Route from \"%s\" to \"%s\"", dst->get_cname(), src->get_cname());
    else
      XBT_DEBUG("Load NetzoneRoute from \"%s(%s)\" to \"%s(%s)\"", dst->get_cname(), gw_src->get_cname(),
                src->get_cname(), gw_dst->get_cname());

    link_table_[dst->id()][src->id()] = std::unique_ptr<Route>(
        new_extended_route(get_hierarchy(), gw_src, gw_dst, get_link_list_impl(link_list, true), false));
    predecessor_table_[dst->id() * blocked_table_size + src->id()] = dst->id();
    cost_table_[dst->id() * blocked_table_size + src->id()] =
        link_table_[dst->id()][src->id()]->link_list_.size(); /* count of links, old model assume 1 */
  }
}

bool WSCZone::init_ckpt_db(const char* db_path)
{

  rocksdb::Options options;
  // Using 16 threads by default
  options.IncreaseParallelism();
  options.OptimizeLevelStyleCompaction();
  options.create_if_missing = false;
  options.compression       = ROCKSDB_NAMESPACE::kLZ4Compression;

  rocksdb::Status s;
  s = rocksdb::DB::OpenForReadOnly(options, db_path, &ckpt_db);
  if(!s.ok() || !ckpt_db) {
    // Open without read-only mode
    options.create_if_missing = true;
    s = rocksdb::DB::Open(options, db_path, &ckpt_db);
  }
  xbt_assert(s.ok(), "Failed to open or create checkpoint database: %s", db_path);

  std::string has_routing_data;
  // Check if the database contains routing data
  s = ckpt_db->Get(rocksdb::ReadOptions(), "sealed", &has_routing_data);
  // We do not check status here since if read failed, we will overwrite the database.

  return (s.ok() && ~s.IsNotFound());
}

// aux functions for floyd routing
uint64_t WSCZone::floyd_warshall_blocked_init(const int n, const int block_size)
{
  [[maybe_unused]] int n_oversized{};
  int block_remainder{n % block_size};
  if (block_remainder == 0) {
    n_oversized = n;
  } else {
    n_oversized = n + block_size - block_remainder;
  }

  // Create the buffer if not existed
  cost_table_.resize(n_oversized * n_oversized, ULONG_MAX);
  predecessor_table_.resize(n_oversized * n_oversized, -1);
  return n_oversized;
}

template <typename COST_T, typename PRED_T>
inline void floyd_warshall_in_place(std::vector<COST_T>& output, const uint64_t c_idx, const uint64_t a_idx,
                                    const uint64_t b_idx, std::vector<PRED_T>& aux, const int ori_n, const int b,
                                    const int n, const int k_base, const int i_base, const int j_base,
                                    const int host_id, const int io_router_id)
{
  for (int k = 0; k < b; k++) {
    if (k + k_base >= ori_n)
      break;
    int kth = k * n;
    for (int i = 0; i < b; i++) {
      if (i + i_base >= ori_n)
        break;
      for (int j = 0; j < b; j++) {
        if (j + j_base >= ori_n)
          break;
        if (k + k_base == io_router_id && i + i_base != host_id && j + j_base != host_id &&
            i + i_base != io_router_id && j + j_base != io_router_id)
          continue;
        COST_T sum = output[a_idx + i * n + k] + output[b_idx + kth + j];
        if (output[a_idx + i * n + k] < ULONG_MAX && output[b_idx + kth + j] < ULONG_MAX &&
            (output[c_idx + i * n + j] > sum || output[c_idx + i * n + j] == ULONG_MAX)) {
          output[c_idx + i * n + j] = sum;

          // Casting is defined in zip_vector::ref
          aux[c_idx + i * n + j] = aux[b_idx + kth + j];
        }
      }
    }
  }
}

void WSCZone::floyd_warshall_blocked(const int ori_n, const int n, const int b, const int io_router_id,
                                     const int host_id)
{
  // for now, assume b divides n
  const int blocks = n / b;

  // note that [i][j] == [i * input_width * block_width + j * block_width]
  for (int k = 0; k < blocks; k++) {
    floyd_warshall_in_place<COST_T, PRED_T>(cost_table_, k * b * n + k * b, k * b * n + k * b, k * b * n + k * b,
                                            predecessor_table_, ori_n, b, n, k * b, k * b, k * b, host_id,
                                            io_router_id);
#pragma omp parallel for
    for (int j = 0; j < blocks; j++) {
      if (j == k)
        continue;
      floyd_warshall_in_place<COST_T, PRED_T>(cost_table_, k * b * n + j * b, k * b * n + k * b, k * b * n + j * b,
                                              predecessor_table_, ori_n, b, n, k * b, k * b, j * b, host_id,
                                              io_router_id);
    }
#pragma omp parallel for
    for (int i = 0; i < blocks; i++) {
      if (i == k)
        continue;
      floyd_warshall_in_place<COST_T, PRED_T>(cost_table_, i * b * n + k * b, i * b * n + k * b, k * b * n + k * b,
                                              predecessor_table_, ori_n, b, n, k * b, i * b, k * b, host_id,
                                              io_router_id);
      for (int j = 0; j < blocks; j++) {
        if (j == k)
          continue;
        floyd_warshall_in_place<COST_T, PRED_T>(cost_table_, i * b * n + j * b, i * b * n + k * b, k * b * n + j * b,
                                                predecessor_table_, ori_n, b, n, k * b, i * b, j * b, host_id,
                                                io_router_id);
      }
    }
  }
}

void WSCZone::do_seal()
{
  auto ckpt_path = this->get_property("ckpt_path");
  /* set the size of table routing */
  auto table_size = get_table_size();
  init_tables(table_size);

  std::cout << "Checkpoint Path in WSCZone: " << ckpt_path << std::endl;
  if (init_ckpt_db(ckpt_path)) {
    bool flag{true};
    // Load checkpoints from database
    for (int row = 0; row < table_size; row++) {
      rocksdb::PinnableSlice cost_row, pred_row;
      flag = true;

      // Prepare the key
      std::ostringstream key_cost;
      std::ostringstream key_pred;
      key_cost << "cost:" << row;
      key_pred << "pred:" << row;

      flag &= (ckpt_db->Get(rocksdb::ReadOptions(), ckpt_db->DefaultColumnFamily(), key_cost.str(), &cost_row)).ok();
      flag &= (ckpt_db->Get(rocksdb::ReadOptions(), ckpt_db->DefaultColumnFamily(), key_pred.str(), &pred_row)).ok();
      if (flag) {
        // Copy memory to cost_table_
        // UNSAFE OPERATION, BE VERY CATIOUS
        memcpy((void*)(cost_table_.data() + row * blocked_table_size), (const void*)cost_row.data_,
               table_size * sizeof(COST_T));
        memcpy((void*)(predecessor_table_.data() + row * blocked_table_size), (const void*)pred_row.data_,
               table_size * sizeof(PRED_T));
        cost_row.Reset();
        pred_row.Reset();
      } else {
        break;
      }
    }

    if (flag) {
      return;
    }
  }

  std::cout << "Routing Checkpoint not found! Generating now." << std::endl;

  unsigned host_id{0};
  unsigned io_router_id{table_size - 1};

  /* Add the loopback if needed */
  if (get_network_model()->loopback_ && get_hierarchy() == RoutingMode::base) {
    for (unsigned int i = 0; i < table_size; i++) {
      auto& route = link_table_[i][i];
      if (not route) {
        route.reset(new Route());
        route->link_list_.push_back(get_network_model()->loopback_.get());
        predecessor_table_[i * blocked_table_size + i] = i;
        cost_table_[i * blocked_table_size + i]        = 1;
      }
    }
  }

  // Blocked Parallel WSC
  // Reference: https://moorejs.github.io/APSP-in-parallel/
  floyd_warshall_blocked(table_size, blocked_table_size, 256, io_router_id, host_id);

  // Save checkpoints to database
  for (int row = 0; row < table_size; row++) {
    // Prepare the key
    std::ostringstream key_cost;
    std::ostringstream key_pred;
    key_cost << "cost:" << row;
    key_pred << "pred:" << row;

    std::string_view cost_row((char*)(cost_table_.data() + row * blocked_table_size), table_size * sizeof(COST_T));
    std::string_view pred_row((char*)(predecessor_table_.data() + row * blocked_table_size),
                              table_size * sizeof(PRED_T));

    auto s = ckpt_db->Put(rocksdb::WriteOptions(), key_cost.str(), cost_row);
    xbt_assert(s.ok(), "Failed to write checkpoint to database %s at row %d", ckpt_path, row);
    s = ckpt_db->Put(rocksdb::WriteOptions(), key_pred.str(), pred_row);
    xbt_assert(s.ok(), "Failed to write checkpoint to database %s at row %d", ckpt_path, row);
  }
  // Seal the database and add metadata
  auto s = ckpt_db->Put(rocksdb::WriteOptions(), "table_size", std::to_string(table_size));
  s      = ckpt_db->Put(rocksdb::WriteOptions(), "sealed", "1");
  xbt_assert(s.ok(), "Failed to seal the database. Checkpoint is corrupted.");
  // Sync WAL before closing
  ckpt_db->Flush(rocksdb::FlushOptions());
  ckpt_db->Close();
}
} // namespace kernel::routing

namespace s4u {
NetZone* create_wsc_zone(const std::string& name)
{
  return (new kernel::routing::WSCZone(name))->get_iface();
}
} // namespace s4u

} // namespace simgrid
